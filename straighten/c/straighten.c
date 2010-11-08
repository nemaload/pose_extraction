#define SD_SAMPLE_SIZE_DEFAULT 1000
#define MST_SAMPLE_SIZE_DEFAULT 150
#define PROGRESS_BAR_WIDTH 100

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>

#include <argtable2.h>

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include "image.h"

/*
 * This struct contains the command-line parameters.
 */
typedef struct {
  const char* filename;
  int width;
  int height;
  int sd_sample_size;
  int mst_sample_size;
} params_t;

/*
 * We use the argtable library to parse command-line arguments.
 */
int parse_args(int argc, char** argv, params_t* params) {
  struct arg_str* file = arg_str1(NULL, NULL, "<file>", "the input image (raw data)");
  struct arg_int* w = arg_int1("w", "width", "<int>", "the width of each image slice");
  struct arg_int* h = arg_int1("h", "height", "<int>", "the height of each image slice");
  struct arg_int* sdss = arg_int0(NULL, "sdss", "<int>", "the sample size to use for computing standard deviation");
  struct arg_int* mstss = arg_int0(NULL, "mstss", "<int>", "the sample size to use for making the MST");
  struct arg_end* end = arg_end(10);
  void* argtable[] = {file, w, h, sdss, mstss, end};
  int nerrors;

  w->ival[0] = h->ival[0] = -1;
  sdss->ival[0] = SD_SAMPLE_SIZE_DEFAULT;
  mstss->ival[0] = MST_SAMPLE_SIZE_DEFAULT;

  nerrors = arg_parse(argc, argv, argtable);
  if(nerrors > 0) {
    arg_print_errors(stderr, end, argv[0]);
    fprintf(stderr,"\nUsage:\n%s",argv[0]);
    arg_print_syntaxv(stderr, argtable, "\n\n");
    arg_print_glossary(stderr, argtable, "\t%-25s %s\n");
    fprintf(stderr,"\n");
    exit(nerrors);
  } else {
    params->filename = file->sval[0];
    params->width = w->ival[0];
    params->height = h->ival[0];
    params->sd_sample_size = sdss->ival[0];
    params->mst_sample_size = mstss->ival[0];
    return 0;
  }
}

/*
 * This function is for a primitive sort of profiling.
 */
double elapsed(struct timespec e_t, struct timespec s_t);
typedef enum {START, END} se_t;
void step(se_t se, const char* desc) {
  static int n = -1;
  static struct timespec s_t_rt;
  static struct timespec s_t_usr;
  static double total_time_rt;
  static double total_time_usr;
  static const char* last_desc;
  if(se == START) {
    n++;
    clock_gettime(CLOCK_REALTIME, &s_t_rt);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &s_t_usr);
    if(desc) {
      printf("[Step %d begin (%s)]\n",n,desc);
      last_desc = desc;
    } else {
      printf("[Step %d begin]\n",n);
      last_desc = NULL;
    }
  } else {
    struct timespec e_t_rt;
    struct timespec e_t_usr;
    double elapsed_rt;
    double elapsed_usr;
    clock_gettime(CLOCK_REALTIME, &e_t_rt);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &e_t_usr);
    if(last_desc && !desc)
      printf("[Step %d end (%s)]\n\n",n,last_desc);
    else if(desc)
      printf("[Step %d end (%s)]\n\n",n,desc);
    else
      printf("[Step %d end]\n\n",n);
    elapsed_rt = elapsed(e_t_rt,s_t_rt);
    elapsed_usr = elapsed(e_t_usr,s_t_usr);

    printf("[Time spent in step %d: %lfms/%lfms]\n",n,elapsed_usr,elapsed_rt);
    total_time_rt += elapsed_rt;
    total_time_usr += elapsed_usr;
    printf("[Running total of time: %lfms/%lfms]\n\n",total_time_usr,total_time_rt);
  }
}
double elapsed(struct timespec e_t, struct timespec s_t) {
    return (e_t.tv_sec-s_t.tv_sec)*1e3+((e_t.tv_nsec-s_t.tv_nsec)*1e-6);
}

/*
 * This function is for displaying progress bars.
 */
void progress(unsigned int i, unsigned int n) {
#ifdef PROGRESS_BAR_WIDTH
  int k;
  if(i<=1) {
    printf("\n[>");
    for(k=0;k<PROGRESS_BAR_WIDTH;k++) {
      printf(" ");
    }
    printf("]");
    fflush(stdout);
  }
  if ((i-1)%(n/PROGRESS_BAR_WIDTH) == 0) {
    if(i/((double)n/PROGRESS_BAR_WIDTH)<=PROGRESS_BAR_WIDTH){
      printf("\033[%dG=>\033[%dG",2+(int)(int)((i-1)/((double)n/PROGRESS_BAR_WIDTH)),PROGRESS_BAR_WIDTH+4);
      printf("%9d/%9d", i, n);
      fflush(stdout);
    }
  }
  if(i>=n) {
    printf("\033[%dG Done!\033[0K\n\n",PROGRESS_BAR_WIDTH+4);
    fflush(stdout);
  }
#endif
}

/*
 * Bundle all the little initializations together
 */
void init(image_t* image, const params_t* params) {
  init_rng(image);
  image->width = params->width;
  image->height = params->height;
  compute_depth(image);
}

/*
 * [Step 0]
 * Here's where the action begins - step 0 of the algorithm, loading
 * the file - well, not loading it per se, but mmaping it.
 */
void* open_mmapped_file(const char* filename, int* length) {
  struct stat fs;
  int fd;
  void* region;

  //First we stat the file to get its length.
  if(stat(filename, &fs)) {
    perror("cannot read file");
    return NULL;
  }
  *length = fs.st_size;
  
  //Now get a file descriptor and mmap!
  fd = open(filename, O_RDONLY);
  region=mmap(NULL, *length, PROT_READ, MAP_SHARED, fd, 0);

  return region;
}

/*
 * [Step 1]
 * Here we compute the mean and standard deviation of the data;
 * but to avoid actually reading in all that data, we take a random
 * sample first.
 */
void compute_sd(const image_t* image, int sample_size, double* mean, double* sd) {
  unsigned short* sample;
  int i;

  sample = malloc(sizeof(unsigned short)*sample_size);
#ifdef DEBUG
  printf("choosing %d samples from %d pixels...\n", sample_size, image->length/2);
#endif

  gsl_ran_sample(image->r, sample, sample_size, image->data, image->length/2, 2);

#ifdef DEBUG
  for(i=0; i<sample_size; i++) {
    printf("Sample %d: %hu\n", i, sample[i]);
  }
#endif

  *mean = gsl_stats_ushort_mean(sample, 1, sample_size);
  *sd = gsl_stats_ushort_sd_m(sample, 1, sample_size, *mean);
}

/*
 * [Step 2]
 * Given a threshhold value, randomly sample points, adding them to
 * a point list if they are above the threshhold, until we have a
 * specified number of points.
 */

point_list_t* sample_bright_points(const image_t* image, double threshhold, int n) {
  point_list_t* list = NULL;
  point_t p;
  int i=0;
  while(i<n) {
    p = random_point(image);
    if(pixel_get(image,p) > threshhold) {
      i++;
      add_point_to_list(&list,p);
      progress(i,n);
    }
  }
  return list;
}

/*
 * [Step 3]
 * Compute the distances between all the points sampled in Step 2.
 */

double* compute_distances(point_list_t* list, int n) {
  point_list_t *i,*j;
  int ix=0;
  int n2=n*n;
  double* table = malloc(n2*sizeof(double));
  for(i=list; i!=NULL; i=i->n) {
    for(j=list; j!=NULL; j=j->n) {
      table[ix++] = distance(i->p,j->p);
      progress(ix,n2);
    }
  }
}

/*
 * Now, we put it all together!
 */
int main(int argc, char** argv) {
  params_t params;
  image_t image;
  double mean;
  double sd;
  double threshhold;

  parse_args(argc, argv, &params);

  printf("Worm straightener v0.0.1\ncreated by David Dalrymple\n============================\n\nLet's straighten this worm!\n\n");

  step(START, "mmapping file");
  image.data=open_mmapped_file(params.filename, &image.length);
  if(image.data==NULL) {
    perror("mmap failed");
    return 1;
  }
  printf("File mmapped successfully...\n");
  init(&image, &params);
  step(END, NULL);
  
  step(START, "computing mean & s.d.");
  compute_sd(&image, params.sd_sample_size, &mean, &sd);
  printf("The standard deviation of %d randomly chosen points is: %lf\nThe mean is: %lf\n", params.sd_sample_size, sd, mean);
  threshhold = mean + sd;
  printf("The threshhold is: %lf\n", threshhold);
  step(END, NULL);

  step(START, "sampling for MST");
  point_list_t* w;
  w = sample_bright_points(&image, threshhold, params.mst_sample_size);
  step(END, NULL);

  step(START, "computing distances for MST");
  double* distances;
  distances = compute_distances(w, params.mst_sample_size);
  step(END, NULL);

  return 0;
}
