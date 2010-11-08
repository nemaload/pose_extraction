#define SAMPLE_SIZE_DEFAULT 1000
#define RANDOM_SEED

#include <time.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>

#include <argtable2.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

/*
 * This struct is passed around by most of the functions,
 * and besides containing a pointer to the region of address space
 * where the input data is found, it also carries housekeeping
 * information like the width, height and depth (number of slices)
 * of the image, the total length of the image, and even more
 * questionable variables like a pointer to a random number
 * generator (so we only have to initialize it once).
 */
typedef struct {
  int width;
  int height;
  int depth;
  int length;
  void* data;
  gsl_rng* r;
} image_t;

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
  struct arg_end* end = arg_end(10);
  void* argtable[] = {file, w, h, sdss, end};
  int nerrors;

  w->ival[0] = h->ival[0] = -1;
  sdss->ival[0] = SAMPLE_SIZE_DEFAULT;

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
    return 0;
  }
}

/*
 * This function is for a primitive sort of profiling.
 */
typedef enum {START, END} se_t;
void step(se_t se, const char* desc) {
  static int n = -1;
  static struct timespec s_t;
  static const char* last_desc;
  if(se == START) {
    n++;
    clock_gettime(CLOCK_REALTIME, &s_t);
    if(desc) {
      printf("[Step %d begin (%s)]\n",n,desc);
      last_desc = desc;
    } else {
      printf("[Step %d begin]\n",n);
      last_desc = NULL;
    }
  } else {
    struct timespec e_t;
    double elapsed;
    clock_gettime(CLOCK_REALTIME, &e_t);
    elapsed = (e_t.tv_sec-s_t.tv_sec)+((e_t.tv_nsec-s_t.tv_nsec)*1e-9);
    if(last_desc && !desc)
      printf("[Step %d end (%s)]\n\n",n,last_desc);
    else if(desc)
      printf("[Step %d end (%s)]\n\n",n,desc);
    else
      printf("[Step %d end]\n\n",n);
    printf("[Time spent in step %d: %lf]\n\n",n,elapsed);
  }
}

/*
 * Taking care of that pesky random-number generator.
 * There's a compile-time flag here to use a random seed,
 * or to use the 0 seed for reproducible results.
 */
void init_rng(image_t* image) {
  image->r = gsl_rng_alloc(gsl_rng_mrg);
#ifdef RANDOM_SEED
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  gsl_rng_set(image->r, t.tv_nsec);
#else
  gsl_rng_set(image->r, 0);
#endif
}

/*
 * Computing the depth of the image given the width,
 * height and length.
 */
void compute_depth(image_t* i) {
  i->depth = i->length / i->width / i->height;
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
 * Now, we put it all together!
 */
int main(int argc, char** argv) {
  params_t params;
  image_t image;
  double mean;
  double sd;
  double threshhold;

  parse_args(argc, argv, &params);

  printf("Worm straightener v0.0.1\ncreated by David Dalrymple\n=========================\n\nLet's straighten this worm!\n\n");

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

  return 0;
}
