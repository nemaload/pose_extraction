#define SAMPLE_SIZE 10000
#define RANDOM_SEED

#include <time.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>

#include <argtable2.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

int parse_args(int argc, char** argv, char** filename, int* width, int* height) {
  struct arg_str *file = arg_str1(NULL, NULL, "<file>", "the input image (raw data)");
  struct arg_int *w = arg_int0("w", "width", "<int>", "the width of each image slice");
  struct arg_int *h = arg_int0("h", "height", "<int>", "the height of each image slice");
  struct arg_end *end = arg_end(10);
  void* argtable[] = {file, w, h, end};
  int nerrors;

  w->ival[0] = h->ival[0] = -1;

  nerrors = arg_parse(argc, argv, argtable);
  if(nerrors > 0) {
    arg_print_errors(stdout, end, argv[0]);
    arg_print_syntax(stdout, argtable, "\n");
    return nerrors;
  } else {
    *filename = file->sval[0];
    *width = w->ival[0];
    *height = h->ival[0];
    return 0;
  }
}

void* open_mmapped_file(char* filename, int* length) {
  struct stat fs;
  int fd;
  void* region;

  if(stat(filename, &fs)) {
    perror("cannot read file");
    return NULL;
  }
  *length = fs.st_size;
  
  fd = open(filename, O_RDONLY);
  region=mmap(NULL, *length, PROT_READ, MAP_SHARED, fd, 0);
  //madvise(region, *length, POSIX_MADV_RANDOM);
  return region;
}

void compute_sd(void* data, int length, double* mean, double* sd) {
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mrg);
#ifdef RANDOM_SEED
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &t);
  gsl_rng_set(r, t.tv_nsec);
#else
  gsl_rng_set(r, 0);
#endif
  unsigned short sample[SAMPLE_SIZE];
  int i;
#ifdef DEBUG
  printf("choosing %d samples from %d pixels...\n", SAMPLE_SIZE, length/2);
#endif
  gsl_ran_sample(r, sample, SAMPLE_SIZE, data, length/2, 2);
#ifdef DEBUG
  for(i=0; i<SAMPLE_SIZE; i++) {
    printf("Sample %d: %hu\n", i, sample[i]);
  }
#endif
  *mean = gsl_stats_ushort_mean(sample, 1, SAMPLE_SIZE);
  *sd = gsl_stats_ushort_sd_m(sample, 1, SAMPLE_SIZE, *mean);
}

int main(int argc, char** argv) {
  char* filename;
  int width;
  int height;
  int length;
  void* data;
  double mean;
  double sd;

  parse_args(argc, argv, &filename, &width, &height);
  data=open_mmapped_file(filename, &length);
  if(data==NULL) {
    perror("mmap failed");
    return 1;
  }

  printf("File mmapped successfully...\n");
  
  compute_sd(data, length, &mean, &sd);
  printf("The standard deviation of %d randomly chosen points is: %lf\nThe mean is: %lf\n", SAMPLE_SIZE, sd, mean);

  return 0;
}
