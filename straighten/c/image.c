#include "image.h"
#include <math.h>
#include <time.h>
#include <unistd.h>
#define RANDOM_SEED

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
  i->depth = i->length / i->width / i->height / sizeof(unsigned short);
}

point_t random_point(const image_t* i) {
  point_t p;
  p.x = gsl_rng_uniform_int(i->r,i->width);
  p.y = gsl_rng_uniform_int(i->r,i->height);
  p.z = gsl_rng_uniform_int(i->r,i->depth);
  return p;
}

unsigned short pixel_get(const image_t* i, point_t p) {
  return ((unsigned short*)i->data)[p.z*i->width*i->height+p.y*i->width+p.x];
}

void add_point_to_list(point_list_t** list, point_t p) {
  point_list_t* old_list = *list;
  *list = malloc(sizeof(point_list_t));
  (*list)->p = p;
  (*list)->n = old_list;
}

double distance(point_t a, point_t b) {
  double dx=(double)a.x-b.x;
  double dy=(double)a.y-b.y;
  double dz=(double)a.z-b.z;
  return sqrt(dx*dx+dy*dy+dz*dz);
}
