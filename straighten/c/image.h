#ifndef IMAGE_UTILS
#define IMAGE_UTILS

#include <gsl/gsl_rng.h>

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

typedef struct {
  int x;
  int y;
  int z;
} point_t;

void init_rng(image_t*);
void compute_depth(image_t*);
point_t random_point(const image_t*);
unsigned short pixel_get(const image_t*, point_t);
double distance(point_t,point_t);

#endif
