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
  int p[3];
  int index;
} point_t;
typedef struct {
  double p[3];
} dpoint_t;

void init_rng(image_t*);
void compute_depth(image_t*);
point_t random_point(const image_t*);
unsigned short pixel_get(const image_t*, point_t);
double distance(point_t,point_t);

point_t* perform_sample(const image_t*, int);
void replace_in_sample(const image_t*, point_t*, int);

typedef struct kd_node {
  point_t location;
  int index;
  struct kd_node* left;
  struct kd_node* right;
} kdtree_t;

kdtree_t* kdtree_build(const point_t*,int);
const kdtree_t* kdtree_search(const kdtree_t*, point_t);

#endif
