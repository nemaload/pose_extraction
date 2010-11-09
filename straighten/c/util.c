#include "util.h"
#include "debug.h"
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#define _GNU_SOURCE
#include <stdlib.h>
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
  p.p[2] = gsl_rng_uniform_int(i->r,i->width);
  p.p[1] = gsl_rng_uniform_int(i->r,i->height);
  p.p[0] = gsl_rng_uniform_int(i->r,i->depth);
  return p;
}

unsigned short pixel_get(const image_t* i, point_t p) {
  return ((unsigned short*)i->data)[p.p[0]*i->width*i->height+p.p[1]*i->width+p.p[2]];
}

double distance(point_t a, point_t b) {
  double dx=(double)a.p[2]-b.p[2];
  double dy=(double)a.p[1]-b.p[1];
  double dz=(double)a.p[0]-b.p[0];
  return sqrt(dx*dx+dy*dy+dz*dz);
}

point_t* perform_sample(const image_t* image, int n) {
  int i;
  point_t* result=malloc(sizeof(point_t)*n);
  for(i=0;i<n;i++) 
    result[i]=random_point(image);
  return result;
}

void replace_in_sample(const image_t* image, point_t* sample, int n) {
  int i;
  for(i=0;i<n;i++)
    sample[i]=random_point(image);
}

#define COMPARISON_FUNCTION(n) \
  int point_compare##n (const void* p1v, const void* p2v) { \
    const point_t* p1 = (const point_t*)p1v; \
    const point_t* p2 = (const point_t*)p2v; \
    if(p1->p[n] > p2->p[n]) \
      return 1; \
    else if (p1->p[n] < p2->p[n]) \
      return -1; \
    else \
      return 0;}

COMPARISON_FUNCTION(0)
COMPARISON_FUNCTION(1)
COMPARISON_FUNCTION(2)

int (*point_compare[3])(const void*,const void*) \
      = {point_compare0,point_compare1,point_compare2};

kdtree_t* kdtree_build_(const point_t* pts, int n, int depth, int* tot_i, int tot_n) {
  if(n==0) return NULL;
  int axis = depth % 3;
  point_t* plist = malloc(sizeof(point_t)*n);
  memcpy(plist,pts,sizeof(point_t)*n);
  qsort(plist,n,sizeof(point_t),point_compare[axis]);
  int median = n/2;
  kdtree_t* node = malloc(sizeof(kdtree_t));
  progress((*tot_i)++,tot_n*2,1,"nodes");
  node->location = plist[median];
  node->left = kdtree_build_(plist,median,depth+1,tot_i,tot_n);
  node->right = kdtree_build_(&(plist[median+1]),n-median-1,depth+1,tot_i,tot_n);
  return node;
}

kdtree_t* kdtree_build(const point_t* pts, int n) {
  int tot_i=0;
  kdtree_build_(pts,n,0,&tot_i,n);
}

const kdtree_t* kdtree_search_(const kdtree_t* here, point_t point, const kdtree_t* best, int axis) {
  if(here==NULL)
    return best;
  if(best==NULL)
    best=here;
  if(distance(here->location,point) < distance(best->location,point))
    best=here;

  int left_nearer;
  kdtree_t* child;
  if(point.p[axis] < here->location.p[axis]) {
    child=here->left;
    left_nearer=1;
  } else {
    child=here->right;
    left_nearer=0;
  }
  best = kdtree_search_(child,point,best,(axis+1)%3);

  if(abs(here->location.p[axis]-point.p[axis]) < distance(best->location,point)) {
    if(left_nearer) {
      child=here->right;
    } else {
      child=here->left;
    }
    best = kdtree_search_(child,point,best,(axis+1)%3);
  }
  return best;
}

const kdtree_t* kdtree_search(const kdtree_t* tree, point_t p) {
  kdtree_search_(tree,p,NULL,0);
}
