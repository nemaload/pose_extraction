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


#define DISTANCE_I(c) double d##c = (double)b.p[c]-a.p[c];
#define DISTANCE_RET return sqrt(d0*d0+d1*d1+d2*d2);
#define DISTANCE_I_ \
  FOREACH3(DISTANCE_I) \
  DISTANCE_RET

inline double distance_id(dpoint_t a, point_t b) {
  DISTANCE_I_
}

inline double distance_i(point_t a, point_t b) {
  DISTANCE_I_
}

#define DISTANCE_D(c) double d##c = b.p[c]-a.p[c];
inline double distance(dpoint_t a, dpoint_t b) {
  FOREACH3(DISTANCE_D)
  DISTANCE_RET
}

#define REFLECT_D(c) b.p[c]=v.p[c]+(v.p[c]-a.p[c]);
dpoint_t reflect(dpoint_t a, dpoint_t v) {
  dpoint_t b;
  FOREACH3(REFLECT_D)
  return b;
}

#define REFLECT2_D(c) b.p[c]=v.p[c]+2*(v.p[c]-a.p[c]);
dpoint_t reflect2(dpoint_t a, dpoint_t v) {
  dpoint_t b;
  FOREACH3(REFLECT2_D)
  return b;
}

#define PRINT_D(c) printf("\t%lf", p.p[c]);
void print_dpoint(dpoint_t p) {
  printf("point:");
  FOREACH3(PRINT_D)
  printf("\n");
}

point_t* perform_sample(const image_t* image, int n, double threshhold) {
  int i=0;
  point_t* result=malloc(sizeof(point_t)*n);
  for(i=0;i<n;i++)
  while(i<n) {
    point_t p = random_point(image);
    if(pixel_get(image,p) > threshhold) {
      result[i++]=p;
      progress(i,n,0,"brights");
    }
  }
  return result;
}

void replace_in_sample(const image_t* image, point_t* sample, int k, int n, double threshhold) {
  int i=0;
  while(i<k) {
    point_t p = random_point(image);
    if(pixel_get(image,p) > threshhold) {
      i++;
      sample[gsl_rng_uniform_int(image->r,n)]=p;
    }
  }
}

#define COMPARISON_FUNCTION(n) \
  int point_compare##n (const void* p1v, const void* p2v) { \
    const dpoint_t* p1 = (const dpoint_t*)p1v; \
    const dpoint_t* p2 = (const dpoint_t*)p2v; \
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

kdtree_t* kdtree_build_(const dpoint_t* pts, int n, int depth, int* tot_i, int tot_n, kdtree_t* parent) {
  if(n==0) return NULL;
  int axis = depth % 3;
  dpoint_t* plist = malloc(sizeof(dpoint_t)*n);
  memcpy(plist,pts,sizeof(dpoint_t)*n);
  qsort(plist,n,sizeof(dpoint_t),point_compare[axis]);
  int median = n/2;
  kdtree_t* node = malloc(sizeof(kdtree_t));
  //progress((*tot_i)++,tot_n*2,1,"nodes");
  node->axis = axis;
  int i;
  if(parent) {
    for(i=0;i<3*2;i++) {
      node->ranges[i]=parent->ranges[i];
    }
  } else {
    for(i=0;i<3;i++) {
      node->ranges[i*2]=-INFINITY;
      node->ranges[i*2+1]=INFINITY;
    }
  }
  node->ranges[node->axis*2]=plist[0].p[axis];
  node->ranges[node->axis*2+1]=plist[n-1].p[axis];
  node->location = plist[median];
  node->up = parent;
  node->left = kdtree_build_(plist,median,depth+1,tot_i,tot_n,node);
  node->right = kdtree_build_(&(plist[median+1]),n-median-1,depth+1,tot_i,tot_n,node);
  return node;
}

kdtree_t* kdtree_build(const dpoint_t* pts, int n) {
  int tot_i=0;
  kdtree_build_(pts,n,0,&tot_i,n,NULL);
}

const kdtree_t* kdtree_search_(const kdtree_t* here, point_t point, const kdtree_t* best, double best_dist, int axis) {
  if(best==NULL)
    best=here;

  double here_dist = distance_id(here->location,point);
  if(here_dist < best_dist) {
    best=here;
    best_dist = here_dist;
  }

  kdtree_t *near_child, *far_child;
  if(point.p[axis] < here->location.p[axis]) {
    near_child=here->left;
    far_child=here->right;
  } else {
    near_child=here->right;
    far_child=here->left;
  }

  if(near_child) {
    best = kdtree_search_(near_child,point,best,best_dist,(axis+1)%3);
    best_dist = distance_id(best->location,point);
  }

  if(far_child) {
    double corner_distance=0;
    int i;
    for(i=0;i<3;i++) {
      if(point.p[i] > far_child->ranges[i*2+1]) {
        corner_distance+=(point.p[i]-far_child->ranges[i*2+1])*(point.p[i]-far_child->ranges[i*2+1]);
      } else if (point.p[i] < far_child->ranges[i*2]) {
        corner_distance+=(point.p[i]-far_child->ranges[i*2])*(point.p[i]-far_child->ranges[i*2]);
      }
    }

    if(sqrt(corner_distance) < best_dist) {
      best = kdtree_search_(far_child,point,best,best_dist,(axis+1)%3);
      best_dist = distance_id(best->location,point);
    }
  }
  return best;
}

const kdtree_t* kdtree_search(const kdtree_t* tree, point_t p) {
  kdtree_search_(tree,p,NULL,0.0,0);
}

void kdtree_free(kdtree_t* tree) {
  if(tree==NULL) return;
  kdtree_free(tree->left);
  kdtree_free(tree->right);
  free(tree);
}
