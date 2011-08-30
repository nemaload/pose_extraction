#ifndef IMAGE_UTILS
#define IMAGE_UTILS

#include <gsl/gsl_rng.h>

#define FOREACH3(M) \
  M(2) \
  M(1) \
  M(0)

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
  int index;
} dpoint_t;

void init_rng(image_t*);
void compute_depth(image_t*);
point_t random_point(const image_t*);
unsigned short pixel_get_(const unsigned short*,int,int,int,int,int);
unsigned short pixel_get(const image_t*, point_t);
double distance_i(point_t,point_t);
double distance(dpoint_t,dpoint_t);

dpoint_t reflect(dpoint_t,dpoint_t);
dpoint_t reflect2(dpoint_t,dpoint_t);

void print_dpoint(dpoint_t);

point_t* perform_sample(const image_t*, int);
void replace_in_sample(const image_t*, point_t*, int, int);

typedef struct kd_node {
  dpoint_t location;
  double ranges[6];
  int axis;
  struct kd_node* up;
  struct kd_node* left;
  struct kd_node* right;
} kdtree_t;

kdtree_t* kdtree_build(const dpoint_t*,int);
const kdtree_t* kdtree_search(const kdtree_t*, point_t);
void kdtree_free(kdtree_t*);

void* open_mmapped_file_read(const char*, int*);
void precache_file(image_t);
void* open_mmapped_file_write(const char*, int);

struct pqn {
  int x;
  int y;
  struct pqn* next;
};
struct pq {
  struct pqn* head;
  struct pqn* tail;
};
void enpq(struct pq*, int, int);
int depq(struct pq*, int*, int*);
void empq(struct pq*);

//When you "test" make sure you haven't already "act"ed!
void floodfill(int, int, int(*test)(int, int), void (*act)(int, int));

#endif
