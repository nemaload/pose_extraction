#include "util.h"
#include "debug.h"
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <png.h>
#define RANDOM_SEED

int _DebugColors = 0;

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

inline unsigned short pixel_get_(const unsigned short* data, int p0, int p1, int p2, int width, int height) {
  return data[p0*width*height+p1*width+p2];
}

inline unsigned short pixel_get(const image_t* i, point_t p) {
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

point_t* perform_sample(const image_t* image, int n) {
  int i=0;
  point_t* result=malloc(sizeof(point_t)*n);
  for(i=0;i<n;i++)
  while(i<n) {
    point_t p = random_point(image);
    if(pixel_get(image,p) > (isnan(image->threshold) ? gsl_rng_uniform_int(image->r, 1<<16) : image->threshold)) {
      result[i++]=p;
      progress(i,n,"brights");
    }
  }
  return result;
}

void replace_in_sample(const image_t* image, point_t* sample, int k, int n) {
  int i=0;
  while(i<k) {
    point_t p = random_point(image);
    if(pixel_get(image,p) > (isnan(image->threshold) ? gsl_rng_uniform_int(image->r, 1<<16) : image->threshold)) {
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

void* open_mmapped_file_read(const char* filename, int* length) {
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

void precache_file(image_t input) {
    int i,j;
    volatile unsigned short foo;
    long pagesize = sysconf(_SC_PAGESIZE);
    for(i=0;i*pagesize<input.length/2;i++) {
      progress(i+1,input.length/2/pagesize,"pages");
      for(;((i+1)%10)!=0 && i*pagesize<input.length/2;i++) {
        foo+=((unsigned short*)input.data)[i*pagesize];
      }
    }
}

void* open_mmapped_file_write(const char* filename, int length) {
  struct stat fs;
  int fd;
  void* region;
  
  //Now get a file descriptor and mmap!
  fd = open(filename, O_RDWR|O_TRUNC|O_CREAT, S_IRUSR|S_IWUSR|S_IROTH|S_IWOTH);
  if(fd<0) {
    perror("couldn't open file");
    printf("file was: %s\n",filename);
  }
  ftruncate(fd,length);
  region=mmap(NULL, length, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
  if(region==MAP_FAILED) {
    perror("couldn't mmap file");
  }

  return region;
}

void copy_file(const char* dest, const char* src) {
  int n;
  void* s = open_mmapped_file_read(src,&n);
  void* d = open_mmapped_file_write(dest,n);
  memcpy(d,s,n);
  munmap(s,n);
  munmap(d,n);
}

/* Convert from L*a*b* doubles to XYZ doubles */
void lab2xyz(double* x, double* y, double* z, double l, double a, double b) {
  double finv(double t) {
    return (t>(6.0/29.0))?(t*t*t):(3*(6.0/29.0)*(6.0/29.0)*(t-4.0/29.0));
  }
  double sl = (l+0.16)/1.16;
  double ill[3] = {0.9643,1.00,0.8251}; //D50
  *y = ill[1] * finv(sl);
  *x = ill[0] * finv(sl + (a/5.0));
  *z = ill[2] * finv(sl - (b/2.0));
  if(_DebugColors) printf("lab2xyz: [ %lf , %lf , %lf ] -> [ %lf , %lf, %lf ]\n",l,a,b,*x,*y,*z);
}

/* Convert from XYZ doubles to sRGB bytes */
void xyz2rgb(unsigned char* r, unsigned char* g, unsigned char* b, double x, double y, double z) {
  double rl =  3.2406*x - 1.5372*y - 0.4986*z;
  double gl = -0.9689*x + 1.8758*y + 0.0415*z;
  double bl =  0.0557*x - 0.2040*y + 1.0570*z;
  int clip = (rl < 0.001 || rl > 0.999 || gl < 0.001 || gl > 0.999 || bl < 0.001 || bl > 0.999);
  if(clip) {
    rl = (rl<0.001)?0.0:((rl>0.999)?1.0:rl);
    gl = (gl<0.001)?0.0:((gl>0.999)?1.0:gl);
    bl = (bl<0.001)?0.0:((bl>0.999)?1.0:bl);
  }
  //if(clip) {rl=1.0;gl=bl=0.0;}
  double correct(double cl) {
    double a = 0.055;
    return (cl<=0.0031308)?(12.92*cl):((1+a)*pow(cl,1/2.4)-a);
  }
  *r = (unsigned char)(255.0*correct(rl));
  *g = (unsigned char)(255.0*correct(gl));
  *b = (unsigned char)(255.0*correct(bl));
  if(_DebugColors) printf("xyz2rgb: [ %lf , %lf , %lf ] -> [ %hhu , %hhu, %hhu ]\n",x,y,z,*r,*g,*b);
}

/* Convert from LAB doubles to sRGB bytes */
void lab2rgb(unsigned char* R, unsigned char* G, unsigned char* B, double l, double a, double b) {
  double x,y,z;
  lab2xyz(&x,&y,&z,l,a,b);
  xyz2rgb(R,G,B,x,y,z);
}

void lab2pix(void* rgb, double l, double a, double b) {
  unsigned char* ptr = (unsigned char*)rgb;
  lab2rgb(ptr,ptr+1,ptr+2,l,a,b);
}

void xyz2pix(void* rgb, double x, double y, double z) {
  unsigned char* ptr = (unsigned char*)rgb;
  xyz2rgb(ptr,ptr+1,ptr+2,x,y,z);
}

void lrl(double* L, double* r, double l) {
  *L = l*0.7;   //L of L*a*b*
  *r = l*0.301+0.125; //chroma
}

/* Convert from a qualitative parameter l and a quantitative parameter c to a 24-bit pixel */
void cl2pix(void* rgb, double c, double l) {
  unsigned char* ptr = (unsigned char*)rgb;
  double L,r;
  lrl(&L,&r,l);
  double angle = TAU/6.0-c*TAU;
  double a = sin(angle)*r;
  double b = cos(angle)*r;
  lab2rgb(ptr,ptr+1,ptr+2,L,a,b);
}

void csl2lab(double* L, double* a, double* b, double c, double s, double l) {
  double r;
  *L=l*0.7;
  r=0.426*s;
  double angle = TAU/6.0-c*TAU;
  *a = sin(angle)*r;
  *b = cos(angle)*r;
  if(_DebugColors) printf("csl2pix: [ %lf , %lf , %lf ] -> L = %lf , r = %lf , angle = %lf , cos(angle)=%lf , a = %lf , b = %lf\n",c,s,l,L,r,angle,cos(angle),a,b);
}

void csl2xyz(double *x, double *y, double *z, double c, double s, double l) {
  double L, a, b;
  csl2lab(&L,&a,&b,c,s,l);
  lab2xyz(x,y,z,L,a,b);
}

void csl2pix(void* rgb, double c, double s, double l) {
  double L, a, b;
  csl2lab(&L,&a,&b,c,s,l);
  lab2pix(rgb,L,a,b);
}

void hsv2pix(void* rgb, double h, double s, double v) {
  double c = v*s;
  double r,g,b;
  h*=6;
  if(h<1) {
    r = c; g = c*h; b = 0;
  } else if(h<2) {
    r = c*(2-h); g = c; b = 0;
  } else if(h<3) {
    r = 0; g = c; b = c*(h-2);
  } else if(h<4) {
    r = 0; g = c*(4-h); b = c;
  } else if(h<5) {
    r = c*(h-4); g = 0; b = c;
  } else {
    r = c; g = 0; b = c*(6-h);
  }
  double m = v-c;
  r+=m; g+=m; b+=m;
  unsigned char* pix = (unsigned char*)rgb;
  *(pix+0) = (unsigned char)(255.0*r);
  *(pix+1) = (unsigned char)(255.0*g);
  *(pix+2) = (unsigned char)(255.0*b);
}

void export_png(char* filename, int width, int height, int bpc, void* data) {
  FILE* fp = fopen(filename,"wb");
  png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
  png_infop info = png_create_info_struct(png);
  png_init_io(png,fp);
  png_byte color_type, bit_depth;
  int bpp;
  switch(bpc) {
    case 9:  bpp=8;  /* Gray:   8+1 */ bit_depth=8; color_type=PNG_COLOR_TYPE_GRAY; break;
    case 17: bpp=16; /* Gray:  16+1 */ bit_depth=16; color_type=PNG_COLOR_TYPE_GRAY; break;
    case 10: bpp=16; /* GA:     8+2 */ bit_depth=8; color_type=PNG_COLOR_TYPE_GRAY_ALPHA; break;
    case 11: bpp=24; /* Color:  8+3 */ bit_depth=8; color_type=PNG_COLOR_TYPE_RGB; break;
    case 12: bpp=32; /* RGBA:   8+4 */ bit_depth=8; color_type=PNG_COLOR_TYPE_RGB_ALPHA; break;
    case 18: bpp=32; /* GA:    16+2 */ bit_depth=16; color_type=PNG_COLOR_TYPE_GRAY_ALPHA; break;
    case 19: bpp=48; /* RGB:   16+3 */ bit_depth=16; color_type=PNG_COLOR_TYPE_RGB; break;
    case 20: bpp=64; /* RGBA:  16+4 */ bit_depth=16; color_type=PNG_COLOR_TYPE_RGB_ALPHA; break;
  }
  png_set_IHDR(png, info, width, height, bit_depth, color_type, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  bpp/=8;
  png_bytep* rows = malloc(height*sizeof(png_bytep));
  int i;
  for(i=0;i<height;i++) rows[i]=data+bpp*width*i;
  png_set_rows(png, info, rows);
  png_write_png(png,info,PNG_TRANSFORM_SWAP_ENDIAN,NULL);
  free(rows);
  fclose(fp);
}

void enpq(struct pq* q, int x, int y) {
  struct pqn *nqn = malloc(sizeof(struct pqn));
  nqn->x=x; nqn->y=y;
  nqn->next=NULL;
  if(q->tail==NULL) {
    q->head=q->tail=nqn;
  } else {
    q->tail->next=nqn;
    q->tail=nqn;
  }
}

int depq(struct pq* q, int* x, int* y) {
  if(q->head == NULL) {
    return 0;
  }
  *x = q->head->x; *y=q->head->y;
  struct pqn* nhead = q->head->next;
  q->head=nhead;
  if(q->head == NULL) {
    q->tail=NULL;
  }
  return 1;
}

void empq(struct pq* q) {
  q->head=q->tail=NULL;
}

void floodfill(int ix, int iy, int (*test)(int x, int y), void (*act)(int x, int y)) {
  struct pq q[0];
  empq(q);
  enpq(q,ix,iy);
  int x,y;
  while(depq(q,&x,&y)) {
    if(test(x-1,y)) { act(x-1,y); enpq(q,x-1,y); }
    if(test(x+1,y)) { act(x+1,y); enpq(q,x+1,y); }
    if(test(x,y-1)) { act(x,y-1); enpq(q,x,y-1); }
    if(test(x,y+1)) { act(x,y+1); enpq(q,x,y+1); }
  }
}
