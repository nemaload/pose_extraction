#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "util.h"
#include "debug.h"

#define X11

#ifdef X11
#include <g2.h>
#include <g2_X11.h>
#endif

/*
 * I wrote a little C macro hackery called ARGBOILER,
 * which uses the argtable2 library to parse arguments.
 * This is the sligtly awkward, but comparatively concise,
 * argument specification list.
 */

#define ARGBOILER(ARG) \
  ARG(ARG_FIL1,input_filename,NULL,NULL,"the input image (as raw data)",NULL) \
  ARG(ARG_INT1,input_width,"w","width","the width of each input image slice",-1) \
  ARG(ARG_INT1,input_height,"h","height","the height of each input image slice",-1) \
  ARG(ARG_INT0,output_width,NULL,"worm-width","the width of the output image",-1) \
  ARG(ARG_INT0,output_height,NULL,"worm-height","the height of the output image",-1) \
  ARG(ARG_INT0,output_extension,NULL,"worm-extension","how many pixels on either edge of the backbone to include (head and tail)",-1) \
  ARG(ARG_FIL0,output_filename,"o","output","the output image (as raw data)",NULL) \
  ARG(ARG_INT0,sd_sample_size,NULL,"sdss","the sample size for computing standard deviation",1000) \
  ARG(ARG_DBL0,thresh_sds,"t","thresh","the number of standard deviations above mean makes a pixel considered part of the worm",0.5) \
  ARG(ARG_INT0,mst_sample_size,NULL,"mstss","the sample size for making the MST",120) \
  ARG(ARG_INT0,refine_sample_size,NULL,"rfss","the sample size of E_image in refining the backbone",3000) \
  ARG(ARG_INT0,refine_refresh_size,NULL,"rfrs","the number of E_image samples to replace each iteration",10) \
  ARG(ARG_DBL0,refine_threshhold,NULL,"rfth","the average distance (in pixels) points must move less than to terminate refinement",1.53) \
  ARG(ARG_INT0,delta_history,NULL,"rfdh","the number of iterations the points must move very little in a row to count",150) \
  ARG(ARG_LIT0,spread_voronoi,"3","spread","adjust control points using the nearest pixels of neighboring control points as well as their own",0) \
  ARG(ARG_LIT0,use_brightness,"r","weight","weight E_image by pixel brightness instead of just threshholding",0) \
  ARG(ARG_LIT0,no_interpolate,"c","no-interpolate","don't interpolate input pixels when restacking",0) \
  ARG(ARG_DBL0,alpha,"a","alpha","weight of E_image; 1 in the original paper",0.15) \
  ARG(ARG_DBL0,beta,"b","beta","weight of E_length; 0.5 in the original paper",1.1) \
  ARG(ARG_DBL0,gamma,"g","gamma","weight of E_smoothness; 0.5 in the original paper",0.6) \
  ARG(ARG_DBL0,delta,"d","delta","'inertia' term (not in the original paper)",2.2) \
  ARG(ARG_DBL0,image_scale,"s","scale","factor to divide image width and height by in display",4.0) \
  ARG(ARG_LITN,precache,"p","precache","mmap the image and read all the pixels; performing no processing if specified twice.",0) \

#include "argboiler.h"

/*
 * Bundle all the little initializations together
 */
void init(image_t* image, const args_t* args) {
  init_rng(image);
  image->width = args->input_width;
  image->height = args->input_height;
  compute_depth(image);
}

/*
 * [Step 0]
 * Here's where the action begins - step 0 of the algorithm, loading
 * the file - well, not loading it per se, but mmaping it.
 */
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

/*
 * [Step 0.5]
 * Precaching option (-p)
 */
void precache_file(image_t input) {
    int i,j;
    volatile unsigned short foo;
    long pagesize = sysconf(_SC_PAGESIZE);
    for(i=0;i*pagesize<input.length/2;i++) {
      progress(i+1,input.length/2/pagesize,0,"pages");
      for(;((i+1)%10)!=0 && i*pagesize<input.length/2;i++) {
        foo+=((unsigned short*)input.data)[i*pagesize];
      }
    }
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

  //gsl_ran_sample(image->r, sample, sample_size, image->data, image->length/2, 2);
  for(i=0; i<sample_size; i++) {
    sample[i]=pixel_get(image,random_point(image));
    progress(i+1,sample_size,0,"samples");
  }

  *mean = gsl_stats_ushort_mean(sample, 1, sample_size);
  *sd = gsl_stats_ushort_sd_m(sample, 1, sample_size, *mean);
  free(sample);
}

/*
 * [Step 2]
 * Given a threshhold value, randomly sample points, adding them to
 * a point list if they are above the threshhold, until we have a
 * specified number of points.
 */

point_t* sample_bright_points(const image_t* image, double threshhold, int n) {
  point_t* list = malloc(sizeof(point_t)*n);
  point_t p;
  int i=0;
  while(i<n) {
    p = random_point(image);
    if(pixel_get(image,p) > threshhold) {
      list[i++]=p;
      progress(i,n,0,"brights");
    }
  }
  return list;
}

/*
 * [Step 3]
 * Compute the distances between all the points sampled in Step 2.
 */

double* compute_distances(point_t* list, int n) {
  int n2=n*n;
  double* table = malloc(n2*sizeof(double));
  int ix=0;
  int i,j;
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      table[ix++] = distance_i(list[i],list[j]);
      progress(ix,n2,0,"distances");
    }
  }
  return table;
}

/*
 * [Step 4]
 * Compute the minimum spanning tree, using Prim's algorithm
 */

int* compute_mst(double* distances, int n) {
  int n2=n*n;
  int* mst = malloc(n2*sizeof(int));
  memset(mst,-1,n2*sizeof(int));
  int* set = malloc(n*sizeof(int));
  int set_b = 1;
  int i,j;
  for(i=0; i<n; i++) {
    set[i]=i;
  }
  while(set_b<n) {
    progress(set_b,n-1,0,"vertices");
    double min_weight = INFINITY;
    int min_edge_1, min_edge_2;
    for(i=0;i<set_b;i++) {
      for(j=set_b;j<n;j++) {
        double d = distances[set[i]*n+set[j]];
        if(d<min_weight) {
          min_weight = d;
          min_edge_1 = i;
          min_edge_2 = j;
        }
      }
    }
    for(i=0;mst[set[min_edge_1]*n+i]!=-1;i++);
    mst[set[min_edge_1]*n+i]=set[min_edge_2];
    for(i=0;mst[set[min_edge_2]*n+i]!=-1;i++);
    mst[set[min_edge_2]*n+i]=set[min_edge_1];
    int tmp;
    tmp = set[set_b];
    set[set_b] = set[min_edge_2];
    set[min_edge_2] = tmp;
    set_b++;
  }
  free(set);
  return mst;
}

void print_mst(int* mst, int n) {
  int i,j;
  for(i=0;i<n;i++) {
    printf("%4d:",i);
    for(j=0;j<n;j++) {
      int x = mst[i*n+j];
      if(x>=0)
        printf("%4d",x);
    }
    printf("\n");
  }
}

/*
 * [Step 5/(6)]
 * Identify the backbone of the MST by dual BFS
 */

void bfs(int k, int old_k, const int* mst, const double* distances, double a, double* output, int* previous, int n) {
  if(output[k] < 0) {
    int i;
    output[k] = a;
    previous[k] = old_k;
    for(i=0;mst[k*n+i]>=0;i++) {
      bfs(mst[k*n+i],k,mst,distances,a+distances[k*n+mst[k*n+i]],output,previous,n);
    }
  }
}

int find_tip(int k, const int* mst, const double* distances, int* previous, int n) {
  double* output = malloc(sizeof(double)*n);
  int malloced_previous = 0;
  if(previous == NULL) {
    previous = malloc(sizeof(int)*n);
    malloced_previous = 1;
  }
  int i;
  for(i=0;i<n;i++) output[i]=-1;
  bfs(k,-1,mst,distances,0,output,previous,n);
  double max=0;
  int argmax=-1;
  for(i=0;i<n;i++) {
    if(output[i]>max) {
      argmax=i;
      max=output[i];
    }
  }
  if(malloced_previous) free(previous);
  free(output);
  return argmax;
}

/*
 * [Step 6]
 * Trace the backbone into a new list of points,
 * discarding the rest of the MST
 */

dpoint_t* trace_backbone(int tip, const int* mst, const double* distances, const point_t* list, int n, int* backbone_length) {
  int* previous = malloc(sizeof(int)*n);
  int tip2 = find_tip(tip,mst,distances,previous,n);
  dpoint_t* backbone = malloc(sizeof(dpoint_t)*n);
  int i,j=0,c=0;
  for(i=tip2;i!=-1;i=previous[i]) {
    for(c=0;c<3;c++) {
      backbone[j].p[c]=list[i].p[c];
    }
    backbone[j].index=j;
    j++;
    progress(j,n,0,"controls");
  }
  *backbone_length = j;
  printf("]  Done tracing backbone! (%d/%d points)\e[K\n\n",j,n);
  free(previous);
  return backbone;
}

/*
 * [Step 7]
 * Refine the backbone, iteratively minimizing energy
 * based on length, smoothness, and correspondence to reality
 */

void draw_image(int g2, const image_t* image, const args_t* args) {
  int width = image->width/args->image_scale;
  int height = image->height/args->image_scale;
  static int* pens = NULL;
  if(pens==NULL) {
    pens = malloc(width*height*sizeof(int));
    int i,j;
    point_t p;
    p.p[0] = image->depth/2;
    for(i=0;i<width;i++) {
      p.p[2] = (int)(i*args->image_scale);
      for(j=0;j<height;j++) {
        printf("(%d,%d)\n",i,j);
        p.p[1] = (int)(j*args->image_scale);
        double brightness = pixel_get(image,p)/(double)(0x1<<10);
        int pen = g2_ink(g2,brightness,brightness,brightness);
        pens[j*width+i]=pen;
      }
    }
  }
  g2_image(g2,0.0,0.0,width,height,pens);
}

void refine_backbone(const image_t* image, point_t* sample, const args_t* args, dpoint_t* backbone, int n) {
  double iter_delta = INFINITY;
  double iter_delta_init = INFINITY;
  dpoint_t* backbone_new = malloc(n*sizeof(dpoint_t));
  double* total_brightness = malloc(n*sizeof(double));
  dpoint_t* weighted_sum = malloc(n*sizeof(dpoint_t));
  int iterations = 0;
  int delta_history = 0;
  char spinner[] = "-\\|/";
  char* dh_str = malloc(args->delta_history+2);
  dh_str[args->delta_history+1]='\0';
  int i;
#ifdef X11
  int g2=g2_open_X11(image->width/args->image_scale,image->height/args->image_scale);
  g2_set_auto_flush(g2,0);
#define G2_TRANSFORM(x) \
  g2,\
                                        x.p[2]/args->image_scale,\
      (image->height/args->image_scale)-x.p[1]/args->image_scale
#endif
  while(delta_history < args->delta_history) {
    kdtree_t* kdtree;

    memset(total_brightness,0,n*sizeof(double));
    memset(weighted_sum,0,n*sizeof(dpoint_t));

    kdtree = kdtree_build(backbone,n);
    for(i=0;i<args->refine_sample_size;i++) {
      unsigned short brightness = 1;
      int nn = kdtree_search(kdtree, sample[i])->location.index;
      if(args->use_brightness) {
        brightness = pixel_get(image,sample[i]);
      }
#define ADD_WSUM(c) weighted_sum[nn].p[c]+=brightness*sample[i].p[c];
#define ACCOUNT_POINT_ FOREACH3(ADD_WSUM) total_brightness[nn]+=brightness
#ifndef X11
#define ACCOUNT_POINT ACCOUNT_POINT_
#else
#define ACCOUNT_POINT ACCOUNT_POINT_;\
      g2_pen(g2,7);\
      g2_move(G2_TRANSFORM(backbone[nn]));\
      g2_line_to(G2_TRANSFORM(sample[i]))
#endif
      ACCOUNT_POINT;
      if(args->spread_voronoi) {
        nn--;
        if(nn>0) {
          ACCOUNT_POINT;
        }
        nn+=2;
        if(nn<=n-1) {
          ACCOUNT_POINT;
        }
      }
    }
    kdtree_free(kdtree);
#ifdef X11
    g2_flush(g2);
#endif

    iter_delta = 0;
    backbone_new[0]=backbone[0];
    for(i=0;i<n;i++) {
      dpoint_t minus2,minus1,here,plus1,plus2;
      if(i==1) {
        minus1=backbone[0];
        minus2=reflect(backbone[1],backbone[0]);
      } else if (i==0) {
        minus1=reflect(backbone[1],backbone[0]);
        minus2=reflect2(backbone[1],backbone[0]);
      } else {
        minus2=backbone[i-2];
        minus1=backbone[i-1];
      }
      if(i==n-2) {
        plus1=backbone[n-1];
        plus2=reflect(backbone[n-2],backbone[n-1]);
      }else if(i==n-1) {
        plus1=reflect(backbone[n-2],backbone[n-1]);
        plus2=reflect2(backbone[n-2],backbone[n-1]);
      } else {
        plus2=backbone[i+2];
        plus1=backbone[i+1];
      }
      here=backbone[i];

      int c;
      if(total_brightness[i]!=0) {
        for(c=0;c<3;c++) {
          backbone_new[i].p[c]=(args->alpha*(weighted_sum[i].p[c]/total_brightness[i]));
        }
      } else {
        for(c=0;c<3;c++) {
          backbone_new[i].p[c]=(args->alpha*here.p[c]);
        }
      }
      for(c=0;c<3;c++) {
        backbone_new[i].p[c]+=(args->beta*(minus1.p[c]+plus1.p[c])\
                             +args->gamma*(1.5*(minus1.p[c]+plus1.p[c])-0.5*(minus2.p[c]+plus2.p[c]))\
                             +args->delta*here.p[c]);
        backbone_new[i].p[c]/=(args->alpha+2*args->beta+2*args->gamma+args->delta);
      }
      iter_delta+=distance(backbone[i],backbone_new[i]);
      backbone_new[i].index=i;
    }
    iter_delta/=n;
    if(iter_delta_init==INFINITY) iter_delta_init=iter_delta;
    if(iter_delta < args->refine_threshhold) {
      delta_history++;
    } else {
      delta_history = 0;
    }
    //printf("iter_delta: %lf\n", iter_delta);

    memcpy(backbone,backbone_new,n*sizeof(dpoint_t));

#ifdef X11
#define G2_DRAW_BACKBONE \
    g2_pen(g2,0);\
    g2_filled_rectangle(g2,0,0,(double)image->width/args->image_scale,(double)image->height/args->image_scale);\
    g2_move(G2_TRANSFORM(backbone[0]));\
    g2_pen(g2,19);\
    for(i=1;i<n;i++) {\
      g2_line_to(G2_TRANSFORM(backbone[i]));\
    }
    G2_DRAW_BACKBONE
#endif

    replace_in_sample(image,sample,args->refine_refresh_size,args->refine_sample_size,image->threshhold);

    memset(dh_str,32,args->delta_history+1);
    memset(dh_str,'=',delta_history);
    dh_str[delta_history]='>';
    printf("\e[0G[\e[1;32;44m%c\e[m] iterations: %d; delta: %s%lf\t\e[32;44m%s\e[m]",\
        spinner[iterations%(sizeof(spinner)-1)],\
        iterations++,\
        (iter_delta>args->refine_threshhold*1.2)?\
          "\e[31m":\
          (iter_delta<args->refine_threshhold)?\
            "\e[32m":\
            "\e[33m",\
        iter_delta,\
        dh_str);
    fflush(stdout);
  }
  printf("\n\n");
#ifdef X11
  G2_DRAW_BACKBONE
  g2_flush(g2);
  usleep(100000);
#endif
  free(total_brightness);
  free(weighted_sum);
}

/*
 * [Step 9]
 * Open the output image with mmap
 */
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

/*
 * [Step 8]
 * Actually straighten the image by restacking along a spline defined by the backbone.
 */
void restack_image(image_t* dst, const image_t* src, const args_t* args, dpoint_t* backbone, int n) {
  /*
   * First, we define the multidimensional spline, by walking along the control
   * points, measuring the total distance walked, and using that as the parameter
   * for three single-dimensional splines.
   */
  double* xa = malloc(sizeof(double)*n);
#define PTS_LIST_ALLOC(c) \
  double* y##c = malloc(sizeof(double)*n);
  FOREACH3(PTS_LIST_ALLOC)
  int i,j,k;
  xa[0]=0.0;
  for(i=1;i<n;i++) {
    xa[i] = xa[i-1]+distance(backbone[i],backbone[i-1]);
#define PTS_SPLIT(c) \
    y##c[i] = backbone[i].p[c];
    FOREACH3(PTS_SPLIT)
  }
  i=0;
  FOREACH3(PTS_SPLIT)
#define GSL_INIT(c) \
  gsl_spline* spl##c = gsl_spline_alloc(gsl_interp_cspline,n); \
  gsl_interp_accel* accel##c = gsl_interp_accel_alloc(); \
  gsl_spline_init(spl##c,xa,y##c,n);
  FOREACH3(GSL_INIT)

  const char* filename;
  const char* suffix=".out";
  if(args->output_filename==NULL) {
    char * filename_=calloc(strlen(args->input_filename)+strlen(suffix),1);
    strcat(filename_,args->input_filename);
    strcat(filename_,suffix);
    filename = filename_;
  } else {
    filename = args->output_filename;
  }
  if(args->output_width==-1) {
    dst->width=src->width/4;
    printf("Output width automatically determined: %d\n",dst->width);
  } else {
    dst->width = args->output_width;
  }
  int extension;
  if(args->output_extension==-1) {
    extension = dst->width;
  } else {
    extension = args->output_extension;
  }
  if(args->output_height==-1) {
    dst->height=(int)xa[n-1]+2*extension;
    printf("Output height automatically determined: %d\n",dst->height);
  } else {
    dst->height = args->output_height;
  }
  if(args->output_width==-1 || args->output_height==-1) {
    printf("\n");
  }
  int length;

#ifndef RESTACK_3D //TODO: do 3D restacking

  dst->depth = 1;
  length = dst->width*dst->height*2;
  dst->data = (open_mmapped_file_write(filename,length));
  unsigned short* data = (unsigned short*)dst->data;

  for(i=0;i<dst->height;i++) {
#define GET_P_D(c) \
    double p##c = gsl_spline_eval(spl##c,i-extension,accel##c); \
    double d##c = gsl_spline_eval_deriv(spl##c,i-extension,accel##c); \
    double dy##c;
    FOREACH3(GET_P_D)
    dy0=0;
    dy1=d2;
    dy2=-d1;
    double magnitude = sqrt(dy1*dy1+dy2*dy2);
    dy1/=magnitude;
    dy2/=magnitude;
    /*if(i<50) {
      printf("spline point:\t%lf\t%lf\t%lf\nd:\t%lf\t%lf\t%lf\ndy:\t%lf\t%lf\t%lf\n",p0,p1,p2,d0,d1,d2,dy0,dy1,dy2);
    }*/
#define INIT_P(c) \
    p##c -= dy##c*dst->width/2;
    FOREACH3(INIT_P)
    for(j=0;j<dst->width;j++) {
      unsigned short pixel = 0;
      if(p0>0&&p1>0&&p2>0&&p0<src->depth-1&&p1<src->height-1&&p2<src->width-1) {
        if(args->no_interpolate) {
          pixel = ((unsigned short*)src->data)[(int)(p0)*src->height*src->width+(int)(p1)*src->width+(int)(p2)];
        } else {
#define TRUNC_P(c) \
          double t##c = p##c - (int)p##c;
          FOREACH3(TRUNC_P)
          double pixel_d=0;
#define ONE_MINUS(a,x) (1-a+(2*a-1)*x)
#define GET_CORNER(z,y,x) \
          pixel_d += (ONE_MINUS(z,t0)*ONE_MINUS(y,t1)*ONE_MINUS(x,t2))*((unsigned short*)src->data)[((int)(p0)+z)*src->height*src->width+((int)(p1)+y)*src->width+((int)(p2)+z)]
          GET_CORNER(0,0,0);
          GET_CORNER(0,0,1);
          GET_CORNER(0,1,0);
          GET_CORNER(0,1,1);
          GET_CORNER(1,0,0);
          GET_CORNER(1,0,1);
          GET_CORNER(1,1,0);
          GET_CORNER(1,1,1);

          pixel = (unsigned short) pixel_d;
        }
      }
      data[i*dst->width+j]=pixel;
#define INC_P(c) \
      p##c += dy##c;
      FOREACH3(INC_P)
    }
    progress(i+1,dst->height,0,"planes");
  }
#endif
}

/*
 * Now, we put it all together!
 */
int main(int argc, char** argv) {
  args_t args;
  image_t input;
  double mean;
  double sd;

  printf("Parsing command line...\n");
  parse_args(argc, argv, &args);

  printf("Worm straightener v0.0.1\ncreated by David Dalrymple\n============================\n\n");

  step_start("mmapping file");
    input.data=open_mmapped_file_read(args.input_filename, &input.length);
    if(input.data==NULL) {
      perror("mmap failed");
      return 1;
    }
    printf("File mmapped successfully...\n");
    init(&input, &args);
  step_end();

  if(args.precache > 0) {
    half_step_start("precaching file");
      precache_file(input);
    step_end();
    if(args.precache > 1) {
      return 0;
    }
  }
  
  step_start("computing mean & s.d.");
    compute_sd(&input, args.sd_sample_size, &mean, &sd);
    printf("The standard deviation of %d randomly chosen points is: %lf\nThe mean is: %lf\n", args.sd_sample_size, sd, mean);
    input.threshhold = mean + sd*args.thresh_sds;
    printf("The threshhold is: %lf\n", input.threshhold);
  step_end();

  step_start("sampling for MST");
    point_t* w;
    w = sample_bright_points(&input, input.threshhold, args.mst_sample_size);
  step_end();

  step_start("computing distances for MST");
    double* distances;
    distances = compute_distances(w, args.mst_sample_size);
  step_end();

  step_start("Prim's algorithm");
    int* mst;
    mst = compute_mst(distances, args.mst_sample_size);
  step_end();

  //print_mst(mst,args.mst_sample_size);

  step_start("Finding tip");
    int tip;
    tip = find_tip(0,mst,distances,NULL,args.mst_sample_size);
  step_end();

  step_start("Tracing backbone");
    dpoint_t* backbone;
    int n;
    backbone=trace_backbone(tip,mst,distances,w,args.mst_sample_size,&n);
  step_end();

  half_step_start("sampling for E_image");
    point_t* refine_sample;
    refine_sample = perform_sample(&input,args.refine_sample_size,input.threshhold);
  step_end();

  step_start("Refining backbone");
    refine_backbone(&input,refine_sample,&args,backbone,n);
  step_end();

  step_start("Restacking to output file");
    image_t output;
    restack_image(&output,&input,&args,backbone,n);
  step_end();

  return 0;
}
