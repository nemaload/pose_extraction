#define QU_(x) #x
#define QU(x) QU_(x)

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>


#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include "util.h"

#define PROGRESS_BAR_WIDTH  120

/*
 * I wrote a little C macro hackery called ARGBOILER,
 * which uses the argtable2 library to parse arguments.
 * This is the sligtly awkward, but comparatively concise,
 * argument specification list.
 */

#define ARGBOILER(ARG) \
  ARG(ARG_STR1,filename,NULL,NULL,"the input image (as raw data)",NULL) \
  ARG(ARG_INT1,width,"w","width","the width of each image slice",-1) \
  ARG(ARG_INT1,height,"h","height","the height of each image slice",-1) \
  ARG(ARG_INT0,sd_sample_size,NULL,"sdss","the sample size for computing standard deviation",1000) \
  ARG(ARG_INT0,mst_sample_size,NULL,"mstss","the sample size for making the MST",150) \
  ARG(ARG_INT0,refine_sample_size,NULL,"rfss","the sample size of E_image in refining the backbone",800) \
  ARG(ARG_INT0,refine_refresh_size,NULL,"rfrs","the number of E_image samples to replace each iteration",50) \
  ARG(ARG_DBL0,refine_threshhold,NULL,"rfth","the average distance (in pixels) points must move less than to terminate refinement",50.0) \
  ARG(ARG_DBL0,alpha,"a","alpha","weight of E_image; 1 in the original paper",1.0) \
  ARG(ARG_DBL0,beta,"b","beta","weight of E_length; 0.5 in the original paper",0.5) \
  ARG(ARG_DBL0,gamma,"g","gamma","weight of E_smoothness; 0.5 in the original paper",0.5) \
  ARG(ARG_DBL0,delta,"d","delta","'inertia' term (not in the original paper)",0.0)

#include "argboiler.h"

/*
 * These functions are for a primitive sort of profiling.
 */
double elapsed(struct timespec e_t, struct timespec s_t);
typedef enum {START, END} se_t;
void step(se_t se, int inc, const char* desc) {
  static int n = -10;
  static struct timespec s_t_rt;
  static struct timespec s_t_usr;
  static double total_time_rt;
  static double total_time_usr;
  static const char* last_desc;
  if(se==START) {
    n+=inc;
    if(inc==10) n-=n%10;
  }
  char num[5];
  if(n%10)
    snprintf(num,sizeof(num),"%d.%d",n/10,n%10);
  else
    snprintf(num,sizeof(num),"%d",n/10);
  if(se == START) {
    clock_gettime(CLOCK_REALTIME, &s_t_rt);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &s_t_usr);

    if(desc) {
      printf("[Step %s begin (%s)]\n\n",num,desc);
      last_desc = desc;
    } else {
      printf("[Step %s begin]\n\n",num);
      last_desc = NULL;
    }
  } else {
    struct timespec e_t_rt;
    struct timespec e_t_usr;
    double elapsed_rt;
    double elapsed_usr;
    clock_gettime(CLOCK_REALTIME, &e_t_rt);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &e_t_usr);
    if(last_desc && !desc)
      printf("[Step %s end (%s)]\n\n",num,last_desc);
    else if(desc)
      printf("[Step %s end (%s)]\n\n",num,desc);
    else
      printf("[Step %s end]\n\n",num);
    elapsed_rt = elapsed(e_t_rt,s_t_rt);
    elapsed_usr = elapsed(e_t_usr,s_t_usr);

    printf("[Time spent in step %s: %lfms/%lfms]\n",num,elapsed_usr,elapsed_rt);
    total_time_rt += elapsed_rt;
    total_time_usr += elapsed_usr;
    printf("[Running total of time: %lfms/%lfms]\n\n",total_time_usr,total_time_rt);
  }
}
double elapsed(struct timespec e_t, struct timespec s_t) {
    return (e_t.tv_sec-s_t.tv_sec)*1e3+((e_t.tv_nsec-s_t.tv_nsec)*1e-6);
}
void step_start(const char* desc) {
  step(START, 10, desc);
}
void step_end(void) {
  step(END, 0, NULL);
}
void half_step_start(const char* desc) {
  step(START, 5, desc);
}

/*
 * This function is for displaying progress bars.
 */
void progress(int i, int n, int l, char* desc) {
#ifdef PROGRESS_BAR_WIDTH
  static int last_l=0;
  int k;
  if(i<=1) {
    if(i==-1) {
      while(l>last_l) {
        printf("\033[u\n\033[s");
        last_l++;
      }
    }
    printf("\033[1G[>");
    for(k=0;k<PROGRESS_BAR_WIDTH;k++) {
      printf(" ");
    }
    printf("]                   %10s\033[s",desc);
    fflush(stdout);
  } else {
    if(l>last_l) {
      printf("\033[%dE",l-last_l);
      last_l=l;
    } else if (last_l>l) {
      printf("\033[%dF",last_l-l);
      last_l=l;
    }
  }
  if(i>=0) {
    if ((n/PROGRESS_BAR_WIDTH) == 0 || (i-1)%(n/PROGRESS_BAR_WIDTH) == 0) {
      if(i/((double)n/PROGRESS_BAR_WIDTH)<=PROGRESS_BAR_WIDTH){
        printf("\033[%dG=>\033[%dG",2+(int)(int)((i-1)/((double)n/PROGRESS_BAR_WIDTH)),PROGRESS_BAR_WIDTH+4);
        printf("%9d/%9d %10s", i, n, desc);
        fflush(stdout);
      }
    }
    if(i>=n) {
      printf("\033[%dG Done with %s!\033[0K",PROGRESS_BAR_WIDTH+4,desc);
      fflush(stdout);
      if(l==0) {
        printf("\033[u\n\n");
      }
    }
  }
#endif
}

/*
 * Bundle all the little initializations together
 */
void init(image_t* image, const args_t* args) {
  init_rng(image);
  image->width = args->width;
  image->height = args->height;
  compute_depth(image);
}

/*
 * [Step 0]
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
      table[ix++] = distance(list[i],list[j]);
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

point_t* trace_backbone(int tip, const int* mst, const double* distances, const point_t* list, int n, int* backbone_length) {
  int* previous = malloc(sizeof(int)*n);
  int tip2 = find_tip(tip,mst,distances,previous,n);
  point_t* backbone = malloc(sizeof(point_t)*n);
  int i,j=0;
  for(i=tip2;i!=-1;i=previous[i]) {
    backbone[j]=list[i];
    backbone[j].index=j;
    j++;
    progress(j,n,0,"controls");
  }
  *backbone_length = j;
  progress(n,n,0,"controls");
  free(previous);
  return backbone;
}

/*
 * [Step 7]
 * Refine the backbone, iteratively minimizing energy
 * based on length, smoothness, and correspondence to reality
 */

void refine_backbone(const image_t* image, point_t* sample, const args_t* args, point_t* backbone, int n) {
  double iter_delta = INFINITY;
  double iter_delta_init = INFINITY;
  point_t* backbone_new = malloc(n*sizeof(point_t));
  double* total_brightness = malloc(n*sizeof(double));
  dpoint_t* weighted_sum = malloc(n*sizeof(dpoint_t));
  int iterations = 1;
  progress(iterations,100,0,"iterations");
  progress(-1,100,1,"nodes");
  progress(-1,100,2,"pixels");
  progress(-1,100,3,"controls");
  while(iter_delta > args->refine_threshhold) {
    kdtree_t* kdtree;

    memset(total_brightness,0,n*sizeof(double));
    memset(weighted_sum,0,n*sizeof(dpoint_t));

    kdtree = kdtree_build(backbone,n);
    int i;
    if(iter_delta!=INFINITY)
      progress(1,args->refine_sample_size,2,"pixels");
    for(i=0;i<args->refine_sample_size;i++) {
      unsigned short brightness = pixel_get(image,sample[i]);
      int nn = kdtree_search(kdtree, sample[i])->location.index;
      total_brightness[nn]+=brightness;
      weighted_sum[nn].p[0]+=brightness*sample[i].p[0];
      weighted_sum[nn].p[1]+=brightness*sample[i].p[1];
      weighted_sum[nn].p[2]+=brightness*sample[i].p[2];
      if(iter_delta!=INFINITY)
        progress(i+1,args->refine_sample_size,2,"pixels");
    }

    iter_delta = 0;
    for(i=0;i<n;i++) {
      point_t minus2,minus1,here,plus1,plus2;
      if(i<2) {
        minus2=backbone[0];
        minus1=backbone[0];
      } else {
        minus2=backbone[i-2];
        minus1=backbone[i-1];
      }
      if(i>=n-2) {
        plus2=backbone[n-1];
        plus1=backbone[n-1];
      } else {
        plus2=backbone[i+2];
        plus1=backbone[i+1];
      }
      here=backbone[i];
      int c;
      for(c=0;c<3;c++) {
        backbone_new[i].p[c]=(args->alpha*(weighted_sum[i].p[c]/total_brightness[i])\
                             +args->beta*(minus1.p[c]+plus1.p[c])\
                             +args->gamma*(1.5*(minus1.p[c]+plus1.p[c])-0.5*(minus2.p[c]+plus2.p[c]))
                             +args->delta*here.p[c])
                             /(args->alpha+args->beta+args->gamma+args->delta);
      }
      iter_delta+=distance(backbone[i],backbone_new[i]);
      progress(i+1,n,3,"controls");
    }
    iter_delta/=n;
    if(iter_delta_init==INFINITY) iter_delta_init=iter_delta;

    memcpy(backbone,backbone_new,n*sizeof(point_t));

    replace_in_sample(image,sample,args->refine_refresh_size);
    
    //printf("\n%lf,%lf\n",(iter_delta_init*4-iter_delta),(iter_delta_init*4-args->refine_threshhold));
    //progress((int)(iter_delta_init*4-iter_delta),(int)(iter_delta_init*4-args->refine_threshhold),0,"iterations");
    if(iterations<99)
      progress(iterations++,100,0,"iterations");
  }
  progress(100,100,0,"iterations");
  free(total_brightness);
  free(weighted_sum);
}

/*
 * Now, we put it all together!
 */
int main(int argc, char** argv) {
  args_t args;
  image_t image;
  double mean;
  double sd;
  double threshhold;

  printf("Parsing command line...\n");
  parse_args(argc, argv, &args);

  printf("Worm straightener v0.0.1\ncreated by David Dalrymple\n============================\n\nLet's straighten this worm!\n\n");

  step_start("mmapping file");
    image.data=open_mmapped_file(args.filename, &image.length);
    if(image.data==NULL) {
      perror("mmap failed");
      return 1;
    }
    printf("File mmapped successfully...\n");
    init(&image, &args);
  step_end();
  
  step_start("computing mean & s.d.");
    compute_sd(&image, args.sd_sample_size, &mean, &sd);
    printf("The standard deviation of %d randomly chosen points is: %lf\nThe mean is: %lf\n", args.sd_sample_size, sd, mean);
    threshhold = mean + sd;
    printf("The threshhold is: %lf\n", threshhold);
  step_end();

  step_start("sampling for MST");
    point_t* w;
    w = sample_bright_points(&image, threshhold, args.mst_sample_size);
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
    point_t* backbone;
    int n;
    backbone=trace_backbone(tip,mst,distances,w,args.mst_sample_size,&n);
  step_end();

  half_step_start("sampling for E_image");
    point_t* refine_sample;
    refine_sample = perform_sample(&image,args.refine_sample_size);
  step_end();

  step_start("Refining backbone");
    refine_backbone(&image,refine_sample,&args,backbone,n);
  step_end();

  return 0;
}

//TODO: Refactor the compute_sd to use perform_sample
//TODO: Rename image.h to util.h
