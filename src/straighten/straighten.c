#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <fenv.h>
#include <math.h>
#include <pthread.h>
#include <sched.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "../common/util.h"
#include "../common/debug.h"

#ifndef NO_X11
#define X11
#endif

#ifdef X11
#include <g2.h>
#include <g2_X11.h>
#endif

/*
 * I wrote a little C macro hackery called ARGBOILER,
 * which uses the argtable2 library to parse arguments.
 */

static const int automatic=-1;

#define ARGBOILER(ARG) \
  ARG(ARG_FIL1,input_filename,NULL,NULL,"the input image (as raw data)","") \
  ARG(ARG_INT1,input_width,"w","width","the width of each input image slice",-1) \
  ARG(ARG_INT1,input_height,"h","height","the height of each input image slice",-1) \
  ARG(ARG_LITN,precache,"p","precache","mmap the image and read all the pixels; performing no processing if specified twice.",0) \
  ARG(ARG_INT0,sd_sample_size,NULL,"sdss","the sample size for computing standard deviation",5000) \
  ARG(ARG_DBL0,thresh_sds,"t","thresh","the number of standard deviations above mean makes a pixel considered part of the worm",0.5) \
  ARG(ARG_INT0,mst_sample_size,NULL,"mstss","the sample size for making the MST",120) \
  ARG(ARG_INT0,refine_sample_size,NULL,"rfss","the sample size of E_image in refining the backbone",3000) \
  ARG(ARG_INT0,refine_refresh_size,NULL,"rfrs","the number of E_image samples to replace each iteration",100) \
  ARG(ARG_DBL0,refine_threshhold,NULL,"rfth","the average distance (in pixels) points must move less than to terminate refinement",1.55) \
  ARG(ARG_INT0,restart_iterations,NULL,"rfri","if the refinement doesn't converge after this many iterations, restart with a new MST",1000) \
  ARG(ARG_INT0,delta_history,NULL,"rfdh","the number of iterations the points must move very little in a row to count",150) \
  ARG(ARG_LIT0,spread_voronoi,"3","spread","adjust control points using the nearest pixels of neighboring control points as well as their own",0) \
  ARG(ARG_LIT0,use_brightness,"r","weight","weight E_image by pixel brightness instead of just threshholding",0) \
  ARG(ARG_DBL0,alpha,"a","alpha","weight of E_image; 1 in the original paper",0.15) \
  ARG(ARG_DBL0,beta,"b","beta","weight of E_length; 0.5 in the original paper",1.1) \
  ARG(ARG_DBL0,gamma,"g","gamma","weight of E_smoothness; 0.5 in the original paper",0.6) \
  ARG(ARG_DBL0,delta,"d","delta","'inertia' term (not in the original paper)",2.2) \
  ARG(ARG_DBL0,image_scale,"s","scale","factor to divide image width and height by in display",4.0) \
  ARG(ARG_INT0,x_delay_ms,"x","xdelay","milliseconds to wait at the end of step 7 before closing the window",0) \
  ARG(ARG_LIT0,no_interpolate,"c","no-interp","don't interpolate input pixels when restacking",0) \
  ARG(ARG_FIL0,dump_spline,NULL,"dump-spline","dump control points to specified TSV file before restack","") \
  ARG(ARG_FIL0,use_spline,NULL,"use-spline","use spline from specified TSV file instead of determining it","") \
  ARG(ARG_LIT0,no_output,NULL,"no-output","do not restack - used in conjunction with --dump-spline",0) \
  ARG(ARG_INT0,n_threads,"n","threads","number of threads to use (for restacking)",2) \
  ARG(ARG_FIL0,output_filename,"o","output","the output image (as raw data) [default: input+'.out']","") \
  ARG(ARG_INT0,output_width,NULL,"out-width","the width of the output image",automatic) \
  ARG(ARG_INT0,output_height,NULL,"out-height","the height of the output image",automatic) \
  ARG(ARG_INT0,output_extension,NULL,"out-extra","how many pixels on either edge of the backbone to include (head and tail)",automatic) \
  ARG(ARG_INT0,output_slice,NULL,"out-slice","output a 2D slice at the specified depth",automatic) \
  ARG(ARG_LIT0,use_rejsampling,NULL,"rejection-sampling","use rejection sampling instead of brightness-based static threshold",0) \

#include "../common/argboiler.h"

/*
 * Bundle all the little initializations together
 */
void init(image_t* image, const args_t* args) {
  init_rng(image);
  if(args->no_interpolate) {
    fesetround(FE_TONEAREST);
  } else {
    fesetround(FE_DOWNWARD);
  }
  image->width = args->input_width;
  image->height = args->input_height;
  compute_depth(image);
  printf("l %d w %d h %d %d\n", image->length, image->width, image->height, image->depth);
}

#include "initial_guess.c"
#include "refine_backbone.c"
#include "restack.c"

/*
 * [Dump/Load Spline]
 */
void dump_spline(const char* filename, dpoint_t* backbone, int n) {
  FILE* fd = fopen(filename,"w");
  int i;
  for(i=0;i<n;i++) {
    fprintf(fd,"%.10lf\t%.10lf\t%.10lf\n",backbone[i].p[0],backbone[i].p[1],backbone[i].p[2]);
  }
  fclose(fd);
}

#define LOAD_SPLINE_BUF 50
dpoint_t* load_spline(const char* filename, int* n) {
  *n=0;
  dpoint_t* backbone = malloc(sizeof(dpoint_t)*LOAD_SPLINE_BUF);
  FILE* fd = fopen(filename,"r");
  while(fscanf(fd,"%lf\t%lf\t%lf\n",&(backbone[*n].p[0]),&(backbone[*n].p[1]),&(backbone[*n].p[2]))!=EOF) {
    (*n)++;
    if((*n)%LOAD_SPLINE_BUF==0) {
      dpoint_t* tmp = realloc(backbone,sizeof(dpoint_t)*((*n)+LOAD_SPLINE_BUF));
      if(!tmp) {
        perror("Could not realloc");
      } else {
        backbone = tmp;
      }
    }
  }
  return backbone;
}

/*
 * Now, we put it all together!
 */
int main(int argc, char** argv) {
  args_t args;
  image_t input;
  point_t* w;
  double* distances;
  int* mst;
  int tip;
  dpoint_t* backbone;
  point_t* refine_sample;
  int n;

  printf("Parsing command line...\n");
  parse_args(argc, argv, &args);

  printf("\e[AWorm straightener v0.0.1\ncreated by David Dalrymple\n============================\n\n");

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
  
  if(!args.use_spline[0]) {
    step_start("computing mean & s.d.");
      if (args.use_rejsampling) {
        printf("rejection sample instead!\n");
        input.threshold = NAN;
      } else {
        double mean, sd;
        compute_sd(&input, args.mst_sample_size, &mean, &sd);
        input.threshold = mean + sd * args.thresh_sds;
        printf("mean %f sd %f threshold %f\n", mean, sd, input.threshold);
      }
    step_end();

    int refine_success=0;
    while(!refine_success) {
      step_start("sampling for MST");
        w = sample_bright_points(&input, args.mst_sample_size);
      step_end();

      step_start("computing distances for MST");
        distances = compute_distances(w, args.mst_sample_size);
      step_end();

      step_start("Prim's algorithm");
        mst = compute_mst(distances, args.mst_sample_size);
      step_end();

      //print_mst(mst,args.mst_sample_size);

      step_start("Finding tip");
        tip = find_tip(0,mst,distances,NULL,args.mst_sample_size);
      step_end();

      step_start("Tracing backbone");
        backbone=trace_backbone(tip,mst,distances,w,args.mst_sample_size,&n);
      step_end();

      half_step_start("sampling for E_image");
        refine_sample = perform_sample(&input, args.refine_sample_size);
      step_end();

      step_start("Refining backbone");
        refine_success=refine_backbone(&input,refine_sample,&args,backbone,n);
      step_end();
    }
  } else {
    step_start("Reading spline");
      backbone=load_spline(args.use_spline,&n);
    step_end();
  }

  if(args.dump_spline[0]) {
    step_start("Dumping spline");
      dump_spline(args.dump_spline,backbone,n);
    step_end();
  }

  if(!args.no_output) {
    step_start("Restacking to output file");
      image_t output;
      restack_image(&output,&input,&args,backbone,n);
    step_end();
  }

  return 0;
}
