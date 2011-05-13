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
  ARG(ARG_DBL0,image_scale,"s","scale","factor to divide image width and height by in display",4.0) \
  ARG(ARG_INT0,x_delay_ms,"x","xdelay","milliseconds to wait at the end of step 7 before closing the window",0) \
  ARG(ARG_INT0,n_threads,"n","threads","number of threads to use (for restacking)",2) \
  ARG(ARG_FIL0,output_filename,"o","output","the output image (as raw data) [default: input+'.out']","") \
  ARG(ARG_INT0,output_width,NULL,"out-width","the width of the output image",automatic) \
  ARG(ARG_INT0,output_height,NULL,"out-height","the height of the output image",automatic) \
  ARG(ARG_INT0,output_slice,NULL,"out-slice","output a 2D slice at the specified depth",automatic) \

#include "../common/argboiler.h"

/*
 * Bundle all the little initializations together
 */
void init(image_t* image, const args_t* args) {
  init_rng(image);
  fesetround(FE_DOWNWARD);
  image->width = args->input_width;
  image->height = args->input_height;
  compute_depth(image);
}

/*
 * Now, we put it all together!
 */
int main(int argc, char** argv) {
  args_t args;
  image_t input;
  int n;

  printf("Parsing command line...\n");
  parse_args(argc, argv, &args);

  printf("\e[AWorm deconvolver [SIPLCA] v0.0.1\ncreated by David Dalrymple\n============================\n\n");

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
  return 0;
}
