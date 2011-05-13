#include <stdint.h>
#include <fenv.h>
#include <svdlib.h>
#include <string.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include "../common/util.h"
#include "../common/debug.h"

#define ARGBOILER(ARG) \
  ARG(ARG_FIL1,input_filename,NULL,NULL,"the input image (as raw data)","") \
  ARG(ARG_INT1,input_width,"w","width","the width of each input image slice",-1) \
  ARG(ARG_INT1,input_height,"h","height","the height of each input image slice",-1) \
  ARG(ARG_INT0,dimensions,"d","dim","the number of SVD axes to extract",100) \
  ARG(ARG_INT0,colors,"c","colors","the number of SVD axes to colorize",20) \
  ARG(ARG_FIL0,output_filename,"o","output","the output image (as raw data) [default: input+'.out']","") \
  ARG(ARG_FIL0,shapes_filename,"s","shapes","the shape stack (as raw data) [default: input+'.shapes']","") \
  ARG(ARG_LIT0,no_recons,"r","no-recons","don't bother writing the output image",0) \

#include "../common/argboiler.h"

/*
 * Bundle all the little initializations together
 */
void init(image_t* image, const args_t* args) {
  init_rng(image);
  fesetround(FE_DOWNWARD);
  image->width = args->input_width;
  image->height = args->input_height;
  image->depth = image->length / image->width / image->height;
}

int main(int argc, char** argv) {
  args_t args;
  image_t input;
  int length;
  int i,j;

  printf("Parsing command line...\n");
  parse_args(argc, argv, &args);

  printf("Singular value decomposer\nDavid Dalrymple\n==================\n");

  step_start("mmapping file");
  input.data = open_mmapped_file_read(args.input_filename, &length);
  input.length = length/2;

  if(input.data==NULL) {
    printf("Couldn't open file\n");
    return 1;
  }
  step_end();

  init(&input, &args);

  step_start("prepping SVD");
  DMat A = svdNewDMat(input.width*input.height,input.depth);
  point_t point; int *p = point.p;
  for(p[1]=0;p[1]<input.height;p[1]++) {
    for(p[2]=0;p[2]<input.width;p[2]++) {
      const int row = p[1]*input.width+p[2];
      for(p[0]=0;p[0]<input.depth;p[0]++) {
        const int col = p[0];
        A->value[row][col]=(double)pixel_get(&input,point);
      }
    }
  }
  step_end();
  half_step_start("converting DMat to SMat");
  SMat As = svdConvertDtoS(A);
  step_end();

  step_start("executing SVD");
  SVDRec result = svdLAS2A(As,args.dimensions);
  step_end();

  const char* ofilename;
  const char* osuffix = ".out";
  if(args.output_filename[0]=='\0') {
    char * filename_ = calloc(strlen(args.input_filename)+strlen(osuffix),1);
    strcat(filename_,args.input_filename);
    strcat(filename_,osuffix);
    ofilename=filename_;
  } else {
    ofilename = args.output_filename;
  }

  step_start("reconstructing image");
  if(args.no_recons) {
    printf("--skipping--\n");
  } else {
    unsigned short* output = open_mmapped_file_write(ofilename,length);
    for(p[1]=0;p[1]<input.height;p[1]++) {
      for(p[2]=0;p[2]<input.width;p[2]++) {
        const int row = p[1]*input.width+p[2];
        for(p[0]=0;p[0]<input.depth;p[0]++) {
          const int col = p[0];
          double v = 0.0;
          unsigned short *pixel = (&(output[p[0]*input.width*input.height+p[1]*input.width+p[2]]));
          for(i=0;i<args.dimensions;i++) {
            v+=(result->S[i])*(result->Ut->value[i][row])*(result->Vt->value[i][col]);
          }
          *pixel = (unsigned short)v;
        }
      }
    }
  }
  step_end();

  const char* sfilename;
  const char* ssuffix = ".shapes";
  if(args.shapes_filename[0]=='\0') {
    char * filename_ = calloc(strlen(args.input_filename)+strlen(ssuffix),1);
    strcat(filename_,args.input_filename);
    strcat(filename_,ssuffix);
    sfilename=filename_;
  } else {
    sfilename = args.output_filename;
  }
  /*
  step_start("reconstructing shapes");
  unsigned short* shapes = open_mmapped_file_write(sfilename,input.width*input.height*args.dimensions*2);
  for(p[1]=0;p[1]<input.height;p[1]++) {
    for(p[2]=0;p[2]<input.width;p[2]++) {
      const int row = p[1]*input.width+p[2];
      for(p[0]=0;p[0]<args.dimensions;p[0]++) {
        const int col = p[0];
        unsigned short *pixel = (&(shapes[p[0]*input.width*input.height+p[1]*input.width+p[2]]));
        *pixel = (unsigned short)((result->Ut->value[col][row])*(result->S[col]));
      }
    }
  }
  step_end();
  */
  step_start("coloring shapes");
  unsigned char* shapes = open_mmapped_file_write(sfilename,input.width*input.height*3);
  double* colors = malloc(args.colors*3*sizeof(double));
  for(i=1;i<args.colors;i++) {
    double r,g,b;
    r=gsl_rng_uniform_pos(input.r);
    g=gsl_rng_uniform_pos(input.r);
    b=gsl_rng_uniform_pos(input.r);
    //double y = 0.2126*r+0.7152*g+0.0722*b;
    double y = 0.33*r+0.33*g+0.33*b;
    double yc = 0.02/y;
    r*=yc; g*=yc; b*=yc;
    colors[i*3]=r;
    colors[i*3+1]=g;
    colors[i*3+2]=b;
  }
  printf("Assigned colors\b\n");
  for(p[1]=0;p[1]<input.height;p[1]++) {
    for(p[2]=0;p[2]<input.width;p[2]++) {
      const int row = p[1]*input.width+p[2];
      printf("Pixel %d,%d:",p[2],p[1]);
      double r=0.33, g=0.33, b=0.33;
      for(p[0]=1;p[0]<args.colors;p[0]++) {
        const int col = p[0];
        double v = (result->Ut->value[col][row])*(result->S[col]);
        v/=(50*args.colors);
        r+=(v*colors[p[0]*3]);
        g+=(v*colors[p[0]*3+1]);
        b+=(v*colors[p[0]*3+2]);
      }
      double mean = (r+g+b)/3.0;
      double var = ((r*r+g*g+b*b)/3.0)-mean*mean;
      double sd = sqrt(var);
      sd*=args.colors;
      printf(" sd(%lf);",sd);
      r*=sd; g*=sd; b*=sd;
      if(r<0) r=0.0;
      if(r>1) r=1.0;
      if(g<0) g=0.0;
      if(g>1) g=1.0;
      if(b<0) b=0.0;
      if(b>1) b=1.0;
      mean = (r+g+b)/3.0;
      if(mean<0.2) {
        r=0; g=0; b=0;
      }
      unsigned char *pixel = (&(shapes[(p[1]*input.width+p[2])*3]));
      printf(" %lf, %lf, %lf -> %d, %d, %d\n",r,g,b,(unsigned char)(r*255.0),(unsigned char)(g*255.0),(unsigned char)(b*255.0));
      pixel[0]=(unsigned char)(r*255.0);
      pixel[1]=(unsigned char)(g*255.0);
      pixel[2]=(unsigned char)(b*255.0);
    }
  }

  return 0;
}
