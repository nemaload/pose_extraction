#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <fenv.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

#include "../common/svdlib.h"

#include "../common/util.h"
#include "../common/debug.h"

#define ARGBOILER(ARG) \
  ARG(ARG_FIL1,input_filename,NULL,NULL,"the input image (as raw data)","") \
  ARG(ARG_INT1,input_width,"w","width","the width of each input image slice",-1) \
  ARG(ARG_INT1,input_height,"h","height","the height of each input image slice",-1) \
  ARG(ARG_FIL0,output_filename,"o","output","name of the output folder [default: truncated input]","") \
  ARG(ARG_INT0,dimensions,"d","dim","the number of SVD axes to extract",100) \
  ARG(ARG_INT0,colors,"c","colors","the number of SVD axes to colorize",20) \
  ARG(ARG_INT0,plot_height,"h","plot-height","height of the plot in each shape layer",50) \
  ARG(ARG_LIT0,no_recons,"r","no-recons","don't bother writing the output image",0) \
  ARG(ARG_INT0,shapes,"s","shapes","how many shapes to guess",100) \
  ARG(ARG_DBL0,scale,"f","scale","the scale factor on SVD output interpretation",60.0) \
  ARG(ARG_DBL0,thresh,"t","threshhold","the similarity threshhold for filling shapes",0.6) \

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
  int i,j,k;

  parse_args(argc, argv, &args);

  
  step_start("L*a*b* test");
  unsigned char* lab_test = malloc(args.input_width*args.input_width*3);
  for(i=0;i<args.input_width;i++) {
    for(j=0;j<args.input_width;j++) {
      unsigned char* ptr = lab_test+(j+i*args.input_width)*3;
      double x,y,l,a,b;
      x=(double)j;
      //x = 20.0*(double)(j/20);
      y=(double)i;
      x/=((double)args.input_width)/1.0;
      y/=((double)args.input_width)/1.0;
      cl2pix(ptr,x,y);
    }
  }
  export_png("lab_test.png",args.input_width,args.input_width,8+3,lab_test);
  step_end();

  step_start("Equivalent RGB test");
  unsigned char* rgb_test = malloc(args.input_width*args.input_width*3);
  for(i=0;i<args.input_width;i++) {
    for(j=0;j<args.input_width;j++) {
      unsigned char* ptr = rgb_test+(j+i*args.input_width)*3;
      double x,y,l,a,b;
      x=(double)j;
      y=(double)i;
      //x = 20.0*(double)(j/20);
      //y = 20.0*(double)(i/20);
      y=(double)i;
      x/=((double)args.input_width)/1.0;
      y/=((double)args.input_width)/1.0;
      hsv2pix(ptr,x,0.5,y);
    }
  }
  export_png("rgb_test.png",args.input_width,args.input_width,8+3,rgb_test);
  step_end();

  return 0;
  //  end */

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
  for(p[0]=0;p[0]<input.depth;p[0]++) {
    const int col = p[0];
    for(p[1]=0;p[1]<input.height;p[1]++) {
      for(p[2]=0;p[2]<input.width;p[2]++) {
        const int row = p[1]*input.width+p[2];
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

  /*
  step_start("reconstructing image");
  if(args.no_recons) {
    printf("--skipping--\n");
  } else {
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
  

  DMat U = svdTransposeD(result->Ut);
  
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
      double r,g,b;
      double* vec = &(U->value[row][1]);
      r=gsl_stats_mean(vec,1,args.colors-1);
      g=gsl_stats_sd_m(vec,1,args.colors-1,r);
      b=g; r=g;
      r*=args.scale;g*=args.scale;b*=args.scale;
      if(r<0) r=0.0;
      if(r>1) r=1.0;
      if(g<0) g=0.0;
      if(g>1) g=1.0;
      if(b<0) b=0.0;
      if(b>1) b=1.0;
      unsigned char *pixel = (&(shapes[(p[1]*input.width+p[2])*3]));
      printf(" %lf, %lf, %lf -> %d, %d, %d\n",r,g,b,(unsigned char)(r*255.0),(unsigned char)(g*255.0),(unsigned char)(b*255.0));
      pixel[0]=(unsigned char)(r*255.0);
      pixel[1]=(unsigned char)(g*255.0);
      pixel[2]=(unsigned char)(b*255.0);
    }
  }
  step_end();*/

  /*
  step_start("isolating shapes");
  unsigned short* shapes = open_mmapped_file_write(sfilename,input.width*(input.height+args.plot_height)*args.shapes*2);
  int* taken = calloc(input.width*input.height,sizeof(int));
  for(j=0;j<args.shapes;j++) {
    double rn;
    int x,y,row;
    double* vec;
    double maxv=0.0;
    do {
      x=gsl_rng_uniform_int(input.r,input.width);
      y=gsl_rng_uniform_int(input.r,input.height);
      rn = gsl_rng_uniform(input.r)/args.scale;
      row = y*input.width+x;
      vec = &(U->value[row][1]);
    } while(rn > gsl_stats_mean(vec,1,args.colors-1) || taken[y*input.width+x]);
    int test(int x, int y) {
      if(x>=input.width) return 0;
      if(y>=input.height) return 0;
      if(x<0) return 0;
      if(y<0) return 0;
      int row = y*input.width+x;
      if(taken[row]==j+1) return 0;
      double* tvec = &(U->value[row][1]);
      double dot=0.0;
      double norm=0.0;
      int i;
      for(i=0;i<args.colors-1;i++) {
        dot+=(tvec[i]*vec[i]);
        norm+=vec[i]*vec[i];
      }
      dot*=args.scale*args.scale;
      if(dot>args.thresh) {
        return 1;
      } else {
        return 0;
      }
    }
    void act(int x, int y) {
      int row = y*input.width+x;
      taken[row]=j+1;
      double* vec = &(U->value[row][1]);
      unsigned short *pixel = (&(shapes[j*input.width*(input.height+args.plot_height)+y*input.width+x]));
      double v = gsl_stats_mean(vec,1,args.colors-1)*args.scale;
      if(v<0) v=0.0; if(v>1) v=1.0;
      if(v>maxv) maxv=v;
      *pixel = (unsigned short)(v*65535.0);
    }
    floodfill(x,y,test,act);
    if(maxv<0.01) {
      j--;
      printf("retrying...\n");
    } else {
      //Make plot
      int px,py;
      double ttymin=100000000, ttymax=-100000000;
      double* tty = malloc(sizeof(double)*input.width);
      for(px=0;px<input.width;px++) {
        double td = ((double)px/(double)input.width)*((double)input.depth);
        int ti = (int)td;
        if(ti<0) ti=0; if(ti>input.depth-1) ti=input.depth-1;
        int d;
        tty[px]=0.0;
        for(d=1;d<args.colors;d++) {
          tty[px]+=(result->S[d])*(U->value[row][d])*(result->Vt->value[d][ti]);
        }
        if(tty[px]<ttymin) ttymin=tty[px];
        if(tty[px]>ttymax) ttymax=tty[px];
      }
      for(px=0;px<input.width;px++) {
        int iy = (int)((tty[px]-ttymin)*(((double)args.plot_height)/(ttymax-ttymin)));
        for(py=0;py<args.plot_height;py++) {
          double v;
          if(py == iy) {
            v = 0.0;
          } else {
            v = 1.0;
          }
          unsigned short *pixel = (&(shapes[j*input.width*(input.height+args.plot_height)+(input.height+args.plot_height-1-iy)*input.width+px]));
          *pixel = (unsigned short)(v*65535.0);
        }
      }
      printf("==========================done with %d!\n",j);
    }
  }
  step_end();
  */

  step_start("creating directory structure");
  const char* dirname;
  if(args.output_filename[0]=='\0') {
    dirname = strdup(args.input_filename);
    *(strchrnul(dirname,'.'))='\0';
  } else {
    dirname = args.output_filename;
  }
  /*
   * dirname/
   *   index.html
   *   dims.json
   *   raw/
   *     0001.png ... NNNN.png
   *   raw_c/
   *     0001.png ... NNNN.png
   *   recons/
   *     0001.png ... NNNN.png
   *   recons_c/
   *     0001.png ... NNNN.png
   *   cell/
   *     0001/ ... MMMM/
   *       bb.json
   *       signal.json
   *       gray.png
   *       color.png
   */
  mkdir(dirname,0777);
  chdir(dirname);
  mkdir("raw",0777);
  mkdir("raw_c",0777);
  mkdir("recons",0777);
  mkdir("recons_c",0777);
  mkdir("cell",0777);
  copy_file("index.html","../../src/svd/index.html");
  step_end();


  step_start("writing raw/ pngs");
  chdir("raw");
  int pngfn_l = strlen("0000.png")+1;
  char* pngfn = malloc(pngfn_l);

  //Find maximum point value to scale by.
  double max=0.0;
  for(i=0;i<input.depth*input.height*input.width;i++) {
    if(max<((unsigned short*)input.data)[i]) {
      max = (double)(((unsigned short*)input.data)[i]);
    }
  }
  double scale = max/255.0;

  unsigned char* pngbuf = malloc(input.width*input.height);
  for(k=0;k<input.depth;k++) {
    for(j=0;j<input.height;j++) {
      for(i=0;i<input.width;i++) {
        unsigned short *pixel = (unsigned short*)(input.data+2*(i+j*input.width+k*input.width*input.height));
        pngbuf[i+j*input.width]=(unsigned char)(((double)(*pixel))/scale);
      }
    }
    snprintf(pngfn,pngfn_l,"%04d.png",k);
    export_png(pngfn,input.width,input.height,8+1,pngbuf);
  }
  chdir("..");
  step_end();


  step_start("writing recons/ pngs");

  step_end();


  step_start("writing raw_c/ pngs");

  step_end();


  step_start("writing recons_c/ pngs");

  step_end();

  
  step_start("making cell/ directories");

  step_end();


  step_start("writing cell/ JSON");

  step_end();


  step_start("writing cell/ pngs");

  step_end();

  return 0;
}
