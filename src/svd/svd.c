#include <stdint.h>
#include <fenv.h>
#include <svdlib.h>
#include <string.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <math.h>
#include "../common/util.h"
#include "../common/debug.h"

#define ARGBOILER(ARG) \
  ARG(ARG_FIL1,input_filename,NULL,NULL,"the input image (as raw data)","") \
  ARG(ARG_INT1,input_width,"w","width","the width of each input image slice",-1) \
  ARG(ARG_INT1,input_height,"h","height","the height of each input image slice",-1) \
  ARG(ARG_INT0,dimensions,"d","dim","the number of SVD axes to extract",100) \
  ARG(ARG_INT0,colors,"c","colors","the number of SVD axes to colorize",20) \
  ARG(ARG_INT0,plot_height,"h","plot-height","height of the plot in each shape layer",50) \
  ARG(ARG_LIT0,no_recons,"r","no-recons","don't bother writing the output image",0) \
  ARG(ARG_FIL0,output_filename,"o","output","the output image (as raw data) [default: input+'.out']","") \
  ARG(ARG_INT0,shapes,"s","shapes","how many shapes to guess",100) \
  ARG(ARG_FIL0,shapes_filename,"","shape-file","the shape stack (as raw data) [default: input+'.shapes']","") \
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

struct pqn {
  int x;
  int y;
  struct pqn* next;
};
struct pq {
  struct pqn* head;
  struct pqn* tail;
};
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
//When you "test" make sure you haven't already "act"ed!
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
    sfilename = args.shapes_filename;
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

  DMat U = svdTransposeD(result->Ut);
  /*
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

  return 0;
}
