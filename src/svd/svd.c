#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <fenv.h>
#include <float.h>
#include <limits.h>
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
  ARG(ARG_INT0,colors,"c","colors","the number of SVD axes to use",20) \
  ARG(ARG_INT0,plot_height,"h","plot-height","height of the plot in each shape layer",50) \
  ARG(ARG_LIT0,no_recons,"r","no-recons","don't bother writing the output image",0) \
  ARG(ARG_LIT0,svd_file,"x","svd-file","check for SVD file; use if exists, write if not",0) \
  ARG(ARG_INT0,masks,"m","masks","how many masks to guess",100) \
  ARG(ARG_DBL0,scale,"f","scale","the scale factor on SVD output interpretation",60.0) \
  ARG(ARG_DBL0,thresh,"t","threshhold","the similarity threshhold for filling masks",0.6) \

#include "../common/argboiler.h"

/*
 * Bundle all the little initializations together
 */
void init(image_t* image, const args_t* args) {
  init_rng(image);
  //fesetround(FE_DOWNWARD);
  image->width = args->input_width;
  image->height = args->input_height;
  image->depth = image->length / image->width / image->height;
}

int main(int argc, char** argv) {
  args_t args;
  image_t input;
  int length;
  int i,j,k;
  point_t point; int *p = point.p;
  SVDRec result;

  parse_args(argc, argv, &args);

  /*
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

  double *U, *Vt, *S;
  char* svdfname;
  FILE* svdf;
  if(args.svd_file) {
    svdfname = malloc(strlen(args.input_filename)+strlen(".svd")+1);
    strcpy(svdfname,args.input_filename);
    strcat(svdfname,".svd");
  }
  if(args.svd_file && (svdf=fopen(svdfname,"r"))) {
    U = malloc(sizeof(double)*args.dimensions*input.width*input.height);
    Vt = malloc(sizeof(double)*args.dimensions*input.width*input.height);
    S = malloc(sizeof(double)*args.dimensions);
    fread(U,sizeof(double),args.dimensions*input.width*input.height,svdf);
    fread(Vt,sizeof(double),args.dimensions*input.width*input.height,svdf);
    fread(S,sizeof(double),args.dimensions,svdf);
  } else {

    step_start("prepping SVD");
    DMat A = svdNewDMat(input.width*input.height,input.depth);
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
    result = svdLAS2A(As,args.dimensions);
    step_end();

    DMat UM = svdTransposeD(result->Ut);
    U = UM->value[0];
    Vt = result->Vt->value[0];
    S = result->S;
    
    if(args.svd_file) {
      svdf=fopen(svdfname,"w");
      fwrite(U,sizeof(double),args.dimensions*input.width*input.height,svdf);
      fwrite(Vt,sizeof(double),args.dimensions*input.width*input.height,svdf);
      fwrite(S,sizeof(double),args.dimensions,svdf);
      fclose(svdf);
    }
  }

  
  step_start("computing cell masks");
  double* masks = malloc(input.width*input.height*args.masks*sizeof(double));
  unsigned short maskbox[args.masks][2][2];
  double* maskcomps = malloc(args.colors*args.masks*sizeof(double));
  double* signal = malloc(input.depth*args.masks*sizeof(double));
  int* taken = malloc(input.width*input.height*sizeof(int));
  for(i=0;i<input.width*input.height;i++) taken[i]=0;
  for(j=0;j<args.masks;j++) {
    double rn;
    int x,y,row;
    double* vec;
    double maxv=0.0;
    double dot=0.0;
    maskbox[j][0][0]=maskbox[j][1][0]=USHRT_MAX;
    maskbox[j][0][1]=maskbox[j][1][1]=0;
    for(i=0;i<input.width*input.height;i++) masks[j*input.width*input.height+i]=0.0;
    do {
      x=gsl_rng_uniform_int(input.r,input.width);
      y=gsl_rng_uniform_int(input.r,input.height);
      rn = gsl_rng_uniform(input.r)/args.scale;
      row = y*input.width+x;
      vec = &(U[row*args.dimensions+1]);
    } while(rn > gsl_stats_mean(vec,1,args.colors-1) || taken[y*input.width+x]);
    int test(int x, int y) {
      if(x>=input.width) return 0;
      if(y>=input.height) return 0;
      if(x<0) return 0;
      if(y<0) return 0;
      int row = y*input.width+x;
      if(taken[row]==j+1) return 0;
      double* tvec = &(U[row*args.dimensions+1]);
      double norm=0.0;
      dot=0.0;
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
      double* vec = &(U[row*args.dimensions+1]);
      double *pixel = (&(masks[j*input.width*input.height+y*input.width+x]));
      double v = dot/args.colors;
      if(v<0) v=0.0; if(v>1) v=1.0;
      if(v>maxv) maxv=v;
      *pixel = v;
      if(x < maskbox[j][0][0]) maskbox[j][0][0]=x;
      if(x > maskbox[j][0][1]) maskbox[j][0][1]=x;
      if(y < maskbox[j][1][0]) maskbox[j][1][0]=y;
      if(y > maskbox[j][1][1]) maskbox[j][1][1]=y;
    }
    floodfill(x,y,test,act);
    if(maxv<0.01) {
      j--;
      printf("."); fflush(stdout);
    } else {
      //Make plot
      int t;
      for(t=0;t<input.depth;t++) {
        int d;
        signal[input.depth*j+t]=0.0;
        for(d=1;d<args.colors;d++) {
          maskcomps[j*args.colors+d] = (S[d])*(U[row*args.dimensions+d]);
          signal[input.depth*j+t]+=maskcomps[j*args.colors+d]*(Vt[d*input.depth+t]);
        }
      }
      printf(" %d ",j); fflush(stdout);
    }
  }
  printf("\n\n");
  step_end();
  

  step_start("creating directory structure");
  const char* outdirname;
  if(args.output_filename[0]=='\0') {
    outdirname = strdup(args.input_filename);
    *(strchrnul(outdirname,'.'))='\0';
  } else {
    outdirname = args.output_filename;
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
  mkdir(outdirname,0777);
  chdir(outdirname);
  mkdir("raw",0777);
  mkdir("raw_c",0777);
  mkdir("recons",0777);
  mkdir("recons_c",0777);
  mkdir("cell",0777);
  FILE *dimsfile = fopen("dims.json","w");
  fprintf(dimsfile,"{width: %d, height: %d, time: %d, cells: %d}\n",input.width,input.height,input.height,args.masks);
  fclose(dimsfile);
  //copy_file("index.html","../../src/svd/index.html");
  step_end();

  step_start("choosing a color scheme");
  //choose a vector for each component
  double* compvects = malloc(sizeof(double)*args.colors*2);
  for(i=1;i<args.colors;i++) {
    compvects[i*2]   = 2*gsl_rng_uniform_pos(input.r)-1.0;
    compvects[i*2+1] = 2*gsl_rng_uniform_pos(input.r)-1.0;
  }
  //choose a color for each pixel
  double* pixcolor = malloc(sizeof(double)*input.height*input.width*2);
  for(p[1]=0;p[1]<input.height;p[1]++) {
    for(p[2]=0;p[2]<input.width;p[2]++) {
      const int row = p[1]*input.width+p[2];
      double x=0.0,y=0.0;
      pixcolor[row*2] = pixcolor[row*2+1] = 0.0;
      for(p[0]=1;p[0]<args.colors;p[0]++) {
        const int col = p[0];
        double v = (U[col+row*args.dimensions])*(S[col]);
        pixcolor[row*2+1] += fabs(v);
        x+=v*compvects[col*2];
        y+=v*compvects[col*2+1];
      }
      pixcolor[row*2] = (atan2(y,x)/TAU)+1/2.0;
      if(isnan(pixcolor[row*2])) pixcolor[row*2]=0.0;
    }
  }
  //choose a color for each mask
  double* maskcolor = malloc(sizeof(double)*args.masks*2);
  for(i=0;i<args.masks;i++) {
    double x=0.0,y=0.0;
    for(j=0;j<args.colors;j++) {
      double v = maskcomps[i*args.colors+j];
      v/=(50*args.colors);
      maskcolor[i*2+1] += fabs(v);
      x+=v*compvects[j*2];
      y+=v*compvects[j*2+1];
    }
    maskcolor[i*2] = (atan2(y,x)/TAU)+1/2.0;
    if(isnan(maskcolor[i*2])) maskcolor[i*2]=0.0;
  }

  double avg_mask_sat=gsl_stats_mean(maskcolor+1,2,args.masks);
  double sd_mask_sat=gsl_stats_sd_m(maskcolor+1,2,args.masks,avg_mask_sat);
  double scale_mask_sat = avg_mask_sat+3*sd_mask_sat;
  for(i=0;i<args.masks;i++) {
    maskcolor[i*2+1] /= scale_mask_sat;
    if(maskcolor[i*2+1]>0.999) maskcolor[i*2+1]=0.999;
  }

  double avg_sat=gsl_stats_mean(pixcolor+1,2,input.width*input.height);
  double sd_sat=gsl_stats_sd_m(pixcolor+1,2,input.width*input.height,avg_sat);
  double scale_sat = avg_sat+3*sd_sat;
  for(i=0;i<input.height*input.width;i++) {
    pixcolor[i*2+1] /= scale_sat;
    if(pixcolor[i*2+1]>0.999) pixcolor[i*2+1]=0.999;
  }

  double avg_mask, var_mask;
  double min_mask=DBL_MAX, max_mask=-DBL_MAX;
  int n_mask;
  for(i=0;i<args.masks*input.width*input.height;i++) {
    double x = masks[i];
    if(x<min_mask) min_mask=x;
    if(max_mask<x) max_mask=x;
    n_mask++;
    double delta = x-avg_mask;
    avg_mask += delta/n_mask;
    var_mask += delta*(x - avg_mask);
  }
  var_mask /= n_mask-1;
  double sd_mask = sqrt(var_mask);
  printf("avg_mask: %lf, sd_mask: %lf\n",avg_mask,sd_mask);
  double scale_mask = max_mask-sd_mask;
  double map_mask(double mask) {
    return (mask-min_mask)/(scale_mask-min_mask);
  }
  for(i=0;i<args.masks*input.width*input.height;i++) {
    masks[i]=map_mask(masks[i]);
    //printf("masks[%d][%d]: %lf\n",i%input.depth,i/input.depth,map_mask(masks[i]));
  }
  for(i=0;i<args.masks;i++) {
    max_mask=-DBL_MAX;
    for(j=0;j<input.width*input.height;j++) {
      double x = masks[i*input.width*input.height+j];
      if(max_mask<x) max_mask=x;
    }
    for(j=0;j<input.width*input.height;j++) {
      masks[i*input.width*input.height+j] /= max_mask;
    }
  }

  double avg_sig, var_sig;
  double min_sig=DBL_MAX, max_sig=-DBL_MAX;
  int n_sig;
  for(i=0;i<args.masks*input.depth;i++) {
    double x = signal[i];
    if(x<min_sig) min_sig=x;
    if(max_sig<x) max_sig=x;
    n_sig++;
    double delta = x-avg_sig;
    avg_sig += delta/n_sig;
    var_sig += delta*(x - avg_sig);
  }
  var_sig /= n_sig-1;
  double sd_sig = sqrt(var_sig);
  printf("avg_sig: %lf, sd_sig: %lf\n",avg_sig,sd_sig);
  double scale_sig = avg_sig;
  double map_sig(double sig) {
    return (sig-min_sig)/(scale_sig-min_sig);
  }
  /*for(i=0;i<args.masks*input.depth;i++) {
    printf("signal[%d][%d]: %lf\n",i%input.depth,i/input.depth,map_sig(signal[i]));
  }*/

  double avg_lum;
  double var_lum;
  int n_lum;
  for(i=0;i<input.depth*input.height*input.width;i++) {
    double x = (double)(((unsigned short*)input.data)[i]);
    n_lum++;
    double delta = x-avg_lum;
    avg_lum += delta/n_lum;
    var_lum += delta*(x - avg_lum);
  }
  var_lum /= n_lum-1;
  double sd_lum = sqrt(var_lum);
  double scale_lum = avg_lum+3*sd_lum;

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
  printf("max point val: %lf\n", max);
  printf("\n");
  step_end();

  unsigned char* pngbuf = malloc(input.width*input.height*4);

  /*
  step_start("writing raw/ pngs");
  chdir("raw");
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
  */

  
  /*
  step_start("writing raw_c/ pngs");
  chdir("raw_c");
  for(k=0;k<input.depth;k++) {
    for(j=0;j<input.height;j++) {
      for(i=0;i<input.width;i++) {
        int n = i+j*input.width;
        unsigned short *pixel = (unsigned short*)(input.data+2*(n+k*input.width*input.height));
        double v = ((double)(*pixel))/scale_lum;
        if(v>1.0) v=1.0;
        else if(v<0.0) v=0.0;
        csl2pix(pngbuf+3*n,pixcolor[n*2],pixcolor[n*2+1],v);
      }
    }
    printf("."); fflush(stdout);
    snprintf(pngfn,pngfn_l,"%04d.png",k);
    export_png(pngfn,input.width,input.height,8+3,pngbuf);
  }
  chdir("..");
  step_end();
  */

  step_start("outputting cell info");
  chdir("cell");
  const int dirname_l=4;
  char dirname[dirname_l+1];
  for(k=0;k<args.masks;k++) {
    snprintf(dirname,dirname_l+1,"%04d",k);
    mkdir(dirname);
    chdir(dirname);

    FILE* bb_json = fopen("bb.json","w");
    fprintf(bb_json,"{xmin: %d, xmax: %d, ymin: %d, ymax: %d}\n",maskbox[k][0][0],maskbox[k][0][1],maskbox[k][1][0],maskbox[k][1][1]);
    fclose(bb_json);

    FILE* signal_json = fopen("signal.json","w");
    fprintf(signal_json,"[\n");
    for(j=0;j<input.depth;j++) {
      fprintf(signal_json,(j==input.depth-1)?"%lf\n":"%lf,\n",map_sig(signal[input.depth*k+j]));
    }
    fprintf(signal_json,"]\n");
    fclose(signal_json);

    for(j=0;j<input.height;j++) {
      for(i=0;i<input.width;i++) {
        int n=i+j*input.width;
        csl2pix(pngbuf+3*n,maskcolor[k*2],maskcolor[k*2+1],masks[n+k*input.height*input.width]);
      }
    }
    printf("."); fflush(stdout);
    export_png("color.png",input.width,input.height,8+3,pngbuf);

    for(j=0;j<input.height;j++) {
      for(i=0;i<input.width;i++) {
        int n=i+j*input.width;
        csl2pix(pngbuf+3*n,0.0,0.0,masks[n+k*input.height*input.width]);
      }
    }
    printf("."); fflush(stdout);
    export_png("gray.png",input.width,input.height,8+3,pngbuf);

    chdir("..");
  }
  chdir("..");
  step_end();


  step_start("writing recons/ pngs");
  step_end();


  step_start("writing recons_c/ pngs");
  chdir("recons_c");
  printf("allocing\n");
  double *recons = malloc(input.depth*input.height*input.width*3*sizeof(double));
#define RECONS(k,j,i,c) recons[k*input.height*input.width*3+j*input.width*3+i*3+c]
  int m;
  printf("starting loop\n");
  for(i=maskbox[m][0][0];i<=maskbox[m][0][1];i++) {
    for(j=maskbox[m][1][0];j<=maskbox[m][1][1];j++) {
      for(k=0;k<input.depth;k++) {
        RECONS(k,j,i,0)=0.0;
        RECONS(k,j,i,1)=0.0;
        RECONS(k,j,i,2)=0.0;
      }
    }
  }
  for(m=0;m<args.masks;m++) {
    for(i=maskbox[m][0][0];i<=maskbox[m][0][1];i++) {
      for(j=maskbox[m][1][0];j<=maskbox[m][1][1];j++) {
        for(k=0;k<input.depth;k++) {
          double x,y,z;
          csl2xyz(&x,&y,&z,maskcolor[m*2],maskcolor[m*2+1],map_sig(signal[m*input.depth+k])*masks[m*input.width*input.height+j*input.width+i]);
          RECONS(k,j,i,0)+=x;
          RECONS(k,j,i,1)+=y;
          RECONS(k,j,i,2)+=z;
        }
      }
    }
  }
  for(k=0;k<input.depth;k++) {
    for(j=0;j<input.height;j++) {
      for(i=0;i<input.width;i++) {
        xyz2pix(pngbuf+3*(i+j*input.width),RECONS(k,j,i,0),RECONS(k,j,i,1),RECONS(k,j,i,2));
      }
    }
    printf("."); fflush(stdout);
    snprintf(pngfn,pngfn_l,"%04d.png",k);
    export_png(pngfn,input.width,input.height,8+3,pngbuf);
  }
  chdir("..");
  step_end();

  return 0;
}
