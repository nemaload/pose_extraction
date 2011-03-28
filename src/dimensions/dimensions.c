#include <stdint.h>
#include <complex.h>
#include <fftw3.h>
#include "../common/util.h"
#include "../common/debug.h"

#define ARGBOILER(ARG) \
  ARG(ARG_FIL1,input_filename,NULL,NULL,"the input image (as raw data)","") \
  ARG(ARG_INT0,max_dim,"m","max-dim","the maximum number of pixels in any dimension",10000) \
  ARG(ARG_INT0,fourier_n,"f","fourier-n","the number of pixels to consider in the FFT",100000) \

#include "../common/argboiler.h"

int main(int argc, char** argv) {
  args_t args;
  uint16_t* data;
  int length;
  int i;

  printf("Parsing command line...\n");
  parse_args(argc, argv, &args);

  printf("Dimension estimator\nDavid Dalrymple\n==================\n");

  step_start("mmapping file");
  data = open_mmapped_file_read(args.input_filename, &length);
  length/=2;
  if(data==NULL) {
    printf("Couldn't open file\n");
    return 1;
  }
  step_end();

  step_start("planning FFT");
  int fourier_n = args.fourier_n;
  if(length<fourier_n) fourier_n=length;
  double* data_d = fftw_malloc(sizeof(double)*fourier_n);
  fftw_complex* fft_data = fftw_malloc(sizeof(fftw_complex)*2*(fourier_n/2+1));
  fftw_set_timelimit(0.5);
  fftw_plan fft = fftw_plan_dft_r2c_1d(fourier_n,data_d,fft_data,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
  step_end();

  step_start("executing FFT");
  for(i=0;i<fourier_n;i++) {
    data_d[i] = (double)data[i];
  }
  fftw_execute(fft);
  step_end();

  step_start("conjugate");
  for(i=0;i<2*(fourier_n/2+1);i++) {
    fft_data[i] = fft_data[i] * conj(fft_data[i]);
  }
  step_end();

  step_start("inverse FFT");
  fftw_plan ifft = fftw_plan_dft_c2r_1d(fourier_n,fft_data,data_d,FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
  fftw_execute(ifft);
  step_end();

  step_start("sorting");
  int nind = args.max_dim;
  if(fourier_n<nind) nind=fourier_n;
  int* indices = malloc(sizeof(int)*nind);
  for(i=0;i<nind;i++)
    indices[i]=i+1;
  int compare_indices_rev(const void *a, const void *b) {
    const int* ia = (const int*)a;
    const int* ib = (const int*)b;
    return (data_d[*ia] < data_d[*ib]) - (data_d[*ia] > data_d[*ib]);
  }
  qsort(indices, nind, sizeof(int), compare_indices_rev);
  step_end();

  for(i=0;i<50;i++) {
    printf("#%d: %d, %lf\n",i+1,indices[i],data_d[indices[i]]/(double)fourier_n);
  }

  printf("\nOkay, I'm going to assume the x dimension is %d.\n",indices[0]);
  int xdim = indices[0];

  int j=1;
  i=0;
  while(j<=nind) {
    int proposed_xy = xdim*j;
    if(proposed_xy > length) {
      break;
    }
    if(length%proposed_xy!=0) {
      j++;
      continue;
    }
    printf("indices[%d]=%d\n",i,proposed_xy);
    indices[i++]=proposed_xy;
    j++;
  }
  qsort(indices, i, sizeof(int), compare_indices_rev);

  for(j=0;j<i;j++) {
    printf("#%d: %dx%dx%d, %lf\n",j+1,xdim,indices[j]/xdim,length/indices[j],data_d[indices[j]]/(double)fourier_n);
  }

  return 0;
}
