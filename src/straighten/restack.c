/*
 * [Step 8]
 * Actually straighten the image by restacking along a spline defined by the backbone.
 */
struct restack_shared {
#define RESTACK_STRUCT(c) \
  gsl_spline* spl##c; \
  gsl_interp_accel* accel##c;
  FOREACH3(RESTACK_STRUCT)
  unsigned short* data;
  int src_height;
  int src_width;
  int src_depth;
  unsigned short* new_data;
  int dst_height;
  int dst_width;
  int dst_depth;
  int extension;
  int no_interpolate;
  int output_slice;
  int three_d_output;
  int* slices_done;
  pthread_mutex_t* slices_done_mutex;
};

typedef struct {
  struct restack_shared sh;
  int j_start;
  int j_end;
} restack_workunit_t;

void* restack_worker(void* workunit);

void restack_image(image_t* dst, const image_t* src, const args_t* args, dpoint_t* backbone, int n) {
  struct restack_shared sh;
  /*
   * First, we define the multidimensional spline, by walking along the control
   * points, measuring the total distance walked, and using that as the parameter
   * for three single-dimensional splines.
   */
  double* xa = malloc(sizeof(double)*n);
#define PTS_LIST_ALLOC(c) \
  double* y##c = malloc(sizeof(double)*n);
  FOREACH3(PTS_LIST_ALLOC)
  int i;
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
  sh.spl##c = gsl_spline_alloc(gsl_interp_cspline,n); \
  sh.accel##c = gsl_interp_accel_alloc(); \
  gsl_spline_init(sh.spl##c,xa,y##c,n);
  FOREACH3(GSL_INIT)

  const char* filename;
  const char* suffix=".out";
  sh.three_d_output=0;
  if(args->output_filename[0]=='\0') {
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
    printf("Output width manually set: %d\n",dst->width);
  }
  if(args->output_extension==-1) {
    sh.extension = dst->width;
  } else {
    sh.extension = args->output_extension;
  }
  if(args->output_height==-1) {
    dst->height=(int)xa[n-1]+2*sh.extension;
    printf("Output height automatically determined: %d\n",dst->height);
  } else {
    dst->height = args->output_height;
    printf("Output height manually set: %d\n",dst->height);
  }
  if(args->output_slice==-1) {
    sh.three_d_output=1;
    if (src->depth > 1)
      dst->depth = (src->depth/2-1)*2;
    else
      dst->depth = src->depth;
    //TODO: output depth option?
    printf("Output depth (same as input depth): %d\n",dst->depth);
  } else {
    dst->depth = 1;
    printf("Output is a single slice at %d (depth is 1)\n",args->output_slice);
  }
  printf("\n");

  int length;
  length = dst->width*dst->height*dst->depth*2;
  dst->data = (open_mmapped_file_write(filename,length));

  sh.data = (unsigned short*)src->data;
  sh.src_height = src->height;
  sh.src_width = src->width;
  sh.src_depth = src->depth;
  
  sh.new_data = (unsigned short*)dst->data;
  sh.dst_height = dst->height;
  sh.dst_width = dst->width;
  sh.dst_depth = dst->depth;

  sh.no_interpolate = args->no_interpolate;
  sh.output_slice = args->output_slice;
  int slices_done=0;
  sh.slices_done=&slices_done;

  int n_threads = args->n_threads;
  pthread_t* thread = malloc(sizeof(pthread_t)*n_threads);
  pthread_mutex_t* slices_done_mutex = malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(slices_done_mutex,NULL);
  sh.slices_done_mutex=slices_done_mutex;

  for(i=0;i<n_threads;i++) {
    restack_workunit_t* wu=malloc(sizeof(restack_workunit_t));
    wu->sh=sh;
    wu->j_start=i*(sh.dst_height/n_threads);
    if(i==n_threads-1) {
      wu->j_end=sh.dst_height;
    } else {
      wu->j_end=(i+1)*(sh.dst_height/n_threads);
    }
    pthread_create(&thread[i],NULL,restack_worker,(void*)wu);
  }
  for(i=0;i<n_threads;i++) {
    void* status;
    pthread_join(thread[i], &status);
    free(status);
  }
  printf("\n\nSyncing to disk...\n");
  msync(sh.new_data,length,MS_SYNC);
}

void* restack_worker(void* workunit) {
  restack_workunit_t* wu=(restack_workunit_t*)workunit;
  struct restack_shared sh=wu->sh;
  int j_start = wu->j_start;
  int j_end = wu->j_end;
  int i,j,k;
  int two_d_input=0;
  if(sh.src_depth==1){two_d_input=1;}

  for(j=j_start;j<j_end;j++) {
    double magnitude;
#define GET_P_D(c) \
    double p##c = gsl_spline_eval(sh.spl##c,j-sh.extension,sh.accel##c); \
    double d##c = gsl_spline_eval_deriv(sh.spl##c,j-sh.extension,sh.accel##c); \
    double dx##c; \
    double dz##c; \
    double pz##c;
    FOREACH3(GET_P_D)

    dx0=0;
    dx1=d2;
    dx2=-d1;
    magnitude = sqrt(dx1*dx1+dx2*dx2);
    dx1/=magnitude;
    dx2/=magnitude;

    dz0=dx1*d2-dx2*d1;
    dz1=dx2*d0;
    dz2=-dx1*d0;
    magnitude = sqrt(dz0*dz0+dz1*dz1+dz2*dz2);
    dz0/=magnitude;
    dz1/=magnitude;
    dz2/=magnitude;

#define COPY_P_Z(c) \
    pz##c = p##c;
    FOREACH3(COPY_P_Z)
    if(!sh.three_d_output) {
#define INC_P_Z_S(c) \
      pz##c += (sh.output_slice-sh.src_depth/2)*dz##c;
      FOREACH3(INC_P_Z_S)
    } else {
#define INIT_P_Z(c) \
      pz##c -= (((double)sh.dst_depth)/2)*dz##c;
      FOREACH3(INIT_P_Z)
    }

    for(k=0;k<sh.dst_depth;k++) {
#define INIT_P_X(c) \
      p##c = pz##c - dx##c*(((double)sh.dst_width)/2);
      FOREACH3(INIT_P_X)
      for(i=0;i<sh.dst_width;i++) {
#define P_INT(c) \
        int pi##c=lrint(p##c);
        FOREACH3(P_INT)
        unsigned short pixel = 0;
        if(pi0>=0&&pi1>=0&&pi2>=0&&(pi0<sh.src_depth-1||two_d_input)&&pi1<sh.src_height-1&&pi2<sh.src_width-1) {
          if(sh.no_interpolate) {
            pixel = pixel_get_(((unsigned short*)sh.data),pi0,pi1,pi2,sh.src_width,sh.src_height);
          } else {
#define TRUNC_P(c) \
            double t##c = p##c - pi##c;
            FOREACH3(TRUNC_P)
            double pixel_d=0;
#define ONE_MINUS(a,x) (1-a+(2*a-1)*x)
#define GET_CORNER(z,y,x) \
            pixel_d += (ONE_MINUS(z,t0)*ONE_MINUS(y,t1)*ONE_MINUS(x,t2))*(sh.data)[((pi0)+z)*sh.src_height*sh.src_width+((pi1)+y)*sh.src_width+((pi2)+z)]
            GET_CORNER(0,0,0);
            GET_CORNER(0,0,1);
            GET_CORNER(0,1,0);
            GET_CORNER(0,1,1);
            if(!two_d_input) {
              GET_CORNER(1,0,0);
              GET_CORNER(1,0,1);
              GET_CORNER(1,1,0);
              GET_CORNER(1,1,1);
            }

            pixel = (unsigned short) pixel_d;
          }
        }
        unsigned short* nd_pixel = (&(sh.new_data[k*sh.dst_width*sh.dst_height+j*sh.dst_width+i]));
        *nd_pixel=pixel;
#define INC_P_X(c) \
        p##c += dx##c;
        FOREACH3(INC_P_X)
      }
      //sched_yield();
#define INC_P_Z(c) \
      pz##c += dz##c;
      FOREACH3(INC_P_Z)
    }
    pthread_mutex_lock(sh.slices_done_mutex);
    (*sh.slices_done)++;
    pthread_mutex_unlock(sh.slices_done_mutex);
    progress_t(j+1,(*sh.slices_done),sh.dst_height,"planes");
  }
  return workunit;
}
