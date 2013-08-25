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

struct refine_shared {
#ifdef X11
  pthread_mutex_t* g2_mutex;
  int g2;
#endif
  pthread_mutex_t* backbone_mutex;
  double* total_brightness;
  dpoint_t* weighted_sum;

  const kdtree_t* kdtree;
  const point_t* sample;
  const image_t* image;
  const args_t* args;
  const dpoint_t* backbone;
  int n;
};

typedef struct {
  struct refine_shared* sh;
  int i_start;
  int i_end;
} account_point_workunit_t;

void* account_point_worker(void* workunit) {
  account_point_workunit_t* wu=(account_point_workunit_t*)workunit;
  struct refine_shared sh=*(wu->sh);
  int i_start = wu->i_start;
  int i_end = wu->i_end;
  int i;
  const image_t* image = sh.image;
  double image_scale = sh.args->image_scale;
  int use_brightness = sh.args->use_brightness;
  int spread_voronoi = sh.args->spread_voronoi;
  int n = sh.n;

  for(i=i_start;i<i_end;i++) {
    unsigned short brightness = 1;
    int nn = kdtree_search(sh.kdtree, sh.sample[i])->location.index;
    if(use_brightness) {
      brightness = pixel_get(sh.image,sh.sample[i]);
    }
#define ADD_WSUM(c) sh.weighted_sum[nn].p[c]+=brightness*sh.sample[i].p[c];
#define ACCOUNT_POINT_ \
    pthread_mutex_lock(&sh.backbone_mutex[nn]); \
    FOREACH3(ADD_WSUM) \
    sh.total_brightness[nn]+=brightness; \
    pthread_mutex_unlock(&sh.backbone_mutex[nn]);
#ifndef X11
#define ACCOUNT_POINT ACCOUNT_POINT_
#else
#define G2_TRANSFORM(x) g2,\
                                x.p[2]/image_scale,\
    (image->height/image_scale)-x.p[1]/image_scale
#define ACCOUNT_POINT ACCOUNT_POINT_;\
    pthread_mutex_lock(sh.g2_mutex);\
    g2_pen(sh.g2,7);\
    g2_move(sh.G2_TRANSFORM(sh.backbone[nn]));\
    g2_line_to(sh.G2_TRANSFORM(sh.sample[i]));\
    pthread_mutex_unlock(sh.g2_mutex);
#endif
    ACCOUNT_POINT;
    if(spread_voronoi) {
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
  return workunit;
}

int refine_backbone(const image_t* image, point_t* sample, const args_t* args, dpoint_t* backbone, int n) {
  int i;
  double iter_delta = INFINITY;
  double iter_delta_init = INFINITY;
  dpoint_t* backbone_new = malloc(n*sizeof(dpoint_t));
  struct refine_shared sh;
  sh.total_brightness = malloc(n*sizeof(double));
  sh.weighted_sum = malloc(n*sizeof(dpoint_t));
  sh.args = args;
  sh.image = image;
  sh.sample = sample;
  sh.backbone = backbone;
  sh.backbone_mutex = malloc(n*sizeof(pthread_mutex_t));
  sh.n = n;
  for(i=0;i<n;i++) {
    pthread_mutex_init(&sh.backbone_mutex[i],NULL);
  }
  int iterations = 0;
  int delta_history = 0;
  char spinner[] = "-\\|/";
  char* dh_str = malloc(args->delta_history+2);
  dh_str[args->delta_history+1]='\0';
#ifdef X11
  sh.g2=g2_open_X11(image->width/args->image_scale,image->height/args->image_scale);
  g2_set_auto_flush(sh.g2,0);
  sh.g2_mutex = malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(sh.g2_mutex,NULL);
#endif
  printf("\e[s");
  int n_threads = args->n_threads;
  int refine_sample_size = args->refine_sample_size;
  double image_scale = args->image_scale;
  while(delta_history < args->delta_history && iterations < args->restart_iterations) {
    memset(sh.total_brightness,0,n*sizeof(double));
    memset(sh.weighted_sum,0,n*sizeof(dpoint_t));

    kdtree_t* kdtree = kdtree_build(backbone,n);
    sh.kdtree = kdtree;

    pthread_t* thread = malloc(n_threads*sizeof(pthread_t));
    for(i=0;i<n_threads;i++) {
      account_point_workunit_t* wu=malloc(sizeof(account_point_workunit_t));
      wu->i_start=(refine_sample_size/n_threads)*i;
      if(i==args->n_threads-1) {
        wu->i_end=refine_sample_size-1;
      } else {
        wu->i_end=(refine_sample_size/n_threads)*(i+1);
      }
      wu->sh=&sh;
      pthread_create(&thread[i],NULL,account_point_worker,(void*)wu);
    }
    for(i=0;i<n_threads;i++) {
      void* status;
      pthread_join(thread[i], &status);
      free(status);
    }

    kdtree_free(kdtree);
#ifdef X11
#define G2_DRAW_BACKBONE \
    g2_pen(sh.g2,0);\
    g2_filled_rectangle(sh.g2,0,0,(double)image->width/image_scale,(double)image->height/image_scale);\
    g2_move(sh.G2_TRANSFORM(backbone[0]));\
    g2_pen(sh.g2,19);\
    for(i=1;i<n;i++) {\
      g2_line_to(sh.G2_TRANSFORM(backbone[i]));\
    }
    G2_DRAW_BACKBONE

    g2_flush(sh.g2);
    usleep(1e3);
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
      if(sh.total_brightness[i]!=0) {
        for(c=0;c<3;c++) {
          backbone_new[i].p[c]=(args->alpha*(sh.weighted_sum[i].p[c]/sh.total_brightness[i]));
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

    memcpy(backbone,backbone_new,n*sizeof(dpoint_t));

    replace_in_sample(image,sample,args->refine_refresh_size,args->refine_sample_size);

    memset(dh_str,32,args->delta_history+1);
    memset(dh_str,'=',delta_history);
    dh_str[delta_history]='>';
    printf("\e[u[\e[1;32;44m%c\e[m] iterations: %d; delta: %s%lf\t\e[32;44m%s\e[m]",\
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
  if(iterations >= args->restart_iterations) return 0;
  printf("\n\n");
#ifdef X11
  G2_DRAW_BACKBONE
  g2_flush(sh.g2);
  usleep(1e3*args->x_delay_ms);
#endif
  free(sh.total_brightness);
  free(sh.weighted_sum);
  return 1;
}
