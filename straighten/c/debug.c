#include "debug.h"
#include <stdlib.h>

/*
 * These functions are for a primitive sort of profiling.
 */
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
    printf("[Running total of time: %lfms/%lfms]\n" \
          "######################################" \
        "\n######################################\n",total_time_usr,total_time_rt);
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
    if(l>last_l) {
      printf("\033[%dE",l-last_l);
      last_l=l;
    } else if (last_l>l) {
      printf("\033[%dF",last_l-l);
      last_l=l;
    }
    printf("\033[1G[>");
    for(k=0;k<PROGRESS_BAR_WIDTH;k++) {
      printf(" ");
    }
    printf("]                    %10s (%d,%d)\033[s",desc,l,last_l);
    fflush(stdout);
  }
  if(i>=0) {
    if(l>last_l) {
      printf("\033[%dE",l-last_l);
      last_l=l;
    } else if (last_l>l) {
      printf("\033[%dF",last_l-l);
      last_l=l;
    }
    if ((n/PROGRESS_BAR_WIDTH) == 0 || (i-1)%(n/PROGRESS_BAR_WIDTH) == 0) {
      if(i/((double)n/PROGRESS_BAR_WIDTH)<=PROGRESS_BAR_WIDTH){
        printf("\033[%dG=>\033[%dG",2+(int)(int)((i-1)/((double)n/PROGRESS_BAR_WIDTH)),PROGRESS_BAR_WIDTH+4);
        printf("%9d/%9d %10s\033[0K", i, n, desc);
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
