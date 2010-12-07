#include <unistd.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifndef NO_PROGRESS_BARS
#define PROGRESS_BAR_WIDTH  120
#endif

double elapsed(struct timespec e_t, struct timespec s_t);
typedef enum {START, END} se_t;
void half_step_start(const char*);
void step_start(const char*);
void step_end(void);
void progress_lt(int, int, int, int, char*);
void progress_l(int, int, int, char*);
void progress_t(int, int, int, char*);
void progress(int, int, char*);
