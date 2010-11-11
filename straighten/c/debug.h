#include <unistd.h>
#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define PROGRESS_BAR_WIDTH  120

double elapsed(struct timespec e_t, struct timespec s_t);
typedef enum {START, END} se_t;
void half_step_start(const char*);
void step_start(const char*);
void step_end(void);
void progress(int, int, int, char*);