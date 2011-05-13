#include <stdint.h>
#include <fenv.h>
#include <SDL/SDL.h>
#include "../common/util.h"
#include "../common/debug.h"

#define ARGBOILER(ARG) \
  ARG(ARG_FIL1,input_filename,NULL,NULL,"the input image (as raw data)","") \
  ARG(ARG_INT1,input_width,"w","width","the width of each input image slice",-1) \
  ARG(ARG_INT1,input_height,"h","height","the height of each input image slice",-1) \
  ARG(ARG_DBL0,image_scale,"s","scale","factor to divide image width and height by in display",4.0) \
  ARG(ARG_INT0,sensitivity,"b","sensitivity","number from 0-255 to count each sample for",50) \

#include "../common/argboiler.h"

void init_video(void) {
	if(SDL_Init(SDL_INIT_VIDEO) < 0) {
		fprintf(stderr, "Unable to init SDL: %s\n", SDL_GetError());
		exit(1);
	}
	atexit(SDL_Quit);
}

/*
 * Bundle all the little initializations together
 */
void init(image_t* image, const args_t* args) {
  init_rng(image);
  fesetround(FE_DOWNWARD);
  image->width = args->input_width;
  image->height = args->input_height;
  compute_depth(image);
  init_video();
}

int main(int argc, char** argv) {
  args_t args;
  image_t input;
  int length;
  int i;
  SDL_Surface *screen;
  SDL_Surface *canvas;

  printf("Parsing command line...\n");
  parse_args(argc, argv, &args);

  printf("Sampling visualizer\nDavid Dalrymple\n==================\n");

  step_start("mmapping file");
  input.data = open_mmapped_file_read(args.input_filename, &length);
  input.length = length/2;

  if(input.data==NULL) {
    printf("Couldn't open file\n");
    return 1;
  }
  step_end();

  init(&input, &args);
  screen = SDL_SetVideoMode(((double)input.width)/args.image_scale+1,((double)input.height)/args.image_scale+2, 32, SDL_HWSURFACE);
  i=0;
  while(1) {
    printf("Samples: %d\n", i++);
    point_t* sample = perform_sample(&input, 1);
    int sx = sample->p[2]/args.image_scale;
    int sy = screen->h - (sample->p[1]/args.image_scale) - 1;
    uint8_t r, g, b;
    uint32_t* pixel = &(((uint32_t*)screen->pixels)[sx+sy*screen->w]);
    SDL_GetRGB(*pixel, screen->format, &r, &g, &b);
    uint8_t update(uint8_t v) {
      int res = (v+5)*2;
      if(res > 255) return 255;
      return (uint8_t)res;
    }
    *pixel= SDL_MapRGB(screen->format, update(r), update(g), update(b));
    SDL_UpdateRect(screen,sx,sy,1,1);
    free(sample);
  }

  return 0;
}
