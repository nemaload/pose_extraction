#include "graphics.h"
#include "util.h"

static xcb_connection_t* X;
static xcb_window_t win;
static xcb_screen_t* screen;

void graphics_init(const image_t* image, const args_t* args) {
  X = xcb_connect(NULL,NULL);
  const xcb_setup_t* setup = xcb_get_setup(X);
  xcb_screen_iterator_t iter = xcb_setup_roots_iterator(setup);
  screen = iter.data;

  win = xcb_generate_id(X);
  xcb_create_window(X,
      XCB_COPY_FROM_PARENT,
      win,
      screen->root,
      0, 0,
      image->width/args->image_scale, image->height/args->image_scale,
      10,
      XCB_WINDOW_CLASS_INPUT_OUTPUT,
      screen->root_visual,
      0,  NULL);

  xcb_map_window(X,win);
  xcb_flush(X);
}



