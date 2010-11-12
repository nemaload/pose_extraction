#ifndef ARGBOILER_ARGS
#define ARGBOILER_ARGS

#define ARGBOILER(ARG) \
  ARG(ARG_FIL1,input_filename,NULL,NULL,"the input image (as raw data)",NULL) \
  ARG(ARG_INT1,input_width,"w","width","the width of each input image slice",-1) \
  ARG(ARG_INT1,input_height,"h","height","the height of each input image slice",-1) \
  ARG(ARG_INT0,output_width,NULL,"worm-width","the width of the output image",-1) \
  ARG(ARG_INT0,output_height,NULL,"worm-height","the height of the output image",-1) \
  ARG(ARG_INT0,output_extension,NULL,"worm-extension","how many pixels on either edge of the backbone to include (head and tail)",-1) \
  ARG(ARG_FIL0,output_filename,"o","output","the output image (as raw data)",NULL) \
  ARG(ARG_INT0,sd_sample_size,NULL,"sdss","the sample size for computing standard deviation",1000) \
  ARG(ARG_DBL0,thresh_sds,"t","thresh","the number of standard deviations above mean makes a pixel considered part of the worm",0.5) \
  ARG(ARG_INT0,mst_sample_size,NULL,"mstss","the sample size for making the MST",120) \
  ARG(ARG_INT0,refine_sample_size,NULL,"rfss","the sample size of E_image in refining the backbone",3000) \
  ARG(ARG_INT0,refine_refresh_size,NULL,"rfrs","the number of E_image samples to replace each iteration",10) \
  ARG(ARG_DBL0,refine_threshhold,NULL,"rfth","the average distance (in pixels) points must move less than to terminate refinement",1.55) \
  ARG(ARG_INT0,restart_iterations,NULL,"rfri","if the refinement doesn't converge after this many iterations, restart with a new MST",1000) \
  ARG(ARG_INT0,delta_history,NULL,"rfdh","the number of iterations the points must move very little in a row to count",150) \
  ARG(ARG_LIT0,spread_voronoi,"3","spread","adjust control points using the nearest pixels of neighboring control points as well as their own",0) \
  ARG(ARG_LIT0,use_brightness,"r","weight","weight E_image by pixel brightness instead of just threshholding",0) \
  ARG(ARG_LIT0,no_interpolate,"c","no-interpolate","don't interpolate input pixels when restacking",0) \
  ARG(ARG_DBL0,alpha,"a","alpha","weight of E_image; 1 in the original paper",0.15) \
  ARG(ARG_DBL0,beta,"b","beta","weight of E_length; 0.5 in the original paper",1.1) \
  ARG(ARG_DBL0,gamma,"g","gamma","weight of E_smoothness; 0.5 in the original paper",0.6) \
  ARG(ARG_DBL0,delta,"d","delta","'inertia' term (not in the original paper)",2.2) \
  ARG(ARG_DBL0,image_scale,"s","scale","factor to divide image width and height by in display",4.0) \
  ARG(ARG_LITN,precache,"p","precache","mmap the image and read all the pixels; performing no processing if specified twice.",0) \

#endif
