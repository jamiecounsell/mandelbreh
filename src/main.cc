/*
   This file is part of the Mandelbox program developed for the course
    CS/SE  Distributed Computer Systems taught by N. Nedialkov in the
    Winter of 2015-2016 at McMaster University.

    Copyright (C) 2015-2016 T. Gwosdz and N. Nedialkov

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "camera.h"
#include "renderer.h"
#include "mandelbulb.h"
#include "mandelbox.h"
#include "walk.h"
#include <unistd.h>

#ifdef BULB

void getParameters(char *filename, CameraParams *camera_params, RenderParams *renderer_params,
		   MandelBulbParams *mandelBulb_paramsP);
void renderFractal(const CameraParams camera_params, const RenderParams renderer_params, 
                  const MandelBulbParams bulb_params, unsigned char* image, int frame);
MandelBulbParams mandelBulb_params;

#else

void getParameters(char *filename, CameraParams *camera_params, RenderParams *renderer_params,
       MandelBoxParams *mandelBox_paramsP);
void renderFractal(const CameraParams camera_params, const RenderParams renderer_params, 
                  const MandelBoxParams box_params, unsigned char* image, int frame);
MandelBoxParams mandelBox_params;

#endif

void init3D       (CameraParams *camera_params, const RenderParams *renderer_params);
void saveBMP      (const char* filename, const unsigned char* image, int width, int height);

int main(int argc, char** argv)
{
  // Get parameters:
  int vflag = 0,
      verbose = 1,
      num_of_iterations = 1,
      c;
  opterr = 0;
  char * fname = argv[1];

  while ((c = getopt (argc, argv, "hvnf:")) != -1)
    switch (c)
    {
      case 'h':
        printf("Usage:");
        #ifdef BULB
          printf(" ./mandelbox ");
        #else
          printf(" ./mandelbulb ");
        #endif
        printf("-hvnf\n");
        printf("   -h   : Display this help message\n");
        printf("   -v   : Generate video from frames\n");
        printf("   -n   : Less verbose output\n");
        printf("   -f n : Generate n frames\n");
        return 0;
      case 'v':
        vflag = 1;
        break;
      case 'n':
        verbose = 0;
        break;
      case 'f':
        if (optarg == NULL){
          printf("Invalid option -f\nInteger required. Ex: -n 100\n");
          return 1;
        }
        num_of_iterations = atoi(optarg);
        break;
      default:
        printf("Unknwon Option. Aborting.\n");
        return 1;
  }

  int frame, i;
  RenderParams renderer_params;
  CameraParams camera_history [num_of_iterations + 1];

  // Get bulb/box params
  for (i=0;i<num_of_iterations;i++){
    #ifdef BULB  
      getParameters(fname, &camera_history[i], &renderer_params, &mandelBulb_params);
    #else
      getParameters(fname, &camera_history[i], &renderer_params, &mandelBox_params);
    #endif
    init3D(&camera_history[i], &renderer_params);
  }
    // Initialize params and image
  int image_size = renderer_params.width * renderer_params.height;
  unsigned char *image = (unsigned char*)malloc(3*image_size*sizeof(unsigned char));

  // Verbose output. Silence with -n
  if (verbose) {
    printf("Image Size:       %dx%d\n", renderer_params.width, renderer_params.height);
    printf("Video:            %s\n", vflag ? "Yes" : "No");
      printf("Number of Frames: %d", num_of_iterations);
    if (num_of_iterations == 1) {
      printf(" (for more frames, use -f)\n\n");
    } else { printf("\n\n"); }
  }

  for (frame = 0; frame < num_of_iterations; frame++){

    // Generate unique image name
    char buf[15];
    sprintf(buf, "../frames/%05d.bmp", frame);

    if (verbose) {
      printf("Rendering frame: %d\n", frame);
    }

    #ifdef BULB
      // Mandelbulb
      walk(camera_history, &renderer_params, &mandelBulb_params, verbose, frame);
      renderFractal(camera_history[frame], renderer_params, mandelBulb_params, image, frame);
    #else
      // Mandelbox
      walk(camera_history, &renderer_params, &mandelBox_params, verbose, frame);
      renderFractal(camera_history[frame], renderer_params, mandelBox_params, image, frame);
    #endif

    // Save image
    saveBMP(buf, image, renderer_params.width, renderer_params.height);
  }
 
  // Cleanup
  free(image);

  // Video shell script
  if (vflag){
    system("./genvideo.sh");
  }

  return 0;
}
