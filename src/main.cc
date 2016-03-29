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
#include "camera.h"
#include "renderer.h"
#include "mandelbulb.h"
#include "mandelbox.h"
#include "walk.h"

void walk(CameraParams *camera_params);

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
  int i;
  CameraParams    camera_params;
  RenderParams    renderer_params;

#ifdef BULB  
  getParameters(argv[1], &camera_params, &renderer_params, &mandelBulb_params);
#else
  getParameters(argv[1], &camera_params, &renderer_params, &mandelBox_params);
#endif

  int num_of_iterations = 1;

  if (argc > 2) {
    num_of_iterations = atoi(argv[2]);
    printf("Running %d iterations.\n", num_of_iterations);
  }

  int image_size = renderer_params.width * renderer_params.height;
  unsigned char *image = (unsigned char*)malloc(3*image_size*sizeof(unsigned char));

  init3D(&camera_params, &renderer_params);

  for (i = 0; i < num_of_iterations; i++){
      char buf[15];

      sprintf(buf, "../frames/%05d.bmp", i);

  #ifdef BULB
      walk(&camera_params, &renderer_params, &mandelBulb_params);
      renderFractal(camera_params, renderer_params, mandelBulb_params, image, i);
  #else
      walk(&camera_params, &renderer_params, &mandelBox_params);
      renderFractal(camera_params, renderer_params, mandelBox_params, image, i);
  #endif
      
      saveBMP(buf, image, renderer_params.width, renderer_params.height);
  }
  free(image);
  #ifdef VIDEO
    system("./genvideo.sh");
  #endif

  return 0;
}
