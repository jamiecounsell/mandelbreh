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

void getParameters(char *filename, CameraParams *camera_params, RenderParams *renderer_params,
		   MandelBulbParams *mandelBulb_paramsP);
void init3D       (CameraParams *camera_params, const RenderParams *renderer_params);
void renderFractal(const CameraParams &camera_params, const RenderParams &renderer_params, unsigned char* image, int frame);
void saveBMP      (const char* filename, const unsigned char* image, int width, int height);

MandelBulbParams mandelBulb_params;

int main(int argc, char** argv)
{
  int i;
  CameraParams    camera_params;
  RenderParams    renderer_params;
  
  getParameters(argv[1], &camera_params, &renderer_params, &mandelBulb_params);

  int image_size = renderer_params.width * renderer_params.height;
  unsigned char *image = (unsigned char*)malloc(3*image_size*sizeof(unsigned char));
  init3D(&camera_params, &renderer_params);
  camera_params.camPos[0] = camera_params.camPos[0];
  camera_params.camPos[1] = camera_params.camPos[1];
  camera_params.camPos[2] = camera_params.camPos[2];
  for (i = 0; i < 1; i++){
    char buf[15];
    printf("Computing frame %d...\n", i);
    sprintf(buf, "../frames/%05d.bmp", i);

    camera_params.camPos[0] = camera_params.camPos[0]-(0.01);
    camera_params.camPos[1] = camera_params.camPos[1]-(0.01);
    camera_params.camPos[2] = camera_params.camPos[2]-(0.01);
    renderFractal(camera_params, renderer_params, image, i);
    saveBMP(buf, image, renderer_params.width, renderer_params.height);
  }
  printf("\n");
  free(image);
  #ifdef VIDEO
    system("./genvideo.sh");
  #endif

  return 0;
}
