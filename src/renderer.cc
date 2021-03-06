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
#include <stdio.h>
#include <stdlib.h>

#include "color.h"
#include "mandelbulb.h"
#include "mandelbox.h"
#include "camera.h"
#include "vector3d.h"
#include "3d.h"

#ifdef _OPENACC
#include <openacc.h>
#endif

// All parallelization is done in renderFractal
// raymarch, getcolor, and unproject are sequential openacc routines
// vec3 has been converted from class to struct with macros (see vector3d.h)
// The cost for distinct Mandelbulb and Mandelbox params is duplicate method signatures
// Specified with -DBULB or -DBOX compiler flag

extern double getTime();
extern void   printProgress( double perc, double time, int frame );

#ifdef BULB
#pragma acc routine seq
extern double rayMarch(const int maxRaySteps, const float maxDistance,
  const float escape_time, const float power, const int num_iter,
  const vec3 &from, const vec3 &direction, double eps, pixelData& pix_data);
#else 
#pragma acc routine seq
extern double rayMarch(const int maxRaySteps, const float maxDistance,
  const int num_iter, const float rMin, const float rFixed, const float escape_time, const float scale,
  const vec3 &from, const vec3 &direction, double eps, pixelData& pix_data);
#endif

#pragma acc routine seq
extern void getcolor(const pixelData &pixData, const int colorType, const float brightness,
		      const vec3 &from, const vec3  &direction, vec3 &result);

#pragma acc routine seq
extern void UnProject(double winX, double winY, const int viewport[4], const double matInvProjModel[16], vec3 &obj);


#ifdef BULB
void renderFractal(const CameraParams camera_params, const RenderParams renderer_params, 
                    const MandelBulbParams bulb_params, unsigned char* image, int frame)
#else
void renderFractal(const CameraParams camera_params, const RenderParams renderer_params, 
                    const MandelBoxParams box_params, unsigned char* image, int frame)
#endif
{
  // DIRECTION, COLOR, PIXEL ARRAYS
  // OpenACC has problems with incorrectly sharing structs
  // Solved by creating struct array (one per pixel, or loop iteration) that is shared in parallel region 
  int size = renderer_params.width * renderer_params.height;
#ifdef _OPENACC
    vec3* direction = (vec3*)acc_malloc(size * sizeof(vec3));
    pixelData* pixel = (pixelData*)acc_malloc(size * sizeof(pixelData));
    vec3* color = (vec3*)acc_malloc(size * sizeof(vec3));
#else
    vec3* direction = (vec3*)malloc(size * sizeof(vec3));
    pixelData* pixel = (pixelData*)malloc(size * sizeof(pixelData));
    vec3* color = (vec3*)malloc(size * sizeof(vec3));
#endif

  // All parameters are explicitly copied into ACC device region
  // parameter structs are no longer passed as function arguments

  // RENDERER PARAMS
  const int colorType = renderer_params.colorType;
  const float brightness = renderer_params.brightness;
  const int height = renderer_params.height;
  const int width  = renderer_params.width;
  const float detail = renderer_params.detail;
  const int maxRaySteps = renderer_params.maxRaySteps;
  const float maxDistance = renderer_params.maxDistance;

  // CAMERA PARAMS
  const double camPos[3] = {camera_params.camPos[0], camera_params.camPos[1], camera_params.camPos[2]} ;
  const double matInvProjModel[16] =
  {
    camera_params.matInvProjModel[0],
    camera_params.matInvProjModel[1],
    camera_params.matInvProjModel[2],
    camera_params.matInvProjModel[3],
    camera_params.matInvProjModel[4],
    camera_params.matInvProjModel[5],
    camera_params.matInvProjModel[6],
    camera_params.matInvProjModel[7],
    camera_params.matInvProjModel[8],
    camera_params.matInvProjModel[9],
    camera_params.matInvProjModel[10],
    camera_params.matInvProjModel[11],
    camera_params.matInvProjModel[12],
    camera_params.matInvProjModel[13],
    camera_params.matInvProjModel[14],
    camera_params.matInvProjModel[15]
  };
  const int viewport[4] =
  {
    camera_params.viewport[0],
    camera_params.viewport[1],
    camera_params.viewport[2],
    camera_params.viewport[3]
  };

  printf("(%lf, %lf, %lf)\n",camera_params.camPos[0], camera_params.camPos[1], camera_params.camPos[2]);

  // copy bulb/box params into device region

  #ifdef BULB
    // MANDELBULB PARAMS
    const  float escape_time =   bulb_params.escape_time;
    const float power =  bulb_params.power;
    const int num_iter =  bulb_params.num_iter; 

    #pragma acc enter data pcopyin(    \
      escape_time, \
      power, \
      num_iter \
    )
  #else
    
    // MANDELBOX PARAMS
    const int num_iter = box_params.num_iter;
    const float rMin = box_params.rMin;
    const float rFixed = box_params.rFixed;
    const float escape_time = box_params.escape_time;
    const float scale = box_params.scale;

    #pragma acc enter data pcopyin(    \
      rMin, \
      rFixed, \
      escape_time, \
      scale, \
      num_iter \
    )

  #endif

  // Copy in image, remaining parameters, and vec3 arrays
  #pragma acc data present_or_copy(image[0:size*3]),   \
  pcopyin(                 \
    camPos[:3],      \
    matInvProjModel[:16], \
    viewport[:4], \
                  \
    colorType, \
    brightness, \
    height, \
    width, \
    detail, \
    maxRaySteps, \
    maxDistance \
  ),            \
    deviceptr(direction, pixel, color)
  {

  // BEGIN DEVICE DATA REGION
    
  #ifndef _OPENACC
  double time = getTime();
  #endif

  const double eps = pow(10.0, detail); 
  const vec3 from = {camPos[0], camPos[1], camPos[2]};

  // for some reason needed for compiler to parallelize loops
  const int cheight = height;
  const int cwidth = width;

  int i,j;
  // total of three parallel pragmas + external routine pragmas
  #pragma acc parallel 
  #pragma acc loop
  for(j = 0; j < cheight; j++)
  {
    #pragma acc loop
    for(i = 0; i < cwidth; i++)
    {
      
      int k, l;
      //vec3 array index
      l = (j * width + i );

  	  // get point on the 'far' plane
      // acc routine
      UnProject(i, j, viewport, matInvProjModel, direction[l]);
  	  
      SUBTRACT_DOUBLE_ARRAY(direction[l], camPos);
      NORMALIZE( direction[l] );	  
  	  
      //render the pixel
      //acc routine
      //difference in box/bulb is the distance estimator used by raymarch
      #ifdef BULB
      rayMarch(maxRaySteps, maxDistance, escape_time, power, num_iter, from, direction[l], eps, pixel[l]); 
      #else
      rayMarch(maxRaySteps, maxDistance, num_iter, rMin, 
        rFixed, escape_time, scale, from, direction[l], eps, pixel[l]);  
      #endif

  	  //get the color at this pixel
      getcolor(pixel[l], colorType, brightness, from, direction[l], color[l]);
        
  	  //save color into texture
  	  k = (j * width + i)*3;
  	  image[k+2] = (unsigned char)(color[l].x * 255);
  	  image[k+1] = (unsigned char)(color[l].y * 255);
  	  image[k]   = (unsigned char)(color[l].z * 255);

    }

      #ifndef _OPENACC
      printProgress((j+1)/(double)height,getTime()-time, frame);
      #endif
  }

  }// END DEVICE DATA REGION

  #ifdef _OPENACC
    // free memory used by acc for vec3 arrays
    acc_free(direction);
    acc_free(pixel);
    acc_free(color);
  #endif
}
