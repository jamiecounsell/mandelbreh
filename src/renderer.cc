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
#include "camera.h"
#include "vector3d.h"
#include "3d.h"

#ifdef _OPENACC
#include <openacc.h>
#endif

extern double getTime();
extern void   printProgress( double perc, double time );

#pragma acc routine seq
extern void rayMarch (const RenderParams render_params, const MandelBulbParams bulb_params,
          const vec3 &from, const vec3  &to, double eps, pixelData &pix_data);

#pragma acc routine seq
extern void getColour(const pixelData &pixData, const RenderParams render_params,
		      const vec3 &from, const vec3  &direction, vec3 &result);

#pragma acc routine seq
extern void UnProject(double winX, double winY, CameraParams camP, double *obj);

void renderFractal(const CameraParams camera_params, const RenderParams renderer_params, 
                    const MandelBulbParams bulb_params, unsigned char* image)
{

  int size = renderer_params.width * renderer_params.height;

  #ifdef _OPENACC
    vec3* to_arr = (vec3*)acc_malloc(size * sizeof(vec3));
    pixelData* pix_arr = (pixelData*)acc_malloc(size * sizeof(pixelData));
    vec3* color_arr = (vec3*)acc_malloc(size * sizeof(vec3));
  #else
    vec3* to_arr = (vec3*)malloc(size * sizeof(vec3));
    pixelData* pix_arr = (pixelData*)malloc(size * sizeof(pixelData));
    vec3* color_arr = (vec3*)malloc(size * sizeof(vec3));
  #endif

 
  #pragma acc data copy(image[0:size*3]),   \
  copyin(                 \
    camera_params, \
    camera_params.camPos[:3],      \
    camera_params.camTarget[:3],   \
    camera_params.camUp[:3],      \
    camera_params.fov,            \
    camera_params.matModelView[:16], \
    camera_params.matProjection[:16], \
    camera_params.matInvProjModel[:16], \
    camera_params.viewport[:4], \
                  \
    renderer_params, \
    renderer_params.fractalType, \
    renderer_params.colourType, \
    renderer_params.brightness, \
    renderer_params.height, \
    renderer_params.width, \
    renderer_params.detail, \
    renderer_params.maxRaySteps, \
    renderer_params.maxDistance, \
                  \
    bulb_params, \
    bulb_params.escape_time, \
    bulb_params.power, \
    bulb_params.num_iter \
  ),            \
    deviceptr(to_arr, pix_arr, color_arr)
  {

  const double eps = pow(10.0, renderer_params.detail); 
  const vec3 from = {camera_params.camPos[0], camera_params.camPos[1], camera_params.camPos[2]};
  double farPoint[3];
  
  //const int height = renderer_params.height;
  //const int width  = renderer_params.width;
    
  #ifndef _OPENACC
  double time = getTime();
  #endif
  
  int i,j;
  #pragma acc parallel //loop collapse(2)
  #pragma acc loop
  for(j = 0; j < renderer_params.height; j++)
    {
      #pragma acc loop
      //for each column pixel in the row
      for(i = 0; i < renderer_params.width; i++)
	    {
        

        int k, l;
        l = (j * renderer_params.width + i );
 
    	  // get point on the 'far' plane
        UnProject(i, j, camera_params, farPoint);
    	  
        SUBTRACT_POINT(to_arr[l], farPoint, camera_params.camPos);
        NORMALIZE( to_arr[l] );	  
    	  
        //render the pixel
    	  
        rayMarch(renderer_params, bulb_params, from, to_arr[l], eps, pix_arr[l]);  	  
      
    	  //get the colour at this pixel
    	  getColour(pix_arr[l], renderer_params, from, to_arr[l], color_arr[l]);
          
    	  //save colour into texture
    	  k = (j * width + i)*3;
    	  image[k+2] = (unsigned char)(color_arr[l].x * 255);
    	  image[k+1] = (unsigned char)(color_arr[l].y * 255);
    	  image[k]   = (unsigned char)(color_arr[l].z * 255);

	   }

      #ifndef _OPENACC
      printProgress((j+1)/(double)height,getTime()-time);
      #endif
      
    }

  }// end device region

  printf("\n rendering done:\n");
}
