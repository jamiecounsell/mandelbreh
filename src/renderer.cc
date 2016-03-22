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
extern void   printProgress( double perc, double time, int frame );

#pragma acc routine seq
extern void rayMarch(const int maxRaySteps, const float maxDistance,
 const float escape_time, const float power, const int num_iter,
 const vec3 &from, const vec3  &direction, double eps, pixelData& pix_data);
//extern void rayMarch (const RenderParams render_params, const MandelBulbParams bulb_params,
//          const vec3 &from, const vec3  &to, double eps, pixelData &pix_data);

#pragma acc routine seq
extern void getColour(const pixelData &pixData, const int colourType, const float brightness,//const RenderParams render_params,
		      const vec3 &from, const vec3  &direction, vec3 &result);

#pragma acc routine seq
extern void UnProject(double winX, double winY, const int viewport[4], const double matInvProjModel[16], vec3 &obj);
//extern void UnProject(double winX, double winY, CameraParams camP, vec3 &obj);//double *obj);

void renderFractal(const CameraParams camera_params, const RenderParams renderer_params, 
                    const MandelBulbParams bulb_params, unsigned char* image, int frame)
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

const int fractalType = renderer_params.fractalType;
  const int colourType = renderer_params.colourType;
  const float brightness = renderer_params.brightness;
  const int height = renderer_params.height;
  const int width  = renderer_params.width;
  const float detail = renderer_params.detail;
  const int maxRaySteps = renderer_params.maxRaySteps;
  const float maxDistance = renderer_params.maxDistance;

  const double camPos[3] = {camera_params.camPos[0], camera_params.camPos[1], camera_params.camPos[2]} ;
  const double camTarget[3] = {camera_params.camTarget[0],camera_params.camTarget[1],camera_params.camTarget[2]};
  const double camUp[3] = { camera_params.camUp[0], camera_params.camUp[1], camera_params.camUp[2] };
  const double fov = camera_params.fov;
  const double matModelView[16] = 
  {
    camera_params.matModelView[0],
    camera_params.matModelView[1],
    camera_params.matModelView[2],
    camera_params.matModelView[3],
    camera_params.matModelView[4],
    camera_params.matModelView[5],
    camera_params.matModelView[6],
    camera_params.matModelView[7],
    camera_params.matModelView[8],
    camera_params.matModelView[9],
    camera_params.matModelView[10],
    camera_params.matModelView[11],
    camera_params.matModelView[12],
    camera_params.matModelView[13],
    camera_params.matModelView[14],
    camera_params.matModelView[15]
  };
  const double matProjection[16] =
  {
    camera_params.matProjection[0],
    camera_params.matProjection[1],
    camera_params.matProjection[2],
    camera_params.matProjection[3],
    camera_params.matProjection[4],
    camera_params.matProjection[5],
    camera_params.matProjection[6],
    camera_params.matProjection[7],
    camera_params.matProjection[8],
    camera_params.matProjection[9],
    camera_params.matProjection[10],
    camera_params.matProjection[11],
    camera_params.matProjection[12],
    camera_params.matProjection[13],
    camera_params.matProjection[14],
    camera_params.matProjection[15]
  };
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

  const  float escape_time =   bulb_params.escape_time;
  const float power =  bulb_params.power;
  const int num_iter =  bulb_params.num_iter; 


 /*
  #pragma acc data copy(image[0:size*3]),   \
  pcopyin(                 \
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
*/

#pragma acc data copy(image[0:size*3]),   \
  pcopyin(                 \
    camPos[:3],      \
    camTarget[:3],   \
    camUp[:3],      \
    fov,            \
    matModelView[:16], \
    matProjection[:16], \
    matInvProjModel[:16], \
    viewport[:4], \
                  \
    fractalType, \
    colourType, \
    brightness, \
    height, \
    width, \
    detail, \
    maxRaySteps, \
    maxDistance, \
                  \
    escape_time, \
    power, \
    num_iter \
  ),            \
    deviceptr(to_arr, pix_arr, color_arr)
  {

  const double eps = pow(10.0, detail); 
  const vec3 from = {camPos[0], camPos[1], camPos[2]};
  

    
  #ifndef _OPENACC
  double time = getTime();
  #endif

  // for some reason needed for compiler to parallelize loops
  const int cheight = height;
  const int cwidth = width;
  
  int i,j;
  #pragma acc parallel 
  #pragma acc loop
  for(j = 0; j < cheight; j++)
    {
      #pragma acc loop
      for(i = 0; i < cwidth; i++)
	    {
        
        int k, l;
        l = (j * width + i );
 
    	  // get point on the 'far' plane
        UnProject(i, j, viewport, matInvProjModel, to_arr[l]);//farPoint);
    	  
        to_arr[l].x = to_arr[l].x - camera_params.camPos[0];
        to_arr[l].y = to_arr[l].y - camera_params.camPos[1];
        to_arr[l].z = to_arr[l].z - camera_params.camPos[2];


        NORMALIZE( to_arr[l] );	  
    	  
        //render the pixel
        rayMarch(maxRaySteps, maxDistance, escape_time, power, num_iter, from, to_arr[l], eps, pix_arr[l]);     

    	  //get the colour at this pixel
        getColour(pix_arr[l], colourType, brightness, from, to_arr[l], color_arr[l]);
          
    	  //save colour into texture
    	  k = (j * width + i)*3;
    	  image[k+2] = (unsigned char)(color_arr[l].x * 255);
    	  image[k+1] = (unsigned char)(color_arr[l].y * 255);
    	  image[k]   = (unsigned char)(color_arr[l].z * 255);

	    }

      #ifndef _OPENACC
      printProgress((j+1)/(double)height,getTime()-time, frame);
      #endif
      
    }

  }// end device data region

  printf("\n rendering done:\n");
}
