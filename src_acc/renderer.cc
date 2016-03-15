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

#include "color.h"
#include "mandelbulb.h"
#include "camera.h"
#include "vector3d.h"
#include "3d.h"

#ifdef _OPENACC
#include <openacc.h>
#endif

extern double getTime();
extern void   printProgress(double perc, double time, int frame);

#pragma acc routine seq
extern double rayMarch (const RenderParams &render_params, const vec3 &from, const vec3  &to, double eps, pixelData &pix_data);

#pragma acc routine seq
extern vec3 getColour(const pixelData &pixData, const RenderParams &render_params,
          const vec3 &from, const vec3  &direction);

void renderFractal(const CameraParams &camera_params, const RenderParams &renderer_params, 
       unsigned char* image, int frame)
{
  
  //double time = getTime();

  const int height = renderer_params.height;
  const int width  = renderer_params.width;
  
  #pragma acc enter data copyin(camera_params)
  #pragma acc enter data copyin(                        \
                          camera_params.camPos[:3],      \
                          camera_params.camTarget[:3],   \
                          camera_params.camUp[:3],      \
                          camera_params.matModelView[:16], \
                          camera_params.matProjection[:16], \
                          camera_params.matInvProjModel[:16], \
                          camera_params.viewport[:4] \
                          )  

  #pragma acc enter data copyin(renderer_params)
  #pragma acc enter data copyin(                  \
                          renderer_params.file_name[:80] \
                          )

  #pragma acc kernels copy(image[0:height*width*3])
  {

  //double testtt = camera_params.matModelView[5];

  const double eps = pow(10.0, renderer_params.detail); 
  double farPoint[3];
  vec3 to, from;
  
  SET_DOUBLE_POINT(from, camera_params.camPos);
  

  // can you use nested structs as long as they are created within?
  pixelData pix_data;
  
  vec3 color;


  int i,j,k;
  for(j = 0; j < height; j++)
  {
      //for each column pixel in the row
    for(i = 0; i <width; i++)
    {
      // get point on the 'far' plane
      // since we render one frame only, we can use the more specialized method
      //UnProject(i, j, camera_params, farPoint);
      UnProject(i, j, camera_params.viewport, camera_params.matInvProjModel, farPoint);

      
      SUBTRACT_POINT(to, farPoint, camera_params.camPos);//SubtractDoubleDouble(farPoint,camera_params.camPos);
      NORMALIZE(to);
      
      //render the pixel
      rayMarch(renderer_params, from, to, eps, pix_data);
      
      //get the colour at this pixel
      // wrong # of arguments error
      color = getColour(pix_data, renderer_params, from, to);
        
      //save colour into texture
      k = (j * width + i)*3;
      image[k+2] = (unsigned char)(color.x * 255);
      image[k+1] = (unsigned char)(color.y * 255);
      image[k]   = (unsigned char)(color.z * 255);
    }
    //printProgress((j+1)/(double)height,getTime()-time, frame);
  }

  }// end pragma
}
