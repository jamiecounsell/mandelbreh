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

// from main.cc
extern MandelBulbParams mandelBulb_params;

#pragma acc routine seq
extern int UnProject(int ix, int iy, const int viewport[4], const double matInvProjModel[16], double* obj);

#pragma acc routine seq
extern void rayMarch (const RenderParams &render_params, const vec3 &from, const vec3  &to, 
  double eps, pixelData &pix_data, const MandelBulbParams &bulb_params);
//void rayMarch(const RenderParams &render_params, const vec3 &from, const vec3  &direction, double eps,
//        vec3 &pd_hit, vec3 &pd_normal, bool pd_escaped, const MandelBulbParams &bulb_params);


#pragma acc routine seq
extern void getColour(const RenderParams &render_params, const vec3 &normal, const vec3 &hit, const bool escaped,
         //const int colourType, const float brightness
         const vec3 &from, const vec3  &direction, vec3 &result
         );
//extern vec3 getColour(const pixelData &pixData, const RenderParams &render_params,
//          const vec3 &from, const vec3  &direction);

void renderFractal(const CameraParams &camera_params, const RenderParams &renderer_params, 
       unsigned char* image, int frame)
{

  const int height = renderer_params.height;
  const int width  = renderer_params.width;


  #pragma acc enter data present_or_copyin(camera_params)
  #pragma acc enter data present_or_copyin(                        \
                          camera_params.camPos[:3],      \
                          camera_params.camTarget[:3],   \
                          camera_params.camUp[:3],      \
                          camera_params.matModelView[:16], \
                          camera_params.matProjection[:16], \
                          camera_params.matInvProjModel[:16], \
                          camera_params.viewport[:4] \
                          )  

  #pragma acc enter data present_or_copyin(renderer_params)
  #pragma acc enter data present_or_copyin(                  \
                          renderer_params.file_name[:80] \
                          )

  #pragma acc kernels copy(image[0:height*width*3]) present_or_copyin(mandelBulb_params, width, height)
  {

  
  #ifndef _OPENACC
  double time = getTime();
  #endif
  
  const double eps = pow(10.0, renderer_params.detail); 
  double farPoint[3];
  vec3 to, from;
  
  //from.SetDoublePoint(camera_params.camPos);
  SET_DOUBLE_POINT(from, camera_params.camPos);


  pixelData pix_data;
  /*
  pix_data.hit.x = 0.0;
  pix_data.hit.y = 0.0;
  pix_data.hit.z = 0.0;
  pix_data.normal.x = 0.0;
  pix_data.normal.y = 0.0;
  pix_data.normal.z = 0.0;
  pix_data.escaped = false;
  */

  vec3 color = {0.0, 50.5, 230.3};
  
  int i,j,k;
  for(j = 0; j < height; j++)
  {
      //for each column pixel in the row
    for(i = 0; i < width; i++)
    {


      // get point on the 'far' plane
      // since we render one frame only, we can use the more specialized method
      UnProject(i, j, camera_params.viewport, camera_params.matInvProjModel, farPoint);
      
      
      // to = farPoint - camera_params.camPos
      SUBTRACT_POINT(to, farPoint, camera_params.camPos);//SubtractDoubleDouble(farPoint,camera_params.camPos);
      
      NORMALIZE(to);
      
      //render the pixel
      rayMarch(renderer_params, from, to, eps, pix_data, mandelBulb_params);

      vec3 pd_norm, pd_hit;//= {pix_data.hit.x, pix_data.normal.x, pix_data.hit.z};
      //double orly = pix_data.normal.y;
      //bool lol = pix_data.escaped;
      VEC(pd_norm, pix_data.normal.x, pix_data.normal.y, pix_data.normal.z); 
      VEC(pd_hit, pix_data.hit.x, pix_data.hit.y, pix_data.hit.z); 

      //get the colour at this pixel
      getColour(renderer_params, pd_norm, pd_hit, pix_data.escaped, 
        //renderer_params.colourType, renderer_params.brightness 
        from, to, color
        );  

      
      //save colour into texture
      k = (j * width + i)*3;
      image[k+2] = (unsigned char)(color.x * 255);
      image[k+1] = (unsigned char)(color.y * 255);
      image[k]   = (unsigned char)(color.z * 255);
    
    }

    #ifndef _OPENACC
    printProgress((j+1)/(double)height,getTime()-time, frame);
    #endif
  }


}//end pragma
}
