#include <stdio.h>

#include "color.h"
#include "mandelbulb.h"
#include "camera.h"
#include "vector3d.h"
#include "3d.h"

#ifdef _OPENACC
#include <openacc.h>
#endif

//extern double getTime();
//extern void   printProgress(double perc, double time, int frame);

// from main.cc
extern MandelBulbParams mandelBulb_params;

#pragma acc routine seq
extern double DE(const vec3 &p, const MandelBulbParams &params);

#pragma acc routine seq
extern double rayMarch (const int maxRaySteps, const float maxDistance, const vec3 &from, 
  const vec3  &to, double eps, pixelData2 &pix_data, const MandelBulbParams &bulb_params);
//extern double rayMarch (const RenderParams &render_params, const vec3 &from, const vec3  &to, double eps, pixelData &pix_data);

#pragma acc routine seq
extern void getColour(const pixelData2 &pixData, const int colourType, const float brightness,
          const vec3 &from, const vec3  &direction, vec3 &result);
//extern vec3 getColour(const pixelData &pixData, const RenderParams &render_params,
//          const vec3 &from, const vec3  &direction);


void renderFractal(const CameraParams &camera_params, const RenderParams &renderer_params, 
       unsigned char* image, int frame)
{
  
  //double time = getTime();

  const int height = renderer_params.height;
  const int width  = renderer_params.width;

  pixelData2 pix;
  
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

  // pow function?
  const double eps = pow(10.0, renderer_params.detail); 
  double farPoint[3];
  vec3 to, from;
  
  SET_DOUBLE_POINT(from, camera_params.camPos);
  


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

      SUBTRACT_POINT(to, farPoint, camera_params.camPos);
      NORMALIZE(to);
      
      //render the pixel
      const int max_rs =  renderer_params.maxRaySteps;
      const float max_d = renderer_params.maxDistance;
      rayMarch(max_rs, max_d, from, to, eps, pix, mandelBulb_params);
      
      const int ct = renderer_params.colourType;
      const float br = renderer_params.brightness;
      getColour(pix, ct, br, from, to, color);
        
      //save colour into texture
      k = (j * width + i)*3;
      image[k+2] = (unsigned char)(color.x * 255);
      image[k+1] = (unsigned char)(color.y * 255);
      image[k]   = (unsigned char)(color.z * 255);
    }
    //printProgress((j+1)/(double)height,getTime()-time, frame);
  }

  }// end pragma

  // pragma acc exit



}
