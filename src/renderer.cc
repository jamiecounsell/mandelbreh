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
#include <stdlib.h>

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
extern double MandelBulbDistanceEstimator(const vec3 &p, const MandelBulbParams &params);



//#pragma acc routine seq
//    extern void rayMarch (const RenderParams &render_params, const vec3 &from, const vec3  &to, 
//      double eps, pixelData &pix_data, const MandelBulbParams &bulb_params);
//void rayMarch(const RenderParams &render_params, const vec3 &from, const vec3  &direction, double eps,
//        vec3 &pd_hit, vec3 &pd_normal, bool pd_escaped, const MandelBulbParams &bulb_params);




//#pragma acc routine seq
//extern void getColour(const RenderParams &render_params, const vec3 &normal, const vec3 &hit, const bool escaped,
         //const int colourType, const float brightness
//         const vec3 &from, const vec3  &direction, vec3* color_ptr
//         );
//extern vec3 getColour(const pixelData &pixData, const RenderParams &render_params,
//          const vec3 &from, const vec3  &direction);

    void renderFractal(const CameraParams &camera_params, const RenderParams &renderer_params, 
     unsigned char* image, int frame)
    {
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


  #ifdef _OPENACC
      vec3* color_ptr = (vec3*)acc_malloc(sizeof(vec3));
  #else
      vec3* color_ptr = (vec3*)malloc(sizeof(vec3));
  #endif

  #pragma acc data present_or_copyin(                        \
      fov,            \
      camPos[:3],      \
      camTarget[:3],   \
      camUp[:3],      \
      matModelView[:16], \
      matProjection[:16], \
      matInvProjModel[:16], \
      viewport[:4] \
      )  

  #pragma acc data present_or_copyin(                  \
    fractalType, \
    colourType, \
    brightness, \
    height, \
    width, \
    detail, \
    maxRaySteps, \
    maxDistance \
  )

  #pragma acc kernels copy(image[0:height*width*3])           \
  present_or_copyin(mandelBulb_params) \
  deviceptr(color_ptr)
  {


  #ifndef _OPENACC
    double time = getTime();
  #endif

    const double eps = pow(10.0, detail); 
    double farPoint[3];
    vec3 to, from;

  //from.SetDoublePoint(camera_params.camPos);
    SET_DOUBLE_POINT(from, camPos);


    pixelData pix_data;
  //vec3 color;

    int i,j,k;
    for(j = 0; j < height; j++)
    {
      //for each column pixel in the row
      for(i = 0; i < width; i++)
      {


      // get point on the 'far' plane
      // since we render one frame only, we can use the more specialized method
        UnProject(i, j, viewport, matInvProjModel, farPoint);


      // to = farPoint - camera_params.camPos
      SUBTRACT_POINT(to, farPoint, camPos);//SubtractDoubleDouble(farPoint,camera_params.camPos);
      NORMALIZE(to);
      
/* INLINE RAYMARCH */
//      rayMarch(renderer_params, from, to, eps, pix_data, mandelBulb_params);
//    extern void rayMarch (const RenderParams &render_params, const vec3 &from, const vec3  &to, 
//      double eps, pixelData &pix_data, const MandelBulbParams &bulb_params);



  double dist = 0.0;
  double totalDist = 0.0;
  
  // We will adjust the minimum distance based on the current zoom
  double epsModified = 0.0;
  
  
  int steps=0;
  vec3 p;
  do 
    {      
      //p = from + direction * totalDist;
      VEC(p,
        from.x + to.x * totalDist, 
        from.y + to.y * totalDist, 
        from.z + to.z * totalDist
      );

      dist = MandelBulbDistanceEstimator(p, mandelBulb_params);
      
      totalDist += .95*dist;
      
      epsModified = totalDist;
      epsModified*=eps;
      steps++;
    }
  while (dist > epsModified && totalDist <= maxDistance && steps < maxRaySteps);
 
  

  if (dist < epsModified) 
    {
      //we didnt escape
      pix_data.escaped = false;
      
      // We hit something, or reached MaxRaySteps
      pix_data.hit = p;
      
      //figure out the normal of the surface at this point
      const vec3 normPos = {
        p.x - to.x * epsModified, 
        p.y - to.y * epsModified, 
        p.z - to.z * epsModified
      };
     

/* INLINE NORMAL */
//      normal(normPos, pix_data.normal, bulb_params);
// inline void normal(const vec3 & p, vec3 & normal, const MandelBulbParams &bulb_params)


 // compute the normal at p
  const double sqrt_mach_eps = 1.4901e-08;

  double eps ; //= std::max( p.Magnitude(), 1.0 )*sqrt_mach_eps;
  MAGNITUDE(eps, normPos);
  eps = MAX(eps, 1.0);
  eps *= sqrt_mach_eps;

  vec3 e1 = {eps, 0, 0}; 
  vec3 e2 = {0, eps, 0}; 
  vec3 e3 = {0, 0, eps}; 

  vec3 vs1, vs2, vs3;
  vec3 vd1, vd2, vd3;
  VECTOR_SUM(vs1, normPos,e1);
  VECTOR_SUM(vs2, normPos,e2);
  VECTOR_SUM(vs3, normPos,e3);

  VECTOR_DIFF(vd1, normPos, e1);
  VECTOR_DIFF(vd2, normPos, e2);
  VECTOR_DIFF(vd3, normPos, e3);

  double nx = MandelBulbDistanceEstimator(vs1, mandelBulb_params) - MandelBulbDistanceEstimator(vd1, mandelBulb_params);
  double ny = MandelBulbDistanceEstimator(vs2, mandelBulb_params) - MandelBulbDistanceEstimator(vd2, mandelBulb_params);
  double nz = MandelBulbDistanceEstimator(vs3, mandelBulb_params) - MandelBulbDistanceEstimator(vd3, mandelBulb_params);
  VEC(pix_data.normal,
    nx,//DE( vector_sum(p,e1)) - DE(vector_diff(p,e1)), 
    ny,//DE(vector_sum(p,e2)) - DE(vector_diff(p,e2)), 
    nz//DE(vector_sum(p,e3)) - DE(vector_diff(p,e3)) 
  );

  NORMALIZE(pix_data.normal);

/* END INLINE NORMAL */


    }
  else {
    //we have the background colour
    pix_data.escaped = true;    
  }
  
  //return dist;




/* END INLINE RAYMARCH */




      //render the pixel
     // rayMarch(renderer_params, from, to, eps, pix_data, mandelBulb_params);

      vec3 pd_norm, pd_hit;//= {pix_data.hit.x, pix_data.normal.x, pix_data.hit.z};
      //double orly = pix_data.normal.y;
      //bool lol = pix_data.escaped;
      VEC(pd_norm, pix_data.normal.x, pix_data.normal.y, pix_data.normal.z); 
     // VEC(pd_hit, pix_data.hit.x, pix_data.hit.y, pix_data.hit.z); 

      color_ptr->x = 124.00;
      color_ptr->y = 33.00;
      color_ptr->z = 111.00;





      /* INLINE GET COLOR */

           //get the colour at this pixel
     // getColour(renderer_params, pd_norm, pd_hit, pix_data.escaped, 
     //   //renderer_params.colourType, renderer_params.brightness 
     //   from, to, color_ptr
     //   );  



      const vec3 baseColor = {1.0, 1.0, 1.0};
      const vec3 backColor = {0.4, 0.4, 0.4};

            //colouring and lightning
            //vec3 hitColor;
            //SET_POINT(hitColor, baseColor);

      if (pix_data.escaped == false) 
      {

                /* MANUALLY INLINED LIGHTING 
                //apply lighting
                lighting(normal, hitColor, direction, result); //, hitColor);
                inline void lighting(const vec3 &n, const vec3 &color, const vec3 &direction,  vec3* outV)

                */

                const vec3 CamLight = {1.0, 1.0, 1.0};

            const double CamLightW = 1.8;// 1.27536;
            const double CamLightMin = 0.3;// 0.48193;

            //vec3 nn = n -1.0;
            vec3 nn = {pd_norm.x - 1.0, pd_norm.y - 1.0, pd_norm.z - 1.0};

            // DOT double, vector
            double dot_res = nn.x*to.x + nn.y*to.y + nn.z*to.z;
            //double ambient = max( CamLightMin, nn.Dot(direction) )*CamLightW;
            double ambient = MAX( CamLightMin, dot_res ) * CamLightW;

            //hitColor = CamLight*ambient*color;
            color_ptr->x = CamLight.x * color_ptr->x * ambient;
            color_ptr->y = CamLight.y * color_ptr->y * ambient;
            color_ptr->z = CamLight.z * color_ptr->z * ambient;



                /*                           */



                //add normal based colouring
            if(colourType == 0 || colourType == 1)
            {

                  //hitColor = hitColor * normal
              color_ptr->x = (color_ptr->x * pd_norm.x + 1.0)/2.0 * brightness;
              color_ptr->y = (color_ptr->y * pd_norm.y + 1.0)/2.0 * brightness;
              color_ptr->z = (color_ptr->z * pd_norm.z + 1.0)/2.0 * brightness;
                  //MULTIPLY_BY_VECTOR(hitColor, normal);
                  //VEC(hitColor, (hitColor.x + 1.0)/2.0, (hitColor.y + 1.0)/2.0, (hitColor.z + 1.0)/2.0 );
                  //MULTIPLY_BY_DOUBLE(hitColor, render_params.brightness);


                  //gamma correction
              v_clamp(color_ptr, 0.0, 1.0);
              color_ptr->x = color_ptr->x * color_ptr->x;
              color_ptr->y = color_ptr->y * color_ptr->y;
              color_ptr->z = color_ptr->z * color_ptr->z;
                  //SQUARE(hitColor);

            }
            
            if(colourType == 1)
            {
              //"swap" colors
              double t = color_ptr->x;//hitColor.x;
              color_ptr->x = color_ptr->z;//hitColor.x = hitColor.z;
              color_ptr->z = t;//hitColor.z = t;
            }
            
          }
          else {
              //we have the background colour
            color_ptr->x = backColor.x;
            color_ptr->y = backColor.y;
            color_ptr->z = backColor.z;


             // SET_POINT(hitColor, backColor);
          }


      /* END GET COLOR */




      //save colour into texture
          k = (j * width + i)*3;
      image[k+2] = (unsigned char)(color_ptr->x *255);//(color.x * 255);
      image[k+1] = (unsigned char)(color_ptr->y *255);//(color.y * 255);
      image[k]   = (unsigned char)(color_ptr->z *255);//(color.z * 255);

    }

    #ifndef _OPENACC
    printProgress((j+1)/(double)height,getTime()-time, frame);
    #endif
  }

}//end pragma

printf("\n rendering done:\n");

}
