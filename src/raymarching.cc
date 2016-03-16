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
#include <assert.h>
#include <algorithm>
#include <stdio.h>

#include "mandelbulb.h"
#include "color.h"
#include "renderer.h"

#ifdef _OPENACC
#include <openacc.h>
#endif

#pragma acc routine seq
extern double MandelBulbDistanceEstimator(const vec3 &p, const MandelBulbParams &params);


inline void normal(const vec3 & p, vec3 & normal, const MandelBulbParams &bulb_params)
{
  // compute the normal at p
  const double sqrt_mach_eps = 1.4901e-08;

  double eps ; //= std::max( p.Magnitude(), 1.0 )*sqrt_mach_eps;
  MAGNITUDE(eps, p);
  eps = MAX(eps, 1.0);
  eps *= sqrt_mach_eps;

  vec3 e1 = {eps, 0, 0}; 
  vec3 e2 = {0, eps, 0}; 
  vec3 e3 = {0, 0, eps}; 
  
  //normal = {DE(p+e1)-DE(p-e1), DE(p+e2)-DE(p-e2), DE(p+e3)-DE(p-e3)};//vec3(DE(p+e1)-DE(p-e1), DE(p+e2)-DE(p-e2), DE(p+e3)-DE(p-e3));
  VEC(normal,
    MandelBulbDistanceEstimator(vector_sum(p,e1), bulb_params) - MandelBulbDistanceEstimator(vector_diff(p,e1), bulb_params), 
    MandelBulbDistanceEstimator(vector_sum(p,e2), bulb_params) - MandelBulbDistanceEstimator(vector_diff(p,e2), bulb_params), 
    MandelBulbDistanceEstimator(vector_sum(p,e3), bulb_params) - MandelBulbDistanceEstimator(vector_diff(p,e3), bulb_params) 
  );
  NORMALIZE(normal);
}

#pragma acc routine seq
void rayMarch(const RenderParams &render_params, const vec3 &from, const vec3  &direction, double eps,
	      pixelData& pix_data, const MandelBulbParams &bulb_params)
{
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
        from.x + direction.x * totalDist, 
        from.y + direction.y * totalDist, 
        from.z + direction.z * totalDist
      );

      dist = MandelBulbDistanceEstimator(p, bulb_params);
      
      totalDist += .95*dist;
      
      epsModified = totalDist;
      epsModified*=eps;
      steps++;
    }
  while (dist > epsModified && totalDist <= render_params.maxDistance && steps < render_params.maxRaySteps);
  
  //vec3 hitNormal; UNUSED???

  if (dist < epsModified) 
    {
      //we didnt escape
      pix_data.escaped = false;
      
      // We hit something, or reached MaxRaySteps
      pix_data.hit = p;
      
      //figure out the normal of the surface at this point
      const vec3 normPos = {
        p.x - direction.x * epsModified, 
        p.y - direction.y * epsModified, 
        p.z - direction.z * epsModified
      };
      normal(normPos, pix_data.normal, bulb_params);
    }
  else {
    //we have the background colour
    pix_data.escaped = true;    
  }

  //return dist;
}
