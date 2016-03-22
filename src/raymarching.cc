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

#include "color.h"
#include "renderer.h"
#include "mandelbulb.h"

#ifdef _OPENACC
#include <openacc.h>
#include <accelmath.h>
#else
#include <math.h>
#endif


inline double MandelBulbDistanceEstimator(const vec3 &p0, 
  const float escape_time, const float power, const int num_iter)
{
  vec3 z;
  //z = p0;
  SET_POINT(z, p0);
  
  double dr = 1.0;
  double r = 0.0;

  double Bailout = escape_time;//params.rMin;
  double Power = power;//params.rFixed;

  for (int i=0; i < num_iter; i++) 
    {
      MAGNITUDE(r,z);
      if(r > Bailout) break; 

      double theta = acos(z.z/r);
      double phi   = atan2(z.y, z.x);
      dr = pow(r, Power - 1.0) * Power * dr + 1.0;

      double zr = pow(r, Power);
      theta     = theta * Power;
      phi       = phi * Power;

      z.x = zr*sin(theta)*cos(phi);
      z.y = zr*sin(phi)*sin(theta);
      z.z = zr*cos(theta);

      z.x = z.x + p0.x;
      z.y = z.y + p0.y;
      z.z = z.z + p0.z;
    }

  return 0.5*log(r)*r/dr;
}

#pragma acc routine seq
void rayMarch(const int maxRaySteps, const float maxDistance,
 const float escape_time, const float power, const int num_iter,
 const vec3 &from, const vec3  &direction, double eps, pixelData& pix_data)
{
  double dist = 0.0;
  double totalDist = 0.0;

  const double sqrt_mach_eps = 1.4901e-08;

  
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

      //dist = MandelBulbDistanceEstimator(p, bulb_params);
      dist = MandelBulbDistanceEstimator(p, escape_time, power, num_iter);

      
      totalDist += .95*dist;
      
      epsModified = totalDist;
      epsModified*=eps;
      steps++;
    }
  while (dist > epsModified && totalDist <= maxDistance && steps < maxRaySteps);
  
  //vec3 hitNormal; unused

  if (dist < epsModified) 
    {
      //we didnt escape
      pix_data.escaped = false;
      
      // We hit something, or reached MaxRaySteps
      pix_data.hit = p;
      
      //figure out the normal of the surface at this point
      //const vec3 normPos = p - direction * epsModified;
      const vec3 normPos = {
        p.x - direction.x * epsModified, 
        p.y - direction.y * epsModified, 
        p.z - direction.z * epsModified
      };
      

      // compute the normal at p
      double eps;
      MAGNITUDE(eps, normPos) ;// std::max( p.Magnitude(), 1.0 )*sqrt_mach_eps;
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
      
      pix_data.normal.x = MandelBulbDistanceEstimator(vs1, escape_time, power, num_iter)-MandelBulbDistanceEstimator(vd1, escape_time, power, num_iter); 
      pix_data.normal.y = MandelBulbDistanceEstimator(vs2, escape_time, power, num_iter)-MandelBulbDistanceEstimator(vd2, escape_time, power, num_iter); 
      pix_data.normal.z = MandelBulbDistanceEstimator(vs3, escape_time, power, num_iter)-MandelBulbDistanceEstimator(vd3, escape_time, power, num_iter);
      
      NORMALIZE(pix_data.normal);
      //normal(bulb_params, normPos, pix_data.normal);

    }
  else 
    //we have the background colour
    pix_data.escaped = true;
}
