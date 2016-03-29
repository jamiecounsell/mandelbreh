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
#include "mandelbox.h"

#ifdef _OPENACC
#include <openacc.h>
#include <accelmath.h>
#else
#include <math.h>
#endif

//copysign(x, y) : magnitude x, sign of y
inline double copysign(double x, double y){

  if(y < -0.00000000000001){
    return -fabs(x);
  }else{
    return fabs(x);
  }


}

#define SQR(x) ((x)*(x))


#define COMPONENT_FOLD(x) { (x) = fabs(x) <= 1? (x) : copysign(2,(x))-(x); }


inline double MandelBoxDE(const vec3 &p0, const int num_iter, const float rMin, 
  const float rFixed, const float escape_time, const float scale, double c1, double c2)
{
  vec3 p = p0;
  double rMin2   = SQR(rMin);
  double rFixed2 = SQR(rFixed);
  double escape  = SQR(escape_time);
  double dfactor = 1; 
  double r2      =-1;
  const double rFixed2rMin2 = rFixed2/rMin2;

  int i = 0;
  while (i< num_iter && r2 < escape)
    {
      COMPONENT_FOLD(p.x);
      COMPONENT_FOLD(p.y);
      COMPONENT_FOLD(p.z);
      
      DOT(r2,p);      

      if (r2<rMin2)
  {
    p.x = p.x * (rFixed2rMin2);
    p.y = p.y * (rFixed2rMin2);
    p.z = p.z * (rFixed2rMin2);

    dfactor *= (rFixed2rMin2);
  }
      else
      if ( r2<rFixed2) 
  {
    const double t = (rFixed2/r2);
    p.x = p.x * (rFixed2/r2);
    p.y = p.y * (rFixed2/r2);
    p.z = p.z * (rFixed2/r2);

    dfactor *= t;
  }
      

      dfactor = dfactor*fabs(scale)+1.0;      
      p.x = p.x * scale + p0.x;
      p.y = p.y * scale + p0.y;
      p.z = p.z * scale + p0.z;

      i++;
    }

  double r = 0.0;
  MAGNITUDE(r, p);
  r -= c1;
  r = r / dfactor;
  r -= c2;
  
  return  r;
}


#pragma acc routine seq
double rayMarch(const int maxRaySteps, const float maxDistance,
 const int num_iter, const float rMin, const float rFixed, const float escape_time, const float scale,
 const vec3 &from, const vec3  &direction, double eps, pixelData& pix_data)
{
  double dist = 0.0;
  double totalDist = 0.0;

  double c1 = fabs(scale - 1.0);
  double c2 = pow( fabs(scale), 1 - num_iter);

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

      dist = MandelBoxDE(p, num_iter, rMin, 
        rFixed, escape_time, scale, c1, c2);

      
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
      
      pix_data.normal.x = MandelBoxDE(vs1, num_iter, rMin, 
        rFixed, escape_time, scale, c1, c2)-MandelBoxDE(vd1, num_iter, rMin, 
        rFixed, escape_time, scale, c1, c2); 
      pix_data.normal.y = MandelBoxDE(vs2, num_iter, rMin, 
        rFixed, escape_time, scale, c1, c2)-MandelBoxDE(vd2, num_iter, rMin, 
        rFixed, escape_time, scale, c1, c2); 
      pix_data.normal.z = MandelBoxDE(vs3, num_iter, rMin, 
        rFixed, escape_time, scale, c1, c2)-MandelBoxDE(vd3, num_iter, rMin, 
        rFixed, escape_time, scale, c1, c2);
      
      NORMALIZE(pix_data.normal);
    }
  else {
    //we have the background color
    pix_data.escaped = true;
    return 0;
  }

  return dist;
}
