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

#include "vector3d.h"
#include "color.h"
#include "mandelbulb.h"
#ifdef _OPENACC
#include <openacc.h>
#include <accelmath.h>
#else
#include <math.h>
#endif


#pragma acc routine seq
double MandelBulbDistanceEstimator(const vec3 &p0, const MandelBulbParams &params)
{
  vec3 z;
  //z = p0;
  SET_POINT(z, p0);

  double dr = 1.0;
  double r = 0.0;

  double Bailout = params.escape_time;//params.rMin;
  double Power = params.power;//params.rFixed;

  for (int i=0; i < params.num_iter; i++) 
    {
      //r = z.Magnitude();//
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

