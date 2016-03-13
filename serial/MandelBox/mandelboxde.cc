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
#include <cmath>
#include <cstdio>
#include <algorithm>
#include "color.h"
#include "mandelbox.h"
#include <math.h>

#define SQR(x) ((x)*(x))

#define COMPONENT_FOLD(x) { (x) = fabs(x) <= 1? (x) : copysign(2,(x))-(x); }


double MandelBoxDE(const vec3 &p0, const MandelBoxParams &params, double c1, double c2)
{
  vec3 p = p0;
  double rMin2   = SQR(params.rMin);
  double rFixed2 = SQR(params.rFixed);
  double escape  = SQR(params.escape_time);
  double dfactor = 1; 
  double r2      =-1;
  const double rFixed2rMin2 = rFixed2/rMin2;

  int i = 0;
  while (i< params.num_iter && r2 < escape)
    {
      COMPONENT_FOLD(p.x);
      COMPONENT_FOLD(p.y);
      COMPONENT_FOLD(p.z);
      
      r2 = p.Dot(p);      

      if (r2<rMin2)
	{
	  p = p*(rFixed2rMin2);
	  dfactor *= (rFixed2rMin2);
	}
      else
      if ( r2<rFixed2) 
	{
	  const double t = (rFixed2/r2);
	  p = p*(rFixed2/r2);;
	  dfactor *= t;
	}
      

      dfactor = dfactor*fabs(params.scale)+1.0;      
      p = p*params.scale+p0;
      i++;
    }
  
  return  (p.Magnitude() - c1) / dfactor - c2;
}


