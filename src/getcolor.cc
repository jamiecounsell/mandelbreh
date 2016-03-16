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
#include "color.h"
#include "renderer.h"
#include "vector3d.h"
#include <cmath>
#include <algorithm>

#ifdef _OPENACC
#include <openacc.h>
#endif

using namespace std;

//---lightning and colouring---------
//static vec3 CamLight(1.0,1.0,1.0);
//static vec3 CamLight = {1.0, 1.0, 1.0};

//static double CamLightW = 1.8;// 1.27536;
//static double CamLightMin = 0.3;// 0.48193;
//-----------------------------------
//static const vec3 baseColor(1.0, 1.0, 1.0);
//static const vec3 backColor(0.4,0.4,0.4);
//static const vec3 baseColor = {1.0, 1.0, 1.0};
//static const vec3 backColor = {0.4, 0.4, 0.4};
//-----------------------------------

inline void lighting(const vec3 &n, const vec3 &color, const vec3 &pos, const vec3 &direction,  vec3 &outV)
{

  const vec3 CamLight = {1.0, 1.0, 1.0};

  const double CamLightW = 1.8;// 1.27536;
  const double CamLightMin = 0.3;// 0.48193;

  //vec3 nn = n -1.0;
  vec3 nn = {n.x - 1.0, n.y - 1.0, n.z - 1.0};

  // DOT double, vector
  double dot_res = nn.x*direction.x + nn.y*direction.y + nn.z*direction.z;
  //double ambient = max( CamLightMin, nn.Dot(direction) )*CamLightW;
  double ambient = max( CamLightMin, dot_res )*CamLightW;

  //outV = CamLight*ambient*color;
  SET_POINT(outV, color);
  MULTIPLY_BY_VECTOR(outV, CamLight);
  MULTIPLY_BY_DOUBLE(outV, ambient);
}

#pragma acc routine seq
vec3 getColour(const pixelData &pixData, const RenderParams &render_params,
	       const vec3 &from, const vec3  &direction)
{

  const vec3 baseColor = {1.0, 1.0, 1.0};
  const vec3 backColor = {0.4, 0.4, 0.4};

  //colouring and lightning
  //vec3 hitColor = baseColor;
  vec3 hitColor;
  SET_POINT(hitColor, baseColor);

  if (pixData.escaped == false) 
    {
      //apply lighting
      lighting(pixData.normal, hitColor, pixData.hit, direction, hitColor);
      
      //add normal based colouring
      if(render_params.colourType == 0 || render_params.colourType == 1)
	{
	  //hitColor = hitColor * pixData.normal;
    VEC(hitColor, hitColor.x * pixData.normal.x, hitColor.y * pixData.normal.y, hitColor.z * pixData.normal.z);
	  //hitColor = (hitColor + 1.0)/2.0;
	  //hitColor = hitColor*render_params.brightness;
    VEC(hitColor, (hitColor.x + 1.0)/2.0, (hitColor.y + 1.0)/2.0, (hitColor.z + 1.0)/2.0 );
	  VEC(hitColor, hitColor.x * render_params.brightness, hitColor.y * render_params.brightness, hitColor.z * render_params.brightness);



	  //gamma correction
	  v_clamp(hitColor, 0.0, 1.0);
	  //hitColor = hitColor*hitColor;
    SQUARE(hitColor);
	}
      if(render_params.colourType == 1)
	{
	  //"swap" colors
	  double t = hitColor.x;
	  hitColor.x = hitColor.z;
	  hitColor.z = t;
	}
    }
  else 
    //we have the background colour
    //hitColor = backColor;
    SET_POINT(hitColor, backColor);

  return hitColor;
}
