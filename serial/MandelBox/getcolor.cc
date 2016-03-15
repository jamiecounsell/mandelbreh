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

using namespace std;

//---lightning and colouring---------
static vec3 CamLight(1.0,1.0,1.0);
static double CamLightW = 1.8;// 1.27536;
static double CamLightMin = 0.3;// 0.48193;
//-----------------------------------
static const vec3 baseColor(1.0, 1.0, 1.0);
static const vec3 backColor(0.4,0.4,0.4);
//-----------------------------------

void lighting(const vec3 &n, const vec3 &color, const vec3 &pos, const vec3 &direction,  vec3 &outV)
{
  vec3 nn = n -1.0;
  double ambient = max( CamLightMin, nn.Dot(direction) )*CamLightW;
  outV = CamLight*ambient*color;
}

vec3 getColour(const pixelData &pixData, const RenderParams &render_params,
	       const vec3 &from, const vec3  &direction)
{
  //colouring and lightning
  vec3 hitColor = baseColor;
  
  if (pixData.escaped == false) 
    {
      //apply lighting
      lighting(pixData.normal, hitColor, pixData.hit, direction, hitColor);
      
      //add normal based colouring
      if(render_params.colourType == 0 || render_params.colourType == 1)
	{
	  hitColor = hitColor * pixData.normal;
	  hitColor = (hitColor + 1.0)/2.0;
	  hitColor = hitColor*render_params.brightness;
	  
	  //gamma correction
	  clamp(hitColor, 0.0, 1.0);
	  hitColor = hitColor*hitColor;
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
    hitColor = backColor;
  
  return hitColor;
}
