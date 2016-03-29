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
#include <accelmath.h>
#endif


inline void lighting(const vec3 &n, const vec3 &color, const vec3 &pos, const vec3 &direction,  vec3 &outV)
{
  vec3 CamLight = { 1.0, 1.0, 1.0};
  double CamLightW = 1.8;// 1.27536;
  double CamLightMin = 0.3;// 0.48193;

  vec3 nn = { n.x -1.0, n.y -1, n.z -1 };
  double dot_res = nn.x * direction.x + nn.y * direction.y + nn.z * direction.z;
  double ambient = MAX( CamLightMin, dot_res ) * CamLightW;
  
  outV.x = CamLight.x * ambient * color.x;
  outV.x = CamLight.y * ambient * color.y;
  outV.x = CamLight.z * ambient * color.z;

}

#pragma acc routine seq
void getcolor(const pixelData &pixData, const int colorType, const float brightness, //const RenderParams render_params,
	       const vec3 &from, const vec3  &direction, vec3 &result)
{

  /* COLOR SCHEMES 

      0, 1: default
      2   : red filter
      3   : grayscale
      4   : less flamboyant rainbow
  */


  vec3 baseColor = {1.0, 1.0, 1.0};
  vec3 backColor = {0.4, 0.4, 0.4};

  //coloring and lightning
  vec3 hitColor =  {baseColor.x, baseColor.y, baseColor.z};
  
  if (pixData.escaped == false) 
  {
      //apply lighting
      lighting(pixData.normal, hitColor, pixData.hit, direction, hitColor);
      
      //add normal based coloring
      if(colorType == 0 || colorType == 1 || colorType == 2 || colorType == 3 || colorType == 4)
    	{
        
        hitColor.x = (hitColor.x * pixData.normal.x + 1.0)/2.0 * brightness;
        hitColor.y = (hitColor.y * pixData.normal.y + 1.0)/2.0 * brightness;
        hitColor.z = (hitColor.z * pixData.normal.z + 1.0)/2.0 * brightness;

    	  //gamma correction
    	  v_clamp(hitColor, 0.0, 1.0);
        SQUARE(hitColor);

    	}
      if(colorType == 1)
    	{
    	  //"swap" colors
    	  double t = hitColor.x;
    	  hitColor.x = hitColor.z;
    	  hitColor.z = t;
    	} else if (colorType == 2)
      {
        // red filter
        hitColor.x = 0.0;
      } else if (colorType == 3)
      { 
        // grayscale
        //weighted average used by GIMP based on human perception
        //double avg = 0.21 * hitColor.x + 0.72 * hitColor.y + 0.07 * hitColor.z;
        double avg = (hitColor.x + hitColor.y + hitColor.z) / 3.0;
        hitColor.x = avg;
        hitColor.y = avg;
        hitColor.z = avg; 
      } else 
      {
        // rainbow
        hitColor.x = 0.85 * hitColor.x;
        hitColor.z = 0.7 * hitColor.z;
      }
      
  }
  else {
    //we have the background color
    hitColor = backColor;
  }
  
  result.x = hitColor.x;
  result.y = hitColor.y;
  result.z = hitColor.z;

}
