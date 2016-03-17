#include "color.h"
#include "renderer.h"
#include "vector3d.h"
//#include <cmath>
//#include <algorithm>

#ifdef _OPENACC
#include <openacc.h>
#endif

//using namespace std;

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

// pos never used
//inline void lighting(const vec3 &n, const vec3 &color, const vec3 &pos, const vec3 &direction,  vec3 &outV)
inline void lighting(const vec3 &n, const vec3 &color, const vec3 &direction,  vec3 &outV)
{

  const vec3 CamLight = {1.0, 1.0, 1.0};

  const double CamLightW = 1.8;// 1.27536;
  const double CamLightMin = 0.3;// 0.48193;

  //vec3 nn = n -1.0;
  vec3 nn = {n.x - 1.0, n.y - 1.0, n.z - 1.0};

  // DOT double, vector
  double dot_res = nn.x*direction.x + nn.y*direction.y + nn.z*direction.z;
  //double ambient = max( CamLightMin, nn.Dot(direction) )*CamLightW;
  double ambient = MAX( CamLightMin, dot_res ) * CamLightW;

  //hitColor = CamLight*ambient*color;
  SET_POINT(outV, color);
  MULTIPLY_BY_VECTOR(outV, CamLight);
  MULTIPLY_BY_DOUBLE(outV, ambient);

}

#pragma acc routine seq
void getColour(const RenderParams &render_params, const vec3 &normal, const vec3 &hit, const bool escaped,
        //const int colourType, const float brightness
        const vec3 &from, const vec3 &direction, vec3 &result
        )
//vec3 getColour(const pixelData &pixData, const RenderParams &render_params,
//	       const vec3 &from, const vec3  &direction)
{

  const vec3 baseColor = {1.0, 1.0, 1.0};
  const vec3 backColor = {0.4, 0.4, 0.4};

  //colouring and lightning
  vec3 hitColor;
  SET_POINT(hitColor, baseColor);

  if (escaped == false) 
  {
   
      //apply lighting
      lighting(normal, hitColor, direction, hitColor);
     
      //add normal based colouring
      if(render_params.colourType == 0 || render_params.colourType == 1)
	   {
	
    // hitColor is corrupted at this point.. fine at end of call to lighting()
    // normal is fine
    //hitColor = hitColor * normal
    //MULTIPLY_BY_VECTOR(hitColor, normal);


    VEC(hitColor, (hitColor.x + 1.0)/2.0, (hitColor.y + 1.0)/2.0, (hitColor.z + 1.0)/2.0 );
	  VEC(hitColor, 
          hitColor.x * render_params.brightness, 
          hitColor.y * render_params.brightness, 
          hitColor.z * render_params.brightness
        );



	  //gamma correction
	  v_clamp(hitColor, 0.0, 1.0);
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
  else {
    //we have the background colour
    SET_POINT(hitColor, backColor);
  }
  
  SET_POINT(result, hitColor);
  //return hitColor;

}
