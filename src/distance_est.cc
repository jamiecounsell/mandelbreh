#include "vector3d.h"
#include "color.h"
#include "mandelbulb.h"
#ifdef _OPENACC
#include <openacc.h>
#include <accelmath.h>
#else
#include <math.h>
#endif

// set by main
extern MandelBulbParams mandelBulb_params;
#pragma acc declare copyin(mandelBulb_params);

//#pragma acc routine seq present_or_copyin(mandelBulb_params)
//extern double MandelBulbDistanceEstimator(const vec3 &p0, const MandelBulbParams &params);

#pragma acc routine seq
inline double MandelBulbDistanceEstimator(const vec3 &p0, const MandelBulbParams &params)
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


//Distance Estimator Field Selector
double DE(const vec3 &p)
{
  double d = 55.0;//MandelBulbDistanceEstimator(p, mandelBulb_params);
  return d;
}

