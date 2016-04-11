#ifndef vec3_h
#define vec3_h

#ifdef _OPENACC
#include <accelmath.h>
#else
#include <math.h>
#endif

// vec3 has been changed to struct with accompanying macros
typedef struct 
{
  double x, y, z;
}  vec3;

// set vec3 p = v
#define SET_POINT(p,v) { p.x=v.x; p.y=v.y; p.z=v.z; }
// set vec3 x, y, z to double[0..2]
#define SET_DOUBLE_POINT(p,v) { p.x=v[0]; p.y=v[1]; p.z=v[2]; }
// set vec3 x,y,z to v[i] - u[i]
#define SUBTRACT_POINT(p,v,u)			\
  {						\
  p.x=(v[0])-(u[0]);				\
  p.y=(v[1])-(u[1]);				\
  p.z=(v[2])-(u[2]);				\
}
// x,y,z = x,y,z - double[0..2]
#define SUBTRACT_DOUBLE_ARRAY(v, d) {v.x = v.x - d[0]; v.y = v.y - d[1]; v.z = v.z - d[2];  }
// (x,y,z)^2
#define SQUARE(p)\
	{\
		p.x = p.x * p.x; \
		p.y = p.y * p.y; \
		p.z = p.z * p.z; \
	}
// normalize vector
#define NORMALIZE(p) {					\
    double fMag = ( p.x*p.x + p.y*p.y + p.z*p.z );	\
    if (fMag != 0)					\
      {							\
	double fMult = 1.0/sqrt(fMag);			\
	p.x *= fMult;					\
	p.y *= fMult;					\
	p.z *= fMult;					\
      }							\
  }
// x*p, y*q, z*r
#define MULTIPLY_BY_VECTOR(v, p) ( { v.x = v.x*p.x; v.y = v.y*p.y; v.z = v.z*p.z; } )
// x,y,z * d
#define MULTIPLY_BY_DOUBLE(v, d) ( { v.x = v.x*d; v.y = v.y*d; v.z = v.z*d; } )
// get vector magnitude
#define MAGNITUDE(m,p) 	({ m=sqrt( p.x*p.x + p.y*p.y + p.z*p.z ); })
// vector dot product
#define DOT(d,p) ({  d= p.x*p.x + p.y*p.y + p.z*p.z ; })
#define MAX(a,b) ( ((a)>(b))? (a):(b) )
// constructor 
#define VEC(v,a,b,c) { v.x = a; v.y = b; v.z = c; }


inline double clamp(double d, double min, double max) 
{
  const double t = d < min ? min : d;
  return t > max ? max : t;
}

// vector addition and subtraction
#define VECTOR_SUM(r, v1, v2) {  r.x = v1.x + v2.x; r.y = v1.y + v2.y; r.z = v1.z + v2.z; }
#define VECTOR_DIFF(r, v1, v2) {  r.x = v1.x - v2.x; r.y = v1.y - v2.y; r.z = v1.z - v2.z; }
// vector clamp
inline void v_clamp(vec3 &v, double min, double max) 
{
  v.x = clamp(v.x,min,max);
  v.y = clamp(v.y,min,max);
  v.z = clamp(v.z,min,max);
}

#endif
