#ifndef vec3_h
#define vec3_h

#ifdef _OPENACC
#include <accelmath.h>
#else
#include <math.h>
#endif


typedef struct 
{
  double x, y, z;
}  vec3;

#define SET_POINT(p,v) { p.x=v.x; p.y=v.y; p.z=v.z; }

#define SET_DOUBLE_POINT(p,v) { p.x=v[0]; p.y=v[1]; p.z=v[2]; }

#define SUBTRACT_POINT(p,v,u)			\
  {						\
  p.x=(v[0])-(u[0]);				\
  p.y=(v[1])-(u[1]);				\
  p.z=(v[2])-(u[2]);				\
}

#define SQUARE(p)\
	{\
		p.x = p.x * p.x; \
		p.y = p.y * p.y; \
		p.z = p.z * p.z; \
	}


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


#define MULTIPLY_BY_VECTOR(v, p) ( { v.x = v.x*p.x; v.y = v.y*p.y; v.z = v.z*p.z; } )
#define MULTIPLY_BY_DOUBLE(v, d) ( { v.x = v.x*d; v.y = v.y*d; v.z = v.z*d; } )


#define MAGNITUDE(m,p) 	({ m=sqrt( p.x*p.x + p.y*p.y + p.z*p.z ); })

#define DOT(d,p) ({  d= p.x*p.x + p.y*p.y + p.z*p.z ; })

#define MAX(a,b) ( ((a)>(b))? (a):(b) )

#define VEC(v,a,b,c) { v.x = a; v.y = b; v.z = c; }


inline double clamp(double d, double min, double max) 
{
  const double t = d < min ? min : d;
  return t > max ? max : t;
}

inline vec3 vector_sum(vec3 v1, vec3 v2){
	vec3 result = {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
	return result;
}

inline vec3 vector_diff(vec3 v1, vec3 v2){
	vec3 result = {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
	return result;
}


//#define VECTOR_CLAMP(v, min, max) { v.x = clamp(v.x,min,max); v.y = clamp(v.y,min,max); v.z = clamp(v.z,min,max); }

inline void clamp(vec3 &v, double min, double max) 
{
  v.x = clamp(v.x,min,max);
  v.y = clamp(v.y,min,max);
  v.z = clamp(v.z,min,max);
}


#endif
