#include "camera.h"
#include "vector3d.h"
#include "mandelbulb.h"
#include "mandelbox.h"
#include "color.h"
#include "camera.h"
#include "renderer.h"
#include "walk.h"

#include <stdio.h>

#ifdef BULB
    #define TOLERANCE (double)0.0001;
    #define STEPSIZE  (double)0.00005;

    extern double rayMarch(const int maxRaySteps, const float maxDistance,
      const float escape_time, const float power, const int num_iter,
      const vec3 &from, const vec3  &direction, double eps, pixelData& pix_data);

    extern double DE(const vec3 &p0, 
        const float escape_time, const float power, const int num_iter);

#else //BOX
    #define TOLERANCE (double)0.000002;
    #define STEPSIZE  (double)0.000001;

    extern double rayMarch(const int maxRaySteps, const float maxDistance,
      const int num_iter, const float rMin, const float rFixed, const float escape_time, const float scale,
      const vec3 &from, const vec3  &direction, double eps, pixelData& pix_data);

    extern double DE(const vec3 &p0, const int num_iter, const float rMin, 
        const float rFixed, const float escape_time, const float scale, double c1, double c2);
#endif

double VECTOR_OPTIONS [4] = {sqrt(1.0/(double)3.0), -sqrt(1.0/(double)3.0), (double)1, (double)-1};
vec3 directions [28];

#ifdef BULB
void walk(CameraParams *camera_history, 
            RenderParams *renderer_params,
            MandelBulbParams *bulb_params,
            int verbose, int frame) 
#else
void walk(CameraParams *camera_history, 
            RenderParams *renderer_params,
            MandelBoxParams *box_params,
            int verbose, int frame)
#endif
{

    double elevation = frame/800.0;
    double inclination = frame/450.0;
    double dist;

    vec3 p;    
    VEC(p,
        camera_history[frame].camPos[0],
        camera_history[frame].camPos[1],
        camera_history[frame].camPos[2]);

    #ifdef BULB
    dist = DE(p, bulb_params->escape_time, bulb_params->power, bulb_params->num_iter);
    #else
    double c1 = fabs(box_params->scale - 1.0);
    double c2 = pow( fabs(box_params->scale), 1 - box_params->num_iter);
    dist = DE(p, box_params->num_iter, box_params->rMin, 
        box_params->rFixed, box_params->escape_time, box_params->scale, c1, c2);
    #endif

    printf("frame %d dist: %lf\n", frame, dist);

    camera_history[frame + 1].camPos[0] = 0.7 * sin(elevation) * cos(inclination);
    camera_history[frame + 1].camPos[1] = 0.7 * sin(elevation) * sin(inclination);
    camera_history[frame + 1].camPos[2] = 0.7 * cos(elevation);
}
