#include "camera.h"
#include "vector3d.h"
#include "mandelbulb.h"
#include "mandelbox.h"
#include "color.h"
#include "3d.h"
#include "camera.h"
#include "renderer.h"
#include "walk.h"

#include <stdio.h>

// Get the next frame to render
// Computes an orbit around the bulb or box and reinitializes camera to point at center
// The computed orbit is a spiral following a circular path in the x-y plane with a decreasing z
// See documentation for full algorithm description

#ifdef BULB
    extern double rayMarch(const int maxRaySteps, const float maxDistance,
      const float escape_time, const float power, const int num_iter,
      const vec3 &from, const vec3  &direction, double eps, pixelData& pix_data);

    extern double DE(const vec3 &p0, 
        const float escape_time, const float power, const int num_iter);

#else //BOX
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
    // full rotation in x,y every 500*pi frames
    double inclination = frame/500.0;
    // normalize z range from 1 to -1 for 7200 frames
    double correction = frame / 3600.0;

    camera_history[frame + 1].camPos[0] = cos(inclination);
    camera_history[frame + 1].camPos[1] = sin(inclination);
    camera_history[frame + 1].camPos[2] = 1 - correction;

    // set camera for the next frame
    init3D(&camera_history[frame + 1], renderer_params);
}
