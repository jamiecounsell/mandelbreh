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
    #define STEPSIZE  0.001;
    #define TOLERANCE 0.005;
#else
    #define STEPSIZE  0.001 * 50;
    #define TOLERANCE 0.005 * 50;
#endif
double VECTOR_OPTIONS [4] = {sqrt(1.0/(double)3.0), -sqrt(1.0/(double)3.0), (double)1, (double)-1};
vec3 directions [24];

#ifdef BULB
    extern double rayMarch(const int maxRaySteps, const float maxDistance,
      const float escape_time, const float power, const int num_iter,
      const vec3 &from, const vec3  &direction, double eps, pixelData& pix_data);
#else //BOX
    extern double rayMarch(const int maxRaySteps, const float maxDistance,
      const int num_iter, const float rMin, const float rFixed, const float escape_time, const float scale,
      const vec3 &from, const vec3  &direction, double eps, pixelData& pix_data);
#endif


#ifdef BULB
void walk(CameraParams *camera_params, RenderParams *renderer_params,
           MandelBulbParams *bulb_params, int verbose) 
#else
void walk(CameraParams *camera_params, RenderParams *renderer_params,
           MandelBoxParams *box_params, int verbose)
#endif
{
    int i, j, k, pos;
    vec3 current;
    double  x = camera_params->camPos[0],
            y = camera_params->camPos[1],
            z = camera_params->camPos[2];
    // Get current position
    VEC(current, x, y, z);
    if (verbose){
        printf("Current position:");PRINTVEC(current, "\n");
    }

    // Get possible movement vectors
    // TODO: move this to main - not necessary to recalculate
    for (i = 0; i < 4; i ++) {
        for (j = 0; j < 4; j++) {
            for (k = 0; k < 4; k++) {

                pos = (i * 4) + (j * 4) + k;
                directions[pos].x= VECTOR_OPTIONS[i];
                directions[pos].y= VECTOR_OPTIONS[j];
                directions[pos].z= VECTOR_OPTIONS[k];

                PRINTVEC(directions[pos], "\n");
            }
        }
    }

    // Make decision for next place
    pixelData pixel;
    const double eps = pow(10.0, renderer_params->detail); 
    vec3 nextdir;
    float tol = TOLERANCE;
    double closed = 1000;
    for (i = 0; i < 125; i++){
        SUBTRACT_DOUBLE_ARRAY(directions[i], camera_params->camPos);
        NORMALIZE( directions[i] );
        #ifdef BULB
        double dist = rayMarch(
            renderer_params->maxRaySteps,
            renderer_params->maxDistance,
            bulb_params->escape_time,
            bulb_params->power,
            bulb_params->num_iter,
            current,
            directions[i],
            eps,
            pixel); 
        #else
        double dist = rayMarch(
            renderer_params->maxRaySteps,
            renderer_params->maxDistance,
            box_params->num_iter, 
            box_params->rMin, 
            box_params->rFixed, 
            box_params->escape_time, 
            box_params->scale,
            current,
            directions[i],
            eps,
            pixel); 

        #endif 
        if (tol < dist < closed) { closed = dist; SET_POINT(nextdir, directions[i]) }
    }
    if (verbose){
        printf("Farthest point: "); PRINTVEC(nextdir, " "); printf("distance: %f\n", closed);
    }

    float step = STEPSIZE;
    // Change params
    camera_params->camPos[0] = x + step * nextdir.x;
    camera_params->camPos[1] = y + step * nextdir.y;
    camera_params->camPos[2] = z + step * nextdir.z;

    printf("\n");
}
