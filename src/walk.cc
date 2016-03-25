#include "camera.h"
#include "vector3d.h"
#include "mandelbulb.h"
#include "color.h"
#include "camera.h"
#include "renderer.h"
#include "walk.h"

#include <stdio.h>

#define STEPSIZE  0.001;
#define TOLERANCE 0.005;

double VECTOR_OPTIONS [5] = {sqrt(1.0/(double)3.0), -sqrt(1.0/(double)3.0), (double)0, (double)1, (double)-1};
vec3 directions [125];

extern double rayMarch(const int maxRaySteps, const float maxDistance,
  const float escape_time, const float power, const int num_iter,
  const vec3 &from, const vec3  &direction, double eps, pixelData& pix_data);

void walk(CameraParams *camera_params, RenderParams *renderer_params,
           MandelBulbParams *bulb_params) {

    int i, j, k, pos;
    vec3 current;

    double  x = camera_params->camPos[0],
            y = camera_params->camPos[1],
            z = camera_params->camPos[2];
    // Get current position
    VEC(current, x, y, z);
    printf("Current position:");PRINTVEC(current, "\n");

    // Measure distance to object along other vectors

    for (i = 0; i < 5; i ++) {
        for (j = 0; j < 5; j++) {
            for (k = 0; k < 5; k++) {
                pos = (i * 25) + (j * 5) + k;
                directions[pos].x= VECTOR_OPTIONS[i];
                directions[pos].y= VECTOR_OPTIONS[j];
                directions[pos].z= VECTOR_OPTIONS[k];
            }
        }
    }

    pixelData pixel;
    const double eps = pow(10.0, renderer_params->detail); 
    vec3 maxdir;
    double maxd = 0;
    for (i = 0; i < 125; i++){
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
        if (dist > maxd) { maxd = dist; SET_POINT(maxdir, directions[i]) }
    }
    printf("Farthest point: "); PRINTVEC(maxdir, ""); printf("distance: %f\n", maxd);
    // Make decision


    // Change params
    camera_params->camPos[0] = x;
    camera_params->camPos[1] = y;
    camera_params->camPos[2] = z;
    printf("\n");
}
