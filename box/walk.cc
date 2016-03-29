#include "camera.h"
#include "vector3d.h"
#include "mandelbox.h"
#include "color.h"
#include "camera.h"
#include "renderer.h"
#include "walk.h"

#include <stdio.h>

#define STEPSIZE  0.001 * 10;
#define TOLERANCE 0.005 * 10;

double VECTOR_OPTIONS [5] = {sqrt(1.0/(double)3.0), -sqrt(1.0/(double)3.0), (double)1, (double)-1};
vec3 directions [125];

extern double rayMarch(const int maxRaySteps, const float maxDistance,
  const int num_iter, const float rMin, const float rFixed, const float escape_time, const float scale,
  const vec3 &from, const vec3  &direction, double eps, pixelData& pix_data);

void walk(CameraParams *camera_params, RenderParams *renderer_params,
           MandelBoxParams *box_params) {

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
    vec3 nextdir;
    float tol = TOLERANCE;
    double closed = 100;
    for (i = 0; i < 125; i++){
        SUBTRACT_DOUBLE_ARRAY(directions[i], camera_params->camPos);
        NORMALIZE( directions[i] );
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
        if (tol < dist < closed) { closed = dist; SET_POINT(nextdir, directions[i]) }
    }
    printf("Farthest point: "); PRINTVEC(nextdir, ""); printf("distance: %f\n", closed);
    // Make decision

    float step = STEPSIZE;
    // Change params
    camera_params->camPos[0] = x + step * nextdir.x;
    camera_params->camPos[1] = y + step * nextdir.y;
    camera_params->camPos[2] = z + step * nextdir.z;
    printf("\n");
}
