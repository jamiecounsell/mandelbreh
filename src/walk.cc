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
#else
    #define TOLERANCE (double)0.000002;
    #define STEPSIZE  (double)0.000001;
#endif

#ifdef BULB
    extern double rayMarch(const int maxRaySteps, const float maxDistance,
      const float escape_time, const float power, const int num_iter,
      const vec3 &from, const vec3  &direction, double eps, pixelData& pix_data);

    extern inline double DE(const vec3 &p0, 
      const float escape_time, const float power, const int num_iter);
#else //BOX
    extern double rayMarch(const int maxRaySteps, const float maxDistance,
      const int num_iter, const float rMin, const float rFixed, const float escape_time, const float scale,
      const vec3 &from, const vec3  &direction, double eps, pixelData& pix_data);
    extern inline double DE(const vec3 &p0, const int num_iter, const float rMin, 
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

    vec3 center, pos;
    double  x = camera_history[frame].camPos[0],
            y = camera_history[frame].camPos[1],
            z = camera_history[frame].camPos[2];
    VEC(center, x, y, z);
    VEC(pos, x, y, z);

    double zer [3] = {0,0,0};

    // Circling algorithm
    SUBTRACT_DOUBLE_ARRAY(pos, zer);
    NORMALIZE(pos);

    pixelData pixel;
    const double eps = pow(10.0, renderer_params->detail); 

    #ifdef BULB
        double dist = DE(
            pos,
            bulb_params->escape_time,
            bulb_params->power,
            bulb_params->num_iter);
    #else
        double c1 = fabs(box_params->scale - 1.0);
        double c2 = pow( fabs(box_params->scale), 1 - bulb_params->num_iter);
        double dist = DE(
            pos,
            bulb_params->num_iter,
            box_params->rMin, 
            box_params->rFixed, 
            box_params->escape_time, 
            box_params->scale,
            c1, c2);
    #endif

    double elevation = frame/800.0;
    double inclination = frame/450.0;

    camera_history[frame + 1].camPos[0] = 0.7 * sin(elevation) * cos(inclination);
    camera_history[frame + 1].camPos[1] = 0.7 * sin(elevation) * sin(inclination);
    camera_history[frame + 1].camPos[2] = 0.7 * cos(elevation);

    printf("walk complete\n");

    // More complex walk algorithm
    // int i, j, k, pos;
    // vec3 current;

    // // Get current position
    // VEC(current, x, y, z);
    // if (verbose){
    //     printf("Current position:"); PRINTVEC(current, "\n");
    // }
    // // Get possible movement vectors
    // // TODO: move this to main - not necessary to recalculate
    // for (i = 0; i < 4; i ++) {
    //     for (j = 0; j < 4; j++) {
    //         for (k = 0; k < 4; k++) {
    //             pos = (i * 4) + (j * 4) + k;
    //             directions[pos].x= VECTOR_OPTIONS[i];
    //             directions[pos].y= VECTOR_OPTIONS[j];
    //             directions[pos].z= VECTOR_OPTIONS[k];
    //         }
    //     }
    // }

    // // Make decision for next place
    // pixelData pixel;
    // const double eps = pow(10.0, renderer_params->detail); 
    // int next;
    // double step = STEPSIZE;
    // double tol = TOLERANCE;
    // double closed = 1000;
    // for (i = 0; i < 28; i++){
    //     SUBTRACT_DOUBLE_ARRAY(directions[i], camera_history[frame].camPos);
    //     NORMALIZE( directions[i] );

    //     #ifdef BULB
    //     double dist = rayMarch(
    //         renderer_params->maxRaySteps,
    //         renderer_params->maxDistance,
    //         bulb_params->escape_time,
    //         bulb_params->power,
    //         bulb_params->num_iter,
    //         current,
    //         directions[i],
    //         eps,
    //         pixel); 
    //     #else
    //     double dist = rayMarch(
    //         renderer_params->maxRaySteps,
    //         renderer_params->maxDistance,
    //         box_params->num_iter, 
    //         box_params->rMin, 
    //         box_params->rFixed, 
    //         box_params->escape_time, 
    //         box_params->scale,
    //         current,
    //         directions[i],
    //         eps,
    //         pixel); 
    //     #endif 
    //     dist = fabs(dist);
    //     PRINTVEC(directions[i],":");printf("%f\n", dist);

    //     if (    // if distance is within tolerance and not background
    //             (
    //                 (tol < dist) && \
    //                 (dist < closed) ) && (

    //                 // and history exists
    //                 frame == 0 || (

    //                     // and next step is different than previous step
    //                     camera_history[frame - 1].camPos[0] != x + step * directions[i].x &&
    //                     camera_history[frame - 1].camPos[1] != y + step * directions[i].y &&
    //                     camera_history[frame - 1].camPos[2] != z + step * directions[i].z
    //                 )
    //             )
    //         ) { closed = dist; next = i; }
    // }
    // if (verbose){
    //     printf("Closest point: "); PRINTVEC(directions[next], " "); printf("distance: %f\n", closed);
    // }

    // // Change params
    // camera_history[frame + 1].camPos[0] = x + step * directions[next].x;
    // camera_history[frame + 1].camPos[1] = y + step * directions[next].y;
    // camera_history[frame + 1].camPos[2] = z + step * directions[next].z;

    printf("\n");
}
