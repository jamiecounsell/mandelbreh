#include "camera.h"
#include "vector3d.h"
#include <stdio.h>
#include <math.h>

#define STEPSIZE  0.001;
#define TOLERANCE 0.005;

double VECTOR_OPTIONS [5] = {sqrt(1.0/(double)3.0), -sqrt(1.0/(double)3.0), (double)0, (double)1, (double)-1};
vec3 directions [100];

void walk(CameraParams *camera_params) {

    int i, j, k, pos;
    double  x = camera_params->camPos[0],
            y = camera_params->camPos[1],
            z = camera_params->camPos[2];
    // Get current position
    printf("Current position: (%f, %f, %f)\n", x, y, z);

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
    //printf("Distance in +x: %f", )

    // Make decision


    // Change params
    camera_params->camPos[0] = x;
    camera_params->camPos[1] = y;
    camera_params->camPos[2] = z;
    printf("\n");
}
