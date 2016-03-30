#ifndef _walk_H
#define _walk_h
#include "camera.h"

#define PRINTVEC(vec, end) ( printf("(%f, %f, %f)%s", vec.x, vec.y, vec.z, end) )

#ifdef BULB
void walk(CameraParams *camera_params, CameraParams *camera_history,
           RenderParams *renderer_params,
           MandelBulbParams *bulb_params,
           int verbose, int frame);
#else
void walk(CameraParams *camera_params, CameraParams *camera_history,
           RenderParams *renderer_params,
           MandelBoxParams *box_params,
           int verbose, int frame);
#endif

#endif
