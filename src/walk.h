#ifndef _walk_H
#define _walk_h
#include "camera.h"

#define PRINTVEC(vec, end) ( printf("(%f, %f, %f)%s", vec.x, vec.y, vec.z, end) )

void walk(CameraParams *camera_params, RenderParams *renderer_params,
           MandelBulbParams *bulb_params);

#endif
