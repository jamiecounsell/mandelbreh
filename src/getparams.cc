/*
   This file is part of the Mandelbox program developed for the course
    CS/SE  Distributed Computer Systems taught by N. Nedialkov in the
    Winter of 2015-2016 at McMaster University.

    Copyright (C) 2015-2016 T. Gwosdz and N. Nedialkov
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "renderer.h"
#include "mandelbulb.h"
#include "mandelbox.h"
#include "camera.h"

#define BUF_SIZE 1024

static char buf[BUF_SIZE];

#ifdef BULB
void getParameters(char *filename, CameraParams *camP, RenderParams *renP, MandelBulbParams *bulbP)
#else
void getParameters(char *filename, CameraParams *camP, RenderParams *renP, MandelBoxParams *boxP)
#endif
{
  FILE *fp;
  int ret;
  double *d;
  
  renP->fractalType = 0;
  renP->maxRaySteps = 8000;
  renP->maxDistance = 1000;

  fp = fopen(filename,"r");
  
  if( !fp ) 
    {
      printf(" *** File %s does not exist\n", filename);
      exit(1);
    }
  
  int count = 0;
  
  while (1)
    {
      memset(buf, 0, BUF_SIZE);
      
      ret = fscanf(fp, "%1023[^\n]\n", buf);
      if (ret == EOF) break;
      
      if(buf[0] == '#') // comment line
	continue;
      
      switch(count)
	{
	  // CAMERA
	  //camera position
	case 0:
	  d = camP->camPos;
	  sscanf(buf, "%lf %lf %lf", d, d+1, d+2);
	  break;
	case 1:
	  //camera target
	  d = camP->camTarget;
	  sscanf(buf, "%lf %lf %lf", d, d+1, d+2);
	  break;
	  //camera up 
	case 2:
	  d = camP->camUp;
	  sscanf(buf, "%lf %lf %lf", d, d+1, d+2);
	  break;
	  //field of view
	case 3:
	  sscanf(buf, "%lf", &camP->fov);
	  break;
	  
	  //IMAGE
	  //width, height
	case 4:
	  sscanf(buf, "%d %d", &renP->width, &renP->height);
	  break;
	  //detail
	case 5:
	  sscanf(buf, "%f", &renP->detail);
	  break;
	  
	  //FRACTAL
	case 6:
    // box: scale, rmin, rfixed
    // bulb: IGNORE, escape time(bailout), power
	  //sscanf(buf, "%f %f %f", &boxP->scale, &boxP->rMin, &boxP->rFixed);
  #ifdef BULB
 	  sscanf(buf, "%*f %f %f", &bulbP->escape_time, &bulbP->power);
  #else
    sscanf(buf, "%f %f %f", &boxP->scale, &boxP->rMin, &boxP->rFixed);
  #endif
    break;
	  
	case 7:
	  //sscanf(buf, "%d %f ", &boxP->num_iter, &boxP->escape_time);
    // bulb: max iterations, IGNORE
  #ifdef BULB
    sscanf(buf, "%d %*f ", &bulbP->num_iter);
  #else	  
    sscanf(buf, "%d %f ", &boxP->num_iter, &boxP->escape_time);
  #endif
    break;
	  
	  //COLORING 
	case 8:
	  sscanf(buf, "%d", &renP->colorType);
	  break;
	case 9:
    // brightness
	  sscanf(buf, "%f ", &renP->brightness);
	  break;
	  //FILENAME
	case 10:
	  strcpy(renP->file_name, buf);
	  break;
	}
      count++;
    }
  fclose(fp);
}

