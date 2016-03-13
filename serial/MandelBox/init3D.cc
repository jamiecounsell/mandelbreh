/*
   This file is part of the Mandelbox program developed for the course
    CS/SE  Distributed Computer Systems taught by N. Nedialkov in the
    Winter of 2015-2016 at McMaster University.

    Copyright (C) 2015-2016 T. Gwosdz and N. Nedialkov

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "camera.h"
#include "renderer.h"
#include "3d.h"

void init3D(CameraParams *camP, const RenderParams *renP)
{
  //set up the viewport for the image
  camP->viewport[0] = 0;
  camP->viewport[1] = 0;
  camP->viewport[2] = renP->width;
  camP->viewport[3] = renP->height;
  
  //init the matricies
  LoadIdentity(camP->matModelView);
  LoadIdentity(camP->matProjection);
  
  //setting up camera lense
  Perspective((65*camP->fov), ((double)renP->width)/((double)renP->height), NEAR, FAR, camP->matProjection);
  
  //setting up model view matrix
  LookAt(camP->camPos, camP->camTarget, camP->camUp, camP->matModelView);

  //setting up the inverse(projection x model) matrix
  double temp[16];
  MultiplyMatrices(temp, camP->matProjection, camP->matModelView);
  //Now compute the inverse of matrix A
  InvertMatrix(temp, camP->matInvProjModel);
}
