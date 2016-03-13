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
#include <stdio.h>
#include <string.h>
#include <math.h>

void printProgress( double perc, double time )
{
  static char delete_space[80];
  static char * OutputString;
  perc *= 100;

  int sec = ceil(time);
  int hr = sec/3600;
  int t = sec%3600;
  int min = t/60;
  sec = t%60;
  
  OutputString = (char*)"*** completed % 5.2f%s ---  cum. time = %02d:%02d:%02d   % e (s)";
  sprintf(delete_space, OutputString, perc, "%%", hr, min, sec, time);

  fprintf( stderr, delete_space);
  for ( unsigned int i = 0; i < strlen(delete_space); i++)
    fputc( 8, stderr);
}
