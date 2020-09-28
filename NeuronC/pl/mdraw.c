/* Segment MDRAW in program MONT */

/*	Latest Mod	20-June-83	R.G.Smith	*/

#include "mdef.h"
#include <stdio.h>
#include "stdplt.h"

#define FSCAL  0.00025

/*int siztab[16] = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32767}; */
int siztab[16] = {1,2,3,4,6,8,12,16,24,32,48,64,96,128,198,256};

#include "gr.h"

/*-------------------------------*/

void drcirc (int size, int fill)

#define CHORD 0.5

{
  double radius,theta;
  double chord;

 radius = siztab[size] * FSCAL;
 radius /= 2.0;
 _drcirc(radius,fill);
/*
 rpmove (radius,-PI/2);
 if (size > 5) chord = CHORD / 2.0;
 else          chord = CHORD;
 for (theta=0.0; theta < PI*2; theta += chord)
   rpdraw (radius*chord,theta);
*/
}

/*-------------------------------*/

void drsqr (int size, int fill)

{
  double side;

 side = siztab[size] * FSCAL;
 rmove (-side/2,side/2);
 rdrrect(0.0,-side,side,0.0,0.0,side,-side,0.0,fill);

/*
 rdraw (0.0,-side);
 rdraw (side,0.0);
 rdraw (0.0,side);
 rdraw (-side,0.0);
*/
}

/*-------------------------------*/

void drrect (int size, int aspect, int fill)
{
  double xside,yside;
  static double aspectab[4][2] = {0.1,1.0,
				  0.3,1.0,
				  1.0,0.1,
				  1.0,0.3};

 aspect &= 3;
 xside = siztab[size] * aspectab[aspect][0] * FSCAL;
 yside = siztab[size] * aspectab[aspect][1] * FSCAL;
 rmove (-xside/2,yside/2);
 rdrrect(0.0,-yside,xside,0.0,0.0,yside,-xside,0.0,fill);
/*
 rdraw (0.0,-yside);
 rdraw (xside,0.0);
 rdraw (0.0,yside);
 rdraw (-xside,0.0);
*/
}

/*-------------------------------*/

void drtri (int size, int fill)
{
  double side;

 side = siztab[size] * FSCAL;
 rmove (0.0, side * 0.5775);
 rdraw (-0.5 * side, -0.866 * side);
 rdraw (side, 0.0);
 rdraw (-0.5 * side, 0.866 * side);
} 

/*-------------------------------*/

void drast (int size)

/* asterisk */

{
  double side;

 side = siztab[size] * FSCAL;
 rmove ( 0.25 * side, .433 * side);
 rdraw (-0.5 * side, -0.866 * side);
 rmove ( 0.5 * side,  0.0);
 rdraw (-0.5 * side,  0.866 * side);
 rmove (-0.25 * side, -.433 * side);
 rdraw ( side, 0.0);
} 

/*-------------------------------*/

void drchar (int dchar, int size)
{
  double side,oldside;

 side = siztab[size] * FSCAL;
 oldside = cwidth(side,".");
 rmove (-.3*side,-.3*side);
 textf ("%c",dchar);
 cwidth(oldside,".");
} 

/*-------------------------------*/

void drnum (ssiz *pnum, int size)
{
  double side,oldside;

 side = siztab[size] * FSCAL;
 oldside = cwidth(side,".");
 rmove (-.3*side,-.3*side);
 textf ("%d",*pnum);
 cwidth(oldside,".");
} 

/*-------------------------------*/

void drhor (int size)
{
  double side;

 side = siztab[size] * FSCAL;
 rmove (-0.5 * side, 0.0);
 rdraw (side, 0.0);
}

/*-------------------------------*/

void drvert (int size)
{
  double side;

 side = siztab[size] * FSCAL;
 rmove (0.0,  0.5 * side);
 rdraw (0.0, -side);
}

/*-------------------------------*/

void drcros (int size)
{
  double side;

 side = siztab[size] * FSCAL;
 rmove (0.0,  0.5 * side);
 rdraw (0.0, -side);
 rmove (-0.5 * side, 0.5 * side);
 rdraw ( side, 0.0);
}

/*-------------------------------*/

