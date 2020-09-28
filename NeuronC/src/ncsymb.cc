/* Segment ncsymb in program NC */

/*	Latest Mod	20-June-83	R.G.Smith	*/

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "stdplt.h"
#include "gr.h"

#ifdef __cplusplus
}
#endif

#include "ncio.h"

/*-------------------------------*/

void drcirc (double size, int fill)

#define CHORD 0.5

{
  double radius,theta;
  double chord;

 radius = size/2.0;
 _drcirc(radius,fill);
}

/*-------------------------------*/

void drsqr (double size, int fill)

{
  double side;

 side = size;
 rmove (-side/2,side/2);
 rdrrect(0.0,-side,side,0.0,0.0,side,-side,0.0,fill);

}

/*-------------------------------*/

void drrect (double size, double aspect, int fill)
{
  double xside,yside;

 if (aspect < 0) {
    ncfprintf (stderr,"drrect: negative aspect %g\n",aspect);
    return;
 }

 if (aspect > 1.0) {
     xside = size; 
     yside = size / aspect;
 }
 else {
     xside = size * aspect;
     yside = size;
 }
 rmove (-xside/2,yside/2);
 rdrrect(0.0,-yside,xside,0.0,0.0,yside,-xside,0.0,fill);
}

/*-------------------------------*/

void drtri (double size, int fill)
{
  double side;

 side = size;
 rmove (0.0, side * 0.5775);
 rdrtri(-0.5*side,-0.866*side,    side,0.0,
		-0.5*side,0.866*side,fill);
/*
 rmove (0.0, side * 0.5775);
 rdraw (-0.5 * side, -0.866 * side);
 rdraw (side, 0.0);
 rdraw (-0.5 * side, 0.866 * side);
*/
} 

/*-------------------------------*/

void drtriu (double size, int fill)
{
  double side;

 side = size;
 rmove (0.0, -side * 0.5775);
 rdrtri(-0.5*side,0.866*side,    side,0.0,
		-0.5*side,-0.866*side,fill);
/*
 rmove (0.0, -side * 0.5775);
 rdraw (-0.5 * side, 0.866 * side);
 rdraw (side, 0.0);
 rdraw (-0.5 * side, -0.866 * side);
*/
} 

/*-------------------------------*/

void drast (double size)

/* asterisk */

{
  double side;

 side = size;
 rmove ( 0.25 * side, .433 * side);
 rdraw (-0.5 * side, -0.866 * side);
 rmove ( 0.5 * side,  0.0);
 rdraw (-0.5 * side,  0.866 * side);
 rmove (-0.25 * side, -.433 * side);
 rdraw ( side, 0.0);
} 

/*-------------------------------*/

void drchar (double size, int dchar)
{
  double side,oldside;

 side = size;
 oldside = cwidth(side,".");
 rmove (-.35*side,-.35*side);
 textf ("%c",(char)dchar);
 cwidth(oldside,".");
} 

/*-------------------------------*/

