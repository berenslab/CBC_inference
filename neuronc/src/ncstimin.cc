/* segment ncstimin in program nc */

/* Low-level routines for making stimuli. */

#include "ncsub.h"
#include "nconst.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif

#ifdef __cplusplus
}
#endif

/*------------------------------------*/

int incirc(double x, double y, double circx, double circy, double rad)

/* calculate whether point is inside a circle
   of given radius.
*/

{
   double tempx, tempy;
   int inside;

   tempx = circx - x;
   tempy = circy - y;
   inside = (sqrt(tempx*tempx + tempy*tempy) <= rad);
   return inside;
}

/*------------------------------------*/

int inrect(double x, double y, double maxx, double minx, 
	double maxy, double miny, double orient)
                                     
/* Calculate whether point is inside a rectangle
   of given limits. Rectangle is rotated about 
   its center by "orient" degrees.
*/

{
   int inside;
   double orientrad, coso, sino, xr, yr, xcent, ycent;

   if (maxx < minx) return 0;
   if (maxy < miny) return 0;
   xcent = (maxx + minx) * 0.5;
   ycent = (maxy + miny) * 0.5;
   orientrad = MPI / 180 * orient;	/* rotate test point opposite way */
   coso = cos(orientrad);
   sino = sin(orientrad);
   xr = (x-xcent) * coso - (y-ycent) * sino + xcent;
   yr = (x-xcent) * sino + (y-ycent) * coso + ycent;
   inside = ((minx <= xr) && (xr <= maxx) && (miny <= yr) && (yr <= maxy));
   return inside;
}

/*------------------------------------*/
