/* line.cc */

/* Computes distance between and points of closest approach between
   two lines in 3 dimensions.
*/

#include "vec.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif

#ifdef __cplusplus
}
#endif

#include "ncio.h"

/*-------------------------------------------*/
vec vcross(vec a, vec b)

/* cross product of vectors a, b */

{
   vec c;

  c.x = a.y * b.z - a.z * b.y;
  c.y = a.z * b.x - a.x * b.z;
  c.z = a.x * b.y - a.y * b.x;
  return c;
}

/*-------------------------------------------*/

double vlen (vec a)

{
  return (sqrt(a.x*a.x + a.y*a.y + a.z*a.z));
}

/*-------------------------------------------*/

double vdist (vec a, vec b)

/* Return dist squared between 2 points (faster than sqrt). */

{
     double x,y,z;

  x = b.x - a.x;
  y = b.y - a.y;
  z = b.z - a.z;
  return (x*x + y*y + z*z);
}

/*-------------------------------------------*/
double vdot(vec a, vec b)

/* dot product of vectors a, b */

{
   return (a.x*b.x + a.y*b.y + a.z*b.z);
}

/*-------------------------------------------*/
vec vadd(vec a, vec b)

/* add vectors a, b */

{
   vec c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  return c;
}
/*-------------------------------------------*/
vec vsub(vec a, vec b)

/* subtract vectors a, b */

{
   vec c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  return c;
}
/*-------------------------------------------*/
void vcopy(vec a, vec b)

/* subtract vectors a, b */

{
   b.x = a.x;
   b.y = a.y;
   b.z = a.z;
}
/*-------------------------------------------*/
void vprint(vec a)

/* print a vector */

{
   vec c;

  ncfprintf (stderr,"%g %g %g\n",a.x,a.y,a.z);
}
/********************************************/

double norm_line (vec a, vec b, vec p, vec q, double *vt, double *vs, 
					vec *lp1, vec *lp2)

/* Compute distance and points of closest approach for lines described by */
/*  a + bt,  p + qs, where a,b,p,q are vectors and s,t are line variables */
/* Store answer in lp1, lp2, return distance */
 
{
   double denom,dist,k,u,t,s;
   double a1,a2,a3;
   double b1,b2,b3;
   double c1,c2,c3;
   double p1,p2,p3;
   double q1,q2,q3;
   vec c;
 

a1 = a.x; a2 = a.y; a3 = a.z;
b1 = b.x; b2 = b.y; b3 = b.z;
p1 = p.x; p2 = p.y; p3 = p.z;
q1 = q.x; q2 = q.y; q3 = q.z;

c = vcross(b,q);
c1 = c.x;
c2 = c.y;
c3 = c.z;

/* vprint (b); */
/* vprint (q); */

if (vlen(c)==0) {		/* if lines are parallel, return neg dist */
   dist = -1;
   return dist;
}

dist = 1 / vlen(c) * vdot(vsub(p,a),c);
if (dist < 0) dist = -dist;

/* ncfprintf (stderr,"dist %g %g %g %g\n",dist,c.x,c.y,c.z); */

/* This solves:  (a + bt) - (p + qs) = kc */
/*  where a,b,p,q are vectors and c is the cross product of b and q. */
/* Solved in Maple */

denom = -c3*q1*b2+q1*c2*b3+c1*q3*b2+b1*c3*q2-b1*c2*q3-c1*b3*q2;

t = -(-c3*q1*a2+c3*q1*p2+q1*c2*a3-q1*c2*p3+c1*q3*a2-c1*q3*p2+a1*c3*q2-
		a1*c2*q3-p1*c3*q2+p1*c2*q3-c1*a3*q2+c1*p3*q2)/ denom;

s = (-c2*a3*b1+c2*p3*b1+c2*b3*a1-c2*b3*p1-p3*c1*b2+b3*c1*p2+a3*c1*b2+
		a2*c3*b1-b3*c1*a2+c3*b2*p1-p2*c3*b1-c3*b2*a1)/ denom;

/* k = (-a3*q1*b2+a3*b1*q2+b3*q1*a2-b3*q1*p2-b3*a1*q2+b3*p1*q2+p3*q1*b2-
		p3*b1*q2-q3*a2*b1-q3*b2*p1+q3*p2*b1+q3*b2*a1)/ denom; */


lp1->x = a1 + b1 * t;
lp1->y = a2 + b2 * t;
lp1->z = a3 + b3 * t;

lp2->x = p1 + q1 * s;
lp2->y = p2 + q2 * s;
lp2->z = p3 + q3 * s;

*vt = t;
*vs = s;

/* ncfprintf (stderr,"t %g s %g %g\n",t,s,denom); */
/* printf (stderr,"%g %g %g   %g %g %g\n",lp1->x,lp1->y,lp1->z,lp2->x,lp2->y,lp2->z);*/

return dist;

}

/* main()

{
   double dist;
   vec a,b,p,q,t,s,lp1,lp2;
*/
/*
  These vectors describe 2 lines (a + bt), (p + qs). Their
  closest points are sqrt(14) apart.  The points of closest
  approach are: (-1,2,1) and (2,1,3).  Use this to check math
  in norm_line above.
*/

/*
a.x = 0;
a.y = 3;
a.z = 0;

b.x = 1;
b.y = 1;
b.z = -1;

p.x = 5;
p.y = 8;
p.z = 2;

q.x = 3;
q.y = 7;
q.z = -1;

 dist = norm_line (a,b,p,q,&t,&s,&lp1,&lp2);
}
*/
