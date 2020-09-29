/* Segment NCROT in program NC */

/* sets up rotation and translation matrices */

#include "nc.h"
#include "ncsub.h"
#include "ncomp.h"

extern double xyscal,zscal;

double ncmat[4][4] =  {0};
double ncmatc[4][4] =  {0};

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

void idmat(double (*mat)[4]);
void revolv (double (*obmat)[4], double theta, double xt, double yt, double zt,
		 double xh, double yh, double zh);
void matmul (double (*m1)[4], double (*m2)[4]);


/****************************************/

void msetup(double xrot, double yrot, double zrot, 
		double centx, double centy, double centz, 
		double tranx, double trany, double tranz, 
		double mag, double xzero, double yzero, double zzero, 
		double scale, double (*mat)[4])
                                       
/* Routine to set up the transformation matrix,
given the rotation, translation, and revolution params.

Parameters needed:
  
	xrot,yrot,zrot:  rotation angle around the axis specified.
	centx,centy,centz: coords of the point of rotation.
	tranx,trany: 	    translation after scaling in the x-y plane.
	xzero,yzero,zzero:  translation before rotation 
				(loc of zero for the axes)
	scale :     preliminary scaling factor
        mag:        user set magnification
*/


{
  int i;
  double x1,x2,y1,y2,z1,z2;
  static double obmat[4][4];
  
  x1 = 1; x2 = 10;			/* the three axes */
  y1 = 1; y2 = 10;
  z1 = 1; z2 = 10;

/*  idmat (mat); 				/* make identity matrix */

  idmat (obmat);
  obmat [3][0] = xzero;			/* add translations to matrix */
  obmat [3][1] = yzero;
  obmat [3][2] = zzero;
  matmul (mat,obmat); 

/*  idmat (obmat);			/* add scale factor */
/*  for (i=0; i<3; i++)
   obmat [i][i] = scale;
  matmul (mat,obmat);  */

  revolv (mat,xrot,x1,centy,centz,x2,centy,centz); 	/* set up xrot */
  revolv (mat,yrot,centx,y1,centz,centx,y2,centz);	/* set up yrot */
  revolv (mat,zrot,centx,centy,z1,centx,centy,z2);	/* set up zrot */

  idmat (obmat);			/* add magnification factor */
  for (i=0; i<3; i++)
   obmat [i][i] = mag*scale;
  matmul (mat,obmat);

  idmat (obmat);
  obmat [3][0] = tranx;			/* add translations to matrix */
  obmat [3][1] = trany;
  obmat [3][2] = tranz;
  matmul (mat,obmat); 
}

/****************************************/

void idmat(double (*mat)[4])

/* set up the identity matrix */

{
  int i,j;

  for (i=0; i<4; i++)			/* Set up identity matrix */
   {
    for (j=0; j<4; j++)
     mat [i][j] = 0.0;
    mat [i][i] = 1.0; 
   } 
}

/****************************************/

void revolv (double (*obmat)[4], double theta, double xt, double yt, double zt,
		 double xh, double yh, double zh)

/* Matrix rotation routine */

/* Courtesy of Neil Friedman and Bill Wood */
/* Moore Sch. of Eng. at U. of PA          */

/* This routine sets up a transformation matrix for 3d rotations,
translations, and scalings. */

{
    int i,j;
    static double mat[4][4];
    double theta2,cosine,sine;
    double a,b,c,L;
    double asin,bsin,csin;
    double asq,bsq,csq;
    double abo,aco,bco;

 for (i=0; i<4; i++)
  for (j=0; j<4; j++)
   mat [i][j] = 0.0;
 mat [3][3] = 1.0;
 a = xh-xt;
 b = yh-yt;
 c = zh-zt;
 L = sqrt (a*a + b*b + c*c);
 if (L < 0.0001) return;
 a /= L;
 b /= L;
 c /= L;
 theta2 = theta * 3.1415926 / 180.0;
 cosine = cos (theta2);
 sine   = -sin (theta2);
 asin = a * sine;
 bsin = b * sine;
 csin = c * sine;
 asq = a*a;
 bsq = b*b;
 csq = c*c;
 abo = a*b*(1-cosine);
 aco = a*c*(1-cosine);
 bco = b*c*(1-cosine);
 mat [0][0] = asq+(1-asq)*cosine;
 mat [0][1] = abo-csin;
 mat [0][2] = aco+bsin;
 mat [1][0] = abo+csin;
 mat [1][1] = bsq+(1-bsq)*cosine;
 mat [1][2] = bco-asin;
 mat [2][0] = aco-bsin;
 mat [2][1] = bco+asin;
 mat [2][2] = csq+(1-csq)*cosine;
 mat [3][0] = xt-xt*mat[0][0]-yt*mat[1][0]-zt*mat[2][0];
 mat [3][1] = yt-xt*mat[0][1]-yt*mat[1][1]-zt*mat[2][1];
 mat [3][2] = zt-xt*mat[0][2]-yt*mat[1][2]-zt*mat[2][2];
 matmul (obmat,mat);
}

/****************************************/

void matmul (double (*m1)[4], double (*m2)[4])
                              
/* multiply 4 x 4 matrices m1 and m2;
leave the result in m1. */

{
    int i,j;
    double t3[4];
  
 for (i=0; i<4; i++)
  {
   for (j=0; j<4; j++)
    {
     t3[j] = m1[i][0] * m2[0][j] + m1[i][1] * m2[1][j] + m1[i][2] * m2[2][j];
     t3[j] = t3[j] + m1[i][3] * m2[3][j];
    }
   for (j=0; j<4; j++)
    m1[i][j] = t3[j];
  }
}

/******************************************/

void dmat (double (*xmat)[4])
{
   int i,j;

   for (i=0; i<4; i++)
    {
     for (j=0; j<4; j++)
      ncfprintf (stderr,"%g ",xmat [i][j]);
     ncfprintf (stderr,"\n");
    }
} 

/****************************************/

void setmat (double (*mat)[4])

#define PSCALE 200.0

{
   double xrot,yrot,zrot;
   double xtran,ytran,mag;
   double xcent,ycent,zcent;

  idmat (mat);

/*  msetup (arr[XROT], arr[YROT], arr[ZROT],
	  0.0, 0.0, 0.0,
	  arr[XTRAN], arr[YTRAN],0.0,
	  arr[MAG],
	  -arr[XCENT]*xyscal, -arr[YCENT]*xyscal, -arr[ZCENT]*zscal,
	  (1.0/PSCALE), mat);
*/

xrot = 0.0;
yrot = 0.0;
zrot = 0.0;
xtran = 0.0;
ytran = 0.0;
xcent= 0.0;
ycent= 0.0;
zcent= 0;
mag = 1.0;
xyscal = 1.0 / PSCALE;

msetup (xrot, yrot, zrot,
        0.0,  0.0,  0.0,
        xtran, ytran, 0.0,
        mag,
	xcent, ycent, zcent,
	xyscal, mat);
}

/****************************************/

void setrot(double xrot, double yrot, double zrot, 
		double xcent, double ycent, double zcent, 
		double rxcent, double rycent, double rzcent, double scal)
{
   double xtran,ytran,mag,(*mat)[4];

if (scal==0.0) scal = .001;
mat = ncmat;

xtran = 0.5;
ytran = 0.5;
if (xcent==LARGENODE && ycent==LARGENODE && zcent==LARGENODE) { 
  xcent = -scal/2.0;
  ycent = -scal/2.0;
  zcent = 0.0;
}
xcent += rxcent;
ycent += rycent;
zcent += rzcent;

mag = 1;			/* set mag = .75, vid -m 1.333 to fill screen*/
xyscal = 1.0/scal;

idmat (mat);			/* set matrix to no rotation, trans or scale */

msetup (xrot, yrot, zrot,
        0.0,  0.0,  0.0,
        xtran, ytran, 0.0,
        mag,
	xcent, ycent, zcent,
	xyscal, mat);
}

/****************************************/

void setray (double xrot, double yrot, double zrot, 
	double xcent, double ycent, double zcent,
	double rxcent, double rycent, double rzcent,
	double scale)

{
idmat (ncmatc);			/* set matrix to no rotation, trans or scale */

msetup (xrot, yrot, zrot,
        0.0,  0.0,  0.0,
        0.0, 0.0, 0.0,
        1.0,
	xcent+rxcent, ycent+rycent, zcent+rzcent,
	scale, ncmatc);
}

/****************************************/
 
void setinv (double (*mat)[4])
{
  idmat (mat);

/*
  msetinv (arr[XROT], arr[YROT], arr[ZROT],
	  0.0, 0.0, 0.0,
	  arr[XTRAN], arr[YTRAN],0.0,
	  arr[MAG],
	  -arr[XCENT]*xyscal, -arr[YCENT]*xyscal, -arr[ZCENT]*zscal,
	  (1.0/PSCALE), mat);
*/
}

/****************************************/

void msetinv(double xrot, double yrot, double zrot, 
	double centx, double centy, double centz, 
	double tranx, double trany, double tranz, 
	double mag, double xzero, double yzero, double zzero, 
	double scale, double (*mat)[4])

/* routine to set up a transformation matrix
   inverse to the matrix created by msetup above. 
*/

{
  int i;
  double x1,x2,y1,y2,z1,z2;
  double obmat[4][4];
  
  x1 = 1; x2 = 10;			/* the three axes */
  y1 = 1; y2 = 10;
  z1 = 1; z2 = 10;

  idmat (mat); 				/* make identity matrix */

  idmat (obmat);
  obmat [3][0] = -tranx;			/* add translations to matrix */
  obmat [3][1] = -trany;
  obmat [3][2] = -tranz;
  matmul (mat,obmat); 
  
  if (mag==0) mag = 1;
  idmat (obmat);			/* add magnification factor */
  for (i=0; i<3; i++)
   obmat [i][i] = 1/mag;
  matmul (mat,obmat);

  revolv (mat,-zrot,centx,centy,z1,centx,centy,z2);	/* set up zrot */
  revolv (mat,-yrot,centx,y1,centz,centx,y2,centz);	/* set up yrot */
  revolv (mat,-xrot,x1,centy,centz,x2,centy,centz); 	/* set up xrot */

  if (scale == 0) scale = 1;
  idmat (obmat);				/* add scale factor */
  for (i=0; i<3; i++)
   obmat [i][i] = 1/scale;
  matmul (mat,obmat);

}

/****************************************/

void transf (double xp, double yp, double zp, double *xv, double *yv, 
	double *zv, double (*tmat)[4])
	                
/* Transform the point xp,yp,zp by matrix
"tmat" to give point xv,yv,zv.
Scaling, rotation, translation, and revolution
are all accomplished in the same step.  */

{
 *xv = xp*tmat[0][0]+yp*tmat[1][0]+zp*tmat[2][0]+tmat[3][0];
 *yv = xp*tmat[0][1]+yp*tmat[1][1]+zp*tmat[2][1]+tmat[3][1];
 *zv = xp*tmat[0][2]+yp*tmat[1][2]+zp*tmat[2][2]+tmat[3][2];
}

/*******************************/
