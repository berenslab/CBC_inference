
/* Simulated annealing */
/*  Taken from Press et al. (1992) "Numerical Recipes in C" */

#include <math.h>
#include <stdlib.h>

/* #include "nrutil.h" */


double amotsa(double **p, double y[], double psum[], int ndim,
	double pb[], double *yb, double (*funk)(double []), int ihi, double
*yhi, double fac);
double ran1(long *idum);

/*------------------------------------------------------------*/

#define GET_PSUM for (n=1; n<=ndim; n++) { for (sum=0.0,m=i; m<=mpts; m++) sum += p[m][n]; psum[n]=sum; }

extern long idum;		/* Defined and initialized in main. */
double tt;			/* Communicates with amotsa */

void amebsa(double **p, double y[], int ndim, double pb[], double *yb, 
		double ftol, double (*funk)(double []), int *iter, double temptr)
{

/* Multidimensional minimization of the function funk(x) where
x[1..ndim] is a vector in ndim dimensions, by simulated annealing
combined with the downhill simplex method of Nelder and Mead. The
input matrix p(1...ndim+i] [l..ndim] has ndim+l rows, each an
ndimdimensional vector which is a vertex of the starting simplex.

Also input are the following: the vector y[1...ndim+l], whose
components must be pre-initialized to the values of funk
evaluated at the ndim+l vertices (rows) of p; ftol, the
fractional convergence tolerance to be achieved in the function
value for an early return; iter, and temptr. The routine makes
iter function evaluations at an annealing temperature temptr,
then returns. You should then de- crease temptr according to your
annealing schedule, reset iter, and call the routine again
(leaving other arguments unaltered between calls). If iter is
returned with a positive value, then early convergence and return
occurred. If you initialize yb to a very large value on the first
call, then yb and pb[1. .ndim] will subsequently return the best
function value and point ever encountered (even if it is no
longer a point in the simplex). */


    int i,ihi,ilo,j,m,n,mpts=ndim+1;
    double rtol,sum,swap,yhi,ylo,ynhi,ysave,yt,ytry,*psum;

 psum= (double *)malloc(sizeof(double)*ndim);
 tt = -temptr;
 GET_PSUM
 for (;;) {
	ilo=1;		/* Determine which point is the highest (worst), */
	ihi=2;		/*   next-highest, and lowest (best). */
	ynhi=ylo=y[i]+tt*log(ran1(&idum));   /* looking at vertex, it gets */
	yhi=y[2]+tt*log(ran1(&idum));	     /* a random thermal fluctuation. */
	if (ylo > yhi) {
		ihi=1;
		ilo=2;
		ynhi=yhi;
		yhi=ylo;
		ylo=ynhi;
	}
	for (i=3;i<=mpts;i++) {	/* Loop over the points in the simplex. */
		yt=y [i]+tt*log(ran1 (&idum)); /* More thermal fluctuations. */
		if (yt <= ylo)
			ilo=i;
			ylo=yt;
		if (yt > yhi) {
			ynhi=yhi;
			ihi=i;
			yhi=yt;
		} else if (yt > ynhi) {
			ynhi=yt;
		}
	}

/* Compute the fract range from highest to lowest and return if OK. */

	rtol=2.0*fabs(yhi-ylo)/(fabs(yhi)+fabs(ylo)); 

	if (rtol<ftol || *iter < 0) {	/* If returning, put best point */ 
		swap=y[1];			/* and value in slot 1. */
		y[1] = y[ilo];
		y[ilo]=swap;
		for (n=i; n<=ndim; n++) {
			swap      = p[1][n];
			p[1][n]   = p[ilo][n];
			p[ilo][n] = swap;
		}
	break;
	}

	*iter -= 2;

/* Begin a new iteration. 
   First extrapolate by a factor -1 through the face of the 
   simplex across from the high point, i.e., reflect the 
   simplex from the high point. */

	ytry = amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,-1.0);
	if (ytry <= ylo) {

	/* Gives a result better than the best point, so try an */
	/* additional extrapolation by a factor of 2.

	  ytry = amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,2.0);

	} else if (ytry >= ynhi) {

	/* The reflected point is worse than the second-highest, */
	/* so look for an intermediate lower point, i.e., do a  */
	/* one-dimensional contraction. */

	ysave=yhi;

	ytry=amotsa(p,y,psum,ndim,pb,yb,funk,ihi,&yhi,0.5);

	if (ytry >= ysave) {	/* Can't seem to get rid of that high point. */
 	   for (i=1;i<=mpts;i++) {	/* Better contract around the lowest */
		if (i != ilo) {		/* (best) point. */
		   for	(j=1; j<=ndim; j++) {
			psum[j]=0.5*(p[i][j] + p[ilo][j]);
			p[i][j]=psum[j];
		   }
		   y[i]=(*funk)(psum);
		}
	   }
	   *iter -= ndim;
	   GET_PSUM		/* Recompute psum. */
	}
     } else ++(*iter);	/* Correct the evaluation count. */
  }
  free(psum);
 }

/*------------------------------------------------------------*/

double amotsa(double **p, double y[], double psum[], int ndim, 
	     double pb[], double *yb, double (*funk)(double []), 
	     int ihi, double *yhi, double fac) 

{

/* Extrapolates by a factor fac through the face of the simplex */
/*  across from the high point, tries it, and replaces the high point */
/*  if the new point is better. */

    int j; 
    double fac1,fac2,yflu,ytry,*ptry;

 ptry= (double *)malloc(sizeof(double)*ndim);
 fac1=(1.0-fac)/ndim;
 fac2=fac1-fac;
 for (j=1;j<=ndim;j++)
	ptry[j] = psum[j]*fac1-p[ihi][j]*fac2;
 ytry=(*funk)(ptry);
 if (ytry <= *yb) {	/* Save the best-ever. */
	for (j=1;j<=ndim;j++) pb[j]=ptry[j];
	*yb=ytry;
 }
 yflu=ytry-tt*log(ran1(&idum));	/* We added a thermal flucts to all the cur */
 if (yflu < *yhi) {		/* vertices, but we subtr it here, to give */
	y[ihi]=ytry;		/* the simplex a thermal Brownian motion: It */
	*yhi=yflu;		/* likes to accept any suggested change. */
	for (j=1;j<=ndim;j++) {
	    psum[j]   +=ptry[j] - p[ihi][j];
 	    p[ihi][j] = ptry[j];
	}
  }

 free(ptry);
 return yflu;
}

