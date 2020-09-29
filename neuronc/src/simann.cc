/* simann.cc */

/* Simulated annealing */

#include <stdio.h>
#include <math.h>

#include "nc.h"
#include "ncsub.h"
#include "ncomp.h"

/* #include "nrutil.h" */

static char sim_ann_state[RNDSIZ];

double drand(void);
char *rsetstate(char *state);
void restorstate(void);
void execerror(char *s, char *t);
char *emalloc(unsigned int n);
void free(void *p);

double callfunc(Symbol *funcp, int npar, double par1, double par2,
                                         double par3, double par4);

double simfunc (double *parms);

double amotsa(double **p, double y[], double psum[], int ndim,
	double pb[], double *yb, double (*simfunc)(double *), int ihi, 
	double *yhi, double fac);


/*------------------------------------------------------------*/

#define GET_PSUM for (n=0; n<ndim; n++) { for (sum=0.0,m=0; m<mpts; m++) sum += p[m][n]; psum[n]=sum; }

extern long idum;		/* Defined and initialized in main. */
double tt;			/* Communicates with amotsa */

void amebsa(double **p, double y[], int ndim, double pb[], double *yb, 
	double ftol, double (*simfunc)(double *), int *iter, double temptr)

/*  Taken from Press et al. (1992) "Numerical Recipes in C", Chap 10, p 452 */

/* Multidimensional minimization of the function simfunc(x) where
x[0..ndim-1] is a vector in ndim dimensions, by simulated annealing
combined with the downhill simplex method of Nelder and Mead. The
input matrix p(0...ndim] [0..ndim-1] has ndim+l rows, each an
ndimdimensional vector which is a vertex of the starting simplex.
Also input are the following: the vector y[1...ndim+l], whose
components must be pre-initialized to the values of simfunc,
evaluated at the ndim+l vertices (rows) of p; ftol, the
fractional convergence tolerance to be achieved in the function
value for an early return; iter, and temptr. 

The routine makes iter function evaluations at an annealing
temperature temptr, then returns. You should then decrease temptr
according to your annealing schedule, reset iter, and call the
routine again (leaving other arguments unaltered between calls). 

If iter is returned with a positive value, then early convergence
and return occurred. If you initialize yb to a very large value
on the first call, then yb and pb[1. .ndim] will subsequently
return the best function value and point ever encountered (even
if it is no longer a point in the simplex). */

{

    int i,ihi,ilo,j,m,n,mpts=ndim+1;
    double rtol,sum,swap,yhi,ylo,ynhi,ysave,yt,ytry,*psum;

 psum= (double *)emalloc(sizeof(double)*ndim);
 tt = -temptr;
 GET_PSUM
 for (;;) {
	ilo=1;		/* Determine which point is the highest (worst), */
	ihi=2;		/*   next-highest, and lowest (best). */
 	rsetstate (sim_ann_state);
	ynhi=ylo=y[0]+tt*log(drand());   /* looking at vertex, it gets */
	yhi=y[1]+tt*log(drand());	     /* a random thermal fluctuation. */
	if (ylo > yhi) {
		ihi=0;
		ilo=1;
		ynhi=yhi;
		yhi=ylo;
		ylo=ynhi;
	}
	for (i=2; i<mpts; i++) { /* Loop over the points in the simplex. */
		yt=y[i]+tt*log(drand()); /* More thermal fluctuations. */
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
	restorstate();

/* Compute the fract range from highest to lowest and return if OK. */

	rtol=2.0*fabs(yhi-ylo)/(fabs(yhi)+fabs(ylo)); 

	if (rtol<ftol || *iter < 0) {	/* If returning, put best point */ 
		swap=y[0];			/* and value in first slot. */
		y[0] = y[ilo];
		y[ilo]=swap;
		for (n=0; n<ndim; n++) {
			swap      = p[0][n];
			p[0][n]   = p[ilo][n];
			p[ilo][n] = swap;
		}
	break;
	}

	*iter -= 2;

/* Begin a new iteration. 
   First extrapolate by a factor -1 through the face of the 
   simplex across from the high point, i.e., reflect the 
   simplex from the high point. */

	ytry = amotsa(p,y,psum,ndim,pb,yb,simfunc,ihi,&yhi,-1.0);
	if (ytry <= ylo) {

	/* Gives a result better than the best point, so try an */
	/* additional extrapolation by a factor of 2.

	  ytry = amotsa(p,y,psum,ndim,pb,yb,simfunc,ihi,&yhi,2.0);

	} else if (ytry >= ynhi) {

	/* The reflected point is worse than the second-highest, */
	/* so look for an intermediate lower point, i.e., do a  */
	/* one-dimensional contraction. */

	ysave=yhi;

	ytry=amotsa(p,y,psum,ndim,pb,yb,simfunc,ihi,&yhi,0.5);

	if (ytry >= ysave) {	/* Can't seem to get rid of that high point. */
 	   for (i=0; i<mpts; i++) {	/* Better contract around the lowest */
		if (i != ilo) {		/* (best) point. */
		   for	(j=0; j<ndim; j++) {
			psum[j]=0.5*(p[i][j] + p[ilo][j]);
			p[i][j]=psum[j];
		   }
		   y[i]=(*simfunc)(psum);
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
	     double pb[], double *yb, double (*simfunc)(double *), 
	     int ihi, double *yhi, double fac) 

/* Extrapolates by a factor fac through the face of the simplex */
/*  across from the high point, tries it, and replaces the high point */
/*  if the new point is better. */

{
    int j; 
    double fac1,fac2,yflu,ytry,*ptry;

 ptry= (double *)emalloc(sizeof(double)*ndim);
 fac1=(1.0-fac)/ndim;
 fac2=fac1-fac;
 for (j=0; j<ndim; j++)
	ptry[j] = psum[j]*fac1-p[ihi][j]*fac2;
 ytry=(*simfunc)(ptry);
 if (ytry <= *yb) {	/* Save the best-ever. */
	for (j=0; j<=ndim; j++) pb[j]=ptry[j];
	*yb=ytry;
 }
 rsetstate (sim_ann_state);
 yflu=ytry-tt*log(drand());	/* We added a thermal flucts to all the cur */
 restorstate();
 if (yflu < *yhi) {		/* vertices, but we subtr it here, to give */
	y[ihi]=ytry;		/* the simplex a thermal Brownian motion: It */
	*yhi=yflu;		/* likes to accept any suggested change. */
	for (j=0;j<ndim; j++) {
	    psum[j]   +=ptry[j] - p[ihi][j];
 	    p[ihi][j] = ptry[j];
	}
  }
 free(ptry);
 return yflu;
}

/*------------------------------------------------------------*/

static Symbol *spap;		/* simulation parameter array */
static int spsize;
static Symbol *simfp;
extern Inst *pc;

void simann(void) 

/* Simulated annealing routine to be called from the interpreter.
This routine is called from the "simann" statement, which has 2
arguments: 1) the simulation function and 2) an array of free
parameters.  The array has 2 dimensions:

    simann simfunc params;

  The first dimension of "params" is always 3, and the second
dimension is the number of free parameters:   

     params[3][nparm]

The first set of parameters are the starting values with which
the simulation is to be run, and the second and third are the
maximum and minimum values.  The starting values get modified and
the "simfunc" function (below) is called with a copy of them in a
1-D array.  This routine in turn copies these parameters back
into the original array (in place of the starting values) and
calls the simulation routine specified on the "simann" statement.  
*/
 
{
     int i;
     array *arrp;
     double **p,*params,*y;

#define EPSILON 1.01

  spap  = (Symbol *)(*pc++);	/* pointer to simulation parameter array */
  simfp = (Symbol *)(*pc++);	/* pointer to simulation function  */
  arrp = spap->arrp;		/* get pointer to address of array */
  if (!arrp) {
     execerror("simann: array space not allocated for",spap->name);
     return;
  }
  spsize = arrp->dim[1];		/* get first dim */
  if (arrp->dim[0] != 3)
    execerror ("simarr: params array doesn't contain 3 values", spap->name);

  if ((p=(double**)emalloc(sizeof (double)*(spsize+1)*spsize))==(double *)NULL){
    execerror ("simarr: can't allocate vertex array","");
  }

  if ((y=(double*)emalloc(sizeof (double)*(spsize+1)))==(double *)NULL){
    execerror ("simarr: can't allocate initial value array","");
  }

   /* first, calculate the initial vertices of the simplex */

  params = arrp->arr;
  for (i=0; i<=spsize; i++) {		/* ndim+1 rows */
    for (j=0; j<spsize; j++) {		/* ndim cols */
     if (i==j) p[i][j] = params[j] * EPSILON;
     else      p[i][j] = params[j];
  }
  for (i=0; i<=spsize; i++) {		/* ndim+1 rows */
    y[i] = simfunc (p[i]); 
  }
  amebsa(p, y, spsize, pb[], *yb, 
	ftol, simfunc, iter, temptr)
}

/*------------------------------------------------------------*/

double simfunc (double *parms)

{
    int i;
    datum d1;
    double *params;

  if (initval=emalloc (sizeof(double) * spsize)==NULL) {
    execerror ("simfunc: can't allocate param array","");
  }
  params = spap->arrp->arr;
  for (i=0; i<spsize; i++) {
    params[i] = parms[i];
  } 
  return (callfunc(simfp,0, 0.0, 0.0, 0.0, 0.0));
}

