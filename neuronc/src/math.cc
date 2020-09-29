/* math.c for "nc" */

#include <errno.h>
#include <string.h>
#include <ctype.h>
#include "nc.h"
#include "ncsub.h"
#include "nconst.h"
#include "y.tab.h"
#include "control.h"

#define MPI 3.14159265358979323846264

extern int errno;
double errcheck(double d, const char *s);

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif
/*
char *strtok(char *s, const char *delim);
int strcmp(const char *s1, const char *s2);
char *strstr(const char *haystak, const char *needle);
*/
double atof(char *s);

#ifdef __cplusplus
}
#endif

#include "ncio.h"

#define max(x, y)       (((x) < (y)) ? (y) : (x))
#define min(x, y)       (((x) < (y)) ? (x) : (y))

double drand();
void execerror(const char *s, const char *t);
char *emalloc(unsigned int n);
void fftstp (double *zinr, double *zini, int after, int now, 
		int before, double *zoutr, double *zouti);
void four1(double *data, int nn, int isign);
int poisdev(double xm);
int binomdev(double x, int n);
double gamdev(double x);
double gasdev();
void setrseed();
void efree(void *ptr);
datum copyarr(datum *src, double **p, int *siz);
void erarr (ncarray *arrp);
double getvarval(const char *name);
double fread_interp(char *str);   /* intepret expressions by fread */
int getarr(datum *d, double **parr, int *siz);

datum domonofunc(datum &x, double (*fv) (double p) );
datum doduofunc(datum &x, datum &y, double (*fv) (double p1, double p2) );

/*-------------------------------------------------------*/

double dolog(double x)
{
	x = errcheck(log(x), "log");
	return x;
}

datum nclog(datum &x)
{
        return domonofunc(x,dolog);
}

double dolog10(double x)
{
	x=errcheck(log10(x), "log10");
	return x;
}
datum nclog10(datum &x)
{
        return domonofunc(x,dolog10);
}

/*-------------------------------------------------------*/

double dosqrt(double x)
{
	x=errcheck(sqrt(x), "sqrt");
	return x;
}

datum ncsqrt(datum &x)
{
        return domonofunc(x,dosqrt);
}

double doexp(double x)
{
	if (x < -745.0) x = 0;
	else x=errcheck(exp(x), "exp");
	return x;
}
datum ncexp(datum &x)
{
        return domonofunc(x,doexp);
}

double doexp10(double x)
{
	if (x < -320.0) x = 0.;
	else x=errcheck(exp(x*LN10), "exp10");
	return x;
}
datum ncexp10(datum &x)
{
        return domonofunc(x,doexp10);
}

/*---------------------------------------------------*/

/*
double doroundup(double x, double precision)
{
	double expon;

	expon = floor(log10(x));
	// fprintf (stderr,"doroundup %.20g %.20g\n",nearbyint(x*exp10(-expon+precision)+0.5), exp10(-precision));
	return floor(x*exp10(-expon+precision) + 1) * exp10(expon-precision);
	 // return int(x * 1e10 + 1) * 1e-10;
}

datum ncroundup(datum &x, datum &precision)
{
        return doduofunc(x,precision,doroundup);
}
/* */

/*---------------------------------------------------*/

double dosin(double x)
{
	x=errcheck(sin(x), "sine");
	return x;
}
datum ncsin(datum &x)
{
        return domonofunc(x,dosin);
}

double doasin(double x)
{
	x=errcheck(asin(x), "asine");
	return x;
}
datum ncasin(datum &x)
{
        return domonofunc(x,doasin);
}
/*---------------------------------------------------*/

double docos(double x)
{
	x=errcheck(cos(x), "cos");
	return x;
}
datum nccos(datum &x)
{
        return domonofunc(x,docos);
}
double doacos(double x)
{
	x=errcheck(acos(x), "acos");
	return x;
}
datum ncacos(datum &x)
{
        return domonofunc(x,doacos);
}
double dotan(double x)
{
	x=errcheck(tan(x), "tan");
	return x;
}
datum nctan(datum &x)
{
        return domonofunc(x,dotan);
}
double doatan(double x)
{
	x=errcheck(atan(x), "atan");
	return x;
}
datum ncatan(datum &x)
{
        return domonofunc(x,doatan);
}
datum ncatan2(datum &x, datum &y)
{
	x.val=errcheck(atan2(y.val,x.val), "atan2");
	x.vtype = NUMBER;
	return x;
}

double dotanh(double x)
{
	x=errcheck(tanh(x), "tanh");
	return x;
}
datum nctanh(datum &x)
{
        return domonofunc(x,dotanh);
}

double dosinh(double x)
{
	x=errcheck(sinh(x), "sinh");
	return x;
}
datum ncsinh(datum &x)
{
        return domonofunc(x,dosinh);
}

double dopow(double x, double y)
{
	x=errcheck(pow(x,y), "pow");
	return x;
}

datum ncpow(datum &x, datum &y)
{
        return doduofunc(x,y,dopow);
}

double xgauss (double &x, double &r)

{
   double val;

  val = -(x*x/(r*r));
  if (val < -745.0) val = 0;
  else val = exp(val);
  return val;
}

datum ncgaus(datum &x, datum &r)
{
	x.val=errcheck(xgauss(x.val,r.val), "gaussian");
	x.vtype = NUMBER;
	return x;
}

datum ncpoisdev(datum &x)
{
	setrseed();
	//x.val=errcheck(((double)poisdev(x.val)), "poisson dev");
	x.val=poisdev(x.val);
	x.vtype = NUMBER;
	return x;
}

datum ncgamdev(datum &x)
{
	setrseed();
	//x.val=errcheck(((double)poisdev(x.val)), "poisson dev");
	x.val=gamdev(x.val);
	x.vtype = NUMBER;
	return x;
}

datum ncbinomdev(datum &x, datum &n)
{
	setrseed();
	//x.val=errcheck((binomdev(x.val,int(n.val))), "binom dev");
	x.val=binomdev(x.val,int(n.val));
	x.vtype = NUMBER;
	return x;
}

datum ncgasdev(void)
{
	datum d;

	setrseed();
	d.val=errcheck(((double)gasdev()), "gaussian dev");
	d.vtype = NUMBER;
	return d;
}

datum ncrgasdev(datum &i)
{
	datum d;
        double rgasdev (int ngen);

	setrseed();
	d.val=errcheck(((double)rgasdev((int)i.val)), "rgaussian dev");
	d.vtype = NUMBER;
	return d;
}

datum ncrand(void)
{
	datum d;

	setrseed();
	d.val=errcheck(drand(), "random number");
	d.vtype = NUMBER;
	return d;
}

datum ncrrand(datum &i)
{
	datum d;
        double rrand (int i);

        
	setrseed();
	d.val=errcheck(rrand((int)i.val), "rrandom number");
	d.vtype = NUMBER;
	return d;
}


/*------------------------------------------------------------*/

datum domonofunc(datum &x, double (*fv) (double p) )
{
	int i;
        datum d2;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"domonofunc\n");
#endif
        if (x.type == ARRAY) {
                double *p1,*p2;
                int siz1,siz2;
          getarr(&x, &p1, &siz1);
          d2 = copyarr(&x, &p2, &siz2);
          if (x.vtype==NUMBER || x.vtype==LITCHAR) {
             for (i=0; i<siz1; i++,p1++,p2++) {
                 *p2 = fv(*p1);
             }
             d2.type = x.type;
             d2.vtype = x.vtype;
	  }
 	  return d2;
        }
        else {
          switch (x.vtype) {
            case NUMBER:
            case LITCHAR:
                  x.val = fv(x.val);
                  break;
            case STRING:
                  break;
          }
          x.type = VAR;
          x.vtype = NUMBER;
 	  return x;
        }
}

/*------------------------------------------------------------*/

datum domonofuncs(datum &x, double (*fv) (char *p) )
{
	int i;
        datum d2;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"domonofunc\n");
#endif
        if (x.type == ARRAY) {
                double *p1,*p2;
                int siz1,siz2;
          getarr(&x, &p1, &siz1);
          d2 = copyarr(&x, &p2, &siz2);
          if (x.vtype==STRING) {
                char **s1;
             s1 = (char **)p1;
             for (i=0; i<siz1; i++,s1++,p2++) {
                 *p2 = fv(*s1);
             }
             d2.type = x.type;
             d2.vtype = NUMBER;
          }
 	  return d2;
        }
        else {
          switch (x.vtype) {
            case NUMBER:
            case LITCHAR:
                  break;
            case STRING:
                  x.val = fv((char *)x.str);
                  break;
          }
          x.type = VAR;
          x.vtype = NUMBER;
 	  return x;
        }
}

/*------------------------------------------------------------*/

datum doduofunc(datum &x, datum &y, double (*fv) (double p1, double p2) )
{
	int i;
        int siz1,siz2,siz3,size;
        datum t;
        datum d3;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"doduofunc\n");
#endif
        if (x.type != ARRAY && y.type ==ARRAY) {
	   t = x; x = y; y = t;
	}
        if (x.type == ARRAY && y.type == ARRAY) {
                double *p1,*p2,*p3;
          getarr(&x, &p1, &siz1);
          getarr(&y, &p2, &siz2);
	  size = min(siz1,siz2);
          d3 = copyarr(&x, &p3, &siz3);
          if (x.vtype!=STRING)
             for (i=0; i<size; i++,p1++,p2++,p3++) {
                 *p3 = fv(*p1,*p2);
             }
          d3.type = x.type;
          d3.vtype = x.vtype;
 	  return d3;
	}
        if (x.type == ARRAY && y.type != ARRAY) {
                double *p1,*p2;
          getarr(&x, &p1, &siz1);
          d3 = copyarr(&x, &p2, &siz2);
          if (x.vtype!=STRING)
             for (i=0; i<siz1; i++,p1++,p2++) {
                 *p2 = fv(*p1,y.val);
             }
          d3.type = x.type;
          d3.vtype = x.vtype;
 	  return d3;
        }
        else {
          switch (x.vtype) {
            case NUMBER:
            case LITCHAR:
                  x.val = fv(x.val,y.val);
                  break;
            case STRING:
                  break;
          }
          x.type = VAR;
          x.vtype = NUMBER;
 	  return x;
        }
}

/*------------------------------------------------------------*/

datum ncinteger(datum &x)

{
   return domonofunc(x,floor);
}

/*------------------------------------------------------------*/

double ncabs(double x)
{
	return (double)((x)<0 ? -(x) : (x));
}

datum ncfabs(datum &a)
{
	return domonofunc(a,ncabs);
}

/*------------------------------------------------------------*/

double fzero (double (*func)(double), double near_val)

/* Finds x-value at a zero of a function of a single variable, */
/*  similar to fzero in matlab. */

{
    double x0, x1, x2;
    double crit = 1e-6;

 for (x0=near_val,x1=x0+1; fabs(x1-x0) > crit; ) {
      // printf ("      x0 %g %g\n",x0,x1);
      x2 = x1 - func(x1) * (x1 - x0) / (func(x1) - func(x0));
      x0 = x1; x1 = x2;
 }
 // printf ("final x0 %g %g\n",x0,x1);
 return x0;
}

/****************************************/

double domax(double a, double b)
{
   return (max(a,b));
}

double domin(double a, double b)
{
   return (min(a,b));
}

datum ncmax(datum &a, datum &b)
{
        return doduofunc(a,b,domax);
}

datum ncmin(datum &a, datum &b)
{
        return doduofunc(a,b,domin);
}

/*------------------------------------------------------------*/

double scalarop (datum d, double(*fa) (double a, double b))

{
	int i,size;
        double *arrp,*p1,val=0;

        if (!fa) return 0;
        getarr(&d, &p1, &size);
        if (d.vtype!=STRING)
             val = d.arrp->arr[0];
             for (i=0; i<size; i++,p1++) {
                 val = fa(*p1,val);
             }
       return val;
}

datum amin(datum &a, datum &b)
{
        double x,y;
	datum t;

        if (a.type != ARRAY && b.type ==ARRAY) {
	   t = a; a = b; b = t;
	}
        if (a.type==ARRAY) {
          a.val = scalarop(a,domin);
        } else {
          x = a.val;
          y = b.val;
          a.val= min(x,y);
        }
        a.type = VAR;
        a.vtype = NUMBER;
        return a;
}

datum amax(datum &a, datum &b)
{
        double x,y;
	datum t;

        if (a.type != ARRAY && b.type ==ARRAY) {
	   t = a; a = b; b = t;
	}
        if (a.type==ARRAY) {
          a.val = scalarop(a,domax);
        } else {
          x = a.val;
          y = b.val;
          a.val= max(x,y);
        }
        a.type = VAR;
        a.vtype = NUMBER;
        return a;
}

/*------------------------------------------------------------*/

datum transpose(datum &m1)
{
   int i,j,ndim,rows,cols,siz;
   datum m2 = {0};
   double *p;
   ncarray *apnt;

  if (m1.type!=ARRAY) return m1;
  if (apnt=m1.arrp) {
    ndim = apnt->ndim;
    if (ndim==2) { 
      rows = apnt->dim[0];
      cols = apnt->dim[1];
      m2 = copyarr(&m1,&p,&siz);
      for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
           m2.arrp->arr[j*rows+i] = apnt->arr[i*cols+j];
        }
      }
      apnt = m2.arrp;
      ndim = apnt->ndim;
      apnt->dim[0] = cols;
      apnt->dim[1] = rows;
    }
    else return m1;
  }
  return m2;
}

datum matmul(datum &m1, datum &m2)

{
   int i,j,k,ndim,siz;
   int rows1,cols1,rows2,cols2,rows3,cols3;
   double *p,val;
   ncarray *apnt1,*apnt2;
   datum m3;

  if (m1.type!=ARRAY) return m1;
  if ((apnt1=m1.arrp) && (apnt2=m2.arrp)) {
    if (apnt1->ndim==2 && apnt2->ndim==2) { 
      rows1 = apnt1->dim[0];
      cols1 = apnt1->dim[1];
      rows2 = apnt2->dim[0];
      cols2 = apnt2->dim[1];
      rows3 = rows1;
      cols3 = cols2;
      if (cols1 != rows2) {
	ncfprintf (stderr,"matmul: ","matrices different size %d %d\n",cols1,rows2);
      }
      m3 = copyarr(&m1,&p,&siz);
      for (i=0; i<rows3; i++) {
        for (j=0; j<cols3; j++) {
            val = 0;
            for (k=0; k<cols1; k++) {
                val += m1.arrp->arr[i*cols1+k] * m2.arrp->arr[k*cols2+j];
            }
            m3.arrp->arr[i*cols1+j] = val;
        }
      }
    }
    else return m1;
  }
  return m3;
}


int ludcmp(double *a, int n, int *indx, double *d)

/* Given a matrix a[0..n-1][0..n-1], this routine replaces it by
   the LU decomposition of a rowwise permutation of itself. a and
   n are input. a is output, arranged as in equation (2.3.14)
   above; indx[0..n-1] is an output vector that records the row
   permutation effected by the partial pivoting; d is output as
   T1 depending on whether the number of row interchanges was
   even or odd, respectively. This routine is used in combination
   with lubksb to solve linear equations or invert a matrix.
 */

{
  int i,imax,imaxn,j,jn,k;
  double big,dum,sum,temp;
  double *vv;           /* vv stores the implicit scaling of each row. */

  vv= (double *)emalloc(sizeof(double)*n);
  *d=1.0;                       /* No row interchanges yet. */
  for (i=0; i<n; i++) {         /* Loop over rows to get the implicit scaling */
    big=0.0;                    /*   information. */
    for (j=0;j<n;j++)
        if ((temp=fabs(a[i*n+j])) > big) big=temp;
           if (big == 0.0) {
              fprintf (stderr,"Singular matrix in routine ludcmp\n");
              return(0);
           }

     /* No nonzero largest element. */
    vv[i]=1.0/big;              /* Save the scaling. */
  }
  for (j=0; j<n; j++) {         /* Loop over columns of Crout's method. */
    for (i=0; i<j; i++) {       /* Equation (2.3.12) except for i = j. */
        sum=a[i*n+j];
        for (k=0; k<i; k++) sum -= a[i*n+k]*a[k*n+j];
        a[i*n+j]=sum;
    }
    big=0.0;                    /* Initialize search for largest pivot elem. */
    for (i=j; i<n; i++) {       /* Is i=j of equation (2.3.12) and i=j+1...N */
        sum=a[i*n+j] ;           /* of equation (2.3.13). */
        for (k=0; k<j; k++)
            sum -= a[i*n+k] * a[k*n+j];
        a[i*n+j]=sum;
        if ( (dum=vv[i]*fabs (sum)) >= big) {
        /* Is the figure of merit for the pivot better than the best so far? */
          big=dum;
          imax=i;
        }
    }
    if (j != imax) {            /* Do we need to interchange rows? */
      imaxn = imax*n;
      jn = j*n;
      for (k=0; k<n; k++) {     /* Yes, do so... */
        dum = a[imaxn+k];
        a[imaxn+k] = a[jn+k];
        a[jn+k] = dum;
      }
     *d = -(*d);                /* ...  and change the parity of d. */
     vv[imax]=vv[j];            /* Also interchange the scale factor. */
    }
    indx[j]=imax;
    if (a[j*n+j] == 0.0) a[j*n+j]=1e-30;

/* If the pivot element is zero the matrix is singular (at least
   to the precision of the algorithm). For some applications on
   singular matrices, it is desirable to substitute TINY for
   zero.
 */

    if (j != n) {               /* Now, finally, divide by the pivot element. */
      dum=1.0 / (a[j*n+j]);
      for (i=j+1; i<n; i++) a[i*n+j] *= dum;
    }
  }     /* Go back for the next column in the reduction. */
  efree(vv);
  return 1;
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void lubksb(double *a, int n, int *indx, double b[])

/* Solves the set of n linear equations A-X = B. Here a[0..n-1]
   [0..n-1] is input, not as the matrix A but rather as its LU
   decomposition, determined by the routine ludcmp. indx[0..n-1]
   is input as the permutation vector returned by ludcmp. b[0..n-1]
   is input as the right-hand side vector B, and returns with
   the solution vector X. a, n, and indx are not modified by this
   routine and can be left in place for successive calls with
   different right-hand sides b. This routine takes into account
   the possibility that b will begin with many zero elements, so
   it is efficient for use in matrix inversion.
 */

{
    int i,ii= -1,ip,j;
    double sum;

  for (i=0; i<n; i++) { /* When ii is set to a positive value, it will become */
    ip=indx[i];         /*  the index of the first nonvanishing element of b. */
    sum=b[ip];          /* We now do the forward substitution, eq(2.3.6). */
    b[ip]=b[i];         /* The only new wrinkle is to unscramble the permut'n */
        if (ii>=0)      /* as we go. */
          for (j=ii; j<i; j++) sum -= a[i*n+j]*b[j];
        else if (sum) ii=i; /* A nonzero element was encountered, so from */
        b[i]=sum;      /* now on we will have to do the sums in the loop above*/
   }
   for (i=n-1; i>=0; i--) { /* Now we do the backsubstitution, eqn (2.3.7).*/
        sum=b[i];
        for (j=i+1; j<n; j++) sum -= a[i*n+j]*b[j];
        b[i]=sum/a[i*n+i];  /* Store a component of the solution vector X. */
   }                       /* All done! */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

datum matinv(datum &m)

{
    int i,j,n,*indx,siz;
    double *a,d,*p,*mp,*y,*col;
    ncarray *apnt;
    datum m2;
 
  if (m.type!=ARRAY) return m;
  if (apnt=m.arrp) {
    if (apnt->ndim==2) {
      n = apnt->dim[0];
      mp = m.arrp->arr;
      m2 = copyarr(&m,&p,&siz);
      for (i=0; i<siz; i++) 
        *(p+i) = *(mp+i);
      a = m2.arrp->arr;
      indx = (int *)emalloc(sizeof(int)*n);
      col = (double *)emalloc(sizeof(double)*n);
      y = (double *)emalloc(sizeof(double)*n*n);
      ludcmp(a,n,indx,&d); 
      for(j=0;j<n;j++) {
        for(i=0;i<n;i++) 
	   col[i]=0.0;
        col[j]=1.0;
        lubksb(a,n,indx,col);
        for(i=0;i<n;i++) 
 	   y[i*n+j]=col[i];
      }
      m2.arrp->arr = y;
      efree(a);
    }
    efree(indx);
    efree(col);
  }
  return m2;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

datum matsolve(datum &m, datum &v)

{
    int i,j,n,*indx,siz;
    double *a,*b, d,*p,*mp,*y;
    ncarray *apnt;
    datum m2;

  if (m.type!=ARRAY) return m;
  if (apnt=m.arrp) {
    if (apnt->ndim==2) {
      n = apnt->dim[0];
      mp = m.arrp->arr;
      m2 = copyarr(&m,&p,&siz);
      for (i=0; i<siz; i++) 
        *(p+i) = *(mp+i);
      a = m2.arrp->arr;
      indx = (int *)emalloc(sizeof(int)*n);
      if (v.arrp->dim[0] == n) {
        b = v.arrp->arr;
        ludcmp(a,n,indx,&d); 
        lubksb(a,n,indx,b);
        m2.arrp->arr = b;
        m2.arrp->ndim = 1;
        m2.arrp->dim[1] = 0;
      }
      efree(a);
      efree(indx);
    }
  }
  return m2;
}

/****************************************/

double errcheck(double d, const char *s)	/* check result of library call */
	         
	        
{
	if (errno == EDOM) {
		errno = 0;
		execerror(s, "argument out of domain");
	} else if (errno == ERANGE) {
		errno = 0;
		execerror(s, "result out of range");
	}
	return d;
}

/****************************************/

void fft(double *real, double *imag, int asiz, int param)

/* Fourier transform of data in "real" and "imag", stored as 
   alternating real and imaginary values in array.  Internal 
   "four()" FFT routine has Fortran index usage and ignores 
   data[0].  First value must be in data[1].  Size of array 
   passed to "four()" must be power of 2, so allocate more 
   memory if necessary. 
*/

{
   int i,j,isign,arrsiz,fdatsiz;
   double *fdata=0,re,im,scale;

 if (!asiz) return;
 for (i=1,arrsiz=1; i<20; i++ ) {    /* find next higher power of 2 */
      arrsiz <<= 1;
      if (arrsiz >= asiz) break;
 }

 fdatsiz = (arrsiz+1) * 2;	   /* allocate space for next power of 2 */
 if ((fdata=(double*)emalloc(fdatsiz*sizeof(double)))==NULL) {
  ncfprintf (stderr,"fft: unable to allocate internal array of size %d.\n",
						fdatsiz);
  execerror ("Memory allocation error.",0);
 }
 else {
       double sum,mean;

    if (param!=IFFT) {			/* if not inverse transform */
      for (i=0; i<asiz; i++) {          /* find mean of real part */
         sum = real[i];
      }
    mean = sum / asiz;
    }
    else mean = 0.0;
    for (i=j=0; i<arrsiz; i++,j+=2) { /* copy  data into new array */
       if (i<asiz) {
         fdata[j] = real[i];
         fdata[j+1] = imag[i];
       }
       else {
         fdata[j] = mean;		/* zero out the rest */
         fdata[j+1] = 0.0;
       }
     }
  }
/* for (i=0; i<arrsiz; i++) 
     ncfprintf(stdout,"%10g %10g\n",real[i],imag[i]);       /* */

/* ncfprintf (stderr,"Running fourier transform of orig %d size %d.\n",
			asiz,arrsiz); /* */

 if (param==IFFT) isign = -1; 		/* inverse FFT */
 else             isign =  1;

 four1(fdata-1,arrsiz,isign); 

  switch (param) {

   case FFT: 
	for (i=j=0; i<asiz; i++,j+=2) { /* copy transf data into orig array */
         real[i] = fdata[j]; 
         imag[i] = fdata[j+1];
        }
	break;
   case IFFT: scale = 1.0/arrsiz;
    	for (i=j=0; i<asiz; i++,j+=2) { /* copy transf data into orig array */
         real[i] = fdata[j] * scale;
         imag[i] = fdata[j+1] * scale;
        }
	break;
   case PFFT: 
	for (i=j=0; i<asiz; i++,j+=2) { /* copy transf data into orig array */
          re = fdata[j];
          im = fdata[j+1];
          real[i] = sqrt (re*re + im*im);
          imag[i] = 0.0;
        }
     break;
   case ACOV: 
        for (i=j=0; i<arrsiz; i++,j+=2) { /* copy transf data into orig array */
          re = fdata[j];
          im = fdata[j+1];
          fdata[j] = (re*re + im*im);	/* mul by complex conjugate */
          fdata[j+1] = 0;
        }
        four1(fdata-1,arrsiz,-1); 	/* reverse transform the product */
        scale = 1.0/arrsiz;
    	for (i=j=0; i<asiz; i++,j+=2) { /* copy transf data into orig array */
          real[i] = fdata[j] * scale;
          imag[i] = fdata[j+1] * scale;
        }
     break;
  }
 if (fdata) efree (fdata);
}

/****************************************/

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(double *data, int nn, int isign)

/* Replaces data by its discrete Fourier transform,
   if isign is input as 1; or replaces data by nn
   times its inverse discrete Fourier transform,
   if isign is input as -1.  data is a complex array
   of length nn, input as a real array data[1..2*nn].
   nn MUST be an integer power of 2 (this is not checked for!).
   Data alternates real and imaginary in data array.
   This routine ignores data[0] because it imitates FORTRAN indexing.
*/

{
   int n, mmax,m,j,istep,i;
   double wtemp,wr,wpr,wpi,wi,theta;    /* double precis for trig recurrences */
   double tempr, tempi;

  n = nn << 1;
  j = 1;
  for (i=1; i<n; i+=2) {                /* the bit-reversal section */
    if (j>i) {
        SWAP(data[j],data[i]);          /* exchange the two complex numbers */
        SWAP(data[j+1],data[i+1]);
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
        j -= m;
        m >>= 1;
    }
    j += m;
  }
  mmax = 2;             /* Here begins the Danielson-Lanczos section */
  while (n>mmax) {      /* outer loop executed log2 nn times */
    istep = mmax << 1;
    theta = isign * (6.28318530717959 / mmax);    /* init for trig recurrence */
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m=1; m<mmax; m+=2) {           /* the two nested loops */
        for (i=m; i<=n; i+= istep) {
           j = i+mmax;                  /* the Danielson-Lanczos formula: */
           tempr = wr*data[j] - wi*data[j+1];
           tempi = wr*data[j+1] + wi*data[j];
           data[j]   = data[i]   - tempr;
           data[j+1] = data[i+1] - tempi;
           data[i] += tempr;
           data[i+1] += tempi;
        }                               /* trigonometric recurrence */
        wr = (wtemp=wr)*wpr - wi*wpi + wr;
        wi = wi*wpr + wtemp*wpi + wi;
    }
    mmax = istep;
  }
}


/****************************************/

double mulreal(double a, double b, double c, double d);
double mulimag(double a, double b, double c, double d);

void dfft(double *z1real, double *z1imag, int n, int param)

/* This routine, while it can transform an array of any size,
   does not do an inverse FFT.  Therefore it is not currently used.  */

/* On input, z1 and z2 are complex n-vectors,
   n is length of z1 and z2,
   inzee indicates whether z1 or z2 is to be transformed.

   Both z1 and z2 are used as work arrays.
   
   On output, z1 or z2 contains the desired transform
   in the correct order. 
   Inzee indicates whether z1 or z2 contains the transform.
   inzee = 1  => z1 contains transform
   inzee = 2  => z2 contains transform

   Method:  

     The integer n is divided into its prime factors (up to
   a point).  For each such factor p, the p-transform of 
   appropriate p-subvectors of z1 (or z2) is calculated
   in fftstp and stored in a suitable way in z2 (or z1).

   Adapted from fortran routine in:

	 Conte and de Boor,
	 'Elementary Numerical Analysis' 1980
		Chapter 6, "Fast Fourier Transforms", p 283
*/

#define FFTSIZ 2048
 
{
    int i;
/*  static double z1real[FFTSIZ] = 0;    */
/*  static double z1imag[FFTSIZ] = 0;	*/
/*  static double z2real[FFTSIZ] = 0;	/* */
/*  static double z2imag[FFTSIZ] = 0;	/* */
    double *z2real, *z2imag, real,imag; 
    static int prime[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    int inzee,after,before,next,nextmx,now;
    char *emalloc(unsigned int n);

 inzee = 1;
 nextmx = 12;
/* n = readin(z1real,z1imag); 		/* read first list */

 z2real = (double *)emalloc(sizeof(double)*n);
 z2imag = (double *)emalloc(sizeof(double)*n);
 if (!z1real || !z2real) {
	ncfprintf (stderr,"fft: can't allocate work space...\n");
	return;
 }

 after = 1;
 before = n;
 next = 1;
 while (before != 1) {
   if ((before / prime[next-1]) * prime[next-1] < before) {
    next++;
     if (next <= nextmx) continue;
     else {
	now = before;
	before = 1;
     }
   }
   else {
     now = prime[next-1];
     before = before / prime[next-1];
   }

   if (inzee == 1)
      fftstp(z1real, z1imag, after, now, before, z2real, z2imag);
   else
      fftstp(z2real, z2imag, after, now, before, z1real, z1imag);
   inzee = 3 - inzee;
   if (before != 1) after *= now;

 }    /* while (before ) */

if (inzee==1) {
   if (param==FFT) {
/*    for (i=0; i<n; i++) 
     ncfprintf(stdout,"%10g %10g\n",z1real[i],z1imag[i]); /* */
   }
   else if (param==PFFT) {
     for (i=0; i<n; i++) {
       real = z1real[i];
       imag = z1imag[i];
       z1real[i] = sqrt (real*real + imag*imag);
       z1imag[i] = 0.0;
/*     ncfprintf(stdout,"%10g %10g\n",z1real[i],z1imag[i]); /* */
    }
  } 	/*  else (param==PFFT) */
}	/* if (inzee==1) */
else {
   if (param==FFT) {
    for (i=0; i<n; i++) {
/*     ncfprintf(stdout,"%10g %10g\n",z2real[i],z2imag[i]); /* */
     z1real[i] = z2real[i];
     z1imag[i] = z2imag[i];
    }
   }
   else if (param==PFFT) {
    for (i=0; i<n; i++) {
       real = z2real[i];
       imag = z2imag[i];
       z1real[i] = sqrt (real*real + imag*imag);
       z1imag[i] = 0.0;
     }
   }
}
 efree (z2real);
 efree (z2imag);

}

/****************************************/

void fftstp (double *zinr, double *zini, int after, int now, 
			int before, double *zoutr, double *zouti)
{
    int ia,ib,in,j;
    static double angle;
    static double argreal,argimag;
    static double omegreal,omegimag;
    static double valreal,valimag;
    static double tvalreal,targreal;
    int pnt;

 angle = 2 * MPI / (now * after);
 omegreal = cos(angle);
 omegimag = -sin(angle); 
 argreal  = 1.0;
 argimag  = 0.0;
 for (j=0; j<now; j++) {
   for (ia=0; ia<after; ia++) {
     for (ib=0; ib<before; ib++) {
       pnt = ((now-1)*before+ib)*after+ia;
       valreal = *(zinr+pnt);
       valimag = *(zini+pnt);
       for (in=now-2; in>=0; in--)  {
	 pnt = (in*before+ib)*after+ia;
	 tvalreal = mulreal(valreal,valimag,argreal,argimag);
	 tvalreal = tvalreal + *(zinr+pnt);
	 valimag = mulimag(valreal,valimag,argreal,argimag);
	 valimag = valimag + *(zini+pnt);
	 valreal = tvalreal;
      }
      pnt = (ib*now+j)*after+ia;
      *(zoutr+pnt) = valreal;
      *(zouti+pnt) = valimag;
    }
    targreal = mulreal(argreal,argimag,omegreal,omegimag);
    argimag = mulimag(argreal,argimag,omegreal,omegimag);
    argreal = targreal;
  }
 }
}

/****************************************/

double mulreal(double a, double b, double c, double d)
{
 return (a*c - b*d);
}

/****************************************/

double mulimag(double a, double b, double c, double d)
{
 return (b*c + a*d);
}

/****************************************/

static char notstr[] = {"variable not a string"};
static char notnum[] = {"variable not a number"};
   
datum xisnum(datum &d)
{
	switch (d.vtype) {
	  case 0:
	  case NUMBER:
		d.val = 1;
		break;
	  case LITCHAR:
		d.val = 0;
		break;
	  case STRING:
		d.val = 0;
		break;
	}
	d.vtype = NUMBER;
	return d;
}


double dostrlen(char *str)
{
	return(strlen (str));
}

datum xstrlen(datum &x)
{
        switch (x.vtype) {
          case 0:
          case NUMBER:
                execerror (notstr,0);
		return x;
                break;
          case LITCHAR:
                x.val = 1;
		return x;
                break;
          case STRING:
                return domonofuncs(x,dostrlen);
                break;
        }
}

datum xstrtok(datum &d1, datum &d2)
{
	if (d1.vtype!=STRING && d1.val!=0.0) execerror (notstr,0);
 	if (d2.vtype!=STRING) execerror ("separator is not a string",0); 
	d1.str = strtok((char *)d1.str, d2.str);
	d1.vtype = STRING;
	return d1;
}

datum xstrcmp(datum &d1, datum &d2)
{
	if (d1.vtype!=STRING && d1.val!=0.0) execerror (notstr,0);
 	if (d2.vtype!=STRING && d2.val!=0.0) execerror (notstr,0); 
	d1.val = strcmp(d1.str, d2.str);
	d1.vtype = NUMBER;
	return d1;
}

datum xstrstr(datum &d1, datum &d2)
{
    const char *p;

	if (d1.vtype!=STRING && d1.val!=0.0) execerror (notstr,0);
 	if (d2.vtype!=STRING && d2.val!=0.0) execerror (notstr,0); 
	if ((p=strstr(d1.str, d2.str))==NULL) 
	   d1.val = -1;
	else d1.val = p - d1.str;
	d1.vtype = NUMBER;
	return d1;
}

datum xindex(datum &d1, datum &d2)
{
    const char *p;
    char ch;

	if (d1.vtype!=STRING && d1.val!=0.0) execerror (notstr,0);
 	if (d2.vtype==STRING) ch = *d2.str;
        else if (d2.vtype==LITCHAR) ch = d2.val;
        else execerror (notstr,0); 
	if ((p=index(d1.str, ch))==NULL) 
	   d1.val = -1;
	else d1.val = p - d1.str;
	d1.vtype = NUMBER;
	return d1;
}

datum xrindex(datum &d1, datum &d2)
{
    const char *p;
    char ch;

	if (d1.vtype!=STRING && d1.val!=0.0) execerror (notstr,0);
 	if (d2.vtype==STRING) ch = *d2.str;
        else if (d2.vtype==LITCHAR) ch = d2.val;
        else execerror (notstr,0); 
	if ((p=rindex(d1.str, ch))==NULL) 
	   d1.val = -1;
	else d1.val = p - d1.str;
	d1.vtype = NUMBER;
	return d1;
}


datum xsubstr(datum &d1, datum &d2, datum &d3)
{
	int i,len,n;
	char *d;
	Symbol *sp;
	const char *s;

	if (d1.vtype!=STRING && d1.val!=0.0) execerror (notstr,0);
 	if (d2.vtype!=NUMBER && d2.val!=0.0) execerror (notnum,0); 
        i = (long int)d2.val;
        n = (long int)d3.val;
        if (n==0) n = 1000;
	s = d1.str;
	len = strlen(s);
	d = emalloc(len+1); 	/* +1 for '\0' */
	if (n > len) n = len;
        strncpy (d,s+i,n);
        d[n] = '\0';
        sp = install (d,STRING,0.0);
	efree (d);
	d1.str = sp->name;
	d1.vtype = STRING;
	return d1;
}

/*
datum xatof(datum &d1)
{
     double val;

	if (d1.vtype!=STRING && d1.val!=0.0) execerror (notstr,0);
	d1.val = atof((char *)d1.str);
	d1.vtype = NUMBER;
	return d1;
}
/* */


datum xatof(datum &d1)
{
     double val;
     char *str;
     Symbol *s;
     char xbuf[150];

	if (d1.vtype!=STRING && d1.val!=0.0) execerror (notstr,0);
	str = (char *)d1.str;
	if (fread_expr) {           // allows string to contain expressions (no spaces)
            val = fread_interp(str);
	} else {
	   if (isalpha(*str) || (*str=='_')) {

	     val = getvarval(str);                     // look in "command line" variable list 
	     if (val==LARGENUM) {
	        s = lookup(str);                       // look in orig symbol table 
	        if (s) val = s->val;
	        else {
		  sprintf (xbuf,
                   "# Missing def for var '%.20s' while converting string.\n", str);
		    execerror ("warning,",xbuf);
	           val = 0;
	        }
	     }
	   }
	   else   {  
              val = atof(str);
	   }
        }
	d1.val = val;
	d1.vtype = NUMBER;
	return d1;
}
/* */

/* - - - - - - - - - - - - - - - - - - - - - */

datum xccstr(datum &d1)
{
	int i,len;
	char *d,*dd;
	Symbol *sp;
	const char *s;

	if (d1.vtype!=STRING && d1.val!=0.0) execerror (notstr,0);
	s = d1.str;
	len = strlen(s);
	dd = d = emalloc(len+1); 		/* +1 for '\0' */
        for (i=0; i<len; i++) {
         if (*(s+i) >= 0x20) *d++ = *(s+i);
        }   
	*d = 0;
        sp = install (dd,STRING,0.0);
	efree (dd);
	d1.str = sp->name;
	d1.vtype = STRING;
	return d1;
}

