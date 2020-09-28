/* Module rndev in Program nc */

/* Generates binomial, Poisson, and gamma deviates (random sequences).
 */

#ifdef __cplusplus
extern "C" {
#endif

#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif
#include <stdio.h>

#ifdef __cplusplus
}
#endif


#include "ncsub.h"
#include "nconst.h"
#include "ncio.h"

void restorstate(void);
char *emalloc(unsigned int n);
void *efree(void *ptr);

/*------------------------------------*/

int binomdev(double pp, int n)
              
/* Taken from Press et al. 1988, "Numerical Recipes in C", p 223.
   Returns as a floating point number (double) an integer value 
   that is a random deviate drawn from a binomial distribution
   of n trials each of probability pp, using drand() as a
   source of uniform random deviates.
*/

{
   register int j;
   register int bnl;
   static int nold = -1;  			     /* must be static */
   static double pold = -1, pc, plog,pclog,en,ex,oldg;  /* must be static */
   static double am,em,g,p,sq,t,y;
   double drand(), gammln2(double xx);
   void restorstate();

#define DIRECTHRESH 25

   //if (pp > 1.0) pp = 1.0;
   //if (pp < 0.0) pp = 0.0;
   //p = (pp <= 0.5? pp : 1.0-pp);
  if (pp <= 0.5) {
   p = pp;
   if (n < DIRECTHRESH) {		/* use direct method if n is small */
      bnl = 0;
      for (j=1; j<=n; j++) 
	 if (drand() < p) bnl++;
   }
   else if ((am=n*p) < 1.0) {		/* if fewer than 1 event is expected */
      g = exp(-am);			/*  in 25 or more trials, then */
      t = 1.0;				/*  use Poisson method */
      for (j=0; j<n; j++) {
         t *= drand();
	 if (t < g) break;
      }
      bnl = j;
   }
   else {				/* use rejection method */
      if (n != nold) {
	 en = n;
         ex = n+1;
	 oldg = gammln2(ex);
	 nold = n;
      }
      if (p != pold) {			/* if p has changed */
	pc = 1.0 - p;
	plog = log(p);
	pclog = log(pc);
	pold = p;
      }
      sq = sqrt(2.0*am*pc);		/* rejection method with Lorentzian */
      do {				/* comparison function */
	do {
	   y = tan(MPI*drand());
	   em = sq*y+am;
	} while (em < 0.0 || em >= ex);	/* reject */
	em = floor(em);
	t = 1.2 * sq * (1.0+y*y) * exp(oldg-gammln2(em+1.0)
		- gammln2(ex-em) + em * plog+ (en-em) * pclog);

      } while (drand() > t);		/* reject. This happens about 1.5 */
      bnl = int(em);			/* times per deviate, on average. */
   }
  }
  else {  /* (pp > 0.5) */
   p = 1.0 - pp;
   if (n < DIRECTHRESH) {		/* use direct method if n is small */
      bnl = 0;
      for (j=1; j<=n; j++) 
	 if (drand() < p) bnl++;
      bnl = n - bnl;
   }
   else if ((am=n*p) < 1.0) {		/* if fewer than 1 event is expected */
      g = exp(-am);			/*  in 25 or more trials, then */
      t = 1.0;				/*  use Poisson method */
      for (j=0; j<n; j++) {
         t *= drand();
	 if (t < g) break;
      }
      bnl = n - j;
   }
   else {				/* use rejection method */
      if (n != nold) {
	 en = n;
         ex = n+1;
	 oldg = gammln2(ex);
	 nold = n;
      }
      if (p != pold) {			/* if p has changed */
	pc = pp;
	plog = log(p);
	pclog = log(pc);
	pold = p;
      }
      sq = sqrt(2.0*am*pc);		/* rejection method with Lorentzian */
      do {				/* comparison function */
	do {
	   y = tan(MPI*drand());
	   em = sq*y+am;
	} while (em < 0.0 || em >= ex);	/* reject */
	em = floor(em);
	t = 1.2 * sq * (1.0+y*y) * exp(oldg-gammln2(em+1.0)
		- gammln2(ex-em) + em * plog+ (en-em) * pclog);

      } while (drand() > t);		/* reject. This happens about 1.5 */
      bnl = n - int(em);		/* times per deviate, on average. */
   }
   //bnl = n -bnl;
  }
//   restorstate();			/* restore default random state  */
   return (bnl);
}

/*------------------------------------*/

double fact(double f) {

  if (f > 1) return (f*fact(f-1));
  else       return 1;

}

/*------------------------------------*/

double gammln2(register double xx)

/* Approximation to the log(gamma) function */
/*  Table lookup is used. Table expands up to 200000 */

#define MAXGAMTABSIZ 200000

{
   register int x;
   double f,val;
   static int gamtabsiz = 0;
   static double *gamtab=(double*)NULL;

#define Y 0.91893853320467267 		    /* 0.5*log(2.0*MPI); */

/* ncfprintf (stderr,"gammln2 %g\n",xx); /* */

 x = int (xx);
 if (x < MAXGAMTABSIZ) {		/* if too big, use approximation */
   if (x >= gamtabsiz) { 		/* make bigger table */
        register int i,j,mingam;
	double *ngamtab;
	static int ngamtabsiz;

     ngamtabsiz = int(xx/1024+1.1) * 1024;	/* make table bigger */
    if (ngamtabsiz > MAXGAMTABSIZ) ngamtabsiz = MAXGAMTABSIZ;

    /* ncfprintf (stderr,"making new gamtab size %d\n",ngamtabsiz); /* */

    ngamtab = (double *)emalloc(ngamtabsiz*sizeof(double));
    if (gamtab) {			
      for (i=0; i<gamtabsiz; i++)	/* copy old table into new */
          ngamtab[i] = gamtab[i]; 
      for (; i<ngamtabsiz; i++)		/* fill rest with zeroes */
          ngamtab[i] = 0.0;
      efree(gamtab);
      mingam = gamtabsiz-1;
    } 
    else {
      ngamtab[0] = 0;
      ngamtab[1] = 0;
      mingam = 1;
    }
    gamtab = ngamtab;
    gamtabsiz = ngamtabsiz;
    for (i=mingam+1,f=j=i-1; i<gamtabsiz; f=j=i++) {
       gamtab[i] = gamtab[j] + log(f);
    }
  /* ncfprintf (stderr,"gamtab done\n"); */
  } 
  return (gamtab[x]);			 /* use table lookup */ 
 }
 else {
   val = (xx - 0.5) * log(xx) - xx + Y + 0.083333754/xx;
 }				/* error corrected < 2e-10 for xx > 2000 */
 return val;
}

#undef Y

/*------------------------------------*/

double gammln(double xx)

/* Taken from Press et al. 1988, "Numerical Recipes in C", 2nd
   ed, p 214.  Returns the value ln(gamma(xx)) where xx > 0.
*/

{
   double x,y,tmp,ser;
   static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,
                        -1.2317395172450155, 0.1208650973866179e-2,
                        -0.5395239384953e-5};
   int j;

   y = x = xx;
   tmp = x + 5.5;
   tmp -= (x+0.5) * log(tmp);
   ser = 1.000000000190015;
   for (j=0; j<=5; j++) {
     ser += cof[j] / ++y;
   }
   return -tmp+log(2.5066282746310005*ser/x);
}

/*------------------------------------*/
#include <errno.h>

int poisdev(double xm)
             
/*
 Return a poisson distribution which on average equals "xm".

 Taken from Press et al. 1988, "Numerical Recipes in C", p 222.
   Returns as a floating point number (double) an integer value 
   that is a random deviate drawn from a Poisson distribution
   of mean xm as a source of uniform random deviates.
*/

{
  static double sq,alxm,g,oldm= -1.0;
  double em,t,y;
  double drand(),gammln2(double xx);
  void restorstate();

#ifdef DEBUG
  if (debug & 1) ncfprintf (stderr,"poisdev\n");
#endif

  if (xm < 12.0) {			/* use direct method */
     if (xm != oldm) {
  	oldm = xm;
	g = exp(-xm);			/* if xm is new, compute exp */
     }
     em = -1;
     t = 1.0;
     do {
       em += 1.0;			/* instead of adding exp deviates */
       t *= drand();			/* we multiply uniform deviates */
     } while (t > g);
  }
  else {				/* use rejection method */
    if (xm != oldm) {
      oldm = xm;
      sq = sqrt(2.0*xm);
      alxm = log(xm);
      g = xm * alxm - gammln2(xm+1.0);
    }
    do {
       do {
	  y = tan(MPI * drand());
	  em = sq*y+xm;
       } while (em < 0.0);		/* reject if zero probability */
       em = (int)(floor(em));
       t = 0.9 * (1.0+y*y) * exp(em*alxm-gammln2(em+1.0)-g);
    } while (drand() > t);
  }
 // restorstate();			/* restore default random state */
 return (int)(em);
}

/*------------------------------------*/

double gasdev(void)

/* Taken from Press et al. 1992, "Numerical Recipes in C", p 289.
   Returns a normally distributed deviate with zero mean and
   unit standard deviation, using drand() as a source of 
   uniform random deviates.
*/

{
    static int iset=0;
    double fac,rsq,v1,v2,retval;
    double drand();

  if (iset==0) {				/* Extra deviate ready ? */
    do {
        v1 = 2.0*drand() - 1.0;			/* Pick 2 uniform numbers */
        v2 = 2.0*drand() - 1.0;			/*  -1 to +1 in each dir */
        rsq = v1*v1 + v2*v2;			/* In unit circle ? */
    } while (rsq >= 1.0 || rsq == 0.0);         /*  not ? try again */

    fac = sqrt(-2.0*log(rsq)/rsq);

			/* now make Box-Muller transformation to get 2
			 * deviates.  Return one of them and save
			 * the other for next time */
/*    gset = v1 * fac;
    iset = 1;
*/
    retval = v2 * fac;
   }
/*
   else {
     iset = 0;
     retval =  gset; 
  }
*/
//  restorstate();			/* restore default random state */
  return retval;
}

/*------------------------------------*/

double rgasdev(int ngen)

/* Taken from Press et al. 1992, "Numerical Recipes in C", p 289.
   Returns a normally distributed deviate with zero mean and
   unit standard deviation, using drand() as a source of 
   uniform random deviates.
*/

{
    static int iset=0;
    double fac,rsq,v1,v2,retval;
    double rrand(int ngen);

  if (iset==0) {				/* Extra deviate ready ? */
    do {
        v1 = 2.0*rrand(ngen) - 1.0;		/* Pick 2 uniform numbers */
	if (v1 < -10) return (0);		/* stop if wrong ngen */
        v2 = 2.0*rrand(ngen) - 1.0;		/*  -1 to +1 in each dir */
        rsq = v1*v1 + v2*v2;			/* In unit circle ? */
    } while (rsq >= 1.0 || rsq == 0.0);         /*  not ? try again */

    fac = sqrt(-2.0*log(rsq)/rsq);

			/* now make Box-Muller transformation to get 2
			 * deviates.  Return one of them and save
			 * the other for next time */
/*    gset = v1 * fac; */
/*    iset = 1; */
    retval = v2 * fac;
   }
/*   else {
     iset = 0;
     retval =  gset; 
  }
*/
  return retval;
}

/*------------------------------------*/

static double gamdev_nr(int ia)

/* Taken from Press et al. 1992, "Numerical Recipes in C", p 292.
   Returns as a floating point number an integer value 
   that is a random deviate drawn from a gamma distribution
   of n trials each of order ia, using drand() as a
   source of uniform random deviates.
*/

{
   int j;
   double am,e,s,v1,v2,x,y;
   double drand();

  if (ia < 1) ia = 1;
  if (ia < 6) {				/* use direct method, adding  */
     x = 1.0;				/* waiting times. */
     for (j=1; j<=ia; j++) x *= drand();
     x = -log(x);
  } else {			/* use rejection method */
   do {
      do {
         do {					/* These 4 lines generate */
	     v1 = 2.0*drand() - 1.0;		/* the tan of a */
	     v2 = 2.0*drand() - 1.0;		/* random angle */
         } while (v1*v1+v2*v2 > 1.0);
         y = v2/v1;
         am = ia - 1;
         s = sqrt(2.0 * am + 1.0);
         x = s * y + am;			/* Decide whether to reject x */
      } while (x <= 0.0);			/* Reject in reg.  of 0 prob */
      e = (1.0+y*y) * exp(am*log(x/am)-s*y);	/* Ratio of prob. fn. to comp */
   } while (drand() > e);			/* Reject on 2nd uniform dev. */
 }
// restorstate();			/* restore default random state */
 return x;
}

/*------------------------------------*/

/* 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Modified from original GSL source code */

#include <math.h>
#include "nconst.h"

static double gamma_large (double a);
static double gamma_int   (unsigned int a);
static double gamma_frac  (double a);

double drand(void);
double drand_pos(void);

/*-----------------------------------------------------------------*/

/* The Gamma distribution of order a>0 is defined by:

   p(x) dx = {1 / \Gamma(a) b^a } x^{a-1} e^{-x/b} dx

   for x>0.  If X and Y are independent gamma-distributed random
   variables of order a1 and a2 with the same scale parameter b, then
   X+Y has gamma distribution of order a1+a2.

   The algorithms below are from Knuth, vol 2, 2nd ed, p. 129. */

double gamdev (double a)
{
  /* assume a > 0 */
  unsigned int na = (int)floor (a);

  if (a == na)
    {
      return gamma_int (na);
    }
  else if (na == 0)
    {
      return gamma_frac (a);
    }
  else
    {
      return (gamma_int (na) + gamma_frac (a - na)) ;
    }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double gamma_int (unsigned int a)
{
  if (a < 12)
    {
      unsigned int i;
      double prod = 1;

      for (i = 0; i < a; i++)
	{
	  prod *= drand();
	}

      /* Note: for 12 iterations we are safe against underflow, since
	 the smallest positive random number is O(2^-32). This means
	 the smallest possible product is 2^(-12*32) = 10^-116 which
	 is within the range of double precision. */

      return -log (prod);
    }
  else
    {
      return gamma_large((double) a);
    }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

static double gamma_large (double a)
{
  /* Works only if a > 1, and is most efficient if a is large

     This algorithm, reported in Knuth, is attributed to Ahrens.  A
     faster one, we are told, can be found in: J. H. Ahrens and
     U. Dieter, Computing 12 (1974) 223-246.  */

  double sqa, x, y, v;
  sqa = sqrt (2 * a - 1);
  do
    {
      do
	{
	  y = tan (MPI * drand());
	  x = sqa * y + a - 1;
	}
      while (x <= 0);
      v = drand ();
    }
  while (v > (1 + y * y) * exp ((a - 1) * log (x / (a - 1)) - sqa * y));

  return x;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

static double gamma_frac (double a)
{
  /* This is exercise 16 from Knuth; see page 135, and the solution is
     on page 551.  */

  double p, q, x, u, v;
  p = M_E / (a + M_E);
  do
    {
      u = drand();
      v = drand_pos();

      if (u < p)
	{
	  x = exp ((1 / a) * log (v));
	  q = exp (-x);
	}
      else
	{
	  x = 1 - log (v);
	  q = exp ((a - 1) * log (x));
	}
    }
  while (drand () >= q);

  return x;
}
