/* Module binom in Program nc */

/* Generates binomial and Poisson deviates (random sequences).
 */

#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif

#include "ncsub.h"
#include "ncio.h"

/*------------------------------------*/

double binomdev(double pp, int n)
              
/* Taken from Press et al. 1988, "Numerical Recipes in C", p 223.
   Returns as a floating point number an integer value 
   that is a random deviate drawn from a binomial distribution
   of n trials each of probability pp, using drand() as a
   source of uniform random deviates.
*/

{
   int j;
   static int nold = -1;  			     /* must be static */
   static double am,em,g,angle,p,bnl,sq,t,y;
   static double pold = -1, pc, plog,pclog,en,oldg;  /* must be static */
   double drand(), gammln(double xx);
   void restorstate();

   if (pp > 1.0) pp = 1.0;
   if (pp < 0.0) pp = 0.0;
   p = (pp <= 0.5? pp : 1-pp);
   am = n * p;				/* am = mean of deviate */
   if (n < 25) {			/* use direct method if n is small */
      bnl = 0.0;
      for (j=1; j<=n; j++) 
	 if (drand() < p) bnl += 1.0;
   }
   else if (am < 1.0) {			/* if fewer than 1 event is expected */
      g = exp(-am);			/*  in 25 or more trials, then */
      t = 1.0;				/*  use Poisson method */
      for (j=0; j<=n; j++) {
         t *= drand();
	 if (t < g) break;
      }
      bnl = (j <= n ? j : n);
   }
   else {				/* use rejection method */
      if (n != nold) {
	 en = n;
	 oldg = gammln(en+1.0);
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
	   angle = drand();
	   if (angle == 0.5) angle = .49999999999;	
	   y = tan(MPI*angle);
	   em = sq*y+am;
	} while (em < 0.0 || em >= (en+1.0));	/* reject */
	em = floor(em);
	t = 1.2 * sq * (1.0+y*y) * exp(oldg-gammln(em+1.0)
		- gammln(en-em+1) + em * plog+ (en-em) * pclog);

      } while (drand() > t);		/* reject. This happens about 1.5 */
      bnl = em;				/* times per deviate, on average. */
   }
   if (p != pp) bnl = n -bnl;
   restorstate();			/* restore default random state  */
   return bnl;
}

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

int poisdev(double xm)
             
/*
 Return a poisson distribution which on average equals "xm".

 Taken from Press et al. 1988, "Numerical Recipes in C", p 222.
   Returns as a floating point number an integer value 
   that is a random deviate drawn from a Poisson distribution
   of mean xm as a source of uniform random deviates.
*/

{
  static double sq,alxm,g,oldm= -1.0;
  double em,t,y;
  double drand(),gammln(double xx);
  void restorstate();

#ifdef DEBUG
  if (debug & 1) ncprintf (stderr,"poisdev\n");
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
      g = xm * alxm - gammln(xm+1.0);
    }
    do {
       do {
	  y = tan(MPI * drand());
	  em = sq*y+xm;
       } while (em < 0.0);		/* reject if zero probability */
       em = (int)(floor(em));
       t = 0.9 * (1.0+y*y) * exp(em*alxm-gammln(em+1.0)-g);
    } while (drand() > t);
  }
 restorstate();				/* restore default random state */
 return (int)(em);
}

/*------------------------------------*/

 
