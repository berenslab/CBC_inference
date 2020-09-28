/* Program bnldev */

/* Finds binomial distribution */

/* prints distribution on standard output */

/* Latest Mod	Apr 84		R.G.Smith */

#define PI 3.141592653589793

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

#include <stdlib.h>

static int N=5;
static int n=10;
static double p=.5;
int cumrand = 0;

FILE *pictin;

double atof();
double log();
double exp();
void run(void);

/****************************************/

main(int argc, char **argv)

{
   char *cptr;
   int i;
   FILE *temp,*freopen();
	 
 pictin = stdin;
 if (argc==1)			/* if user needs help */
   {
    printf ("\nBinomial deviation program Dec 94\r\n");
    run();
   }
 else
 do					/* if there are any switches */
  {
   argc--; argv++;
   cptr = *argv;
   if (argc)
    {
     if (*cptr == '-')
      {
       cptr++;
       switch (*cptr)
        {

      	  case 'N': 
		argv++; argc--;
		N = atoi(*argv);
		break;

      	  case 'n': 
		argv++; argc--;
		n = atoi(*argv);
		break;

      	  case 'p': 
		argv++; argc--;
		p = atof(*argv);
		break;

        }  /* switch */
      }	   /* if */
     else
      {
       if((pictin=fopen(cptr,"r"))==NULL)
         {
           fprintf(stderr,"Binom: cannot open %s\n",cptr);
	   fflush (stderr);
           continue;
         }
       run();
       if (argc <= 1) break;
     }
    }
   else run();
  }
 while (argc > 0);
}

/*------------------------------------------*/

void run(void)
{

   int i;
   double val, tot;
   double bnldev(double pp, int n);

  tot = 0.0;  
  for (i=0; i<n; i++) {
   val = bnldev(p,N);
   printf ("%g\n",val);

  }
/* printf ("avg %g p %g n %d\n",tot/k,p,n); /* */
}

/*------------------------------------------*/

double bnldev(double pp, int n)

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
   double drand(void), gammln(double xx); 

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
	   if (angle == 0.5) angle = 0.4999999999;
	   y = tan(PI*angle);
	   em = sq*y+am;
	} while (em < 0.0 || em >= (en+1.0));	/* reject */
	em = floor(em);
	t = 1.2 * sq * (1.0+y*y) * exp(oldg-gammln(em+1.0)
		- gammln(en-em+1) + em * plog+ (en-em) * pclog);

      } while (drand() > t);		/* reject. This happens about 1.5 */
      bnl = em;				/* times per deviate, on average. */
   }
   if (p != pp) bnl = n -bnl;
   return bnl;
}

/*------------------------------------*/

double gammln(double xx)
{
   double x,y,tmp,ser;
   static double
   cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,
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

