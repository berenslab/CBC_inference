/* segment initchan in program nc */

/* sets up channel parameters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncomp.h"
#include "control.h"
#include "ncio.h"

static double qfgj;			/* q10 constant for gap junctions */

#define EXPSIZ 1000
#define EXPEXPSIZ 300

double exptab [EXPSIZ+1]={0};    
double expexptab [EXPEXPSIZ+1]={0};    

#define LOG3  1.09861228866

chantype *makchantype(int ctype, int cnum, int nstates, int nparm, 
				int nq, double bt);

#ifdef __cplusplus
extern "C" {
#endif

double exp(double);
double log(double);

#ifdef __cplusplus
}
#endif

chantype *makna0(void);
chantype *makna1(void);
chantype *makna2(void);
chantype *makna3(void);
chantype *makna4(void);
chantype *makna5(void);
chantype *makna6(void);
chantype *makna8(void);
chantype *makna20(void);
chantype *makif(void);  /* type 21 */

chantype *makk0(void);
chantype *makk1(void);
chantype *makk2(void);
chantype *makk3(void);
chantype *makk4(void);
chantype *makk5(void);
chantype *makk6(void);
chantype *makk7(void);
chantype *makk8(void);
chantype *makk9(void);
// begin joesterle
chantype *makk10(void);
chantype *makk11(void);
chantype *makk12(void);
chantype *makk13(void);
chantype *makk14(void);
// end joesterle

chantype *makkca0(void);
chantype *makkca1(void);
chantype *makkca2(void);
chantype *makkca3(void);
chantype *makkca4(void);
chantype *makkca5(void);
chantype *makkca6(void);
chantype *makkca8(void);

chantype *makclca1(void);
chantype *makclca2(void);

chantype *makca0(void);
chantype *makca1(void);
chantype *makca2(void);
chantype *makca3(void);
chantype *makca4(void);
chantype *makca5(void);
chantype *makca6(void);
chantype *makca7(void);

chantype *makcgmp1(void);
chantype *makcgmp2(void);
chantype *makcgmp3(void);
chantype *makcgmp4(void);
chantype *makcgmp5(void);
chantype *makcgmp6(void);
chantype *makcgmp7(void);
chantype *makcgmp8(void);
chantype *makcgmp9(void);
chantype *makcgmp10(void);
chantype *makcgmp11(void);

chantype *makampa1(void);
chantype *makampa2(void);
chantype *makampa3(void);
chantype *makampa4(void);
chantype *makampa5(void);

chantype *maksyn1(void);

chantype *maknmda1(void);
chantype *maknmda2(void);

chantype *makgaba1(void);
chantype *makgaba2(void);
chantype *makgaba3(void);
chantype *makgaba4(void);

chantype *makgly1(void);

chantype *makchr2(void);

double qrate(chanparm *chp);
void getntresp(chan *cpnt, comp *comp1);
void isntresp(chan *cpnt, comp *comp1);
double rnt(chan *cpnt);

double ncabs(double x);

double lkexp(double val);

/*----------------------------------------*/

void makexp(void)
/* set up exponential function lookup table.  */
{
   int i;
   static int done=0;

 if (done) return;

 done = 1;
 for (i=0; i<EXPSIZ+1; i++)
   exptab[i] = exp((double)i/EXPSIZ); 
 for (i=0; i<EXPEXPSIZ+1; i++)
   expexptab[i] = exp((double)(i-EXPEXPSIZ/2)); 

/*  Test lkexp(). It's accurate to factor of 1.25e-7 */

/*  for (i=0; i<100; i++) {
      double scal;

       scal = 0.0001;
   ncfprintf (stdout,
  "%g %20.20g %20.20g %20.20g\n",(double)i*scal,exp((double)i*scal), 
     lkexp((double)i*scal), (exp((double)i*scal) - lkexp((double)i*scal)) /
			       exp((double)i*scal));
 } /* */
}

/*----------------------------------------*/

double lkexp(double val)

/* lookup an exp() value in the range -250 to 250. */

/* On some machines, it's faster to use this 
    interpolation scheme than use the exp() 
    routine supplied with the library.
*/

{
  int i,r;
  double rr,x;
  
  if      (val > 250.0) val = 250.0; 
  else if (val <-250.0) val = -250.0; 
  i = int(val);				/* integer for lookup in expexptab[] */
  if (val<0) i--;
  rr = (val - i) * EXPSIZ;		/* remainder for lookup in exptab[] */
  r = int(rr + 0.5);
  rr -=  r;
  x = exptab[r] * (1.0+(1.0/EXPSIZ)*rr); /* use linear interpolation */
  return (x*expexptab[i+EXPEXPSIZ/2]); 
}

/*--------------------------------------------*/

double qrategj (void)

/* Return rate constant multiplier for gap junctions 
   which depends on temperature. Computed at a constant time step.
   Assume normal temperature of 22 deg C (just for simplicity).
*/

{
 return stiminc * exp(LOG3 * (tempcel - 22.0) / 10.0); /* Q10 of 3 */
}

/*--------------------------------------------*/

#define GJLAMBDA 0.0013
#define GJVOFF 0.0147
#define Aa 0.077
#define Ab 0.14

double gjalpha (double v, gj *gjpnt)

/* Rate function for gap junction conductance change */

/* v is in mV, return rate in 1/sec */

/* For details, see Harris, Spray, Bennett, */
/* J.  Gen Physiol. 77:95-117, 1981.  */

{
    double vj;                          /* junction voltage */
    double alpha, n;

   if (v<0 && (gjpnt->rev)) {		/* reverse parameters */
     vj = ncabs(v);
     alpha = qfgj * gjpnt->rtaun * GJLAMBDA * 
		exp(-Aa * (vj - gjpnt->rvoff)*VTOMV/gjpnt->rvgain);
   }
   else {				/* forward parameters */
     vj = ncabs(v);
     alpha = qfgj * gjpnt->taun * GJLAMBDA * 
		exp(-Aa * (vj - gjpnt->voff)*VTOMV/gjpnt->vgain);
   }
   return alpha*MSSEC;
}

/*--------------------------------------------*/

double gjbeta (double v, gj *gjpnt)

/* Rate function for gap junction conductance change */

/* v is in mV */

/* For details, see Harris, Spray, Bennett, */
/* J.  Gen Physiol. 77:95-117, 1981.  */

{
    double vj;                          /* junction voltage */
    double beta, n;

   if (v<0 && (gjpnt->rev)) {		/* reverse parameters */
     vj = ncabs(v);
     beta  = gjpnt->rtaun * GJLAMBDA * exp( Ab * (vj - gjpnt->rvoff) *
		VTOMV/gjpnt->rvgain);
   }
   else {				/* forward parameters */
     vj = ncabs(v);
     beta  = gjpnt->taun * GJLAMBDA * exp( Ab * (vj - gjpnt->voff) *
		VTOMV/gjpnt->vgain);
   }
   beta = qfgj * beta / (1.0 + 50 * beta); 
   return beta*MSSEC;
}

/*----------------------------------------*/
 
void initchan(void)

/* This subroutine sets up all sequential-state channel tables
    at the beginning of a "run".  To set up additional types,
    copy the "makna1()" or "makk1()" routine and modify it, then
    call it from this routine.  
*/

{
   static int done = 0;

 if (done) return;
 else done = 1;

 makexp();			/* set up exp() table for synapses */
 qfgj = qrategj();

 makna0();			/* set up Na tables */
 makna1();			/* set up Na sequential states */
 makna2();
 makna3();
 makna4();
 makna5();
 makna6();
 makna8();

 makna20();
 makif ();			/* set int-and-fire states=Na type 21 */


 makk0();			/* set up Khh */
 makk1();			/* set up K 1 sequential states */
 makk2();			/* set up KA hh  sequential states */
 makk3();			/* set up KA markov */
 makk4();			/* set up K 4 sequential states */
 makk5();			/* set up K 5 sequential states */
 makk6();			/* set up K 6 sequential states */
 makk7();			/* set up K 7 sequential states */
 //makk8();			/* set up K 8 sequential states */
 makk9();			/* set up K 9 sequential states */
 // begin joesterle
 makk10();		/* set up K 10 sequential states */
 makk11();		/* set up K 10 sequential states */
 makk12();		/* set up K 10 sequential states */
 makk13();		/* set up K 10 sequential states */
 makk14();		/* set up K 10 sequential states */
 // end joesterle

 makkca0();			/* set up Kca sequential states */
 makkca1();			/* set up Kca sequential states */
 makkca2();			/* set up Kca sequential states */
 makkca3();			/* set up Kca sequential states */
 makkca4();			/* set up Kca sequential states */
 makkca5();			/* set up Kca sequential states */
 makkca6();			/* set up Kca sequential states */
 makkca8();			/* set up Kca sequential states */

 makclca1();			/* set up Clca sequential states */
 makclca2();			/* set up Clca sequential states */

 makca0();			/* set up Ca tables states */
 makca1();
 makca2();
 makca3();
 makca4();
 makca5();
 makca6();
 makca7();

 makampa1();
 makampa2();
 makampa3();
 makampa4();
 makampa5();

 makcgmp1();
 makcgmp2();
 makcgmp3();
 makcgmp4();
 makcgmp5();
 makcgmp6();
 makcgmp7();
 makcgmp8();
 makcgmp9();
 makcgmp10();
 makcgmp11();

 maknmda1();
 maknmda2();

 makgaba1();
 makgaba2();
 makgaba3();
 makgaba4();
 
 makgly1();
 
 makchr2();

 maksyn1();			/* 2 state markov synapse from conduct */
}


