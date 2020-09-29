/* segment channa6 in program nc */

/* sets up parameters for channel na type 6 */

/* Model of NaV1.6 channel (SCN8a) from Raman & Bean, 2001 */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncomp.h"
#include "control.h"

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

chantype *makchantype(int ctype, int cnum, int nstates, int nparm, 
					int nq, double bt);
void chanrate(chan *cpnt, chanparm *chp);
double rt(chan *spnt);

/*---------------------------------------------------------------*/

double calcna6m (double v, int func)

/* Calculate na type 5 rate functions alpha, beta, gamma, delta, epsilon, zeta
   from Raman IM and Bean BP (2001) Biophysical J, 80: 729-737.
   */

/* Calculate na6 rate funcs, */
/* Calibrated in mV, per sec */
/* Temperature 22 deg C */

{
   double val;

   if (func > 6) func = 6;
   switch (func) {

     case 1: val = 150. * exp ( v /  20.0);  break; /* alpha */
     case 2: val = 3.   * exp ( v / -20.0);  break; /* beta */
     case 3: val = 150.;                     break; /* gamma */
     case 4: val = 40.;                      break; /* delta */
     case 5: val = 1.75;                     break; /* epsilon */
     case 6: val = 0.03 * exp ( v / -25.0);  break; /* zeta */
     }
   return val * MSSEC;
}

/*--------------------------------------------*/

double calcna6h (double v, int func)

/* Calculate na type 5 rate functions Coff, Con, Ooff, Oon
   from Raman IM and Bean BP (2001) Biophysical J, 80: 729-737.
   */

/* Calculate na6 rate funcs, */
/* Calibrated in mV, per sec */
/* Temperature 22 deg C */

{
   double val;

   if (func > 4) func = 4;
   switch (func) {

     case 1: val = 0.5;    break; /* Coff */
     case 2: val = 0.005;  break; /* Con  */
     case 3: val = 0.005;  break; /* Ooff */
     case 4: val = 0.75;   break; /* Oon  */
     }
   return val * MSSEC;
}

/*--------------------------------------------*/

double na6a(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na6b(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na6g(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na6d(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[3];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na6e(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[4];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na6z(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[5];
}

/*--------------------------------------------*/

double na6c1(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na6c2(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na6o1(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na6o2(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[3];
}


/*----------------------------------------*/

chantype *makna6(void)

/* Discrete-State-Markov description of 
   Na+ channel that very closely mimics 
   resurgent current in mouse purkinje cells.

   from Raman IM and Bean BP (2001) Biophysical J, 80: 729-737.

     (Similar to Na type 3 from Kuo CC and Bean BP (1994) Neuron, 12: 819-829.)
     (  but has extra "blocked" state. )

   Number states from 0 to n; 
   Set numstate = n;
   Set numtrans = number of transitions.
   Set cond     = conductance of state;

   For each state, set transitions:
    Set trans =   State to go to on this transition. 
    Set trate =   Function that returns basic rate for transition. 
    Set ratemul = Multiplier for rate function.
    Set rateo  =  Sets additional tau multiplier:
                        0->am, 1->bm, 2->ah, 3->bh.
*/

/*
	                                                       open     blocked
          4a          3a          2a           a           g         e

     0    <->    1    <->    2    <->    3    <->    4    <->    5  <->  12

           b          2b          3b           4b          d         z 

  C1        C1/A       C1/A^2      C1/A^3      C1/A^4          O1
     |           |           |           |           |           | 
       C2         C2*A        C2*A^2      C2*A^3      C2*A^4      O2  

          4aA         3aA         2aA         aA            g 

     6    <->    7    <->    8    <->    9    <->    10    <->   11

          b/A         2b/A        3b/A        4b/A          d 


C1 = Coff = 0.5e3/sec
C2 = Con  = 0.005e3/sec
O1 = Ooff = 0.005e3/sec
O2 = Oon  = 0.75e3/sec
A  = [(Coff/Con) / (Ooff/Oon) ] ^ 1/8 =   3.3266829 

*/


#define A   3.3266829 
#define AI  0.30059974

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 13;				/* 12 Markov states */
   nparm = 2;				/* make 2 set of seq. state states */
   nq = 6;				/* make 6 Q10 values per parm */
   ch = makchantype(NA, 6, nstate, nparm, nq, dbasetc); /* make chan info */
   ch->unitary = dnau;
   ch->ions->ionp[PNA] = 1.0;			/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = dpkna;		/* permeability to K+ ions / Na perm */
   ch->vrev = vna;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 6;		/* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;		/* don't use implicit vals (yet) */
   parm[0].chancalc =  calcna6m;        /* default rate function */
   parm[0].funcname = (char *)"calcna6m"; /* user rate function */
   parm[0].voff = dnaoffsm;             /* voltage offset from user */

   parm[0].dq[0] = dqm;                 /* Q10 for alpha rate function */
   parm[0].dq[1] = dqm;                 /* Q10 for beta rate function */
   parm[0].dq[2] = dqm;                 /* Q10 for gamma rate function */
   parm[0].dq[3] = dqm;                 /* Q10 for delta rate function */
   parm[0].dq[4] = dqm;                 /* Q10 for epsilon rate function */
   parm[0].dq[5] = dqm;                 /* Q10 for zeta rate function */

   parm[1].nfval    = 4;		/* num of func vals (ah, bh, etc.) */
   parm[1].nival    = 0;		/* don't use implicit vals (yet) */
   parm[1].chancalc =  calcna6h;        /* default rate function */
   parm[1].funcname = (char *)"calcna6h"; /* user rate function */
   parm[1].voff = dnaoffsh;             /* voltage offset from user */

   parm[1].dq[0] = dqh;                 /* Q10 for Coff */
   parm[1].dq[1] = dqh;                 /* Q10 for Con */
   parm[1].dq[2] = dqh;                 /* Q10 for Ooff */
   parm[1].dq[3] = dqh;                 /* Q10 for Oon */
   parm[1].dq[4] = 1.0;
   parm[1].dq[5] = 1.0;

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = na6a;
   spnt[0].ratemul[0] = 4.0 * m;
   spnt[0].rateo  [0] = 0;  /* taum */
   spnt[0].trans  [1] = 6;
   spnt[0].trate  [1] = na6c2;
   spnt[0].ratemul[1] = 1.0 * h;
   spnt[0].rateo  [1] = 2; 	/* tauh */

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = na6a;
   spnt[1].ratemul[0] = 3.0 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = na6b;
   spnt[1].ratemul[1] = 1.0 * m;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 7;	
   spnt[1].trate  [2] = na6c2;
   spnt[1].ratemul[2] = A * h;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = na6a;
   spnt[2].ratemul[0] = 2.0 * m;
   spnt[2].rateo  [0] = 0;	/* taum */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = na6b;
   spnt[2].ratemul[1] = 2.0 * m;
   spnt[2].rateo  [1] = 1;	/* taum */
   spnt[2].trans  [2] = 8;
   spnt[2].trate  [2] = na6c2;
   spnt[2].ratemul[2] = A*A * h;
   spnt[2].rateo  [2] = 2;	/* tauh, tauc */

   spnt[3].numtrans   = 3;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = na6a;
   spnt[3].ratemul[0] = 1.0 * m;
   spnt[3].rateo  [0] = 0;	/* taum */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = na6b;
   spnt[3].ratemul[1] = 3.0 * m;
   spnt[3].rateo  [1] = 1;	/* taum */
   spnt[3].trans  [2] = 9;
   spnt[3].trate  [2] = na6c2;
   spnt[3].ratemul[2] = A*A*A * h;
   spnt[3].rateo  [2] = 2;	/* tauh, tauc */

   spnt[4].numtrans   = 3;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = na6g;
   spnt[4].ratemul[0] = 1.0 * m;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 3;
   spnt[4].trate  [1] = na6b;
   spnt[4].ratemul[1] = 4.0 * m;
   spnt[4].rateo  [1] = 1;
   spnt[4].trans  [2] = 10;
   spnt[4].trate  [2] = na6c2;
   spnt[4].ratemul[2] = A*A*A*A * h;
   spnt[4].rateo  [2] = 2;	/* tauh, tauc */

   spnt[5].numtrans   = 3;
   spnt[5].cond       = 1.0;		/* state 5 = the open state */
   spnt[5].trans  [0] = 11;
   spnt[5].trate  [0] = na6o2;
   spnt[5].ratemul[0] = 1.0 * h;
   spnt[5].rateo  [0] = 2;	/* tauh, tauc */
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = na6d;
   spnt[5].ratemul[1] = 1.0 * m;
   spnt[5].rateo  [1] = 1;
   spnt[5].trans  [2] = 12;		/* goes to blocked state */
   spnt[5].trate  [2] = na6e;
   spnt[5].ratemul[2] = 1.0 * m;
   spnt[5].rateo  [2] = 1;	/* taum */

   spnt[6].numtrans   = 2;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = na6a;
   spnt[6].ratemul[0] = 4.0 * A * m;	
   spnt[6].rateo  [0] = 0;
   spnt[6].trans  [1] = 0;
   spnt[6].trate  [1] = na6c1;
   spnt[6].ratemul[1] = 1.0 * h;
   spnt[6].rateo  [1] = 3;	/* tauh, taud */

   spnt[7].numtrans   = 3;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 8;
   spnt[7].trate  [0] = na6a;
   spnt[7].ratemul[0] = 3.0 * A * m;
   spnt[7].rateo  [0] = 0;
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = na6b;
   spnt[7].ratemul[1] = 1.0 * AI * m;
   spnt[7].rateo  [1] = 1;
   spnt[7].trans  [2] = 1;
   spnt[7].trate  [2] = na6c1;
   spnt[7].ratemul[2] = AI * h;
   spnt[7].rateo  [2] = 3;	/* tauh, taud */

   spnt[8].numtrans   = 3;
   spnt[8].cond       = 0;
   spnt[8].trans  [0] = 9;
   spnt[8].trate  [0] = na6a;
   spnt[8].ratemul[0] = 2.0 * A * m;
   spnt[8].rateo  [0] = 0;
   spnt[8].trans  [1] = 7;
   spnt[8].trate  [1] = na6b;
   spnt[8].ratemul[1] = 2.0 * AI * m;
   spnt[8].rateo  [1] = 1;
   spnt[8].trans  [2] = 2;
   spnt[8].trate  [2] = na6c1;
   spnt[8].ratemul[2] = AI*AI * h;
   spnt[8].rateo  [2] = 3;	/* tauh, taud */

   spnt[9].numtrans   = 3;
   spnt[9].cond       = 0;
   spnt[9].trans  [0] = 10;
   spnt[9].trate  [0] = na6a;
   spnt[9].ratemul[0] = 1.0 * A * m;
   spnt[9].rateo  [0] = 0;
   spnt[9].trans  [1] = 8;
   spnt[9].trate  [1] = na6b;
   spnt[9].ratemul[1] = 3.0 * AI * m;
   spnt[9].rateo  [1] = 1;
   spnt[9].trans  [2] = 3;
   spnt[9].trate  [2] = na6c1;
   spnt[9].ratemul[2] = AI*AI*AI * h;
   spnt[9].rateo  [2] = 3;	/* tauh, taud */

   spnt[10].numtrans   = 3;
   spnt[10].cond       = 0;
   spnt[10].trans  [0] = 11;
   spnt[10].trate  [0] = na6g;
   spnt[10].ratemul[0] = 1.0 * m;
   spnt[10].rateo  [0] = 0;
   spnt[10].trans  [1] = 9;
   spnt[10].trate  [1] = na6b;
   spnt[10].ratemul[1] = 4.0 * AI * m;
   spnt[10].rateo  [1] = 1;
   spnt[10].trans  [2] = 4;
   spnt[10].trate  [2] = na6c1;
   spnt[10].ratemul[2] = AI*AI*AI*AI * h;
   spnt[10].rateo  [2] = 3;	/* tauh, taud */

   spnt[11].numtrans   = 2;
   spnt[11].cond       = 0;
   spnt[11].trans  [0] = 5;
   spnt[11].trate  [0] = na6o1;
   spnt[11].ratemul[0] = 1.0 * h;
   spnt[11].rateo  [0] = 3;	/* tauh, taud */
   spnt[11].trans  [1] = 10;
   spnt[11].trate  [1] = na6d;
   spnt[11].ratemul[1] = 1.0 * m;
   spnt[11].rateo  [1] = 1;

   spnt[12].numtrans   = 1;
   spnt[12].cond       = 0;
   spnt[12].trans  [0] = 5;
   spnt[12].trate  [0] = na6z;
   spnt[12].ratemul[0] = 1.0;
   spnt[12].rateo  [0] = 0;	/* taum, taua */

   return ch;
}
/*----------------------------------------*/


