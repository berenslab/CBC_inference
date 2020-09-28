/* segment channa3 in program nc */

/* sets up channel na 3 parameters */

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

double calcna3m (double v, int func)

/* Calculate na type 3 rate functions alpha, beta, gamma, delta
   from Kuo CC and Bean BP (1994) Neuron, 12: 819-829.
   */

/* Calculate na3 rate funcs, */
/* Calibrated in mV, per sec */
/* Temperature 22 deg C */

{
   double val;

   if (func > 4) func = 4;
   switch (func) {

     case 1: val = 20. * exp ( v /  40.0);  break; /* alpha */
     case 2: val = 3.  * exp ( v / -40.0);  break; /* beta */
     case 3: val = 50. * exp ( v / 100.0);  break; /* gamma */
     case 4: val = 0.2 * exp ( v / -25.0);  break; /* delta */
     }
   return val * MSSEC;
}

/*--------------------------------------------*/

double calcna3h (double v, int func)

/* Calculate na type 3 rate functions Coff, Con, Ooff, Oon
   from Kuo CC and Bean BP (1994) Neuron, 12: 819-829.
   */

/* Calculate na5 rate funcs, */
/* Calibrated in mV, per sec */
/* Temperature 22 deg C */

{
   double val;

   if (func > 4) func = 4;
   switch (func) {

     case 1: val = 4.5;    break; /* Coff */
     case 2: val = 0.004;  break; /* Con  */
     case 3: val = 0.008;  break; /* Ooff */
     case 4: val = 4.0;    break; /* Oon  */
     }
   return val * MSSEC;
}

/*--------------------------------------------*/

double na3a(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na3b(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na3g(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na3d(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[3];
}

/*--------------------------------------------*/

double na3c1(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na3c2(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na3o1(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na3o2(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[3];
}

/*----------------------------------------*/

chantype *makna3(void)

/* Discrete-State-Markov description of 
   Na+ channel that very closely mimics 
   kinetics of patch clamp data from squid giant axon.

   See Kuo CC and Bean BP (1994) Neuron, 12: 819-829.

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
	                                                           open 
              4a          3a          2a           a           g      

         0    <->    1    <->    2    <->    3    <->    4    <->    5 

               b          2b          3b           4b          d 

      C1        C1/A       C1/A^2      C1/A^3      C1/A^4          O1
         |           |           |           |           |           | 
           C2         C2*A        C2*A^2      C2*A^3      C2*A^4      O2  

              4aA         3aA         2aA         aA            g 

         6    <->    7    <->    8    <->    9    <->    10    <->   11

              b/A         2b/A        3b/A        4b/A          d 


C1 = Coff = 4.5e3/sec
C2 = Con  = 0.004e3/sec
O1 = Ooff = 0.008e3/sec
O2 = Oon  = 4e3/sec
A = [(Coff/Con) / (Ooff/Oon) ] ^ 1/8 =  5.2331757

*/

#define A  5.2331757
#define AI 0.19108856

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 12;				/* 12 Markov states */
   nparm = 2;				/* make 2 sets of seq. state states */
   nq = 4;				/* make 4 Q10 value per parm */
   ch = makchantype(NA,3, nstate, nparm, nq, dbasetc); /* make chan info */
   ch->unitary = dnau;
   ch->ions->ionp[PNA] = 1.0;			/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = dpkna;		/* permeability to K+ ions / Na perm */
   ch->vrev = vna;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 4;		/* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;		/* don't use implicit vals (yet) */
   parm[0].chancalc =  calcna3m;        /* default rate function */
   parm[0].funcname = (char *)"calcna3m"; /* user rate function */
   parm[0].voff = dnaoffsm;             /* voltage offset from user */

   parm[0].dq[0] = dqm;                 /* Q10 for alpha rate function */
   parm[0].dq[1] = dqm;                 /* Q10 for beta rate function */
   parm[0].dq[2] = dqm;                 /* Q10 for gamma rate function */
   parm[0].dq[3] = dqm;                 /* Q10 for delta rate function */

   parm[1].nfval    = 4;                /* num of func vals (ah, bh, etc.) */
   parm[1].nival    = 0;                /* don't use implicit vals (yet) */
   parm[1].chancalc =  calcna3h;        /* default rate function */
   parm[1].funcname = (char *)"calcna3h"; /* user rate function */
   parm[1].voff = dnaoffsh;             /* voltage offset from user */

   parm[1].dq[0] = dqh;                 /* Q10 for Coff */
   parm[1].dq[1] = dqh;                 /* Q10 for Con */
   parm[1].dq[2] = dqh;                 /* Q10 for Ooff */
   parm[1].dq[3] = dqh;                 /* Q10 for Oon */

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = na3a;
   spnt[0].ratemul[0] = 4.0 * m;
   spnt[0].rateo  [0] = 0;  /* taum */
   spnt[0].trans  [1] = 6;
   spnt[0].trate  [1] = na3c2;
   spnt[0].ratemul[1] = 1.0 * h;
   spnt[0].rateo  [1] = 2; 

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = na3a;
   spnt[1].ratemul[0] = 3.0 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = na3b;
   spnt[1].ratemul[1] = 1.0 * m;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 7;	
   spnt[1].trate  [2] = na3c2;
   spnt[1].ratemul[2] = A * h;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = na3a;
   spnt[2].ratemul[0] = 2.0 * m;
   spnt[2].rateo  [0] = 0;	/* taum */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = na3b;
   spnt[2].ratemul[1] = 2.0 * m;
   spnt[2].rateo  [1] = 1;	/* taum */
   spnt[2].trans  [2] = 8;
   spnt[2].trate  [2] = na3c2;
   spnt[2].ratemul[2] = A*A * h;
   spnt[2].rateo  [2] = 2;	/* taum */

   spnt[3].numtrans   = 3;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = na3a;
   spnt[3].ratemul[0] = 1.0 * m;
   spnt[3].rateo  [0] = 0;	/* taum */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = na3b;
   spnt[3].ratemul[1] = 3.0 * m;
   spnt[3].rateo  [1] = 1;	/* taum */
   spnt[3].trans  [2] = 9;
   spnt[3].trate  [2] = na3c2;
   spnt[3].ratemul[2] = A*A*A * h;
   spnt[3].rateo  [2] = 2;	/* tauh, tauc */

   spnt[4].numtrans   = 3;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = na3g;
   spnt[4].ratemul[0] = 1.0 * m;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 3;
   spnt[4].trate  [1] = na3b;
   spnt[4].ratemul[1] = 4.0 * m;
   spnt[4].rateo  [1] = 1;
   spnt[4].trans  [2] = 10;
   spnt[4].trate  [2] = na3c2;
   spnt[4].ratemul[2] = A*A*A*A * h;
   spnt[4].rateo  [2] = 2;	/* tauh, tauc */

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 1.0;		/* state 5 = the open state */
   spnt[5].trans  [0] = 11;
   spnt[5].trate  [0] = na3o2;
   spnt[5].ratemul[0] = 1.0 * h;
   spnt[5].rateo  [0] = 2;
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = na3d;
   spnt[5].ratemul[1] = 1.0 * m;
   spnt[5].rateo  [1] = 1;

   spnt[6].numtrans   = 2;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = na3a;
   spnt[6].ratemul[0] = 4.0 * A * m;	
   spnt[6].rateo  [0] = 0;
   spnt[6].trans  [1] = 0;
   spnt[6].trate  [1] = na3c1;
   spnt[6].ratemul[1] = 1.0 * h;
   spnt[6].rateo  [1] = 3;	/* tauh, taud */

   spnt[7].numtrans   = 3;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 8;
   spnt[7].trate  [0] = na3a;
   spnt[7].ratemul[0] = 3.0 * A * m;
   spnt[7].rateo  [0] = 0;
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = na3b;
   spnt[7].ratemul[1] = 1.0 * AI * m;
   spnt[7].rateo  [1] = 1;
   spnt[7].trans  [2] = 1;
   spnt[7].trate  [2] = na3c1;
   spnt[7].ratemul[2] = AI * h;
   spnt[7].rateo  [2] = 3;	/* tauh, taud */

   spnt[8].numtrans   = 3;
   spnt[8].cond       = 0;
   spnt[8].trans  [0] = 9;
   spnt[8].trate  [0] = na3a;
   spnt[8].ratemul[0] = 2.0 * A * m;
   spnt[8].rateo  [0] = 0;
   spnt[8].trans  [1] = 7;
   spnt[8].trate  [1] = na3b;
   spnt[8].ratemul[1] = 2.0 * AI * m;
   spnt[8].rateo  [1] = 1;
   spnt[8].trans  [2] = 2;
   spnt[8].trate  [2] = na3c1;
   spnt[8].ratemul[2] = AI*AI * h;
   spnt[8].rateo  [2] = 3;	/* tauh, taud */

   spnt[9].numtrans   = 3;
   spnt[9].cond       = 0;
   spnt[9].trans  [0] = 10;
   spnt[9].trate  [0] = na3a;
   spnt[9].ratemul[0] = 1.0 * A * m;
   spnt[9].rateo  [0] = 0;
   spnt[9].trans  [1] = 8;
   spnt[9].trate  [1] = na3b;
   spnt[9].ratemul[1] = 3.0 * AI * m;
   spnt[9].rateo  [1] = 1;
   spnt[9].trans  [2] = 3;
   spnt[9].trate  [2] = na3c1;
   spnt[9].ratemul[2] = AI*AI*AI * h;
   spnt[9].rateo  [2] = 3;	/* tauh, taud */

   spnt[10].numtrans   = 3;
   spnt[10].cond       = 0;
   spnt[10].trans  [0] = 11;
   spnt[10].trate  [0] = na3g;
   spnt[10].ratemul[0] = 1.0 * m;
   spnt[10].rateo  [0] = 0;
   spnt[10].trans  [1] = 9;
   spnt[10].trate  [1] = na3b;
   spnt[10].ratemul[1] = 4.0 * AI * m;
   spnt[10].rateo  [1] = 1;
   spnt[10].trans  [2] = 4;
   spnt[10].trate  [2] = na3c1;
   spnt[10].ratemul[2] = AI*AI*AI*AI * h;
   spnt[10].rateo  [2] = 3;	/* tauh, taud */

   spnt[11].numtrans   = 2;
   spnt[11].cond       = 0;
   spnt[11].trans  [0] = 5;
   spnt[11].trate  [0] = na3o1;
   spnt[11].ratemul[0] = 1.0 * h;
   spnt[11].rateo  [0] = 3;	/* tauh, taud */
   spnt[11].trans  [1] = 10;
   spnt[11].trate  [1] = na3d;
   spnt[11].ratemul[1] = 1.0 * m;
   spnt[11].rateo  [1] = 1;

   return ch;
}
#undef A
#undef AI

/*---------------------------------------------------------------*/

double calcna4m (double v, int func)

/* Calculate Na type 4 rate functions alpha, beta, gamma, delta
   from Taddese A and Bean BP (2002) Neuron, 33: 587-600.
   */

/* Calculate na4 rate funcs, */
/* Calibrated in mV, per sec */
/* Temperature 22 deg C */

{
   double val;

   if (func > 4) func = 4;
   switch (func) {

     case 1: val = 140. * exp ( v /  27.0);  break; /* alpha */
     case 2: val = 10.  * exp ( v / -27.0);  break; /* beta */
     case 3: val = 150.;  break;		    /* gamma */
     case 4: val = 40.; break;			    /* delta */
     }
   return val * MSSEC;
}

/*--------------------------------------------*/

double calcna4h (double v, int func)

/* Calculate Na type 4 rate functions Coff, Con, Ooff, Oon
   from Taddese A and Bean BP (2002) Neuron, 33: 587-600.
   */

/* Calculate na5 rate funcs, */
/* Calibrated in mV, per sec */
/* Temperature 22 deg C */

{
   double val;

   if (func > 4) func = 4;
   switch (func) {

     case 1: val = 9;      break; /* Coff */
     case 2: val = 0.004;  break; /* Con  */
     case 3: val = 0.008;  break; /* Ooff */
     case 4: val = 0.8;    break; /* Oon  */
     }
   return val * MSSEC;
}

/*--------------------------------------------*/

chantype *makna4(void)

/* Discrete-State-Markov description of 
   Na2v Na+ channel that gives a persistent current.

   See Taddese A and Bean BP (2002) Neuron, 33: 587-600.

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
	                                                           open 
              4a          3a          2a           a           g      

         0    <->    1    <->    2    <->    3    <->    4    <->    5 

               b          2b          3b           4b          d 

      C1        C1/B       C1/B^2      C1/B^3      C1/B^4          O1
         |           |           |           |           |           | 
           C2         C2*A        C2*A^2      C2*A^3      C2*A^4      O2  

              4aA         3aA         2aA         aA            g 

         6    <->    7    <->    8    <->    9    <->    10    <->   11

              b/B         2b/B        3b/B        4b/B          d 


C1 = Coff =     9e3/sec
C2 = Con  = 0.004e3/sec
O1 = Ooff = 0.004e3/sec
O2 = Oon  =   0.8e3/sec
A = [(.05/0.004) ] ^ 1/4 =  1.88
B = [(9/0.00025) ] ^ 1/4 =  13.8

*/

#define A  1.88
#define B  13.8
#define BI 0.072463768

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 12;				/* 12 Markov states */
   nparm = 2;				/* make 2 sets of seq. state states */
   nq = 4;				/* make 4 Q10 value per parm */
   ch = makchantype(NA,4, nstate, nparm, nq, dbasetc); /* make chan info */
   ch->unitary = dnau;
   ch->ions->ionp[PNA] = 1.0;			/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = dpkna;		/* permeability to K+ ions / Na perm */
   ch->vrev = vna;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 4;		/* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;		/* don't use implicit vals (yet) */
   parm[0].chancalc =  calcna4m;        /* default rate function */
   parm[0].funcname = (char *)"calcna4m"; /* user rate function */
   parm[0].voff = dnaoffsm;             /* voltage offset from user */

   parm[0].dq[0] = dqm;                 /* Q10 for alpha rate function */
   parm[0].dq[1] = dqm;                 /* Q10 for beta rate function */
   parm[0].dq[2] = dqm;                 /* Q10 for gamma rate function */
   parm[0].dq[3] = dqm;                 /* Q10 for delta rate function */

   parm[1].nfval    = 4;                /* num of func vals (ah, bh, etc.) */
   parm[1].nival    = 0;                /* don't use implicit vals (yet) */
   parm[1].chancalc =  calcna4h;        /* default rate function */
   parm[1].funcname = (char *)"calcna4h"; /* user rate function */
   parm[1].voff = dnaoffsh;             /* voltage offset from user */

   parm[1].dq[0] = dqh;                 /* Q10 for Coff */
   parm[1].dq[1] = dqh;                 /* Q10 for Con */
   parm[1].dq[2] = dqh;                 /* Q10 for Ooff */
   parm[1].dq[3] = dqh;                 /* Q10 for Oon */

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = na3a;
   spnt[0].ratemul[0] = 4.0 * m;
   spnt[0].rateo  [0] = 0;  /* taum */
   spnt[0].trans  [1] = 6;
   spnt[0].trate  [1] = na3c2;
   spnt[0].ratemul[1] = 1.0 * h;
   spnt[0].rateo  [1] = 2; 

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = na3a;
   spnt[1].ratemul[0] = 3.0 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = na3b;
   spnt[1].ratemul[1] = 1.0 * m;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 7;	
   spnt[1].trate  [2] = na3c2;
   spnt[1].ratemul[2] = A * h;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = na3a;
   spnt[2].ratemul[0] = 2.0 * m;
   spnt[2].rateo  [0] = 0;	/* taum */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = na3b;
   spnt[2].ratemul[1] = 2.0 * m;
   spnt[2].rateo  [1] = 1;	/* taum */
   spnt[2].trans  [2] = 8;
   spnt[2].trate  [2] = na3c2;
   spnt[2].ratemul[2] = A*A * h;
   spnt[2].rateo  [2] = 2;	/* taum */

   spnt[3].numtrans   = 3;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = na3a;
   spnt[3].ratemul[0] = 1.0 * m;
   spnt[3].rateo  [0] = 0;	/* taum */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = na3b;
   spnt[3].ratemul[1] = 3.0 * m;
   spnt[3].rateo  [1] = 1;	/* taum */
   spnt[3].trans  [2] = 9;
   spnt[3].trate  [2] = na3c2;
   spnt[3].ratemul[2] = A*A*A * h;
   spnt[3].rateo  [2] = 2;	/* tauh, tauc */

   spnt[4].numtrans   = 3;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = na3g;
   spnt[4].ratemul[0] = 1.0 * m;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 3;
   spnt[4].trate  [1] = na3b;
   spnt[4].ratemul[1] = 4.0 * m;
   spnt[4].rateo  [1] = 1;
   spnt[4].trans  [2] = 10;
   spnt[4].trate  [2] = na3c2;
   spnt[4].ratemul[2] = A*A*A*A * h;
   spnt[4].rateo  [2] = 2;	/* tauh, tauc */

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 1.0;		/* state 5 = the open state */
   spnt[5].trans  [0] = 11;
   spnt[5].trate  [0] = na3o2;
   spnt[5].ratemul[0] = 1.0 * h;
   spnt[5].rateo  [0] = 2;
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = na3d;
   spnt[5].ratemul[1] = 1.0 * m;
   spnt[5].rateo  [1] = 1;

   spnt[6].numtrans   = 2;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = na3a;
   spnt[6].ratemul[0] = 4.0 * A * m;	
   spnt[6].rateo  [0] = 0;
   spnt[6].trans  [1] = 0;
   spnt[6].trate  [1] = na3c1;
   spnt[6].ratemul[1] = 1.0 * h;
   spnt[6].rateo  [1] = 3;	/* tauh, taud */

   spnt[7].numtrans   = 3;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 8;
   spnt[7].trate  [0] = na3a;
   spnt[7].ratemul[0] = 3.0 * A * m;
   spnt[7].rateo  [0] = 0;
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = na3b;
   spnt[7].ratemul[1] = 1.0 * BI * m;
   spnt[7].rateo  [1] = 1;
   spnt[7].trans  [2] = 1;
   spnt[7].trate  [2] = na3c1;
   spnt[7].ratemul[2] = BI * h;
   spnt[7].rateo  [2] = 3;	/* tauh, taud */

   spnt[8].numtrans   = 3;
   spnt[8].cond       = 0;
   spnt[8].trans  [0] = 9;
   spnt[8].trate  [0] = na3a;
   spnt[8].ratemul[0] = 2.0 * A * m;
   spnt[8].rateo  [0] = 0;
   spnt[8].trans  [1] = 7;
   spnt[8].trate  [1] = na3b;
   spnt[8].ratemul[1] = 2.0 * BI * m;
   spnt[8].rateo  [1] = 1;
   spnt[8].trans  [2] = 2;
   spnt[8].trate  [2] = na3c1;
   spnt[8].ratemul[2] = BI*BI * h;
   spnt[8].rateo  [2] = 3;	/* tauh, taud */

   spnt[9].numtrans   = 3;
   spnt[9].cond       = 0;
   spnt[9].trans  [0] = 10;
   spnt[9].trate  [0] = na3a;
   spnt[9].ratemul[0] = 1.0 * A * m;
   spnt[9].rateo  [0] = 0;
   spnt[9].trans  [1] = 8;
   spnt[9].trate  [1] = na3b;
   spnt[9].ratemul[1] = 3.0 * BI * m;
   spnt[9].rateo  [1] = 1;
   spnt[9].trans  [2] = 3;
   spnt[9].trate  [2] = na3c1;
   spnt[9].ratemul[2] = BI*BI*BI * h;
   spnt[9].rateo  [2] = 3;	/* tauh, taud */

   spnt[10].numtrans   = 3;
   spnt[10].cond       = 0;
   spnt[10].trans  [0] = 11;
   spnt[10].trate  [0] = na3g;
   spnt[10].ratemul[0] = 1.0 * m;
   spnt[10].rateo  [0] = 0;
   spnt[10].trans  [1] = 9;
   spnt[10].trate  [1] = na3b;
   spnt[10].ratemul[1] = 4.0 * BI * m;
   spnt[10].rateo  [1] = 1;
   spnt[10].trans  [2] = 4;
   spnt[10].trate  [2] = na3c1;
   spnt[10].ratemul[2] = BI*BI*BI*BI * h;
   spnt[10].rateo  [2] = 3;	/* tauh, taud */

   spnt[11].numtrans   = 2;
   spnt[11].cond       = 0;
   spnt[11].trans  [0] = 5;
   spnt[11].trate  [0] = na3o1;
   spnt[11].ratemul[0] = 1.0 * h;
   spnt[11].rateo  [0] = 3;	/* tauh, taud */
   spnt[11].trans  [1] = 10;
   spnt[11].trate  [1] = na3d;
   spnt[11].ratemul[1] = 1.0 * m;
   spnt[11].rateo  [1] = 1;

   return ch;
}

#undef A
#undef B
#undef BI

/*----------------------------------------*/


