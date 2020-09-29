/* segment channa5 in program nc */

/* sets up parameters for channel na type 5 */

/* Model of NaV1.1 channel (SCN1A) from */
/* Clancy C and Kass RS (2004) Biophysical J, 86: 2606-2614. */

/* This channel has a slower recovery from inactivation than the
   standard Hodgkin-Huxley Na channel. */

/*  Also included in that paper is a closely related mutated
    channel that is thought to be involved in epilepsy. */

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

double calcna5m (double v, int func)

/* Calculate na type 5 (Nav1.1 or SCN1A) rate constant functions 
     from:  Clancy C and Kass RS (2004) Biophysical J, 86: 2606-2614.
             Rate functions from original source code 
             provided by Colleen E Clancy, Columbia University.
   */

/* Calculate na5 rate funcs, */
/* Calibrated in mV, per sec */
/* Temperature 22 deg C */

{
   double val;

   if (func > 6) func = 6;
   switch (func) {

     case 1: val = 2.802/(0.21*exp(v/-17.0)+0.23*exp(v/-150)); break; /* a11 */
     case 2: val = 2.802/(0.23*exp(v/-15.0)+0.25*exp(v/-150)); break; /* a12 */
     case 3: val = 2.802/(0.25*exp(v/-12.0)+0.27*exp(v/-150)); break; /* a13 */
     case 4: val = 0.4*exp(v/-20.3);			       break; /* b11 */
     case 5: val = 0.4*exp((v-5)/-20.3);		       break; /* b12 */
     case 6: val = 0.4*exp((v-10)/-20.3)/4.5;		       break; /* b13 */
     }
   return val * MSSEC;
}

/*--------------------------------------------*/

double na5_a2 (double v)
{
    return ((9.178*exp(v/29.68))/4.5);
}

/*--------------------------------------------*/

double na5_a3 (double v)
{
    return (3.7933e-7*exp(v/-7.6))*3; 
}
/*--------------------------------------------*/

double na5_a13 (double v)
{
     return (2.802/(0.25*exp(v/-12.0)+0.27*exp(v/-150)));
}
 
/*--------------------------------------------*/

double na5_b3 (double v)
{
     return (0.0084+.00002*v);
}

/*--------------------------------------------*/

double na5_b13 (double v)
{
     return (0.4*exp((v-10)/-20.3)/4.5);
}
 
/*--------------------------------------------*/

double calcna5h (double v, int func)

/* Calculate na type 5 rate functions a2,a3 and b2,b3
   from  Clancy C and Kass RS (2004) Biophysical J, 86: 2606-2614.
   Rate functions from original source code 
   provided by Colleen E Clancy, Columbia University.
   */

/* Calculate na5 rate funcs, */
/* Calibrated in mV, per sec */
/* Temperature 22 deg C */

{
   double val;

   if (func > 4) func = 4;
   switch (func) {
 
     case 1: val = ((9.178*exp(v/29.68))/4.5);					break; /* a2 */
     case 2: val = (3.7933e-7*exp(v/-7.6))*3;   				break; /* a3 */
     case 3: val = ((na5_a13(v)*na5_a2(v)*na5_a3(v))/(na5_b13(v)*na5_b3(v)));	break; /* b2 */
     case 4: val = (0.0084 + 0.00002*v);					break; /* b3 */

     }
   return val * MSSEC;
}

/*--------------------------------------------*/

double calcna5i (double v, int func)

/* Calculate na type 5 slow inactivating rate functions a4,a5 and b4,b5
   from  Clancy C and Kass RS (2004) Biophysical J, 86: 2606-2614.
   Rate functions from original source code 
   provided by Colleen E Clancy, Columbia University.
   */

/* Calculate na5 rate funcs, */
/* Calibrated in mV, per sec */
/* Temperature 22 deg C */

{
   double val;

   if (func > 4) func = 4;
   switch (func) {
 
     case 1: val = (na5_a2(v)/100)*1.5;		break; /* a4 */
     case 2: val = (na5_a2(v)/95000)*80.0;	break; /* a5 */
     case 3: val = na5_a3(v)/5;			break; /* b4 */
     case 4: val = (na5_a3(v)/30)/10.0;		break; /* b5 */

     }
   return val * MSSEC;
}

/*--------------------------------------------*/

double na5a11(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na5a12(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na5a13(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na5b11(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[3];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na5b12(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[4];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na5b13(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[5];
}

/*--------------------------------------------*/

double na5a2(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na5a3(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na5b2(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na5b3(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[3];
}


/*--------------------------------------------*/

double na5a4(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[2];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na5a5(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[2];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na5b4(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[2];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double na5b5(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[2];
  chanrate(cpnt, chp);
  return chp->fval[3];
}


/*----------------------------------------*/

chantype *makna5(void)

/* Discrete-State-Markov description of SCN1A (Nav1.1) channel.

   Markov diagram from  
      Clancy C and Kass RS (2004) Biophysical J, 86: 2606-2614.
   original recording from  
      Lossin C, Wang DW, Rhodes TH, Vanoye CG, George AL (2002) Neuron 34: 877-884.

     (Similar to cardiac Na  channel from Clancy & Rudy, 1999.)

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
State equivalents from Clancy source code
to numbering scheme used here:

UC3  = Na_MC3  -> state 0
UC2  = Na_MC2  -> state 1
UC1  = Na_MC   -> state 2
UO   = Na_MO   -> state 3
UIC3 = IC3     -> state 4
UIC2 = IC2     -> state 5
UIF  = NA_MI   -> state 6
UIM1 = Na_MIs  -> state 7
UIM2 = Na_MIs2 -> state 8

	                                  open
          a11 >       a12 >        a13 >

    Na_MC3      Na_MC2      Na_MC         Na_MO
    UC3         UC2         UC1           UO

     0    <->    1    <->    2     <->     3   

         < b11      < b12        < b13      

  ^           ^           ^        ^
  a3 | b3     a3 | b3     a3 | b3  b2 / a2 (sic)


          a11 >       a12 >       a4 >        a5 >

    IC3         IC2         Na_MI       Na_MIs      NA_MIs2
    UIC3        UIC2        UIF         UIM1        UIM2
     4    <->    5    <->    6    <->    7    <->    8            inactivated states

        < b11       < b12       < b4         < b5


*/

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 9;				/* 9 Markov states */
   nparm = 3;				/* make 3 set of seq. state states */
   nq = 6;				/* make 6 Q10 values per parm */
   ch = makchantype(NA, 5, nstate, nparm, nq, dbasetc); /* make chan info */
   ch->unitary = dnau;
   ch->ions->ionp[PNA] = 1.0;			/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = dpkna;		/* permeability to K+ ions / Na perm */
   ch->vrev = vna;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 6;		/* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;		/* don't use implicit vals (yet) */
   parm[0].chancalc =  calcna5m;        /* default rate function */
   parm[0].funcname = (char *)"calcna5m";       /* user rate function */
   parm[0].voff = dnaoffsm;             /* voltage offset from user */

   parm[0].dq[0] = dqm;                 /* Q10 for a11 */
   parm[0].dq[1] = dqm;                 /* Q10 for a12 */
   parm[0].dq[2] = dqm;                 /* Q10 for a13 */
   parm[0].dq[3] = dqm;                 /* Q10 for b13 */
   parm[0].dq[4] = dqm;                 /* Q10 for b14 */
   parm[0].dq[5] = dqm;                 /* Q10 for b15 */

   parm[1].nfval    = 4;		/* num of func vals (ah, bh, etc.) */
   parm[1].nival    = 0;		/* don't use implicit vals (yet) */
   parm[1].chancalc =  calcna5h;        /* default rate function */
   parm[1].funcname = (char *)"calcna5h";       /* user rate function */
   parm[1].voff = dnaoffsh;             /* voltage offset from user */

   parm[1].dq[0] = dqh;                 /* Q10 for a2 */
   parm[1].dq[1] = dqh;                 /* Q10 for b2 */
   parm[1].dq[2] = dqh;                 /* Q10 for a3 */
   parm[1].dq[3] = dqh;                 /* Q10 for b3 */
   parm[1].dq[4] = 1.0;
   parm[1].dq[5] = 1.0;

   parm[2].nfval    = 4;		/* num of func vals (am, bm, etc.) */
   parm[2].nival    = 0;		/* don't use implicit vals (yet) */
   parm[2].chancalc =  calcna5i;        /* default rate function */
   parm[2].funcname = (char *)"calcna5i";       /* user rate function */
   parm[2].voff = dnaoffsm;             /* voltage offset from user */

   parm[2].dq[0] = dqm;                 /* Q10 for a4 */
   parm[2].dq[1] = dqm;                 /* Q10 for a5 */
   parm[2].dq[2] = dqm;                 /* Q10 for b4 */
   parm[2].dq[3] = dqm;                 /* Q10 for b5 */
   parm[2].dq[4] = 1.0;
   parm[2].dq[5] = 1.0;

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = na5a11;
   spnt[0].ratemul[0] = 1.0 * m;
   spnt[0].rateo  [0] = 0;  /* taum */
   spnt[0].trans  [1] = 4;
   spnt[0].trate  [1] = na5b3;
   spnt[0].ratemul[1] = 1.0 * h;
   spnt[0].rateo  [1] = 2; 	/* tauh */

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = na5a12;
   spnt[1].ratemul[0] = 1.0 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = na5b11;
   spnt[1].ratemul[1] = 1.0 * m;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 5;	
   spnt[1].trate  [2] = na5b3;
   spnt[1].ratemul[2] = 1.0 * h;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = na5a13;
   spnt[2].ratemul[0] = 1.0 * m;
   spnt[2].rateo  [0] = 0;	/* taum */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = na5b12;
   spnt[2].ratemul[1] = 1.0 * m;
   spnt[2].rateo  [1] = 1;	/* taum */
   spnt[2].trans  [2] = 6;
   spnt[2].trate  [2] = na5b3;
   spnt[2].ratemul[2] = 1.0 * h;
   spnt[2].rateo  [2] = 2;	/* tauh, tauc */

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;		/* state 3 = the open state */
   spnt[3].trans  [0] = 2;
   spnt[3].trate  [0] = na5b13;
   spnt[3].ratemul[0] = 1.0 * m;
   spnt[3].rateo  [0] = 1;	/* taum, taub */
   spnt[3].trans  [1] = 6;
   spnt[3].trate  [1] = na5a2;
   spnt[3].ratemul[1] = 1.0 * m;
   spnt[3].rateo  [1] = 2;	/* tauh */

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = na5a11;
   spnt[4].ratemul[0] = 1.0 * m;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 0;
   spnt[4].trate  [1] = na5a3;
   spnt[4].ratemul[1] = 1.0 * m;
   spnt[4].rateo  [1] = 3;

   spnt[5].numtrans   = 3;
   spnt[5].cond       = 0;		
   spnt[5].trans  [0] = 6;
   spnt[5].trate  [0] = na5a12;
   spnt[5].ratemul[0] = 1.0 * m;
   spnt[5].rateo  [0] = 0;	/* taum, taua */
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = na5b11;
   spnt[5].ratemul[1] = 1.0 * m;
   spnt[5].rateo  [1] = 1;	/* taum, taub */
   spnt[5].trans  [2] = 1;
   spnt[5].trate  [2] = na5a3;
   spnt[5].ratemul[2] = 1.0 * h;
   spnt[5].rateo  [2] = 3;	/* tauh, tauc */

   spnt[6].numtrans   = 4;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = na5a4;
   spnt[6].ratemul[0] = 1.0 * m;	
   spnt[6].rateo  [0] = 0;
   spnt[6].trans  [1] = 5;
   spnt[6].trate  [1] = na5b12;
   spnt[6].ratemul[1] = 1.0 * h;
   spnt[6].rateo  [1] = 2;	/* tauh, tauc */
   spnt[6].trans  [2] = 2;
   spnt[6].trate  [2] = na5a3;
   spnt[6].ratemul[2] = 1.0 * h;
   spnt[6].rateo  [2] = 3;	/* tauh, taud */
   spnt[6].trans  [3] = 3;
   spnt[6].trate  [3] = na5b2;
   spnt[6].ratemul[3] = 1.0 * h;
   spnt[6].rateo  [3] = 3;	/* tauh, taud */

   spnt[7].numtrans   = 2;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 8;
   spnt[7].trate  [0] = na5a5;
   spnt[7].ratemul[0] = 1.0 * m;
   spnt[7].rateo  [0] = 0;
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = na5b4;
   spnt[7].ratemul[1] = 1.0 * m;
   spnt[7].rateo  [1] = 1;

   spnt[8].numtrans   = 1;
   spnt[8].cond       = 0;
   spnt[8].trans  [0] = 7;
   spnt[8].trate  [0] = na5b5;
   spnt[8].ratemul[0] = 1.0 * m;
   spnt[8].rateo  [0] = 1;

   return ch;
}
/*----------------------------------------*/


