/* segment channa2 in program nc */

/* sets up channel na 2 parameters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
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

/*--------------------------------------------*/

/* For the purpose of normalizing kinetics to 22 deg C, */
/*  we assume Q10=2 for m and Q10=3 for h. */

/* #define HHM5RATE  (exp(log(2.0) * (BASETC- 5.0)/10.0); */
#define    HHM5RATE  3.2490096

/* #define HHH5RATE  (exp(log(3.0) * (BASETC- 5.0)/10.0); */
#define    HHH5RATE  6.4730078

/* - - - - - - - - - - - - - - - - - - - - - -*/

double calcna2m (double v, int func)

/* Calculate na type 2 rate functions x, y, a, b, c, d.
   Rate functions are taken from Vandenberg and Bezanilla, 1991,
   Biophys J. 60:1511-1533. */

/* Calculate na2 rate funcs, from Vandenberg and Bezanilla, 1991. */
/* Calibrated in mV, per sec */
/* Temperature 5 deg C */

/* 
   Changed the divisor in exponent for y and z to 18 from the
   original 24.  This shortens toe of activation curve.
   RGS, Sept 27, 2001
*/

/*   The original:                                                     */

/*     case 1: val = 16609. * exp ( 1.50 * 0.22/24. * v);  break; /* y */
/*     case 2: val = 971.   * exp (-1.50 * 0.78/24. * v);  break; /* z */

{
   double val;

   if (func > 6) func = 6;
   switch (func) {

     case 1: val = 16609. * exp ( 1.50 * 0.22/18. * v);  break; /* y */
     case 2: val = 971.   * exp (-1.50 * 0.78/18. * v);  break; /* z */
     case 3: val = 5750.  * exp ( 0.42 * 0.99/24. * v);  break; /* a */
     case 4: val = 4325.  * exp (-0.42 * 0.01/24. * v);  break; /* b */
     case 5: val = 15669. * exp ( 1.91 * 0.75/24. * v);  break; /* c */
     case 6: val = 1361.  * exp (-1.91 * 0.25/24. * v);  break; /* d */
      
     }
   return val * HHM5RATE;
}

/*--------------------------------------------*/

double calcna2h (double v, int func)

/* Calculate na2 rate inactivation funcs g, f, i, j, 
   from Vandenberg and Bezanilla, 1991,
   Biophys J. 60:1511-1533. */

/* Calibrated in mV, per sec */
/* Temperature 5 deg C */
/* This function is separate from calcna2m because the Q10 for
   inactivation functions is assumed to possibly differ from
   that for activation functions.  Therefore, we place these rate
   function values into a different "chanparm" which has a separate
   Q10 value.
*/

{
   double val;

   if (func > 4) func = 4;
   switch (func) {

     case 1: val = 4.     * exp (-0.91 * 0.999/24. * v); break; /* i */
     case 2: val = 432.   * exp ( 0.91 * 0.001/24. * v); break; /* f */
     case 3: val = 770. / 432.  *                             /* g*i/f */
                   (4.    * exp (-0.91 * 0.999/24. * v)); break; /* j */
     case 4: val = 770.   * exp ( 0.91 * 0.001/24. * v); break; /* g */
      
     }
   return val * HHH5RATE;
}

/*--------------------------------------------*/

double nay(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double naz(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double naa(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double nab(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[3];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double nac(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[4];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double nad(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[5];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double naf(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double nag(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[3];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double nai(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/*- - - - - - - - - - - - - - - - - - - - - - */

double naj(chan *cpnt)

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[2];
}

/*--------------------------------------------*/

chantype *makna2(void)

/* Discrete-State-Markov description of 
   Na+ channel that very closely mimics 
   kinetics of patch clamp data from squid giant axon.

   See Vandenberg and Bezanilla (1991) Biophysical Journal 60:1511-1533.

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

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 9;				/* 9 Markov states */
   nparm = 2;				/* make 2 sets of seq. state states */
   nq = 1;				/* make 1 Q10 value per parm */
   ch = makchantype(NA,2, nstate, nparm, nq, dbasetc); /* make chan info */
   ch->unitary = dnau;
   ch->ions->ionp[PNA] = 1.0;			/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = dpkna;		/* permeability to K+ ions / Na perm */
   ch->ions->ionp[PCA] = dpcana;		/* permeability to Ca++ ions / Na perm */
   ch->vrev = vna;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 6;		/* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;		/* don't use implicit vals (yet) */
   parm[0].chancalc =  calcna2m;        /* default rate function */
   parm[0].funcname = (char *)"calcna2m"; /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */
   parm[0].voff = dnaoffsm;             /* voltage offset from user */

   parm[1].nfval    = 4;		/* num of func vals (ah, bh, etc.) */
   parm[1].nival    = 0;		/* don't use implicit vals (yet) */
   parm[1].chancalc =  calcna2h;        /* default rate function */
   parm[1].funcname = (char *)"calcna2h"; /* user rate function */
   parm[1].dq[0] = dqh;                 /* Q10 for h rate function */
   parm[1].voff = dnaoffsh;             /* voltage offset from user */

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = nay;
   spnt[0].ratemul[0] = 1.0 * m;
   spnt[0].rateo  [0] = 0;  /* taum */

			     /*                               open */

                             /*   0 <-> 1 <-> 2 <-> 3 <-> 4 <-> 5  */

                             /*                     |           |  */

                             /*                     |           |  */

                             /*                     6 <-> 7 <-> 8  */
   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = nay;
   spnt[1].ratemul[0] = 1.0 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = naz;
   spnt[1].ratemul[1] = 1.0 * m;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = nay;
   spnt[2].ratemul[0] = 1.0 * m;
   spnt[2].rateo  [0] = 0;	/* taum */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = naz;
   spnt[2].ratemul[1] = 1.0 * m;
   spnt[2].rateo  [1] = 1;	/* taum */

   spnt[3].numtrans   = 3;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = naa;
   spnt[3].ratemul[0] = 1.0 * m;
   spnt[3].rateo  [0] = 0;	/* taum */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = naz;
   spnt[3].ratemul[1] = 1.0 * m;
   spnt[3].rateo  [1] = 1;	/* taum */
   spnt[3].trans  [2] = 6;
   spnt[3].trate  [2] = nag;
   spnt[3].ratemul[2] = 1.0 * h;
   spnt[3].rateo  [2] = 3;	/* tauh, taud */

			     /*      y>    y>    y>    a>    c> open */

                             /*   0 <-> 1 <-> 2 <-> 3 <-> 4 <-> 5   */

                             /*     <z    <z    <z  | <b    <d  |   */
			     
			     /*                       ^           ^ */
                             /*                   g | j       f | i */

			     /*                        a>    c>     */

                             /*                     6 <-> 7 <-> 8   */

			     /*                       <b    <d      */
   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = nac;
   spnt[4].ratemul[0] = 1.0 * m;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 3;
   spnt[4].trate  [1] = nab;
   spnt[4].ratemul[1] = 1.0 * m;
   spnt[4].rateo  [1] = 1;

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 1.0;		/* state 5 = the open state */
   spnt[5].trans  [0] = 8;
   spnt[5].trate  [0] = naf;
   spnt[5].ratemul[0] = 1.0 * h;
   spnt[5].rateo  [0] = 3;		/* tauh, taud */
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = nad;
   spnt[5].ratemul[1] = 1.0 * m;
   spnt[5].rateo  [1] = 1;		/* taum, taub */

   spnt[6].numtrans   = 2;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = naa;
   spnt[6].ratemul[0] = 1.0 * m;	
   spnt[6].rateo  [0] = 0;
   spnt[6].trans  [1] = 3;
   spnt[6].trate  [1] = naj;
   spnt[6].ratemul[1] = 1.0 * h;
   spnt[6].rateo  [1] = 2;		/* tauh, tauc */

   spnt[7].numtrans   = 2;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 8;
   spnt[7].trate  [0] = nac;
   spnt[7].ratemul[0] = 1.0 * m;
   spnt[7].rateo  [0] = 0;		/* taua, taum */
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = nab;
   spnt[7].ratemul[1] = 1.0 * m;
   spnt[7].rateo  [1] = 1;		/* taub, taum */

   spnt[8].numtrans   = 2;
   spnt[8].cond       = 0;		/* state 8 = the inactivated state */
   spnt[8].trans  [0] = 5;
   spnt[8].trate  [0] = nai;
   spnt[8].ratemul[0] = 1.0 * h;
   spnt[8].rateo  [0] = 2;		/* tauc */
   spnt[8].trans  [1] = 7;
   spnt[8].trate  [1] = nad;
   spnt[8].ratemul[1] = 1.0 * m;
   spnt[8].rateo  [1] = 1;

   return ch;
}

/*--------------------------------------------*/

chantype *makna8(void)

/* Discrete-State-Markov description of 
   Na+ channel that mimics behavior of slow TTX-resistant Na
   channel (Nav1.8).

   Same as Na type 2 above except different multipliers:

   taua = 10;    slows activation
   taub = 10;    slows activation
   tauc = 3;     makes plateau about 50% of peak
   taud = 80;    makes longer inactivation time constant

   See Vandenberg and Bezanilla (1991) Biophysical Journal 60:1511-1533.

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

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double m,h;
     int nstate, nparm, nq;

   nstate = 9;				/* 9 Markov states */
   nparm = 2;				/* make 2 sets of seq. state states */
   nq = 1;				/* make 1 Q10 value per parm */
   ch = makchantype(NA,8, nstate, nparm, nq, dbasetc); /* make chan info */
   ch->unitary = dnau;
   ch->ions->ionp[PNA] = 1.0;			/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = dpkna;		/* permeability to K+ ions / Na perm */
   ch->ions->ionp[PCA] = dpcana;		/* permeability to Ca++ ions / Na perm */
   ch->vrev = vna;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 6;		/* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;		/* don't use implicit vals (yet) */
   parm[0].chancalc =  calcna2m;        /* default rate function */
   parm[0].funcname = (char *)"calcna2m"; /* user rate function */
   parm[0].dq[0] = dqm;                 /* Q10 for m rate function */
   parm[0].voff = dnaoffsm;             /* voltage offset from user */

   parm[1].nfval    = 4;		/* num of func vals (ah, bh, etc.) */
   parm[1].nival    = 0;		/* don't use implicit vals (yet) */
   parm[1].chancalc =  calcna2h;        /* default rate function */
   parm[1].funcname = (char *)"calcna2h"; /* user rate function */
   parm[1].dq[0] = dqh;                 /* Q10 for h rate function */
   parm[1].voff = dnaoffsh;             /* voltage offset from user */

   m = 1.0;
   h = 1.0;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = nay;
   spnt[0].ratemul[0] = 0.1 * m;
   spnt[0].rateo  [0] = 0;  /* taum */

			     /*                               open */

                             /*   0 <-> 1 <-> 2 <-> 3 <-> 4 <-> 5  */

                             /*                     |           |  */

                             /*                     |           |  */

                             /*                     6 <-> 7 <-> 8  */
   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = nay;
   spnt[1].ratemul[0] = 0.1 * m;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;	
   spnt[1].trate  [1] = naz;
   spnt[1].ratemul[1] = 0.1 * m;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = nay;
   spnt[2].ratemul[0] = 0.1 * m;
   spnt[2].rateo  [0] = 0;	/* taum */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = naz;
   spnt[2].ratemul[1] = 0.1 * m;
   spnt[2].rateo  [1] = 1;	/* taum */

   spnt[3].numtrans   = 3;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = naa;
   spnt[3].ratemul[0] = 0.1 * m;
   spnt[3].rateo  [0] = 0;	/* taum */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = naz;
   spnt[3].ratemul[1] = 0.1 * m;
   spnt[3].rateo  [1] = 1;	/* taum */
   spnt[3].trans  [2] = 6;
   spnt[3].trate  [2] = nag;
   spnt[3].ratemul[2] = 0.0125 * h;
   spnt[3].rateo  [2] = 3;	/* tauh, taud */

			     /*      y>    y>    y>    a>    c> open */

                             /*   0 <-> 1 <-> 2 <-> 3 <-> 4 <-> 5   */

                             /*     <z    <z    <z  | <b    <d  |   */
			     
			     /*                       ^           ^ */
                             /*                   g | j       f | i */

			     /*                        a>    c>     */

                             /*                     6 <-> 7 <-> 8   */

			     /*                       <b    <d      */
   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = nac;
   spnt[4].ratemul[0] = 0.1 * m;
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 3;
   spnt[4].trate  [1] = nab;
   spnt[4].ratemul[1] = 0.1 * m;
   spnt[4].rateo  [1] = 1;

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 1.0;		/* state 5 = the open state */
   spnt[5].trans  [0] = 8;
   spnt[5].trate  [0] = naf;
   spnt[5].ratemul[0] = 0.0125 * h;
   spnt[5].rateo  [0] = 3;		/* tauh, taud */
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = nad;
   spnt[5].ratemul[1] = 0.1 * m;
   spnt[5].rateo  [1] = 1;		/* taum, taub */

   spnt[6].numtrans   = 2;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 7;
   spnt[6].trate  [0] = naa;
   spnt[6].ratemul[0] = 0.1 * m;	
   spnt[6].rateo  [0] = 0;
   spnt[6].trans  [1] = 3;
   spnt[6].trate  [1] = naj;
   spnt[6].ratemul[1] = 0.333 * h;
   spnt[6].rateo  [1] = 2;		/* tauh, tauc */

   spnt[7].numtrans   = 2;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 8;
   spnt[7].trate  [0] = nac;
   spnt[7].ratemul[0] = 0.1 * m;
   spnt[7].rateo  [0] = 0;		/* taua, taum */
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = nab;
   spnt[7].ratemul[1] = 0.1 * m;
   spnt[7].rateo  [1] = 1;		/* taub, taum */

   spnt[8].numtrans   = 2;
   spnt[8].cond       = 0;		/* state 8 = the inactivated state */
   spnt[8].trans  [0] = 5;
   spnt[8].trate  [0] = nai;
   spnt[8].ratemul[0] = 0.333 * h;
   spnt[8].rateo  [0] = 2;		/* tauc */
   spnt[8].trans  [1] = 7;
   spnt[8].trate  [1] = nad;
   spnt[8].ratemul[1] = 0.1 * m;
   spnt[8].rateo  [1] = 1;

   return ch;
}

