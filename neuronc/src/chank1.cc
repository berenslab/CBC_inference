/* segment chank1 in program nc */

/* sets up channel parameters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncomp.h"
#include "control.h"

chantype *makchantype(int ctype, int cnum, int nstates, int nparm, 
					int nq, double bt);
chantype *getchantype(int ctype, int cnum);
void mak2state(chanstate *spnt, double (*frate)(chan *cpnt),
                                double (*rrate)(chan *cpnt), double rate);
double rchanf(chan *cpnt);
double rchanr(chan *cpnt);

double al0 (chan *cpnt);
double bet0(chan *cpnt);

static double (*aln) (chan*) = al0;
static double (*betn)(chan*) = bet0;

#ifdef __cplusplus
extern "C" {
#endif

double exp(double);

#ifdef __cplusplus
}
#endif

double ncabs(double x);

/*----------------------------------------*/

/* Vo is used here to define the resting potential
   that Hodgkin and Huxley (1952) used as a base 
   for their rate constant functions.  Vo is only
   used to normalize the modern definition of 
   membrane potential to their definition.
   "Vo" should not be modified as it refers only to
   the antiquated definition of membrane voltage: */

#define Vo (-65.)

/* The transformation from HH voltage to the standard
   used here is:

   new V =   - (oldV-Vo)
or
   new V =   (Vo - oldV)

*/

/*----------------------------------------*/

double calck1n (double v, int func)

/* Calculate K rate functions given voltage in mv.
   All rates are calculated exactly as in HH (1952) paper.
   The "func" parameter defines:

    1	alpha n 
    2	beta  n
*/

{
   double val,x,y;

  switch (func) {
						/* alpha functions */

  case 1:					/* alpha n */
/*    y = 0.1 * (v+10.);                        /* the old way */
    y = -0.1 * (v+55.);                         /* the modern way */
    x = exp (y) - 1.;
    if (ncabs(x) > 1e-5)                        /* singularity when x==0 */
       val = MSSEC * 0.1 * y / x;
    else
       val = MSSEC * 0.1;
    val *= dratehhna;				/* normalize to dbasetc */ 
    break;

  case 2:					/* beta n */
/*    val = 0.125 * exp (v / 80.);              /* the old way */

    val = MSSEC * 0.125 * exp ((v-Vo) / -80.);  /* the modern way */
    val *= dratehhnb;				/* normalize to dbasetc */ 
    break;

  }
  return val;
}						/*  i.e. 22 deg C */

/*--------------------------------------------*/

chantype *makk0(void)
{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 2;                          /* 2 Markov states */
   nparm = 1;                           /* make 1 set of params, "n" */
   nq = 2;				/* make 2 Q10 values */
   ch=makchantype(K,0,nstate,nparm,nq,dbasetc); /* make chan state const, parm*/
   ch->hh = 1;				/* HH channel = uses "n" = "m" */
   ch->unitary = dku;			/* 6 pS @ 6.3 -> 33 pS @ 22 deg */
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions /perm K */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 2;                /* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calck1n;         /* default rate function */
   parm[0].funcname = (char *)"calck0n"; /* user rate function */
   parm[0].dq[0] = dqna;                /* Q10 for n rate function */
   parm[0].dq[1] = dqnb;                /* Q10 for n rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

   mak2state(spnt,rchanf,rchanr,CHANRATE / dktauf * 1.5);

   return ch;
}

/*--------------------------------------------*/

chantype *makk1(void)
{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 5;                          /* 5 Markov states */
   nparm = 1;                           /* make 1 set of params, "n" */
   nq = 2;				/* make 1 Q10 value */
   ch=makchantype(K,1,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dku;
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions /perm K */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calck1n;         /* default rate function */
   parm[0].funcname = (char *)"calck1n"; /* user rate function */
   parm[0].dq[0] = dqna;                /* Q10 for n rate function */
   parm[0].dq[1] = dqnb;                /* Q10 for n rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

   n = 1.0;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = aln;
   spnt[0].ratemul[0] = 4.0*n;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;              /*   0 <-> 1 <-> 2 <-> 3 <-> 4 */
   spnt[1].trate  [0] = aln;
   spnt[1].ratemul[0] = 3.0*n;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = betn;
   spnt[1].ratemul[1] = 1.0*n;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = aln;
   spnt[2].ratemul[0] = 2.0*n;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = betn;
   spnt[2].ratemul[1] = 2.0*n;
   spnt[2].rateo  [1] = 1;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = aln;
   spnt[3].ratemul[0] = 1.0*n;
   spnt[3].rateo  [0] = 0;
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = betn;
   spnt[3].ratemul[1] = 3.0*n;
   spnt[3].rateo  [1] = 1;

   spnt[4].numtrans   = 1;                      /* state 4 = the open state */
   spnt[4].cond       = 1.0;
   spnt[4].trans  [0] = 3;
   spnt[4].trate  [0] = betn;
   spnt[4].ratemul[0] = 4.0*n;
   spnt[4].rateo  [0] = 1;
 
   return ch;
}

/*----------------------------------------*/

chantype *makk6(void)
{
 /* Kv3 channel, duplicate version of K1 above, used for slower version */

     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 5;                          /* 5 Markov states */
   nparm = 1;                           /* make 1 set of params, "n" */
   nq = 2;				/* make 1 Q10 value */
   ch=makchantype(K,6,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dku;
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions /perm K */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calck1n;         /* default rate function */
   parm[0].funcname = (char *)"calck1n"; /* user rate function */
   parm[0].dq[0] = dqna;                /* Q10 for n rate function */
   parm[0].dq[1] = dqnb;                /* Q10 for n rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

   n = 1.0;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = aln;
   spnt[0].ratemul[0] = 4.0*n;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;              /*   0 <-> 1 <-> 2 <-> 3 <-> 4 */
   spnt[1].trate  [0] = aln;
   spnt[1].ratemul[0] = 3.0*n;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = betn;
   spnt[1].ratemul[1] = 1.0*n;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = aln;
   spnt[2].ratemul[0] = 2.0*n;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = betn;
   spnt[2].ratemul[1] = 2.0*n;
   spnt[2].rateo  [1] = 1;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = aln;
   spnt[3].ratemul[0] = 1.0*n;
   spnt[3].rateo  [0] = 0;
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = betn;
   spnt[3].ratemul[1] = 3.0*n;
   spnt[3].rateo  [1] = 1;

   spnt[4].numtrans   = 1;                      /* state 4 = the open state */
   spnt[4].cond       = 1.0;
   spnt[4].trans  [0] = 3;
   spnt[4].trate  [0] = betn;
   spnt[4].ratemul[0] = 4.0*n;
   spnt[4].rateo  [0] = 1;
 
   return ch;
}

/*----------------------------------------*/

chantype *makk7(void)
{
 /* Kv3 channel, duplicate version of K1 above, used for slower version */

     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 5;                          /* 5 Markov states */
   nparm = 1;                           /* make 1 set of params, "n" */
   nq = 2;				/* make 1 Q10 value */
   ch=makchantype(K,7,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = dku;
   ch->ions->ionp[PK]  = 1.0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = dpnak;		/* permeability to Na+ ions / perm K */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions /perm K */
   ch->vrev = vk;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 2;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* number of implicit vals (m1, m2) */
   parm[0].chancalc =  calck1n;         /* default rate function */
   parm[0].funcname = (char *)"calck1n"; /* user rate function */
   parm[0].dq[0] = dqna;                /* Q10 for n rate function */
   parm[0].dq[1] = dqnb;                /* Q10 for n rate function */
   parm[0].voff = dkoffsn;              /* voltage offset from user */

   n = 1.0;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = aln;
   spnt[0].ratemul[0] = 4.0*n;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;              /*   0 <-> 1 <-> 2 <-> 3 <-> 4 */
   spnt[1].trate  [0] = aln;
   spnt[1].ratemul[0] = 3.0*n;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = betn;
   spnt[1].ratemul[1] = 1.0*n;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = aln;
   spnt[2].ratemul[0] = 2.0*n;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = betn;
   spnt[2].ratemul[1] = 2.0*n;
   spnt[2].rateo  [1] = 1;

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;
   spnt[3].trate  [0] = aln;
   spnt[3].ratemul[0] = 1.0*n;
   spnt[3].rateo  [0] = 0;
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = betn;
   spnt[3].ratemul[1] = 3.0*n;
   spnt[3].rateo  [1] = 1;

   spnt[4].numtrans   = 1;                      /* state 4 = the open state */
   spnt[4].cond       = 1.0;
   spnt[4].trans  [0] = 3;
   spnt[4].trate  [0] = betn;
   spnt[4].ratemul[0] = 4.0*n;
   spnt[4].rateo  [0] = 1;
 
   return ch;
}

/*----------------------------------------*/


