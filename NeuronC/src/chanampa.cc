/* segment chanampa in program nc */

/* sets up channel parameters */

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

double exp(double);
double log(double);

#ifdef __cplusplus
}
#endif

chantype *makchantype(int ctype, int cnum, int nstates, int nparm, 
				int nq, double bt);

double rnt(chan *spnt);
double rt(chan *spnt);
int set_perms (chantype *chtyp, iontab *ions, double svrev);

/*--------------------------------------------*/

double calcampa1(double v, int func) 

{
return 1.0;
}

/*--------------------------------------------*/

chantype *makampa1(void)

/* number states from 0 to n; 
   set numstate = n;
   set numtrans = number of transitions.
   set cond     = conductance of state;

   for each state, set transitions:
    set trans = state to go to on this transition. 
    set trate = function that returns basic rate for transition. 
    set ratemul = multiplier for rate function.
    set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.
*/

/* Taken from:
                Jonas, P., Major, G., and Sakmann, B. (1993) 
                J. Physiol. 472:615-5663.

   in which a 7-state AMPA markov diagram is provided.  Two sets
of rate constants were derived, 1) which gave a good fit, and 2)
which was derived assuming independence of binding between
binding of glutamate molecules.  The one defined here is the "good
fit", not assuming independence.

     k1f       k2f        a
  0  <->   1   <->   2   <->   3
     k1r       k2r        b
                                    
             ^        ^         ^
        a1 | b1   a2 | b2   a3 | b3
        \/        \/        \/

               k3f        a4
           4   <->   5   <->   6
               k3r        b4         

State 0 is unbound
      1 is partially bound
      2 is fully bound
      3 is bound and open
      4,5,6 are deactivated
 
*/

  {
     double rnt(chan *cp),rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k, pca;
     int nstate, nparm, nq;

#define ampa1f1 4.59e6
#define ampa1r1 4.26e3
#define ampa1f2 28.4e6
#define ampa1r2 3.26e3
#define ampa1f3 1.27e6
#define ampa1r3 45.7
#define ampa1a  4.24e3
#define ampa1b  900
#define ampa1a1 2.89e3
#define ampa1b1 39.2
#define ampa1a2 172
#define ampa1b2 0.727
#define ampa1a3 17.7
#define ampa1b3 4.0
#define ampa1a4 16.8 
#define ampa1b4 190.4 

   nstate = 7;                          /* 7 Markov states */
   nparm = 1;                           /* make 1 set of params */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(AMPA,1,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   ch->unitary = dampau;
   ch->trconc = dstr;
   pca = dpcaampa;
   //ch->ions->ionp[PCA] = pca;		/* permeability to Ca++ ions */
   ch->ions->ionp[PCA] = 0;		/* permeability to Ca++ ions */
   ch->ions->ionp[PNA] = 1.0;		/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0;		/* permeability to K+ ions */

   ch->vrev = vglu;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcampa1;       /* default rate function */
   parm[0].funcname = (char *)"calcampa1";      /* user rate function */
   parm[0].dq[0] = dqsyn;               /* Q10 for m rate function */

   respamp[GLU-GLU] = 1.0;		/* response to glutamate */
   respamp[AMPA-GLU] = 0.5;		/* response to AMPA */
   respamp[NMDA-GLU] = 1e-3;		/* response to NMDA */
   respamp[CNQX-GLU] = -1.0;		/* response to CNQX */

   k=1;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = ampa1f1 * k;
   spnt[0].rateo  [0] = 0;	/* offsm */

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0; 
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = ampa1f2 * k;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = ampa1r1 * k;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 4;
   spnt[1].trate  [2] = rt;
   spnt[1].ratemul[2] = ampa1a1 * k;
   spnt[1].rateo  [2] = 2;      /* offsh */

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rt;		/* function of time */
   spnt[2].ratemul[0] = ampa1a * k;
   spnt[2].rateo  [0] = 0;	/* offsm */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;		/* function of nt conc, time */
   spnt[2].ratemul[1] = ampa1r2 * k;
   spnt[2].rateo  [1] = 1;	/* offsm */
   spnt[2].trans  [2] = 5;
   spnt[2].trate  [2] = rt;
   spnt[2].ratemul[2] = ampa1a2 * k;
   spnt[2].rateo  [2] = 2;	/* offsh */

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 6;			/* state 3 = the open state */
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = ampa1a3 * k;
   spnt[3].rateo  [0] = 2;	/* offsh */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = ampa1b * k;
   spnt[3].rateo  [1] = 1;	/* offsm */

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rnt;
   spnt[4].ratemul[0] = ampa1f3;	
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 1;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = ampa1b1 * k;
   spnt[4].rateo  [1] = 3;

   spnt[5].numtrans   = 3;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 6;			/* state 5 = inactivated state*/
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = ampa1a4 * k;
   spnt[5].rateo  [0] = 0;
   spnt[5].trans  [1] = 2;
   spnt[5].trate  [1] = rt;
   spnt[5].ratemul[1] = ampa1b2 * k;
   spnt[5].rateo  [1] = 3;
   spnt[5].trans  [2] = 4;
   spnt[5].trate  [2] = rt;
   spnt[5].ratemul[2] = ampa1r3;
   spnt[5].rateo  [2] = 1;

   spnt[6].numtrans   = 2;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 3;			/* state 6 = inactivated state*/
   spnt[6].trate  [0] = rt;
   spnt[6].ratemul[0] = ampa1b3 * k;
   spnt[6].rateo  [0] = 3;	/* offsh */
   spnt[6].trans  [1] = 5;
   spnt[6].trate  [1] = rt;
   spnt[6].ratemul[1] = ampa1b4 * k;
   spnt[6].rateo  [1] = 1;	/* offsm */

   return ch;

#undef ampa1f1 
#undef ampa1r1
#undef ampa1f2
#undef ampa1r2
#undef ampa1f3
#undef ampa1r3
#undef ampa1a 
#undef ampa1b 
#undef ampa1a1
#undef ampa1b1
#undef ampa1a2
#undef ampa1b2
#undef ampa1a3
#undef ampa1b3
#undef ampa1a4
#undef ampa1b4
}

/*--------------------------------------------*/

chantype *makampa2(void)

/* number states from 0 to n; 
   set numstate = n;
   set numtrans = number of transitions.
   set cond     = conductance of state;

   for each state, set transitions:
    set trans = state to go to on this transition. 
    set trate = function that returns basic rate for transition. 
    set ratemul = multiplier for rate function.
    set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.
*/

/* Taken from:
                Jonas, P., Major, G., and Sakmann, B. (1993) 
                J. Physiol. 472:615-5663.

   in which a 7-state AMPA markov diagram is provided.  Two sets
of rate constants were derived, 1) which gave a good fit, and 2)
which was derived assuming independence of binding between
binding of glutamate molecules.  The one defined here is the 
"assuming independence".  States are same as ampa1 above but with
different rate functions.

State 0 is unbound
      1 is partially bound
      2 is fully bound
      3 is bound and open
      4,5,6 are deactivated
 
*/

  {
     double rnt(chan *cp),rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double a,b,c,d;
     double pca;
     int nstate, nparm, nq;

#define ampa1f1 26.6e6
#define ampa1r1 6.24e3
#define ampa1f2 13.3e6
#define ampa1r2 12.5e3
#define ampa1f3 2.41e6
#define ampa1r3 57.1 
#define ampa1a  4.20e3
#define ampa1b  302 
#define ampa1a1 513 
#define ampa1b1 28.1 
#define ampa1a2 395 
#define ampa1b2 0.546
#define ampa1a3 109 
#define ampa1b3 33.4 
#define ampa1a4 8.15 
#define ampa1b4 130 

   nstate = 7;                          /* 7 Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(AMPA,2,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   pca = dpcaampa;
   ch->ions->ionp[PCA] = pca;		/* permeability to Ca++ ions */
   ch->ions->ionp[PNA] = 1.0;		/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0;		/* permeability to K+ ions */

   ch->vrev = vglu;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   ch->unitary = dampau;
   ch->trconc = dstr;
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcampa1;       /* default rate function */
   parm[0].funcname = (char *)"calcampa2";      /* user rate function */
   parm[0].dq[0] = dqsyn;               /* Q10 for m rate function */

   respamp[GLU-GLU] = 1.0;		/* response to glutamate */
   respamp[AMPA-GLU] = 0.5;		/* response to AMPA */
   respamp[NMDA-GLU] = 1e-3;		/* response to NMDA */
   respamp[CNQX-GLU] = -1.0;		/* response to CNQX */

   a = 1.0;	/* forward: alpha */
   b = 1.0;     /* reverse: beta */
   c = 1.0;     /* deactivation */
   d = 1.0;     /* re-activation */
 
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = ampa1f1 * a;
   spnt[0].rateo  [0] = 0;	/* offsm */

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = ampa1f2 * a;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = ampa1r1 * b;
   spnt[1].rateo  [1] = 0;
   spnt[1].trans  [2] = 4;
   spnt[1].trate  [2] = rt;
   spnt[1].ratemul[2] = ampa1a1 * c;
   spnt[1].rateo  [2] = 1;      /* offsh */

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rt;		/* function of time */
   spnt[2].ratemul[0] = ampa1a * a;
   spnt[2].rateo  [0] = 0;	/* offsm */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;		/* function of nt conc, time */
   spnt[2].ratemul[1] = ampa1r2 * b;
   spnt[2].rateo  [1] = 0;	/* offsm */
   spnt[2].trans  [2] = 5;
   spnt[2].trate  [2] = rt;
   spnt[2].ratemul[2] = ampa1a2 * c;
   spnt[2].rateo  [2] = 1;	/* offsh */

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 6;			/* state 3 = the open state */
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = ampa1a3 * c;
   spnt[3].rateo  [0] = 1;	/* offsh */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = ampa1b * b;
   spnt[3].rateo  [1] = 0;	/* offsm */

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rnt;
   spnt[4].ratemul[0] = ampa1f3 * a;	
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 1;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = ampa1b1 * d;
   spnt[4].rateo  [1] = 1;

   spnt[5].numtrans   = 3;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 6;			/* state 5 = inactivated state*/
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = ampa1a4 * a;
   spnt[5].rateo  [0] = 1;	/* offsh */
   spnt[5].trans  [1] = 2;
   spnt[5].trate  [1] = rt;
   spnt[5].ratemul[1] = ampa1b2 * d;
   spnt[5].rateo  [1] = 1;	/* offsh */
   spnt[5].trans  [2] = 4;
   spnt[5].trate  [2] = rt;
   spnt[5].ratemul[2] = ampa1r3 * b;
   spnt[5].rateo  [2] = 0;	/* offsm */

   spnt[6].numtrans   = 2;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 3;			/* state 6 = inactivated state*/
   spnt[6].trate  [0] = rt;
   spnt[6].ratemul[0] = ampa1b3 * d;
   spnt[6].rateo  [0] = 1;	/* offsh */
   spnt[6].trans  [1] = 5;
   spnt[6].trate  [1] = rt;
   spnt[6].ratemul[1] = ampa1b4 * b;
   spnt[6].rateo  [1] = 0;	/* offsm */

   return ch;

#undef ampa1f1 
#undef ampa1r1
#undef ampa1f2
#undef ampa1r2
#undef ampa1f3
#undef ampa1r3
#undef ampa1a 
#undef ampa1b 
#undef ampa1a1
#undef ampa1b1
#undef ampa1a2
#undef ampa1b2
#undef ampa1a3
#undef ampa1b3
#undef ampa1a4
#undef ampa1b4
}

/*--------------------------------------------*/

chantype *makampa3(void)

/* number states from 0 to n; 
   set numstate = n;
   set numtrans = number of transitions.
   set cond     = conductance of state;

   for each state, set transitions:
    set trans = state to go to on this transition. 
    set trate = function that returns basic rate for transition. 
    set ratemul = multiplier for rate function.
    set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.
*/
/* Taken from:
                Raman, I.M, and Trussel, L.O. (1995) 
                J. Physiol. 472:615-5663.

   in which a 9-state AMPA markov diagram is provided:

                          4
                          |
                          |
              0 <-> 1 <-> 2 <-> 3
                    |     |
                    |     |
                    5 <-> 6 <-> 7 <-> 8

     State 0 is unbound
           1, 2 are bound
           3, 4, and 8 are bound and open
           5, 6 are bound and desensitized 
           7 is bound but closed
           8 is open
*/

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double a,b,c,d;
     double pca;
     int nstate, nparm, nq;

   nstate = 9;                          /* 9 Markov states */
   nparm = 1;                           /* make 1 set of params */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(AMPA,3,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   pca = dpcaampa;
   ch->ions->ionp[PCA] = pca;		/* permeability to Ca++ ions */
   ch->ions->ionp[PNA] = 1.0;		/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0;		/* permeability to K+ ions */

   ch->vrev = vglu;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   ch->unitary = dampau;
   ch->trconc = 10e-3;			/* default max neurotrans conc */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcampa1;       /* default rate function */
   parm[0].funcname = (char *)"calcampa3";      /* user rate function */
   parm[0].dq[0] = dqsyn;                  /* Q10 for m rate function */

   respamp[GLU-GLU] = 0.8;		/* response to glutamate */
   respamp[AMPA-GLU] = 1.0;		/* response to AMPA */
   respamp[NMDA-GLU] = 1e-3;		/* response to NMDA */
   respamp[CNQX-GLU] = -1.0;		/* response to CNQX */

   a=1e7;	/* Kf alpha: forward */
   b=1.0;	/* Kr beta: reverse  */
   c=1.0;	/* desensitize f,r */
   d=1.0;	/* noise flicker f,r */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 3 * a;	/* Kf */
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 1.0 * a;	/* K'f */
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 300 * b;  /* Kr */
   spnt[1].rateo  [1] = 0;
   spnt[1].trans  [2] = 5;
   spnt[1].trate  [2] = rt;
   spnt[1].ratemul[2] = 1000.0 * c;	/* Lf */
   spnt[1].rateo  [2] = 1;

   spnt[2].numtrans   = 4;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rt;
   spnt[2].ratemul[0] = 3000.0 * d;	/* IIf */
   spnt[2].rateo  [0] = 3;		/* taud */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 3e5 * b;	/* K'r */
   spnt[2].rateo  [1] = 0;
   spnt[2].trans  [2] = 4;
   spnt[2].trate  [2] = rt;
   spnt[2].ratemul[2] = 6e4 * d;	/* If */
   spnt[2].rateo  [2] = 3;		/* taud */
   spnt[2].trans  [3] = 6;
   spnt[2].trate  [3] = rt;
   spnt[2].ratemul[3] = 2.7e4 * c;	/* Mf */
   spnt[2].rateo  [3] = 1;

   spnt[3].numtrans   = 1;
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 2;			/* state 3 = open state */
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 350 * d;	/* IIr */
   spnt[3].rateo  [0] = 3;

   spnt[4].numtrans   = 1;
   spnt[4].cond       = 1.0;
   spnt[4].trans  [0] = 2;			/* state 4 = open state */
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = 3000 * d;	/* Ir */
   spnt[4].rateo  [0] = 3;

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 1;			/* state 5 = inactivated state*/
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = 300 * c;	/* Lr */
   spnt[5].rateo  [0] = 1;
   spnt[5].trans  [1] = 6;
   spnt[5].trate  [1] = rnt;
   spnt[5].ratemul[1] = 1.0 * a;	/* K''f */
   spnt[5].rateo  [1] = 0;

   spnt[6].numtrans   = 3;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 7;	
   spnt[6].trate  [0] = rnt;
   spnt[6].ratemul[0] = 1.0 * a;	/* K'''f */
   spnt[6].rateo  [0] = 0;
   spnt[6].trans  [1] = 5;
   spnt[6].trate  [1] = rt;
   spnt[6].ratemul[1] = 519 * b;	/* K''r */
   spnt[6].rateo  [1] = 0;
   spnt[6].trans  [2] = 2;
   spnt[6].trate  [2] = rt;
   spnt[6].ratemul[2] = 14 * c;		/* Mr */
   spnt[6].rateo  [2] = 1;

   spnt[7].numtrans   = 2;
   spnt[7].cond       = 0;
   spnt[7].trans  [0] = 8;
   spnt[7].trate  [0] = rt;
   spnt[7].ratemul[0] = 6.0 * d;	/* IIIf */
   spnt[7].rateo  [0] = 3;
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = rt;
   spnt[7].ratemul[1] = 220 * 3 * b;	/* K'''r */
   spnt[7].rateo  [1] = 0;

   spnt[8].numtrans   = 1;
   spnt[8].cond       = 1.0;
   spnt[8].trans  [0] = 7;			/* state 8 = open state */
   spnt[8].trate  [0] = rt;
   spnt[8].ratemul[0] = 2000 * d;	/* IIIr */
   spnt[8].rateo  [0] = 3;

   return ch;
}

/*--------------------------------------------*/

chantype *makampa4(void)

/* number states from 0 to n; 
   set numstate = n;
   set numtrans = number of transitions.
   set cond     = conductance of state;

   for each state, set transitions:
    set trans = state to go to on this transition. 
    set trate = function that returns basic rate for transition. 
    set ratemul = multiplier for rate function.
    set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.

     0 <-> 1 <-> 2 <-> 3
                 |     |
                 4 <-> 5

*/

  {
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double a,b,c,d;
     double pca;
     int nstate, nparm, nq;

   nstate = 6;                          /* 6 Markov states */
   nparm = 1;                           /* make 2 sets of params, m, h */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(AMPA,4,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   pca = dpcaampa;
   ch->ions->ionp[PCA] = pca;		/* permeability to Ca++ ions */
   ch->ions->ionp[PNA] = 1.0;		/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0;		/* permeability to K+ ions */

   ch->vrev = vglu;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   ch->unitary = dampau;
   ch->trconc = dstr;
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcampa1;       /* default rate function */
   parm[0].funcname = (char *)"calcampa4";      /* user rate function */
   parm[0].dq[0] = dqsyn;                  /* Q10 for m rate function */

   respamp[GLU-GLU] = 0.8;		/* response to glutamate */
   respamp[AMPA-GLU] = 1.0;		/* response to AMPA */
   respamp[NMDA-GLU] = 1e-3;		/* response to NMDA */
   respamp[CNQX-GLU] = -1.0;		/* response to CNQX */

   a=9e6;	/* alpha: forward */
   b=6e3;	/* beta: reverse  */
   c=1.0;
   d=1.0;

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 3.0 * a;
   spnt[0].rateo  [0] = 0;	/* offsm */

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 2.0 * a;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 1.0 * b;
   spnt[1].rateo  [1] = 0;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rnt;
   spnt[2].ratemul[0] = 1.0 * a;
   spnt[2].rateo  [0] = 0;	/* offsm */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 2.0 * b;
   spnt[2].rateo  [1] = 0;	/* offsm */
   spnt[2].trans  [2] = 4;
   spnt[2].trate  [2] = rt;
   spnt[2].ratemul[2] = 2.0 * c;
   spnt[2].rateo  [2] = 1;	/* offsh */

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 5;			/* state 3 = the open state */
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 1e3 * c;
   spnt[3].rateo  [0] = 1;	/* offsh */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = 300 * b;
   spnt[3].rateo  [1] = 0;	/* offsm */

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rnt;
   spnt[4].ratemul[0] = 8.0 * a;	
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 2;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = 0.5 * d;
   spnt[4].rateo  [1] = 1;

   spnt[5].numtrans   = 2;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 3;			/* state 5 = inactivated state*/
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = 4.0 * c;
   spnt[5].rateo  [0] = 1;	/* offsh */
   spnt[5].trans  [1] = 4;
   spnt[5].trate  [1] = rt;
   spnt[5].ratemul[1] = 100 * c;
   spnt[5].rateo  [1] = 0;	/* offsm */

   return ch;
}

/*--------------------------------------------*/

chantype *makampa5(void)

/* number states from 0 to n; 
   set numstate = n;
   set numtrans = number of transitions.
   set cond     = conductance of state;

   for each state, set transitions:
    set trans = state to go to on this transition. 
    set trate = function that returns basic rate for transition. 
    set ratemul = multiplier for rate function.
    set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.
*/

/* Same as "ampa type 1" above, but without desensitization */

/* This channel's behavior is duplicated by running ampa type 1 */
/*  with tauc=100. */
 
/* Taken from:
                Jonas, P., Major, G., and Sakmann, B. (1993) 
                J. Physiol. 472:615-5663.

   in which a 7-state AMPA markov diagram is provided.  Two sets
of rate constants were derived, 1) which gave a good fit, and 2)
which was derived assuming independence of binding between
binding of glutamate molecules.  The one defined here is the "good
fit", not assuming independence.

     k1f       k2f        a
  0  <->   1   <->   2   <->   3
     k1r       k2r        b
                                    
State 0 is unbound
      1 is partially bound
      2 is fully bound
      3 is bound and open
      4,5,6 are deactivated
 
*/

  {
     double rnt(chan *cp),rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double k, pca;
     int nstate, nparm, nq;

#define ampa1f1 4.59e6
#define ampa1r1 4.26e3
#define ampa1f2 28.4e6
#define ampa1r2 3.26e3
#define ampa1f3 1.27e6
#define ampa1r3 45.7
#define ampa1a  4.24e3
#define ampa1b  900

   nstate = 4;                          /* 7 Markov states */
   nparm = 1;                           /* make 1 set of params */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(AMPA,5,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   ch->unitary = dampau;
   ch->trconc = dstr;
   pca = dpcaampa;
   ch->ions->ionp[PCA] = pca;		/* permeability to Ca++ ions */
   ch->ions->ionp[PNA] = 1.0;		/* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0;		/* permeability to K+ ions */

   ch->vrev = vglu;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcampa1;       /* default rate function */
   parm[0].funcname = (char *)"calcampa1";      /* user rate function */
   parm[0].dq[0] = dqsyn;               /* Q10 for m rate function */

   respamp[GLU-GLU] = 1.0;		/* response to glutamate */
   respamp[AMPA-GLU] = 0.5;		/* response to AMPA */
   respamp[NMDA-GLU] = 1e-3;		/* response to NMDA */
   respamp[CNQX-GLU] = -1.0;		/* response to CNQX */

   k=1;
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = ampa1f1 * k;
   spnt[0].rateo  [0] = 0;	/* offsm */

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 0; 
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = ampa1f2 * k;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = ampa1r1 * k;
   spnt[1].rateo  [1] = 1;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rt;		/* function of time */
   spnt[2].ratemul[0] = ampa1a * k;
   spnt[2].rateo  [0] = 0;	/* offsm */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;		/* function of nt conc, time */
   spnt[2].ratemul[1] = ampa1r2 * k;
   spnt[2].rateo  [1] = 1;	/* offsm */

   spnt[3].numtrans   = 1;
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 2;
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = ampa1b * k;
   spnt[3].rateo  [0] = 1;	/* offsm */

   return ch;

#undef ampa1f1 
#undef ampa1r1
#undef ampa1f2
#undef ampa1r2
#undef ampa1f3
#undef ampa1r3
#undef ampa1a 
#undef ampa1b 
#undef ampa1a1
#undef ampa1b1
#undef ampa1a2
#undef ampa1b2
#undef ampa1a3
#undef ampa1b3
#undef ampa1a4
#undef ampa1b4
}

/*--------------------------------------------*/

