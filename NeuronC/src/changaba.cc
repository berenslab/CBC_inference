/* segment changaba in program nc */

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

/*--------------------------------------------*/

double calcgaba1(double v, int func) {}

/*--------------------------------------------*/

chantype *makgaba1(void)

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

		Busch and Sakmann (1990) 
			[cited by Destexhe, Mainen and Sejnowski (1998),
                              Kinetic Models of Synaptic Transmission,
                              in Methods in Neuronal Modeling, 
                              Ed. by Koch and Segev ]

   in which a 5-state GABA-A Markov diagram is provided:

                  gaba      gaba
                0  <->  1   <->   2

                        |         |

                        3         4     Open states

*/

  {
     double rnt(chan *cp),rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double a,b,c,d;
     int nstate, nparm, nq;


   nstate = 5;                          /* 5 Markov states */
   nparm = 1;                           /* make 1 set of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(GABA,1,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   ch->ions->ionp[PCL] = 1.0;			/* permeability to Cl- ions */
   ch->ions->ionp[PK]  = 0.0;			/* permeability to K+ ions */
   ch->unitary = dgabau;
   ch->trconc = 100e-6;			/* default max trans conc */
   ch->vrev = vcl;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcgaba1;       /* default rate function */
   parm[0].funcname = (char *)"calcgaba1";      /* user rate function */
   parm[0].dq[0] = dqsyn;               /* Q10 for m rate function */

   respamp[GABA-GLU] = 1.0;		/* response to GABA */
   respamp[BIC-GLU] = -1.0;		/* response to bicucculine */
   respamp[PTX-GLU] = -1.0;		/* response to picrotoxin */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 20e6;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 10e6;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 4.6e3;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 3;
   spnt[1].trate  [2] = rt;
   spnt[1].ratemul[2] = 3.3e3;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 4;
   spnt[2].trate  [0] = rt;
   spnt[2].ratemul[0] = 10.6e3;
   spnt[2].rateo  [0] = 2;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 9.2e3;
   spnt[2].rateo  [1] = 1;	

   spnt[3].numtrans   = 1;
   spnt[3].cond       = 1.0;			/* state 3 = open state */
   spnt[3].trans  [0] = 1;
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 9.8e3;
   spnt[3].rateo  [0] = 3;

   spnt[4].numtrans   = 1;
   spnt[4].cond       = 1.0;			/* state 4 = open state */
   spnt[4].trans  [0] = 2;
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = 410;	
   spnt[4].rateo  [0] = 3;

   return ch;

}
/*--------------------------------------------*/

chantype *makgaba2(void)

/* number states from 0 to n; 
   set numstate = n;
   set numtrans = number of transitions.
   set cond     = conductance of state;

   for each state, set transitions:
    set trans = state to go to on this transition. 
    set trate = function that returns basic rate for transition. 
    set ratemul = multiplier for rate function.
    set rateo  =  Set to 0 for m, set to 1 for h:  sets voffset and tau.

       k1f       k2f        a
    0  <->   1   <->   2   <->   3
       k1r       k2r        b
                                   
               ^        ^         ^
          a1 | b1   a2 | b2   a3 | b3
          \/        \/        \/

                  k3f        a4
             4   <->   5   <->   6
                  k3r        b4         

*/

/* Taken from:
                Jonas, P., Major, G., and Sakmann, B. (1993) 
                J. Physiol. 472:615-5663.

   in which a 7-state AMPA markov diagram is provided. 

   Simple GABA channel based on AMPA state diagram. Need to play
with kinetics (a,b,c,d).

*/

  {
     double rnt(chan *cp),rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double a,b,c,d;
     int nstate, nparm, nq;


   nstate = 7;                          /* 7 Markov states */
   nparm = 1;                           /* make 1 set of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(GABA,2,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   ch->ions->ionp[PCL] = 1.0;			/* permeability to Cl- ions */
   ch->ions->ionp[PK]  = 0;			/* permeability to K+ ions */
   ch->ions->ionp[PNA] = 0;			/* permeability to Na+ ions */
   ch->unitary = dgabau;
   ch->trconc = 100e-6;			/* default max trans conc */
   ch->vrev = vcl;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcgaba1;       /* default rate function */
   parm[0].funcname = (char *)"calcgaba2";      /* user rate function */
   parm[0].dq[0] = dqsyn;               /* Q10 for m rate function */

   respamp[GABA-GLU] = 1.0;		/* response to GABA */
   respamp[BIC-GLU] = -1.0;		/* response to bicucculine */
   respamp[PTX-GLU] = -1.0;		/* response to picrotoxin */

   a = 4.0e4;	/* forward: alpha */
   b = 2.0e2;   /* reverse: beta */
   c = 1.0e2;   /* deactivation */
   d = 5.0e2;   /* re-activation */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 1.0 * a;
   spnt[0].rateo  [0] = 0;	/* offsm */

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 1.0 * a;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 1.0 * b;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 4;
   spnt[1].trate  [2] = rt;
   spnt[1].ratemul[2] = 1.0 * c;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = rt;		/* function of time */
   spnt[2].ratemul[0] = 1.0 * a;
   spnt[2].rateo  [0] = 0;	/* offsm */
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;		/* function of nt conc, time */
   spnt[2].ratemul[1] = 1.0 * b;
   spnt[2].rateo  [1] = 1;	/* offsm */
   spnt[2].trans  [2] = 5;
   spnt[2].trate  [2] = rt;
   spnt[2].ratemul[2] = 1.0 * c;
   spnt[2].rateo  [2] = 2;	/* offsh */

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 1.0;
   spnt[3].trans  [0] = 6;			/* state 3 = the open state */
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 1.0 * c;
   spnt[3].rateo  [0] = 2;	/* offsh */
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = 1.0 * b;
   spnt[3].rateo  [1] = 1;	/* offsm */

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 5;
   spnt[4].trate  [0] = rnt;
   spnt[4].ratemul[0] = 1.0 * a;	
   spnt[4].rateo  [0] = 0;
   spnt[4].trans  [1] = 1;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = 1.0 * d;
   spnt[4].rateo  [1] = 3;

   spnt[5].numtrans   = 3;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 6;			/* state 5 = inactivated state*/
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = 1.0 * a;
   spnt[5].rateo  [0] = 0;	/* offsh */
   spnt[5].trans  [1] = 2;
   spnt[5].trate  [1] = rt;
   spnt[5].ratemul[1] = 1.0 * d;
   spnt[5].rateo  [1] = 3;	/* offsh */
   spnt[5].trans  [2] = 4;
   spnt[5].trate  [2] = rt;
   spnt[5].ratemul[2] = 1.0 * b;
   spnt[5].rateo  [2] = 1;	/* offsm */

   spnt[6].numtrans   = 2;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 3;			/* state 6 = inactivated state*/
   spnt[6].trate  [0] = rt;
   spnt[6].ratemul[0] = 1.0 * d;
   spnt[6].rateo  [0] = 3;	/* offsh */
   spnt[6].trans  [1] = 5;
   spnt[6].trate  [1] = rt;
   spnt[6].ratemul[1] = 1.0 * b;
   spnt[6].rateo  [1] = 1;	/* offsm */

   return ch;

}


/*--------------------------------------------*/

chantype *makgaba3(void)

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

    Jones MV, Westbrook GL (1995), Desensitized states prolong GABA-A channel responses to brief agonist pulses. 
      Neuron 15:181-191.

in wich a 7-state GABA-A Markov diagram is provided:

                Dslow         Dfast
              d1 /| r1      d2 /| r2
                  |             |
        2kon      |/    kon     |/
Unbound <--->   Bound1 <--->  Bound2
        kOff     /|    2kOff   /|
               a1 | b1       a2 | b2
                  |/            |/
                Open1         Open2


                  5             6

                 /|            /| 
              d1  | r1      d2  | r2
        2KON      |/    KON     |/

    0   <--->     1    <--->    2

        KOFF     /|    2KOFF   /|
               a1 | b1       a2 | b2
                  |/            |/

                  3             4

 kinetic rate constants:

 kon (/uM/s)        3
 koff (s-1)       150
 Bl (s-1)         200
 a1 (s-1)        1111
 B2 (s-1)        2500
 a2 (s-l)         142
 d1 (s-l)          13
 r1 (s-l)        0.13
 d2 (s-1)    750-1000
 r2 (s-1)       15-25
 Po (max)     0.6-0.8

*/
     

#define KON      3
#define KOFF   150
#define d1      13
#define d2     800
#define r1    0.13
#define r2      20
#define a1    1111
#define a2     142
#define b1     200
#define b2    2500


/*
#define KON      3
#define KOFF   150
#define d1      0
#define d2     0
#define r1    0.0
#define r2      0
#define a1    0
#define a2     0
#define b1     0
#define b2    0
*/

  {
     double rnt(chan *cp),rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     int nstate, nparm, nq;

   nstate = 7;                          /* 5 Markov states */
   nparm = 1;                           /* make 1 sets of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(GABA,3,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   ch->ions->ionp[PCL] = 1.0;			/* permeability to Cl- ions */
   ch->ions->ionp[PK]  = 0.0;			/* permeability to K+ ions */
   ch->unitary = dgabau;
   ch->trconc = 100e-6;			/* default max trans conc */
   ch->vrev = vcl;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcgaba1;       /* default rate function */
   parm[0].funcname = (char *)"calcgaba3";      /* user rate function */
   parm[0].dq[0] = dqsyn;               /* Q10 for m rate function */

   respamp[GABA-GLU] = 1.0;		/* response to GABA */
   respamp[BIC-GLU] = -1.0;		/* response to bicucculine */
   respamp[PTX-GLU] = -1.0;		/* response to picrotoxin */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = KON * 2e6;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 4;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = KON * 1e6;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = KOFF;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 3;
   spnt[1].trate  [2] = rt;
   spnt[1].ratemul[2] = b1;
   spnt[1].rateo  [2] = 2;
   spnt[1].trans  [3] = 5;
   spnt[1].trate  [3] = rt;
   spnt[1].ratemul[3] = d1;
   spnt[1].rateo  [3] = 3;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 1.0;
   spnt[2].trans  [0] = 1;
   spnt[2].trate  [0] = rt;
   spnt[2].ratemul[0] = KOFF;
   spnt[2].rateo  [0] = 1;
   spnt[2].trans  [1] = 4;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = b2;
   spnt[2].rateo  [1] = 2;
   spnt[2].trans  [2] = 6;
   spnt[2].trate  [2] = rt;
   spnt[2].ratemul[2] = d2;
   spnt[2].rateo  [2] = 3;

   spnt[3].numtrans   = 1;
   spnt[3].cond       = 0.0;
   spnt[3].trans  [0] = 1;
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = a1;
   spnt[3].rateo  [0] = 2;

   spnt[4].numtrans   = 1;
   spnt[4].cond       = 0.0;
   spnt[4].trans  [0] = 2;
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = a2;
   spnt[4].rateo  [0] = 2;

   spnt[5].numtrans   = 1;
   spnt[5].cond       = 0;
   spnt[5].trans  [0] = 1;
   spnt[5].trate  [0] = rt;
   spnt[5].ratemul[0] = r1;
   spnt[5].rateo  [0] = 2;

   spnt[6].numtrans   = 1;
   spnt[6].cond       = 0;
   spnt[6].trans  [0] = 2;
   spnt[6].trate  [0] = rt;
   spnt[6].ratemul[0] = r2;
   spnt[6].rateo  [0] = 2;

#undef KON
#undef KOFF
#undef d1
#undef d2
#undef r1
#undef r2
#undef a1
#undef a2
#undef b1
#undef b2

}

/*--------------------------------------------*/

chantype *makgaba4(void)

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

		Busch and Sakmann (1990) 
			[cited by Destexhe, Mainen and Sejnowski (1998),
                              Kinetic Models of Synaptic Transmission,
                              in Methods in Neuronal Modeling, 
                              Ed. by Koch and Segev ]

   in which a 5-state GABA-A Markov diagram is provided:

                  gaba      gaba
                0  <->  1   <->   2

                        |         |

                        3         4     Open states


 Like GABA1 channel (above) except with slower deactivation for GABA-C 

 Test with nc/tcomp66b4

*/

  {
     double rnt(chan *cp),rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double *respamp;
     double a,b,c,d;
     int nstate, nparm, nq;

#define GABAC_drate 0.14

   nstate = 5;                          /* 5 Markov states */
   nparm = 1;                           /* make 1 set of params, m, h */
   nq = 1;				/* number of Q10 values */
   ch=makchantype(GABA,4,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   ch->ions->ionp[PCL] = 1.0;			/* permeability to Cl- ions */
   ch->ions->ionp[PK]  = 0.0;			/* permeability to K+ ions */
   ch->unitary = dgabau;
   ch->trconc = 100e-6;			/* default max trans conc */
   ch->vrev = vcl;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;
   respamp = ch->respamp;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcgaba1;       /* default rate function */
   parm[0].funcname = (char *)"calcgaba1";      /* user rate function */
   parm[0].dq[0] = dqsyn;               /* Q10 for m rate function */

   respamp[GABA-GLU] = 1.0;		/* response to GABA */
   respamp[BIC-GLU] = -1.0;		/* response to bicucculine */
   respamp[PTX-GLU] = -1.0;		/* response to picrotoxin */

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = 20e6;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = 10e6;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = 4.6e3 * GABAC_drate;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 3;
   spnt[1].trate  [2] = rt;
   spnt[1].ratemul[2] = 3.3e3;
   spnt[1].rateo  [2] = 2;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 4;
   spnt[2].trate  [0] = rt;
   spnt[2].ratemul[0] = 10.6e3;
   spnt[2].rateo  [0] = 2;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = 9.2e3 * GABAC_drate;
   spnt[2].rateo  [1] = 1;	

   spnt[3].numtrans   = 1;
   spnt[3].cond       = 1.0;			/* state 3 = open state */
   spnt[3].trans  [0] = 1;
   spnt[3].trate  [0] = rt;
   spnt[3].ratemul[0] = 9.8e3 * GABAC_drate;
   spnt[3].rateo  [0] = 3;

   spnt[4].numtrans   = 1;
   spnt[4].cond       = 1.0;			/* state 4 = open state */
   spnt[4].trans  [0] = 2;
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = 410 * GABAC_drate;	
   spnt[4].rateo  [0] = 3;

   return ch;

}

