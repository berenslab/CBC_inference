/* segment chank4 in program nc */
/* Ih, the hyperpolarization-activated inward rectifier */

/* sets up channel parameters */

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "ncsub.h"
#include "ncomp.h"
#include "nconst.h"
#include "control.h"

#ifdef __cplusplus
extern "C" {
#endif

double log(double);
double exp(double);

#ifdef __cplusplus
}
#endif

chantype *makchantype(int ctype, int cnum, int nstates, int nparm, 
				int nq, double bt);
void chanrate(chan *cpnt, chanparm *chp);

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
int set_perms (chantype *chtyp, iontab *ions, double svrev);

double rnt(chan *spnt);
double rt(chan *spnt);

/*----------------------------------------*/

/* HCN2 Markov model from Kusch et al., 2011 */
/* kc2o2 from Benndorf et al., 2012 */

/*  preactivated by step to -130 mV */

#define kc0c1 5.4e6  // s-1M-1
#define kc1c0 3.7e0  // s-1
#define kc1c2 8.4e5  // s-1M-1
#define kc2c1 9.3e-2 // s-1
#define kc2c3 9.9e4  // s-1M-1
#define kc3c2 8.5e0  // s -1
#define kc3c4 2.2e7  // s-1M-1
#define kc4c3 8.5e0  // s -1
#define ko0o1 8.4e6  // s-1M-1
#define ko1o0 2.5e0  // s -1
#define ko1o2 1.4e6  // s-1M-1
#define ko2o1 1.5e-1 // s -1
#define ko2o3 1.0e6  // s-1M-1
#define ko3o2 4.6e0  // s -1
#define ko3o4 2.2e7  // s-1M-1
#define ko4o3 8.5e0  // s -1
#define kc0o0 4.0e-2 // s -1
#define ko0c0 3.0e-2 // s -1
#define kc1o1 2.2e0  // s -1
#define ko1c1 7.1e-1 // s -1
#define kc2o2 1.6e-2 // s -1
#define ko2c2 4.6e-3 // s -1
#define kc3o3 4.0e1  //  > 2.0e1 -> 2.0×102 s -1
#define ko3c3 6.0e-1 //  > 3.0e-1 -> 3.0×100 s -1
#define kc4o4 4.0e1  //  > 2.0e1 -> 2.0×102 s -1
#define ko4c4 6.0e-1 //  > 3.0e-1 -> 3.0×10  s -1


/*--------------------------------------------*/

double calck10 (double v, int func) {}

/*----------------------------------------*/

chantype *makk10(void)

/* The HCN channel: Ih */

/*   Kinetic scheme and rate functions are taken from         */

/*   Kusch J1, Thon S, Schulz E, Biskup C, Nache V, Zimmer T, Seifert R, Schwede F, Benndorf K. (2011) */
/*      Nat Chem Biol. 2011 Dec 18;8(2):162-9. doi: 10.1038/nchembio.747. */

/*  Also see: */
 
/*   Benndorf K1, Thon S, Schulz E. (2012)                                   */
/*      Biophys J. 2012 Nov 7;103(9):1860-9. doi: 10.1016/j.bpj.2012.09.024. */

/*   Thon S1, Schmauder R, Benndorf K. (2013)                                */
/*      Biophys J. 2013 Oct 1;105(7):1581-9. doi: 10.1016/j.bpj.2013.08.027. */

/*   Young (2013) Biophys J. 105(7):1549-50. doi: 10.1016/j.bpj.2013.08.028. */


{
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double n;
     int nstate, nparm, nq;

   nstate = 10;                         /* 10 Markov states */
   nparm = 2;                           /* make 1 set of params, "n" */
   nq = 1;				/* make 1 Q10 value */
   ch=makchantype(K,10,nstate,nparm,nq,dbasetc); /* make chan state info */
   ch->unitary = 1.67e-12;		/* from Thon et al., 2013 (above) */
   ch->ions->ionp[PK]  = 1.0;		/* permeability to K+ ions */
   ch->ions->ionp[PNA] = 0.48;		/* permeability to Na+ ions */
   ch->ions->ionp[PCA] = dpcak;		/* permeability to Ca++ ions */
   ch->vrev = -0.03;			/* default reversal potential */
   set_perms (ch, ch->ions, ch->vrev);  /* fine tune permeabilities from vrev */

   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 0;                /* number of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calck10;         /* default rate function */
   parm[0].funcname = (char *)"calck10"; /* user rate function */
   parm[0].dq[0]    = dqm;              /* Q10 for m rate function */
   parm[0].voff = 0;                    /* voltage offset from user */

   parm[1].nfval    = 0;                /* number of func vals (am, bm, etc.) */
   parm[1].nival    = 0;                /* don't use implicit vals (yet) */
   parm[1].chancalc =  calck10;         /* default rate function */
   parm[1].funcname = (char *)"calck10";/* user rate function */
   parm[1].dq[0]    = dqm;              /* Q10 for m rate function */
   parm[1].voff = 0;                    /* voltage offset from user */

   n = 1.0;
   spnt[0].numtrans   = 2;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = rnt;
   spnt[0].ratemul[0] = kc0c1*n;
   spnt[0].rateo  [0] = 1;
   spnt[0].trans  [1] = 5;
   spnt[0].trate  [1] = rt;
   spnt[0].ratemul[1] = kc0o0*n;
   spnt[0].rateo  [1] = 0;

			     /*          kc0c1     kc1c2     kc2c3     kc3c4                     */
                             /*       0 <-----> 1 <-----> 2 <-----> 3 <-----> 4   closed states  */
			     /*          kc1c0     kc2c1     kc3c2     kc4c3                     */
			     /*      /|        /|        /|        /|        /|                  */
			     /*  ko0c0|    ko1c1|    ko2c2|    ko3c3|    ko4c4|                  */
			     /*       |kc0o0    |kc1o1    |kc2o2    |kc3o3    |kc4o4             */
			     /*       |         |         |         |         |                  */
			     /*       |/        |/        |/        |/        |/                 */
			     /*          ko0o1     ko1o2      ko2o3     ko3o4                    */
                             /*       5 <-----> 6 <-----> 7  <-----> 8 <-----> 9   open states   */
			     /*          ko1oO     ko2o1      ko3o2     ko4o3                    */
   spnt[1].numtrans   = 3;
   spnt[1].cond       = 0;
   spnt[1].trans  [0] = 2;	  	
   spnt[1].trate  [0] = rnt;
   spnt[1].ratemul[0] = kc1c2*n;
   spnt[1].rateo  [0] = 1;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = rt;
   spnt[1].ratemul[1] = kc1c0*n;
   spnt[1].rateo  [1] = 1;
   spnt[1].trans  [2] = 6;
   spnt[1].trate  [2] = rt;
   spnt[1].ratemul[2] = kc1o1*n;
   spnt[1].rateo  [2] = 0;

   spnt[2].numtrans   = 3;
   spnt[2].cond       = 0;
   spnt[2].trans  [0] = 3;	  	
   spnt[2].trate  [0] = rnt;
   spnt[2].ratemul[0] = kc2c3*n;
   spnt[2].rateo  [0] = 1;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = rt;
   spnt[2].ratemul[1] = kc2c1*n;
   spnt[2].rateo  [1] = 1;
   spnt[2].trans  [2] = 7;
   spnt[2].trate  [2] = rt;
   spnt[2].ratemul[2] = kc2o2*n;
   spnt[2].rateo  [2] = 0;

   spnt[3].numtrans   = 3;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 4;	  	
   spnt[3].trate  [0] = rnt;
   spnt[3].ratemul[0] = kc3c4*n;
   spnt[3].rateo  [0] = 1;
   spnt[3].trans  [1] = 2;
   spnt[3].trate  [1] = rt;
   spnt[3].ratemul[1] = kc3c2*n;
   spnt[3].rateo  [1] = 1;
   spnt[3].trans  [2] = 8;
   spnt[3].trate  [2] = rt;
   spnt[3].ratemul[2] = kc3o3*n;
   spnt[3].rateo  [2] = 0;

   spnt[4].numtrans   = 2;
   spnt[4].cond       = 0;
   spnt[4].trans  [0] = 3;	  	
   spnt[4].trate  [0] = rt;
   spnt[4].ratemul[0] = kc4c3*n;
   spnt[4].rateo  [0] = 1;
   spnt[4].trans  [1] = 9;
   spnt[4].trate  [1] = rt;
   spnt[4].ratemul[1] = kc4o4*n;
   spnt[4].rateo  [1] = 0;
			     /*          kc0c1     kc1c2     kc2c3     kc3c4                     */
                             /*       0 <-----> 1 <-----> 2 <-----> 3 <-----> 4   closed states  */
			     /*          kc1c0     kc2c1     kc3c2     kc4c3                     */
			     /*      /|        /|        /|        /|        /|                  */
			     /*  ko0c0|    ko1c1|    ko2c2|    ko3c3|    ko4c4|                  */
			     /*       |kc0o0    |kc1o1    |kc2o2    |kc3o3    |kc4o4             */
			     /*       |         |         |         |         |                  */
			     /*       |/        |/        |/        |/        |/                 */
			     /*          ko0o1     ko1o2      ko2o3     ko3o4                    */
                             /*       5 <-----> 6 <-----> 7  <-----> 8 <-----> 9   open states   */
			     /*          ko1oO     ko2o1      ko3o2     ko4o3                    */

   spnt[5].numtrans   = 2;		/* open state */
   spnt[5].cond       = 1;
   spnt[5].trans  [0] = 6;
   spnt[5].trate  [0] = rnt;
   spnt[5].ratemul[0] = ko0o1*n;
   spnt[5].rateo  [0] = 1;
   spnt[5].trans  [1] = 0;
   spnt[5].trate  [1] = rt;
   spnt[5].ratemul[1] = ko0c0*n;
   spnt[5].rateo  [1] = 0;

   spnt[6].numtrans   = 3;		/* open state */
   spnt[6].cond       = 1;
   spnt[6].trans  [0] = 7;	  	
   spnt[6].trate  [0] = rnt;
   spnt[6].ratemul[0] = ko1o2*n;
   spnt[6].rateo  [0] = 1;
   spnt[6].trans  [1] = 5;
   spnt[6].trate  [1] = rt;
   spnt[6].ratemul[1] = ko1o0*n;
   spnt[6].rateo  [1] = 1;
   spnt[6].trans  [2] = 1;
   spnt[6].trate  [2] = rt;
   spnt[6].ratemul[2] = ko1c1*n;
   spnt[6].rateo  [2] = 0;

   spnt[7].numtrans   = 3;		/* open state */
   spnt[7].cond       = 1;
   spnt[7].trans  [0] = 8;	  	
   spnt[7].trate  [0] = rnt;
   spnt[7].ratemul[0] = ko2o3*n;
   spnt[7].rateo  [0] = 1;
   spnt[7].trans  [1] = 6;
   spnt[7].trate  [1] = rt;
   spnt[7].ratemul[1] = ko2o1*n;
   spnt[7].rateo  [1] = 1;
   spnt[7].trans  [2] = 2;
   spnt[7].trate  [2] = rt;
   spnt[7].ratemul[2] = ko2c2*n;
   spnt[7].rateo  [2] = 0;

   spnt[8].numtrans   = 3;		/* open state */
   spnt[8].cond       = 1;
   spnt[8].trans  [0] = 9;	  	
   spnt[8].trate  [0] = rnt;
   spnt[8].ratemul[0] = ko3o4*n;
   spnt[8].rateo  [0] = 1;
   spnt[8].trans  [1] = 7;
   spnt[8].trate  [1] = rt;
   spnt[8].ratemul[1] = ko3o2*n;
   spnt[8].rateo  [1] = 1;
   spnt[8].trans  [2] = 3;
   spnt[8].trate  [2] = rt;
   spnt[8].ratemul[2] = ko3c3*n;
   spnt[8].rateo  [2] = 0;

   spnt[9].numtrans   = 2;		/* open state */
   spnt[9].cond       = 1;
   spnt[9].trans  [0] = 8;	  	
   spnt[9].trate  [0] = rt;
   spnt[9].ratemul[0] = ko4o3*n;
   spnt[9].rateo  [0] = 1;
   spnt[9].trans  [1] = 4;
   spnt[9].trate  [1] = rt;
   spnt[9].ratemul[1] = ko4c4*n;
   spnt[9].rateo  [1] = 0;

   return ch;
}

#undef kc0c1 
#undef kc1c0 
#undef kc1c2 
#undef kc2c1 
#undef kc2c3
#undef kc3c2
#undef kc3c4
#undef kc4c3
#undef ko0o1
#undef kO1o0
#undef ko1o2
#undef ko2o1
#undef ko2o3
#undef ko3o2
#undef ko3o4
#undef ko4o3
#undef kc0o0
#undef ko0c0
#undef kc1o1
#undef ko1c1
#undef kc2o2
#undef ko2c2
#undef kc3o3
#undef ko3c3
#undef kc4o4
#undef ko4c4

