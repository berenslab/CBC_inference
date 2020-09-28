/* segment chanchr2 in program nc */

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
double tanh(double);

#ifdef __cplusplus
}
#endif

chantype *makchantype(int ctype, int cnum, int nstates, 
		int nparm, int nq, double bt);
double qraten (chanparm *chp, int t);

double rnt(chan *spnt);
double rt(chan *spnt);

/*--------------------------------------------*/

double calcchr2(double v, int func) {}

/*--------------------------------------------*/

#define e1 0.8535
#define e2 0.14

// #define e12d 0.011
#define e12d 0.051
#define c1 0.005
#define c2 0.024

// #define e21d 0.008
#define e21d 0.001
#define c3 0.004
#define c4 0.024

#define tchr2 1.3	/* time const msec */
#define o2gamma 0.1	/* weight for O2 conductance, see makphotrec() in ncmak.cc */

/* note that qraten() multiplies by timinc */
/* see changes (below) to gd1, gd2, and gr that reduce gr rate and increase initial nose of response */
 
/* - - - - - - - - - - - - - - - - - - - - - - */

double chr2_p(chan *cpnt) 

/* calculate integrated flux "p", called once per timinc step */

{
   double p,s;

   p = ((chrchan*)cpnt)->p;
   s = 0.5 * (1 + tanh(120 * (((chrchan*)cpnt)->iflux*100.0 - 0.1)));
   p += (s - p) / tchr2 * MSSEC * timinc;
   ((chrchan*)cpnt)->p = p;
   return p;
}

/* - - - - - - - - - - - - - - - - - - - - - - */

double k1,k2,gd1,gd2,e12, e21,gr;

double chr2_k1(chan *cpnt) 
{
	double p;
   p = chr2_p(cpnt);		/* calculate p, once per time step */
// fprintf (stderr,"c k1 %g k2 %g gd1 %g gd2 %g e12 %g e21 %g gr %g iflux %g q %g cond %g C1 %g O1 %g O2 %g C2 %g\n",k1,k2,gd1,gd2,e12,e21,gr,((chrchan*)cpnt)->iflux,qraten(cpnt->chtyp->parm,0)*MSSEC,cpnt->conduct,cpnt->conc[0].cest,cpnt->conc[1].cest,cpnt->conc[2].cest,cpnt->conc[3].cest);
   k1 = e1 * ((chrchan*)cpnt)->iflux * ((chrchan*)cpnt)->p * MSSEC * qraten(cpnt->chtyp->parm,0);
   return k1;
}

/* - - - - - - - - - - - - - - - - - - - - - - */

double chr2_k2(chan *cpnt) 
{
   k2 = e2 * ((chrchan*)cpnt)->iflux * ((chrchan*)cpnt)->p * MSSEC  * qraten(cpnt->chtyp->parm,1);
   return k2;
}

/* - - - - - - - - - - - - - - - - - - - - - - */

double chr2_gd1(chan *cpnt) 
{
   // gd1 = (0.075 + 0.043 * tanh((cpnt->comp1->v*VTOMV + 20) / -20)) * MSSEC * qraten(cpnt->chtyp->parm,2);
   gd1 = (0.025 + 0.043 * tanh((cpnt->comp1->v*VTOMV + 20) / -20)) * MSSEC * qraten(cpnt->chtyp->parm,2);
   return gd1;
}

/* - - - - - - - - - - - - - - - - - - - - - - */

double chr2_gd2(chan *cpnt) 
{
// #define gd2_val 0.05
#define gd2_val 0.001
   gd2 = gd2_val * MSSEC * qraten(cpnt->chtyp->parm,3);
   return gd2;
}

/* - - - - - - - - - - - - - - - - - - - - - - */

double chr2_e12(chan *cpnt) 
{
   e12 =  (e12d + c1 * log(1.0 + ((chrchan*)cpnt)->iflux / c2)) * MSSEC * qraten(cpnt->chtyp->parm,4);
   return e12;
}

/* - - - - - - - - - - - - - - - - - - - - - - */

double chr2_e21(chan *cpnt) 
{
   e21 =  (e21d + c3 * log(1.0 + ((chrchan*)cpnt)->iflux / c4)) * MSSEC * qraten(cpnt->chtyp->parm,5);
   return e21;
}

/* - - - - - - - - - - - - - - - - - - - - - - */

double chr2_gr(chan *cpnt) 
{
   // gr = (4.34587e5 * exp (-0.0211539274*cpnt->comp1->v*VTOMV)) * MSSEC * qraten(cpnt->chtyp->parm,6);
   gr = (4.34587e2 * exp (-0.0211539274*cpnt->comp1->v*VTOMV)) * MSSEC * qraten(cpnt->chtyp->parm,6);
   return gr;
}

/*--------------------------------------------*/

chantype *makchr2(void)

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

/* channel-rhodopsin phototransduction from Williams et al (2013) */

/* Williams JC, Xu J, Lu Z, Klimas A, Chen X, Ambrosi CM, Cohen IS, Entcheva E.
   (2013) Computational optogenetics: empirically-derived voltage- and
   light-sensitive channelrhodopsin-2 model. PLoS Comput Biol. 2013;9(9):e1003220.
   doi: 10.1371/journal.pcbi.1003220. https://www.ncbi.nlm.nih.gov/pubmed/24068903
*/

/*  Channel-Rhodopsin 2 (ChR2) channel 

    Four states:

                e12
        O1     <--->     O2
                e21
      ^  |             ^ |
      k1 | Gd1        k2 | Gd2
           \/              \/
                Gr
        C1     <---     C2


        1     <--->     2

      ^ |             ^ |
        | \/            | \/

        0     <---      3             4 = "p parameter"
*/

{
     double rnt(chan *cp),rt(chan *cp);     /* rate functions */
     chantype *ch;
     chanstate *spnt;
     chanparm *parm;
     double a,b,c,d;
     int nstate, nparm, nq;


   nstate = 4;                          /* 5 Markov states */
   nparm = 1;                           /* make 1 set of params, m, h */
   nq = 7;				/* number of Q10 values */
   ch=makchantype(CHRC,ChR2,nstate,nparm,nq,dbasetsyn);/* make chan const, parm */
   ch->ions->ionp[PCL] = 0;			/* permeability to Cl- ions */
   ch->ions->ionp[PCA] = 0;             /* permeability to Ca++ ions */
   ch->ions->ionp[PNA] = 1.0;           /* permeability to Na+ ions */
   ch->ions->ionp[PK]  = 1.0;           /* permeability to K+ ions */

   ch->unitary = dpru;
   ch->gamma = 0.1;			/* weight for O2 cond */
   ch->vrev = 0;			/* default reversal potential */
   spnt = ch->state;
   parm = ch->parm;

   parm[0].nfval    = 0;                /* num of func vals (am, bm, etc.) */
   parm[0].nival    = 0;                /* don't use implicit vals (yet) */
   parm[0].chancalc =  calcchr2;       /* default rate function */
   parm[0].funcname = (char *)"calcchr2";   /* user rate function */

   /*
   cpt->q10_e1  = qcalc(1.46);
   cpt->q10_e2  = qcalc(2.77);
   cpt->q10_gd1 = qcalc(1.97); 
   cpt->q10_gd2 = qcalc(1.77);
   cpt->q10_e12 = qcalc(1.1); 
   cpt->q10_e21 = qcalc(1.95);
   cpt->q10_gr  = qcalc(2.56);
   */
   parm[0].dq[0] = 1.46;		/* Q10 for C->O rate functions */
   parm[0].dq[1] = 2.77;
   parm[0].dq[2] = 1.97;		/* Q10 for O->rate functions */
   parm[0].dq[3] = 1.77;
   parm[0].dq[4] = 1.1; 		/* Q10 for forward-reverse functions */
   parm[0].dq[5] = 1.95;
   parm[0].dq[6] = 2.56;

   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = chr2_k1;
   spnt[0].ratemul[0] = 1;
   spnt[0].rateo  [0] = 0;

   spnt[1].numtrans   = 2;
   spnt[1].cond       = 1.0;
   spnt[1].trans  [0] = 2;
   spnt[1].trate  [0] = chr2_e12;
   spnt[1].ratemul[0] = 1;
   spnt[1].rateo  [0] = 0;
   spnt[1].trans  [1] = 0;
   spnt[1].trate  [1] = chr2_gd1;
   spnt[1].ratemul[1] = 1;
   spnt[1].rateo  [1] = 0;

   spnt[2].numtrans   = 2;
   spnt[2].cond       = ch->gamma;
   spnt[2].trans  [0] = 3;
   spnt[2].trate  [0] = chr2_gd2;
   spnt[2].ratemul[0] = 1;
   spnt[2].rateo  [0] = 0;
   spnt[2].trans  [1] = 1;
   spnt[2].trate  [1] = chr2_e21;
   spnt[2].ratemul[1] = 1;
   spnt[2].rateo  [1] = 0;	

   spnt[3].numtrans   = 2;
   spnt[3].cond       = 0;
   spnt[3].trans  [0] = 2;
   spnt[3].trate  [0] = chr2_k2;
   spnt[3].ratemul[0] = 1;
   spnt[3].rateo  [0] = 0;
   spnt[3].trans  [1] = 0;	
   spnt[3].trate  [1] = chr2_gr;
   spnt[3].ratemul[1] = 1;
   spnt[3].rateo  [1] = 0;

   return ch;

}

/*--------------------------------------------*/

