/* Segment ncomp in Program nc */

/* Simulates neuronal circuits */
/*  by numerically integrating difference equations */
/*  using iterative relaxation technique */
/*  to implement the "implicit" integration method */ 

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "adef.h"
#include "nc.h"
#include "ytab.h"
#include "control.h"

#ifdef __cplusplus
}
#endif

#include "ncsub.h"
#include "ncomp.h"
#include "ncelem.h" 

/* #define DEBUG			/* */
#define NOSKIP 128

#define SQRT10 3.16227766
#define SSQRT10 1.77827942

extern comp *compnt,*compend;

#ifdef __cplusplus
extern "C" {
#endif
    double log(double);
    double log10(double);
    double exp(double);
    double sqrt(double);
    void exit(int err);
#ifdef __cplusplus
}
#endif

double ncabs(double x);

static double maxerr;
static int ncomp,tcomp;

double akcacalc(double v, double ca, double rate, double d1, double k1);
double bkcacalc(double v, double ca, double rate, double d2, double k2);

void docomp(comp *p);
void runsyn(int t);
void chanimpl(chanparm *chp, chan* cpnt);	
void kcatab(double voffs, double cais, int stype, double taua,
	     double d1, double d2, double k1, double k2, chanparm *chp);
char *rsetstate (char *state);
double binomdev(double pp, int n);
chan *alloc_chan(short ctype, short stype);
char *emalloc(unsigned int n);
double drand(void);
double gasdev(void);
chantype *getchantype(int ctype, int cnum);
double ghkv (chantype *chtyp,double cao, double cai);

/*------------------------------------*/

double dochan2(chan *chpnt)
              
/* Second order calculation of channel states. */
/*  Assumes "cest" contains old state concentration. */
/*  Leaves new estimate in "cest". */
/*  Returns conductance. */
/*  If "cest" goes out of range, assume that timestep */
/*   is too large.  Reduce timestep and try again. */

{
   register int i,j;
   int t,nstate,numtr,errcnt;
   static chanstate *stpnt;
   static stconc *conc,*mconc;
   static double (**f_rate)(chan *cpnt), trrate;
   static double conduct,critc,totconc,*tcrp;
   static float *trmul;
   static char *rateo;
   double trm;			/* multiplier for rate func */
   char *trans;
   double cerr,cmaxerr,cest;
   static dstconc oldconc[NUMSTATE];
   static dstconc *nconc,newconc[NUMSTATE];
   static double tcr[NUMSTATE][NUMTRANS]; /* temp store for trrate in calc */
   static int k=0;

#define MAXERRCNT 100

  t = chpnt->stype;
  nstate = chpnt->numstate;

   /* Set up rate constant multiplier */

  stpnt = chpnt->chtyp->state;
  for (i=0; i<nstate; i++,stpnt++) {  		/* zero deltas, set tcr */
    tcrp = tcr[i];				/* pointer to precalc values */
    f_rate = stpnt->trate;
    trmul = stpnt->ratemul;
    rateo = stpnt->rateo;
    numtr = stpnt->numtrans;			/*per timinc */
    trans = stpnt->trans;	     		/* get parameters */
    for (j=0; j<numtr; j++,trans++,f_rate++,trmul++,rateo++,tcrp++){ 
	switch (*rateo) {
	default:
	case 0: trm = chpnt->arate; break; 
	case 1: trm = chpnt->brate; break; 
	case 2: trm = chpnt->crate; break; 
	case 3: trm = chpnt->drate; break; 
        }
        *tcrp = (*f_rate)(chpnt) * *trmul * trm; /* save precalc val */
    }
  }	

  if (chpnt->nchan > 0) {			/* Noise on */

        int delch, nchano, totnchan,r,s,ii,jj,t;
	static int rstate1[NUMSTATE];
	static int rstate2[NUMSTATE];

    conc = chpnt->conc;
    stpnt = chpnt->chtyp->state;
    nconc = newconc;
    for (r=0; r<nstate; r++) { 		/* state order */
       rstate1[r] = r; 
    }
    for (r=0; r<nstate; r++) { 		/* make random state order */
       t = rstate1[r]; 
       s = int(drand()*nstate); 
       rstate1[r] = rstate1[s];
       rstate1[s] = t;
    }
    for (i=0; i<nstate; i++) { /* find change in conc */
      ii = rstate1[i];
      numtr = stpnt[ii].numtrans;
      trans = stpnt[ii].trans;                     /* get parameters */
      tcrp = tcr[ii];                            /* pointer to precalc values */
      for (r=0; r<numtr; r++) { 		/* state order */
         rstate2[r] = r; 
      }
      for (r=0; r<numtr; r++) { 		/* make random state order */
         t = rstate2[r]; 
         s = int(drand()*nstate); 
         rstate2[r] = rstate2[s];
         rstate2[s] = t;
      }
      for (j=0; j<numtr; j++){
        jj = rstate2[j];
        if (chpnt->cstate) rsetstate (chpnt->cstate);
        delch = int(binomdev(tcrp[jj],conc[ii].nchan)); 
	  conc[trans[jj]].nchan += delch;
	  conc[ii].nchan -= delch;
//fprintf (stderr,"state %d tran %d pop %d delch %d\n",i,j,conc->nchan,delch);
      }
//fprintf (stderr,"\n");
    }

/*
    for (i=0; i<nstate; i++,stpnt++,conc++,nconc++) { /* find change in conc */
/*      numtr = stpnt->numtrans;
      trans = stpnt->trans;                     /* get parameters */
/*      tcrp = tcr[i];                            /* pointer to precalc values */
/*      for (j=0; j<numtr; j++,trans++,tcrp++){     
        if (chpnt->cstate) rsetstate (chpnt->cstate);
        delch = int(binomdev(*tcrp,conc->nchan)); 
//        delch = int(binomdev(*tcrp, int(chpnt->nchan*conc->cest))); 
//        newconc[*trans].dchan += delch;
//        nconc->dchan -= delch;
	chpnt->conc[*trans].nchan += delch;
	  conc->nchan -= delch;
//fprintf (stderr,"state %d tran %d pop %d delch %d\n",i,j,conc->nchan,delch);
      }
//fprintf (stderr,"\n");
    }
*/

//fprintf (stderr,"\n");
    conc = chpnt->conc;
    nconc = newconc;
    stpnt = chpnt->chtyp->state;
    nchano = 0;
    totnchan = 0;
    for (i=0; i<0;/*nstate;*/ i++,stpnt++,conc++,nconc++) {     /* Total */
       conc->nchan += nconc->dchan;
       if (conc->nchan < 0) {		/* if nchan is negative */	
             int maxtr;
             double maxval;
//fprintf (stderr,"nchan neg state %d val %d\n",i,conc->nchan);
          maxtr = 0;
          maxval = 0.0;
          numtr = stpnt->numtrans;
          trans = stpnt->trans;
          tcrp = tcr[i]; 
          for (j=0; j<numtr; j++,trans++,tcrp++) {     
            if (*tcrp > maxval) {
                 maxval = *tcrp;	/* find the largest transition */
		 maxtr = j;
	    }
          }
	  mconc = &chpnt->conc[maxtr];
	  if (mconc->nchan + conc->nchan < 0) {  /* decr the largest trans */
             for (k=0; k<nstate; k++) {
                if (chpnt->conc[k].nchan + conc->nchan > 0) {
                   chpnt->conc[k].nchan += conc->nchan; /* decr first one */
		   break;
		}
	     }
	  }
	  else mconc->nchan += conc->nchan;
	  conc->nchan = 0;
       }
       nconc->dchan = 0;
       nchano += int(conc->nchan*stpnt->cond+0.5);
//       totnchan += conc->nchan;
    }
    if (nchano > chpnt->nchan) nchano = int(chpnt->nchan+0.5);
    chpnt->cno = nchano;
    conduct = ((double)nchano)/chpnt->nchan;
//fprintf (stderr,"totchan %d\n",totnchan);

  } /* noise */

  else {

   conc = chpnt->conc;
   nconc = newconc;
   for (i=0; i<nstate; i++,conc++,nconc++) {  /* zero deltas */
     conc->cval = conc->cest;
     nconc->dcon = 0;
     nconc->icon = 0;
   }
   conc = chpnt->conc;
   nconc = newconc;
   stpnt = chpnt->chtyp->state;
   for (i=0; i<nstate; i++,conc++,nconc++,stpnt++) {  /* zero deltas, set tcr */
     tcrp = tcr[i];				/* pointer to precalc values */
     numtr = stpnt->numtrans;			/*per timinc */
     trans = stpnt->trans;	     		/* get parameters */
     for (j=0; j<numtr; j++,trans++,f_rate++,trmul++,rateo++,tcrp++){ 
          trrate = conc->cval * *tcrp;
          newconc[*trans].dcon += trrate;
          nconc->dcon -= trrate;
          nconc->icon += *tcrp;			/* save for implicit factor */
     }
   }	
	/* after calculating deltas, find first-order estimate: */

   conc = chpnt->conc;
   nconc = newconc;
   stpnt = chpnt->chtyp->state;
   for (i=0; i<nstate; i++,stpnt++,conc++,nconc++) {
       conc->cest = conc->cval + nconc->dcon;
   }

	/* compute second-order estimate: */

   nconc = newconc;
   for (i=0; i<nstate; i++,nconc++) {  	     /* zero deltas, set tcr */
     oldconc[i].dcon = nconc->dcon;  	     /* save first order estimate */
     nconc->dcon = 0;
   }

   critc = 1e-8;
   for (cmaxerr=1.0,errcnt=0; cmaxerr>critc && errcnt<MAXERRCNT; errcnt++) {

     nconc = newconc;
     for (i=0; i<nstate; i++,nconc++) {	/* find change in conc */
       nconc->dcon = 0;			/* erase delta for next iteration */
     }
     cmaxerr = 0.0;
//     totconc = 0.0;
     conc = chpnt->conc;
     nconc = newconc;
     stpnt = chpnt->chtyp->state;
     for (i=0; i<nstate; i++,stpnt++,conc++,nconc++) {	/* find change in conc */
        trans = stpnt->trans;	     		/* get parameters */
        tcrp = tcr[i];				/* pointer to precalc values */
        numtr = stpnt->numtrans;
        for (j=0; j<numtr; j++,trans++,tcrp++){		/*per timinc */
          trrate = conc->cest * *tcrp;

/* if (chpnt->ctype==NA) fprintf
(stderr,
  "i %d j %d t %d trrat %9.4g cest %9.4g v %g tcrp %g\n",
	i,j,t,trrate,conc->cest,chpnt->comp2->v, *tcrp); /* */

         newconc[*trans].dcon += trrate;

        }  /* for (j=0; j<numtr; ...) */
      }   /*  for (i=0; i<nstate; ...) */

	/* after calculating all deltas, find new conc estimates: */

    conc = chpnt->conc;
    nconc = newconc;
    stpnt = chpnt->chtyp->state;
    for (i=0; i<nstate; i++,stpnt++,conc++,nconc++) {

    if (implicit) cest = (conc->cval + nconc->dcon) / (1 + nconc->icon);

    else  cest = (conc->cval + 0.5*(nconc->dcon+oldconc[i].dcon )) / 
			(1 + 0.5 * nconc->icon);
    cerr = conc->cest;
    conc->cest = cest;
    cerr -= cest;
//    totconc += cest;

#ifdef DEBUG
  if (debug & 16 && debugz && 1)
 fprintf(stderr,"2 typ %d conc %10.4g dconc %10.4g err %10.4g cond %g  i %d\n",
		chpnt->ctype,conc->cest,nconc->dcon,cerr,stpnt->cond,i); /* */
#endif
     if (cerr < 0.0) cerr = -cerr; 	/* absolute value */
     if (cerr > cmaxerr) cmaxerr = cerr;

 /*fprintf (stderr,"cmax %10.4g\n",cmaxerr);   /* */

   }  /*  for (i=0; i<nstate; ...) */

/*  fprintf (stderr,"\n");   /* */


#ifdef DEBUG
  if (debug & 16 && debugz && 1)
     fprintf (stderr,"totconc %g\n",totconc);   /* */
#endif
/* fprintf (stderr,"cmax %9.4g\n",cmaxerr);   /* */
  }    /* for (cmaxerr;;) */
/* fprintf (stderr,"count %d\n",errcnt);   /* */

    conduct = 0.0;
    conc = chpnt->conc;
    stpnt = chpnt->chtyp->state;
    for (i=0; i<nstate; i++,stpnt++,conc++) { /* integrate and total cond */
      conduct += conc->cest * stpnt->cond;
  /*    if (chpnt->ctype==NA) fprintf (stderr,"%8.4g ",conc->cest); /* */
    }
  /* if (chpnt->ctype==NA) fprintf (stderr,"\n"); /* */

    chpnt->cno = int(conduct*chpnt->nchan);
  }
  return (conduct);
}

/*------------------------------------*/

#define MAXOLDCHAN 50

static chan *oldchans = (chan *)NULL;

void save_oldchan(chan *cp, double conduct, double critc) 

/* save a computed equilibrium state for a channel */

{
    int i,nstate;
    chan *op,*np;
    static stconc *conc,*oconc;
				/* count up saved states */	
  for (op=oldchans,i=0; op; op=(chan*)op->next,i++) 
	{ }
 				/* check to see if state exists */ 
  if (i < MAXOLDCHAN) {
        np = alloc_chan(cp->ctype,cp->stype);
        np->next = oldchans;
	oldchans = np;
        op = np;

        nstate = cp->numstate;
        conc = cp->conc;
        oconc = op->conc;
        for (i=0; i<nstate; i++,conc++,oconc++) { /* copy the old equilib vals */
           oconc->cest = conc->cest;
        }

        op->ctype   = cp->ctype;
        op->stype   = cp->stype;
        op->voffsm  = cp->voffsm;
        op->voffsh  = cp->voffsh;
        op->arate    = cp->arate;
        op->brate    = cp->brate;
        op->crate    = cp->crate;
        op->drate    = cp->drate;
        op->conduct = conduct;
	op->vrev    = cp->comp1->v;
	op->maxcond = critc;		/* save accuracy crit */
  }
}

/*------------------------------------*/

int find_oldchan(chan *cp, double *conduct, double critc) 

/* find a previously computed equilibrium state for a channel */

{
    int i,found,nstate;
    chan *op;
    static stconc *conc,*oconc;

  found = 0;
  for (op=oldchans; op; op=(chan*)op->next) {/* check to see if state exists */
    if (op->vrev  == cp->comp1->v &&
	op->ctype == cp->ctype &&
        op->stype == cp->stype &&
        op->voffsm == cp->voffsm &&
        op->voffsh == cp->voffsh &&
        op->arate == cp->arate &&
        op->brate == cp->brate &&
        op->crate == cp->crate &&
        op->drate == cp->drate &&
	op->maxcond <= crit)
     { found = 1; break; }
  }
/* fprintf (stderr,"found %d %d %d op %d\n",found,cp->ctype,cp->stype,op); */
  if (found) {
      nstate = cp->numstate;
      conc = cp->conc;
      oconc = op->conc;
      for (i=0; i<nstate; i++,conc++,oconc++) {  /* copy the old equilib vals */
         conc->cest = oconc->cest;
      }
      *conduct = op->conduct;
       return 1;
  }
  else return 0;
}

/*------------------------------------*/

double dochani(chan *chpnt, double critc)

/* Implicit calculation of equilibrium state concentrations */

/*  Try to save time by returning old value when channel and
    voltage is the same as in previous calculation. */

{
   register int i,j;
   int t,nstate,numtr,oldfl,n,nt;
   static chanstate *stpnt;
   static stconc *conc;
   static double (**f_rate)(chan *cpnt), trrate, conduct, *tcrp;
   static float *trmul;
   static char *rateo;
   double trm;			/* multiplier for rate func */
   char *trans;
   double critcx,cerr,cmaxerr,tcomp,totconc;
   static double tcr[NUMSTATE][NUMTRANS]; /* temp store for trrate in calc */
   static dstconc *nconc,newconc[NUMSTATE];

#ifdef DEBUG
  if (debug & 1 && debugz && 16)
	 fprintf (stderr,"dochani type %d %d\n",chpnt->ctype,chpnt->stype);
#endif

  if (find_oldchan(chpnt,&conduct,critc)) {

/*  fprintf (stderr,"dochani: using old values %d\n",chpnt->ctype); /* */
     return conduct;
  }

 nt = 0;
 tcomp = 10.0;
 for (t=0; t<7; t++) {
  tcomp *= 10;
  nstate = chpnt->numstate;
  nconc = newconc;
  conc = chpnt->conc;
  stpnt = chpnt->chtyp->state;
  for (i=0; i<nstate; i++,conc++,nconc++,stpnt++) {  /* zero deltas, set tcr */
    conc->cval = conc->cest;
    nconc->dcon = 0;
    nconc->icon = 0;
    tcrp = tcr[i];                          /* pointer to precalc values */
    f_rate = stpnt->trate;
    trmul = stpnt->ratemul;
    numtr = stpnt->numtrans;
    rateo = stpnt->rateo;
    for (j=0; j<numtr; j++,trans++,f_rate++,trmul++,rateo++,tcrp++){ 
	switch (*rateo) {
	default:
	case 0: trm = chpnt->arate; break; 
	case 1: trm = chpnt->brate; break; 
	case 2: trm = chpnt->crate; break; 
	case 3: trm = chpnt->drate; break; 
        }
        *tcrp = (*f_rate)(chpnt) * *trmul * trm * tcomp;
        nconc->icon += *tcrp;

/*  fprintf (stderr,"i %d j %d tcrp %9.4g rate %9.4g mul %9.4g trm %g\n",
	i,j, *tcrp,(*f_rate)(chpnt),*trmul,trm); /* */

    } 						/* save precalc val */
  }	
	/* calculate delta concentrations from rate functions: */

  critcx = critc / tcomp * 1e6;

/* fprintf (stderr,"critcx %g critc %g tcomp %g\n",critcx,critc,tcomp); /* */

  for (n=0,cmaxerr=1.0; cmaxerr>critcx; n++) {
    cmaxerr = 0.0;
    conc = chpnt->conc;
    stpnt = chpnt->chtyp->state;
    for (i=0; i<nstate; i++,stpnt++,conc++) { /* find change in conc */
       trans = stpnt->trans;	     		/* get parameters */
       tcrp = tcr[i];				/* pointer to precalc values */
       numtr = stpnt->numtrans;
       for (j=0; j<numtr; j++,trans++,tcrp++){	/* per timinc */
         trrate = conc->cest * *tcrp;
         newconc[*trans].dcon += trrate;

/* if (chpnt->ctype==NA)
  fprintf (stderr,"i %d j %d trrate %10.4g cest %10.4g tcrp %g\n",
	i,j,trrate,conc->cest,*tcrp); /* */
       }
    }
	/* after calculating all deltas, find new conc estimates: */

    conc = chpnt->conc;
    nconc = newconc;
    stpnt = chpnt->chtyp->state;
    for (i=0; i<nstate; i++,stpnt++,conc++,nconc++) {
       cerr = -conc->cest;
       conc->cest = (conc->cval + nconc->dcon) / (1 + nconc->icon);
       cerr += conc->cest; 
       conc->cest -= cerr * .017;	/* relaxation factor */
					/* Tune this for fastest init. */ 
/* fprintf			
 (stderr,"init typ %d conc %9.4g dconc %9.4g err %9.4g mxerr %9.4g c %g st %d\n",
	chpnt->ctype,conc->cest,nconc->dcon,cerr,cmaxerr,stpnt->cond,i); /* */

       if (cerr < 0.0) cerr = -cerr; 	/* absolute value */
       if (cerr > cmaxerr) cmaxerr = cerr;
       nconc->dcon = 0;			/* erase delta for next iteration */
    }
/* fprintf (stderr,"\n");   /* */
  }    /* for (cmaxerr;;) */
/* fprintf (stderr,"n %d\n",n);   /* */
  nt += n;

 }  /* for (t;;) */
/* fprintf (stderr,"nt %d\n\n",nt);   /* */

  totconc = 0.0;
  conduct = 0.0;
  conc = chpnt->conc;
  stpnt = chpnt->chtyp->state;
  for (i=0; i<nstate; i++,stpnt++,conc++) { /* integrate and total cond */
    conduct += conc->cest * stpnt->cond;
/*  totconc += conc->cest; /* */
  }

 /* fprintf (stderr,"totconc %g\n",totconc);   /* */

  save_oldchan(chpnt,conduct,critc);

#ifdef DEBUG
  if (debug & 1 && debugz && 16) fprintf (stderr,"dochani end\n");
#endif

  return (conduct);
}

/*------------------------------------*/

void docacomp (cacomp *capnt)
                
/* Calculate calcium concentration for next time step. 
   Concentration is computed for thin shell at inner surface
   of membrane.  Calcium diffuses through ten shells and then
   into the inner volume of the cell.  Calcium flux is converted
   into an equivalent concentration at a shell outside the membrane.
   This allows the diffusion equation to be run identically for
   all the compartments (except the inner core).  

   Based on Yamada, Koch, Adams (1998) and de Schutter & Smolen (1998) 

   Use equation		caflux = D * time * surf/dist * (C1 - C2).

*/

{
  register int i;
  static int j,k,n;
  double *cais,*cab,cacrit;
  int cashell,cashellp,nshell,z;
  comp *pnt;
  conlst *lpnt;
  chan *chpnt;
  hhchan *hhpnt;
  chanparm *chp;
  static double caratio,cadf,gfrac,chdf,chgfrac,chgvrev;
  static double vf,vf2,vfi,vfi2,cvfi,cavfi;
  static double nao, nai, ko, ki, cao;
  static double navfi, kvfi;
  static double c,h,v,vm,conduct,ica,eica,pica,ina;
  static double caflux,eex,vmax,kexo,ca1,cai,tcai,camaxerr,err;
  static double casf0,casfn,casfc,camin,camind,implfn,implfc,implf0;
  static double cabt, cabti, cabf, cabr,ea,eb,ebca,dvm,vm2,ea2,eb2,eica2;
  static int caisiz = 0;		/* number of ca shell estimates */
  static double *ccaise=(double*)NULL;	/* calcium internal shell estimates */
  static double *caiso=(double*)NULL;	/* calcium shell old estimates */
  static double *caisn=(double*)NULL;	/* calcium shell neighbors */
  static double *cabe=(double*)NULL;	/* calcium buffer new estimates */
  static double *caba=(double*)NULL;	/* adder to Ca on rhs */
  static double *cabm=(double*)NULL;	/* multiplier for Best */
  static double *cabea=(double*)NULL;	/* adder for rhs of Best */
  static double *implf=(double*)NULL;	/* implicit factor for whole rhs */
  static double *implfb=(double*)NULL;	/* implicit factor for whole rhs */
  register double *caise;		/* calcium internal shell estimates */
  
#ifdef DEBUG
  if (debug & 1 && debugz && 16) fprintf (stderr,"doca\n");
#endif

  if (!(pnt=capnt->comp1)) {		/* voltage compartment */
     fprintf (stderr,"docacomp: can't find voltage compartment\n");
     return;
  }
  vm = pnt->v - capnt->cavoff;

	/* First, calculate driving force from cao, cai, nao, nai, ko, ki */
        /*  using GHK current eqn.    */
	/*  see Hille 2nd ed, p 341-342. */

	/* Find channel's conductance for computation of Ca channel current */ 
	/*  by calculating the slope of the driving force */
	/*  See de Schutter & Smolen (1998) in "Methods..." p 218 */

	/* Note that this calculation is done for mono- and  divalent ions */

    if (ncabs(vm) < 1e-8) vm=1e-8;
    vf = exp(vm*frt);			/* for monovalent ions */
    vf2 = vf*vf;			/* for divalent ions, exp (vm*f2rt) */
    vfi  = 1.0 / (1.0 - vf);
    vfi2 = 1.0 / (1.0 - vf2);
    cao = capnt->cao;
    caratio = capnt->cais[0] / cao;
    cavfi = (1.0-caratio*vf2)*vfi2;
    if(caratio < 10.0) { 
			/* Calculate Ca current from GHK current equation */
			/*  for calculation of Ica for Ca diffusion below */
			/* See Hille, 1992, p 341 */
			/* de Schutter & Smolen, (1998) for deriv of "gfrac" */
       cadf = -vm * cavfi;
     }
     else {
      cadf = capnt->vrev - vm;		/* driving force */
     } 
		/* Next, calculate ca channel conductances. */
		/* Could be several ca channels in this ca compartment */

  ica = 0.0;
  for (lpnt=capnt->clst; lpnt; lpnt=lpnt->next) {
     chpnt = (chan *)lpnt->conpnt;
     if (! chpnt) continue;
/* fprintf (stderr,"doca: chan %d %d\n",chpnt->ctype,chpnt->stype); /* */
     hhpnt = (hhchan *)chpnt; 
    switch (chpnt->ctype) {

     case CA:					/* Ca L-type */
       switch (chpnt->stype) {
       case 0:
	  chp = &hhpnt->chtyp->parm[0];
      	  chanimpl(chp,hhpnt); /* set Ca rates */
	  c = hhpnt->m * chp->ival[0] + chp->ival[1];
          hhpnt->m = c;
          conduct = c*c*c;
	  if (hhpnt->nchan) {
	    hhpnt->conduct = conduct;
  	    conduct = dochan2(hhpnt);		/* do noise only */
	  }
          break;

       case 2:					/* Ca T-type */
	  chp = &hhpnt->chtyp->parm[0];
      	  chanimpl(chp,hhpnt);
	  c = hhpnt->m * chp->ival[0] + chp->ival[1];
	  chp++;
      	  chanimpl(chp,hhpnt);
	  h = hhpnt->h * chp->ival[0] + chp->ival[1];
	  hhpnt->m = c;
	  hhpnt->h = h;
  	  conduct = c*c*h;
	  if (hhpnt->nchan) {
	    hhpnt->conduct = conduct;
  	    conduct = dochan2(hhpnt);		/* do noise only */
	  }
	  break;
       case 1:
       case 3:
  	  conduct = dochan2(chpnt);
	  break;
       }	  /* switch stype */

       if (chpnt->setvrev) {	/* if chan vrev has been individually set */
         caratio = 1.0/exp(chpnt->vrev * f2rt);
         cvfi = (1.0-caratio*vf2)*vfi2;
         chgfrac = cvfi + f2rt*vm*(1.0-caratio)*vf2*vfi2*vfi2;
       }
       else {
		/* Nernst potential modified by GHK for [K]i */
		/* see Hille 2nd ed, p 107 & p 344           */
		/* and Lewis, 1979 J. Physiol 286:417-665    */

		/* For each Ca channel find its "slope" (chgfrac) and  */
		/*  "slope  reversal potential" (chgvrev) */
		/*  for computation of total channel current */ 
		/*  below in "runcomp()" */

		/* Calculate chgrvrev and chgfrac for this chan only */
		/* Note that this calculation is based on [Ca] and [K] */

	     static chantype *chtyp;
	     static double kratio, naratio, caperm, kperm, naperm;
	     static double ko, nao, kca, naca; 

         nao = chpnt->nao;
	 ko  = chpnt->ko;
	 if (cao==0) {
	   naca = 10.0;
	   kca = 2.0;
	 } 
	 else {
	   naca = nao/cao;
	   kca = ko/cao;
	 }
         naratio = chpnt->naratio;
	 kratio  = chpnt->kratio;
    /*	 navfi = (1.0-naratio*vf)*vfi*naca; */ /* leave out for time being */
	 kvfi  = (1.0-kratio*vf) *vfi*kca;
	 if (!(chtyp=chpnt->chtyp)) {
	    naperm = dpnaca;
     	    kperm  = dpkca;
	 }
	 else {
	    naperm = chpnt->chtyp->perm[PNA];
	    kperm  = chpnt->chtyp->perm[PK];
	 }
         /* cvfi = (1.0-caratio*vf2)*vfi2; */	/* calc based only on [Ca] */
         /* chgfrac = cvfi + f2rt*vm*(1.0-caratio)*vf2*vfi2*vfi2; */

	 caperm = chpnt->caperm;
         cvfi = (caperm*cavfi + 
             /*  naperm*navfi*0.25 + */		/* leave out Na from Ca chans */
	         kperm*kvfi*0.25);		/* for time being */

         chgfrac = cvfi + frt*vm*(2 *   caperm*(1.0-caratio)*vf2*vfi2*vfi2 +
			  /*  naca*0.25*naperm*(1.0-naratio)*vf*vfi*vfi) +  */
				kca*0.25*kperm*(1.0- kratio)*vf*vfi*vfi);  
       }
       if(caratio < 10.0) { 
          chdf = -vm * cvfi;
          chgvrev = chdf/chgfrac + vm;	/* "slope vrev" * frac of g */
       }
       else {
          chgvrev = capnt->vrev;	/* vrev for total channel current */
          chgfrac = 1.0;		/* conductance frac for tot current */
       } 

      conduct = conduct * chpnt->maxcond;
      ica += conduct * cadf;		/* calcium current, used here only */
      chpnt->gvrev = chgvrev;		/* Set "slope" vrev and cond */
      chpnt->conduct = conduct * chgfrac; /*  for use in runcomp() below. */
    break;

      default: 
    case SYN2:
    case AMPA:
    case NMDA:
    case CGMP:				/* retrieve prev. ica estimate */
        conduct = chpnt->conduct * chpnt->caperm;
        ica += conduct * cadf;		/* Ca current, used here only */
	break;
    }	/* switch ctype */

  }

/*  fprintf (stderr,"ica %g\n",ica);   	/* */
  capnt->ica = ica;			/* current through ca-selective chan */

  casf0 = capnt->casf0;			/* flux factor=timinc*F2/(4PI*r*r*dr) */
  cais  = capnt->cais;
  casfn = capnt->casfn;			/* shell factor = 1 / (dr*dr) */
  casfc = capnt->casfc;			/* shell factor for core */

  ca1 = capnt->cais[0];
  camin = capnt->cai;			/* min for pump = starting cai */
  camind = 10e-9;			/* min for diffusion in shells */ 
  tcai = ca1-camind;
  tcai = max (tcai,0);

  if ((kexo=capnt->kexo) > 0) {		/* sodium-calcium exchanger */
					/* From de Schutter & Smolen 1998 */
					/* after DiFrancesco & Noble 1985 */
				        /* Gabbiani,Midtgaard and Knopfel 1994*/

		/* orig calculation */
/*     eica = kexo*exp(capnt->ego*vm) - tcai*capnt->kexi*exp(capnt->egi*vm); */

     ea = kexo*exp(capnt->ego*vm);
     eb = 				     capnt->kexi*exp(capnt->egi*vm);
     eica = ea - capnt->cais[0] * eb;
     ebca = eb * casf0;

     dvm = 1e-6;
     vm2 = vm + dvm;
     ea2 = kexo*exp(capnt->ego*vm2);
     eb2 = 				     capnt->kexi*exp(capnt->egi*vm2);

     if (!implicit) {		/* CN */

       eica2 = ea2 - capnt->cais[0] * eb2;
       capnt->slope = (1-NAEXCHR/CAZ)*(eica2-eica) / dvm;
       capnt->eica = eica*(1-NAEXCHR/CAZ) - vm*capnt->slope;
     }

     ina = -NAEXCHR * eica / CAZ;	/* sodium current */
     ica += ea;				/* additional calcium current */
/*   ica += eica;			/* additional calcium current */
  }
  else eica = ina = ea = eb = 0.0;

  if ((vmax=capnt->vmax) > 0) { 	/* calcium pump */
     pica = -vmax * tcai / (capnt->pkm + tcai);
     ica += pica;			/* pumped calcium */
  }
  else pica = 0.0;

		 /* Compute Ca concentration at inner mem surface. */

/* fprintf (stderr,"ica %g pica %g eica %g ina %g\n",ica,pica,eica,ina); /* */

  capnt->ipump = pica;			/* Ca pump current */
  capnt->iexch = eica + ina;		/* total exchanger current */

					/*       = timestep * M/sec/A/vol(dm) */
  caflux = ica * casf0;   		/* flux in M (moles/liter) for 1st sh */

/* fprintf (stderr,"vm %g ica %g caflux %g eb*cai %g\n",
			vm,ica,caflux,ebca*cais[0]); /* */

  cab   = capnt->cab;
  cashell = capnt->cashell;		/* number of shells */
  cashellp = cashell - 1;		/* loop iteration limit */
  nshell  = cashell - 1;
  if (nshell <= 1) nshell=1;
					/* make temporary space */
					/* to compute diffusion */
  if (cashell > caisiz) {		/* if bigger than previous */
   caisiz = cashell;
   if (ccaise) free(ccaise);		/* free old space */
   if ((ccaise=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
    fprintf(stderr,"docacomp: no space left for %d calcium shells\n",nshell);
    return;
   }
   if (caiso) free(caiso);		/* free old space */
   if ((caiso=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
      fprintf(stderr,"docacomp: no space left for %d calcium shells\n",nshell);
      return;
   }
   if (caisn) free(caisn);
   if ((caisn=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
     fprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
      return;
   }
  }
  if (cab) {			/* temporary calc arrays for Ca buffer */
	static int cabsiz;
    if(cashell > cabsiz) {
     cabsiz = cashell;
     if (cabe) free(cabe);
     if ((cabe=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       fprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (implf) free(implf);
     if ((implf=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       fprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (caba) free(caba);
     if ((caba=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       fprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (cabm) free(cabm);
     if ((cabm=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       fprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (cabea) free(cabea);
     if ((cabea=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       fprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (implfb) free(implfb);
     if ((implfb=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       fprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
    }
  }
 
  caise = ccaise;

  for (i=0; i<cashell; i++) {		/* reset old estimates */
    caise[i] = cais[i];
    caiso[i] = 0.0;
  }

 implf0 = 1.0 + casfn + ebca;
 implfn = 1.0 + 2.0*casfn;
 implfc = 1.0 + casfc;

 if (cab) {				/* set up implicit factors */
   cabt  = capnt->cabt;
   cabti = capnt->cabti;
   cabf  = capnt->cabf;
   cabr  = capnt->cabr;
   for (i=0; i<cashell; i++) {		/* reset old estimates */
     cabe[i] = cab[i];
   }
   if (implicit) {			/* precompute for iteration below */
    for (i=0; i<cashell; i++) {	
     if (i==0) {
       caba[i] = cais[i]+ caflux + cabr*cabti;	/* adder to [Ca] */
       cabea[i] = cab[i] + cabr*cabti;	/* multipler for caest */
     }
     else {
       caba[i] = cais[i] + cabr*cabt;	/* adder to [Ca] */
       cabea[i] = cab[i] + cabr*cabt;	/* multipler for caest */
     }
    }
   }
   else {	/* CN */
    for (i=0; i<cashell; i++) {	
     if (i==0) 
        implf[i] = implf0 + cabf*cab[i]*0.5; 
     else
        implf[i] = implfn + cabf*cab[i]*0.5; 
     implfb[i] = 1.0 + cabr*0.5; /* buf: divider for whole rhs*/
     if (i==0) {
        caba[i] = cais[i]+caflux+(cabti- cab[i]*0.5)*cabr;   /* adder to rhs */
        cabea[i] = (cab[i]*(1.0-(cabr+cabf*cais[i])*0.5) + cabr*cabti);
     }
     else {
        caba[i] = cais[i]+(cabt - cab[i]*0.5)*cabr;
        cabea[i] = (cab[i]*(1.0-(cabr+cabf*cais[i])*0.5) + cabr*cabt);
     }
     cabm[i]  = (cais[i]*cabf + cabr)*0.5;		/* multipler for cabe */
    }
   implf[i-1] = implfc+cabf*cab[i-1]*0.5; /* core: divider for whole rhs */
  } /* CN */
 } /* cab */

if (!implicit) {
  for (i=0; i<cashell; i++) {	
     if (i==0)           caisn[i] = cais[i+1] - cais[i];
     else if (i<cashellp)caisn[i] = cais[i-1] + cais[i+1] - 2.0*cais[i];
     else                caisn[i] = cais[i-1] - cais[i];
  }
}

 cacrit = 1e-10;
 for (camaxerr=1e8,n=0; camaxerr>cacrit; n++) {	/* relaxation */

 if (camaxerr > 1e8 || n>1000) {
    fprintf (stderr,"Cacomp: panic, no convergence, continuing...\n");
    break;
 }
 camaxerr = 0.0;
					/* solve ca concentration */

 if (cab) {				/* active buffer */
   if (implicit) {

     caise[0]=(caba[0]+caise[1]*casfn - cabr*cabe[0]) /
		(implf0 + cabf*cabe[0]); 		/* first shell */
     if (caise[0] < camind)  caise[0] = camind;
     cabe[0] = cabea[0] / (1.0+cabr+cabf*caise[0]);	

    for (i=1; j=i-1,k=i+1,i<cashellp; i=k) {
     caise[i]=(caba[i]+(caise[j]+caise[k])*casfn - cabr*cabe[i]) /
		(implfn + cabf*cabe[i]); /* normal shells */
     if (caise[i] < camind)  caise[i] = camind;
     cabe[i] = cabea[i] / (1.0+cabr+cabf*caise[i]);	
    }
    caise[i]= (caba[i] + caise[j]*casfc - cabr*cabe[i]) /
		 (implfc+cabf*cabe[i]);   		/* core */
    if (caise[i] < camind)  caise[i] = camind;
    cabe[i] = cabea[i] / (1.0+cabr+cabf*caise[i]);	
   }
 else {		/* CN */

      caise[0] = (caba[0]+ (caise[1]+caisn[0]) * casfn -
	cabe[0]*cabm[0]) / implf[0];			/* first shell */
      if (caise[0] < camind)  caise[0] = camind;
      cabe[0] = cabea[0] / (implfb[0] + cabf*caise[0]*0.5);

    for (i=1; j=i-1,k=i+1,i<cashellp; i=k) {
      caise[i] = (caba[i]+ (caise[j]+caise[k]+caisn[i]) * casfn -
	cabe[i]*cabm[i]) / implf[i];
      if (caise[i] < camind)  caise[i] = camind;
      cabe[i] = cabea[i] / (implfb[i] + cabf*caise[i]*0.5);
      }
							/* core */
     caise[i]=(caba[i] + (caise[j]+caisn[i])*casfc - cabe[i]*cabm[i]) /
			implf[i]; 
      if (caise[i] < camind)  caise[i] = camind;
     cabe[i] = cabea[i] / (implfb[i] + cabf*caise[i]*0.5);
    }
 }
 else {					/* passive buffer */
   if (implicit) {
							/* first shell */
     caise[0]=(cais[0]+caflux+(caise[1])*casfn)/implf0; 
     if (caise[0] < camind)  caise[0] = camind;

    for (i=1; j=i-1,k=i+1,i<cashellp; i=k) {	  
     caise[i]=(cais[i]+(caise[j]+caise[k])*casfn)/implfn; /* normal shells */
    }
    caise[i]= (cais[i] + caise[j]*casfc) / implfc;   /* core */
/* fprintf (stderr,"caise0 %g caise1 %g\n",caise[0],caise[1]); */

   }
 else {		/* CN */
//fprintf (stderr,"111 t %g caflux %g caise0 %20.16g caise1 %20.16g caisn0 %20.16g casfn %20.16g implf0 %20.16g\n",
//		xtime,caflux,caise[0],caise[1],caisn[0],casfn,implf0);
     caise[0]= (cais[0]+caflux+(caise[1]+caisn[0])*casfn) / implf0;
     if (caise[0] < camind)  caise[0] = camind;

//fprintf (stderr,"222 t %g caflux %g caise0 %20.16g caise1 %20.16g caisn0 %20.16g casfn %20.16g implf0 %20.16g\n",
//		xtime,caflux,caise[0],caise[1],caisn[0],casfn,implf0);

    for (i=1; j=i-1,k=i+1,i<cashellp; i=k) {
     caise[i]= (cais[i]+
       (caise[j]+caise[k]+caisn[i]) * casfn) / implfn;
     }
					/* core */
      caise[i]= (cais[i] + (caise[j]+caisn[i])*casfc) / implfc;
/* fprintf (stderr,"caise0 %g caise1 %g\n",caise[0],caise[1]); /* */
    }
  }
				/* find error */
  for (i=0; i<cashell; i++) {
    err = caise[i] - caiso[i];
    if (ncabs(err) > camaxerr) camaxerr = err;
//    caise[i] += err * .001;
    caiso[i] = caise[i];		/* save old estimate */
  }

/*  fprintf (stderr,"doca: camaxerr %g\n",camaxerr);  /* */

 }	/* for (camaxerr;;) */

/* fprintf (stderr,"n %d\n",n);  /* */

  for (i=0; i<cashell; i++) {
       cais[i] = caise[i];			/* save estimate */
  }
  if (cab) 
    for (i=0; i<cashell; i++) {
       cab[i] = cabe[i];
    }
		 /* Compute calcium reversal potential */
		 /* Use Nernst equation */

/*  if (xtime > .014) {
    i = i+1;
  } 			/* */

/* compute new exchr current */

 if (kexo > 0) {
     eica  = ea  - capnt->cais[0] * eb;		
     eica2 = ea2 - capnt->cais[0] * eb2;
     if (implicit) {
        capnt->slope = (1-NAEXCHR/CAZ)*(eica2-eica) / dvm;
        capnt->eica  = (1-NAEXCHR/CAZ)*eica - vm*capnt->slope;
     }
     else {
        capnt->slope = (capnt->slope +(1-NAEXCHR/CAZ)*(eica2-eica)/dvm)*0.5;
        capnt->eica  = (capnt->eica  +(1-NAEXCHR/CAZ)*eica-vm*capnt->slope)*0.5;
     }
 }
			/* Nernst potential for Ca. */
			/* Not correct for Ca chans, as they need *
			/* modified GHK voltage equation for [K]i */
			/* see Hille 2nd ed, p 107 & p 344        */

  capnt->vrev = r2ft*log(capnt->cao/capnt->cais[0]); 

#ifdef DEBUG
  if (debug & 1 && debugz && 16) fprintf (stderr,"doca end\n");
#endif

  return;
}

/*------------------------------------*/

void docacompi (cacomp *capnt)
                
/* Initialize calcium concentration in shell for buffers.
   Use implicit mode with long time step.

   Concentration is computed for thin shell at inner surface
   of membrane.  Calcium diffuses through ten shells and then
   into the inner volume of the cell.  Calcium flux is converted
   into an equivalent concentration at a shell outside the membrane.
   This allows the diffusion equation to be run identically for
   all the compartments (except the inner core).  

   Based on Yamada, Koch, Adams (1998) and de Schutter & Smolen (1998) 

   Use equation		caflux = D * time * surf/dist * (C1 - C2).

*/

{
  register int i;
  int j,k,n;
  double *cais,*cab,timeimpl,cacrit;
  int cashell,cashellp,nshell,z;
  comp *pnt;
  conlst *lpnt;
  chan *chpnt;
  static double caflux,camaxerr,err;
  static double casf0,casfn,casfc,camin,camind,implfn,implfc,implf0;
  static double cabt, cabti, cabf, cabr;
  static int caisiz = 0;		/* number of ca shell estimates */
  static double *ccaise=(double*)NULL;	/* calcium internal shell estimates */
  static double *caiso=(double*)NULL;	/* calcium shell old estimates */
  static double *cabe=(double*)NULL;	/* calcium buffer new estimates */
  static double *caba=(double*)NULL;	/* adder to Ca on rhs */
  static double *cabea=(double*)NULL;	/* adder for rhs of Best */
  register double *caise;		/* calcium internal shell estimates */
  
#ifdef DEBUG
  if (debug & 1 && debugz && 16) fprintf (stderr,"doca\n");
#endif

  if (!(pnt=capnt->comp1)) {		/* voltage compartment */
     fprintf (stderr,"docacompi: can't find voltage compartment\n");
     return;
  } 

  caflux = 0;   			/* flux in moles */

/* fprintf (stderr,"ica %g caflux %g\n",ica,caflux);   /* */

  cab   = capnt->cab;
  cashell = capnt->cashell;		/* number of shells */
  cashellp = cashell - 1;		/* loop iteration limit */
  nshell  = cashell - 1;
  if (nshell <= 1) nshell=1;
					/* make temporary space */
					/* to compute diffusion */
  if (cashell > caisiz) {		/* if bigger than previous */
   caisiz = cashell;
   if (ccaise) free(ccaise);		/* free old space */
   if ((ccaise=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
    fprintf(stderr,"docacomp: no space left for %d calcium shells\n",nshell);
    return;
   }
   if (caiso) free(caiso);		/* free old space */
   if ((caiso=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
      fprintf(stderr,"docacomp: no space left for %d calcium shells\n",nshell);
      return;
   }
   if (cab) {			/* temporary calc arrays for Ca buffer */
     if (cabe) free(cabe);
     if ((cabe=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       fprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (caba) free(caba);
     if ((caba=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       fprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (cabea) free(cabea);
     if ((cabea=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       fprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
   }
  }

  timeimpl = 1e6/timinc;		/* timinc = 1e6 seconds */ 
  cais  = capnt->cais;
  casf0 = capnt->casf0 * timeimpl;	/* caflux factor = dt*2F/(4PI*r*r*dr) */
  casfn = capnt->casfn * timeimpl;	/* shell factor = 1 / (dr*dr) */
  casfc = capnt->casfc * timeimpl;	/* shell factor for core */

  caise = ccaise;

  for (i=0; i<cashell; i++) {		/* reset old estimates */
    caise[i] = cais[i];
    caiso[i] = 0.0;
  }

 implf0 = 1.0 + casfn;
 implfn = 1.0 + 2.0*casfn;
 implfc = 1.0 + casfc;

 if (cab) {				/* set up implicit factors */
   cabt  = capnt->cabt;
   cabti = capnt->cabti;
   cabf  = capnt->cabf * timeimpl;
   cabr  = capnt->cabr * timeimpl;
   for (i=0; i<cashell; i++) {		/* reset old estimates */
     cabe[i] = cab[i];
   }
   for (i=0; i<cashell; i++) {	
     if (i==0) {
       caba[i] = cais[i] + cabr*cabti;	/* adder to [Ca] */
       cabea[i] = cab[i] + cabr*cabti;	/* multipler for caest */
     }
     else {
       caba[i] = cais[i] + cabr*cabt;	/* adder to [Ca] */
       cabea[i] = cab[i] + cabr*cabt;	/* multipler for caest */
     }
   }
  }

 cacrit = 1e-10;
 for (camaxerr=1e8,n=0; camaxerr>cacrit; n++) {	/* relaxation */

 if (camaxerr > 1e8 || n>1000) {
    fprintf (stderr,"Cacomp: panic, no convergence, continuing...\n");
    break;
 }
 camaxerr = 0.0;
					/* solve ca concentration */

 if (cab) {				/* active buffer */
     caise[0]=(caba[0]+(caise[1])*casfn - cabr*cabe[0]) /
		(implf0 + cabf*cabe[0]); /* normal shells */
     if (caise[0] < camind)  caise[0] = camind;
     cabe[0] = cabea[0] / (1.0+cabr+cabf*caise[0]);	

    for (i=1; j=i-1,k=i+1,i<cashellp; i=k) {
     caise[i]=(caba[i]+(caise[j]+caise[k])*casfn - cabr*cabe[i]) /
		(implfn + cabf*cabe[i]); /* normal shells */
     if (caise[i] < camind)  caise[i] = camind;
     cabe[i] = cabea[i] / (1.0+cabr+cabf*caise[i]);	
    }
    caise[i]= (caba[i] + caise[j]*casfc - cabr*cabe[i]) /
		 (implfc+cabf*cabe[i]);   /* core */
    if (caise[i] < camind)  caise[i] = camind;
    cabe[i] = cabea[i] / (1.0+cabr+cabf*caise[i]);	
 }
  for (i=0; i<cashell; i++) {
    err = caise[i] - caiso[i];
    if (ncabs(err) > camaxerr) camaxerr = err;
    caise[i] += err * .6;
    caiso[i] = caise[i];		/* save old estimate */
  }

/*  fprintf (stderr,"doca: camaxerr %g\n",camaxerr);  /* */

 }	/* for (camaxerr;;) */

/* fprintf (stderr,"n %d\n",n);  /* */

  for (i=0; i<cashell; i++) {
       cais[i] = caise[i];			/* save estimate */
  }
  if (cab) 
    for (i=0; i<cashell; i++) {
       cab[i] = cabe[i];
    }
			 /* Compute calcium reversal potential */
			 /* Use Nernst equation */

			/* Nernst potential modified by GHK for [K]i */
			/* see Hille 2nd ed, p 107 & p 344           */

  capnt->vrev = r2ft*log(capnt->cao/capnt->cais[0]); 

		/* Compute reversal potential for calcium channels */

  for (lpnt=capnt->clst; lpnt; lpnt=lpnt->next) {
     chpnt = (chan *)lpnt->conpnt;
     if (! chpnt) continue;
     if (chpnt->ctype == CA && !chpnt->setvrev)
		 chpnt->vrev = ghkv(getchantype(chpnt->ctype,chpnt->stype),
					capnt->cao,capnt->cais[0]);
  }

#ifdef DEBUG
  if (debug & 1 && debugz && 16) fprintf (stderr,"doca end\n");
#endif

  return;
}

/*------------------------------------*/

#ifdef DOCHNOISE

double dochnoise(double cconduct, chan *chpnt)

/* Run noise on a channel. Use a two-state machine to calculate
number of events going from closed -> open, and open -> closed.
Number of events is a binomial function of the originating state
and the probability for the transition.  Transition probabilities
are calculated from the tau set by the user (cdur) and the conductance.
*/
 
{
   double cdur,p,chanev;
   double tco, toc, pa, pb;
   int nchan, no, nc, nco, noc;

  if (nchan=chpnt->nchan) {
     if (cconduct > 0 && cconduct < 1.0) {
       cdur = chpnt->cdur / timinc;		/* ch event duration / timinc */
       if (cdur <= 0) cdur = 1e-6;
       no = chpnt->cno;				/* num of prev. open chans */
       nc = nchan - no;				/* num of prev. closed chans */
       tco = cdur / cconduct;			/* tau for close -> open */
       toc = cdur / (1.0-cconduct); /* tau for open -> close */
       pa = 1 - exp (-1.0/tco);			/* rate of closed -> open */
       pb = 1 - exp (-1.0/toc);			/* rate of open -> closed */
       if (chpnt->cstate) rsetstate (chpnt->cstate);
       nco = int(binomdev(pa,nc)); 		/* delta closed to open */
       if (chpnt->cstate) rsetstate (chpnt->cstate);
       noc = int(binomdev(pb,no)); 		/* delta open to closed */
       no += nco - noc;
       chpnt->cno = no;
       cconduct = ((double)no) / nchan;
     }
  }
  return cconduct;
}

#endif

/*------------------------------------*/

void runcomp(void)

/* Advance the network of compartments one step in time,
using a second-order implicit numerical integration.
Iterate the approximation until the error is
reduced below criterion.
*/

{
  register comp *pnt;
  register conlst *lpnt;
  static int left,niter,titer,less,noskip=0;
  static int stype;
  static double nodcur;                    /* total node current */
  static double nodcr;                     /* total unvarying node curr */
  static double ocais;                     /* original cais[0] */
  static conn *conpnt;
  static comp *ocpnt;
  static double delcrit,tdelcrit;
  static double relmult,critx,critxx,critk;
  static double tcond;			   /* total conductance in comp */
  static double tcondn;			   /* tot cond not varying w/vest */
  static double tcondm;			   /* tot membrane cond */
  static double conduct;		   /* conductance of channel */
  static double vest;		   	   /* conductance of channel */
  static recpar *cpt;                      /* receptor constants */
  static hhchan *hhpnt;
  static kcachan *kcapnt;
  static kcachan *kcapnts;
  static chan *chpnt;
  static cacomp *capnt;
  static chanparm *chp;
  static double a,m,n,h,dm,dn,dh;
  static double mmaxerr;		   /* max maxerr for level of critx */
  static double oldmaxerr;		   /* looking for maxerr increases */

#ifdef DEBUG
  if (debug & 1) fprintf (stderr,"runcomp\n");
#endif

 
 if (timinc > 1.0) runsyn(0);			/* if static, calc synapses */

 if (euler) implicit = 0;

 for (pnt=compnt; pnt; pnt=pnt->next)                   /* First estimate */
  {                                                     /* for all comps */
   nodcur = 0.0;		/* current that is also calc in docomp() */
   nodcr = 0.0;			/* current that is precalc for dcomp() */
   tcond = tcondn = 0.0;
   if (capnt=pnt->capnt) {		/* Compute calcium concentration. */
     ocais = capnt->cais[0];
     docacomp(capnt);			/* Must do first because other chans */
     nodcr += capnt->ipump;		/* pump current */
     nodcur+= capnt->iexch;		/* exchanger current */
   }					/*  may be dependent on ca. */

   for (lpnt=pnt->clst; lpnt; lpnt=lpnt->next) 
    {                                           /* check all connections */
     conpnt = lpnt->conpnt;
     if (! conpnt) break;
     if (conpnt->comp1 == pnt) ocpnt = conpnt->comp2; /* get other compartment*/
     else                      ocpnt = conpnt->comp1;

     switch (conpnt->ctype) { 	/* if maxcond is zero, don't compute */

        case GJ:      if (((gj*)conpnt)->maxcond == 0) continue; break;
        case SYNAPSE: if (((synap*)conpnt)->maxcond == 0) continue; break;
        case ROD:
        case CONE:    if (((recep*)conpnt)->maxcond == 0) continue; break;
        case SYN2:
        case AMPA:
        case CGMP:
        case GABA:
        case NMDA:
        case CA:
        case NA:
        case K:       if (((chan*)conpnt)->maxcond == 0) continue; break;
       default:   break;
     }

     switch (conpnt->ctype)
      {                         /* add connection conductances which change */

       case GJ:
                nodcur += conpnt->conduct * ocpnt->v;   /* current thru g j */
                tcond += conpnt->conduct;              /* total cond */
               break;
       case SYNAPSE:					 /* cur thru syn chan */
		if (!((synap*)conpnt)->resp) {		/* regular synapse */
		if (pnt==((synap *)conpnt)->comp2) {
                  nodcr += ((synap *)conpnt)->conduct *((synap *)conpnt)->vrev;
                  tcondn += ((synap *)conpnt)->conduct;   /* total cond */
  	/* fprintf (stderr,"vrev %10.3g conduct %10.3g\n",
		((synap *)conpnt)->vrev,((synap *)conpnt)->conduct); /* */
		}
	       } break;
      case ROD:
      case CONE:                                /* find receptor instantan */
                                                /* conductance */
          cpt = ((recep *)conpnt)->chtyp;
          nodcr += ((recep *)conpnt)->conduct * cpt->vrev;/* cur thru rec chan*/
          tcondn  += ((recep *)conpnt)->conduct;               /* total cond */
/*   fprintf (stderr,"vrev %10.3g conduct %10.3g\n",
		cpt->vrev,((recep *)conpnt)->conduct); /* */
         break;

      case LOAD:
                nodcr += ((load *)conpnt)->conduct * 
			  ((load *)conpnt)->vrev;     /* cur to "ground" */
                break;              /* conductance is collected in pnt->tcond */

      case AXIALRES:
      case RESISTOR:
                nodcur += conpnt->conduct * ocpnt->v;  /* cur to next comp */
                break;              /* conductance is collected in pnt->tcond */
      case BATT:
     		if (ocpnt == conpnt->comp1)
                  nodcur += (ocpnt->v - conpnt->conduct) * BATTCOND;
		else 
                  nodcur += (ocpnt->v + conpnt->conduct) * BATTCOND;
                break;		/* cur to next comp through batt */

      case BUF:	if (pnt == ((dbuf *)conpnt)->comp2) {
		   if (((dbuf *)conpnt)->delay) {
		     pnt->extv = *((dbuf *)conpnt)->delpnt * DELCONV;
		   }
		   else pnt->extv = ((dbuf *)conpnt)->comp1->v;
		  }
		break;

      case CAP: nodcr += (pnt->oldv - ocpnt->oldv) * conpnt->conduct/timinc;
                break;			/* -(v3 - v2) * c / t */

      case CGMP:		/* Channels with Ca current */
      case NMDA:
      case AMPA:		/* Chans without Ca current */
      case GABA:
      case SYN2:
		chpnt = (chan *)conpnt;
  		conduct = dochan2(chpnt);
  		conduct *= chpnt->maxcond;
		chpnt->conduct = conduct;
                nodcr += conduct * chpnt->vrev;
                tcondn += conduct;    			/* total cond */
		break;

      case NA:  hhpnt = (hhchan *)conpnt;
		stype = hhpnt->stype;
		switch (stype) {

		case 0:				/* Hodgkin-Huxley */

		chp = &hhpnt->chtyp->parm[0];
      		chanimpl(chp,hhpnt);
		m = hhpnt->m * chp->ival[0] + chp->ival[1];
		chp++;
      		chanimpl(chp,hhpnt);
		h = hhpnt->h * chp->ival[0] + chp->ival[1];

		hhpnt->m = m;
		hhpnt->h = h;
  		conduct = m*m*m*h;
		if (hhpnt->nchan) {
		  hhpnt->conduct = conduct;
  		  conduct = dochan2(hhpnt);		/* do noise only */
		}
  		conduct *= hhpnt->maxcond;
		hhpnt->conduct = conduct;
                nodcr += conduct * hhpnt->vrev;
                tcondn += conduct;    			/* total cond */
		break;

		case 1:				/* Sequential state */
		case 2:
		case 3:
		case 4:

	 	chpnt = (chan *)conpnt;
  		conduct = dochan2(chpnt);
  		conduct *= chpnt->maxcond;
		chpnt->conduct = conduct;
                nodcr += conduct * chpnt->vrev;
                tcondn += conduct;    			/* total cond */
		break;

		}	/* switch (stype) */
                break;

      case K:   
	 	hhpnt = (hhchan *)conpnt;
		stype = hhpnt->stype;
		switch (stype) {

		case 0:				/* Hodgkin-Huxley */

		chp = &hhpnt->chtyp->parm[0];
      		chanimpl(chp,hhpnt);
		n = hhpnt->m * chp->ival[0] + chp->ival[1];
		hhpnt->m = n;
  		conduct = (n*n)*(n*n);
		if (hhpnt->nchan) {
		  hhpnt->conduct = conduct;
  		  conduct = dochan2(hhpnt);
		}
  		conduct *= hhpnt->maxcond;
		hhpnt->conduct = conduct;
                nodcr += conduct * hhpnt->vrev;
                tcondn += conduct;    			/* total cond */
                break;

		case 1:				/* Sequential state */
		case 4:
		case 5:
	 	chpnt = (chan *)conpnt;
  		conduct = dochan2(chpnt);
  		conduct *= chpnt->maxcond;
		chpnt->conduct = conduct;
                nodcr += conduct * chpnt->vrev;
                tcondn += conduct;    		/* total cond */
		break;

		case 2:				/* Calcium-activated K chan */
	 	kcapnt = (kcachan *)conpnt;
		chp = &kcapnt->chtyp->parm[0];
		if (pnt->capnt) {  /* calcium compartment exists */
		  if (!kcapnt->initfl){	/* Calcium not set up before K chan. */
			double ca;	/* Find equilibrium value for "n". */ 
			double alph,bet;
		    ca = pnt->capnt->cais[0];
		    alph = akcacalc((pnt->v-kcapnt->voffsm),ca,
			kcapnt->arate,kcapnt->d1,kcapnt->k1);
		    bet = bkcacalc((pnt->v-kcapnt->voffsm),ca,
			kcapnt->arate,kcapnt->d2,kcapnt->k2);
		    kcapnt->m = alph / (alph + bet);
		    kcapnt->initfl = 1;
		    ocais = ca;
		  }
 		  kcatab((pnt->v-kcapnt->voffsm),    /* set K rate chtyp */ 
		     (ocais+pnt->capnt->cais[0])*.5,stype,kcapnt->arate,
		     kcapnt->d1,kcapnt->d2,kcapnt->k1,kcapnt->k2,chp);
		}
		else {		/* Can't find ca comp for kca. */
				/* Use default cai */
     		  kcatab((pnt->v-kcapnt->voffsm),    /* set K rate chtyp */ 
		     dcai,stype,kcapnt->arate,
		     kcapnt->d1,kcapnt->d2,kcapnt->k1,kcapnt->k2,chp);
		}
		n = kcapnt->m * chp->ival[0] + chp->ival[1];
		kcapnt->m = n;
  		conduct = n;
		if (kcapnt->nchan) {
		  kcapnt->conduct = conduct;
  		  conduct = dochan2(kcapnt);		/* do noise */
		}
  		conduct *= kcapnt->maxcond;
		kcapnt->conduct = conduct;
                nodcr += conduct * kcapnt->vrev;
                tcondn += conduct;    			/* total cond */
		break;

		case 6:
	        kcapnts = (kcachan *)conpnt;
		if (pnt->capnt) {  /* calcium compartment exists */
		  if (!kcapnts->initfl){/* Calcium not set up before K chan. */
		     dochani((chan *)kcapnts,1e-9);
		     kcapnts->initfl = 1;
		  }
		}
  		conduct = dochan2(kcapnts);
  		conduct *= kcapnts->maxcond;
		kcapnts->conduct = conduct;
                nodcr += conduct * kcapnts->vrev;
                tcondn += conduct;    			/* total cond */
		break;

		case 3:				/* Type A channel */
						/* fast inactivating */
		chp = &hhpnt->chtyp->parm[0];
      		chanimpl(chp,hhpnt);
		a = hhpnt->m * chp->ival[0] + chp->ival[1];
		chp++;
      		chanimpl(chp,hhpnt);
		h = hhpnt->h * chp->ival[0] + chp->ival[1];

		hhpnt->m = a;
		hhpnt->h = h;
  		conduct = a*a*a*h;
		if (hhpnt->nchan) {
		  hhpnt->conduct = conduct;
  		  conduct = dochan2(hhpnt);		/* do noise */
		}
  		conduct *= hhpnt->maxcond;
		hhpnt->conduct = conduct;
                nodcr += conduct * hhpnt->vrev;
                tcondn += conduct;    			/* total cond */
                break;

		}	/* switch (stype) */
		break;

       case CA:   
	 	chpnt = (chan *)conpnt;
		stype = chpnt->stype;
		switch (stype) {	/* calculated by docacomp() above */

		case 0:			/* Hines Ca type c3 "L-type" */
		case 1:			/* Markov version of type 0 */
		case 2:			/* G. Matthews Ca c3h "T-type" */
		case 3:			/* Markov version of type 2 */
		case 4:
		conduct = chpnt->conduct; 
                nodcr += conduct * chpnt->gvrev;
                tcondn += conduct;	/* total cond */
/* fprintf (stderr,"conduct %g ipump %g\n",conduct,capnt->ipump); /* */
		break;

	    default:
	 	chpnt = (chan *)conpnt;
		conduct = chpnt->conduct; 
                nodcr += conduct * chpnt->vrev;
                tcondn += conduct;    			/* total cond */
		break;

		}     /* switch (stype) */
		break;

      default:  break;

      } /* switch */
    }      /* for (lpnt=; ;) */

   tcondm = tcond + tcondn + pnt->tcondm;      /* total membrane conductance */
   tcond  += tcondn + pnt->tcond;              /* total conductance */
   nodcur -= pnt->v * tcond;                   /* current due to comp volts */
   nodcr  += pnt->vrev * pnt->rm;              /* membrane conductance */

   if (pnt->jnoise) {			       /* Johnson noise */
  	if (pnt->jstate) rsetstate (pnt->jstate); /* if set individually */
	nodcr += gasdev() * ktrb * sqrt(tcondm) * pnt->jnoise;
                        /* See Hille 1992 p 325 on Johnson noise */
   }


/*  fprintf (stderr,"pnt->vrev %10.3g\n",pnt->vrev); /* */
/*  fprintf (stderr,"nodcr %10.3g nodcur %10.3g pntcond %g tcond %g\n",
			nodcr,nodcur,pnt->tcond,tcond); /* */

   if (pnt->miscfl & IEXT)                    /* comp has a current src */
      nodcr += pnt->exti;                     /* external current */
   nodcur += nodcr;			      /* total comp current */

 if (euler)	{     			      /* forward-Euler */
     vest = pnt->v + nodcur * pnt->k;	      /* simple first-order est */
     if (pnt->miscfl & VEXT) {                /* comp is a battery */
       pnt->extvi = (vest - pnt->extv) / (pnt->k);
       vest = pnt->extv;		      /*  battery voltage  */
     }
     else {
       pnt->extvi = 0.0;
     }
     pnt->vest = vest;				/* just save the estimate */
 }						/* (not quite) end of Euler */
						/*   (see below) */

 else {						/* Implicit modes */

   pnt->tcondn = tcondn;			/* cond not varying w/vest */
   pnt->implf = 1. / (1. + tcond * pnt->k);     /* implicit factor   */

   if (implicit) {			/* Backward Euler, fully implicit */
     pnt->nodc = nodcr;				/* curr for estimate */
     nodcur *=  pnt->k;                      	/* first-ord only */
     pnt->vest = pnt->v + nodcur;               /* simple first-order est */
     if (pnt->miscfl & VEXT) {                   /* comp is a battery */
          pnt->extvi = (pnt->vest - pnt->extv) / (pnt->k);
          pnt->vest = pnt->extv;
       }
     else {
       pnt->extvi = 0.0;
     }
     pnt->verr = 1.;
   }						/* end of implicit */

   else	{     /* Crank-Nicolson */
     pnt->nodc = nodcur + nodcr;	      /* curr for estimate */
     nodcur *=  pnt->k * 2.;                  /* first-ord only */
     pnt->vest = pnt->v + nodcur;             /* simple first-order est */
     if (pnt->miscfl & VEXT) {                /* comp is a battery */
       pnt->extvi = (pnt->vest - pnt->extv) / (pnt->k*2.);
       pnt->vest = pnt->extv;		      /*  battery voltage  */
     }
     else {
       pnt->extvi = 0.0;
     }
     pnt->verr = 1.;
   }    /* crank-nicholson */
  }  /* CN or implicit */
 }   /* for (pnt=;;) */

 if (euler) {					/* after everything else, */
   for (pnt=compnt; pnt; pnt=pnt->next) {
	pnt->v = pnt->vest;			/* save as new voltage */
   }
   return;
 }

 relmult = relincr * timinc / 1e-4;
 mmaxerr=oldmaxerr=maxerr = 1.0;                 /* set maxerr so loop runs */ 
 ncomp = tcomp = 0; 
 delcrit =   SQRT10; /* */
 left = 1;
 noskip = 0;

for (titer=0,critx=1.000001e-2;
		 mmaxerr >= crit;
		 critx = min(mmaxerr,critx),
		  critx /= delcrit, maxerr = 1.0, titer+=niter, left++) /* */
 {
 critxx = critx;
 oldmaxerr=1.0;
 mmaxerr = 0.0;

 for (niter=0; maxerr>critx; niter++)	/* Iterate until convergence */
  {
   if (timinc > 1.0) runsyn(0);		/* if static, calc synapses */

/* fprintf (stderr,"maxerr %g critx %g critxx %g\n",maxerr,critx,critxx); /* */


/*   if (mmaxerr > critxx*delcrit) critxx = mmaxerr/delcrit; /* */

   maxerr = 0.0;
   if (left&1)
     for (ncomp=0,pnt=compend; pnt; pnt=pnt->last)     /* left estimate */
      {
#ifdef DEBUG
        if (!(debug&NOSKIP))
#endif
        if (!noskip)
          if (pnt->verr<critxx) continue; /* */
        docomp(pnt);
      }
   else
     for (ncomp=0,pnt=compnt; pnt; pnt=pnt->next)      /* right estimate */
      {
#ifdef DEBUG
        if (!(debug&NOSKIP))
#endif
        if (!noskip)
          if (pnt->verr<critxx) continue; /* */
        docomp(pnt);

      }           /* for (ncomp=0,pnt=compnt;) */


   /* Reduce the criterion if error increases, and also at least
       once every 10 iterations.  This causes compartments
       outside the region of iteration (i.e. otherwise skipped over)
       to be recomputed and reduces the total number of iterations
       necessary for complete convergence.
   */
       
  if ((maxerr>oldmaxerr && oldmaxerr!=0) || ((niter>0) && (niter%10==0))) {
    tdelcrit = sqrt(delcrit);
    if (tdelcrit>1.05) {
      delcrit=tdelcrit;			/* */
#ifdef DEBUG
      if (debug & 4 && debugz & 1) {
	 if (maxerr>oldmaxerr)
	 fprintf (stderr,"Maxerr incr:   reducing delcrit to %8.4g\n",delcrit);
	 else
	 fprintf (stderr,"10 iterations: reducing delcrit to %8.4g\n",delcrit);
      }
#endif
    } 
      oldmaxerr = 1.0; 
      noskip = 1; 
    critxx /= delcrit;

#ifdef DEBUG
  if (debug & 4 && debugz & 1)
	 fprintf (stderr,"maxerr %8.3g, reducing critxx to %8.3g\n",
					maxerr,critxx);
#endif

   }  /* if (maxerr>oldmaxerr || ) */

   else {
      oldmaxerr = maxerr; 
   }

   if (noskip) noskip--;
   if (maxerr==0) noskip = 2; 			/* don't skip when unstable */
   if (critxx < critx*1e-2) noskip = 10; 	/* don't skip when unstable */

   if (mmaxerr<maxerr) mmaxerr = maxerr;   /* remember largest maxerr */

#ifdef DEBUG
  if (debug & 4 && debugz & 1)
   fprintf 
    (stderr,"ncomp %4d mmaxerr %-8.3g maxerr %-8.3g critxx %8.3g\n",
					ncomp,mmaxerr,maxerr,critxx);
#endif

   if (critxx < 1e-30) {
      fprintf (stderr,"Ncomp: panic, no convergence, continuing...");
      return;
   }
 
 }            /* for (niter=0; maxerr>critx; niter++)  */


/* Reset error after each level of critx */
/*  This has the effect of spreading the region of iteration,
    and causes compartments just outside to be recomputed.
*/

 critk = .95 * critx / delcrit;
 for (pnt=compnt; pnt; pnt=pnt->next)
  {
   pnt->verr += critk;
  }			/* */


#ifdef DEBUG
  if (debug & 2 && debugz & 1)
   fprintf 
 (stderr,"niter %2d   mmaxerr %-8.3g maxerr %-8.3g critx  %8.3g critk %7.3g\n",
				niter,mmaxerr,maxerr,critx,critk);
#endif

  }         /* for (critx ; ; ) */

#ifdef DEBUG
  if (debug & 2 && debugz & 1)
	fprintf (stderr,"*tcomp %4d titer %d\n",tcomp,titer);
#endif

/* Just save it and calculate relax for next time step. */

 for (pnt=compnt; pnt; pnt=pnt->next)
  {
   pnt->oldv = pnt->v = pnt->vest;		/* save it */
  }

/* This code runs only once per time step, */
/* so it adds negligible time to computation */
/* if model is large and requires many iterations */
/* per time step. */

if (relincr && titer > 10)			/* if lots of iterations */
 for (pnt=compnt; pnt; pnt=pnt->next) {
   if (!(pnt->miscfl & VEXT)) {			/* if comp is not a battery */
      less = pnt->t - pnt->g;			/* then compute oscillations */
      m = ((pnt->g < less) ? less : pnt->g);	/* max of less, greater */
      if (m==pnt->t) pnt->relax += relmult;	/* adjust comp relax */
      else           pnt->relax -= (1.0 - ((double) m)/pnt->t) * 10.0 * relmult;
      if      (pnt->relax > 1.0) pnt->relax = 1.0;
      else if (pnt->relax < 0.0) pnt->relax = 0.0;
    }

#ifdef DEBUG
  if (debug & 8 && debugz & 1)
     fprintf (stderr,"c %-4d relax %-8.4g g %-4d t %-4d misc %d\n",
			   pnt->num,pnt->relax,pnt->g,pnt->t,pnt->miscfl); /* */
#endif
   pnt->g = pnt->t = 0;				/* reset greater and total */
  }

#ifdef DEBUG
  if (debug & 1) fprintf (stderr,"runcomp end.\n");
#endif
}

/*------------------------------------*/

void docomp(register comp *pnt)
{
  register conlst *lpnt;
  register conn *conpnt;
  static cacomp *capnt;                 /* pntr to calcium compartment */
  static double nodcur;                 /* total comp current */
  static double tcondv;                 /* total comp conductance */
  static double conduct;                /* channel conductance */
  static double implf;			/* new implicit factor */
  static comp *ocpnt;
  static double err;
  static double icomp;
  static hhchan *hhpnt;
  static double m,n,h,dm,dn,dh;
  int i;

#ifdef DEBUG
     tcomp++;
     ncomp++;
#endif
     nodcur = pnt->nodc;                        /* start with old node cur */
     tcondv = 0.0;				/* non-varying conductance */
     if (capnt=pnt->capnt) {			/* exchanger current. */
	 static double eica;
       if (capnt->kexo) {
	     double eica;
	 if (implicit) eica = capnt->slope*pnt->v + capnt->eica;
 	 else          eica = capnt->slope*(pnt->v+pnt->vest)*0.5+capnt->eica;
         nodcur += eica;
       }
     }					/*  may be dependent on ca. */
     for (lpnt=pnt->clst; lpnt; lpnt=lpnt->next) 
      {
       conpnt = lpnt->conpnt;                   /* check connection */
       if (! conpnt) break;
       if (conpnt->comp1 == pnt) ocpnt = conpnt->comp2; /* get other comp */
       else                      ocpnt = conpnt->comp1;
       switch (conpnt->ctype)
        {
         case AXIALRES:
         case RESISTOR:
         case GJ:
                nodcur += ocpnt->vest * conpnt->conduct;  /* curr thru resist */
                break;
         case BATT:
     		if (ocpnt == conpnt->comp1) 
                  nodcur += (ocpnt->vest - conpnt->conduct) * BATTCOND;
		else
                  nodcur += (ocpnt->vest + conpnt->conduct) * BATTCOND;
                break;		     /* curr thru batt small resistor */

         case BUF: if (pnt == ((dbuf *)conpnt)->comp2) {
		     if (((dbuf *)conpnt)->delay) {
		       pnt->extv = *((dbuf *)conpnt)->delpnt * DELCONV;
		     }
		     else pnt->extv = ((dbuf *)conpnt)->comp1->vest;
		   }
		break;
         case CAP:
		nodcur += ocpnt->vest * conpnt->conduct / timinc;	
                break; 
         case NA:
      	 case K:   
         case CA:
		break;

         case SYNAPSE:                            /* currents don't need est */
         case ROD:
         case CONE:
         case LOAD:
         default:
                break;
        }  /* switch */
       }  /* for (lpnt= ; ; ) */

     err = pnt->vest;
     if (tcondv) {				/* if any new conductances */
        implf = 1. / (1. + (tcondv+pnt->tcondn+pnt->tcond) * pnt->k);
        pnt->vest = (pnt->v + nodcur * pnt->k) * implf;
     }
     else  {
 	implf = pnt->implf;
        pnt->vest = (pnt->v + nodcur * pnt->k) * implf;
     }
     if (pnt->miscfl & VEXT) {			/* comp is a battery */
	icomp = (pnt->vest - pnt->extv) / (pnt->k * implf);
	if (implicit) {
           pnt->extvi = icomp;
        }
	else {    			/* crank-nicholson */
           pnt->extvi = icomp * 0.5;
	}   
       pnt->vest = pnt->extv;
       pnt->verr = 1.0;			/* set verr so vclamp always runs */
     }
     else {
        err -= pnt->vest;                         /* find convergence deriv */
        pnt->vest -= err * pnt->relax;
        pnt->t++;
        if (err>=0.0) pnt->g++;
	else err = -err;
        pnt->verr = err;
        if (err > maxerr) maxerr = err;
     }
}

/*------------------------------------*/

