/* Segment ncomp in Program nc */

/* Simulates neuronal circuits */
/*  by numerically integrating difference equations */
/*  using iterative relaxation technique */
/*  to implement the "implicit" integration method */ 

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <math.h>
#include <unistd.h>

#ifdef __cplusplus
}
#endif

#include "ndef.h"
#include "nc.h"
#include "y.tab.h"
#include "control.h"

#include "ncsub.h"
#include "ncomp.h"
#include "ncelem.h" 
#include "ncio.h"
#include "nconst.h"

#define DEBUG			/* */
#ifdef DEBUG
#include "ncdebug.h"
#endif

#define NOSKIP 256

#define SQRT10 3.16227766

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

double akcacalc(double v, double ca, double rate, double d1, double k1);
double bkcacalc(double v, double ca, double rate, double d2, double k2);

void docomp(comp *p);
void runsyn(int t);
void chanimpl(chanparm *chp, chan* cpnt);	
void kcatab(double voffs, double cais, int stype, double taua,
	     double d1, double d2, double k1, double k2, chanparm *chp);
void drand_setstate (char *state);
int binomdev(double pp, int n);
chan *alloc_chan(short ctype, short stype);
char *emalloc(unsigned int n);
double drand(void);
double gasdev(void);
chantype *getchantype(int ctype, int cnum);
double ghkv (chantype *chtyp,double cao, double cai);
void restorstate(void);
void efree(void *ptr);
double hillcalc(double ligand, double kd, double p);
double ccavoff(chan *ch);
double compcavoff(comp *pnt);
double vext(conn *chpnt);
double getnt(comp *cpnt, int nttype);
double setnt(comp *cpnt, int nttype, double val);

double get_gvrev(chan *ch, double v);			// in chanfuncs.cc
double get_gfrac(chan *ch, double v);
double get_ioni(chan *ch, int ion);
double get_iono(chan *ch, int ion);
double get_ionp(chan *ch, int ion);
iontab *init_ions(iontab *itab);

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
   static double (**f_rate)(chan *cpnt), trrate, implfh;
   static double conduct,critc,totconc,*tcrp,icon,dcon;
   static double *trmul;
   static char *rateo;
   double trm;			/* multiplier for rate func */
   char *trans;
   double cerr,cmaxerr,cest;
   static dstconc oldconc[NUMSTATE];
   static dstconc *nconc,newconc[NUMSTATE];
   static double tcr[NUMSTATE][NUMTRANS]; /* temp store for trrate in calc */

#define MAXERRCNT 100
#define SINGTHRESH 50
#define MINIMPLFH 0.58
//#define MINIMPLFH 1.0

/* Rationale: implfh sets the mode of integration. 
    implfh = 0.5 -> Crank-Nicolson (second order)
    implfh = 1.0 -> full implicit   (first order)
    implfh = 0.58 -> more stable than C-N
    implfk = overall control of implicit mode, set in "runcomp()" below
*/
  if ((implfh=implfk) < MINIMPLFH) implfh = MINIMPLFH;

  t = chpnt->stype;
  nstate = chpnt->numstate;

   /* Set up rate constant multiplier */

  stpnt = chpnt->chtyp->state;
  for (i=0; i<nstate; i++,stpnt++) {  		/* zero deltas, set tcr */
    tcrp = tcr[i];				/* pointer to precalc values */
    f_rate = stpnt->trate;
    trmul = stpnt->ratemul;
    rateo = stpnt->rateo;
    numtr = stpnt->numtrans;			/* per timinc */
    for (j=0; j<numtr; j++,f_rate++,trmul++,rateo++,tcrp++){ 
	switch (*rateo) {
	default:
	case 0: trm = chpnt->arate; break; 
	case 1: trm = chpnt->brate; break; 
	case 2: trm = chpnt->crate; break; 
	case 3: trm = chpnt->drate; break; 
	case 4: trm = chpnt->erate; break; 
	case 5: trm = chpnt->frate; break; 
        }
        *tcrp = (*f_rate)(chpnt) * *trmul * trm; /* save precalc val */
    }
  }

  if (chpnt->nchan > 0) {			/* Noise on */

        int k,r,dchan,delch,nchano,totnchan,multp;
	double maxp, divp, totp;

		/* When noise is specified, compute number of channels */
		/*  in the states. */

		/* Total up all exit rates from a state and find */
		/*  the maximum exit rate for all the states.*/
		/*  Then divide each rate by the maximum rate */ 
		/*  to make sure that the maximum rate (probability) */
		/*  is less than 1 so it doesn't get truncated in */
		/*  binomdev(p,n). The transitions are iterated to keep the */
		/* total rate constant. */

   stpnt = chpnt->chtyp->state;
   for (maxp=i=0; i<nstate; i++,stpnt++) {
     numtr = stpnt->numtrans;
     tcrp = tcr[i];				/* pointer to precalc rates */
     totp = 0.0;
     for (j=0; j<numtr; j++,tcrp++){ 
       totp += *tcrp;				/* total rate from state */
     }
     if (maxp < totp) maxp = totp;		/* find max state rate */
   }
   multp = int(maxp*1.9 + 1.0);			/* max acceptable p = 0.5 */
   divp = 1.0 / multp;
   /* ncfprintf (stderr,"time %g multp %d\n",time,multp); /* */
 
   if (chpnt->cstate) drand_setstate (chpnt->cstate);
   for (r=0; r<multp; r++) {			/* repeat if p reduced */
    conc = chpnt->conc;
    stpnt = chpnt->chtyp->state;
    nconc = newconc;

    if (chpnt->nchan < SINGTHRESH) {		/* single channel method */

	/* If the number of channels is small, it is faster to  */
	/*  compute 1 random number for each channel to specify */
	/*  which transition to make */

	double stay_trans;
       	double cum_trans;
       	double cum_tr[NUMTRANS];
       	double p;
	int k;

	/* first, find cumulative transition probability */

     for (i=0; i<nstate; i++,stpnt++,conc++,nconc++){/* find change in conc */
      numtr = stpnt->numtrans;
      tcrp = tcr[i];                            /* pointer to precalc values */
      cum_trans = 0;
      for (j=0; j<numtr; j++,tcrp++) {		/* sum total transitions */
        cum_trans += *tcrp*divp; 
        cum_tr[j] = cum_trans;
      }
      /* if (cum_trans > 0.5) {
	stay_trans = 1.0 - exp (-cum_trans);	
      }						/* not nec unless total p>1 */

      dchan = 0;
      for (k=0; k<conc->nchan; k++) {		/* for each chan in state */
        p = drand();				/* prob for transition */
        trans = stpnt->trans;
        for (j=0; j<numtr; j++,trans++) {	/* sum total transitions */
	  if (p < cum_tr[j]) {			/* find transition */
	    newconc[*trans].dchan++;
	    dchan++;
	    break;
	  }
        }
      }
      nconc->dchan -= dchan;
     }  /* for (;i<nstate;) */

    }  /* single channel method */

    else {
     for (i=0; i<nstate; i++,stpnt++,conc++,nconc++){/* find change in conc */
      numtr = stpnt->numtrans;
      trans = stpnt->trans;                     /* get parameters */
      tcrp = tcr[i];                            /* pointer to precalc values */
      dchan = 0;
      for (j=0; j<numtr; j++,trans++,tcrp++){     
        delch = binomdev(*tcrp*divp,conc->nchan); 
/*  ncfprintf (stderr,"multp %d tcrp %g nchan %d del %d state %d tr %d\n",
		multp,*tcrp*divp,conc->nchan,delch,i,j); /**/
        newconc[*trans].dchan += delch;
        dchan += delch;
//        chpnt->conc[*trans].nchan += delch;	/* avoids making nchan neg */
//        conc->nchan -= delch;
      }
      nconc->dchan -= dchan;
     }
    }
    conc = chpnt->conc;
    nconc = newconc;
    stpnt = chpnt->chtyp->state;
    nchano = 0;
    totnchan = 0;
    for (i=0; i<nstate; i++,stpnt++,conc++,nconc++) {     /* Total */
	conc->nchan += nconc->dchan;
	if (conc->nchan < 0) {           /* if nchan is negative */
	     int maxtr;
	     double maxval;
	  maxtr = 0;
	  maxval = 0.0;
	  numtr = stpnt->numtrans;
	  trans = stpnt->trans;
	  tcrp = tcr[i];
	  for (j=0; j<numtr; j++,trans++,tcrp++) {
	    if (*tcrp > maxval) {
		 maxval = *tcrp;        /* find the largest transition */
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
    }  /* for (i;;) */
   }  /* for (r;;) */
   restorstate();		/* restore orig random sequence */

    if (nchano > chpnt->nchan) nchano = int(chpnt->nchan+0.5);
    chpnt->cno = nchano;
    conduct = ((double)nchano)/chpnt->nchan;
//ncfprintf (stderr,"totchan %d\n",totnchan);

  } /* noise */

  else {

   conc = chpnt->conc;
   nconc = newconc;
   for (i=0; i<nstate; i++,conc++,nconc++) {  /* zero deltas */
     conc->cval = conc->cest;		      /* save previous values */
     nconc->dcon = 0;
   }
   conc = chpnt->conc;
   nconc = newconc;
   stpnt = chpnt->chtyp->state;

   for (i=0; i<nstate; i++,conc++,nconc++,stpnt++) { /* zero deltas, set tcr */
      dcon = icon = 0.0;
      tcrp = tcr[i];				/* pointer to precalc values */
      numtr = stpnt->numtrans;			/*per timinc */
      trans = stpnt->trans;	     		/* get parameters */
      for (j=0; j<numtr; j++,trans++,f_rate++,trmul++,tcrp++){ 
           trrate = conc->cval * *tcrp;
           newconc[*trans].dcon += trrate;
           dcon += trrate;
           icon += *tcrp;			/* save for implicit factor */
      }
      nconc->dcon -= dcon;
      nconc->icon = 1.0 / (1.0 + implfh * icon); /* implfh=1 -> implicit */
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

   critc = 1e-10;
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
     for (i=0; i<nstate; i++,stpnt++,conc++,nconc++) {	/* find change in conc*/
        trans = stpnt->trans;	     		/* get parameters */
        tcrp = tcr[i];				/* pointer to precalc values */
        numtr = stpnt->numtrans;
        for (j=0; j<numtr; j++,trans++,tcrp++){		/*per timinc */

#ifdef DEBUG
          trrate = conc->cest * *tcrp;

/* if (chpnt->ctype==NA) ncfprintf
(stderr,
  "i %d j %d t %d trrat %9.4g cest %9.4g v %g tcrp %g\n",
	i,j,t,trrate,conc->cest,chpnt->comp2->v, *tcrp); /* */
#endif

         newconc[*trans].dcon += conc->cest * *tcrp;

        }  /* for (j=0; j<numtr; ...) */
      }   /*  for (i=0; i<nstate; ...) */

	/* after calculating all deltas, find new conc estimates: */

    conc = chpnt->conc;
    nconc = newconc;
    stpnt = chpnt->chtyp->state;
    for (i=0; i<nstate; i++,stpnt++,conc++,nconc++) {

     cest = (conc->cval + implfh*nconc->dcon + (1.-implfh)*oldconc[i].dcon) 
								* nconc->icon;
     cerr = conc->cest;
     conc->cest = cest;
     cerr -= cest;
//    totconc += cest;

#ifdef DEBUG
  if (debug & NCOMP && debugz & 64)
 ncfprintf(stderr,"2 typ %d conc %10.4g dconc %10.4g err %10.4g cond %g  i %d\n",
		chpnt->ctype,conc->cest,nconc->dcon,cerr,stpnt->cond,i); /* */
#endif
     if (cerr < 0.0) cerr = -cerr; 	/* absolute value */
     if (cerr > cmaxerr) cmaxerr = cerr;

 /*ncfprintf (stderr,"cmax %10.4g\n",cmaxerr);   /* */

   }  /*  for (i=0; i<nstate; ...) */

/*  ncfprintf (stderr,"\n");   /* */


#ifdef DEBUG
  if (debug & NCOMP && debugz & 64)
     ncfprintf (stderr,"totconc %g\n",totconc);   /* */
#endif
/* ncfprintf (stderr,"cmax %9.4g\n",cmaxerr);   /* */
  }    /* for (cmaxerr;;) */
/* ncfprintf (stderr,"count %d\n",errcnt);   /* */

    conduct = 0.0;
    conc = chpnt->conc;
    stpnt = chpnt->chtyp->state;
    for (i=0; i<nstate; i++,stpnt++,conc++) { /* integrate and total cond */
      conduct += conc->cest * stpnt->cond;
  /*    if (chpnt->ctype==NA) ncfprintf (stderr,"%8.4g ",conc->cest); /* */
    }
  /* if (chpnt->ctype==NA) ncfprintf (stderr,"\n"); /* */

    chpnt->cno = int(conduct*chpnt->nchan);
  }
  if (chpnt->comp1) {			/* calc inhibition from ca binding */ 
	cacomp *capnt;

     if (chpnt->chtyp->cabnd && (capnt=chpnt->comp1->capnt)) {
       conduct *= (1-hillcalc(capnt->cais[0],chpnt->cakd,chpnt->cahc)); 
     }
  }
  return (conduct);
}

/*------------------------------------*/

#define MAXOLDCHAN 50

static chan *oldchans = (chan *)NULL;

chan *find_oldchan(chan *cp, double *conduct, double critc);

void save_oldchan(chan *cp, double conduct, double critc) 

/* save a computed equilibrium state for a channel */

{
    int i,nstate;
    chan *op,*np;
    static stconc *conc,*oconc;
    double tconduct;

				/* count up saved states */	
  for (op=oldchans,i=0; op; op=(chan*)op->next,i++) 
	{ }

//ncfprintf (stderr,"old chan state  %d\n",i);

 				/* check to see if state exists */ 
  if (i < MAXOLDCHAN) {
     if (!find_oldchan(cp, &tconduct, critc)) {
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

  /* ncfprintf(stderr,"chan  ctype stype %d %d\n",cp->ctype,cp->stype); /* */

        op->ctype   = cp->ctype;
        op->stype   = cp->stype;
        op->voffsm  = cp->voffsm;
        op->voffsh  = cp->voffsh;
        op->arate    = cp->arate;
        op->brate    = cp->brate;
        op->crate    = cp->crate;
        op->drate    = cp->drate;
        op->erate    = cp->erate;
        op->frate    = cp->frate;
        op->conduct = conduct;
	op->vrev    = cp->comp1->v;
	op->maxcond = critc;		/* save accuracy crit */
     }
  }
}

/*------------------------------------*/

chan *find_oldchan(chan *cp, double *conduct, double critc) 

/* Find a previously computed equilibrium state for a channel. */
/*  Check to make sure that all controlling params are the same, */
/*   except Ca++ or ligands. */

{
    int found;
    chan *op;

  found = 0;
  for (op=oldchans; op; op=(chan*)op->next) {/* check to see if state exists */
    if (op->vrev  == cp->comp1->v &&
	op->ctype == cp->ctype &&
        op->stype == cp->stype &&
        op->voffsm == cp->voffsm &&
        op->voffsh == cp->voffsh &&
       !op->nchan == !cp->nchan  &&       /* both same noise status ? */
        op->arate == cp->arate &&
        op->brate == cp->brate &&
        op->crate == cp->crate &&
        op->drate == cp->drate &&
        op->erate == cp->erate &&
        op->frate == cp->frate &&
	op->maxcond <= critc)
     { found = 1; break; }
  }
/* ncfprintf (stderr,"found %d %d %d op %d\n",found,cp->ctype,cp->stype,op); */

  if (found) return op;
  else       return NULL;
}

/*------------------------------------*/

int get_oldchan(chan *cp, double *conduct, double critc) 

/* Look for previously computed equilibrium state for a channel, */
/*  and if it is found, get its markov states. */

{
    int i,nstate;
    static stconc *conc,*oconc;

    chan *op;

  if (op=find_oldchan(cp, conduct, critc)) {
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
   static double (**f_rate)(chan *cpnt), trrate, conduct, *tcrp, icon;
   static double *trmul;
   static char *rateo;
   double trm;			/* multiplier for rate func */
   char *trans;
   double critcx,cerr,cmaxerr,tcomp,totconc;
   static double tcr[NUMSTATE][NUMTRANS]; /* temp store for trrate in calc */
   static dstconc *nconc,newconc[NUMSTATE];

#ifdef DEBUG
  if (debug & NCOMP && debugz & 16)
	 ncfprintf (stderr,"dochani type %d %d\n",chpnt->ctype,chpnt->stype);
#endif

  if (get_oldchan(chpnt,&conduct,critc)) {

/*  ncfprintf (stderr,"dochani: using old values %d\n",chpnt->ctype); /* */
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
    icon = 0;
    tcrp = tcr[i];                          /* pointer to precalc values */
    f_rate = stpnt->trate;
    trmul = stpnt->ratemul;
    numtr = stpnt->numtrans;
    rateo = stpnt->rateo;
    for (j=0; j<numtr; j++,f_rate++,trmul++,rateo++,tcrp++){ 
	switch (*rateo) {
	default:
	case 0: trm = chpnt->arate; break; 
	case 1: trm = chpnt->brate; break; 
	case 2: trm = chpnt->crate; break; 
	case 3: trm = chpnt->drate; break; 
	case 4: trm = chpnt->erate; break; 
	case 5: trm = chpnt->frate; break; 
        }

        *tcrp = (*f_rate)(chpnt) * *trmul * trm * tcomp;
        icon += *tcrp;

/*  ncfprintf (stderr,"i %d j %d tcrp %9.4g rate %9.4g mul %9.4g trm %g\n",
	i,j, *tcrp,(*f_rate)(chpnt),*trmul,trm); /* */

    } 						/* save precalc val */
    nconc->icon = 1.0 / (1.0 + icon);
  }	
	/* calculate delta concentrations from rate functions: */

  critcx = critc / tcomp * 1e6;

/* ncfprintf (stderr,"critcx %g critc %g tcomp %g\n",critcx,critc,tcomp); /* */

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
  ncfprintf (stderr,"i %d j %d trrate %10.4g cest %10.4g tcrp %g\n",
	i,j,trrate,conc->cest,*tcrp); /* */
       }
    }
	/* after calculating all deltas, find new conc estimates: */

    conc = chpnt->conc;
    nconc = newconc;
    stpnt = chpnt->chtyp->state;
    for (i=0; i<nstate; i++,stpnt++,conc++,nconc++) {
       cerr = -conc->cest;
       conc->cest = (conc->cval + nconc->dcon) * nconc->icon;
       cerr += conc->cest; 
       conc->cest -= cerr * .017;	/* relaxation factor */
					/* Tune this for fastest init. */ 
/* ncfprintf			
 (stderr,"init typ %d conc %9.4g dconc %9.4g err %9.4g mxerr %9.4g c %g st %d\n",
	chpnt->ctype,conc->cest,nconc->dcon,cerr,cmaxerr,stpnt->cond,i); /* */

       if (cerr < 0.0) cerr = -cerr; 	/* absolute value */
       if (cerr > cmaxerr) cmaxerr = cerr;
       nconc->dcon = 0;			/* erase delta for next iteration */
    }
/* ncfprintf (stderr,"\n");   /* */
  }    /* for (cmaxerr;;) */
/* ncfprintf (stderr,"n %d\n",n);   /* */
  nt += n;

 }  /* for (t;;) */
/* ncfprintf (stderr,"nt %d\n\n",nt);   /* */

  totconc = 0.0;
  conduct = 0.0;
  conc = chpnt->conc;
  stpnt = chpnt->chtyp->state;
  for (i=0; i<nstate; i++,stpnt++,conc++) { /* integrate and total cond */
    conduct += conc->cest * stpnt->cond;
/*  totconc += conc->cest; /* */
  }

 /* ncfprintf (stderr,"totconc %g\n",totconc);   /* */

  save_oldchan(chpnt,conduct,critc);

#ifdef DEBUG
  if (debug & NCOMP && debugz & 16) ncfprintf (stderr,"dochani end\n");
#endif

    if (chpnt->comp1) {			/* calc inhibition from ca binding */ 
          cacomp *capnt;
          double cakd;
       if (capnt=chpnt->comp1->capnt) {
         cakd = chpnt->cakd;   /* Ca binding Kd (Ca closes channel ) */
         conduct *= cakd / (capnt->cais[0] + cakd);
       }
    }
  return (conduct);
}

/*------------------------------------*/

double cicrfluxu(cacomp *capnt, double cais)

/* calculate cicr uptake flux */

{  
         double flux;

  flux = capnt->vm2 * pow(cais,capnt->ncicr) / (pow(capnt->k2cicr,capnt->ncicr) + pow(cais,capnt->ncicr)); 
  return flux;
}

/* - - - - - - - - - - - - - - - - - - */

double cicrflux(cacomp *capnt, double cas, double cais)

/* calculate cicr release flux */

{ 
     double flux;

  flux = capnt->vm3 * pow(cas, capnt->mcicr) / (pow(capnt->krcicr,capnt->mcicr) + pow(cas,capnt->mcicr)) * 
		      pow(cais,capnt->pcicr) / (pow(cais,capnt->pcicr) + pow(capnt->kacicr,capnt->pcicr));
  return flux;
} 

/* - - - - - - - - - - - - - - - - - - - -*/

double ip3Taum(cacomp * capnt, double cais)
{
       double Taum;
        Taum = 1/(capnt->b2ip3 + (capnt->a2ip3*cais));
        return Taum;
}

 /* - - - - - - - - - - - - - - - - - - - */

double ip3M(cacomp * capnt, double cais)
 
{
       double mest, dmest;
       double ip3Minf;
       double m;
 
         ip3Minf = (capnt->ip3i/(capnt->ip3i +
         capnt->d1ip3))*(cais/(cais + capnt->d4ip3));
 
       if (capnt->mtypeip3 == 0){
               m = ip3Minf;
       } else {
               mest = capnt->mip3 + timinc*(ip3Minf - capnt->mip3)/ip3Taum(capnt,cais);
               dmest = timinc*(ip3Minf - mest)/ip3Taum(capnt,cais);
               m = capnt->mip3 + (1- implfk)*(mest -
               capnt->mip3) + implfk*dmest;
       }
       capnt->mip3 = m;
       return m;
}

/* - - - - - - - - - - - - - - - - - - - -*/

double ip3TauH(cacomp * capnt, double cais, double Qfactor)
{
      double TauH = 1/(capnt->a3ip3*(Qfactor + cais));
      return TauH;
}


/* - - - - - - - - - - - - - - - - - - - -*/

double ip3H(cacomp * capnt, double cais)

{
      double h, hest, dhest, ip3Hinf, Qfactor;

      Qfactor = capnt->d2ip3*(capnt->ip3i+capnt->d1ip3)/(capnt->ip3i + capnt->d3ip3);
      ip3Hinf = Qfactor/(Qfactor + cais);
      
      hest = capnt->hip3 + timinc*(ip3Hinf -
      capnt->hip3)/(ip3TauH(capnt,cais,Qfactor));
      dhest = timinc*(ip3Hinf-hest) / (ip3TauH(capnt,cais,Qfactor));
      
      h = capnt->hip3 + (1-implfk)*(hest - capnt->hip3) +
      implfk*dhest;
      capnt->hip3 = h;
      return h;
}

/*- - - - - - - - - - - - - - - - - - -*/

 double ip3flux(cacomp *capnt, double cas2, double cais)

 /* calculate ip3 release flux */
 
{
       double flux;
 
  flux = capnt->v4ip3 * pow(ip3M(capnt,cais),capnt->mip3) * pow(ip3H(capnt,cais),capnt->hip3)*(cas2-cais) + 
		capnt->v2ip3*(cas2-cais);
  return flux;
}

/* - - - - - - - - - - - - - - - - - - -*/
 
double ip3fluxu(cacomp *capnt, double cais)
 
 /*calculate ip3 uptake flux */
 
{
       double flux;
      
  flux = capnt->v3ip3 * pow(cais,capnt->oip3) / (pow(capnt->k3ip3,capnt->oip3) + pow(cais,capnt->oip3));
  return flux;
}       


/*------------------------------------*/

void docacomp (cacomp *capnt)
                
/* Calculate calcium concentration for next time step. 
   Concentration is computed for thin shell at inner surface
   of membrane.  Calcium diffuses through several (ten) shells and then
   into the inner volume of the cell.  Calcium flux is converted
   into an equivalent concentration at a shell outside the membrane.
   This allows the diffusion equation to be run identically for
   all the compartments (except the inner and outer cores).  

   Based on Yamada, Koch, Adams (1998) and de Schutter & Smolen (1998) 

   Use equation		caflux = D * time * surf/dist * (C1 - C2).

*/

{
  register int i;
  static int j,k,n;
  double *cais,*caos,*cab,cacrit;
  int cashell,cashellp,nshell,z;
  int caoshell,caoshellp,nshello;
  comp *pnt;
  conlst *lpnt;
  chan *chpnt;
  hhchan *hhpnt;
  chanparm *chp;
  static double caratio,cadf,gfrac,chdf,chgfrac,chgvrev;
  static double vf,vf2,vfi,vfi2,cvfi,cavfi;
  static double ph,vph,ca_relg;
  static double nao, nai, ko, ki, cao;
  static double navfi, kvfi;
  static double c,h,v,vm,conduct,ica,eica,pica,ina;
  static double caflux,eex,vmax,kexo,ca1,cai,tcai,camaxerr,err;
  static double cicr_r,cicr_u,ip3_r,ip3_u,dcas,dcas2,dcais,tempcais,tempcas,tempcas2,implfh;
  static double casf0,casfn,casfc,camind,implfn,implfc,implf0;
  static double casfno,casfco,implfno,implfco,implf0o;
  static double cabt, cabti, cabf, cabr,ea,eb,ebca,dvm,vm2,ea2,eb2,eica2;
  static int caisiz = 0;		/* number of ca shell estimates */
  static double *ccaise=(double*)NULL;	/* calcium internal shell estimates */
  static double *caiso=(double*)NULL;	/* calcium shell old estimates */
  static double *caisn=(double*)NULL;	/* calcium shell neighbors */
  static double *ccaose=(double*)NULL;	/* calcium external shell estimates */
  static double *caoso=(double*)NULL;	/* calcium shell old estimates */
  static double *caosn=(double*)NULL;	/* calcium shell neighbors */
  static double *cabe=(double*)NULL;	/* calcium buffer new estimates */
  static double *caba=(double*)NULL;	/* adder to Ca on rhs */
  static double *cabm=(double*)NULL;	/* multiplier for Best */
  static double *cabea=(double*)NULL;	/* adder for rhs of Best */
  static double *implf=(double*)NULL;	/* implicit factor for whole rhs */
  static double *implfb=(double*)NULL;	/* implicit factor for whole rhs */
  register double *caise;		/* calcium internal shell estimates */
  register double *caose;		/* calcium external shell estimates */
  static double implfca;
  static int implicitca;
 
#ifdef DEBUG
  if (debug & NCOMPCA && debugz & 16) ncfprintf (stderr,"doca\n");
#endif

  if (!(pnt=capnt->comp1)) {		/* voltage compartment */
     ncfprintf (stderr,"docacomp: can't find voltage compartment\n");
     return;
  }
	// calculate offset for vrev due to surface charge

  vm = pnt->v - pnt->vext - compcavoff(pnt);
  if (pnt->pvext!=NULL) {
      ph = getnt(pnt->pvext, PH);		   /* get pH from external compartment */
      if (ph==0) ph = 7.4;
      if (ph < 6) ph = 6.0;			   /* limit to reasonable values */
      if (ph > 8) ph = 8.0;
      vph = 0.317 - 0.071*ph + 0.00381*ph*ph;  /* shift in gating voltage from ext pH, see Barnes & Bui, 1991 */
  } else {
      ph = 7.4;
      vph = 0;  			   /* shift in gating voltage from ext pH, see Barnes & Bui, 1991 */
  }					   /*  only used for Ca chans */
  pnt->vph = vph;
                                           /* from Master.m by Mike van Rijssel */
  ca_relg = (ph - 7.4) * dcaphg + 1;       /* relative conductance from ext pH, see Barnes & Bui, 1991 */
                                           /* from Master.m by Mike van Rijssel */

  implicitca = implicit;
  implfca = implfk;
  //implicitca = 1;
  //implfca = 1;
 
	/* First, calculate driving force from cao, cai, nao, nai, ko, ki */
        /*  using GHK current eqn.    */
	/*  see Hille 2nd ed, p 341-342. */

	/* Find channel's conductance for computation of Ca channel current */ 
	/*  by calculating the slope of the driving force */
	/*  See de Schutter & Smolen (1998) in "Methods in Neuronal Modeling." p 218 */

	/* Note that this calculation is done for mono- and  divalent ions */

    if (ncabs(vm) < 1e-8) vm=1e-8;
    vf = exp(vm*frt);			/* for monovalent ions */
    vf2 = vf*vf;			/* for divalent ions, exp (vm*f2rt) */
    vfi  = 1.0 / (1.0 - vf);
    vfi2 = 1.0 / (1.0 - vf2);
    cao = capnt->cao;
    if (capnt->caos!=NULL && capnt->caos[0] > 0) caratio = capnt->cais[0] / capnt->caos[0];
    else if       (cao > 0) caratio = capnt->cais[0] / cao;
    else	            caratio = capnt->cais[0] / 1e-9;

    //cavfi = capnt->caos[0]/dcao*(1.0-caratio*vf2)*vfi2;
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
  for (lpnt=pnt->clst; lpnt; lpnt=lpnt->next) {
     chpnt = (chan *)lpnt->conpnt;
     if (! chpnt) continue;
/* ncfprintf (stderr,"doca: chan %d %d\n",chpnt->ctype,chpnt->stype); /* */
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
	    hhpnt->conduct = conduct*ca_relg;
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
	    hhpnt->conduct = conduct*ca_relg;
  	    conduct = dochan2(hhpnt);		/* do noise only */
	  }
	  break;
       case 4:					/* Ca T-type */
	  chp = &hhpnt->chtyp->parm[0];
      	  chanimpl(chp,hhpnt);
	  c = hhpnt->m * chp->ival[0] + chp->ival[1];
	  chp++;
      	  chanimpl(chp,hhpnt);
	  h = hhpnt->h * chp->ival[0] + chp->ival[1];
	  hhpnt->m = c;
	  hhpnt->h = h;
  	  conduct = c*c*c*h;
	  if (hhpnt->nchan) {
	    hhpnt->conduct = conduct*ca_relg;
  	    conduct = dochan2(hhpnt);		/* do noise only */
	  }
	  break;
       case 1:
       case 3:
       case 5:
       case 6:
       case 7:
      default:
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

	 if (chpnt->ions!=NULL) {
            nao = chpnt->ions->iono[PNA];
	    ko  = chpnt->ions->iono[PK];
	 } else {
            nao = chpnt->chtyp->ions->iono[PNA];
	    ko  = chpnt->chtyp->ions->iono[PK];
	 }
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
	 if ((chtyp=chpnt->chtyp)==0) {
	    naperm = dpnaca;
     	    kperm  = dpkca;
	 }
	 else {
	    naperm = chpnt->chtyp->ions->ionp[PNA];
	    kperm  = chpnt->chtyp->ions->ionp[PK];
	 }
         /* cvfi = (1.0-caratio*vf2)*vfi2; */	/* calc based only on [Ca] */
         /* chgfrac = cvfi + f2rt*vm*(1.0-caratio)*vf2*vfi2*vfi2; */

	 caperm = get_ionp(chpnt, PCA);
         cvfi = (caperm*cavfi 
             /*  + naperm*navfi*0.25 + */		/* leave out Na from Ca chans */
	         +  kperm*kvfi*0.25   
		 );		/* for time being */

         chgfrac = cvfi + frt*vm*(2 *   caperm*(1.0-caratio)*vf2*vfi2*vfi2
			  /* +  naca*0.25*naperm*(1.0-naratio)*vf*vfi*vfi   */
		 	     +  kca*0.25*kperm*(1.0- kratio)*vf*vfi*vfi      
			         );  
       }
       if(caratio < 10.0) { 
          chdf = -vm * cvfi;
          chgvrev = chdf/chgfrac + vm;	/* "slope vrev" * frac of g */
       }
       else {
          chgvrev = capnt->vrev;	/* vrev for total channel current */
          chgfrac = 1.0;		/* conductance frac for tot current */
       } 

      conduct = conduct * chpnt->maxcond * ca_relg;
      ica += conduct * cadf;		            /* calcium current, used here only */
      chpnt->conduct = conduct * chgfrac;           /*  for use in runcomp() below. */
      chpnt->gvrev = chgvrev;		            /* Set "slope" vrev and cond */
    break;

    case NA:
    case K:
    case ClCa:
    case KCa:
    case SYN2:
    case AMPA:
    case NMDA:
    case GABA:
    case GLY:
    case CGMP:				/* retrieve prev. ica estimate */
        conduct = chpnt->conduct * get_ionp(chpnt, PCA);
        ica += conduct * cadf;		/* Ca current from other chan types, used here only */
	break;

    case AXIALRES:
    case RESISTOR:
    case DIODE:
      default: break;
	
    }	/* switch ctype */

  }

/*  ncfprintf (stderr,"ica %g\n",ica);   	/* */
  capnt->ica = ica;			/* current through ca-selective chan */

  casf0 = capnt->casf0;			/* flux factor=timinc*F2/(4PI*r*r*dr) */
  cais  = capnt->cais;
  casfn = capnt->casfn * implfca;	/* shell factor = 1 / (dr*dr) */
  casfc = capnt->casfc * implfca;	/* shell factor for core */

  caos   = capnt->caos;
  casfno = capnt->casfno * implfca;	/* shell factor = 1 / (dr*dr) */
  casfco = capnt->casfco * implfca;	/* shell factor for core */

  ca1 = capnt->cais[0];
  camind = dtcai;			/* min for pump & diffusion */ 
  tcai = ca1-camind;
  tcai = max (tcai,0);

  if ((kexo=capnt->kexo) > 0) {	/* sodium-calcium exchanger */
					/* From de Schutter & Smolen 1998 */
					/* after DiFrancesco & Noble 1985 */
				        /* Gabbiani,Midtgaard and Knopfel 1994*/
     if (caos!=NULL)
           kexo *= caos[0];
     else
           kexo *= cao;

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

     if (!implicitca) {		/* CN */

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
/*     ncfprintf (stderr,"ica %g vmax %g tcai %g pkm %g\n",ica,vmax,tcai,capnt->pkm); /* */
  }
  else pica = 0.0;

		 /* Compute Ca concentration at inner mem surface. */

/* ncfprintf (stderr,"ica %g pica %g eica %g ina %g\n",ica,pica,eica,ina); /* */

  capnt->ipump = pica * dicafrac;	    /* Ca pump current */
  capnt->iexch = (eica + ina) * diexchfrac; /* total exchanger current */

					/*       = timestep * M/sec/A/vol(dm) */
  caflux = ica * casf0;   		/* flux in M (moles/liter) for 1st sh */


		/* Dan Emerson's CICR code */
		/* Averaged implicit + explicit CICR */

  if (capnt->vm3 > 0) { 
        double cais0, implfcas, dcais_est, dcas_est, dcas2_est, dcais_err;

    implfcas = implfca;

    cais0 = capnt->cais[0]; 					/* [Ca]i in outermost shell, M/L */
    cicr_u = cicrfluxu(capnt, cais0);				/* uptake from outer Ca shell, M/s */
    cicr_r = cicrflux (capnt, capnt->cas, cais0);		/* release into outer Ca shell, M/s */
    ip3_u = ip3fluxu(capnt,cais0);				/* uptake from outer Ca shell, M/s */
    ip3_r = ip3flux(capnt,capnt->cas2, cais0);			/* release into outer Ca Shell, M/s */

      /* make initial estimate for beginning of time step */
      /* dcais calib in M/L */

    dcais = timinc * (capnt->c1cicr * 
			(cicr_r - cicr_u + capnt->kfcicr *  capnt->cas + 
			  ip3_r - ip3_u) - capnt->k1cicr*cais0 + capnt->vip3*capnt->bip3);

    dcas  = timinc *    (cicr_u - cicr_r - capnt->kfcicr*capnt->cas);
    dcas2 = timinc *    ( ip3_u - ip3_r );

    dcais_est = dcais;		/* save for implicit est */
    dcas_est  = dcas;
    dcas2_est  = dcas2;

      /* make new estimate for end of time step, iterate for implicit */

    cacrit = 1e-10;
    for (n=0; n<100 && dcais_err > cacrit; n++) {

        tempcais = cais0      + dcais_est;
        tempcas  = capnt->cas + dcas_est;
        tempcas2  = capnt->cas2 + dcas2_est;

        cicr_u = cicrfluxu(capnt, tempcais);
        cicr_r = cicrflux (capnt, tempcas, tempcais);

        ip3_u = ip3fluxu(capnt,tempcais);
        ip3_r = ip3flux(capnt,tempcas2, tempcais);
 
        dcais_est = timinc * (capnt->c1cicr * 
                 (cicr_r - cicr_u + capnt->kfcicr*tempcas + 
		  ip3_r -  ip3_u) - capnt->k1cicr*tempcais + capnt->vip3*capnt->bip3);

	dcas_est  = timinc * (cicr_u - cicr_r - capnt->kfcicr*tempcas);
	dcas2_est = timinc * (ip3_u -  ip3_r);

        dcais_err = ncabs(dcais_est - dcais);
    }

	// Add averaged cicr flux to total caflux

    capnt->cicrflux = implfcas*dcais_est + (1-implfcas)*dcais; 
    caflux += capnt->cicrflux;

	// Add averaged cas flux to original cas

    capnt->cas  += implfcas*dcas_est  + (1-implfcas)*dcas; 
    capnt->cas2 += implfcas*dcas2_est + (1-implfcas)*dcas2;

  /* ncfprintf (stderr,"vm %g ica %g caflux %g eb*cai %g\n",
			vm,ica,caflux,ebca*cais[0]); /* */

  } /* end of CICR calculations */

  cab   = capnt->cab;
  cashell = capnt->cashell;		/* number of shells */
  cashellp = cashell - 1;		/* loop iteration limit */
  nshell  = cashell - 1;
  if (nshell <= 1) nshell=1;

  caoshell = capnt->caoshell;		/* number of shells */
  caoshellp = caoshell - 1;		/* loop iteration limit */
  nshello  = caoshell - 1;
  if (nshello <= 1) nshello=1;
					/* make temporary space */
					/* to compute diffusion */
  if (cashell > caisiz) {		/* if bigger than previous */
   caisiz = cashell;
   if (ccaise) efree(ccaise);		/* free old space */
   if ((ccaise=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
    ncfprintf(stderr,"docacomp: no space left for %d calcium shells\n",nshell);
    return;
   }
   if (caiso) efree(caiso);		/* free old space */
   if ((caiso=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
      ncfprintf(stderr,"docacomp: no space left for %d calcium shells\n",nshell);
      return;
   }
   if (caisn) efree(caisn);
   if ((caisn=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
     ncfprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
      return;
   }
					/* external shells */
   if (ccaose) efree(ccaose);		/* free old space */
   if ((ccaose=(double *)emalloc((caoshell)*sizeof(double))) == (double*)NULL) {
    ncfprintf(stderr,"docacomp: no space left for %d calcium shells\n",nshello);
    return;
   }
   if (caoso) efree(caoso);		/* free old space */
   if ((caoso=(double *)emalloc((caoshell)*sizeof(double))) == (double*)NULL) {
      ncfprintf(stderr,"docacomp: no space left for %d calcium shells\n",nshello);
      return;
   }
   if (caosn) efree(caosn);
   if ((caosn=(double *)emalloc((caoshell)*sizeof(double))) == (double*)NULL) {
     ncfprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshello);
      return;
   }
  }
  if (cab) {			/* temporary calc arrays for Ca buffer */
	static int cabsiz;
    if(cashell > cabsiz) {
     cabsiz = cashell;
     if (cabe) efree(cabe);
     if ((cabe=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       ncfprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (implf) efree(implf);
     if ((implf=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       ncfprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (caba) efree(caba);
     if ((caba=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       ncfprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (cabm) efree(cabm);
     if ((cabm=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       ncfprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (cabea) efree(cabea);
     if ((cabea=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       ncfprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (implfb) efree(implfb);
     if ((implfb=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       ncfprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
    }
  }
 
  caise = ccaise;
  for (i=0; i<cashell; i++) {		/* reset old estimates */
    caise[i] = cais[i];
    caiso[i] = 0.0;
  }

  if (caos!=NULL) {
    caose = ccaose;
    for (i=0; i<caoshell; i++) {		/* reset old estimates */
      caose[i] = caos[i];
      caoso[i] = 0.0;
    }
  }

 implf0 = 1.0 + casfn + ebca;
 implfn = 1.0 + 2.0*casfn;
 implfc = 1.0 + casfc;

 implf0o = 1.0 + casfno + ebca;
 implfno = 1.0 + 2.0*casfno;
 implfco = 1.0 + casfco;

 if (cab) {				/* set up implicit factors */
   cabt  = capnt->cabt;
   cabti = capnt->cabti;
   cabf  = capnt->cabf;
   cabr  = capnt->cabr;
   for (i=0; i<cashell; i++) {		/* reset old estimates */
     cabe[i] = cab[i];
   }
   if (implicitca) {			/* precompute for iteration below */
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

if (!implicitca) {
  for (i=0; i<cashell; i++) {	
     if (i==0)           caisn[i] = cais[i+1] - cais[i];
     else if (i<cashellp)caisn[i] = cais[i-1] + cais[i+1] - 2.0*cais[i];
     else                caisn[i] = cais[i-1] - cais[i];
  }
  if (caos!=NULL) {
    for (i=0; i<caoshell; i++) {	
      if (i==0)            caosn[i] = caos[i+1] - caos[i];
      else if (i<caoshellp)caosn[i] = caos[i-1] + caos[i+1] - 2.0*caos[i];
      else                 caosn[i] = caos[i-1] - caos[i];
    }
  }
}

 cacrit = 1e-10;
 for (camaxerr=1e8,n=0; camaxerr>cacrit; n++) {	/* relaxation */

 if (camaxerr > 1e8 || n>1000) {
    ncfprintf (stderr,"# Cacomp: error, no convergence, continuing, pid %d\n",getpid());
    break;
 }
 camaxerr = 0.0;
					/* solve ca concentration */

 if (cab) {				/* active buffer */
   if (implicitca) {

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
 } /* cab */
 else {					/* passive buffer */
   if (implicitca) {
							/* first shell */
     caise[0]=(cais[0]+caflux+(caise[1])*casfn)/implf0; 
     if (caise[0] < camind)  caise[0] = camind;

    for (i=1; j=i-1,k=i+1,i<cashellp; i=k) {	  
     caise[i]=(cais[i]+(caise[j]+caise[k])*casfn)/implfn; /* normal shells */
    }
    caise[i]= (cais[i] + caise[j]*casfc) / implfc;   /* core */
/* ncfprintf (stderr,"caise0 %g caise1 %g\n",caise[0],caise[1]); */

   }
 else {		/* CN */
//ncfprintf (stderr,"111 t %5.3g caflux %8.3g caise0 %20.14g caise1 %20.14g caisn0 %20.14g casfn %20.14g implf0  %20.14g\n",
//		simtime,caflux,caise[0],caise[1],caisn[0],casfn,implf0);
     caise[0]= (cais[0]+caflux+(caise[1]+caisn[0])*casfn) / implf0;
     if (caise[0] < camind)  caise[0] = camind;

//ncfprintf (stderr,"222 t %5.3g caflux %8.3g caise0 %20.14g caise1 %20.14g caisn0 %20.14g casfn %20.14g implf0  %20.14g\n",
//		simtime,caflux,caise[0],caise[1],caisn[0],casfn,implf0);

    for (i=1; j=i-1,k=i+1,i<cashellp; i=k) {
     caise[i]= (cais[i]+
       (caise[j]+caise[k]+caisn[i]) * casfn) / implfn;
     }
					/* core */
      caise[i]= (cais[i] + (caise[j]+caisn[i])*casfc) / implfc;
/* ncfprintf (stderr,"caise0 %g caise1 %g\n",caise[0],caise[1]); /* */
    }
  }

 if (nshello>0 && caos!=NULL) {					/* passive buffer */
   if (implicitca) {
							/* first shell */
     caose[0]=(caos[0]-caflux+(caose[1])*casfno)/implf0; 
     if (caose[0] < camind)  caose[0] = camind;

    for (i=1; j=i-1,k=i+1,i<caoshellp; i=k) {	  
     caose[i]=(caos[i]+(caose[j]+caose[k])*casfno)/implfno; /* normal shells */
    }
    caose[i]= (caos[i] + caise[j]*casfco) / implfco;   /* core */
/* ncfprintf (stderr,"caose0 %g caose1 %g\n",caose[0],caose[1]); */

   }
 else {		/* CN */
//ncfprintf (stderr,"111 o %5.3g caflux %8.3g caose0 %20.14g caose1 %20.14g caosn0 %20.14g casfno %20.14g implf0o %20.14g\n",
//		simtime,caflux,caose[0],caose[1],caosn[0],casfno,implf0o);
     caose[0]= (caos[0]-caflux+(caose[1]+caosn[0])*casfno) / implf0o;
     if (caose[0] < camind)  caose[0] = camind;

//ncfprintf (stderr," %5.3g caflux %8.3g\n",simtime,caflux);
//ncfprintf (stderr,"222 o %5.3g caflux %8.3g caose0 %20.14g caose1 %20.14g caosn0 %20.14g casfno %20.14g implf0o %20.14g\n",
//		simtime,caflux,caose[0],caose[1],caosn[0],casfno,implf0o);

    for (i=1; j=i-1,k=i+1,i<caoshellp; i=k) {
     caose[i]= (caos[i]+
       (caose[j]+caose[k]+caosn[i]) * casfno) / implfno;
     }
					/* core */
      caose[i]= (caos[i] + (caose[j]+caosn[i])*casfco) / implfco;
/* ncfprintf (stderr,"caose0 %g caose1 %g\n",caose[0],caose[1]); /* */
    }
  }
				/* find error */
  for (i=0; i<cashell; i++) {
    err = caise[i] - caiso[i];
    if (ncabs(err) > camaxerr) camaxerr = err;
//    caise[i] += err * .001;
    caiso[i] = caise[i];		/* save old estimate */
  }
  if (caos!=NULL) {
    for (i=0; i<caoshell; i++) {
      err = caose[i] - caoso[i];
//fprintf (stderr,"erro %g\n",err);
      if (ncabs(err) > camaxerr) camaxerr = err;
  //    caose[i] += err * .001;
      caoso[i] = caose[i];		/* save old estimate */
    }
  }

/*  ncfprintf (stderr,"doca: camaxerr %g\n",camaxerr);  /* */

 }	/* for (camaxerr;;) */

/* ncfprintf (stderr,"n %d\n",n);  /* */

  for (i=0; i<cashell; i++) {
       cais[i] = caise[i];			/* save estimate */
  }
  if (caos!=NULL) {
    for (i=0; i<caoshell; i++) {
       caos[i] = caose[i];			/* save estimate */
    }
  }
  if (cab) 
    for (i=0; i<cashell; i++) {
       cab[i] = cabe[i];
    }
		 /* Compute calcium reversal potential */
		 /* Use Nernst equation */

/*  if (time > .014) {
    i = i+1;
  } 			/* */

/* compute new exchr current */

 if (kexo > 0) {
     eica  = ea  - capnt->cais[0] * eb;		
     eica2 = ea2 - capnt->cais[0] * eb2;
     if (implicitca) {
        capnt->slope = (1-NAEXCHR/CAZ)*(eica2-eica) / dvm;
        capnt->eica  = (1-NAEXCHR/CAZ)*eica - vm*capnt->slope;
     }
     else {
        capnt->slope = (capnt->slope +(1-NAEXCHR/CAZ)*(eica2-eica)/dvm)*0.5;
        capnt->eica  = (capnt->eica  +(1-NAEXCHR/CAZ)*eica-vm*capnt->slope)*0.5;
     }
 }
			/* Nernst potential for Ca. */
			/* Not completely correct for Ca chans, as they need *
			/* modified GHK voltage equation for [K]i */
			/* see Hille 2nd ed, p 107 & p 344        */

  if (capnt->caos!=NULL)
     capnt->vrev = r2ft*log(capnt->caos[0]/capnt->cais[0]); 
  else
     capnt->vrev = r2ft*log(capnt->cao/capnt->cais[0]); 

#ifdef DEBUG
  if (debug & NCOMPCA && debugz & 16) ncfprintf (stderr,"doca end\n");
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
   all the compartments (except the inner and outer cores).  

   Based on Yamada, Koch, Adams (1998) and de Schutter & Smolen (1998) 

   Use equation		caflux = D * time * surf/dist * (C1 - C2).

*/

{
  register int i;
  int j,k,n;
  double *cais,*caos,*cab,timeimpl,cacrit;
  int cashell,cashellp,nshell,z;
  int caoshell,caoshellp,nshello;
  comp *pnt;
  conlst *lpnt;
  chan *chpnt;
  static double caflux,camaxerr,err;
  static double casf0,casfn,casfc,camind,implfn,implfc,implf0;
  static double casfno,casfco,implfno,implfco;
  static double cabt, cabti, cabf, cabr;
  static int caisiz = 0;		/* number of ca shell estimates */
  static double *ccaise=(double*)NULL;	/* calcium internal shell estimates */
  static double *caiso=(double*)NULL;	/* calcium int shell old estimates */
  static double *ccaose=(double*)NULL;	/* calcium external shell estimates */
  static double *caoso=(double*)NULL;	/* calcium ext shell old estimates */
  static double *cabe=(double*)NULL;	/* calcium buffer new estimates */
  static double *caba=(double*)NULL;	/* adder to Ca on rhs */
  static double *cabea=(double*)NULL;	/* adder for rhs of Best */
  register double *caise;		/* calcium internal shell estimates */
  register double *caose;		/* calcium external shell estimates */
  
#ifdef DEBUG
  if (debug & NCOMPCA && debugz & 16) ncfprintf (stderr,"doca\n");
#endif

  if (!(pnt=capnt->comp1)) {		/* voltage compartment */
     ncfprintf (stderr,"docacompi: can't find voltage compartment\n");
     return;
  } 

  caflux = 0;   			/* flux in moles */
  camind = dtcai;			/* min for pump & diffusion */ 

/* ncfprintf (stderr,"ica %g caflux %g\n",ica,caflux);   /* */

  cab   = capnt->cab;
  cashell = capnt->cashell;		/* number of shells */
  caoshell = capnt->caoshell;		/* number of shells outside */
  cashellp = cashell - 1;		/* loop iteration limit */
  caoshellp = caoshell - 1;		/* loop iteration limit */
  nshell  = cashell - 1;
  nshello  = caoshell - 1;
  if (nshell <= 1) nshell=1;
  if (nshello <= 1) nshello=1;
					/* make temporary space */
					/* to compute diffusion */
  if (cashell > caisiz) {		/* if bigger than previous */
					/* internal shells */
   caisiz = cashell;
   if (ccaise) efree(ccaise);		/* free old space */
   if ((ccaise=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
    ncfprintf(stderr,"docacomp: no space left for %d calcium shells\n",nshell);
    return;
   }
   if (caiso) efree(caiso);		/* free old space */
   if ((caiso=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
      ncfprintf(stderr,"docacomp: no space left for %d calcium shells\n",nshell);
      return;
   }
					/* external shells, no buffer */
   if (ccaose) efree(ccaose);		/* free old space */
   if ((ccaose=(double *)emalloc((caoshell)*sizeof(double))) == (double*)NULL) {
    ncfprintf(stderr,"docacomp: no space left for %d calcium shells\n",nshello);
    return;
   }
   if (caoso) efree(caoso);		/* free old space */
   if ((caoso=(double *)emalloc((caoshell)*sizeof(double))) == (double*)NULL) {
      ncfprintf(stderr,"docacomp: no space left for %d calcium shells\n",nshello);
      return;
   }

   if (cab) {			/* temporary calc arrays for Ca buffer */
     if (cabe) efree(cabe);
     if ((cabe=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       ncfprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (caba) efree(caba);
     if ((caba=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       ncfprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
     if (cabea) efree(cabea);
     if ((cabea=(double *)emalloc((cashell)*sizeof(double))) == (double*)NULL) {
       ncfprintf(stderr,"docacomp: no space left for %d calcium buffers\n",nshell);
        return;
     }
   }
  }

  timeimpl = 1e6/timinc;		/* timinc = 1e6 seconds */ 
  cais  = capnt->cais;
  caos  = capnt->caos;
  casf0 = capnt->casf0 * timeimpl;	/* caflux factor = dt*2F/(4PI*r*r*dr)*/
  casfn = capnt->casfn * timeimpl; 	/* shell factor = 1 / (dr*dr) */
  casfc = capnt->casfc * timeimpl; 	/* shell factor for core */

  casfno = capnt->casfno * timeimpl; 	/* shell factor = 1 / (dr*dr) */
  casfco = capnt->casfco * timeimpl; 	/* shell factor for core */

  caise = ccaise;
  caose = ccaose;

  for (i=0; i<cashell; i++) {		/* reset old estimates */
    caise[i] = cais[i];
    caiso[i] = 0.0;
  }
  if (caos!=NULL) {
    for (i=0; i<caoshell; i++) {	/* reset old estimates */
      caose[i] = caos[i];
      caoso[i] = 0.0;
    }
  }

 implf0 = 1.0 + casfn;
 implfn = 1.0 + 2.0*casfn;
 implfc = 1.0 + casfc;
 implfno = 1.0 + 2.0*casfno;
 implfco = 1.0 + casfco;

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
    ncfprintf (stderr,"# Cacompi: error, no convergence, continuing, pid %d\n",getpid());
    break;
 }
 camaxerr = 0.0;
					/* solve ca concentration */

 if (cab) {				/* active buffer */
     caise[0]=(caba[0]+(caise[1])*casfn - cabr*cabe[0]) /
		(implf0 + cabf*cabe[0]); /* first shell */
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
 if (caos!=NULL) {
   for (i=0; i<caoshell; i++) {
    err = caose[i] - caoso[i];
    if (ncabs(err) > camaxerr) camaxerr = err;
    caose[i] += err * .6;
    caoso[i] = caose[i];		/* save old estimate */
   }
 }

/*  ncfprintf (stderr,"doca: camaxerr %g\n",camaxerr);  /* */

 }	/* for (camaxerr;;) */

/* ncfprintf (stderr,"n %d\n",n);  /* */

  for (i=0; i<cashell; i++) {
       cais[i] = caise[i];			/* save estimate */
  }
  if (caos!=NULL) {
    for (i=0; i<caoshell; i++) {
       caos[i] = caose[i];			/* save estimate */
    }
  }
  if (cab) 
    for (i=0; i<cashell; i++) {
       cab[i] = cabe[i];
    }
			 /* Compute calcium reversal potential */
			 /* Use Nernst equation */

			/* Nernst potential modified by GHK for [K]i */
			/* see Hille 2nd ed, p 107 & p 344           */

  if (caos!=NULL) 
    capnt->vrev = r2ft*log(capnt->caos[0]/capnt->cais[0]); 
  else
    capnt->vrev = r2ft*log(capnt->cao/capnt->cais[0]); 

		/* Compute reversal potential for calcium channels */

  for (lpnt=pnt->clst; lpnt; lpnt=lpnt->next) {
     chpnt = (chan *)lpnt->conpnt;
     if (! chpnt) continue;
     if (chpnt->ctype == CA && !chpnt->setvrev) {
        if (caos!=NULL) 
		 chpnt->vrev = ghkv(getchantype(chpnt->ctype,chpnt->stype),
					capnt->caos[0],capnt->cais[0]);
        else
		 chpnt->vrev = ghkv(getchantype(chpnt->ctype,chpnt->stype),
					capnt->cao,capnt->cais[0]);
      }
  }

#ifdef DEBUG
  if (debug & NCOMPCA && debugz & 16) ncfprintf (stderr,"doca end\n");
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
       if (chpnt->cstate) drand_setstate (chpnt->cstate);
       nco = binomdev(pa,nc); 			/* delta closed to open */
       noc = binomdev(pb,no); 			/* delta open to closed */
       no += nco - noc;
       chpnt->cno = no;
       cconduct = ((double)no) / nchan;
       restorstate();
     }
  }
  return cconduct;
}

#endif

/*------------------------------------*/

void doverr(register comp *pnt, int noskip, double critxx);

void runcomp(void)

/* Advance the network of compartments one step in time,
using a second-order implicit numerical integration.
Iterate the approximation until the error is
reduced below criterion.
*/

{
  register comp *pnt;
  register conlst *lpnt;
  static unsigned int left=1;
  static int ncompl,ncompr,tcomp;
  static int niter,titer,less,noskip=0,ocpntf=0;
  static int stype;
  static double nodcur;                    /* total node current */
  static double nodcr;                     /* total unvarying node curr */
  static double ocais;                     /* original cais[0] */
  static conn *conpnt;
  static comp *ocpnt;
  static double delcrit;
  static double relmult,critx,critxx,critk;
  static double tcond;			   /* total conductance in comp */
  static double tcondn;			   /* tot cond not varying w/vest */
  static double conduct;		   /* conductance of channel */
  static double vest;		   	   /* conductance of channel */
  static recpar *cpt;                      /* photrec constants */
  static hhchan *hhpnt;
  static kcachan *kcapnt;
  static kcachan *kcapnts;
  static chan *chpnt;
  static cacomp *capnt;
  static chanparm *chp;
  static double a,m,n,h,dm,dn,dh;
  static double mmaxerr;		   /* max maxerr for level of critx */
  static double oldmaxerr;		   /* looking for maxerr increases */
  static double oimplfk;		   /* save for implfk */
  static int oimplicit;		   /* save for implicit */
  double avgrelax;
  int nrelax;

#ifdef DEBUG
  if (debug & NCOMP && debugz & 1) ncfprintf (stderr,"runcomp\n");
#endif

 if (timinc > 1.0) runsyn(0);		/* if static, calc synapses */
 if (euler) implicit = 0;

// if (timinc >= 5e-5) {
//   implicit = 1;
//   implfk = 1;
// }
// implfk = .5;
 oimplfk = implfk;			/* save implicit mode for change */
 oimplicit = implicit;			

 for (pnt=compnt; pnt; pnt=pnt->next) {
   pnt->oldv = pnt->v;		/* save old voltage, allows save/restore */
 }

 for (pnt=compnt; pnt; pnt=pnt->next)                   /* First estimate */
  {                                                     /* for all comps */
   nodcur = 0.0;		/* current that is also calc in docomp() */
   nodcr = 0.0;			/* current that is precalc for docomp() */
   tcond = 0.0;			/* conductance that changes to another node */
   tcondn = 0.0;		/* conductance that changes in channels to outside */

   if (pnt->pvext!=NULL) pnt->vext = pnt->pvext->v; /* set external voltage for cacomp */
   else                  pnt->vext = 0;

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

        case GJ:   
        case PNX:     if (((gj*)conpnt)->maxcond == 0) continue; break;
        case SYNAPSE: if (((synap*)conpnt)->maxcond == 0) continue; break;

        case ROD:
        case CONE:    if (((photrec*)conpnt)->maxcond == 0) continue; break;

        case TRANSDUCER: break;

        case SYN2:
        case AMPA:
        case CGMP:
        case GABA:
        case GLY:
        case NMDA:
        case CA:
        case NA:
        case ClCa:
        case KCa:
        case K:       if (((chan*)conpnt)->maxcond == 0) continue; break;
       default:   break;
     }

     ocpntf = 0;
     switch (conpnt->ctype) {
	    case SYN2:
	    case AMPA:
	    case CGMP:
	    case GABA:
	    case GLY:
	    case NMDA:
            case NA: 
            case K: 
            case KCa: 
            case ClCa:   
            case CA:   
		   if (pnt==conpnt->compe) {	/* this comp is vext, other comp computes conductance */
		     chpnt = (chan *)conpnt;
		     conduct = chpnt->conduct * get_gfrac(chpnt,chpnt->comp1->v);
                     // nodcr -= conduct * (chpnt->vrev - (chpnt->comp1->v - pnt->v));  /* sub instead of add */
                     //  add GHK current equation here 
                     nodcr -= conduct * (get_gvrev(chpnt,chpnt->comp1->v) - 
				                        (chpnt->comp1->v - pnt->v));  /* sub instead of add */
                     tcondn += conduct;    		/* total cond */
		     ocpntf = 1;
		   }
		   break;
       case SYNAPSE:
		    if (pnt==conpnt->compe) { /* this comp is vext, other comp computes conductance */
		      chpnt = (chan *)conpnt;
		      conduct = chpnt->conduct;
                      nodcr -= conduct * (chpnt->vrev-(chpnt->comp2->v-pnt->v)); /* sub instead of add */
                      tcondn += conduct;    		/* total cond */
		      ocpntf = 1;
		    }
		    break;

            default:   
		    ocpntf = 0;
		    break;
     }

     if (!ocpntf)
      switch (conpnt->ctype)
      {                         /* add connection conductances which change */

       case GJ:
       case PNX:
                nodcur += conpnt->conduct * ocpnt->v;   /* current thru g j */
                tcond  += conpnt->conduct;              /* total cond */
               break;
       case DIODE:
		if ((pnt->v - ocpnt->v) < 0) {
                  nodcur += conpnt->conduct * ocpnt->v;  /* cur to next comp */
                  tcond  += conpnt->conduct;              /* total cond */
		}
                break;              /* conductance is collected in pnt->tcond */
       case SYNAPSE:					 /* cur thru syn chan */
		if (!(((synap*)conpnt)->resp1 || ((synap*)conpnt)->resp2)) {	/* regular synapse */
		if (pnt==((synap *)conpnt)->comp2) {
                  nodcr += ((synap *)conpnt)->conduct * (((synap *)conpnt)->vrev-compcavoff(pnt));
          	  nodcur -= ((synap *)conpnt)->conduct * (pnt->v - vext(conpnt));  /* current due to Vm */
                  tcondn += ((synap *)conpnt)->conduct;   /* total cond */
  	/* ncfprintf (stderr,"vrev %10.3g conduct %10.3g\n",
		((synap *)conpnt)->vrev,((synap *)conpnt)->conduct); /* */
		}
	       } break;

      case TRANSDUCER: break;

      case ROD:
      case CONE:                                /* find photrec instantan */
                                                /* conductance */
          cpt = ((photrec *)conpnt)->chtyp;
          nodcr += ((photrec *)conpnt)->conduct * (cpt->vrev-compcavoff(pnt));  /* cur thru rec chan*/
          nodcur -= ((photrec *)conpnt)->conduct * (pnt->v - vext(conpnt)); 	/* current due to Vm */
          tcondn  += ((photrec *)conpnt)->conduct;               /* total cond */
/*   ncfprintf (stderr,"vrev %10.3g conduct %10.3g\n",
		cpt->vrev,((photrec *)conpnt)->conduct); /* */
         break;

      case LOAD:
		conduct = ((load *)conpnt)->conduct;
                nodcr += conduct * (((load *)conpnt)->vrev-compcavoff(pnt));	/* cur to "ground" */
                tcond += conduct;    				/* total cond */
                break;              /* conductance is *not* collected in pnt->tcond */

      case RESISTOR:
		conduct = conpnt->conduct;
                nodcur += conduct * ocpnt->v;		/* cur to next comp */
                tcond += conduct;    			/* total cond */
                break;              /* conductance is *not* collected in pnt->tcond */

      case AXIALRES:
                nodcur += conpnt->conduct * ocpnt->v;  /* cur to next comp */
                break;              /* conductance is collected in pnt->tcond */

      case BATT:
     		if (ocpnt == conpnt->comp1)
                  nodcur += (ocpnt->v - conpnt->conduct) * BATTCOND;
		else 
                  nodcur += (ocpnt->v + conpnt->conduct) * BATTCOND;
                break;		/* cur to next comp through batt */

      case BUF:	if (pnt == ((dbuf *)conpnt)->comp2) {
                     if (((dbuf *)conpnt)->delay || ((dbuf *)conpnt)->filt) {
		          pnt->extv = ((dbuf*)conpnt)->v;
		     }
		     else pnt->extv = (((dbuf *)conpnt)->comp1->v - ((dbuf *)conpnt)->offset) *
		     		       ((dbuf *)conpnt)->gain;
		  }
		break;

      case NBUF: if (pnt == ((ndbuf *)conpnt)->comp2) {
		    setnt (((ndbuf *)conpnt)->comp2, ((ndbuf *)conpnt)->ntrans, 
			( ((ndbuf *)conpnt)->comp1->v - ((ndbuf *)conpnt)->offset )  *
			  ((ndbuf *)conpnt)->gain +     ((ndbuf *)conpnt)->ntoffset ); 
		 }
		break;

      case CAP: nodcr += (pnt->oldv - ocpnt->oldv) * conpnt->conduct/timinc;
                break;			/* -(v3 - v2) * c / t */

      case CGMP:		/* Channels possibly with Ca current */
      case NMDA:
      case AMPA:
      case GABA:
      case GLY:
      case SYN2:
		chpnt = (chan *)conpnt;
  		conduct = dochan2(chpnt);
  		conduct *= chpnt->maxcond * get_gfrac(chpnt,pnt->v);
		chpnt->conduct = conduct;
                nodcr += conduct * (get_gvrev(chpnt,pnt->v) - compcavoff(pnt));
                nodcur -= conduct * (pnt->v - vext(conpnt)); 		/* current due to Vm */
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
  		conduct *= hhpnt->maxcond * get_gfrac(hhpnt,pnt->v);
		hhpnt->conduct = conduct;
                nodcr += conduct * (get_gvrev(hhpnt,pnt->v) - compcavoff(pnt));
                nodcur -= conduct * (pnt->v - vext(conpnt));  		/* current due to Vm */
                tcondn += conduct;    			/* total cond */
		break;

		case 1:				/* Sequential state */
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		case 8:
		default:

	 	chpnt = (chan *)conpnt;
  		conduct = dochan2(chpnt);
  		conduct *= chpnt->maxcond * get_gfrac(chpnt,pnt->v);
		chpnt->conduct = conduct;
                nodcr += conduct * (get_gvrev(chpnt,pnt->v) - compcavoff(pnt));
                nodcur -= conduct * (pnt->v - vext(conpnt)); 		/* current due to Vm */
                tcondn += conduct;    					/* total cond */
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
  		conduct *= hhpnt->maxcond * get_gfrac(hhpnt,pnt->v);
		hhpnt->conduct = conduct;
                nodcr += conduct * (get_gvrev(hhpnt,pnt->v) - compcavoff(pnt));
                nodcur -= conduct * (pnt->v - vext(conpnt));  			/* current due to Vm */
                tcondn += conduct;    			/* total cond */
                break;

		case 1:				/* Sequential state */
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
		case 8:
		case 9:
		default:
	 	chpnt = (chan *)conpnt;
  		conduct = dochan2(chpnt);
  		conduct *= chpnt->maxcond * get_gfrac(chpnt,pnt->v);
		chpnt->conduct = conduct;
                nodcr += conduct * (get_gvrev(chpnt,pnt->v) - compcavoff(pnt));
                nodcur -= conduct * (pnt->v - vext(conpnt));  			/* current due to Vm */
                tcondn += conduct;    		/* total cond */
		break;

		case 2:				/* Type A channel */
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
  		conduct *= hhpnt->maxcond * get_gfrac(hhpnt,pnt->v);
		hhpnt->conduct = conduct;
                nodcr += conduct * (get_gvrev(hhpnt,pnt->v) - compcavoff(pnt));
                nodcur -= conduct * (pnt->v - vext(conpnt));  			/* current due to Vm */
                tcondn += conduct;    			/* total cond */
                break;

		}	/* switch (stype) */
		break;

      case ClCa:   
	        kcapnts = (kcachan *)conpnt;
		if (pnt->capnt) {  /* calcium compartment exists */
		  if (!kcapnts->initfl){/* Calcium not set up before K chan. */
		     dochani((chan *)kcapnts,1e-9);
		     kcapnts->initfl = 1;
		  }
		}
  		conduct = dochan2(kcapnts);
  		conduct *= kcapnts->maxcond * get_gfrac(kcapnts,pnt->v);
		kcapnts->conduct = conduct;
                nodcr += conduct * (get_gvrev(kcapnts,pnt->v) - compcavoff(pnt));
                nodcur -= conduct * (pnt->v - vext(conpnt));  			/* current due to Vm */
                tcondn += conduct;    			/* total cond */
		break;
      case KCa:   
	 	hhpnt = (hhchan *)conpnt;
		stype = hhpnt->stype;
		switch (stype) {

		case 0:				/* Calcium-activated K chan */
		case 2:	
	 	kcapnt = (kcachan *)conpnt;
		chp = &kcapnt->chtyp->parm[0];
		if (pnt->capnt) {  /* calcium compartment exists */
		  if (!kcapnt->initfl){	/* Calcium not set up before K chan. */
			double ca;	/* Find equilibrium value for "n". */ 
			double alph,bet;
		    ca = pnt->capnt->cais[0];
		    alph = akcacalc((pnt->v-vext(conpnt)-ccavoff(kcapnt)-kcapnt->voffsm),ca,
			kcapnt->arate,kcapnt->d1,kcapnt->k1);
		    bet = bkcacalc((pnt->v-vext(conpnt)-ccavoff(kcapnt)-kcapnt->voffsm),ca,
			kcapnt->arate,kcapnt->d2,kcapnt->k2);
		    kcapnt->m = alph / (alph + bet);
		    kcapnt->initfl = 1;
		    ocais = ca;
		  }
 		  kcatab((pnt->v-vext(conpnt)-ccavoff(kcapnt)-kcapnt->voffsm),    /* set K rate chtyp */ 
		     (ocais+pnt->capnt->cais[0])*.5,stype,kcapnt->arate,
		     kcapnt->d1,kcapnt->d2,kcapnt->k1,kcapnt->k2,chp);
		}
		else {		/* Can't find ca comp for kca. */
				/* Use default cai */
     		  kcatab((pnt->v-vext(kcapnt)-ccavoff(kcapnt)-kcapnt->voffsm),    /* set K rate chtyp */ 
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
  		conduct *= kcapnt->maxcond * get_gfrac(kcapnt,pnt->v);
		kcapnt->conduct = conduct;
                nodcr += conduct * (get_gvrev(kcapnt,pnt->v) - compcavoff(pnt));
                nodcur -= conduct * (pnt->v - vext(conpnt));  			/* current due to Vm */
                tcondn += conduct;    			/* total cond */
		break;

		case 1:					/* markov KCa chans */
		case 3:
		case 4:
		case 5:
		case 6:
		case 7:
		case 8:
		default:
	        kcapnts = (kcachan *)conpnt;
		if (pnt->capnt) {  /* calcium compartment exists */
		  if (!kcapnts->initfl){/* Calcium not set up before K chan. */
		     dochani((chan *)kcapnts,1e-9);
		     kcapnts->initfl = 1;
		  }
		}
  		conduct = dochan2(kcapnts);
  		conduct *= kcapnts->maxcond * get_gfrac(kcapnts,pnt->v);
		kcapnts->conduct = conduct;
                nodcr += conduct * (get_gvrev(kcapnts,pnt->v) - compcavoff(pnt));
                nodcur -= conduct * (pnt->v - vext(conpnt));  		/* current due to Vm */
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
		case 5:
		case 6:
		case 7:
		conduct = chpnt->conduct; /* already has gfrac */
                nodcr += conduct * (chpnt->gvrev - compcavoff(pnt));
                nodcur -= conduct * (pnt->v - vext(conpnt));   		/* current due to Vm */
                tcondn += conduct;	/* total cond */
/* ncfprintf (stderr,"conduct %g v %g vext %g gvrev %g\n",conduct,pnt->v,pnt->vext,chpnt->gvrev); /* */
/* ncfprintf (stderr,"conduct %g ipump %g\n",conduct,capnt->ipump); /* */
		break;

	    default:
	 	chpnt = (chan *)conpnt;
		conduct = chpnt->conduct * get_gfrac(chpnt,pnt->v); 
                nodcr += conduct * (get_gfrac(chpnt,pnt->v) - compcavoff(pnt));
                nodcur -= conduct * (pnt->v - vext(conpnt)); 		/* current due to Vm */
                tcondn += conduct;	/* total cond */
		break;

		}     /* switch (stype) */
		break;

      default:  break;

      } /* switch */
    }      /* for (lpnt=; ;) */

   if (pnt->jnoise) {				/* Johnson noise */
         static double tcondm;			/* tot membrane cond */

  	if (pnt->jstate) drand_setstate (pnt->jstate); /* if set individually */
        tcondm = tcond + tcondn + pnt->tcondm;      /* total membrane conductance */
	nodcr += gasdev() * ktrb * sqrt(tcondm) * pnt->jnoise;
                        /* See Hille 1992 p 325 on Johnson noise */
	restorstate();
   }
   tcond  += pnt->tcond;              	       /* total conductance */
   nodcur -= pnt->v * tcond;		       /* current due to comp volts through resistors */
   //nodcur -= (pnt->v-pnt->vext) * tcondn;      /* current due to comp volts through chans */
   //						/* now distributed to indiv chans above */
   //						/* to allow different vext for each chan */
   nodcr  += (pnt->vrev-compcavoff(pnt)) * pnt->rm;  /* membrane conductance */

/*  ncfprintf (stderr,"pnt->vrev %10.3g\n",pnt->vrev); /* */
/*  ncfprintf (stderr,"nodcr %10.3g nodcur %10.3g pntcond %g tcond %g\n",
			nodcr,nodcur,pnt->tcond,tcond); /* */

   // if (pnt->miscfl & IEXT)                    /* comp has a current src */
   nodcr += pnt->exti;                           /* external current */
   nodcur += nodcr;			      /* total comp current */

 if (euler)	{     			      /* forward-Euler */
     vest = pnt->v + nodcur * pnt->k;	      /* simple first-order est */
     if (pnt->miscfl & VEXT) {                /* comp is a battery */
       // pnt->extvi = (vest - pnt->extv) / (pnt->k);
       pnt->extvi = -nodcur;
       vest = pnt->extv;		      /*  battery voltage  */
     }
     else {
       pnt->extvi = 0.0;
     }
     pnt->vest = vest;				/* just save the estimate */
 }						/* (not quite) end of Euler */
						/*   (see below) */

 else {						/* Implicit modes */

   pnt->tcondt = tcond + tcondn;		/* total cond incl. new */
 //pnt->implf = 1. / (1. + pnt->tcondt * pnt->k * implfk); /* implicit factor   */

   /* Crank-Nicolson,  or backward Euler fully implicit */

    pnt->nodc = nodcr;				/* comp current precalc */
    pnt->nodcur = nodcur;			/* comp current recalc in docomp() */

    nodcur *=  pnt->k;                      	/* first-ord only */
    pnt->vest = pnt->v + nodcur;               /* simple first-order est */
    if (pnt->miscfl & VEXT) {                   /* comp is a battery */
         //pnt->extvi = (pnt->vest - pnt->extv) / (pnt->k);
         pnt->extvi = -pnt->nodcur;
         pnt->vest = pnt->extv;
      }
    else {
      pnt->extvi = 0.0;
    }
    pnt->verr = 1.;

  }  /* CN or implicit */

  pnt->vesto = pnt->vest;		/* save initial voltage for err est */

 }   /* for (pnt=;;) */

 if (euler) {					/* after everything else, */
   for (pnt=compnt; pnt; pnt=pnt->next) {
	pnt->v = pnt->vest;			/* save as new voltage */
   }
   return;
 }

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 relmult = relincr * timinc / 1e-4;
 mmaxerr=maxerr = 0.0;                 /* set maxerr so loop runs */ 
 ncompl = ncompr = tcomp = 0; 
 delcrit =   SQRT10; /* */
 /* delcrit =   10; /* */

/*-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

  /* always do one time through */

  for (pnt=compend; pnt; pnt=pnt->last)      /* left estimate */
   {
       docomp(pnt);
#ifdef DEBUG
       ncompl++;
#endif
   } 

  for (pnt=compnt; pnt; pnt=pnt->next)      /* right estimate */
   {
       docomp(pnt);
#ifdef DEBUG
       ncompr++;
#endif
   } 

  for (pnt=compnt; pnt; pnt=pnt->next)
   {
    pnt->g = pnt->t = 0;			/* start by resetting greater and total */
    doverr(pnt,noskip,critxx);			/* compute error estimate */
   }

  maxerr = mmaxerr = 1;

#ifdef DEBUG
  tcomp += ncompl + ncompr;
#endif

#ifdef DEBUG
  if (debug & NCOMP && debugz & 2)
   ncfprintf 
    (stderr,"ncomp before iter %-4d %-4d mmaxerr %-8.3g maxerr %-8.3g\n",
				ncompl,ncompr,mmaxerr,maxerr);
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

if (maxerr>crit) noskip = 1;	/* if needs more convergence */

for (titer=1, critxx=1, critx = min (maxerr*0.5, crit);
		 maxerr>crit;
		  titer+=niter) /* */
 {
 critx = min(critx,(mmaxerr+critx)/2),
 critx /= delcrit;
 critxx = critx;
 oldmaxerr=1.0;
 mmaxerr = 0.0;
 
 for (niter=0; maxerr>critx && maxerr>crit; niter++)	/* Iterate until convergence */
  {
   if (timinc > 1.0) runsyn(0);		/* if static, calc synapses */

/* ncfprintf (stderr,"maxerr %g critx %g critxx %g\n",maxerr,critx,critxx); /* */

   maxerr = 0.0;
   ncompl = ncompr = 0;

/*-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

   for (pnt=compend; pnt; pnt=pnt->last)     /* left estimate */
   {

// ncfprintf (stderr,"verr %g critxx %g ncompl %d\n",pnt->verr,critxx,ncompl);

#ifdef DEBUG
  if (debug & NCOMP && debugz & 32) {
     if (pnt->verr >= critxx) 
      ncfprintf (stderr,"c %-4d v %-10.6g relax %-8.4g g %-4d t %-4d misc %d\n",
			   pnt->num,pnt->v,pnt->relax,pnt->g,pnt->t,pnt->miscfl); /* */
   }
#endif
#ifdef DEBUG
       if (!(debug&NCOMP && debugz & NOSKIP))
#endif
       if (noskip==0) {
         if (pnt->verr<critxx) continue; /* */
       }
       docomp(pnt);
#ifdef DEBUG
       ncompl++;
#endif
   }

/*-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */
  //if (niter>0 || titer>1)
   for (pnt=compnt; pnt; pnt=pnt->next)      /* right estimate */
   {
#ifdef DEBUG
       if (!(debug&NCOMP && debugz & NOSKIP))
#endif
       if (noskip==0)
         if (pnt->verr<critxx) continue; /* */
       docomp(pnt);
#ifdef DEBUG
       ncompr++;
#endif
   }           /* for (ncomp=0,pnt=compnt;) */

/*-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

   for (pnt=compnt; pnt; pnt=pnt->next)
    {
     doverr(pnt,noskip,critxx);			/* compute error estimate */
    }

/*-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

#ifdef DEBUG
     tcomp += ncompl + ncompr;
#endif

   /* Reduce the criterion if error increases, and also at least
       once every 10 iterations.  This causes compartments
       outside the region of iteration (i.e. otherwise skipped over)
       to be recomputed and reduces the total number of iterations
       necessary for complete convergence.
   */
       
   if ((maxerr>(oldmaxerr*1.5) && oldmaxerr!=0 && critxx > crit*0.5)) {
    noskip += 1; 
    if (noskip > 100) noskip = 100;
    critxx /= delcrit/2;

#ifdef DEBUG
  if (debug & NCOMP && debugz & 16)
 ncfprintf (stderr,"maxerr   %8.3g, oldmaxerr %8.3g noskip %d reducing critxx to %8.3g\n",
					maxerr,oldmaxerr,noskip,critxx);
#endif
    oldmaxerr = 1.0; 

   }  /* if (maxerr>oldmaxerr || ) */

   else if ((niter>0) && (niter%20==0)) {
     noskip++;
#ifdef DEBUG
     if (debug & NCOMP && debugz & 16)
	 ncfprintf (stderr,"maxerrok %8.3g, noskip %d\n",
					maxerr,noskip);
#endif
     if ((niter%40)==0) {
       critxx /= delcrit;

#ifdef DEBUG
  if (debug & NCOMP && debugz & 16)
	 ncfprintf (stderr,"maxerrok %8.3g, noskip %d reducing critxx to %8.3g\n",
					maxerr,noskip,critxx);
#endif
     }
   }

   else {
      oldmaxerr = maxerr; 
      if (noskip>0) noskip--;		/* allow skip after first iteration */
   }

   if (maxerr==0) noskip++;		/* don't skip when no comps done */

   if (critxx < critx*1e-4) { 		/* don't skip when unstable */
	noskip = 10; 	
   }
   if (noskip > 0) {
	oimplfk = implfk;
	implfk = 1;
   }
   else {
	implfk = oimplfk;
   }
   if (mmaxerr<maxerr) mmaxerr = maxerr;   /* remember largest maxerr */

#ifdef DEBUG
  if (debug & NCOMP && debugz & 16)
   ncfprintf 
    (stderr,"ncomp %-4d %-4d mmaxerr %-8.3g maxerr %-8.3g critx %-8.3g critxx %8.3g\n",
				ncompl,ncompr,mmaxerr,maxerr,critx,critxx);
#endif

   if (critxx < 1e-30) {
    ncfprintf (stderr,"# Ncomp: error, no convergence, continuing, pid %d\n",getpid());
      return;
   }
 
 }            /* for (niter=0; maxerr>critx; niter++)  */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

/* Reset error after each level of critx */
/*  This has the effect of spreading the region of iteration,
    and causes compartments just outside to be recomputed.
*/

 critk = .9 * critx / delcrit;
 for (pnt=compnt; pnt; pnt=pnt->next)
  {
   pnt->verr += critk;
  }			/* */

#ifdef DEBUG
  if (debug & NCOMP && debugz & 8)
   ncfprintf 
(stderr,"niter %2d        mmaxerr %-8.3g maxerr %-8.3g critx %-8.3g\n",
				niter,mmaxerr,maxerr,critx);
#endif

  }         /* for (critx ; ; ) */

/* Just save it and calculate relax for next time step. */

#ifdef DEBUG
 avgrelax = 0;
 nrelax = 0; 
#endif
 for (pnt=compnt; pnt; pnt=pnt->next)
  {
   pnt->v = pnt->vest;		/* save it */
#ifdef DEBUG
   avgrelax += pnt->relax;
   nrelax++;
#endif
  }
#ifdef DEBUG
  if (nrelax==0) nrelax = 1;
  avgrelax /= nrelax;
#endif

/* This code runs only once per time step, */
/* so it adds negligible time to computation */
/* if model is large and requires many iterations */
/* per time step. */

if (relincr>0 && titer>10) {			/* if lots of iterations */
 for (pnt=compnt; pnt; pnt=pnt->next) {
   if (!(pnt->miscfl & VEXT)) {			/* if comp is not a battery */
      less = pnt->t - pnt->g;			/* then compute oscillations */
      m = ((pnt->g < less) ? less : pnt->g);	/* max of less, greater */
      if (m==pnt->t) pnt->relax += relmult;	/* adjust comp relax */
      else           pnt->relax -= (1.0 - ((double) m)/pnt->t) * 10.0 * relmult;
      if      (pnt->relax > 1.0) pnt->relax = 1.0;
      else if (pnt->relax < 0.0) pnt->relax = 0.0;
    }
  }
 }

#ifdef DEBUG
  if (debug & NCOMP && debugz & 2) ncfprintf (stderr,"time %-8.6g ", simtime);
  if (debug & NCOMP && debugz & 2 && !(debugz & 4)) ncfprintf (stderr,"\n");

  if (debug & NCOMP && debugz & 4)
	ncfprintf (stderr,"tcomp %4d titer %d avgrelax %-9.4g\n",
			tcomp,titer,avgrelax);
#endif


#ifdef DEBUG
  if (debug & NCOMP && debugz & 1) ncfprintf (stderr,"runcomp end.\n");
#endif
}

/*------------------------------------*/

void docomp(register comp *pnt)
{
  register conlst *lpnt;
  register conn *conpnt;
  static cacomp *capnt;                 /* pntr to calcium compartment */
  static double nodcur;                 /* total comp current */
  static double conduct;                /* channel conductance */
  static double implf;			/* new implicit factor */
  static comp *ocpnt;
  static double icomp;
  static hhchan *hhpnt;
  static double m,n,h,k,dm,dn,dh;
  int i;

     /* start with old node cur */
     /* Proportion of nodc and nodcur sets relative amount of implicit. */
     /* implfk = 0.5 -> C-N, implfk = 1 -> fully implicit. */
     /*  Can change implfk during interation if necessary */

     nodcur = (implfk*pnt->nodc + (1-implfk)*pnt->nodcur)/implfk;  

     if (capnt=pnt->capnt) {			/* exchanger current. */
	 static double eica;
       if (capnt->kexo) {
	     double eica;
 	 eica = capnt->slope*(pnt->v*implfk + pnt->vest*(1.-implfk))+capnt->eica;
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
         case PNX:
                nodcur += ocpnt->vest * conpnt->conduct;  /* curr thru resist */
                break;
	 case DIODE:
		if ((pnt->vest - ocpnt->vest) < 0)
                   nodcur += ocpnt->vest * conpnt->conduct;  /* curr thru resist */
                break;
         case BATT:
     		if (ocpnt == conpnt->comp1) 
                  nodcur += (ocpnt->vest - conpnt->conduct) * BATTCOND;
		else
                  nodcur += (ocpnt->vest + conpnt->conduct) * BATTCOND;
                break;		     /* curr thru batt small resistor */

         case BUF: if (pnt == ((dbuf *)conpnt)->comp2) {
                     if (((dbuf *)conpnt)->delay || ((dbuf *)conpnt)->filt) {
		          pnt->extv = ((dbuf*)conpnt)->v;
		     }
		     else pnt->extv = (((dbuf *)conpnt)->comp1->vest - ((dbuf *)conpnt)->offset) *
		     		       ((dbuf *)conpnt)->gain;
		   }
		break;

         // case NBUF: if (pnt == ((ndbuf *)conpnt)->comp2) {
	 //	    setnt (((ndbuf *)conpnt)->comp2, ((ndbuf *)conpnt)->ntrans, ((ndbuf *)conpnt)->comp1->v);
 	 //	 }
	 // 	break;

         case CAP:
		nodcur += ocpnt->vest * conpnt->conduct / timinc;	
                break; 
         case NA:
      	 case K:   
      	 case ClCa:   
      	 case KCa:   
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

     k = pnt->k * implfk;
     implf = 1. / (1. + pnt->tcondt * k);  	/* implicit factor */
     pnt->vest = (pnt->v + nodcur * k) * implf;
     if (pnt->miscfl & VEXT) {			/* comp is a battery */
	  // icomp = (pnt->vest - pnt->extv) / (k * implf);
          // pnt->extvi = -icomp;
//fprintf (stderr,"icomp %g nodcur %g\n",icomp,nodcur*implf);
        pnt->vest = pnt->extv;
     }
}

/*------------------------------------*/

void doverr(register comp *pnt, int noskip, double critxx) 
{
   static double err;

  if (pnt->miscfl & VEXT) {		/* comp is a battery */
     pnt->verr = 1.0;			/* set verr so vclamp always runs */
  }
  else {
    err = pnt->vesto - pnt->vest;   /* find convergence deriv */
    pnt->vest -= err * pnt->relax;
    if (noskip || pnt->verr>=critxx) {
      pnt->t++;
      if (err>=0.0) pnt->g++;
    }
    else err = -err;
    pnt->verr = err;
    if (err > maxerr) maxerr = err;
  }
  pnt->vesto = pnt->vest;   		/* save volt est for convergence deriv */
}

/*------------------------------------*/

