/* segment chanfunc in program nc */

/* routines to set up channel parameters */

#include <stdio.h>
#include <string.h>
#include "nc.h"
#include "y.tab.h"
#include "nconst.h"
#include "ncsub.h"
#include "ncomp.h"
#include "control.h"
#include "ncelem.h"
#include "ndef.h"

static Symbol *fpnt;

static int nchantype = 0;
static chantype chtypes[NCHANCONST];
chantype *findchantype(int ctype, int stype);
chantype *getchantype(int ctype, int stype);
void initchan(void);

// extern FILE *stdout;
// extern FILE *stderr;
extern int ncfprintf(FILE *stream, const char *format, ...);
extern int interp;			/* running interpeter */

extern iontab *ions;		// ion concentrations and valence in init.cc
extern int     ionz[]; 		// ion valence in init.cc

#ifdef __cplusplus
extern "C" {
#endif

double exp(double);
double log(double);
double log10(double);

#ifdef __cplusplus
}
#endif

/* #define DEBUG        /* */

double ncabs(double x);
double callfunc(Symbol *funcp, int npar, double par1, double par2,
					 double par3, double par4);
char *emalloc(unsigned int n);
char *findsym(int);
elem *findelem (int num);
void execerror (const char *s, const char *t);
void varcopyu(void);
attrib *getattrib(int elnum);
Symbol *lookup(const char *s);
double bkcainf (double v, double ca, mattrib *mpnt, double arate, double brate);
double bkcatau (double v, double ca, mattrib *mpnt, double arate, double brate, 
		chanparm *chp);
double vext(conn *ch);
void makghk_tab (chantype *ch);

double (rt)(chan *c);
double (rnt)(chan *c);
double (rca)(chan *c);
double (al0)(chan *c);
double (bet0)(chan *c);
double (al1)(chan *c);
double (bet1)(chan *c);

iontab *init_ions(iontab *itab);		// in init.cc
double nernstc (double vrev, double co, int z);
double ghkv (chantype *chtyp, iontab *ions, double cao, double cai);
double ghkv (chantype *chtyp, double cao, double cai);
double ghkv (chan *ch, double cao, double cai);

/*----------------------------------------*/

chantype *makchantype(int ctype, int stype, int nstates, int nparm, 
				int nq, double basetempc)

/* Information describing the state variables of a HH-type channel
  (e.g. m or h, etc.) or a Markov sequential-state channel.

  Contains table to hold values of alpha and beta calculated by
  the same function call, and with the same Q10, table to hold
  value of m1, m2 for implicit integration of the variable, pointer
  to default function to calculate rates, name of user-defined
  function to override calculation of rates, and storage space for
  the Q10 value, v, timinc, and tempcel variables used to
  determine at run time if table calculation is still valid.
*/

{
   chantype *ch;
   chanparm *chp;
   int i, t;

/*
  if (!(ch=(chanstate *)emalloc(sizeof(chantype)))) {
     ncfprintf (stderr,"makchantype: can't allocate table of size %d\n",nstates);
     return ((chantype*)NULL);
  }
*/

  if  (!(ch=findchantype(ctype,stype))) {
      ch = &chtypes[nchantype++];		/* allocate new chan type */
  }
  else {
     ncfprintf (stderr,"makchantype: chan %s type %d is already defined.\n",
		findsym(ctype),stype);
  }
  if (nchantype > NCHANCONST) {
	nchantype = NCHANCONST-1;
	ncfprintf (stderr,"makchantype: too many chan types\n");
  }
 
  ch->ctype = ctype;
  ch->stype  = stype;

  ch->numstate = nstates;
  ch->numparm = nparm;		 
  ch->hh = 0;		 	/* = 1 -> HH implicit tables used for rates */
  ch->qc = dqc;			/* Q10 factor for conductance */
  ch->qt = basetempc;		/* base temperature for channel rate & cond */
  ch->qcond = 0;		/* Q10 for unitary conductance */
  ch->vrev = NULLVAL;
  ch->unitary = 22e-12;
  ch->trconc = dstr;
  ch->state = (chanstate *)NULL;
  ch->parm = (chanparm *)NULL;
  ch->respamp = (double *)NULL;
  ch->cabnd   = 0;		/* needs Ca binding */
  ch->cgbnd   = 0;		/* needs cGMP binding and saturation */
  ch->gamma   = 0;		/* weight for O2 cond in ChR2 */
  ch->gfractab = (double *)NULL;
  ch->gvrevtab = (double *)NULL;
  ch->ions = (iontab *)NULL;

  if (nstates) {		/* allocate space for state definitions */
    if (nstates > NUMSTATE) {
	ncfprintf (stderr,"makchantype: %d too many states\n",nstates);
        nstates = NUMSTATE;
    }
    if (!(ch->state=(chanstate *)emalloc(nstates*sizeof(chanstate)))) {
      ncfprintf(stderr,"makchantype: can't allocate table of size %d\n",nstates);
       return ((chantype*)NULL);
    }
  }
				/* allocate space for param table */
  if (!(ch->parm=(chanparm *)emalloc(nparm*sizeof(chanparm)))) {
    ncfprintf (stderr,"Makechantype: can't allocate table of size %d\n",nparm);
    return ((chantype*)NULL);
  }
  for (i=0; i<nparm; i++) {
    chp = &ch->parm[i];
    chp->pn  = i;
    chp->nfval  = 0;
    chp->nival  = 0;
    chp->fval  = (double *)NULL;
    chp->ival  = (double *)NULL;
    chp->chanftab = (double *)NULL;
    chp->chanitab = (double *)NULL;
    chp->chancalc = (double (*)(double,int))NULL;
    chp->funcname = (char *)NULL;
    chp->nq = nq;			/* number of Q10 values */
    if (!(chp->qrate=(double *)emalloc(nq*sizeof(double)))) {
      ncfprintf (stderr,"Makechantype: can't alloc qrate table of size %d\n",nq);
      return ((chantype*)NULL);
    }
    if (!(chp->dq=(double *)emalloc(nq*sizeof(double)))) {
      ncfprintf (stderr,"Makechantype: can't alloc dq table of size %d\n",nq);
      return ((chantype*)NULL);
    }
    for (t=0; t<nq; t++) {
      chp->qrate[t] = 1.0;	/* Q val calculated from qt, dq */
      chp->dq[t] = 1.0;		/* Q10 temp factor for param alpha */
    }
    chp->qt = ch->qt;		/* base temperature for rate calc */
    chp->dqca = dcavoff;	/* Q10 ca voffset factor for param */
    chp->qcao = DCAO + DMGO;	/* base cao, mgo for voffset calc */
    chp->qcavoff = 0;		/* voffset calc from dqca, dcao,dcoo,dmgo */
    chp->voff = 0;		/* voffset from user */
    chp->oldvf = -1e10;
    chp->oldvi = -1e10;
    chp->oldar = -1e10;
    chp->oldbr = -1e10;
    chp->oldtimincf = 0;
    chp->oldtiminci = 0;
    chp->oldtempcel = 0;
    chp->oldcao = 0;
  }
  if (ctype >= GLU && ctype <= CGMP) {	   /* if receptor binds ligand */
    if (!(ch->respamp=(double *)emalloc(NNTR*sizeof(double)))) {
      ncfprintf(stderr,"makchantype: can't allocate respamp of size %d\n",NNTR);
      return ((chantype*)NULL);
    }
    for (i=0; i<NNTR; i++) {
      ch->respamp[i] = 0;
    }
  }
  ch->ions = init_ions(ch->ions);
  for (i=0; i<NIONS; i++) {
     ch->ions->ioni[i] = ions->ioni[i];
     ch->ions->iono[i] = ions->iono[i];
     ch->ions->ionp[i] = 0;
  }
  return ch;
}

/*--------------------------------------------*/

chantype *copychantyp(int sctype, int sstype, int dctype,
					int dstype, int nstates)

/* Copy a channel type to allow creating a new type. */

{
     int i, j, n, t;
     chantype *chs, *chd;	/* source, destination channel types */
     chanparm *chps,*chpd;	/* source, destination params */

  if  (!(chs=getchantype(sctype,sstype))) {
      ncfprintf (stderr,"copychantype: can't find original type %d, %d\n",
					sctype,sstype);
      return NULL;
  }
  if ((sctype == dctype) && (sstype == dstype)) {
      ncfprintf 
(stderr,"copychantype: can't make same type: '%s' type %d = '%s' type %d\n",
					findsym(sctype),sstype,
					findsym(dctype),dstype);
      return NULL;
  }
  if (nstates < 2) nstates = chs->numstate;
  if ((chd=makchantype(dctype,dstype,nstates,chs->numparm,
		chs->parm[0].nq,chs->qt))==NULL) {
    ncfprintf (stderr,"copychantype: can't make new channel type %d %d\n", 
		dctype, dstype);
    return NULL;
  }
  chd->unitary = chs->unitary;		/* unitary current */
  chd->trconc =  chs->trconc;		/* max conc of nt */
  chd->vrev   = chs->vrev;		/* reversal potential */
  chd->cabnd  = chs->cabnd;		/* needs Ca binding */
  chd->cgbnd  = chs->cgbnd;		/* needs cGMP binding */

  for (i=0; i<chs->numparm; i++) {	/* copy parameter tables */
    chps = &chs->parm[i];
    chpd = &chd->parm[i]; 
    chpd->pn = chps->pn;
    chpd->nfval = chps->nfval;
    chpd->nival = chps->nival;
    chpd->nq = chps->nq;
    chpd->qt = chps->qt;
    chpd->dqca = chps->dqca;
    chpd->qcao = chps->qcao;
    chpd->qcavoff = chps->qcavoff;
    chpd->voff = chps->voff;
    chpd->oldvf = chps->oldvf;
    chpd->oldvi = chps->oldvi;
    chpd->oldar = chps->oldar;
    chpd->oldbr = chps->oldbr;
    chpd->oldtimincf = chps->oldtimincf;
    chpd->oldtiminci = chps->oldtiminci;
    chpd->oldtempcel = chps->oldtempcel;
    chpd->oldcao = chps->oldcao;
    if (chps->fval)
      for (j=0; j<chps->nfval; j++) {	/* copy values */
        chpd->fval[j] = chps->fval[j];
      }
    if (chps->ival)
      for (j=0; j<chps->nival; j++) {	/* copy values */
        chpd->ival[j] = chps->ival[j];
      }
    if (chps->chanftab) {
      n = chps->nfval * (RATESIZ+1);	/* copy rate tables */
      for (j=0; j<n; j++) {
        chpd->chanftab[j] = chps->chanftab[j];
      }
    }
    if (chps->chanitab) {
      n = chps->nival * (RATESIZ+1);	/* copy rate tables */
      for (j=0; j<n; j++) {
        chpd->chanitab[j] = chps->chanitab[j];
      }
    }
    chpd->chancalc = chps->chancalc;	/* copy func */
    n = strlen(chps->funcname);		/* copy funcname */
    chpd->funcname = (char*)emalloc(n);
    for (j=0; j<n; j++) {
      chpd->funcname[j] = chps->funcname[j];
    }
    for (t=0; t<chps->nq; t++) {	/* copy temp rates */
      chpd->qrate[t] = chps->qrate[t];
      chpd->dq[t]    = chps->dq[t];
    }
  }
  for (i=0; i<nstates; i++) {		/* copy states */
    chd->state[i] = chs->state[i];
  } 
  for (i=0; i<NIONS; i++) {		/* copy permeabilities */
     chd->ions->ioni[i] = chs->ions->ioni[i];
     chd->ions->iono[i] = chs->ions->iono[i];
     chd->ions->ionp[i] = chs->ions->ionp[i];
  }
  for (i=0; i<NNTR; i++) {		/* copy resp amplitudes */
    chd->respamp[i] = chs->respamp[i];
  }
  return chd;
}

/*--------------------------------------------*/

datum copychan(datum &sctype, datum &sstype, datum &dctype,
					datum &dstype, datum &nstates)

/* make a new channel type by copying an old one */

{
  int sctyp, dctyp;

  if (!(sctyp=lookup(sctype.str)->type))
    ncfprintf (stderr,"copychan: can't find %s\n",sctype.str);
  if (!(dctyp=lookup(dctype.str)->type))
    ncfprintf (stderr,"copychan: can't find %s\n",dctype.str);

   copychantyp(sctyp, (int)sstype.val, dctyp, (int)dstype.val, 
						(int)nstates.val);
   return dstype;
}

/*--------------------------------------------*/

datum setchan_ntrans (datum &ctype, datum &stype, datum &state, datum &val)

/* Set "cond" for a channel state. */

{
   chantype *ch;
   double tval;

  if  (!(ch=getchantype(lookup(ctype.str)->type,(int)stype.val))) {
      ncfprintf (stderr,"setchan_ntrans: can't find chan type %s, %d\n",
					ctype.str,(int)stype.val);
      return ctype;
  }
  if (state.val < 0) state.val = 0;
  if (state.val >= ch->numstate) state.val = ch->numstate-1;
  tval = ch->state[int(state.val)].numtrans;
  ch->state[int(state.val)].numtrans = int(val.val);
  ctype.val = tval; 
  return ctype;
}

/*--------------------------------------------*/

datum setchan_cond (datum &ctype, datum &stype, datum &state, datum &val)

/* Set "cond" for a channel state. */

{
   chantype *ch;
   double tval;

  if  (!(ch=getchantype(lookup(ctype.str)->type,(int)stype.val))) {
      ncfprintf (stderr,"setchan_cond: can't find chan type %s, %d\n",
					ctype.str,(int)stype.val);
      return ctype;
  }
  if (state.val < 0) state.val = 0;
  if (state.val >= ch->numstate) state.val = ch->numstate-1;
  tval = ch->state[int(state.val)].cond;
  ch->state[int(state.val)].cond = val.val;
  ctype.val = tval; 
  return ctype;
}

/*--------------------------------------------*/

datum setchan_trans (datum &ctype, datum &stype, 
				datum &state, datum &trans, datum &val)

/* Set "trans" for a channel state transition, */
/* i.e. the state that this transition goes to. */
{
   chantype *ch;
   double tval;

  if  (!(ch=getchantype(lookup(ctype.str)->type,(int)stype.val))) {
      ncfprintf (stderr,"setchan_trans: can't find chan type %s, %d\n",
					ctype.str,(int)stype.val);
      return ctype;
  }
  if (trans.val < 0) trans.val = 0;
  if (trans.val >= NUMTRANS) trans.val = NUMTRANS-1;
  if (state.val < 0) state.val = 0;
  if (state.val >= ch->numstate) state.val = ch->numstate-1;
  if (trans.val < 0) trans.val = 0;
  if (trans.val >= NUMTRANS) trans.val = NUMTRANS-1;
  tval = ch->state[int(state.val)].trans[int(trans.val)];
  ch->state[int(state.val)].trans[int(trans.val)] = int(val.val);
  ctype.val = tval; 
  ctype.vtype = NUMBER; 
  return ctype;
}

/*--------------------------------------------*/

datum setchan_trate (datum &ctype, datum &stype, 
				datum &state, datum &trans, datum &val)

/* Set "trate" for a channel state transition, */
/* i.e. the rate function. */

{
   chantype *ch;
   double (*f)(chan *c), (*oldf)(chan *c);

  if  (!(ch=getchantype(lookup(ctype.str)->type,(int)stype.val))) {
      ncfprintf (stderr,"setchanmul: can't find chan type %s, %d\n",
					ctype.str,(int)stype.val);
      return ctype;
  }
  if (trans.val < 0) trans.val = 0;
  if (trans.val >= NUMTRANS) trans.val = NUMTRANS-1;
  if (state.val < 0) state.val = 0;
  if (state.val >= ch->numstate) state.val = ch->numstate-1;
  if (trans.val < 0) trans.val = 0;
  if (trans.val >= NUMTRANS) trans.val = NUMTRANS-1;
  if     (strcmp(val.str,"rt"))     f = rt;
  else if(strcmp(val.str,"rnt"))    f = rnt;
  else if(strcmp(val.str,"rca"))    f = rca;
  else if(strcmp(val.str,"alpha0")) f = al0;    /* alpha, beta for parm 0 */
  else if(strcmp(val.str,"beta0"))  f = bet0;
  else if(strcmp(val.str,"alpha1")) f = al1;    /* alpha, beta for parm 1 */
  else if(strcmp(val.str,"beta1"))  f = bet1;
  else f = NULL;
  oldf = ch->state[int(state.val)].trate[int(trans.val)];
  ch->state[int(state.val)].trate[int(trans.val)] = f;
  ctype.str = (char*)oldf; 
  ctype.vtype = NUMBER; 
  return ctype;
}

/*--------------------------------------------*/

datum setchan_mul (datum &ctype, datum &stype, 
				datum &state, datum &trans, datum &val)

/* Set "ratemul" for a channel state transition, */
/* i.e. the rate multiplier. */

{
   chantype *ch;
   double tval;

  if  (!(ch=getchantype(lookup(ctype.str)->type,(int)stype.val))) {
      ncfprintf (stderr,"setchanmul: can't find chan type %s, %d\n",
					ctype.str,(int)stype.val);
      return ctype;
  }
  if (trans.val < 0) trans.val = 0;
  if (trans.val >= NUMTRANS) trans.val = NUMTRANS-1;
  if (state.val < 0) state.val = 0;
  if (state.val >= ch->numstate) state.val = ch->numstate-1;
  if (trans.val < 0) trans.val = 0;
  if (trans.val >= NUMTRANS) trans.val = NUMTRANS-1;
  tval = ch->state[int(state.val)].ratemul[int(trans.val)];
  ch->state[int(state.val)].ratemul[int(trans.val)] = val.val;
  ctype.val = tval; 
  ctype.vtype = NUMBER; 
  return ctype;
}

/*--------------------------------------------*/

datum setchan_rateo (datum &ctype, datum &stype, 
				datum &state, datum &trans, datum &val)

/* Set "rateo" for a channel state transition, */
/* i.e. the user-specified multiplier for the rate function */
{
   chantype *ch;
   double tval;

  if  (!(ch=getchantype(lookup(ctype.str)->type,(int)stype.val))) {
      ncfprintf (stderr,"setchanrateo: can't find chan type %s, %d\n",
					ctype.str,(int)stype.val);
      return ctype;
  }
  if (trans.val < 0) trans.val = 0;
  if (trans.val >= NUMTRANS) trans.val = NUMTRANS-1;
  if (state.val < 0) state.val = 0;
  if (state.val >= ch->numstate) state.val = ch->numstate-1;
  if (trans.val < 0) trans.val = 0;
  if (trans.val >= NUMTRANS) trans.val = NUMTRANS-1;
  tval = ch->state[int(state.val)].rateo[int(trans.val)];
  ch->state[int(state.val)].rateo[int(trans.val)] = int(val.val);
  ctype.val = tval; 
  ctype.vtype = NUMBER; 
  return ctype;
}

/*--------------------------------------------*/

double qrate (chanparm *chp)

/* Compute Q10 multiplier for alpha for channel rate functions */

/* Return rate function multiplier which depends on temperature.
   Base for HH rate functions = 6.3 deg C.

   Base temp set to 22 deg because many channels are now studied
   at room temperature.

  Assumes that timinc and/or tempcel may change during a run,
    but that the Q10 value and base temperature do not.
    Function of time step.

     chp->dq = Q10 value
     chp->qt = base temperature for calculation
*/

{
   int i;
   double dq;

  if (tempcel != chp->oldtempcel) {

    for (i=0; i<chp->nq; i++) {
      if ((dq=chp->dq[i]) <= 0) dq = 2;
      chp->qrate[i] = exp(log(dq) * (tempcel - chp->qt) / 10.0);
    }
    chp->oldtempcel = tempcel;
  }
  return (chp->qrate[0] * timinc);
}

/*--------------------------------------------*/

double qraten (chanparm *chp, int t)

/* Compute Q10 multiplier for alpha for channel rate functions */

/* Return rate function multiplier which depends on temperature.
   Base for HH rate functions = 6.3 deg C.

   Base temp set to 22 deg because many channels are now studied
   at room temperature.

  Assumes that timinc and/or tempcel may change during a run,
    but that the Q10 value and base temperature do not.
    Function of time step.

     chp->dq = Q10 value
     chp->qt = base temperature for calculation

  The value of "nq" sets the behavior of the Q10 values.  If it
  is set to 1, then the same Q10 is used for all the rate
  functions (defined in dq[0].  If nq > 1, then each Q10 value
  must be set individually in the "dq" array.

*/

{
   int i, nq;
   double dq;

  if ((nq=chp->nq) < t) t = nq;

  if (tempcel != chp->oldtempcel) {

    for (i=0; i<nq; i++) {
      if ((dq=chp->dq[i]) <= 0) dq = 2;
      chp->qrate[i] = exp(log(dq) * (tempcel - chp->qt) / 10.0);
    }
    chp->oldtempcel = tempcel;
  }
  return (chp->qrate[t] * timinc);
}

/*--------------------------------------------*/

double qcond (chantype *ch)

/* Compute Q10 multiplier for channel unitary conductance */

/* Return conductance multiplier which depends on temperature.
   Original base for HH rate functions = 6.3 deg C.

   Base temp set to 22 deg because many channels are now studied
   at room temperature.

     ch->qc    = Q10 value
     ch->qt    = base temperature for calculation
     ch->qcond = multiplier calculated from qc, qt.
*/

{
  ch->qcond = exp(log(ch->qc) * (tempcel - ch->qt) / 10.0);
  return ch->qcond;
}

/*--------------------------------------------*/

double qcavoff (chanparm *chp)

/* Compute initial channel gating voltage offset due to surface charge from Cao */

/* Return voffset that depends on cao. */

{
  return log10( (dcao+dmgo+dcoo) / chp->qcao);
}

/*--------------------------------------------*/

double ccavoff (chan *ch)

/* Compute runtime channel gating voltage offset due to surface charge from Cao */

/* Return voffset which depends on cao. */

{
    double cao, cavoff;
   chanparm *chp;
   cacomp *capnt;

  chp = ch->chtyp->parm;
  if ((ch->comp1!=NULL) && (capnt=ch->comp1->capnt) != NULL) {
     if (capnt->caos!=NULL)
        cao = capnt->caos[0];
     else cao = capnt->cao;
  } else {
    cao = dcao;
  }
  cavoff = chp->dqca * log10( (cao+dmgo+dcoo) / chp->qcao);
  // if (ch->compp !=NULL) cavoff += ch->compp->v;		/* add pH effect on gating voltage */
  // fprintf (stderr,"cavoff %g\n",cavoff);
  return cavoff;
}

/*----------------------------------------*/

double compcavoff (comp *pnt)

/* Compute runtime membrane Vrev offset due to surface charge from Cao */

{
    double cao, cavoff;
    cacomp *capnt;

   if ((capnt=pnt->capnt)!=NULL) {
      if (capnt->caos!=NULL) 
         cao = capnt->caos[0]; 
      else cao = capnt->cao;
   } else {
      cao = dcao;
   }
   cavoff = dcaspvrev * log10( (cao+dmgo+dcoo) / (DCAO+DMGO));
   return cavoff;
}

/*----------------------------------------*/

chantype *findchantype(int ctype, int stype)

{
   int i,found;
   chantype *ch;

  for (found=i=0; !found && i<nchantype; i++) {
    if ((chtypes[i].ctype == ctype) && (chtypes[i].stype == stype)) {
      ch = &chtypes[i];
      found = 1;
      break;
    }
  }
  if (!found) ch = ((chantype *)NULL);
  return ch;
}

/*----------------------------------------*/

chantype *getchantype(int ctype, int stype)

{
     chantype *ch;
     int found;

  found = 1;
  if (!(ch=findchantype(ctype,stype))) {
     initchan();
     if (!(ch=findchantype(ctype,stype))) {
	found = 0;
	if (stype==0) {		/* default type not found, try type 1 */
	   found = 1;
           if (!(ch=findchantype(ctype,1))) found=0; 
	}
     }
  }
  if (!found) {
   ncfprintf(stderr,"## getchantype: don't know how to make chan '%s type %d'\n",
		 findsym(ctype),stype);
    ch  = ((chantype *)NULL); 
  }
  return ch;
}

/*----------------------------------------*/

int find_pri_ion (iontab *ions)

/* check ion permeabilities in local ions[] to find the primary ion */

#define IONTHR 0.2
{
    int i,pri_ion;
    double maxperm,perm;
 
  for (maxperm=i=0; i<NIONS; i++) {		// copy ion concs from chan type
    if ((perm=ions->ionp[i]) > maxperm) {
        maxperm = perm; 
        pri_ion = i;
    }
  } 
  return pri_ion;
}

/*----------------------------------------*/

int find_sec_ion (iontab *ions, int pri_ion)

/* check ion permeabilities in local ions[] to find the second ion */

#define IONTHR 0.2
{
    int i,sec_ion;
    double maxperm, perm;
    
  sec_ion = 1;
  for (maxperm=-1,i=0; i<NIONS; i++) {		// copy ion concs from chan type
    if (i==pri_ion) continue;			// skip primary ion
    if ((perm=ions->ionp[i]) > maxperm) {
        maxperm = perm; 
        sec_ion = i;
    }
  } 
  return (sec_ion);
}

/*----------------------------------------*/

double calc_cvfi (chan *ch, double v)

/* Calculate most of the GHK current equation for a channel */
/* See cadfxx.n series */

{
    int ion, chion;
    double chioni,chiono,chionm,perm;
    double vf,cvfi,z,zm,zmi;
    iontab *ions;


  if (ch==NULL) {
      ncfprintf(stderr,"## calc_cvfi: missing channel\n");
      return 0;
  }
  if ((ions=ch->ions)==NULL) ions = ch->chtyp->ions;

  chion = find_pri_ion(ions);
  chionm = ions->iono[chion];		// The ext conc of the channel's "main ion".
  				        // This normalizes the GHK curr eqn to that ion.
				 	// See "cadf_7.n"
  if (chionm < ions->ioni[chion]) chionm = ions->ioni[chion]; // set to whichever one is greatest

  zm = ionz[chion];
  zmi = zm*zm;				// Normalizes other chan types to Ca, see cadf_7.n
					// For Na, K, Cl ions, it mults their effect by 0.25 for Ca chans
  for (cvfi=ion=0; ion<NIONS; ion++) {
    z = ionz[ion];
    chioni = ions->ioni[ion];
    chiono = ions->iono[ion];
    perm   = ions->ionp[ion];
    vf = exp(-z*v*frt);
    cvfi += perm * z*z/zmi * (chioni-chiono*vf)/(1-vf)/chionm;
    // fprintf (stderr,"cvfi %g z %g zmi %g chioni %g chiono %g perm %g chion %d chionm %g\n",cvfi,z,zmi,chioni,chiono,perm,chion,chionm);
  }
  return cvfi;
}

/*----------------------------------------*/

double calc_gslope (chan *ch, double v)

/* Calculate the slope conductance for a channel */
/* See cadfxx.n series */

{
    int ion, chion;
    double chioni,chiono,chionm,perm;
    double vf,vfi,gfrac,z,zm,zmi;
    iontab *ions;

  if (ch==NULL) {
      ncfprintf(stderr,"## calc_cglsope: missing channel\n");
      return 0;
  }
  if ((ions=ch->ions)==NULL) ions = ch->chtyp->ions;

  chion = find_pri_ion(ions);
  chionm = ions->iono[chion];		// The ext conc of the channel's "main ion".
  				        // This normalizes the GHK curr eqn to that ion.
				 	// See "cadf_7.n"	
  if (chionm < ions->ioni[chion]) chionm = ions->ioni[chion];// set to whichever one is greatest
 
  zm = ionz[chion];
  zmi = zm*zm;				// Normalizes other chan types to Ca, see cadf_7.n
					// For Na, K, Cl ions, it mults their effect by 0.25 for Ca chans
  for (gfrac=ion=0; ion<NIONS; ion++) {
    z    = ionz[ion];
    chioni = ions->ioni[ion];
    chiono = ions->iono[ion];
    perm   = ions->ionp[ion];
    vf = exp(-z*v*frt);
    vfi = 1 / (1 - vf);
    gfrac += perm * frt * v *z*z*z/zmi * (chioni-chiono)*vf*vfi*vfi/chionm;
  }
  return gfrac;
}

/*----------------------------------------*/

void set_perms (chantype *chtyp, iontab *ions, double svrev, int pri_ion, int sec_ion)

/* Set permeabilities of secondary ion to match vrev */
/*  Does not set primary or tertiary ions */
/*  Can set chantype perms or chan perms */

{
  int j, temp;
  double chratio, chratio2, perm, perm2;
  double sign, sign2;
  double dvrev, vrev, ovrev;
#define MAXPERMITER 50

  chratio  = ions->ioni[pri_ion] / ions->iono[pri_ion];  			
  chratio2 = ions->ioni[sec_ion] / ions->iono[sec_ion];  			
  sign  = (chratio >1 ? -1: 1) * ionz[pri_ion];
  sign2 = (chratio2>1 ? -1: 1) * ionz[sec_ion];
  perm  = ions->ionp[pri_ion];
  perm2 = ions->ionp[sec_ion];
  for (ovrev=-1,dvrev=0.01,j=0; j<MAXPERMITER && abs(dvrev) > 1e-8; j++) {
      vrev = ghkv(chtyp,ions,ions->iono[PCA],ions->ioni[PCA]);   /* vrev calculated from ions */
      dvrev = svrev - vrev;
      dvrev = limit (dvrev,0.03,-0.03);
      perm2 *= 1+sign2*dvrev*40;
      perm2 = limit (perm2,0,1);
      ions->ionp[sec_ion] = perm2;
      if (ovrev==vrev) { 		// If no convergence, swap secondary ion 
	 j = 0;				//  This may indicate wrong starting perms
         temp = pri_ion;
         pri_ion = sec_ion;
         sec_ion = temp;
	 chratio  = ions->ioni[pri_ion] / ions->iono[pri_ion];  			
	 chratio2 = ions->ioni[sec_ion] / ions->iono[sec_ion];  			
	 sign  = (chratio >1 ? -1: 1) * ionz[pri_ion];
	 sign2 = (chratio2>1 ? -1: 1) * ionz[sec_ion];
	 perm  = ions->ionp[pri_ion];
	 perm2 = ions->ionp[sec_ion];
         ovrev = -1;
      }
      else ovrev = vrev;
      // fprintf (stderr,"%s %d chratio %g chratio2 %g perm %g, perm2 %g vrev %g dvrev %g\n",
      //	findsym(chtyp->ctype),j,chratio,chratio2,perm,perm2,vrev,dvrev);
  }
}

/*----------------------------------------*/

int set_perms (chantype *chtyp, iontab *ions, double svrev)

{
    int pri_ion, sec_ion;

  pri_ion = find_pri_ion(ions);
  sec_ion = find_sec_ion(ions,pri_ion);
  set_perms (chtyp, ions, svrev, pri_ion, sec_ion);
}

/*----------------------------------------*/

void makghk_tab (chan *ch)

/* Make interpolation tables for ghk current equation for "slope conductance" and "slope vrev" */
/*  See de Schutter & Smolen (1998) in "Methods in Neuronal Modeling." p 218 */
/*  Also see docacomp() in ncomp.cc, and cadfxx.n series */

{
   int i, vn, pri_ion, sec_ion;
   double **gvrevtab, **gfractab, *gp, *vp;
   double v; 
   double cvfi, gslope;

#define PERMTHR 0.1

  if (ch==NULL) {
      ncfprintf(stderr,"## calc_cvfi: missing channel\n");
      return;
  }
  if (ch->setvrev) {			// if user has set vrev for this chan
     ch->ions = init_ions(ch->ions);
     for (i=0; i<NIONS; i++) {		// copy ion concs from chan type to indiv chan
        ch->ions->ioni[i] = ch->chtyp->ions->ioni[i];
        ch->ions->iono[i] = ch->chtyp->ions->iono[i];
        ch->ions->ionp[i] = ch->chtyp->ions->ionp[i];
     }

     // if 2 major ions, set permeabilities from user-set vrev
     // otherise, set internal concentration from user-set vrev

     pri_ion = find_pri_ion(ch->ions);
     sec_ion = find_sec_ion(ch->ions,pri_ion);
     if (ch->ions->ionp[sec_ion] > PERMTHR) { // If major secondary ion
        set_perms (ch->chtyp, ch->ions, ch->vrev, pri_ion, sec_ion);
     }
     else ch->ions->ioni[pri_ion] = nernstc(ch->vrev, ch->ions->iono[pri_ion],ionz[pri_ion]); // only pri ion
  }

  // fprintf (stderr,"chan type %s %d\n", findsym(ch->ctype),ch->stype);

  if (ch->ions!=NULL) { 		// if indiv channel has ions set
     gfractab = &ch->gfractab;
     gvrevtab = &ch->gvrevtab;
  } else {
     gfractab = &ch->chtyp->gfractab;
     gvrevtab = &ch->chtyp->gvrevtab;
  }
  if (*gfractab==NULL) {		
        if (!(*gfractab=(double *)emalloc((RATESIZ+1) * sizeof(double)))) {
        ncfprintf (stderr, "Makghk_tab: can't allocate gfrac table\n");
        return;
     }
  }
  if (*gvrevtab==NULL) {		
        if (!(*gvrevtab=(double *)emalloc((RATESIZ+1) * sizeof(double)))) {
        ncfprintf (stderr, "Makghk_tab: can't allocate gvrev table\n");
        return;
     }
  }
  gp = *gfractab;
  vp = *gvrevtab;
  for (vn=0; vn<=RATESIZ; vn++) {   /* make tables from -250 to 250 mV, */
                                    /*  i.e. vn==0 -> Vm= -250 mV */
                                    /* v = Vm, calib in V */
    v = (vn - (RATESIZ/2)) / MVOLT; 
    if (v==0) v=1e-6;
    cvfi     = calc_cvfi (ch,v);		// calculate most of GHK current eqn
    gslope = calc_gslope (ch,v);		// calculate rest of g'
    *(gp + vn) = cvfi - gslope;			// calc gfrac
    *(vp + vn) = v - v*cvfi/(cvfi-gslope); 	// calc gvrev

    //  ncfprintf (stderr,"%g %g %g\n", v, v*cvfi, v-v*cvfi/(cvfi-gslope)); /* compare to cadf_7.n */
    //  ncfprintf (stderr,"%g %g %g\n", v, (double)*(gp+vn), (double)*(vp+vn)); /* */

  }  /* for (vn;;) */

}

/*----------------------------------------*/

double get_ioni (chan *ch, int ion) 
{
     double ion_p;

  if (ch->ions!=NULL) 
     ion_p = ch->ions->ioni[ion];
  else if (ch->chtyp->ions!=NULL) 
     ion_p = ch->chtyp->ions->ioni[ion];
  else 
     ion_p = ions->ioni[ion];
  return ion_p;
}

/*----------------------------------------*/

double get_iono (chan *ch, int ion) 
{
     double ion_p;

  if (ch->ions!=NULL) 
     ion_p = ch->ions->iono[ion];
  else if (ch->chtyp->ions!=NULL) 
     ion_p = ch->chtyp->ions->iono[ion];
  else 
     ion_p = ions->iono[ion];
  return ion_p;
}

/*----------------------------------------*/

double get_ionp (chan *ch, int ion) 
{
     double ion_p;

  if (ch->ions!=NULL) 
     ion_p = ch->ions->ionp[ion];
  else if (ch->chtyp->ions!=NULL) 
     ion_p = ch->chtyp->ions->ionp[ion];
  else 
     ion_p = ions->ionp[ion];
  return ion_p;
}

/*----------------------------------------*/

double get_gfrac(chan *ch, double v)

{
    int vn, vnn;
    double indx, r, rr, gfrac;
    double *vp,**gfractab;

  if (use_ghki) {
    if (ch->ions!=NULL) 			/* if ions set for indiv chan */
       gfractab = &ch->gfractab;		/* use table in chan */
    else
       gfractab = &ch->chtyp->gfractab;		/* use table in chantyp */

    if (abs(v) < 1e-20) v = 0.0;
    indx = (v * MVOLT) + (RATESIZ/2.0);		/* convert to mv + 250 */
    if (indx < 0) indx = 0;
    if (indx > RATESIZ) indx = RATESIZ;
    vn = (int)indx;				/* table starts at -250 mv */
    vnn = vn + 1;

    r = indx - vn;				/* remainder for interpolation */
    rr = 1.0 - r;
    if (*gfractab==NULL) makghk_tab(ch);
    if (ch->ions!=NULL) 			/* if ions set for indiv chan */
       gfractab = &ch->gfractab;		/* use table in chan */
    vp = *gfractab;				/* use whichever table was initialized */
    gfrac  =  vp[vn]*rr + vp[vnn]*r;		/* interpolate gvrev table */ 
  }
  else gfrac = 1.0;
  return gfrac;
}

/*----------------------------------------*/

double get_gvrev(chan *ch, double v)

{
    int vn, vnn;
    double indx, r, rr, gvrev;
    double *vp,**gvrevtab;

  if (use_ghki) {
    if (ch->ions!=NULL) 				/* if ions set for indiv chan */
       gvrevtab = &ch->gvrevtab;
    else
       gvrevtab = &ch->chtyp->gvrevtab;

    if (abs(v) < 1e-20) v = 0.0;
    indx = (v * MVOLT) + (RATESIZ/2.0);		/* convert to mv + 250 */
    if (indx < 0) indx = 0;
    if (indx > RATESIZ) indx = RATESIZ;
    vn = (int)indx;				/* table starts at -250 mv */
    vnn = vn + 1;

    r = indx - vn;				/* remainder for interpolation */
    rr = 1.0 - r;
    if (*gvrevtab==NULL) makghk_tab(ch);
    if (ch->ions!=NULL) 				/* if ions set for indiv chan */
       gvrevtab = &ch->gvrevtab;
    vp = *gvrevtab;
    gvrev  =  vp[vn]*rr + vp[vnn]*r;		/* interpolate gvrev table */ 
  }
  else gvrev = ch->vrev;
  return gvrev;
}

/*----------------------------------------*/

double callchancalc (double v, int func)

/* Call a user-defined function for calculating chan rate functions. */

{
   return (callfunc (fpnt,2,v,(double)func,0,0));
}

/*----------------------------------------*/

void set_chancalc (int ctype, int stype, int np, double (*func)(double, int))

/* set a user-defined C/C++ rate function. */

{
   chantype *ch;

  if ((ch=findchantype(ctype,stype))==NULL) initchan();
  if ((ch=findchantype(ctype,stype))==NULL) {
       ncfprintf (stderr,"set_chantype: can't set rate function: chan %s type %d not found.\n",
		findsym(ctype),stype);
  }
  else {
    if (np > ch->numparm) {
       ncfprintf (stderr,"set_chantype: param too high %d.\n", np);
    }
    else ch->parm[np].chancalc = func;
  }
}


/*----------------------------------------*/

void makchanratetab(chanparm *chp)

/* Make tables for rate functions (e.g. alpha and beta) for any 
   standard ("HH") or Markov sequential-state channel type.
   Makes "ntab" tables (e.g. am, bm, etc.).

   Hodgkin-Huxley equations for rate constants are calculated
   from the negative of the difference between the voltage and the
   resting potential.  The tables are calculated here from -250 to
   250 mv (0 means -250 mv, 500 means 250 mv), and are calibrated in
   units of per second.  A final multiplication by "timinc"
   calibrates the tables in terms of the time step.  Therefore, the
   tables need to be recalculated every time the time step changes.

   Note that the rate function called here starts with the first 
   rate func = 1, but that the rate table starts at 0.

   Note that am and bm rate functions are calculated with a
   temperature dependence factor of q.  Normally this represents a
   Q10 of 3, but am and bm are calculated with a Q10 of 2
   (Frankenhaeuser and Moore, 1963, J Physiol 169:431-467).  This
   temperature dependence is implemented in the function qrate();

*/

{
  double v,tval;
  double q;
  double *npt,*chantab;
  double *fval, voff;
  int vn,t,nfval,ratesiz;
  Symbol *f;
  double (*chanf)(double,int);
  extern double calck5n(double, int);

  static double *am,*bm;

  static double *m1tab;
  static double *m2tab;

/* ncfprintf (stderr,"maktabk begin\n"); /* */

  nfval = (int)chp->nfval;
  if (nfval < 2) nfval = 2;	/* always make space for at least 2 values */

  fval = chp->fval;
  if (!fval) {		
        if (!(fval=(double *) emalloc(nfval * sizeof(double)))) {
          ncfprintf (stderr, "Makchanstatetab: can't allocate table for %s\n",
	 		chp->funcname);
        return;
     }
    chp->fval = fval;
    if (chp->nfval < 1) {	/* make dummy values in case no rate function */
      chp->fval[0] = 1.0;
      chp->fval[1] = 1.0;
    }
  }
			/* Don't make voltage table if none specified */

  if ((nfval=(int) chp->nfval) < 1) return; 

  if (!(f=lookup(chp->funcname))) chanf = chp->chancalc;
  else {
     fpnt = f;
     chanf = callchancalc;
     if (!vidmode) ncfprintf (stdout,"## substituting the %s function\n",chp->funcname);
  }

  if (!chanf) return;	/* if no rate function defined, return */
 
  chantab = chp->chanftab;
  if (!chantab) {		
        if (!(chantab=(double *)		
                emalloc(nfval * (RATESIZ+1) * sizeof(double)))) {
        ncfprintf (stderr, "Makchanstatetab: can't allocate table for %s\n",
	 		chp->funcname);
        return;
     }
    chp->chanftab = chantab;
  }
 ratesiz = RATESIZ+1;
 npt = chantab;
 for (t=1; t<=nfval; t++) {  /* calc tables for alpha, beta rate funcs */ 
				 /*  and their m1tab, m2tab lookups */

  if (chp->nq==1) 
     q = qrate(chp);	/* find Q10 for rate func, func of V and temperature */
  else
     q = qraten(chp,t-1);	

  for (vn=0; vn<=RATESIZ; vn++) {   /* make tables from -250 to 250 mV, */
                                    /*  i.e. vn==0 -> Vm= -250 mV */
                                    /* v = Vm, calib in mV */
    v = vn - (RATESIZ/2) - (chp->voff) * MVOLT; 
    *(npt + vn) = q * (*chanf)(v,t);

      // ncfprintf (stderr,"func %d V %g rate %g %g\n", t, v, *(npt+vn), calcna5m(v,t)); /* */
      // ncfprintf (stderr,"func %d V %g rate %g %g\n", t, v, (*chanf)(v,t), calckna5h(v,t)); /* */

  }  /* for (vn;;) */
  npt += ratesiz;

 }  /* for (t;;) */

 chp->oldvf = -1e10;
 chp->oldtimincf = timinc;
}

/*----------------------------------------*/

void makchanimpltab(chanparm *chp)

/* Make implicit tables for rate functions alpha and beta for any 
   standard ("HH") channel type.  Makes 2 tables (nival=2), m1tab, m2tab.

   Offset table voltage by default offset value.
*/

{
  double *npt,*chanftab,*chanitab;
  double *am,*bm,*m1tab,*m2tab;
  double *ival,tval;
  int nival,ratesiz,vn;

  nival = (int)chp->nival;
  if (nival < 2) nival = 2;	/* always make space for at least 2 values */


  ival = chp->ival;
  if (!ival) {		
        if (!(ival=(double *) emalloc(nival * sizeof(double)))) {
          ncfprintf (stderr, "Makchanstatetab: can't allocate table for %s\n",
	 		chp->funcname);
        return;
     }
    chp->ival = ival;
  }
			/* Don't make voltage table if none specified */
  if ((nival=(int) chp->nival) < 2) return; 

  chanitab = chp->chanitab;
  if (!chanitab) {		
        if (!(chanitab=(double *)		
                emalloc(nival * (RATESIZ+1) * sizeof(double)))) {
        ncfprintf (stderr, "Makchanimpltab: can't allocate table for %s\n",
	 		chp->funcname);
        return;
     }
     chp->chanitab = chanitab;
  }

  ratesiz = RATESIZ+1;
  chanftab = chp->chanftab;
  am    = chanftab;			/* alpha, beta tables, already done */
  bm    = chanftab + ratesiz;

  m1tab = chanitab;			/* implicit tables */
  m2tab = chanitab + ratesiz;

  if (am && m1tab) { 

		/* next, calculate implicit tables */

   for (vn=0; vn<=RATESIZ; vn++) {   /* make tables from -250 to 250 mV, */

    tval  =	  (1.0 - (am[vn] + bm[vn])*0.5) /
    	    	  (1.0 + (am[vn] + bm[vn])*0.5); 

    if (tval>0) {
      m1tab[vn] = tval;				/* C-N implicit */
      m2tab[vn] = am[vn] / (1.0+(am[vn]+bm[vn])*0.5); 
    }
    else {					/* purely implicit */
      m1tab[vn] = 1.0    / (1.0+am[vn]+bm[vn]);
      m2tab[vn] = am[vn] / (1.0+am[vn]+bm[vn]); 
    }

/*  ncfprintf (stderr,"t %d %4d %10.5g %10.5g\n",
		t,vn-RATESIZ/2,am[vn],bm[vn]);  /* */

/* ncfprintf (stderr,"t %d %4d %10.5g %10.5g\n",
		t,vn-RATESIZ/2,m1tab[vn],m2tab[vn]);  /* */

   }  /* for (vn;;) */
 }
 chp->oldtiminci = timinc;
 chp->oldvi = -1e10;
 chp->oldar = -1e10;
 chp->oldbr = -1e10;
}


/*----------------------------------------*/

void makchanimpl (chanparm *chp)

/* Make alpha, beta rate tables and implicit tables for 
   one state variable (either HH type or sequential-state).
   If chp->nitab is nonzero, implicit tables will be made
   for an HH channel.
*/

{
  makchanratetab (chp);		/* make rate function tables */
  makchanimpltab (chp);		/* make implicit tables */
}

/*----------------------------------------*/

void calcchanimpl(double v, double arate, double brate, chanparm *chp)
                       
/* Used at runtime, for HH chans only. Given input voltage in
   volts, lookup HH rate constants in tables that were constructed
   as a function of millivolts above -250 mv. Do interpolation
   between adjacent table entries and calculate m1, m2 values
   from implicit lookup tables.  Tables give second-order equation
   for HH rate constants.  Create tables if they have not been
   created yet.

*/

#define MAXVOLTS 0.25
{
  int vn,vnn,ratesiz;
  static double x,r,rr,tval,ar,br;
  double a,b,ab,ab2, val;
  double *am,*bm;
  double *m1tab,*m2tab;
  double m1, m2;
  char cbuf[100];
  extern char *progname;

  if (v < -MAXVOLTS || v > MAXVOLTS) {
	  //execerror ("calcchanimpl: t %9.6g membrane voltage %5.3g out of range","stopping...");
	  // execerror (cbuf," limiting...");
	  // return;
          // sprintf (cbuf,"calcchanimpl: t %9.6g membrane voltage %5.3g out of range",simtime,v);
          // execerror (cbuf,"limiting...");

     // ncfprintf (stderr,"%s: calcchanimpl: t %9.6g membrane voltage %5.3g out of range, limiting...\n",
     // 		 			progname,simtime,v);
 
     if (v < -MAXVOLTS) v = -MAXVOLTS;
     if (v >  MAXVOLTS) v =  MAXVOLTS;
  }
  if ((!chp->chanitab) || (timinc!=chp->oldtiminci)) makchanimpl(chp);
  if (v != chp->oldvi || arate != chp->oldar || brate != chp->oldbr) {

    ratesiz = RATESIZ+1;
    m1tab = chp->chanitab;
    m2tab = chp->chanitab + ratesiz;
 
    if (abs(v) < 1e-20) v = 0.0;
    x = (v * 1000.0) + (RATESIZ/2.0);		/* convert to mv + 250 */
    if (x < 0) x = 0;
    if (x > RATESIZ) x = RATESIZ;
    vn = (int)x;				/* table starts at -250 mv */
    vnn = vn + 1;

    r = x - vn;				/* remainder for interpolation */
    rr = 1.0 - r;

    if (arate==1.0 && brate==1.0) {		/* find activation rates */
      m1  =  m1tab[vn]*rr + m1tab[vnn]*r;	/* interpolate m1 table */ 
      m2  =  m2tab[vn]*rr + m2tab[vnn]*r;	/* interpolate m2 table */ 
    }
    else {		/* non-standard tau */
      am = chp->chanftab;
      bm = chp->chanftab + ratesiz;

      if (arate==0) arate = 1e-4;
      if (brate==0) brate = 1e-4;
      ar = arate;
      br = brate;

      a = (am[vn]*rr + am[vnn]*r) * ar;
      b = (bm[vn]*rr + bm[vnn]*r) * br;
      ab  = a+b;
      ab2 = ab*0.5;

      tval  =	  (1.0 - ab2) /
    	    	  (1.0 + ab2); 

      if (tval>0) {
        m1 = tval;				/* C-N implicit */
        m2 = a / (1.0+ab2); 
      }
      else {					/* purely implicit */
        m1 = 1.0 / (1.0+ab);
        m2 = a   / (1.0+ab); 
      }
    }
    chp->oldvi = v;
    chp->oldar = arate;
    chp->oldbr = brate;
    chp->oldtiminci = timinc;
    chp->ival[0] = m1;				/* save values for impl calc */
    chp->ival[1] = m2;
  }
}
/*--------------------------------------------*/

void calcchanrate(double v, chanparm *chp)

/* Calculate rates from tables. Note that rate tables start at 0. */

{
  int i,vn,vnn,ratesiz;
  double x,r;
  double *am,*bm;
  double *chanftab;
  double *fval;


  if (v != chp->oldvf) {

    if (!chp->chanftab || (timinc!=chp->oldtimincf)) makchanratetab (chp);

   if (abs(v) < 1e-10) v = 0.0;
   x = (v * 1000.0) + (RATESIZ/2.0);		/* convert to mv + 250 */
   if (x < 0) x = 0;
   if (x > RATESIZ) x = RATESIZ;
   vn = (int)x;			     /* table starts at -250 mv */
   r = x - vn;			     /* remainder for interpolation */
 
   fval = chp->fval;
   ratesiz = RATESIZ+1;
   chanftab = chp->chanftab;
   for (i=0; i<chp->nfval; i++, chanftab+=ratesiz) {
        register double chftab;
     chftab = chanftab[vn];		     /* alpha, beta, etc. */
     fval[i] = chftab + (chanftab[vn+1] - chftab) * r; /* interpolate table */
   } 

   chp->oldvf = v;

/* ncfprintf (stderr,"v %g vn %d am %g bm %g\n",
		v,vn,fval[0],fval[1]); /* */
  }
}

/*--------------------------------------------*/

void chanrate(chan *cpnt, chanparm *chp)
             
/* Given chan, find its voltage in volts, lookup alpha and 
   beta rate functions in tables that were constructed as a 
   function of millivolts above -250 mv. Do interpolation 
   between adjacent table entries. Used by markov sequential-state 
   integration method. Note that "vext" is the voltage on
   the external side of the membrane from ephaptic effects,
   and "vph" is the shift in gating for Ca chans from external pH.
*/

{
#ifdef DEBUG
  if (!cpnt || !cpnt->comp1) {
    ncfprintf (stderr, "nc: error in chanrate\n");
    return;
  }
#endif
  calcchanrate(cpnt->comp1->v - cpnt->comp1->vph - vext(cpnt) - ccavoff(cpnt) -
		  (chp->pn ? cpnt->voffsh: cpnt->voffsm), chp);
}

/*--------------------------------------------*/

void chanimpl(chanparm *chp, chan *cpnt)
             
/* Given chan, find its voltage in volts, lookup alpha and 
   beta rate functions in tables that were constructed as a 
   function of millivolts above -250 mv. Do interpolation 
   between adjacent table entries. Used by HH  
   integration method.
*/

{
  comp *comp1;

//  if (!cpnt || !cpnt->comp1) return;
  comp1 = cpnt->comp1;
  if (chp->pn == 0) 
       calcchanimpl(comp1->v - comp1->vph - vext(cpnt) - ccavoff(cpnt) - cpnt->voffsm, 
			cpnt->arate, cpnt->brate, chp);
  else calcchanimpl(comp1->v - comp1->vph - vext(cpnt) - ccavoff(cpnt) - cpnt->voffsh, 
			cpnt->crate, cpnt->drate, chp);
}

/*----------------------------------------*/

void makchan_dftab(chantype *chtyp, int ndftab)

/* Make tables for driving force functions (e.g. GHK current equation)
   for any channel type.
   Makes "ntab" tables (e.g. am, bm, etc.).
*/

{

}

/*--------------------------------------------*/

void calcchan_dftab(double v, chantype *chtyp)

{

}

/*--------------------------------------------*/

double calcchaninf(double v, chanparm *chp, double arate, double brate)

{
  double alpha, beta;

  if (chp->nfval < 2) {		/* must calc rate functions here */

  }
  calcchanrate(v,chp);
  alpha = chp->fval[0] * arate;
  beta  = chp->fval[1] * brate;
  return (alpha / (alpha + beta));
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

double calcchantau(double v, chanparm *chp, double arate, double brate)

{
  double alpha, beta;

 calcchanrate(v,chp);
 alpha = chp->fval[0] * arate;
 beta = chp->fval[1] * brate;
 return (timinc / (alpha + beta));
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

double calcchantaua(double v, chanparm *chp)

{
  double alpham;

 calcchanrate(v,chp);
 alpham = chp->fval[0];
 return (timinc / alpham);
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

double calcchantaub(double v, chanparm *chp)

{
  double betam;

 calcchanrate(v,chp);
 betam = chp->fval[1];
 return (timinc / betam);
}

/*--------------------------------------------*/

double al0(chan *cpnt)

/* Return the rate constant calculated as an implicit value */
/*  by "calcchanrate()" and placed in fval[0] (normally alpha). */

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

double bet0(chan *cpnt)

/* Return the rate constant calculated as an implicit value */
/*  by "calcchanrate()" and placed in fval[0] (normally beta). */

{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[0];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

double al1(chan *cpnt)
{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[0];
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

double bet1(chan *cpnt)
{
    chanparm *chp;

  chp = &cpnt->chtyp->parm[1];
  chanrate(cpnt, chp);
  return chp->fval[1];
}

/*--------------------------------------------*/

double pinf (double v, int elnum, int p)

/* Return the equilibrium activation for a state variable 
   param (e.g m, h) at a given voltage v. */

{
  attrib *apnt;
  mattrib *mpnt;
  chantype *cht;
  chanparm *chp;
  double voff, arate, brate, val;

  if (interp) varcopyu();
  if (!(apnt=getattrib(elnum))) {
    ncfprintf (stderr,"pinf: can't find element %d\n",elnum);
    return 0;
  }
  if (!(cht=getchantype(apnt->ctype, apnt->stype))) {
     ncfprintf (stderr,"pinf: incorrect channel type for element %d\n",elnum);
     return 0;
  }
  if (p < cht->numparm) {
    chp = &cht->parm[p];
    switch (apnt->ctype) {
     case CA:
     case K:
     case ClCa:
     case KCa:
     case NA:
	mpnt = (mattrib *)apnt; break;
     default: return 0; break;
    }
    if (chp->pn) { 		/* h */
      voff =  mpnt->voffsh;
      arate = mpnt->tauc;
      brate = mpnt->taud;
    }
    else {			/* m */
      voff =  mpnt->voffsm;
      arate = mpnt->taua;
      brate = mpnt->taub;
    } 
    if (voff==NULLVAL)  voff = 0;
    if (arate==NULLVAL) arate = 1;
    if (brate==NULLVAL) brate = 1;
    arate = 1/arate;		/* rate is inverse tau */
    brate =  1/brate;
    if (apnt->ctype==KCa) {	/* special case for Kca chans */
      val = bkcainf(v-voff, dcai, mpnt, arate, brate);
    }
    else val = calcchaninf(v-voff,chp,arate,brate);
  }
  else val = 0;
  return val;
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

double ptau (double v, int elnum, int p)

/* Return the absolute tau (in sec) for activation at a given voltage.
*/

{
  elem *epnt;
  attrib *apnt;
  mattrib *mpnt;
  chantype *cht;
  chanparm *chp;
  chan *ch;
  double voff, arate, brate, val;

  if (interp) varcopyu();
  if (!(apnt=getattrib(elnum))) {
    ncfprintf (stderr,"pinf: can't find element %d\n",elnum);
    return 0;
  }
  if (!(cht=getchantype(apnt->ctype, apnt->stype))) {
     ncfprintf (stderr,"pinf: incorrect channel type for element %d\n",elnum);
     return 0;
  }
  if (p < cht->numparm) {
    chp = &cht->parm[p];
    switch (apnt->ctype) {
     case CA:
     case K:
     case ClCa:
     case KCa:
     case NA:
	  mpnt = (mattrib *)apnt; break;
     default: return 0; break;
    }
    if (chp->pn) { 		/* h */
      voff =  mpnt->voffsh;
      arate = mpnt->tauc;
      brate =  mpnt->taud;
    }
    else {			/* m */
      voff =  mpnt->voffsm;
      arate = mpnt->taua;
      brate =  mpnt->taub;
    } 
    if (voff==NULLVAL)  voff = 0;
    if (arate==NULLVAL) arate = 1;
    if (brate==NULLVAL) brate = 1;
    arate = 1/arate;		/* rate is inverse tau */
    brate = 1/brate;
    if (apnt->ctype==KCa) {	/* special case for Kca chans */
      val = bkcatau(v-voff, dcai, mpnt, arate, brate, chp);
    }
    else val = calcchantau(v-voff,chp,arate,brate);
  }
  else val = 0;
  return val;
}

/*--------------------------------------------*/

double ptaua (double v, int elnum, int p)

/* taua with user-defined tau for alpha applied */

{
  attrib *apnt;
  mattrib *mpnt;
  chantype *cht;
  chanparm *chp;
  chan *ch;
  double voff, alpha, val;

  if (interp) varcopyu();
  if (!(apnt=getattrib(elnum))) {
    ncfprintf (stderr,"pinf: can't find element %d\n",elnum);
    return 0;
  }
  if (!(cht=getchantype(apnt->ctype, apnt->stype))) {
     ncfprintf (stderr,"pinf: incorrect channel type for element %d\n",elnum);
     return 0;
  }
  if (p < cht->numparm) {
    chp = &cht->parm[p];
    switch (apnt->ctype) {
     case CA:
     case K:
     case NA:
	  mpnt = (mattrib *)apnt; break;
     default: return 0; break;
    }
    if (chp->pn) { 		/* h */
      voff =  mpnt->voffsh;
      alpha = mpnt->tauc;
    }
    else {			/* m */
      voff =  mpnt->voffsm;
      alpha = mpnt->taua;
    } 
    if (voff==NULLVAL)  voff = 0;
    if (alpha==NULLVAL) alpha = 1;
     val = (calcchantaua(v-voff,chp)*alpha);
  }
  else val = 0;
  return val;
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

double ptaub (double v, int elnum, int p)

/* taub with user-defined tau for beta applied */

{
  attrib *apnt;
  mattrib *mpnt;
  chantype *cht;
  chanparm *chp;
  chan *ch;
  double voff, beta, val;

  if (interp) varcopyu();
  if (!(apnt=getattrib(elnum))) {
    ncfprintf (stderr,"pinf: can't find element %d\n",elnum);
    return 0;
  }
  if (!(cht=getchantype(apnt->ctype, apnt->stype))) {
     ncfprintf (stderr,"pinf: incorrect channel type for element %d\n",elnum);
     return 0;
  }
  if (p < cht->numparm) {
    chp = &cht->parm[p];
    switch (apnt->ctype) {
     case CA:
     case K:
     case NA:
	  mpnt = (mattrib *)apnt; break;
     default: return 0; break;
    }
    if (chp->pn) { 		/* h */
      voff =  mpnt->voffsh;
      beta = mpnt->taud;
    }
    else {			/* m */
      voff =  mpnt->voffsm;
      beta = mpnt->taub;
    } 
    if (voff==NULLVAL)  voff = 0;
    if (beta==NULLVAL)  beta = 1;
    val = (calcchantaub(v-voff,chp)*beta);
  }
  else val = 0;
  return val;
}

/*--------------------------------------------*/

datum ctaua (datum &v, datum &elname)

{
 v.val = ptaua (v.val, (int)elname.val, 0);
 return v;
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

datum ctaub (datum &v, datum &elname)

{
 v.val = ptaub (v.val, (int)elname.val, 0);
 return v;
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

datum ctauc (datum &v, datum &elname)

{
 v.val = ptaua (v.val, (int)elname.val, 1);
 return v;
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

datum ctaud (datum &v, datum &elname)

{
 v.val = ptaub (v.val, (int)elname.val, 1);
 return v;
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

datum ctaue (datum &v, datum &elname)

{
 v.val = ptaua (v.val, (int)elname.val, 2);
 return v;
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

datum ctauf (datum &v, datum &elname)

{
 v.val = ptaub (v.val, (int)elname.val, 2);
 return v;
}

/*--------------------------------------------*/

datum minf (datum &v, datum &elname)

{
 v.val = pinf (v.val, (int)elname.val, 0);
 return v;
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

datum mtau (datum &v, datum &elname)
{
 v.val = ptau (v.val, (int)elname.val, 0);
 return v;
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

datum hinf (datum &v, datum &elname)
{
 v.val = pinf (v.val, (int)elname.val, 1);
 return v;
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

datum htau (datum &v, datum &elname)
{
 v.val = ptau (v.val, (int)elname.val, 1);
 return v;
}

/* - - - - - - - - - - - - - - - - - - - - - -*/

datum (*ninf) (datum &v, datum &elname) = minf;
datum (*ntau) (datum &v, datum &elname) = mtau;
datum (*dinf) (datum &v, datum &elname) = hinf;
datum (*dtau) (datum &v, datum &elname) = htau;
datum (*cinf) (datum &v, datum &elname) = minf;
datum (*ctau) (datum &v, datum &elname) = mtau;
datum (*chinf) (datum &v, datum &elname) = hinf;
datum (*chtau) (datum &v, datum &elname) = htau;

/*--------------------------------------------*/

double getntresp (chan *ch, comp *cpnt)

/* Calculate a channel's response to ligands present. 
   Agonist ligands have positive weight, and
   antagonist ligands have negative weight.

   Not totally correct, but simple; otherwise,
   need to have a different binding func for 
   each ligand (see "rantag()" and "ragonist()") below.
*/
   
{
   int type;
   ntcomp *npnt;
   double a,ampl;

  ampl = 0.0;
  for (npnt=cpnt->ntpnt; npnt; npnt=npnt->ntpnt) { /* search list */
	type = npnt->ctype - GLU;
        if (type <= 0 || type >= NNTR) continue;
        a = ch->chtyp->respamp[type];
        ampl += a * npnt->val / npnt->n; 
  }
  if (ampl < 0) ampl = 0;
  return ampl; 
}

/*--------------------------------------------*/

double getagonresp (chan *ch, comp *cpnt)

/* Calculate a channel's response to agonist ligands present */ 
/* Agonist ligands have positive weight */
   
{
   int type;
   ntcomp *npnt;
   double a,ampl;

  ampl = 0.0;
  for (npnt=cpnt->ntpnt; npnt; npnt=npnt->ntpnt) { /* search list */
	type = npnt->ctype - GLU;
        if (type <= 0 || type >= NNTR) continue;
        a = ch->chtyp->respamp[type];
	if (a > 0) ampl += a * npnt->val / npnt->n; 
  }
  if (ampl < 0) ampl = 0;
  return ampl; 
}

/*--------------------------------------------*/

double getantresp (chan *ch, comp *cpnt)

/* Calculate a channel's response to competitive antagonists present.  
   Antagonists have negative weight in the respamp[] table.
    Sum all the negative weights and invert.
*/
   
{
   int type;
   ntcomp *npnt;
   double a,ampl;

  ampl = 0.0;
  for (npnt=cpnt->ntpnt; npnt; npnt=npnt->ntpnt) { /* search list */
	type = npnt->ctype - GLU;
        if (type <= 0 || type >= NNTR) continue;
        a = ch->chtyp->respamp[type];
	if (a < 0) ampl -= a * npnt->val / npnt->n; 
  }
  if (ampl < 0) ampl = 0;
  return ampl; 
}

/*--------------------------------------------*/

int isntresp (chan *ch, comp *cpnt)

/* Find if a channel can respond to ligands present */ 
   
{
   int type;
   ntcomp *npnt;
   double a,ampl;

  ampl = 0.0;
  for (npnt=cpnt->ntpnt; npnt; npnt=npnt->ntpnt) { /* search list */
	type = npnt->ctype - GLU;
        if (type <= 0 || type >= NNTR) continue;
        a = ch->chtyp->respamp[type];
	if (a > 0) return 1;
  }
  return 0;
}

/*--------------------------------------------*/

int isantag (chan *ch, comp *cpnt)

/* Find if a channel can respond to antagonists present */ 
   
{
   int type;
   ntcomp *npnt;
   double a,ampl;

  ampl = 0.0;
  for (npnt=cpnt->ntpnt; npnt; npnt=npnt->ntpnt) { /* search list */
	type = npnt->ctype - GLU;
        if (type <= 0 || type >= NNTR) continue;
	a = ch->chtyp->respamp[type];
	if (a < 0) return (1);
  }
  return 0;
}

/*--------------------------------------------*/

double rchanf(chan *cpnt)

/* Forward binding rate constant for channel opening.
   Calculated from channel conductance.
   For use in dochan2().
 */

{
   double conduct;

 if (!cpnt) return (0.0);
 if (cpnt->comp2)
   conduct=((synap *)cpnt->comp2)->conduct;	/* find synap conductance */
 else
   conduct=cpnt->conduct;			/* find chan conductance */
 return (conduct*qrate(cpnt->chtyp->parm));	/* rate per second */
}

/*--------------------------------------------*/

double rchanr(chan *cpnt)

/* Reverse binding rate constant for channel opening.
   Calculated from channel conductance.
   For use in dochan2().
 */

{
   double conduct;

 if (!cpnt) return (0.0);
 if (cpnt->comp2)
   conduct=((synap *)cpnt->comp2)->conduct;	/* find synap conductance */
 else
   conduct=cpnt->conduct;			/* find chan conductance */
 return ((1.0-conduct)*qrate(cpnt->chtyp->parm));	/* rate per second */
}

/*--------------------------------------------*/

double rt(chan *cpnt)

/* Return a constant rate scaled by time increment. */
/* Uses first q10 value. */

{
  return (qrate(cpnt->chtyp->parm));
}

/*--------------------------------------------*/

double rt2(chan *cpnt)

/* Return a constant rate scaled by time increment. */
/* Uses first q10 value. */

{
  return (qraten(cpnt->chtyp->parm,1));
}

/*--------------------------------------------*/

double rtq(chan *cpnt, int p, int n)

/* Return a constant rate scaled by time increment.
   Q10 value is from table of q values for chanparm p, indexed by
   n, i.e. it can be different from the Q10 value for the "first
   param", for example, a different rate function such as one
   sensitive to voltage.  First n value is n=0. */

{
   chanparm *parm;

  if (!(parm= &cpnt->chtyp->parm[p])) {
        return (0);
  } 
  return (qraten(parm,n));
}

/*--------------------------------------------*/

double rnt(chan *cpnt)

/* Forward binding rate constant for neurotransmitter.
   For use in Markov schemes where agonist and antagonist
   concentrations are weighted positive and negative,
   respectively and are summed to get one effect.
   (Not strictly very correct, but simple.)

   Multiply neurotrans concentration by timinc, 
   Q10 value, and scaling factor "trconc", then
   return this for use in dochan2().
 */

{
   double nt;
   synap *spnt;

  if (!cpnt) return (0.0);
  if (isntresp(cpnt,cpnt->comp1)) {      /* find postsyn ligand conc */
     nt = getntresp (cpnt,cpnt->comp1);  /* if it exists, use it */
     					 /* (already mult by trconc/mesgconc). */
  }
  else {                                 /* otherwise use direct signal */
    if (spnt=(synap*)cpnt->comp2){ 	 /* set in "rchan()" */
      nt = spnt->conduct;   		 /* synaptic nt conc */
    }
    else nt = 0;
  }
  return (nt*qrate(cpnt->chtyp->parm));    /* rate prop. to nt conc.*/
}                                       /*  per second */

/*--------------------------------------------*/

double rca(chan *cpnt)

/* Forward binding rate constant for Ca in KCa and ClCa chans.
   For use in Markov schemes where Ca gates channel.

   Return Ca concentration, multiplied by timinc and Q10 value.
 */

{
   double ca;
   cacomp *capnt;

  if (!cpnt) return (0.0);
  if (!(capnt=cpnt->comp1->capnt)) {
        ca = dcai;
  }
  else ca = capnt->cais[0];
  return (ca*qrate(cpnt->chtyp->parm));    /* rate prop. to ca conc.*/
}                                       /*  per second */

/*--------------------------------------------*/

double rca2(chan *cpnt)

/* Forward binding rate constant for Ca in KCa and ClCa chans.
   For use in Markov schemes where Ca gates channel.
   Q10 value is from "second param", i.e. it can be different
   from the Q10 value for the "first param", e.g. a different
   rate function such as one sensitive to voltage.

   Return Ca concentration, multiplied by timinc and Q10 value.
 */

{
   double ca;
   cacomp *capnt;
   chanparm *parm;

  if (!cpnt) return (0.0);
  if (!(capnt=cpnt->comp1->capnt)) {
        ca = dcai;
  }
  else ca = capnt->cais[0];
  return (ca*qraten(cpnt->chtyp->parm,1));    /* rate prop. to ca conc.*/
}                                 /*  per second */

/*--------------------------------------------*/

double rclca(chan *cpnt)

/* Forward binding rate constant for Ca in KCa and ClCa chans.
   For use in Markov schemes where Ca and voltage gate channel.
   Q10 is 8 for a 120 mV change in voltage.
   Yang et al. (2008) TMEM16A confers ...  Nature 455:1210-1216.

   Return Ca concentration, multiplied by timinc, and voltage Q10 value.
 */

{
   double ca,vr;
   cacomp *capnt;
   chanparm *parm;

  if (!cpnt) return (0.0);
  if (!(capnt=cpnt->comp1->capnt)) {
        ca = dcai;
  }
  else ca = capnt->cais[0];
  vr = exp(log(dclcavs)*(cpnt->comp1->v - vext(cpnt) - 
			      compcavoff(cpnt->comp1) - -0.04)/0.120);
  return (ca*vr*qrate(cpnt->chtyp->parm));    /* rate prop. to ca conc.*/
}                                 /*  per second */

/*--------------------------------------------*/

double rclcac(chan *cpnt)

/* Forward binding rate constant for Ca in KCa and ClCa chans.
   For use in Markov schemes where Ca and voltage gate channel.
   Q10 is 8 for a 120 mV change in voltage.
   Yang et al. (2008) TMEM16A confers ...  Nature 455:1210-1216.

   Return Ca concentration, multiplied by timinc, and voltage Q10 value.
 */

{
   double ca,vr;
   cacomp *capnt;
   chanparm *parm;

  if (!cpnt) return (0.0);
  if (!(capnt=cpnt->comp1->capnt)) {
        ca = dcai;
  }
  else ca = capnt->cais[int(dsclcac>capnt->cashell? capnt->cashell:dsclcac)];	// calcium from core
  vr = exp(log(dclcavsc)*(cpnt->comp1->v - vext(cpnt) - 
			      compcavoff(cpnt->comp1) - -0.04)/0.120);
  return (ca*vr*qrate(cpnt->chtyp->parm));    /* rate prop. to ca conc.*/
}                                 /*  per second */

/*--------------------------------------------*/

double ragonist(chan *cpnt)

/* Forward binding rate constant for neurotransmitter.
   For use in Markov schemes where separate functions
   are used for agonist and antagonist concentrations.
 
   Return neurotrans concentration, multiplied by
   timinc and scaling factor "trconc" for use in dochan2().
 */

{
   double nt;
   synap *spnt;

  if (!cpnt) return (0.0);
  if (isntresp(cpnt,cpnt->comp1)) {      /* find postsyn ligand conc */
     nt = getagonresp (cpnt,cpnt->comp1);  /* if it exists, use it */
  }
  else {                                 /* otherwise use direct signal */
    if (spnt=(synap*)cpnt->comp2) {
      nt = spnt->conduct * spnt->trconc;   /* synaptic nt conc */
    }
    else nt = 0;
  }
  return (nt*qrate(cpnt->chtyp->parm));    /* rate prop. to nt conc.*/
}                                       /*  per second */

/*--------------------------------------------*/

double rantag(chan *cpnt)

/* Forward binding rate constant for competitive antagonist 
   ligands.  Return neurotrans concentration, multiplied by
   timinc and scaling factor "trconc" for use in dochan2().
 */

{
   double nt;
   synap *spnt;

  if (!cpnt) return (0.0);
  if (isantag(cpnt,cpnt->comp1)) {      /* find postsyn ligand conc */
     nt = getantresp (cpnt,cpnt->comp1);  /* if it exists, use it */
     return (nt*qrate(cpnt->chtyp->parm));    /* rate prop. to nt conc.*/
  }
  else return 0;
}                                       /*  per second */

/*--------------------------------------------*/

void mak2state(chanstate *spnt, double (*frate)(chan *cpnt),
				double (*rrate)(chan *cpnt), double rate)

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

  {
     double a,b;

   if (!spnt) {
	ncfprintf (stderr,"mak2state: no chanstate\n");
	return;
   }

   a = rate;
   b = rate;
					/*      states     */

   					/*    0   <->  1   */ 
   spnt[0].numtrans   = 1;
   spnt[0].cond       = 0;
   spnt[0].trans  [0] = 1;
   spnt[0].trate  [0] = frate;
   spnt[0].ratemul[0] = a;
   spnt[0].rateo  [0] = 5;		/* use frate = for noise */

   spnt[1].numtrans   = 1;
   spnt[1].cond       = 1.0;		/* 1 = the open state */
   spnt[1].trans  [0] = 0;
   spnt[1].trate  [0] = rrate;
   spnt[1].ratemul[0] = b;
   spnt[1].rateo  [0] = 5;
}

