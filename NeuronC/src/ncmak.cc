/* Module ncmak in program nc */

/* Makes cable branches, nodes, compartments,
   connections, and synapses. 
*/


#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "ncsub.h"
#include "ncomp.h"
#include "ncelem.h"
#include "nconst.h"
#include "control.h"
#include "ncio.h"

/* #define MALLOPT		/* */

#ifdef MALLOPT
#include <malloc.h>
#endif

#include "nclist.h"

#define DEBUG 		/* */

#ifdef DEBUG
#include "ncdebug.h"
#endif

node *nodepnt=NULL;
node *nodend=NULL;
elem *elempnt=NULL;
elem *elemend=NULL;
elem *oelemend=NULL;

extern double iono[];
extern double ioni[];
extern int    ionz[];
extern int numplots;
void initpl(int c);

extern double pigmlen[];	/* path length through o.s. (in ncstim) */

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

double ncabs(double x);
char *emalloc(unsigned int n);
double calctau(double tau, double timestep);
double akcacalc(double v, double ca, double tau, double d1, double k1);
double bkcacalc(double v, double ca, double tau, double d2, double k2);
void execerror (const char *s, const char *t);
void freelst (conlst *clst);
void maklst(conlst **head, conn *cpnt);
void addlst(conlst **head, conlst *h2);
void dellst(conlst **head, conn *cpnt);
void makdyad(dyadlst **head, synap *cpnt);
chan *makchan(attrib *apnt, comp *comp1, comp *comp2, int saved);
chattrib *make_chan (elem *epnt, int ctype, int stype);
void coneset(int n);
void rodiset(int n);
void chr2set(int n);
void savephotrec(photrec *rpnt);
void restorphotrec(photrec *rpnt);
void runrec(photrec *rpnt, int dtimestep);
double rsens(photrec *rpnt, double wavel, int filt);
double calcchaninf(double v, chanparm *chp, double arate, double brate);
double dochani(chan *chpnt, double crit);
void drand_setstate (char *state);
void initstate (unsigned long int seed, char *state, int n);
void delcomp  (comp *cpnt);
void delchan  (chan *chpnt);
void delsynap (synap *spnt);
void delphotrec (photrec *rpnt);
void delelem  (elem *epnt);
void einstall (elem *nelem);
elem *findelem (int num);
elem *findelemne (int num);
void unsetnodlst (node *nodp, elem *epnt);
node *findnode(nodeint node1a, nodeint node1b, 
		nodeint node1c, nodeint node1d, const char *s);
void setphotrectim(photrec *rpnt, double timec, double loopgain);
chantype *getchantype(int ctype, int cnum);
double setbind(synap *spnt, double chligand);
double setcond(synap *spnt, double bound);
double setbind2 (synap *spnt, double active);

double gjalpha (double v, gj *gjpnt);
double gjbeta  (double v, gj *gjpnt);
double gamdev(double a);
double qcond (chantype *ch);
double nernstv(double cratio);
void docacompi (cacomp *capnt);
double *makfiltarr (int newsiz, int oldsiz, double *oarr, double val);
lpfilt *maksynfilt (int stype, int nfilt, double offset, double tfall, double val,
                double *timec, double timestep);
char *findsym(int num);
char *prnode (int n1, int n2, int n3, int n4);
double qcavoff (chanparm *chp);
void restorstate(void);
void efree (void *ptr);
cacomp *makcacomp(elem *epnt, cattrib *apnt, comp *comp1, double area);
double fzero (double (*func)(double), double near_val);
double ccavoff(chan *ch);
double vext(conn *ch);
iontab *init_ions(iontab *itab);

void ninithash(void);
void einithash(void);

/*------------------------------------------------*/

void initmk(void)

/* reset all lists */

{
#ifdef MALLOPT
  static int called=0;
if (!called) {
   called = 1;
   mallopt(M_MXFAST,sizeof(comp)+1);	/* allocate comps, nodes, etc. faster */
   mallopt(M_NLBLKS,100);		/* set size of holding blocks  */
   mallopt(M_GRAIN,8);			/* set grain size  */
}
#endif

cumelem=0;				/* cumulative lengths of lists */
cumattr=0;
cumnattr=0;
cumcattr=0;
cumnode=0;
cumcomp=0;
cumconn=0;
cumchan=0;
cumsynap=0;
cumload=0;
cumphotrec=0;
cumrecst=0;
cumcacomp=0;
cumcycacomp=0;
cumcycgcomp=0;
cumclst=0;
cumdlst=0;
cumrand=0;
cumlpfilt=0;
cumntcomp=0;
cumgj=0;
cumdbuf=0;
cumnbuf=0;
delelmt=0;


elempnt=NULL;				/* pointers to lists */
elemend=NULL;
oelemend=NULL;
nodepnt=NULL;
nodend=NULL;
compnt=NULL;
compend=NULL;
connpnt=NULL;
connend=NULL;
synpnt=NULL;
synend=NULL;
loadpnt=NULL;
loadend=NULL;
chanpnt=NULL;
chanend=NULL;
recpnt=NULL;
recend=NULL;
recspnt=NULL;
recsend=NULL;

  if (makestim) {
   reccum = 0;                     /* reset number of receptors */
   reclist = NULL;
   reclend = NULL;
  }
}

/*------------------------------------------------*/

void ncleanup(void)

/* erase all lists */

{
   elem *epnt,*tepnt;
   comp *cpnt,*tcpnt;
   conn *chpnt,*tchpnt;
   conn *cnpnt,*tcnpnt;
   node *npnt, *tnpnt;
   synap *spnt,*tspnt;
   load  *lpnt,*tlpnt;
   photrec *rpnt,*trpnt;
   recstim *rspnt, *trspnt;
   recnod *npt, *tnpt;

 for (epnt=elempnt; epnt; ) {		/* elements */
    tepnt = epnt;
    epnt = epnt->next;
    delelem(tepnt);			/* delete the element */
 }
 einithash();				/* erase the element hash table */

 for (npnt=nodepnt; npnt; ) {		/* nodes */
    freelst (npnt->elemlst);		/* free the node's element list */
    tnpnt = npnt;
    npnt = npnt->next;
    efree (tnpnt);
 }
 ninithash();				/* erase the node hash table */

 for (cpnt=compnt; cpnt; ) {		/* compartments */
    tcpnt = cpnt;
    cpnt = cpnt->next;
    delcomp (tcpnt);
 }

 for (cnpnt=connpnt; cnpnt; ) {		/* connections */
    tcnpnt = cnpnt;
    cnpnt = cnpnt->next;
    efree (tcnpnt);
 }

 for (lpnt=loadpnt; lpnt; ) {		/* loads */
    tlpnt = lpnt;
    lpnt = lpnt->next;
    efree (tlpnt);
 }

 for (spnt=synpnt; spnt; ) {		/* synapses */
    tspnt=spnt;
    spnt= (synap*)spnt->next;
    delsynap (tspnt);
 }

 for (chpnt=chanpnt; chpnt; ) {		/* channels */
    tchpnt=(chan*)chpnt;
    chpnt= chpnt->next;
    delchan((chan*)tchpnt);
 }

 for (rpnt=recpnt; rpnt; ) {		/* photoreceptors */
    trpnt=rpnt;
    rpnt=(photrec*)rpnt->next;
    delphotrec(trpnt);
 }

 for (rspnt=recspnt; rspnt; ) {		/* photoreceptor stimuli */
    trspnt = rspnt;
    rspnt = (recstim *)rspnt->next;
    efree (trspnt);
 }

 initpl(0);				/* delete all plots */

 if (makestim) {
   for (npt=reclist; npt; ) {
    tnpt = npt;
    npt = npt->next;
    efree (tnpt);
   }
 }

initmk();				/* reset all lists and counts */
}

/*------------------------------------------------*/

void freelst (conlst *pnt)

/* Free a conlst (list of connections), part of either
   a node, compartment, or cacomp, but not its connection (conpnt).
*/

{
   conlst *lpnt,*tlpnt;

   for (lpnt=pnt; lpnt; ) {
        tlpnt = lpnt;
        lpnt = lpnt->next;
	efree (tlpnt);
    }
}

/*------------------------------------------------*/

void freendp (node *nodp, elem *epnt)
             
/* Reverse of "setnodlst()" in modcode.c.
   Subtracts a "conlst" from a node's element list,
    but maintains the list by patching pointers.
   Frees the list pointer from a node to an element.
   Used when the element is to be freed. 
*/

{
  unsetnodlst (nodp,epnt);
}

/*------------------------------------------------*/

void copyattr(attrib *src, attrib *dest)

/* Copy params but not pointers from one attrib to another. */

{
   chattrib *chsrc,*chdest;
   cattrib *casrc,*cadest;
   nattrib *nsrc,*ndest;

   if (!src || !dest) {
      ncfprintf (stderr,"copyattr: Missing pointer.\n");
      return;
   }

   dest->ctype   = src->ctype;
   dest->stype   = src->stype;
   dest->lptr    = src->lptr;

   switch (src->ctype) {

      case GABA:
      case GLY:
      case KAINATE:
      case NA:
      case K:
      case ClCa:
      case KCa:
	  chsrc = (chattrib *)src;
	  chdest = (chattrib *)dest;
	  chdest->vrev    = chsrc->vrev;
	  chdest->voffsm  = chsrc->voffsm;
	  chdest->voffsh  = chsrc->voffsh;
	  chdest->maxcond = chsrc->maxcond;
	  chdest->density = chsrc->density;
	  chdest->ndensity = chsrc->ndensity;
	  chdest->unitary = chsrc->unitary;
	  chdest->caperm  = chsrc->caperm;
	  chdest->cakd    = chsrc->cakd;
	  chdest->cahc    = chsrc->cahc;
	  chdest->taua    = chsrc->taua;
	  chdest->taub    = chsrc->taub;
	  chdest->tauc    = chsrc->tauc;
	  chdest->taud    = chsrc->taud;
	  chdest->taue    = chsrc->taue;
	  chdest->tauf    = chsrc->tauf;
	  chdest->d1      = chsrc->d1;
	  chdest->d2      = chsrc->d2;
	  chdest->k1      = chsrc->k1;
	  chdest->k2      = chsrc->k2;
	  chdest->n       = chsrc->n;
          break;

      case NMDA:
      case CGMP:
      case SYN2:
      case AMPA:
      case CA:
      case CACOMP:
	  casrc = (cattrib *)src;
	  cadest = (cattrib *)dest;

	  cadest->vrev    = casrc->vrev;
	  cadest->voffsm  = casrc->voffsm;
	  cadest->voffsh  = casrc->voffsh;
	  cadest->maxcond = casrc->maxcond;
	  cadest->density = casrc->density;
	  cadest->ndensity = casrc->ndensity;
	  cadest->unitary = casrc->unitary;
	  cadest->caperm  = casrc->caperm;
	  cadest->cakd    = casrc->cakd;
	  cadest->cahc    = casrc->cahc;
	  cadest->taua    = casrc->taua;
	  cadest->taub    = casrc->taub;
	  cadest->tauc    = casrc->tauc;
	  cadest->taud    = casrc->taud;
	  cadest->taue    = casrc->taue;
	  cadest->tauf    = casrc->tauf;
	  cadest->n      = casrc->n;

	  cadest->exch   = casrc->exch;
	  cadest->pump   = casrc->pump;
	  cadest->cabuf  = casrc->cabuf;
	  cadest->cao    = casrc->cao;
	  cadest->cai    = casrc->cai;
	  cadest->tcai   = casrc->tcai;
	  cadest->kex    = casrc->kex;
	  cadest->ekm    = casrc->ekm;
	  cadest->vmax   = casrc->vmax;
	  cadest->pkm    = casrc->pkm;


          cadest->cas	 = casrc->cas;                     /* CICR parameters */
	  cadest->cas2	 = casrc->cas2;
	  cadest->vm2	 = casrc->vm2;
	  cadest->vm3	 = casrc->vm3;
	  cadest->ncicr	 = casrc->ncicr;
	  cadest->mcicr	 = casrc->mcicr;
	  cadest->pcicr	 = casrc->pcicr;
	  cadest->kacicr = casrc->kacicr;
	  cadest->kfcicr = casrc->kfcicr;
	  cadest->krcicr = casrc->krcicr;
	  cadest->k1cicr = casrc->k1cicr;
	  cadest->k2cicr = casrc->k2cicr;
	  cadest->c1cicr = casrc->c1cicr;

	  cadest->ip3i	 = casrc->ip3i;                    /* IP3 parameters */;
	  cadest->bip3	 = casrc->bip3;
	  cadest->b2ip3	 = casrc->b2ip3;
	  cadest->vip3	 = casrc->vip3;
	  cadest->v2ip3	 = casrc->v2ip3;
	  cadest->v3ip3	 = casrc->v3ip3;
	  cadest->v4ip3	 = casrc->v4ip3;
	  cadest->mip3	 = casrc->mip3;
	  cadest->hip3	 = casrc->hip3;
	  cadest->oip3	 = casrc->oip3;
	  cadest->a2ip3	 = casrc->a2ip3;
	  cadest->a3ip3	 = casrc->a3ip3;
	  cadest->mtypeip3 = casrc->mtypeip3;
	  cadest->d1ip3	 = casrc->d1ip3;
	  cadest->d2ip3	 = casrc->d2ip3;
	  cadest->d3ip3	 = casrc->d3ip3;
	  cadest->d4ip3	 = casrc->d4ip3;
	  cadest->k3ip3	 = casrc->k3ip3;

	  cadest->bkd    = casrc->bkd;
	  cadest->bmax   = casrc->bmax;
	  cadest->btot   = casrc->btot;
	  cadest->btoti  = casrc->btoti;
	  cadest->cabnd  = casrc->cabnd;
	  cadest->sarea  = casrc->sarea;
	  cadest->mg     = casrc->mg;
	  cadest->cashell= casrc->cashell;
	  cadest->caoshell= casrc->caoshell;
	  cadest->caflg  = casrc->caflg;

          break;

   case VESNOISE:
   case CCHNOISE:
	  nsrc  = (nattrib *)src;
	  ndest = (nattrib *)dest;
	  ndest->n     = nsrc->n;   
	  ndest->vsize = nsrc->vsize;   
	  ndest->vcov  = nsrc->vcov;   
	  ndest->cov   = nsrc->cov;   
	  ndest->refr  = nsrc->refr;   
	  ndest->unitary = nsrc->unitary;   
	  ndest->tauf  = nsrc->tauf;   
	  ndest->rseed = nsrc->rseed;   
          break;
   }
}

/*------------------------------------------------*/

cattrib *makcattr(void)

/* make a new neural element to contain calcium channel attributes.
   Return a pointer to the new element. */

{
    cattrib *apnt;

#ifdef DEBUG 
  if (debug & NCMAK && debugz & 8)  ncfprintf (stderr,"makcattrib %d\n",cumelem);
#endif

  if ((apnt=(cattrib *)emalloc(sizeof(cattrib))) == (cattrib *)NULL) {
     ncfprintf (stderr,"no space left for cattrib %d\n",cumelem);
     return ((cattrib *)NULL);  
  }
  apnt->stype   = NULLVAL;
  apnt->vrev    = NULLVAL;
  apnt->voffsm  = NULLVAL;
  apnt->voffsh  = NULLVAL;
  apnt->taua    = NULLVAL;
  apnt->taub    = NULLVAL;
  apnt->tauc    = NULLVAL;
  apnt->taud    = NULLVAL;
  apnt->taue    = NULLVAL;
  apnt->tauf    = NULLVAL;
  apnt->maxcond = NULLVAL;
  apnt->density = NULLVAL;
  apnt->ndensity= NULLVAL;
  apnt->caperm  = NULLVAL;
  apnt->n       = NULLVAL;
  apnt->unitary = NULLVAL;
  apnt->cakd    = NULLVAL;
  apnt->cahc    = NULLVAL;
  apnt->trconc  = NULLVAL;

  apnt->cao     = NULLVAL;
  apnt->cai     = NULLVAL;
  apnt->tcai    = NULLVAL;
  apnt->kex     = NULLVAL;
  apnt->vmax    = NULLVAL;

  apnt->cas     = NULLVAL;
  apnt->cas2    = NULLVAL;
  apnt->ip3i    = NULLVAL;
  apnt->vm2     = NULLVAL;
  apnt->vm3     = NULLVAL;
  apnt->ncicr   = NULLVAL;
  apnt->mcicr   = NULLVAL;
  apnt->pcicr   = NULLVAL;
  apnt->kacicr  = NULLVAL;
  apnt->kfcicr  = NULLVAL;
  apnt->krcicr  = NULLVAL;
  apnt->k1cicr  = NULLVAL;
  apnt->k2cicr  = NULLVAL;
  apnt->bip3    = NULLVAL;
  apnt->vip3    = NULLVAL;
  apnt->c1cicr  = NULLVAL;
  apnt->b2ip3   = NULLVAL;
  apnt->a2ip3   = NULLVAL;
  apnt->a3ip3   = NULLVAL;
  apnt->v2ip3   = NULLVAL;
  apnt->v3ip3   = NULLVAL;
  apnt->v4ip3   = NULLVAL;
  apnt->k3ip3   = NULLVAL;
  apnt->d1ip3   = NULLVAL;
  apnt->d2ip3   = NULLVAL;
  apnt->d3ip3   = NULLVAL;
  apnt->d4ip3   = NULLVAL;
  apnt->hip3    = NULLVAL;
  apnt->mip3    = NULLVAL;
  apnt->oip3    = NULLVAL;
  apnt->mtypeip3 = NULLVAL;

  apnt->pkm     = NULLVAL;
  apnt->ekm     = NULLVAL;
  apnt->btot    = NULLVAL;
  apnt->btoti   = NULLVAL;
  apnt->bkd     = NULLVAL;
  apnt->bmax    = NULLVAL;
  apnt->cabnd   = NULLVAL;
  apnt->cashell = NULLVAL;
  apnt->caoshell = NULLVAL;
  apnt->sarea   = NULLVAL;
  apnt->mg      = NULLVAL;
  apnt->pump    = 0;
  apnt->caflg   = 0;			/* set to 1 if cacomp in synapse */
  apnt->cabuf   = 0;
  apnt->exch    = 0;
  apnt->cicr    = 0;
  apnt->ip3     = 0;
  apnt->attpnt  = (attrib *)NULL;
  cumcattr++;
  return (apnt); 
}

/*------------------------------------------------*/

chattrib *makchattr(void)

/* make a new neural element to contain channel attributes.
   Return a pointer to the new element. */

{
    chattrib *apnt;

#ifdef DEBUG 
  if (debug & NCMAK && debugz & 8)  ncfprintf (stderr,"makattrib %d\n",cumelem);
#endif

  if ((apnt=(chattrib *)emalloc(sizeof(chattrib))) == (chattrib *)NULL) {
     ncfprintf (stderr,"no space left for attrib %d\n",cumelem);
     return ((chattrib *)NULL);  
  }
  apnt->ctype   = NULLVAL;
  apnt->stype   = NULLVAL;
  apnt->lptr    = (elem *)NULL;
  apnt->vrev    = NULLVAL;
  apnt->voffsm  = NULLVAL;
  apnt->voffsh  = NULLVAL;
  apnt->maxcond = NULLVAL;
  apnt->rextern = NULLVAL;
  apnt->density = NULLVAL;
  apnt->ndensity= NULLVAL;
  apnt->unitary	= NULLVAL;
  apnt->caperm	= NULLVAL;
  apnt->cakd    = NULLVAL;
  apnt->cahc    = NULLVAL;
  apnt->taua    = NULLVAL;
  apnt->taub    = NULLVAL;
  apnt->tauc    = NULLVAL;
  apnt->taud    = NULLVAL;
  apnt->taue    = NULLVAL;
  apnt->tauf    = NULLVAL;
  apnt->n	= NULLVAL;
  apnt->d1 	= NULLVAL;
  apnt->d2 	= NULLVAL;
  apnt->k1 	= NULLVAL;
  apnt->k2 	= NULLVAL;
  apnt->trconc 	= NULLVAL;
  apnt->attpnt  = (attrib *)NULL;
  cumattr++;
  return (apnt); 
}

/*------------------------------------------------*/

nattrib *maknattr(void)

/* make a new neural element to contain noise attributes.
   Return a pointer to the new element. */

{
    nattrib *apnt;

#ifdef DEBUG 
  if (debug & NCMAK && debugz & 8)  ncfprintf (stderr,"maknattrib %d\n",cumelem);
#endif

  if ((apnt=(nattrib *)emalloc(sizeof(nattrib))) == (nattrib *)NULL) {
     ncfprintf (stderr,"no space left for nattrib %d\n",cumelem);
     return ((nattrib *)NULL);  
  }
  apnt->ctype   = NULLVAL;
  apnt->stype   = NULLVAL;
  apnt->lptr    = (elem *)NULL;
  apnt->n       = NULLVAL;
  apnt->vsize   = NULLVAL;
  apnt->unitary = NULLVAL;
  apnt->tauf    = NULLVAL;
  apnt->rseed   = 0;
  apnt->vcov    = NULLVAL;
  apnt->cov     = NULLVAL;
  apnt->refr    = NULLVAL;
  apnt->attpnt  = (attrib *)NULL;
  cumnattr++;
  return (apnt); 
}

/*------------------------------------*/

load *makload(loadelem *epnt, comp *comp1)
               
/* make a new resistive load and link it to the load list. */

{
    load *lpnt;
    double rp, r;
 
    if ((lpnt=(load *)emalloc(sizeof(load))) == (load*)NULL) {
      ncfprintf (stderr,"no space left for load %d\n",cumload+1);
      return ((load *)NULL);  
    }
    lpnt->next = (load *)NULL;
    if (!loadpnt) loadpnt = lpnt;  	/* save head if first load */
    if (loadend)
      loadend->next = lpnt;
    loadend = lpnt;

    loadend->comp1  = comp1; 		 /* compartment we're connected to */
    loadend->comp2  = NULL; 		 /* */
    loadend->compe  = NULL; 		 /* ext compartment */
    loadend->compp  = NULL; 		 /* ext compartment */
    maklst (&comp1->clst,(conn*)loadend); /* link presynaptic cell conn */
    loadend->ctype  = epnt->ctype;	 /* set load type */
    if (loadend->ctype == ELEMENT)
		 loadend->ctype = LOAD;  /* default load type */
 					 /* use default values if necess */
   
    if ((r=epnt->r)==NULLVAL)   r = 1e6; /* resistance of load */
    if ((rp=epnt->vrev)==NULLVAL) rp = 0; /* reversal potential */
    if (r==0) r = 1.0;
    loadend->conduct = 1.0 / r;
    loadend->vrev = rp;                 /* reversal potential */
    cumload++; 				/* increment total */
    return (loadend); 
}

/*------------------------------------*/

void modload(loadelem *ldpnt)
               
/* modify an existing load */

{
    loadelem *epnt;
    load *lpnt;
    double r;


    if (!(epnt=(loadelem *)findelem(ldpnt->modif))) {
      ncfprintf (stderr,"modload: can't find element %d\n",ldpnt->modif);
      return;
    }

    lpnt = (load *)epnt->lptr;	   /* get pointer to synapse from elem */
    if (lpnt==(load *)NULL) {
       ncfprintf (stderr,"modload: can't find load for elem %d\n",ldpnt->elnum);
       return;  
    }
    if (lpnt->ctype != LOAD) {
       ncfprintf (stderr,"modload: element %d is not a load\n",ldpnt->modif);
       return;  
    }

    if (ldpnt->r != NULLVAL) {
       epnt->r = ldpnt->r; 
       r = ldpnt->r;
       if (r==0.0) r = 1e8;
       lpnt->conduct = 1.0/r;  				/* conductance */
    }
    if (ldpnt->vrev!=NULLVAL)  lpnt->vrev = ldpnt->vrev; /* reversal pot */
}

/*------------------------------------*/

char *makrand (int rndsiz, const char *mesg, int num)

{
      char *stpnt;

 if ((stpnt=(char *)emalloc(rndsiz)) == (char *)NULL) {
   ncfprintf (stderr,"upk_chan: No space left for %s rand state %d.\n", mesg,num);
 }
 cumrand++;
 return stpnt;
}

/*------------------------------------*/

void linksyn (conn *cpnt)

/* Add a connection to the list of synapses. */

{
 if (!synpnt) synpnt = (synap *)cpnt;  	/* save head if first synap */
 if (synend)
   synend->next = (synap *)cpnt;
 synend = (synap *)cpnt;
}

/*------------------------------------*/

ntcomp *findntcomp (comp *cpnt, int type)

/* Find a nt concentration given its type. */

{
     int found;
     ntcomp *npnt;

   if (!cpnt) {
	// ncfprintf (stderr,"findntcomp: can't find compartment for %d\n",type);
	return (ntcomp*)NULL;
   }
   for (found=0,npnt=cpnt->ntpnt; npnt; npnt=npnt->ntpnt) { /* search list */
      if (npnt->ctype==type) {
	 found = 1;
	 break;
      }
   }
   if (!found) return (ntcomp*)NULL;
   else        return npnt;
}

/*------------------------------------*/

ntcomp *makntcomp (comp *cpnt, int type)

/* Make nt concentration for a compartment. 
   Each time, increment count.  This means that
   there are several sources that should be averaged.
*/

{
     ntcomp *npnt;

   if (!cpnt) return (ntcomp*)NULL;
   if (!(npnt=findntcomp(cpnt,type))) {
      if ((npnt=(ntcomp *)emalloc(sizeof(ntcomp))) == (ntcomp*)NULL) {
        ncfprintf (stderr,"makntcomp: no space left for ntcomp\n");
        return ((ntcomp *)NULL);  
      }
      npnt->ntpnt = cpnt->ntpnt;
      cpnt->ntpnt = npnt;
      npnt->ctype = type;
      npnt->n = 0;
      cumntcomp++;
   }
   npnt->n++;
   return (npnt);
}

/*------------------------------------*/

ntcomp *makaddntcomp (comp *cpnt, int type, double val, int n)

/* Add a nt concentration to a compartment. 
   Each time, add count "n".  This means that
   there are several sources that should be averaged.
*/

{
     ntcomp *npnt;

   npnt = makntcomp(cpnt,type);
   npnt->val += val;
   npnt->n += n-1;
   return (npnt);
}

/*------------------------------------*/

void delntcomp (ntcomp *npnt)

/* Delete a list of ntcomps. */

{
   ntcomp *npt,*tpt;

  if (!npnt) return;

  for (npt=npnt; npt; npt=tpt) {	/* erase list */
    tpt = npt->ntpnt;
    efree (npt);
  }
}

/*------------------------------------*/

int setnt (comp *cpnt, int type, double val)

/* Set concentration of nt in a compartment */

{
     ntcomp *npnt;

  if (!(npnt=findntcomp(cpnt,type))) {
	npnt = makntcomp (cpnt, type);
  }
  if (npnt) {
    npnt->val = val;
    return 1;
  }
  else return 0;
}

/*------------------------------------*/

double addnt (comp *cpnt, int type, double val)

/* Add a concentration of nt to a compartment. */
/* Allows several release sites to sum */

{
     ntcomp *npnt;

  if (!(npnt=findntcomp(cpnt,type))) {
	npnt = makntcomp (cpnt, type);
  }
  if (npnt) {
    npnt->val += val;
    return npnt->val;
  }
  else return 0;
}

/*------------------------------------*/

int setntfromcomp (comp *cpnt1, comp *cpnt2, double offset, double gain, int type)

/* Set concentration of nt of type in a compartment cpnt2 */
/*     from voltage in another compartment cpnt1*/

{
     ntcomp *npnt;

  if (cpnt1==NULL) return 0;

  if (!(npnt=findntcomp(cpnt2,type))) {
	npnt = makntcomp (cpnt2, type);
  }
  if (npnt) {
    npnt->val = (cpnt1->v - offset) * gain;
    return 1;
  }
  else return 0;
}

/*------------------------------------*/

double getnt (comp *cpnt, int type)

/* Return neurotransmitter concentration for a compartment. */
/*  Idea is that nt can be puffed onto a node and the    */
/*  compartment keeps track of each such substance. */

/* If a nt comp has not already been made, make a new one. */

{
     ntcomp *npnt;

   if (!(npnt=findntcomp(cpnt,type))) return 0;
   else return (npnt->val / npnt->n);
}

/*------------------------------------*/

int isnt (comp *cpnt, int type)

/* Return 1 if nt comp has been found. */

/* If a nt comp has not already been made, make a new one. */

{
     ntcomp *npnt;

   if (!(npnt=findntcomp(cpnt,type))) return 0;
   else return 1;
}

/*------------------------------------*/

double ghkvi (chantype *chtyp, iontab *ions, double vrev)

/* calculate denominator of GHK potential for a channel type */
/*  but leave out Ca terms as these will be added later */

/*  Use const field eqn from Lewis 1979, J.Physiol 286:417.*/

{
    int i, ionv, ionv2;
    double perm, ci, df, nf, cavoff;

   cavoff = qcavoff(chtyp->parm) * dcaspvrev;
   df = exp ( FR * (vrev-cavoff) / ktemp); /* multiplier for Ca perm */
   nf = 1.0 / ( 1.0 + df);
   for (ci=i=0; i<NIONS; i++) {
       if (i==PCA) continue;	      /* don't add Ca terms here */
       perm = ions->ionp[i];
       ionv = ionz[i];
       ionv2 = ionv*ionv;             /* perm multiplier of 4 for Ca chans */
       switch(ionv) {
	 case -1:					/* Cl, etc */
          	ci += perm * ions->iono[i]; break;
	 case 1:					/* Na, K, etc */
          	ci += perm * ions->ioni[i] * ionv; break;
	 case 2:					/* Mg++ etc. */
          	ci += perm * ions->ioni[i] * ionv2 * nf * df; break;
	}
     }
   return (ci);
}
/*------------------------------------*/

double ghkvo (chantype *chtyp, iontab *ions, double vrev)

/* calculate numerator of GHK potential for a channel type */
/*  but leave out Ca terms as these will be added later */

/*  Use const field eqn from Lewis 1979, J.Physiol 286:417.*/

{
    int i, ionv, ionv2;
    double perm, co, df, nf, cavoff;

   cavoff = qcavoff(chtyp->parm) * dcaspvrev;
   df = exp ( FR * (vrev-cavoff) / ktemp); /* multiplier for Ca perm */
   nf = 1.0 / ( 1.0 + df);
   for (co=i=0; i<NIONS; i++) {
       if (i==PCA) continue;	      /* don't add Ca terms here */
       perm = ions->ionp[i];
       ionv = ionz[i];
       ionv2 = ionv*ionv;             /* perm multiplier of 4 for Ca chans */
       switch(ionv) {
	 case -1:					/* Cl, etc */
          	co += perm * ions->ioni[i]; break;
	 case 1:					/* Na, K, etc */
          	co += perm * ions->iono[i]; break;
	 case 2:
          	co += perm * ions->iono[i] * ionv2 * nf; break; /* allow Mg++, etc */
	}
     }
   return (co);
}

/*------------------------------------*/

double ghkv (chantype *chtyp, iontab *ions, double cao, double cai)

/* Calculate GHK potential for a channel type. */

/* Na, K, and Cl concentrations are in ioni[] and iono[]. */
/* Ca concentration is given in cao, cai and is assumed to vary. */
/* The terms "df" and "nf" allow ions with valence of 2 to be */
/*  combined with Na, K, and Cl ions. */

/*  See derivation of const field eqn in Lewis 1979, J.Physiol 286:417. */

/* "cavoff" is surface potential, empirically fit to data in Lewis (1979) */ 
/*  It is derived from a gating shift of 18 mV/factor of 10 change in cao */
/*  Then it is multiplied by "dcasurf" (.18) to give surface potential here */

{
    int j, ctype, ion, ionv, ionv2;
    double perm, ci, co, df, nf;
    double ionox, ionix, cavoff, err, vrev;

   ctype = chtyp->ctype;
   vrev = chtyp->vrev;               /* start by assuming vrev = 0 */
   cavoff = qcavoff(chtyp->parm) * dcaspvrev;
   for (err=1.0,j=0; j<20 && ncabs(err) > 1e-4; j++) {
     df = exp ( FR * (vrev-cavoff) / ktemp); /* multiplier for Ca perm */
     nf = 1.0 / ( 1.0 + df);
     ci = ghkvi(chtyp,ions,vrev);
     co = ghkvo(chtyp,ions,vrev);
     ion = PCA;
     perm = ions->ionp[ion];
     ionv = ionz[ion];
     ionv2 = ionv*ionv;               /* perm multiplier of 4 for Ca chans */
     switch(ionv) {
	 case 2:
		if (ion==PCA) {
		  ionox = cao;
		  ionix = cai;
	  	}
	  	else {
	 	  ionox = ions->iono[ion];
		  ionix = ions->ioni[ion];
	  	}
          	co += perm * ionox * ionv2 * nf;
          	ci += perm * ionix * ionv2 * nf * df;
	  break;
      }
      err = vrev;
      vrev = nernstv (co/ci);           /* find reversal potential */
      vrev += cavoff;
      err -= vrev;
      switch (ctype) {	/* empirically derived relaxation multipliers */
        case CA:   vrev += err * .45; break;
        case NMDA:
        case AMPA: 
        case CHRC: 
        case CGMP: vrev += err * .06; break;
        case NA:   vrev += err * .0015; break;
        case GABA: 
        case GLY: 
        case KAINATE: 
        case K: 
        case ClCa:    
        case KCa:    
         default:
                  vrev -= err * .015; break;
		break;
      };
      /* ncfprintf (stderr,"%d %s %d vrev %g err * %g\n",
		j, findsym(chtyp->ctype),chtyp->stype,vrev,err); /* */
   }
   return (vrev);
}

/*------------------------------------*/

double ghkv (chantype *chtyp, double cao, double cai)

{
   return ghkv(chtyp,chtyp->ions,cao,cai); 
}

/*------------------------------------*/

double ghkv (chan *ch, double cao, double cai)

{
   return ghkv(ch->chtyp,ch->ions,cao,cai); 
}

/*------------------------------------*/

unsigned short crc (unsigned short crc, unsigned char val)

/* Given a remainder up to now, return the new CRC after one
   char is added. */

{
    int i;
    unsigned int ans;

    ans = crc ^ val << 8;
    for (i=0; i<8; i++) {
       if (ans & 0x8000)
	   ans = (ans <<= 1) ^ 4129;
       else
	   ans <<= 1; 
    }
    return ans;
}

/*------------------------------------*/

int noduniq(nodeint nodea, nodeint nodeb, 
		nodeint nodec, nodeint noded, int seed)

/* Make unique int from 4 integer (short int) values (e.g. node * number). */
/* Output val needs to be 4 byte integer for rsrandom(). */

{
    unsigned int x;
    unsigned short crcval;

   x = crcval = 0;
   crcval = crc (crcval,  noded & 0x00ff);
   crcval = crc (crcval, (noded & 0xff00) >> 8);
   crcval = crc (crcval,  nodec & 0x00ff);
   crcval = crc (crcval, (nodec & 0xff00) >> 8);
   x = crcval;
   x <<= 16;
   crcval = 0;
   crcval = crc (crcval,  nodeb & 0x00ff);
   crcval = crc (crcval, (nodeb & 0xff00) >> 8);
   crcval = crc (crcval,  nodea & 0x00ff);
   crcval = crc (crcval, (nodea & 0xff00) >> 8);
   x += crcval;

   if (seed < 0) seed = rseed;
   x ^= seed;			/* random modifier for phot noise */

   return (x & 0x7fffffff);
}

/*------------------------------------*/

void maksseed (nattrib *napnt, synap *spnt)

/* make a unique rseed for each synapse */

{
    unsigned int nrseed;

  if ((nrseed=napnt->rseed) >= 0) {
    if (nrseed==0) 
      nrseed = noduniq(27091, spnt->num, 12491, spnt->num,srseed^rseed^4323);
    spnt->vstate = makrand(RNDSIZ,"vesicle interval",cumsynap);
    if (spnt->vstate) initstate (nrseed,spnt->vstate,RNDSIZ);
    if (spnt->vstate) drand_setstate (spnt->vstate);
  }
  if (spnt->vcov > 0) {
    if ((nrseed=napnt->rseed) >= 0) {
      if (nrseed==0) 
        nrseed=noduniq(18097, spnt->num+1, 11699, spnt->num+1,srseed^rseed^6779);
      spnt->gstate = makrand(RNDSIZ,"vesicle size",cumsynap);
      if (spnt->gstate) initstate (nrseed+11,spnt->gstate,RNDSIZ);
    }
  }
  spnt->vtint = gamdev(1);
  restorstate();
}

/*------------------------------------*/

synap *maksynap(synapse *epnt, comp *comp1, comp *comp2)
               
/* Make a new synapse and link it to the synapse list. */

{
    int sf1,sf1h,sf2,sf3,nt,stop,secmsg,sens;
    synap *spnt;
    attrib *apnt;
    chattrib *ntpnt;
    nattrib *napnt;
    double st, sn, sg, so, nkd, npow, ckd, chc, sc, rp, msgc, mx, re, trc, eg, vg, rrp, mrrp, ms, rrpg;
    double sf1hg, sf1ho, tf2, tf3, vsize;
    synapse *oepnt;
    void makexp(void);
 
    if ((spnt=(synap *)emalloc(sizeof(synap))) == (synap *)NULL) {
      ncfprintf (stderr,"no space left for synapse %d\n",cumsynap+1);
      return ((synap *)NULL);  
    }
    spnt->next = (synap *)NULL;
    linksyn (spnt);			/* add to list of synapses */

    cumsynap++; 			/* increment total */
    spnt->num = cumsynap;

    spnt->comp1  = comp1; 		/* presynaptic cell */
    spnt->comp2  = comp2; 		/* postsynaptic cell */
    spnt->compe  = NULL; 		/* external compartment */
    spnt->compp  = NULL; 		/* external pH compartment */
    spnt->num1   = comp1->num;	/* comp # for presynap cell */
    spnt->num2   = comp2->num;	/* comp # for postsynap cell */
    maklst (&comp1->clst,(conn*)spnt);/* link presynaptic cell conn */
    maklst (&comp2->clst,(conn*)spnt);/* link postsynaptic cell conn */
    spnt->ctype  = epnt->ctype;	/* set synapse type */
    if (spnt->ctype == ELEMENT)
		 spnt->ctype = SYNAPSE;  /* default synapse type */
    					/* use default values if necess */
    if ((sn=epnt->ngain)==NULLVAL) sn = dsn; 	/* neurotrans gain */
    if ((vg=epnt->vgain)==NULLVAL) vg = dvg; /* neurotrans gain */
    if ((eg=epnt->caegain)==NULLVAL) eg = dscaeg; /* ca power ves release gain */
    if ((sc=epnt->curve)==NULLVAL) sc = dsc; 	/* synap curve: lin or expon */
    if ((sg=epnt->cgain)==NULLVAL)  sg = dsg; 	/* default cycG gain */
    if ((so=epnt->coff)==NULLVAL)   so = 1.0; 	/* default cycG offset */
    if ((ckd=epnt->ckd)==NULLVAL)   ckd = dckd;	/* default cycG Kd */
    if ((chc=epnt->chc)==NULLVAL)   chc = dchc;	/* default cycG Hill coeff */
    if ((sens=epnt->sens)==NULLVAL)  sens = V; 	/* presyn sensitivity */
    if ((nt=epnt->ntact)==NULLVAL)   nt = OPEN; /* neurotransmitter action */
    if ((secmsg=epnt->secmsg)==NULLVAL)secmsg = 0;  /* second messenger */
    if ((nkd=epnt->nkd)==NULLVAL)    nkd = dskd; /* half-max saturation point */
    if ((npow=epnt->npow)==NULLVAL) npow = dshc; /* hill coeff at postsyn rec */
    if ((st=epnt->thresh)==NULLVAL)  st = dst;     /* synaptic threshold */
    if ((rrp=epnt->rrpool)==NULLVAL) rrp = dsrrp;  /* readily releasible pool */
    if ((rrpg=epnt->rrpoolg)==NULLVAL) rrpg = dsrrg; /* readil releas pool gain */
    if ((mrrp=epnt->mrrpool)==NULLVAL) mrrp = dsmrrp;  /* max readily releasible pool */
    if ((ms=epnt->maxsrate)==NULLVAL) ms = dsms;   /* maximum sustained rate */

    if ((sf1=(int)(epnt->nfilt1))==NULLVAL)sf1=1;  /* number of filters */
    if (!epnt->timec1) epnt->timec1=makfiltarr(sf1, 0, (double*)NULL, 1.0);
    if (epnt->timec1)
      if (epnt->timec1[0]==0) epnt->timec1[0]= dfta; /* tau of filters (msec) */

    if ((sf1h=(int)(epnt->nfilt1h))==NULLVAL)sf1h=0; /* number of h filters */
    if ((sf1hg=(epnt->filt1hg))==NULLVAL)sf1hg=0;    /* gain of h filter */
    if ((sf1ho=(epnt->filt1ho))==NULLVAL)sf1ho=0;    /* dc offset of h filter */
    if (epnt->timec1h)
      if (epnt->timec1h[0]==0) epnt->timec1h[0]= dftah; /* high pass tau (ms) */

    if ((sf2=(int)(epnt->nfilt2))==NULLVAL)sf2=1;/* number of filters */
    if (!epnt->timec2) epnt->timec2=makfiltarr(sf2, 0, (double*)NULL, 1.0);
    if (epnt->timec2)
      if (epnt->timec2[0]==0) epnt->timec2[0]=dftb; /* tau of filters (msec) */
    if ((tf2=epnt->tfall2)==NULLVAL) tf2 = 0.0;	 /* falling tau of filt (msec)*/

    if ((sf3=(int)(epnt->nfilt3))==NULLVAL)sf3=0;/* number of filters */
    if ((tf3=epnt->tfall3)==NULLVAL) tf3 = 0.0;	 /* falling tau of filt (msec)*/
    if ((trc=epnt->trconc)==NULLVAL) {
	if (dstr==0) dstr = 1e-4;
	    trc = dstr;
	// trc = 1.0;
    }
    if ((msgc=epnt->mesgconc)==NULLVAL) msgc=dsmsgc;  /* 2nd mesg conc, def overridden by chan def below */
    if ((rp=epnt->vrev)==NULLVAL) {	/* reversal potential */
       rp =  ghkv (getchantype(AMPA,0),dcao,dcai);
       spnt->setvrev = 0;
    }
    else spnt->setvrev = 1;

    spnt->trconc = trc;
    spnt->thresh = st;
    if (mrrp < 1) mrrp = 1;
    if (rrp < 1)   rrp = 1;
    if (rrp > mrrp) rrp = mrrp;
    spnt->rrpool = rrp;		/* actual size of readily releasible pool */
    spnt->rrpoolg = rrpg;	/* readily releasible pool gain mult */
    spnt->maxrrpool = mrrp;	/* maximum size of readily releasible pool */
    spnt->maxsrate = ms;	/* maximum sustained rate */
    spnt->nkd    = nkd;
    spnt->npow   = npow;
    spnt->ckd    = ckd;
    spnt->chc    = chc;
    spnt->vrev = rp;			/* reversal potential */
    spnt->maxcond = NULLVAL;		/* unset to allow makchan to set */
    spnt->transrate = 0;
    if ((mx=epnt->maxcond)!=NULLVAL) spnt->maxcond = mx;
    if ((re=epnt->rextern)!=NULLVAL) spnt->rextern = re;
    else			     spnt->rextern = 0;
    spnt->oldc1   = 0.0;
    spnt->oldc2   = 0.0;
    spnt->ntact = nt;
    if (nt==CLOSE) secmsg = 1;	/* use second-messenger pathway */
    spnt->secmsg = secmsg;
				/* allocate filter space if filters specified */

    if (epnt->attpnt) {			/* if synapse has chan attributes */

      nattrib  *unitsyn=(nattrib*)NULL;
      cattrib *chsyn=(cattrib*)NULL;
      for (stop=0,apnt=epnt->attpnt; apnt && !stop; apnt=apnt->attpnt) {
           ntchan   *chpnt;

       switch (apnt->ctype) {

     /* The idea here is to make a synaptic channel if a
	specific postsynaptic channel type has been defined,
	or if N, tauf, or unit has been defined, or if channel 
	noise has been defined. If channel noise has been defined, 
	but a channel type (and N, tauf, and unit) have not been 
	defined, then make a new chattrib and link it before the
	chnoise nattrib, then use the normal method for defining a
	postsynaptic channel with noise (channel chattrib comes first,
	then nattrib comes after it in the attrib list). Always stop
	after the first channel definition, look at noise attribs 
	further below.
      */

         case NUNIT:  unitsyn = (nattrib *)apnt;	/* remember N, unit */
			      break;
						/* default:        */ 
         case CCHNOISE: if (((nattrib*)apnt)->n!=0){ /* stop if if nchan = 0 */
         		  chsyn = makcattr();	/* make dummy chan */
			  chsyn->ctype = SYN2; 
			  chsyn->attpnt = apnt;
			  chsyn->tauf = ((nattrib*)apnt)->tauf;
			  apnt = chsyn;
			}
			else break;	/* stop if nchan == 0 */
         case AMPA: 
         case KAINATE: 
         case CGMP: 
         case NMDA: 
         case GABA: 
         case GLY: 
         case SYN2: 
	   if (unitsyn) {	/* if N or unit was defined before resp */
		chattrib *chapnt= (chattrib *)apnt;

		if (unitsyn->n !=NULLVAL) 
		  if (chapnt->n==NULLVAL) chapnt->n = unitsyn->n;
		
		if (unitsyn->unitary !=NULLVAL) {
		  if (chapnt->unitary==NULLVAL) 
			chapnt->unitary = unitsyn->unitary;
		}
	   }
           chpnt = (ntchan*)makchan(apnt,spnt->comp2,(comp*)spnt, epnt->nocondens); 

	   if (chpnt) {
            switch (apnt->ctype) {		/* possibly make cacomp */
	         cacomp *capnt;
		 double area, lcap;

	      case SYN2:
	      case AMPA:
	      case CGMP:
	      case NMDA:
	       if (((cattrib *)apnt)->caflg)	/* user has specified cacomp */
  	        if (!(capnt=comp2->capnt)) { 

    		 lcap = dcm;
    		 if (lcap==0) lcap = 1e-6;
		 area = comp2->cap/lcap;

    		 if ((capnt=makcacomp(epnt,(cattrib*)apnt,comp2,area))==NULL) { 
      	        ncfprintf (stderr,"makca: no space for ca comp %d\n",comp2->num);
      		   return ((synap*)NULL);
 		 }
		 else comp2->capnt = capnt;
  	       }
	     break;
	     default: break;
            }
           }  /* if (chpnt) */

	   if (secmsg) spnt->resp2 = chpnt;	/* Give 2nd mesg chan props */
           else        spnt->resp1 = chpnt;	/* Give postsyn chan props */
						/*  to replace filt3. */

          ntpnt = (chattrib *)apnt;
	  if (chpnt) {
	     if (epnt->trconc==NULLVAL) spnt->trconc = chpnt->chtyp->trconc;
    	     msgc = chpnt->chtyp->trconc;
             spnt->maxcond = chpnt->maxcond;
             spnt->rextern = chpnt->rextern;
	     if ((rp=ntpnt->vrev)!=NULLVAL) {	/* channel rev potential */
	        chpnt->vrev = rp;
		if (ncabs(rp-chpnt->chtyp->vrev) >= 0.001) chpnt->setvrev = 1;
             }
             if ((rp=epnt->vrev)!=NULLVAL) {	/* synaptic rev potential */
	        chpnt->vrev = rp;		/* if explicitly set */
		if (ncabs(rp-chpnt->chtyp->vrev) >= 0.001) chpnt->setvrev = 1;
             }
	     if (info>=3 && chpnt->setvrev==1) 
		     ncfprintf (stderr,"setting vrev of %s %d to %g\n",findsym(chpnt->ctype),chpnt->stype,rp);
	  }
	 if (chsyn) { efree(chsyn); }
	 stop = 1;
	 break;

	default: break;

       }     /* switch (ctype) */

                
      }  /* for (apnt= ; ; ) */
         if (unitsyn) unitsyn = (nattrib *)NULL;
    }  /* if (epnt->attpnt) */

    if (spnt->maxcond == NULLVAL) spnt->maxcond=dmaxsyn *
					 qcond(getchantype(SYN2,0));

    spnt->vstate = (char *)NULL;
    spnt->vsites = 0;
    spnt->cov = 1.0;
    spnt->vcov = 0.0;
    spnt->vrefr = 0;
    spnt->vtint = 1.0;		/* max vesicle interval */
    spnt->vtime = 0;
    spnt->vflag = 0;
    spnt->mesg1 = 0;
    spnt->mesg2 = 0;
    if ((vsize=epnt->vsize)==NULLVAL)	/* user set ves size when no noise */ 
	 vsize = dvsz;		
    spnt->vsize = vsize;
     
    for (apnt=epnt->attpnt; apnt; apnt=apnt->attpnt) {/* set up ves noise parms*/
       napnt = (nattrib *)apnt;
       switch (napnt->ctype) {

 	  case NUNIT:
 	  case CCHNOISE:
	  break;

 	  case VESNOISE:
	  if (napnt->n==NULLVAL)     	/* use defaults if not defined */	
	      spnt->vsites = (int)(dsvn);	
	  else
	      spnt->vsites = (int)(napnt->n);	

	  if (napnt->vsize==NULLVAL)
	      spnt->vsize = vsize;	
	  else
	      spnt->vsize = napnt->vsize;	

	  if (napnt->vcov==NULLVAL)
	      spnt->vcov = 0;	
	  else
	      spnt->vcov = napnt->vcov;	
	  if (spnt->vcov  > 1) spnt->vcov = 1.0;
	  if (spnt->vcov <= 0) spnt->vcov = 0.0;
	
	  if (napnt->refr==NULLVAL)
	       spnt->vrefr = 0;
	  else if (stiminc) spnt->vrefr = int(napnt->refr/stiminc+0.5);
	  if (spnt->vrefr < 0) spnt->vrefr = 0;
	  if (spnt->vrefr > MAXVINT) spnt->vrefr = MAXVINT;

	  if (napnt->cov==NULLVAL)
	       spnt->cov = 1.0;
	  else spnt->cov = napnt->cov;
	  if (spnt->cov >  1) spnt->cov = 1.0;
	  if (spnt->cov <= 0) spnt->cov = 0.0;

	  maksseed (napnt, spnt);
 
	  break;

         }  /* switch */
      }  /* for (apnt;;) */

      spnt->sens = sens;
      spnt->vgain = vg;
      spnt->caegain = eg;
      spnt->ngain = sn;
      spnt->cgain = sg;
      spnt->coff = so;
      spnt->curve = (int)(sc);
      if (epnt->dyadelem > 0) {		/* if this synapse is 2nd in dyad */
        if (!(oepnt=(synapse*)findelemne(epnt->dyadelem))) {  /* find saved synap elem */
           if (info>=4) ncfprintf (stderr,"maksynap: can't find dyad element %d, using orig synapse.\n",epnt->dyadelem);
	   epnt->dyadelem = 0;				    /* if 1st synap doesn't exist */
           spnt->spost = NULL;
           spnt->spre  = NULL;
	   oepnt = epnt;
        }
         else { /* making a dyad synapse, connect back to original */
          if (oepnt->lptr) { 			  /* if 1st synap exists */ 
           makdyad((dyadlst **)&((synap *)oepnt->lptr)->spost,spnt);  
           makdyad((dyadlst **)&spnt->spre,(synap *)oepnt->lptr);
          }
        } 
      }
      else {	/* making a normal synapse, can connect dyad to it later */
         spnt->spost = NULL;
         spnt->spre  = NULL;
         oepnt = epnt;
      }
      if (epnt->spost!=NULLVAL) {	/* If this synapse is presynaptic */
					/*   to another postsyn elem */

        if (!(oepnt=(synapse*)findelemne((int)epnt->spost))){/* find saved elem */
           if (info>=4) ncfprintf (stderr,"maksynap: spost: can't find presyn element %d, using orig synapse.\n", (int)epnt->spost);
	   spnt->spost = NULL;
        }
        else {	/* making a dyad synapse, connect back to original */
          if (oepnt->lptr) { 			  /* if 1st synap exists */ 
           makdyad((dyadlst **)&((synap *)oepnt->lptr)->spre,spnt);  
           makdyad((dyadlst **)&spnt->spost,(synap *)oepnt->lptr);
          }
        }
      }
      else spnt->spost = NULL;

		/* make nt comps, connect to postsynaptic comp2 */

    if (epnt->mesg1!=NULLVAL) {
       spnt->mesg1 = epnt->mesg1;
       makntcomp (comp2,spnt->mesg1); 
    }
    if (epnt->mesg2!=NULLVAL) {
       spnt->mesg2 = epnt->mesg2;
       makntcomp (comp2,spnt->mesg2); 
       spnt->mesgconc = msgc;
    }


    if (tf2) tf2 = 1.0 - calctau(tf2,stiminc/MSEC);
    if (tf3) tf3 = 1.0 - calctau(tf3,stiminc/MSEC);
    
    spnt->filt1  = maksynfilt (LP,sf1, 0,    0,    0,oepnt->timec1, stiminc/MSEC);
    spnt->filt1h = maksynfilt (HP,sf1h,sf1ho,sf1hg,0,oepnt->timec1h,stiminc/MSEC);
    spnt->filt2  = maksynfilt (LP,sf2, 0,    tf2,  0, epnt->timec2, stiminc/MSEC);
    spnt->filt3  = maksynfilt (LP,sf3, 0,    tf3,  0, epnt->timec3, stiminc/MSEC);


/* ncfprintf (stderr,"sf1 %d sf1h %d ft1 %g sf2 %d ft2 %g sf3 %d ft3 %g\n",
				sf1,sf1h,ft1,sf2,ft2,sf3,tf3); /* */

/* ncfprintf (stderr,"nt %d se %g si %g sc %g st %g nkd %g mx %g rp %g\n",
				nt,se,si,sc,st,nkd,mx,rp); /* */

    return (spnt); 
}

/*------------------------------------*/

void setfiltval (lpfilt *fpnt, double val)

/* set a temporal filter for a synapse */

{
     int i;

  if (!fpnt) return;
  switch (fpnt->ctype) {
    case LPFILT:
        for (i=0; i<fpnt->nfilt; i++) {
           fpnt->lfilt[i] = val;	 /* init low-pass filt array */
        }
	break;

    case MFILT:  break;

  }
}

/*------------------------------------*/

void setfilttau (lpfilt *fpnt, int nfilt, double *timec, double tfall)

/* set a temporal filter for a synapse */

{
     int i;

  if (!fpnt) return;
  switch (fpnt->ctype) {
    case LPFILT:
      if (timec) {
        if (nfilt > fpnt->nfilt) {
          fpnt->ftau = makfiltarr(nfilt, fpnt->nfilt,fpnt->ftau,0.0);
        }
        fpnt->nfilt = nfilt;
        for (i=0; i<nfilt; i++) {
           fpnt->ftau[i] = calctau(timec[i],stiminc/MSEC); /* set time const */
        }
      }
      switch (fpnt->stype) {

           case LP: if (tfall!=NULLVAL) 
			fpnt->tfall = 1.0-calctau(tfall,stiminc/MSEC);
		    break;

           case HP: if (tfall != NULLVAL) fpnt->tfall = tfall;   
		    break;
      }
     break;

    case MFILT:  break;

  }
}

/*------------------------------------*/

void modfiltval (lpfilt *fpnt, int nfilt, double val)

/* change the number of filters in a synaptic temporal filter */

{
     double *lfilt;

  if (!fpnt) return;
  switch (fpnt->ctype) {
    case LPFILT:
	if (nfilt != NULLVAL) {
	  lfilt = fpnt->lfilt;
	  if (nfilt > fpnt->nfilt) {
             fpnt->lfilt = makfiltarr(nfilt, fpnt->nfilt,lfilt,val);
	  }
    	  fpnt->nfilt = nfilt;
	}
	break;

    case MFILT:  break;
  }
}

/*------------------------------------*/

void setsynv(synap *spnt, double v)
              
/* Set a synapse's internal activation from 
   its presynaptic voltage.
   If "synapse" is a voltage buffer, then
   set "destination" voltage equal to the
   "source" voltage.
 */

{
   int i;
   double cond,sa,st,tr,x;
   dbuf *dpnt;
   void maknois(chan *chpnt);
   double getcurv(synap *spnt, double val);

switch (spnt->ctype) {

case SYNAPSE:

if (makestim) {
    tr = 0;        			/* find resting transmitter rel */
    cond = 0;        			/* resting conductance */
}
else {
   st = spnt->thresh;			/* synaptic threshold */
   sa = v - st;				/* synaptic activation in mv */
   if (!oldsynapz)
     tr = getcurv(spnt,sa);        		/* find resting transmitter rel */
   else
    tr = sa;

   setfiltval(spnt->filt1, tr);
   setfiltval(spnt->filt1h,0);
   if (spnt->filt1h) spnt->filt1h->lfilt[0] = tr;   /* initialize first val */
   if (spnt->spre != NULL) {
         synap *tspnt;

       if (spnt->spre) {
         tspnt = spnt->spre->sdyad;
         if (tspnt) tr = getcurv (tspnt,v - tspnt->thresh);
       }
    }
   else {
     if (oldsynapz) tr = getcurv(spnt,sa);         /* find resting transmitter rel */
   }
   tr *= spnt->trconc;
   setfiltval(spnt->filt2,tr);		 /* initialize filter 2 */

   if (spnt->resp1) {
      switch (spnt->resp1->ctype){       /* initialize postsynaptic chan */
      
         case AMPA: 			/* set states & conductance */
         case KAINATE: 
         case NMDA: 
         case GABA: 
         case GLY: 
	   cond = dochani((chan*)spnt->resp1,1e-8);  
	   maknois((chan*)spnt->resp1);	/* set chans in states */
	   break;
      }
   }
   if (spnt->resp2) {
      switch (spnt->resp2->ctype){       /* initialize postsynaptic chan */
         case CGMP: 
         case SYN2:
           cond = setbind(spnt, tr);	/* calc binding */
           if (spnt->filt3) setfiltval(spnt->filt3,cond); /* init filter 3 */
           cond = setcond(spnt, cond);	/* calc static conductance */
	   cond = setbind2(spnt,cond);   /* calc 2nd mesg binding */
	   spnt->conduct = cond * spnt->maxcond;
	   cond = dochani((chan*)spnt->resp2,1e-8);  
	   maknois((chan*)spnt->resp2);	/* set chans in states */
	   ((chan*)spnt->resp2)->conduct = cond * spnt->maxcond;
	   break;

      }   /* switch (spnt->ctype) */
    }
    if (!(spnt->resp1 || spnt->resp2)) {      /* standard postsynaptic apparatus or mesg */
       cond = setbind(spnt, tr/spnt->trconc); /* calc binding */
       setfiltval(spnt->filt3,cond);	 /* initialize filter 3 */
       cond = setcond(spnt, cond);	 /* calc static conductance */
       cond = setbind2(spnt,cond);       /* calc 2nd mesg binding */
    }
}

    if (spnt->mesg1){
        x = addnt(spnt->comp2,spnt->mesg1, (cond-spnt->oldc1));
        spnt->oldc1 = cond;
    }
    if (spnt->mesg2){
        x = addnt(spnt->comp2,spnt->mesg2, (cond-spnt->oldc2)*spnt->mesgconc);
        spnt->oldc2 = cond;
    }
    break;

  case BUF:

    dpnt = (dbuf *)spnt;
    if (dpnt->delay) {
      if (dpnt->delbuf)
        for (i=0; i<dpnt->delay; i++) 
           dpnt->delbuf[i] = v;
      dpnt->delpnt = dpnt->delbuf;
    }
   else {
     dpnt->comp2->v = v; 
   }
   break;
  }  /* end switch */
}

/*------------------------------------*/

synap *modsynap(synapse *sepnt)
               
/* modify an old low-level synapse */

{
    int i, sfi;
    double dval;
    synapse *epnt;
    synap *spnt;
    attrib *apnt;
    nattrib *napnt;
    double getcurv(synap *spnt, double val);

    if (!(epnt=(synapse *)findelem(sepnt->modif))) {
      ncfprintf (stderr,"modsyn: can't find element %d\n",sepnt->modif);
      return ((synap*)NULL);
    }

    spnt = (synap *)epnt->lptr;		/* get pointer to synapse from elem */
    if (spnt==(synap*)NULL) return ((synap *)NULL); 
    if (spnt->ctype != SYNAPSE) {
       ncfprintf (stderr,"modsyn: element %d is not a synapse\n",sepnt->modif);
       return ((synap*)NULL);  
    }

 if (!makestim) {
   if (sepnt->ngain!=NULLVAL) {				/* neurotrans gain */
	epnt->ngain = spnt->ngain= sepnt->ngain;
   }
   if (sepnt->vgain!=NULLVAL) {				/* nt linear rel gain */
	epnt->vgain = spnt->vgain= sepnt->vgain; 
   }
   if (sepnt->curve!=NULLVAL) epnt->curve=spnt->curve = sepnt->curve; /* synap curv:lin,exp*/
   if (sepnt->cgain!=NULLVAL) {				/* cycG gain */
	epnt->cgain = spnt->cgain = sepnt->cgain;
   }
   if (sepnt->coff!=NULLVAL) {				/* cycG offset */
	epnt->coff = spnt->coff = sepnt->coff;
   }
   if (sepnt->ntact!=NULLVAL) epnt->ntact = spnt->ntact = sepnt->ntact; /* neurotrans action*/
   if (sepnt->nkd!=NULLVAL)   epnt->nkd   = spnt->nkd   = sepnt->nkd;   /* half-max sat point*/
   if (sepnt->npow!=NULLVAL)  epnt->npow  = spnt->npow  = sepnt->npow;  /* Hill coeff */
   if (sepnt->ckd!=NULLVAL)   epnt->ckd   = spnt->ckd   = sepnt->ckd;   /* half-max sat point*/
   if (sepnt->chc!=NULLVAL)   epnt->chc   = spnt->chc   = sepnt->chc;   /* cycG Hill coeff */
   if (sepnt->sens!=NULLVAL)  epnt->sens  = spnt->sens  = sepnt->sens;  /* presyn sensitivity */
   if (sepnt->thresh!=NULLVAL)  {
	epnt->thresh = spnt->thresh = sepnt->thresh;     		/* synaptic threshold */
   }
   if (sepnt->maxcond!=NULLVAL) {			  /* max conduct*/
	epnt->maxcond = spnt->maxcond = sepnt->maxcond; 
	if (spnt->resp1!=NULL) spnt->resp1->maxcond = sepnt->maxcond;
	if (spnt->resp2!=NULL) spnt->resp2->maxcond = sepnt->maxcond;
   }
   if (sepnt->vrev!=NULLVAL)  {				  /* reversal potential*/
	epnt->vrev = spnt->vrev = sepnt->vrev;  
   }
   if (sepnt->trconc!=NULLVAL) epnt->trconc = spnt->trconc = sepnt->trconc;
   if (sepnt->mesgconc!=NULLVAL) epnt->mesgconc = spnt->mesgconc = sepnt->mesgconc;
   if (sepnt->mesg1!=NULLVAL) epnt->mesg1 = spnt->mesg1 = sepnt->mesg1;
   if (sepnt->mesg2!=NULLVAL) epnt->mesg2 = spnt->mesg2 = sepnt->mesg2;
   if (sepnt->rrpool!=NULLVAL) epnt->rrpool = spnt->rrpool = sepnt->rrpool; /* readily releaseable pool*/
   if (sepnt->maxsrate!=NULLVAL) epnt->maxsrate = spnt->maxsrate = sepnt->maxsrate; /* max sus rate */

   if (sepnt->nfilt1 > 0) epnt->nfilt1  = sepnt->nfilt1;
   if (sepnt->nfilt1 > 0) epnt->nfilt1h = sepnt->nfilt1h;
   if (sepnt->nfilt2 > 0) epnt->nfilt2  = sepnt->nfilt2;
   if (sepnt->nfilt3 > 0) epnt->nfilt3  = sepnt->nfilt3;
   setfilttau(spnt->filt1, sepnt->nfilt1, sepnt->timec1,NULLVAL);
   setfilttau(spnt->filt1h,sepnt->nfilt1h,sepnt->timec1h,sepnt->filt1hg);
   setfilttau(spnt->filt2, sepnt->nfilt2, sepnt->timec2,sepnt->tfall2);
   setfilttau(spnt->filt3, sepnt->nfilt3, sepnt->timec3,sepnt->tfall3);

	/* If we're changing number of filters, */
        /*  we initialize them to the value of filter before them. */
        /* But if number of filters is 0, set them from comp */

   dval = spnt->comp1->v - spnt->comp1->vext - spnt->thresh;
   dval = getcurv(spnt, dval);
   modfiltval(spnt->filt1, sepnt->nfilt1, dval);
   modfiltval(spnt->filt1h, sepnt->nfilt1h, dval);
   modfiltval(spnt->filt2, sepnt->nfilt2, dval);

   if ((sfi=(sepnt->nfilt3)!=NULLVAL)) {	/* modify nfilt in third filter */
        double tr,bnd,cond;
        double k,kk,ll,pp;

      tr = dval;
      pp = spnt->npow;
      k = spnt->nkd;
      if (floor(pp) == pp) {         /*  integer power */
      	kk = k;
      	ll = tr;
              for (i=int(pp); i>1; i--) {
         	   kk *= k;
        	   ll *= tr;
              }
      }
      else {				/* non-integer power */
       	kk = pow(k,pp);
       	ll = pow(tr,pp);
      }
      bnd = ll / (ll + kk);		/* fraction bound */

      switch (spnt->ntact) {		/* find resting conductance */
      	case OPEN:
 		cond = bnd;
		break;
	case CLOSE:
		cond = (spnt->coff - bnd * spnt->cgain);
		break;
      }	
      modfiltval(spnt->filt3, sepnt->nfilt3, cond);
   }

	/* if access to old synapse's attribs elem is needed, */
	/* we must change this code to modify the old attrib */ 

    for (apnt=sepnt->attpnt; apnt; apnt=apnt->attpnt) {/* set up noise params */
      napnt = (nattrib *)apnt;
      switch (napnt->ctype) {
 	case NUNIT:
 	case CCHNOISE:
	break;

 	case VESNOISE:
	if (napnt->n!=NULLVAL)
	    spnt->vsites = (int)(napnt->n);	
	if (napnt->vsize!=NULLVAL)
	    spnt->vsize = napnt->vsize;	
        if (napnt->rseed > 0) initstate (napnt->rseed,spnt->vstate,RNDSIZ);
	break;

        }  /* switch */
      }  /* for (apnt;;) */
 }
 return (spnt); 
}

/*------------------------------------*/

void delsynap (synap *spnt)

{
  synap *tspnt;
  dyadlst *lpnt, *nlpnt;
  void deldyad(dyadlst **head);

switch (spnt->ctype) {

 case SYNAPSE:
  tspnt = spnt;
  if (spnt->vstate) efree (spnt->vstate);
  if (spnt->filt1)  efree (spnt->filt1);
  if (spnt->filt2)  efree (spnt->filt2);
  if (spnt->filt3)  efree (spnt->filt3);
  if (spnt->spost) {
   for (lpnt=spnt->spost; lpnt; lpnt=nlpnt) {
     nlpnt = lpnt->next;
     if (lpnt->sdyad) {
        deldyad(&lpnt->sdyad->spost);   /* delete pointer in other synapse */
     }
     efree (lpnt);
   }
  }
  if (spnt->spre) {
   for (lpnt=spnt->spre; lpnt; lpnt=nlpnt) {
     nlpnt = lpnt->next;
     if (lpnt->sdyad) {
        deldyad(&lpnt->sdyad->spre);   /* delete pointer in other synapse */
     }
     efree (lpnt);
   }
  }
  spnt = (synap*)spnt->next;
  efree (tspnt);
  cumsynap--;
  break;

  case GJ:
   efree (spnt);
  break;

  case BUF:
   if (((dbuf *)spnt)->filt!=NULL) efree (((dbuf*)spnt)->filt);
   efree (((dbuf *)spnt)->delbuf);
   efree (spnt);
  break;

 }
}

/*------------------------------------*/
recpar *rectypes[TOTREC];

void initrec(void)

/* initialize receptor transduction cascade constants */

/* Calibration procedure:

 Read:

    Tranchina, D. (1998) J. Gen. Physiol. 111: 3-6.

    Nikonov, S., Engheta, N., and Pugh, E.N. Jr. (1998) 
	J.  Gen. Physiol. 111: 7-37.

    Calvert, P.D., Ho, T.W., LeFebvre, Y.M, and Arshavsky, V.Y (1998)   
	J.  Gen. Physiol. 111: 39-51.


    1) adjust delay times for rhod, rhod2, rstar. Time constant
for rstar should be 0.5 msec @ 37 deg (Nikonov et al, 1998).

    2) adjust rise and fall of gpr (G), it is fast here, and its
time course (duration) sets rise time of pde.

    3) adjust rise and fall of pde (E). Its rise detemines rise
time of response. Its fall determines fall time constant of light response 
(Nikonov et al, 1998).

    4) adjust starting cycg, cond, ca, cax to give flat trace
with no response (tcomp21c).

    5) adjust sfactor (speed factor) to give peak at 200 msec.

    6) adjust cax rise and fall to give response proper asymmetry.

    7) adjust starting cycg, cond, ca, cax to give flat trace
with no response (tcomp21c).

    8) readjust sfactor (speed factor) to give peak at 200 msec.

    9) Adjust dmaxrod to give 34 pA max response.

    10) Adjust lgain to give single photon response of 0.7 pA. 
*/

{
    int i;
    recpar *cpt;
    double rrate, qrate;
    double cycgbuf, loopgain;
    int pigm;
    double sfactor;
    phvars *vpt;

    for (i=0; i<RODI; i++) {		/* allocate space for receptor constants */
       if ((rectypes[i]=(recpar*)emalloc(sizeof(recpar))) == (recpar*)NULL) {
           ncfprintf (stderr,"no space left for receptor const %d\n",i);
  	   return;
       }
    }
    if ((rectypes[RODI]=(recpar*)emalloc(sizeof(recpari))) == (recpar*)NULL) {      // invergo rod constants
        ncfprintf (stderr,"no space left for receptor const %d\n",RODI);
  	   return;
    }
    if ((rectypes[ChR2]=(recpar*)emalloc(sizeof(recparc))) == (recpar*)NULL) {      // invergo rod constants
        ncfprintf (stderr,"no space left for receptor const %d\n",ChR2);
  	   return;
    }

    cpt = rectypes[0];		/* set up rod cascade constants */

    cpt->ctype = ROD;
    cpt->stype = PORIG;	// original receptor type

    cpt->recslow = 1;		/* rod is slow by this factor */
    cpt->loopgain = 1;		/* gain for ca loop, sets stability */

    cpt->qc = dqrec;		/* Q10 factor for conductance shape */
    cpt->qcond = dqcrec;	/* Q10 factor for conductance cond */
				/* not currently used */

    cpt->qt = 37.0;		/* base temperature for channel rate & cond */
    cpt->oldtempcel = tempcel;
    qrate = exp(log(cpt->qc) * (tempcel - cpt->qt) / 10.0);

    sfactor = 1.0;		/* speedup of cascade & intensity */
    rrate = stiminc * qrate /sfactor;    /* stiminc = timinc for syn/stimuli */

    pigm = 0;
    cpt->pigm = pigm;
    cpt->pathl = pigmlen[pigm];	/* path length through o.s. */

    cycgbuf  = 1.0;
    loopgain = cpt->loopgain;

    cpt->powcag  = 2.0;
    cpt->powcak  = 2.0;
    cpt->powcycg = 3.0;

    cpt->pdebase = .04;
    cpt->kdgcyc  = .0575;	/* compensates for change from cax to ca -> gcyc */
    cpt->kdrhodk = .002;
    cpt->cat     = 100;
    cpt->totcond = 20.0;

    vpt = &cpt->vars;

    vpt->rhod  = 0.0;
    vpt->rhodk = 0.0;
    vpt->rhod2 = 0.0;
    vpt->rstar = 0.0;
    vpt->gpr   = 0.0;
    vpt->pde   = cpt->pdebase;
    vpt->resp_b  = 0;

    vpt->cycg  = 1.05605; 	/* values taken from tcomp21b */
    vpt->cond  = .994444; 
    vpt->ca    =  0.364626; 
    vpt->cax   =  0.152104; 
    vpt->cab   =  0.18*cpt->cat; 

    cpt->lgain    = rrate * 2.5;	/* adjust to give R* resp of 0.7 pA */
    cpt->rgain1   = calctau (0.0004,rrate);
    cpt->rgain2   = calctau (0.0004,rrate); /* 3 rhod steps give .5 ms tau */
    cpt->rgain3   = rrate * 800; 
    cpt->rkgain   = rrate * 5;	 /* adjust this for amount of adaptation */

    cpt->cn    = 4e-9/rrate; 	 // fwd rate for rod bleaching (cone: 4.1e-9 in coneset if tau_r=3.4 ms)
    cpt->rk_b  = 0.2;            // parameter of bleaching recovery
    cpt->tau_b = 100/rrate;      // time constant of rod bleaching recovery (cone: 10)

					     /*  for R* activation */
    cpt->ggain    = rrate * 8.0;
    cpt->gcygain  = rrate * 11.2 / cycgbuf * loopgain;
    cpt->pdegain  = rrate * 80.0 / cycgbuf;
    cpt->condgain = rrate * 448   * loopgain;
    cpt->gca      = rrate * 270.0 * loopgain;
    cpt->gcab     = rrate * 0.001 * loopgain;
    cpt->gcabr    = rrate * 0.001 * loopgain;
    cpt->cxgain   = rrate * 128.0 * loopgain;

    cpt->decrstar = 1 - (rrate * 20);
    cpt->pdedec   = 1 - (rrate * 16);	/* controls tail of cond */
    cpt->deccond  = 1 - (rrate * 480);
    cpt->capump   =      rrate * 1000;
    cpt->deccax   = 1 - (rrate * 320);

					/* voltage activation of dark curr */
					/* assume linear, but offset */
    cpt->vrev  =  0.72;			/* 0.1 activation / -80 mv */

    for (pigm=1; pigm<RODI; pigm++) {
      coneset(pigm);
    }
    rodiset(RODI);			/* set mouse rod, pigm == RODI, in rodi.cc */
    chr2set(ChR2);			/* set channel-rhodopsin, pigm=ChR2 in chr2.cc */
}

/*------------------------------------*/

void coneset(int pigm)
{
    recpar *cpt;
    double rrate, qrate, sfactor;
    phvars *vpt;

  if (pigm>=NUMREC) pigm = NUMREC-1;
  cpt = rectypes[pigm];		/* set up cone cascade constants */
				/* last pigm type (pigm=NUMREC) is ksens */
  cpt->ctype = CONE;

  if (pigm<CONEVH1) 
      cpt->stype = PORIG;	// original receptor types
  else if (pigm<RODI)
      cpt->stype = PVH;	        // van hateren receptor types

  cpt->pigm = pigm;
  cpt->recslow = 1.0;		/* cone is slow by this factor */
  cpt->loopgain = 1;		/* gain for ca loop, sets stability */

  cpt->qc = dqrec;		/* Q10 factor for conductance shape */
  cpt->qcond = dqcrec;	/* Q10 factor for conductance cond */
				/* not currently used */

  cpt->qt = 37.0;		/* base temperature for channel rate & cond */
  cpt->oldtempcel = tempcel;
  qrate = exp(log(cpt->qc) * (tempcel - cpt->qt) / 10.0);

  sfactor = 1.1;		/* speedup of cascade & intensity */
  rrate = stiminc * qrate /sfactor;    /* stiminc = timinc for syn/stimuli */

  cpt->pathl = pigmlen[pigm];	/* path length through o.s. */
 
  if (cpt->stype==PORIG) {	// original cone type
       double cycgbuf, loopgain;

    loopgain = cpt->loopgain;
    cycgbuf  = 1.0;

    cpt->powcag  = 2.3;
    cpt->powcak  = 2.0;
    cpt->powcycg = 3.0;

    cpt->pdebase = .04;
    cpt->kdgcyc  = .020;
    cpt->kdrhodk = .01;
    cpt->cat     = 250;
    cpt->totcond = 20.0;

    vpt = &cpt->vars;

    vpt->rhod  = 0.0;
    vpt->rhodk = 0.0;
    vpt->rhod2 = 0.0;
    vpt->rstar = 0.0;
    vpt->gpr   = 0.0;
    vpt->pde   = cpt->pdebase;
    vpt->resp_b  = 0;

    vpt->cycg  = 1.03805; 	/* values taken from tcomp21b */
    vpt->cond  = .902444; 
    vpt->ca    =  0.405; 
    vpt->cax   =  0.23119; 
    vpt->cab   =  9.7458; 

    cpt->lgain    = rrate * 20.0/loopgain;  /* adjust to give R* resp ~ 0.02 pA */
    cpt->rgain1   = calctau (0.0004,rrate);
    cpt->rgain2   = calctau (0.0004,rrate); /* 3 rhod steps give .5 ms tau */
    cpt->rgain3   = rrate * 800; 
    cpt->rkgain   = rrate * 20000;	/* adjust this for amount of adaptation */

    cpt->cn    = 4e-9/rrate; 	 // fwd rate for cone-bleaching (4.1e-9 for van Hat below if tau_r=3.4 ms)
    cpt->rk_b  = 0.2;            // parameter of bleaching recovery
    cpt->tau_b = 10/rrate;      // time constant of bleaching recovery

    cpt->ggain    = rrate * 40.0;
    cpt->gcygain  = rrate * 30.0  /cycgbuf * loopgain;  /* controls tail of cond */
    cpt->pdegain  = rrate * 100.0 /cycgbuf; 
    cpt->condgain = rrate * 448   * loopgain;
    cpt->gca      = rrate * 325.0 * loopgain;
    cpt->gcab     = rrate * 20  * loopgain;	/* phase delay for Ca++ feedback */
    cpt->gcabr    = rrate * 200 * loopgain;
    cpt->cxgain   = rrate * 2.0  * loopgain;	/* phase delay for rhodk */

    cpt->decrstar = 1 - (rrate * 20);
    cpt->pdedec   = 1 - (rrate * 100);
    cpt->deccond  = 1 - (rrate * 480);
    cpt->capump   =      rrate * 1000;
    cpt->deccax   = 1 - (rrate * 4);
					/* voltage activation of dark curr */
					/* assume linear, but offset */
    cpt->vrev  =  -0.008;		/* reversal pot for TRPM1 channel */
 }

//   orig:
//    sfactor = .5;		/* speedup of cascade & intensity */
//    cpt->kdrhodk = .002;
//    cpt->lgain    = rrate * .25;	/* adjust to give R* resp ~ 0.02 pA */
//    cpt->rkgain   = rrate * 50;	     /* adjust this for amount of adaptation */
//    cpt->gcab     = rrate * 1e-6 * loopgain;
//    cpt->gcabr    = rrate * 1e-6 * loopgain;
//    cpt->cxgain   = rrate * 5.0 * loopgain;  
//    cpt->deccax   = 1 - (rrate * 10);  
//
//    sfactor = 1.2;		/* speedup of cascade & intensity */
//    cpt->kdrhodk = .0005;
//    cpt->lgain    = rrate * .25*3;	/* adjust to give R* resp ~ 0.02 pA */
//    cpt->rkgain   = rrate * 8000;	     /* adjust this for amount of adaptation */
//    cpt->gcab     = rrate * 30  * loopgain;
//    cpt->gcabr    = rrate * 100 * loopgain;
//    cpt->cxgain   = rrate * 25.0 * loopgain;  
//    cpt->deccax   = 1 - (rrate * 50);  

  else if (cpt->stype==PVH) {	// van Hateren receptor type
       double tau=10.0;

   cpt->lgain  = rrate * 6e5; // rrate * 4.0;	// adjust to give R* resp ~ 0.02 pA
   cpt->tau_r=3.4*tau;          // R* lifetime (ms)
   cpt->cn=4.1e-9*cpt->tau_r/3.4; // norm const for cone-bleaching (4.1e-9 if tau_r=3.4 ms)
   cpt->tau_b=25000*tau;        // time constant of bleaching recovery
   cpt->tau_e=8.7*tau;          // E* lifetime
   cpt->tau_c=3*tau;            // time constant of Ca2+ extrusion
   cpt->tau_vc=10*tau;          // membrane time constant
   cpt->tau_is=90*tau;          // time constant of membrane nonlinearity
   cpt->rk_b=0.2;               // parameter of bleaching recovery
   cpt->c_beta=2.8e-3/tau;      // dark PDE activity (=1/tau_D)
   cpt->rk_beta=1.4e-4/tau;     // E* dependence of PDE activity
   cpt->beta_e_max=4/tau;       // parameter diffusion-limited cGMP hydrolysis
   cpt->rnx=1;                  // apparent Hill coefficient CNG channels
   cpt->rnc=4;                  // Hill coefficient GC activation
   cpt->a_c=0.23;               // scaling constant of GC activation
   cpt->gamma_is=0.7;           // parameter of membrane nonlinearity
   cpt->a_is=2.9e-2;            // parameter of membrane nonlinearity
   cpt->vrev  =  -0.008;	// reversal pot for TRPM1 channel

  }
}

/*------------------------------------*/

static recpar *rpt;	// pointer to van hateren cone params for steady_ode()
static double stimh;	// stimulus for equilib

// steady_ode
//
// finds steady-state of van hateren cone
//  returns 0 if x=steady-state of resp_c
//
//     Function accompanying 'Simulating human cones from mid-mesopic
//     up to high-photopic luminances' (J.H. van Hateren, 2006)
//
//	 used to initialize van Hateren cones
//
//     J.H. van Hateren, 28/8/06
//     translated into C by R.G Smith, 9/2011

  double steady_ode (double x)
  
{
     double rval, beta, beta_e, stimn;
     double resp_r, resp_b, resp_e;
     phvars *vpt;

  vpt = &rpt->vars;
  if (stimh==0) vpt->resp_e=0;
  else {
    stimn=stimh*rpt->cn;
    resp_b=0.5 * (1-rpt->rk_b-rpt->tau_r/rpt->tau_b*rpt->rk_b*(stimn+1)/stimn +
         sqrt(pow(1-rpt->rk_b-rpt->tau_r/rpt->tau_b*rpt->rk_b*(stimn+1)/stimn,2)+4*rpt->rk_b));
    resp_r=(1-resp_b)*stimh/(1+stimn);
    resp_e=resp_r;
  }
  beta=rpt->c_beta+rpt->rk_beta*resp_e;
  beta_e=beta/(1+beta/rpt->beta_e_max);
  rval = x - pow((1/(1+pow(rpt->a_c*x,rpt->rnc)))/beta_e,rpt->rnx);
  return rval;
}
  
/*------------------------------------*/

void setdark_vhcone (recpar *rpnt)

/* Adapt a van hateren cone type to dark */
/* Call to get dark response for normalizing response output */ 

{
    double beta, beta_e;
    double resp_q, gain_x; 
    phvars *vpnt;

  rpt = rpnt;		// external variables for steady_ode called by fzero()
  stimh = 0;

  vpnt = &rpnt->vars;
  vpnt->resp_c = fzero(steady_ode,1);
  // printf ("# dark resp_c %g\n",resp_c);

  beta    = rpnt->c_beta;
  beta_e  = beta / (1+beta/rpnt->beta_e_max);
  resp_q  = 1/beta_e;
  gain_x  = 1/(1 + pow(rpnt->a_c*vpnt->resp_c,rpnt->rnc));
  vpnt->resp_x  = gain_x*resp_q;
  rpnt->dark_resp_os = pow(vpnt->resp_x,rpnt->rnx);
  rpnt->dark_resp_is= pow(rpnt->dark_resp_os/rpnt->a_is,(1/(1+rpnt->gamma_is)));
}

/*------------------------------------*/

void adapt_vhcone (photrec *rpnt, double stim) 

/* adapt a van hateren cone type to a given intensity */

// see cone_ode.cc 
//
//  van Hateren JH, Snippe HP (2006) Simulating human cones from
//     mid-mesopic up to high-photopic luminances. Journal of Vision 7(4):1-11.

{
    double beta, beta_e, resp_q, gain_x, resp_im, stimn;
    phvars *vpnt;

  rpt = rpnt->chtyp;		// external variables for steady_ode called by fzero()
  stimh = stim;
  vpnt = &rpnt->vars;

  vpnt->resp_c = fzero(steady_ode,1);

  if (stim==0) {
    vpnt->resp_b  = 0;
    vpnt->resp_r  = 0;
    vpnt->resp_e  = 0;
  }
  else {
    stimn  = stim * rpt->cn;
    vpnt->resp_b = 0.5 * (1-rpt->rk_b-rpt->tau_r/rpt->tau_b*rpt->rk_b*(stimn+1)/stimn +
                 sqrt(pow(1-rpt->rk_b-rpt->tau_r/rpt->tau_b*rpt->rk_b*(stimn+1)/stimn,2)+
			 4*rpt->rk_b));
    vpnt->resp_r = (1-vpnt->resp_b)*stim/(1+stimn);
    vpnt->resp_e = vpnt->resp_r;
  }
  beta    = rpt->c_beta + rpt->rk_beta*vpnt->resp_e;
  beta_e  = beta/(1 + beta/rpt->beta_e_max);
  resp_q  = 1/beta_e;
  gain_x  = 1/(1+pow(rpt->a_c*vpnt->resp_c,rpt->rnc));
  vpnt->resp_x  = gain_x*resp_q;
  vpnt->resp_os = pow(vpnt->resp_x,rpt->rnx);
  vpnt->resp_is = pow(vpnt->resp_os/rpt->a_is,(1/(1+rpt->gamma_is)));
  resp_im = pow(vpnt->resp_is,rpt->gamma_is);
  vpnt->atten_i = rpt->a_is*resp_im;
}

/*------------------------------------*/

photrec *makr(photorec *epnt, int ctype)

/* make a dummy receptor and link it to the list. */

{
    photrec *rpnt;
    int i, errflag=0;

  if ((rpnt=(photrec *)emalloc(sizeof(photrec))) == (photrec*)NULL) {
     ncfprintf (stderr,"no space left for receptor %d\n",cumphotrec+1);
      return ((photrec*)NULL);  
  }
  if ((rpnt->aflux=(double *)emalloc(numstimchan*sizeof(double))) == (double*)NULL) errflag=1;
  if ((rpnt->mflux=(double *)emalloc(numstimchan*sizeof(double))) == (double*)NULL) errflag=1;
  if ((rpnt->mask =(double *)emalloc(numstimchan*sizeof(double))) == (double*)NULL) errflag=1;
  if ((rpnt->chanw =(double *)emalloc(numstimchan*sizeof(double))) == (double*)NULL) errflag=1;
  if (errflag) {
     ncfprintf (stderr,"no space left for photoreceptor %d\n",cumphotrec+1);
     return ((photrec*)NULL);  
  }
  rpnt->next = (photrec*)NULL;
  if (!recpnt) recpnt = rpnt;	  	/* save head if first photrec */
  if (recend)
    recend->next = rpnt;
  rpnt->last = recend;
  recend = rpnt;

  rpnt->ctype = ctype;			/* set channel to dummy type */
  rpnt->xloc  = epnt->xpos;		/* x loc of receptor */
  rpnt->yloc  = epnt->ypos;		/* y loc of receptor */
  rpnt->recnm1 = epnt->node1a;		/* receptor number from node */
  rpnt->recnm2 = epnt->node1b;
  rpnt->recnm3 = epnt->node1c;
  rpnt->recnm4 = epnt->node1d;
  for (i=0; i<numstimchan; i++) {
      rpnt->aflux[i] = 0.0;
      rpnt->mflux[i] = 0.0;
      rpnt->mask[i]  = 0.0;
      rpnt->chanw[i]  = 0.0;
  }
  rpnt->stimchan = 0;			/* can set to give private stimulus channel */
  rpnt->chanw[0] = 1.0;
  rpnt->iflux = 0.0;
  cumphotrec++;			/* increment total photoreceptors */
  return (rpnt); 
}

/*---------------------------------------------------*/

photrec *makphotrec(photorec *epnt, comp *cpnt)
               
/* make a new photoreceptor and link it to the list. */

{
    int i,pigm,filt,cnois,pnois, tspeed, errflag=0;;
    int nrseed;
    double dnois,unitary;
    double linit,sens;
    photrec *rpnt;
    photreci *rpnti;
    photrecc *rpntc;
    recpar *cpt;
    double maxcond, dia, pathl, attf;
    double loopgain, recslow, timec, timefactor;
    phvars *vpnt;
    char sbuf[100];

  if (epnt->pigm==RODI) {			// Invergo et al (2014) mouse rod
    if ((rpnti=(photreci *)emalloc(sizeof(photreci))) == (photreci*)NULL) {
       ncfprintf (stderr,"no space left for receptor %d\n",cumphotrec+1);
        return ((photrec*)NULL);  
    }
    rpnt = (photrec*)rpnti;
  } else if (epnt->pigm==ChR2){			// Williams et al (2013) channel rhodopsin
    if ((rpntc=(photrecc *)emalloc(sizeof(photrecc))) == (photrecc*)NULL) {
       ncfprintf (stderr,"no space left for receptor %d\n",cumphotrec+1);
        return ((photrec*)NULL);  
    }
    rpnt = (photrec*)rpntc;
  } else {
    if ((rpnt=(photrec *)emalloc(sizeof(photrec))) == (photrec*)NULL) {
       ncfprintf (stderr,"no space left for receptor %d\n",cumphotrec+1);
        return ((photrec*)NULL);  
    }
  }
  if ((rpnt->aflux=(double *)emalloc(numstimchan*sizeof(double))) == (double*)NULL) errflag=1;
  if ((rpnt->mflux=(double *)emalloc(numstimchan*sizeof(double))) == (double*)NULL) errflag=1;
  if ((rpnt->mask =(double *)emalloc(numstimchan*sizeof(double))) == (double*)NULL) errflag=1;
  if ((rpnt->chanw=(double *)emalloc(numstimchan*sizeof(double))) == (double*)NULL) errflag=1;
  if (errflag) {
     ncfprintf (stderr,"no space left for photoreceptor %d\n",cumphotrec+1);
     return ((photrec*)NULL);  
  }
  rpnt->next = (photrec*)NULL;
  if (!recpnt) recpnt = rpnt;  	/* save head if first recep */
  if (recend)
    recend->next = rpnt;
  rpnt->last = recend;
  recend = rpnt;
  cumphotrec++;	 			/* increment total receptors */
  rpnt->comp1 = cpnt; 			/* no presynaptic cell */
  maklst(&cpnt->clst,(conn*)rpnt);	/* compartment assoc with recep */
  rpnt->ctype = epnt->ctype;		/* set channel type: ROD or CONE */
  rpnt->xloc  = epnt->xpos;		/* x loc of receptor */
  rpnt->yloc  = epnt->ypos;		/* y loc of receptor */

  vpnt = &rpnt->vars;

  if (epnt->ctype == ROD) {		/* rod */

    if ((maxcond=epnt->maxcond)==NULLVAL) maxcond=dmaxrod;/* max conductance*/
    if ((dia=epnt->dia)==NULLVAL) dia = 1.674; /* 2.2 um2 area */
    epnt->dia = dia;
    rpnt->area = dia * dia * MPI / 4;	/* photon collecting area of rod */
    rpnt->maxcond = maxcond;		/* maximum conductance of channel */
    if ((attf=epnt->attf)==NULLVAL) attf = .9;  /* misc attenuation factor */

    if ((pigm=(int)(epnt->pigm))==NULLVAL) pigm = 0; /* rod default pigment */
    if (pigm>=NUMREC) pigm = NUMREC-1;
    else if (pigm<0)  pigm = 0;

    rpnt->chtyp = rectypes[pigm];		/* rod constants */
    if (rpnt->chtyp->stype == PORIG) {
	rpnt->stype = PORIG;
        vpnt->ksens = rpnt->chtyp->vars.ksens = 500;
    } else if (rpnt->chtyp->stype==PINVG) {
	rpnt->stype = PINVG;		/* invergo rod type */
    } else {
	ncfprintf (stderr,"makphotrec: wrong pigm type %d for rod\n",pigm);
	return (NULL);
    }
  }

  else if (epnt->ctype == CONE || epnt->ctype == CHR ||
		  epnt->ctype==VTRANSDUCER ||
		  epnt->ctype==ITRANSDUCER) {	/* cone */

    if ((pigm=(int)(epnt->pigm))==NULLVAL) {
	if (epnt->ctype==CONE) pigm = 1; /* cone default pigment */
	if (epnt->ctype==CHR) pigm = ChR2; /* CHR default pigment */
    }
    if (pigm>=NUMREC) pigm = NUMREC-1;
    else if (pigm<1)  pigm = 1;

    rpnt->chtyp = rectypes[pigm];		/* cone constants */
    if (pigm>=19 && pigm<=21) 
         rpnt->stype = PVH;			/* set VH type */ 
    else if (pigm==0 || pigm==RODI) {
	if      (epnt->ctype==CONE) 
		sprintf (sbuf,"makphotrec: wrong pigm type %d for cone\n",pigm);
	else if (epnt->ctype==CHR)
		sprintf (sbuf,"makphotrec: wrong pigm type %d for ChR\n",pigm);
	execerror(sbuf,"");
	return (NULL);
    }
    else if (pigm==ChR2) {
         rpnt->stype = ChR2;
    }
    else rpnt->stype = PORIG;
    if (epnt->ctype==CONE && pigm==ChR2) {
	sprintf (sbuf,"makphotrec: wrong pigm type %d for cone\n",pigm);
	execerror(sbuf,"");
	return (NULL);
    }
    if (epnt->ctype==CHR && pigm!=ChR2) {
	sprintf (sbuf,"makphotrec: wrong pigm type %d for ChR\n",pigm);
	execerror(sbuf,"");
	return (NULL);
    }

    if ((maxcond=epnt->maxcond)==NULLVAL) maxcond=dmaxcon; /* max conductance */
    if ((dia=epnt->dia)==NULLVAL) dia = 1.674; /* 2.2 um2 area */
    epnt->dia = dia;
    rpnt->area = dia * dia * MPI / 4;	/* photon collecting area of cone */
    rpnt->maxcond = maxcond;		/* maximum conductance of channel */
    if ((attf=epnt->attf)==NULLVAL) attf = .9;  /* misc attenuation factor */
    vpnt->ksens = rpnt->chtyp->vars.ksens = 50000.0;
    if ((unitary=epnt->unitary)==NULLVAL) unitary=dpru;		/* unitary conductance */

    if (epnt->ctype==CHR) {		/* make channel for ChR2 */
	    chattrib *apnt;
	    chan *chpnt;
	    int chr2_channel = 1;	/* set to 0 to use runchr() in chr2.cc */
					/* test with tcomp23c */
	if (chr2_channel > 0) {
	   apnt = make_chan(epnt,CHRC,ChR2);
           apnt->maxcond = maxcond;	
	   if (epnt->channoise>0) {
               apnt->attpnt = maknattr();
	       apnt->attpnt->ctype=CCHNOISE;
	       if (unitary<=0) unitary = 2-12;
	       ((nattrib*)(apnt->attpnt))->unitary=unitary;
	       ((nattrib*)(apnt->attpnt))->n = maxcond/unitary;
	   }
	   if ((chpnt=makchan(apnt, cpnt, NULL, epnt->nocondens))==NULL) 
		fprintf (stderr,"# makphotrec, can't make ChR2 channel\n"); 

 	   rpnt->phchan = chpnt;	/* use channel defined in chanchr2() */
	} else
	   rpnt->phchan = NULL;	/* don't use ChR2 channel, use runchr() in chr2.cc */
    } else {
	rpnt->phchan = NULL;	/* no channel for rods,cones */
    }
  }

  for (i=0; i<numstimchan; i++) {
        rpnt->aflux[i] = 0.0;
        rpnt->mflux[i] = 0.0;
        rpnt->mask[i]  = 0.0;
        rpnt->chanw[i] = 0.0;
  }
  rpnt->chanw[0]  = 0.0;
  if (epnt->stimchan==NULLVAL) rpnt->stimchan = 0;	/* set to give private stimulus channel */
  else rpnt->stimchan = epnt->stimchan;	
  rpnt->iflux = 0.0;

  cpt = rpnt->chtyp;			/* pointer to receptor constants */

  if (pigm>=19 && pigm<=21) setdark_vhcone(cpt);  /* set dark resting output level */
  restorphotrec(rpnt);			/* set receptor state starting vals */

  if ((pathl=epnt->pathl)==NULLVAL) pathl=cpt->pathl; /*default path length*/
  rpnt->pathl = pathl;
  if ((filt=epnt->filt)==NULLVAL) filt = 0; 	/* filter over pigment */
  rpnt->filt = filt;
  
  if ((timec=epnt->timec1)==NULLVAL) timec = 1.0; /* time speed factor */
  if (timec > 10.0) timec = 10.0;
  if (timec <= 0.0) timec = 0.001;
  recslow = cpt->recslow * stiminc/HMICSEC;
  if (timec > recslow) {
    tspeed = 1;				/* must have at least 1 iteration */
    timefactor = recslow / timec;	/* but can slow it down a lot */
  } 
  else {
    tspeed = (int) (recslow / timec);		/* number of time steps */
    timefactor = recslow / (tspeed * timec);
  }
  rpnt->tspeed = tspeed; 	/* number of time steps to run photrec */
  rpnt->intenadj = 1.0 / (tspeed * timefactor); /* norm intensity */

  if ((loopgain=epnt->loopg)==NULLVAL) loopgain = 1.0; /* gain for ca loop */
  if (loopgain < 0.2) loopgain = 0.2;
  if (loopgain > 1.0) loopgain = 1.0;
  rpnt->loopgain = loopgain;

  if (timefactor != 1.0 || loopgain != 1.0) {
    savephotrec(rpnt);		/* if timing is nonstandard, save chtyp */
    setphotrectim(rpnt,timefactor,loopgain);
  }

#define IMASK 0x7fffffff
  if ((dnois=epnt->darknoise)==NULLVAL) dnois = 0;	/* dark noise */
  if (dnois>0) {
     rpnt->dnois = dnois;
     rpnt->dstate = (char*)NULL;
     if ((nrseed=epnt->dkrseed) >= 0) {		/* use default random if neg */
        if (nrseed==0) /* use nodenum if rseed not set */
	  nrseed = noduniq(epnt->node1a, epnt->node1b,
			 epnt->node1c, epnt->node1d, srseed^rseed^dkrseed^1111);
        rpnt->dstate  = makrand (RNDSIZ,"dark",cumphotrec);
        if(rpnt->dstate)
            initstate(((nrseed&IMASK)+1818213)&IMASK,rpnt->dstate,RNDSIZ);
     }
  }
  if ((pnois=epnt->photnoise)==NULLVAL) pnois = 0;	/* photon noise */
  if (pnois) {
     rpnt->pnois = pnois;
     rpnt->pstate = (char*)NULL;
     if ((nrseed=epnt->phrseed) >= 0) {		/* use default random if neg */
        if (nrseed==0) /* use nodenum if rseed not set */
	  nrseed = noduniq(epnt->node1a, epnt->node1b,
			   epnt->node1c, epnt->node1d, srseed^rseed^phrseed^2222);
        rpnt->pstate = makrand (RNDSIZ,"phot",cumphotrec);
        if(rpnt->pstate)
            initstate(((nrseed&IMASK)+1098667)&IMASK,rpnt->pstate,RNDSIZ);
#undef IMASK
     }
  }
  rpnt->attf = attf;

  if ((linit=epnt->linit)==NULLVAL) linit = 0.0;  /* initial light value */
  rpnt->iflux = 0.0;

  if (!makestim) {

#define NITER 10000

  if (linit >= 0) { 		/* initialize photoreceptor activity */
       int i;
       double cond, ocond, err;

  if (epnt->ctype == CONE || epnt->ctype==ROD ||	/* don't initialize vclamp or iclamp */ 
		  	     epnt->ctype==CHR) {
     sens = rsens(rpnt,1,rpnt->filt);
     rpnt->aflux[0] = rpnt->area * linit * sens * rpnt->attf * dqeff * stiminc;
     rpnt->chanw[0] = 1;

     switch (pigm) {			/* initialize van hateren cones */
	case 19:
	case 20:
	case 21:
		adapt_vhcone (rpnt, rpnt->aflux[0]);
		break;
     }
  }
  // fprintf (stderr,"RK %g %g\n",rpnti->vars.RK,rpnti->chtyp->vars.RK);

  if (epnt->ctype == CONE || epnt->ctype==ROD ||	/* don't equilibrate vclamp or iclamp */
		  	     epnt->ctype==CHR) {
	  double dnois;
	  int equil;

     dnois = rpnt->dnois;
     rpnt->dnois = 0;
     runrec(rpnt,50);
     ocond = vpnt->cond;
     for (equil=i=0; i<NITER; i++) {		/* equilibrate photoreceptor */
	 runrec(rpnt,20);		/*  until convergence */
         cond = vpnt->cond;
	 err = cond - ocond;
	 err = ncabs(err);
         if (err<1e-7) {
           //  fprintf (stderr,"equil_phtorec: i %d err %g cond %g\n",i,err,cond);
	    equil = 1;
	    break;
         }
         ocond = cond; 
         // fprintf (stderr,"equil_phtorec: i %d err %g\n",i,err);
     }
     // if (!equil) fprintf (stderr,"# equil_phtorec: no equil i %d err %g cond %g\n",i,err,cond);
     // rpnt->aflux[0] = 0.0;
     rpnt->dnois = dnois;
   }
#undef NITER
  }
 }

  rpnt->recnm1 = epnt->node1a;		/* receptor number from node */
  rpnt->recnm2 = epnt->node1b;
  rpnt->recnm3 = epnt->node1c;
  rpnt->recnm4 = epnt->node1d;
  return (rpnt); 
}

/*------------------------------------*/

void savephotrec(photrec *rpnt)

/* Save a photoreceptor's kinetic states for a later restore */

{
    recpar *rp,*cpt;
    recpari *rpi;
    photreci *rpnti;
    char *emalloc(unsigned int n);

  cpt = rpnt->chtyp;			/* pointer to receptor constants */
  if (cpt== rectypes[cpt->pigm]) {      /* if orig values, make space */

    if (cpt->pigm==RODI) { 		/* invergo mouse rod */
      if ((rpi=(recpari *)emalloc(sizeof(recpari))) == (recpari *)NULL) {
         ncfprintf (stderr,"no space left for receptor save\n");
         return;  
      }
      rpnti = (photreci*)rpnt;
      *rpi = *(rpnti->chtyp);		/* copy the old constants */

    } else {				/* orig rod, cones */
      if ((rp=(recpar *)emalloc(sizeof(recpar))) == (recpar *)NULL) {
         ncfprintf (stderr,"no space left for receptor save\n");
         return;  
      }
      *rp = *cpt;			/* copy the old constants */
    }
  }
  else {				/* else, use the existing space */
     if (cpt->pigm==RODI) rpi = (recpari*)cpt;
     else                  rp = cpt;
  }

 if (cpt->pigm < RODI) {
   *(&rp->vars) = *(&rpnt->vars);       /* set the new values */
    rpnt->chtyp = rp;			/* pointer to receptor variables */
 } else {
   *(&rpi->vars) = *(&rpnti->vars);     /* set the new values */
   rpnti->chtyp = rpi;			/* save pointer to receptor variables */
 }
}

/*------------------------------------*/

void restorphotrec(photrec *rpnt)

/* restore photoreceptor kinetic states */

{
    recpar *cpt;
    recpari *cpti;
    photreci *rpnti;
    recparc *cptc;
    photrecc *rpntc;

  if (rpnt->chtyp->stype < PINVG) {     /* if original type */
    cpt = rpnt->chtyp;			/* pointer to receptor constants */
    *(&rpnt->vars) = *(&cpt->vars);	/* copy the receptor state variables */

  } else if (rpnt->chtyp->stype == PINVG) { 
    rpnti = (photreci*)rpnt;		/* invergo mouse rod type */
    cpti = rpnti->chtyp;		/* pointer to receptor constants */
    *(&rpnti->vars)  = *(&cpti->vars);	/* copy the receptor state variables */

  } else if (rpnt->chtyp->stype == ChR2) { 
    rpntc = (photrecc*)rpnt;		/* invergo mouse rod type */
    cptc = rpntc->chtyp;		/* pointer to receptor constants */
    *(&rpntc->vars)  = *(&cptc->vars);	/* copy the receptor state variables */
  }
}

/*------------------------------------*/

void setphotrectim(photrec *rpnt, double timec, double loopgain)
                 

/* set time constant of photoreceptor kinetic states */

{
    recpar *cpt;

  cpt = rpnt->chtyp;			/* pointer to receptor constants */
  if (cpt== rectypes[cpt->pigm]) { /* if orig estimates, make space */
     ncfprintf (stderr,"Error, can't change original receptor constants\n");
     return;
  }
   cpt->loopgain = loopgain;

   cpt->lgain    *= timec / loopgain;
   cpt->rgain1   *= timec;
   cpt->rgain2   *= timec;
   cpt->rgain3   *= timec;
   cpt->rkgain   *= timec;
   cpt->ggain    *= timec; 
   cpt->gcygain  *= timec * loopgain;
   cpt->pdegain  *= timec;
   cpt->condgain *= timec * loopgain;
   cpt->gca      *= timec * loopgain;
   cpt->gcab     *= timec * loopgain;
   cpt->gcabr    *= timec * loopgain;
   cpt->cxgain   *= timec * loopgain;

   cpt->decrstar = 1.0 - ((1.0 - cpt->decrstar) * timec);
   cpt->pdedec   = 1.0 - ((1.0 - cpt->pdedec) * timec);
   cpt->deccond  = 1.0 - ((1.0 - cpt->deccond) * timec);
   cpt->capump  *=    timec; 
   /* cpt->capump =  1.0 - ((1.0 - cpt->capump) * timec);   ??  new */
   cpt->deccax   = 1.0 - ((1.0 - cpt->deccax) * timec);

   cpt->tau_r 	*= timec;      // R* lifetime (ms)
   cpt->cn 	*= timec;      // normalization constant for cone- bleaching (4.1e-9 if tau_r=3.4 ms)
   cpt->tau_b 	*= timec;      // time constant of bleaching recovery
   cpt->tau_e 	*= timec;      // E* lifetime
   cpt->tau_c 	*= timec;      // time constant of Ca2+ extrusion
   cpt->tau_vc 	*= timec;      // membrane time constant
   cpt->tau_is 	*= timec;      // time constant of membrane nonlinearity
   cpt->rk_b 	*= timec;      // parameter of bleaching recovery
   cpt->c_beta 	*= timec;      // dark PDE activity (=1/tau_D)
   cpt->rk_beta *= timec;      // E* dependence of PDE activity
   cpt->beta_e_max *= timec;   // parameter diffusion-limited cGMP hydrolysis

}

/*------------------------------------*/

photrec *modphotrec(photorec *pepnt)
               
/* modify a photoreceptor */

{
    photrec *rpnt;
    photorec *epnt;
    recpar *cpt;
    double maxcond, dia, attf, pathl;
    int filt, tspeed;
    double dnois;
    double loopgain, recslow, timec, timefactor;

   if (!(epnt=(photorec*)findelem(pepnt->modif))) {
      ncfprintf (stderr,"modphotrec: can't find element %d\n",pepnt->modif);
      return ((photrec*)NULL);
   }

   rpnt = (photrec *)epnt->lptr;	/* get pointer to recept from elem */
   if (rpnt==(photrec*)NULL) return ((photrec*)NULL); 
   if (rpnt->ctype != ROD && rpnt->ctype != CONE && rpnt->ctype != CHR) {
    ncfprintf (stderr,"modphotrec: element %d is not a photoreceptor\n",pepnt->modif);
       return ((photrec*)NULL);  
   }
					/* don't want to change type or loc: */

/*  rpnt->ctype = pepnt->ctype;		/* set channel type: ROD or CONE */
/*  rpnt->xloc  = pepnt->xpos;		/* x loc of receptor */
/*  rpnt->yloc  = pepnt->ypos;		/* y loc of receptor */

  if ((maxcond=pepnt->maxcond)!=NULLVAL)	/* maximum conductance */
     epnt->maxcond = rpnt->maxcond = maxcond;	/* maximum conductance of channel */
  if ((dia=pepnt->dia)!=NULLVAL) {
     epnt->dia = pepnt->dia;
     rpnt->area = dia * dia * MPI / 4;		/* photon collecting area of rod */
  }
  if ((attf=pepnt->attf)!=NULLVAL)		/* misc attenuation factor */
     epnt->attf = rpnt->attf = attf;
  if ((pathl=pepnt->pathl)!=NULLVAL)		/* path length */
      epnt->pathl = rpnt->pathl = pathl;
  if ((filt=pepnt->filt)!=NULLVAL)		/* filter over pigment */
      epnt->filt = rpnt->filt = filt;
  if ((dnois=pepnt->darknoise)!=NULLVAL)	/* dark noise */
      epnt->darknoise = rpnt->dnois = dnois;
/*  if ((pnois=pepnt->photnoise)!=NULLVAL)	/* photon noise */
/*      epnt->photnoise = rpnt->pnois = pnois;   */       /* never modify photnoise */

  if ((timec=pepnt->timec1)!=NULLVAL) {  /* time speed factor */
    if (timec > 10.0) timec = 10.0;
    if (timec <= 0.0) timec = 0.001;
    epnt->timec1 = timec;
    cpt = rpnt->chtyp;
    recslow = cpt->recslow;
    if (timec > recslow) {
      tspeed = 1;			/* must have at least 1 iteration */
      timefactor = recslow / timec;	/* but can slow it down a lot */
    } 
    else {
      tspeed = (int) (recslow / timec);		/* number of time steps */
      timefactor = recslow / (tspeed * timec);
    }
    rpnt->tspeed = tspeed; 	/* number of time steps to run photrec */
    rpnt->intenadj = 1.0 / (tspeed * timefactor); /* norm intensity */
  }
  else timefactor = 1.0;

  if ((loopgain=epnt->loopg)!=NULLVAL) {   /* gain for ca loop */
    if (loopgain < 0.2) loopgain = 0.2;
    if (loopgain > 1.0) loopgain = 1.0;
    rpnt->loopgain = loopgain;
  }
  else loopgain = 1.0;

  if (timefactor != 1.0 || loopgain != 1.0) {
    savephotrec(rpnt);	/* if timec is nonstandard, save chtyp */
    setphotrectim(rpnt,timefactor,loopgain);
  }

  if (pepnt->save) {			/* save equilibrium kinetic states */
     savephotrec(rpnt);
  }
  else 
    if (pepnt->restore) {		/* restore kinetic states */
     restorphotrec(rpnt);
  }

  return (rpnt);
}

/*------------------------------------*/

void delphotrec (photrec *rpnt)

{
     recpar *cpt;
     photrec *next;
     photrec *last;

  if (!rpnt) return;
  if ((cpt=rpnt->chtyp))		  /* pointer to receptor constants */
    if (cpt!= rectypes[cpt->pigm]) {     /* if saved, delete new chtyp */
     efree (cpt);
    }
  if (rpnt->dstate) efree (rpnt->dstate);
  if (rpnt->pstate) efree (rpnt->pstate);
  if (rpnt->aflux) efree (rpnt->aflux);
  if (rpnt->mflux) efree (rpnt->mflux);
  if (rpnt->mask)  efree (rpnt->mask);
  if (rpnt->chanw) efree (rpnt->chanw);
  next = (photrec *)rpnt->next;
  last = (photrec *)rpnt->last;
  if (next) next->last = last;
  if (last) last->next = next; 
  if (rpnt==recpnt) recpnt = next;
  if (rpnt==recend) recend = last;
  efree (rpnt);
  cumphotrec--;
}

/*------------------------------------*/

conn *makconn(comp *comp1, comp *comp2, double s, int type)

/* make a new connection betw compartments and link
    it to both compartments' conn lists. */
/* connect comp1 to comp2 with conductance s; */

{
    conn *cpnt;

#ifdef DEBUG
 if (debug & NCMAK && debugz & 16)
   ncfprintf (stderr,"makconn: comp %d comp %d\n",comp1->num, comp2->num);
#endif

    if ((cpnt=(conn *)emalloc(sizeof(conn))) == (conn*)NULL) {
      ncfprintf (stderr,"no space left for conn %d\n",cumconn+1);
      return ((conn*)NULL);  
    }
					/* add gj or op amp to synapse list */
					/*  so it can be made time-varying */ 
/*    spnt->next = (conn*)NULL;
    if (!synpnt) synpnt = cpnt;  	/* save head if first synap */
/*    if (synend)
      synend->next = cpnt;
    synend = cpnt;
*/
    cpnt->next = (conn*)NULL;			/* zero list pointer */
    if (!connpnt) connpnt = cpnt;	/* increment connection list */
    else if (connend) connend->next = cpnt;
    cpnt->last = connend;		/* pointer to last connection */
    connend = cpnt;

    maklst(&comp1->clst, cpnt);
    maklst(&comp2->clst, cpnt);
    cpnt->comp1 = comp1; 		/* set pointers to comps */
    cpnt->comp2 = comp2; 
    cpnt->compe = NULL; 
    cpnt->compp = NULL; 
    cpnt->num1 = comp1->num;		/* set comp numbers */
    cpnt->num2 = comp2->num;		/*  for parallel implementation */
    cpnt->conduct = s;			/* conductance of connection */
    cpnt->ctype = type;
    cumconn++; 				/* increment total */
    return (cpnt); 
}

/*------------------------------------*/

conn *modconn(elem *cepnt, double s)
               
/* modify a connection betw compartments */

{
    conn *cpnt;
    resistor *epnt;

#ifdef DEBUG
 if (debug & NCMAK && debugz & 8)
   ncfprintf (stderr,"modconn: elem %d\n",cepnt->elnum);
#endif

  if (!(epnt=(resistor *)findelem(cepnt->modif))) {
      ncfprintf (stderr,"modconn: can't find element %d\n",cepnt->modif);
      execerror ("Missing element: "," stopping... ");
      return ((conn*)NULL);
  }

  cpnt = (conn *)epnt->lptr;	/* get pointer to recept from elem */
  if (cpnt==(conn*)NULL) return ((conn*)NULL); 
  /* if (cpnt->ctype != GJ && cpnt->ctype != RESISTOR) {
     ncfprintf (stderr,"modconn: element %d is not a gj %d\n",cepnt->modif);
     return ((conn*)NULL);  
  } */

    cpnt->conduct = s;			/* new conductance of connection */
    epnt->r = 1/s;
    return (cpnt); 
}

/*------------------------------------*/

double getmesg(gj *gjpnt, int mesg)

/* get cA/GMP concentration from remote node (compartment)
   if one is defined, or from neighboring compartments to
   element (normally a gap junction). */

{
   double val;
   node *npnt;

  if (!gjpnt) return 0;
  val = 0;
  if ((npnt=gjpnt->nodec)) {	/* modulation is from defined node */
        val = getnt(npnt->comptr,mesg);
  } /* if nodec */

  else {                /* modulation is from neighboring compartments */
    val =  getnt(gjpnt->comp1,mesg);
    val = (getnt(gjpnt->comp2,mesg) + val) / 2.0;
  }
  return val;
}

/*------------------------------------*/

gj *makgj(gapjunc *gepnt, comp *comp1, comp *comp2, int type)

/* Make a new gj betw compartments and link
    it to both compartments' conn lists. */
/* Connect comp1 to comp2 with conductance s; */

{
    gj *gjpnt;
    node *npnt;
    double  gmax, gnv, n, tau, vg, vj, voff;
    double alpha, beta, specres, area, cyc, atp_decr;

#ifdef DEBUG
 if (debug & NCMAK && debugz & 16)
   ncfprintf (stderr,"makgj: comp %d comp %d\n",comp1->num, comp2->num);
#endif

    if (type != GJ && type != PNX) {
      ncfprintf (stderr,"makgj: incorrect type %d\n",type);
      return ((gj*)NULL);
    }
    if (type==GJ) {
       if ((gjpnt=(gj *)emalloc(sizeof(gj))) == (gj*)NULL) {
          ncfprintf (stderr,"no space left for gj %d\n",cumconn+1);
          return ((gj*)NULL);  
       }
    } else {   // PNX
       if ((gjpnt=(gj *)emalloc(sizeof(pnx))) == (gj*)NULL) {
         ncfprintf (stderr,"no space left for pnx %d\n",cumconn+1);
         return ((gj*)NULL);  
       }
    }
    gjpnt->next = (gj*)NULL;		/* zero list pointer */
    linksyn (gjpnt);			/* add to list of synapses */

    maklst(&comp1->clst, gjpnt);	/* set the list pointers */
    maklst(&comp2->clst, gjpnt);
    gjpnt->comp1 = comp1; 		/* set pointers to comps */
    gjpnt->comp2 = comp2; 
    gjpnt->num1 = comp1->num;		/* set comp numbers */
    gjpnt->num2 = comp2->num;		/*  for parallel implementation */
    gjpnt->ctype = type;
    cumgj++; 				/* increment total */

			/* Set gj conductance. */
			/* The idea is to use what info is given, */
			/* letting "gmax" override "area", */
			/* and "specres" override "drg". */
 
    if ((gmax=gepnt->gmax)==NULLVAL) { 
      if ((specres=gepnt->specres) == NULLVAL) specres = drg;
      if (specres <= 0.0) {
	 ncfprintf (stderr,"makgj: specres is zero, using default 5e6\n"); 
	 specres = 5e6;
      } 
      if ((area=gepnt->area) == NULLVAL) area = 0.0004;  /* um2 */
      if (area <= 0.0) {
	 ncfprintf (stderr,"makgj: area is zero, using default 0.04\n"); 
	 area = 0.0004;
      } 
      gmax = area / specres;          /* conductance = gj area / spec res */
    }
    if (gmax < 0.0) gmax = 100e-12;  /* absolute default */

    if ((vg=gepnt->vgain)==NULLVAL) vg = 1.0;
    gjpnt->vgain = vg; 			/* voltage gain */
    if ((voff=gepnt->voff)==NULLVAL) voff = dgjoff;
    gjpnt->voff = voff; 		/* voltage offset */
    if ((tau=gepnt->taun)==NULLVAL) tau = 1.0;
    gjpnt->taun = tau;			/* extra rate factor */

	/* reverse parameters */

    if ((vg=gepnt->rvgain)!=NULLVAL) gjpnt->rev=1;
    else vg = gjpnt->vgain;
    if ((voff=gepnt->rvoff)!=NULLVAL) gjpnt->rev=1;
    else voff = gjpnt->voff;
    if ((tau=gepnt->rtaun)!=NULLVAL) gjpnt->rev=1;
    else tau = gjpnt->taun;  
    gjpnt->rvgain = vg; 	
    gjpnt->rvoff = voff; 	
    gjpnt->rtaun = tau;			

    gjpnt->rect = gepnt->rect;			/* gj is pure rect */
    if ((gnv=gepnt->gnv)==NULLVAL) gnv = dgjnv;
    if (gnv < 0) gnv = 0;
    if (gnv > 1) gnv = 1;
    gjpnt->maxcond = gmax; 		/* max conductance, never varies */
    gjpnt->gnv = gnv;			/* frac of conductance not modulated */


   cyc = 0;
   gjpnt->modtyp=gepnt->modtyp;
   if (gjpnt->modtyp != 0) {   		/* if cycA/G modulation */
     gjpnt->sign = gepnt->sign;		/* sign of modulation */
     if (gepnt->nodeca != NULLVAL) {
       if (!(npnt=findnode(gepnt->nodeca,gepnt->nodecb,
			   gepnt->nodecc,gepnt->nodecd,"makgj"))){
		execerror ("makgj: can't find mod node for gj","Stopping...");
	}
	gjpnt->nodec=npnt;	/* save node pointer for future ref */
     }
     cyc = getmesg (gjpnt, gepnt->modtyp);
     if (gjpnt->sign) cyc = 1.0 - cyc;
   }
					/* calculate equilib conductance */
    vj = ncabs(gjpnt->comp1->v - gjpnt->comp2->v);
    alpha = gjalpha(vj, gjpnt);
    beta =  gjbeta (vj, gjpnt);
    n = alpha / (alpha + beta);
    gjpnt->n = n;
    if (gjpnt->rect != 0) {
         if ((gjpnt->comp1->v - gjpnt->comp2->v) > 0)
              gjpnt->conduct = 1;
         else
              gjpnt->conduct = 0;
         if (gjpnt->rect < 0) gjpnt->conduct = !gjpnt->conduct;
             gjpnt->conduct *= gjpnt->maxcond;
    }
    else gjpnt->conduct = ( gnv + n * (1 - gnv) ) * gmax * (1.0-cyc);
    return (gjpnt); 
}

/*------------------------------------*/

pnx *makpnx(pannex *pepnt, comp *comp1, comp *comp2, int type)

/* make a pannexin channel that uses gj structure along with
 * additional state variables and parameters */

{
    double atp_decr, atp, atp_gain, adp, amp, h;
    pnx *pnxpnt;

    pnxpnt = (pnx*)makgj( (gapjunc *)pepnt, comp1, comp2, type);

    if ((atp_decr=pepnt->atp_decr)==NULLVAL) atp_decr = dgjatpdec;
    pnxpnt->atp_decr = atp_decr;	/* decrement for atp in cleft per time step */
    pnxpnt->pnxfunc=(double(*)(pnx*,double))pepnt->pnxfunc;	/* function to process atp signal */

    if ((atp=pepnt->atp)!=NULLVAL) pnxpnt->atp   = atp;
    else                           pnxpnt->atp   = 0.2;

    if ((atp_gain=pepnt->atp_gain)!=NULLVAL) pnxpnt->atp_gain   = atp_gain;
    else                           pnxpnt->atp_gain   = 2e14;
    pnxpnt->atpx  = 0;
    pnxpnt->atpo  = 0;
    if ((adp=pepnt->adp)!=NULLVAL) pnxpnt->adp   = adp;
    else                           pnxpnt->adp   = 0.1;
    pnxpnt->adpx  = 0;
    pnxpnt->adpo  = 0;
    if ((amp=pepnt->amp)!=NULLVAL) pnxpnt->amp   = amp;
    else                           pnxpnt->amp   = 0.1;
    pnxpnt->ampx  = 0;
    pnxpnt->ad    = 0;
    if ((h=pepnt->h)!=NULLVAL)     pnxpnt->h   = h;
    else                           pnxpnt->h   = exp(-7.4*LN10);
    pnxpnt->pi    = 0;
    pnxpnt->pio   = 0;
    pnxpnt->atpase = 0;
}

/*------------------------------------*/

gj *modgj(gapjunc *gepnt)
               
/* Modify a gap junction. */
/* Currently only maxcond can be modified. */

{
    gj *gjpnt;
    elem *epnt;
    gapjunc *gopnt;
    node *npnt;
    double gmax, n, vj, tau, alpha, beta;
    double area, specres, gnv, cyc, gjrate, atp_decr;
    int recalc;

#ifdef DEBUG
 if (debug & NCMAK && debugz & 16)
   ncfprintf (stderr,"modconn: elem %d\n",gepnt->elnum);
#endif

  if (!(epnt=findelem(gepnt->modif))) {
      ncfprintf (stderr,"modgj: can't find element %d\n",gepnt->modif);
      execerror ("Missing element: "," stopping... ");
      return ((gj*)NULL);
  }

  gjpnt = (gj *)epnt->lptr;	/* get pointer to recept from elem */
  if (gjpnt==NULL) return ((gj*)NULL); 
  /* if (gjpnt->ctype != GJ) {
     ncfprintf (stderr,"modgj: element %d is not a gj %d\n",gepnt->modif);
     return ((gj*)NULL);  
  } */


    if ((gmax=gepnt->gmax)!=NULLVAL) {
	recalc = 1;
    }
    else {
      gopnt = (gapjunc *)epnt;  /* pointer to orig gap junction spec */

      recalc = 0;		/* first, find if new numbers */
      if ((specres = gepnt->specres) != NULLVAL) recalc = 1;
      if ((area=gepnt->area) != NULLVAL) recalc = 1;
      if (recalc) {		/* if so, find correct vals */
        if ((specres = gepnt->specres) == NULLVAL) {
          if ((specres = gopnt->specres) == NULLVAL) {
           specres = drg;
	  }
	  gopnt->specres = specres;
	}
        if (specres <= 0.0) {		/* check specres for legal value */
	   execerror ("modgj: specific res is zero.","stopping...\n");
	   return ((gj*)NULL);
        } 
        if ((area=gepnt->area) == NULLVAL) {		/* um2 */
          if ((area=gopnt->area) == NULLVAL) {		/* um2 */
	   area = 0.0004;		/* small default value */
	  }
	  gopnt->area = area;
	}
        if (area <= 0.0) {		/* check area for legal value */
	   execerror ("modgj: area is zero.","stopping..."); 
	   return ((gj*)NULL);
        } 
      }  /* if recalc */
      gmax = area / specres; 
    }
    if (gmax < 0.0) gmax = 100e-12;	/* absolute default */
    if (recalc) gjpnt->maxcond = gmax;
    else        gmax = gjpnt->maxcond;

    if (gepnt->vgain !=NULLVAL) gopnt->vgain  = gjpnt->vgain  = gepnt->vgain; 
    if (gepnt->voff  !=NULLVAL) gopnt->voff   = gjpnt->voff   = gepnt->voff; 
    if (gepnt->taun  !=NULLVAL) gopnt->taun   = gjpnt->taun   = gepnt->taun; 
    if (gepnt->rvgain!=NULLVAL) gopnt->rvgain = gjpnt->rvgain = gepnt->rvgain; 
    if (gepnt->rvoff !=NULLVAL) gopnt->rvoff  = gjpnt->rvoff  = gepnt->rvoff; 
    if (gepnt->rtaun !=NULLVAL) gopnt->rtaun  = gjpnt->rtaun  = gepnt->rtaun; 
    // if (gepnt->atp_decr!=NULLVAL) gopnt->atp_decr  = gjpnt->atp_decr  = gepnt->atp_decr; 
    if (gepnt->gnv   !=NULLVAL) gopnt->gnv    = gjpnt->gnv    = gepnt->gnv;
    if (gepnt->modtyp) gopnt->modtyp = gjpnt->modtyp = gepnt->modtyp;   

    if ((tau=gepnt->taun)==NULLVAL) gjrate = 1.0;
    else {
         if (tau < 0) tau = 1;
         if (tau == 0) gjrate = 0;
         else gjrate = dgjtau / tau;
    }
    gjpnt->taun = gjrate;		/* extra rate factor */

    if ((tau=gepnt->rtaun)==NULLVAL) gjrate = 1.0;
    else {
         if (tau < 0) tau = 1;
         if (tau == 0) gjrate = 0;
         else gjrate = dgjtau / tau;
    }
    gjpnt->rtaun = gjrate;		/* extra rate factor */
    // if ((atp_decr=gepnt->atp_decr)==NULLVAL) atp_decr = 1.0;
    // gjpnt->atp_decr = atp_decr;                               //add this to modpannex()

    cyc = 0;
    if (gjpnt->modtyp != 0) {   /* if cycA/G modulation */
      if (gepnt->nodeca != NULLVAL) {
        if (!(npnt=findnode(gepnt->nodeca,gepnt->nodecb,
			    gepnt->nodecc,gepnt->nodecd,"makgj"))){
 		execerror ("modgj: can't find mod node for gj","Stopping...");
 	}
	gjpnt->nodec=npnt;	/* save node pointer for future ref */
      }
      cyc = getmesg (gjpnt, gepnt->modtyp);
      if (gjpnt->sign) cyc = 1.0 - cyc;
    }
					/* calculate equilib conductance */
    vj = ncabs(gjpnt->comp1->v - gjpnt->comp2->v);
    alpha = gjalpha(vj, gjpnt);
    beta =  gjbeta (vj, gjpnt);
    n = alpha / (alpha + beta);
    gjpnt->n = n;
    if (gjpnt->rect != 0) {
         if ((gjpnt->comp1->v - gjpnt->comp2->v) > 0)
              gjpnt->conduct = 1;
         else
              gjpnt->conduct = 0;
         if (gjpnt->rect < 0) gjpnt->conduct = !gjpnt->conduct;
             gjpnt->conduct *= gmax;
    }
    gnv = gjpnt->gnv;			/* conductance not modulated */
    gjpnt->conduct = ( gnv + n * (1 - gnv) ) * gmax * (1.0-cyc);
    return (gjpnt); 
}

/*------------------------------------*/

dbuf *makdelbuf(elem *epnt, comp *comp1, comp *comp2)
                       
/* make an amplifier with gain, offset, and delay connecting comp1 to comp2 and link
    it to both compartments' conn lists. */

{
    dbuf *cpnt;
    int i,idelay,lphp;
    double *dpnt, tau;
    double v, delay;
    vbuf *vpnt;

#ifdef DEBUG
 if (debug & NCMAK && debugz & 16)
   ncfprintf (stderr,"makdelbuf: comp %d comp %d\n",comp1->num, comp2->num);
#endif

    vpnt = (vbuf *)epnt;

    if ((cpnt=(dbuf *)emalloc(sizeof(dbuf))) == (dbuf*)NULL) {
      ncfprintf (stderr,"no space left for conn %d\n",cumconn+1);
      return ((dbuf*)NULL);  
    }
#define MAXDELAY 10000
   if (vpnt->delay != NULLVAL) delay = vpnt->delay;
   else                        delay = 0;
					/* delay in msec, */
   idelay = (int)(delay / stiminc);	/*  time steps of 100 usec (==stiminc) */
   if (idelay < 0) {
      idelay = 0; ncfprintf (stderr,"Error: negative delay in buffer, setting to zero.\n");
   }   
   else if (idelay > 10000) {
      ncfprintf (stderr,"Warning: possible error; excessive delay %d in buffer\n", 
				idelay);
      if (idelay > MAXDELAY) idelay = MAXDELAY;
   }

   if (idelay)  {
     if ((dpnt=(double *)emalloc(idelay * sizeof(double))) == 
					(double*)NULL) {
       ncfprintf (stderr,"no space left for delay buffer %d\n",cumconn+1);
       return ((dbuf*)NULL);  
     }
   }
   else dpnt = (double*)NULL;
					/* add delay buffer to synapse list */
					/*  so it can be made time-varying */ 
   cpnt->next = (dbuf*)NULL;
   linksyn(cpnt);			/* add to list of synapses */

   maklst(&comp1->clst, (conn*)cpnt);
   maklst(&comp2->clst, (conn*)cpnt);
   cpnt->comp1 = comp1; 
   cpnt->comp2 = comp2; 
   cpnt->conduct = 0;			/* conductance of connection */
   cpnt->delay = idelay;		/* length of delay buffer (in stiminc==1e-4sec) */
   if (vpnt->offset  !=NULLVAL) cpnt->offset = vpnt->offset;	/* voltage offset (volts) */
   else                         cpnt->offset = 0;
   if (vpnt->gain  !=NULLVAL)   cpnt->gain = vpnt->gain;	/* voltage gain */
   else                         cpnt->gain = 1.0;
   if (vpnt->tau   !=NULLVAL)   tau = vpnt->tau;		/* low, high pass time const, sec */
   else                         tau = 0;
   if (vpnt->lphp   !=NULLVAL)  lphp = vpnt->lphp;		/* low pass or high pass */
   else                         lphp = LP;

   cpnt->ctype = vpnt->ctype;
   cpnt->delbuf = dpnt;		/* pointer to circular buffer */
   cpnt->delpnt = dpnt;		/* input/output pointer */

   v = (comp1->v - cpnt->offset) * cpnt->gain;
   if (dpnt)
     for (i=0; i<idelay; i++) 
       *dpnt++ = v;			/* set buffer to constant voltage */
   cpnt->v = v;	
   if (tau>0) 
	   cpnt->filt  = maksynfilt (lphp, 1, 0, 0, v, &tau, stiminc);
   cumdbuf++;				/* increment total */
   return (cpnt); 
}

/*------------------------------------*/

ndbuf *makndbuf(elem *epnt, comp *comp1, comp *comp2)
                       
/* make a ntrans buffer connecting comp1 to ntrans comp in comp2 and link
    it to comp1's conn list. */

{
    ndbuf *cpnt;
    int ntrans;
    nbuf *npnt;

#ifdef DEBUG
 if (debug & NCMAK && debugz & 16)
   ncfprintf (stderr,"makndbuf: comp %d comp %d\n",comp1->num, comp2->num);
#endif

   npnt = (nbuf *)epnt;

   if ((cpnt=(ndbuf *)emalloc(sizeof(ndbuf))) == (ndbuf*)NULL) {
      ncfprintf (stderr,"no space left for conn %d\n",cumconn+1);
      return ((ndbuf*)NULL);  
   }
   maklst(&comp1->clst, (conn*)cpnt);
   maklst(&comp2->clst, (conn*)cpnt);
   cpnt->comp1 = comp1;
   cpnt->comp2 = comp2;

   cpnt->ctype = npnt->ctype;
   if (npnt->offset !=NULLVAL)   cpnt->offset = npnt->offset;	/* voltage offset */
   else				 cpnt->offset = 0;
   if (npnt->ntoffset !=NULLVAL) cpnt->ntoffset = npnt->ntoffset; /* nt offset */
   else				 cpnt->ntoffset = 0;
   if (npnt->gain !=NULLVAL)     cpnt->gain = npnt->gain;	/* voltage gain */
   else				 cpnt->gain = 1;
   if (npnt->ntrans !=NULLVAL)   cpnt->ntrans = npnt->ntrans;	/* ntrans */
   else				 cpnt->ntrans = 0;

   cumnbuf++;				/* increment total */
   return (cpnt);
}

/*------------------------------------*/

chan *alloc_chan(short ctype, short stype)

/* Allocate space for channels according to their type */

{
   int chansiz,hh,ns;
   chan *chpnt;
   chantype *chtyp;
   stconc *conc;

  if (stype < 0) stype = 0;
  if (stype > NCHANCONST) stype = NCHANCONST;

  chpnt = (chan *)NULL;

  chtyp = getchantype(ctype,stype);
  ns = chtyp->numstate;
  hh = chtyp->hh;
  switch (ctype) {
    case NA:
 	 if (hh) chansiz = sizeof(hhchan);
         else    chansiz = sizeof(chan);
       break;	/* case NA */
 
    case K:
       switch (stype) {

	case 0:					/* HH */
	case 2:					/* KA */
	 chansiz = sizeof(hhchan);
         break;

	case 1:					/* Markov */
	case 3:					/* Markov */
	case 4:					/* Markov */
	case 5:					/* Markov */
	case 6:					/* Markov */
	case 7:					/* Markov */
	case 8:					/* Markov */
	case 9:					/* Markov */
         chansiz = sizeof(chan);
         break;

	default:
 	  if (hh) chansiz = sizeof(hhchan);
          else    chansiz = sizeof(chan);
         break;

       }  /* switch (stype) */
       break;

    case ClCa:
    case KCa:
       switch (stype) {

	case 6:					/* bKCa, Markov */
	case 5:					/* sKCa, Markov */
	case 4:					/* sKCa, Markov */
	case 3:					/* bKCa, Markov */
	case 2:					/* bKCa */
	case 1:					/* sKCa, Markov */
	case 0:					/* sKCa */
       default:
         chansiz = sizeof(kcachan);
         break;

       }  /* switch (stype) */
       break;

    case CA:
       switch (stype) {

	case 0:
	   chansiz = sizeof(hhchan);	/* simple H-H */
	  break;
        default:			/* channel defined by sequen states */
 	  if (hh) chansiz = sizeof(hhchan);
          else    chansiz = sizeof(chan);
	  break;
        }  /* switch (stype) */
       break;

    case AMPA:
    case KAINATE:
    case CGMP:
    case NMDA:
    case GABA:
    case GLY:
    case SYN2:
         chansiz = sizeof(ntchan);
        break;

    case CHRC:
         chansiz = sizeof(chrchan);
        break;

    }  /* switch (ctype) */

  if ((chpnt=(chan *)emalloc(chansiz)) == (chan*)NULL) {
      ncfprintf (stderr,"no space left for chan %d\n",cumchan+1);
      return ((chan*)NULL);  
  }

  chpnt->numstate = ns;
  if (ns) {
     if ((conc=(stconc *)emalloc(ns*sizeof(stconc))) == (stconc*)NULL) {
      ncfprintf (stderr,"no space left for %d stconc in chan %d\n",ns,cumchan+1);
      return ((chan*)NULL);  
    } 
   chpnt->conc = conc; 
  }
  chpnt->num = cumchan;
  return chpnt;
}

/*------------------------------------*/

void linkchan (chan *chpnt)

/* link a channel into chan list */

{
  chpnt->next = (chan*)NULL;
  if (!chanpnt) {
    chpnt->last = (chan*)NULL;
    chanpnt = chpnt;  	/* save head if first chan */
  }
  if (chanend) {
    chanend->next = chpnt;
    chpnt->last = chanend;
  }
  chanend = (chan*)chpnt;
}

/*------------------------------------*/

int chaniss (int ctype, int stype)

/* return a 1 if channel is sschan */

{
    int typ;
    chantype *chtyp;

  chtyp = getchantype(ctype,stype);
  if (chtyp->numstate==0) typ = 0;
  else                    typ = 1;
  return typ;
}

/*------------------------------------*/

void initsschan(chan *chpnt, comp *comp1, comp *comp2, chantype *chtyp,
		double maxcond, double vrev, 
		double arate, double brate, double crate, 
		double drate, double erate, double frate)

/* initialize a sequential-state channel's states */

{
    int i, nstate;

  chpnt->comp1 = comp1;			/* set loc of voltage sens. */
  chpnt->comp2 = comp2;
  chpnt->compe = NULL;
  chpnt->compp = NULL;
  nstate = chpnt->numstate;  
  if (nstate > 0) {
    for (i=0; i<nstate; i++) {		/* set initial concs */
       chpnt->conc[i].cval = 0.0;
       chpnt->conc[i].cest = 0.0;
       chpnt->conc[i].nchan = 0;
    }
    chpnt->conc[0].cval = 1.0;
    chpnt->conc[0].cest = 1.0;
  }
  chpnt->maxcond  = maxcond;
  chpnt->vrev = vrev;
  chpnt->arate = arate;
  chpnt->brate = brate;
  chpnt->crate = crate;
  chpnt->drate = drate;
  chpnt->erate = erate;
  chpnt->frate = frate;
}

/*------------------------------------*/

void equilchans(void)

/* Equilibrate all channels so they start in correct states. */
/*  Normally this routine should be called only after condense().  */

{
   conn *cpnt;
  for (cpnt=(conn*)chanpnt; cpnt; cpnt=cpnt->next) {
    dochani((chan*)cpnt,1e-8);
  }
}

/*------------------------------------*/

void setnchan(chan *chpnt)

/* Set the nchan (number of channels) in a chan to be integer. */
/* Use after all chans in a compartment have been condensed.   */

/* The nchan param of a channel is floating point (double) to allow */
/*  several channels to be added together during the condensation */
/*  process.  This procedure enforces an integral number and */
/*  modifies "maxcond" so its value is consistent with "nchan". */

{
     double unitary, nchan;

   if (!chpnt) return;
   switch (chpnt->ctype) {

         case CA:
         case  K:
        case ClCa:
        case KCa:
         case NA:
       case GABA:
       case GLY:
       case CHRC:
       case SYN2:
       case NMDA:
       case CGMP: break;

         case GJ:
        case CAP:
        case ROD:
       case CONE:
       case CHR:
       case VTRANSDUCER:
       case ITRANSDUCER:
       case BATT:
    case SYNAPSE:
   case AXIALRES:
   case RESISTOR:
       case LOAD:
         default: return;    /* ignore non-channel connections */
		  break;
   }
   if (chpnt->nchan > 0) {
      unitary = chpnt->maxcond / chpnt->nchan;
      nchan = int(chpnt->nchan + 0.5);
      chpnt->nchan = nchan;  
      chpnt->maxcond = unitary * nchan;  
   }
}

/*------------------------------------*/

void maknois(chan *chpnt)

/* Calculate the population of an sschan's states when it is first defined.  */

{
        int i,nstate, nchano, totchan;
        chanstate *stpnt;
        stconc *conc;

   if (chpnt->nchan > 0) {		/* set up noisy states */
      nstate = chpnt->numstate;
      conc = chpnt->conc;
      stpnt = chpnt->chtyp->state;
      nchano = totchan = 0;
      for (i=0; i<nstate; i++,stpnt++,conc++) {            /* Total */
         conc->nchan = int(conc->cest * chpnt->nchan + 0.5);
         totchan += conc->nchan;
         nchano += int(conc->nchan*stpnt->cond+0.5);
      }
      chpnt->conc[0].nchan += int(chpnt->nchan-totchan); /* Make correct tot */
      if (nchano > chpnt->nchan) nchano = int(chpnt->nchan);
      chpnt->cno = nchano;
      /* conduct = ((double)nchano)/chpnt->nchan; */
   }
}

/*------------------------------------*/

void sethhinf (chan *chpnt, comp *comp1)

/* set HH m and h to equilibrium values */

{
   chanparm *chp;

  chp = &chpnt->chtyp->parm[0];
  ((hhchan *)chpnt)->m  = calcchaninf(comp1->v-vext(chpnt)-ccavoff(chpnt)-chpnt->voffsm,chp,
					chpnt->arate,chpnt->brate);
  if (chpnt->chtyp->numparm > 1) {
    chp++;
    ((hhchan *)chpnt)->h = calcchaninf(comp1->v-vext(chpnt)-ccavoff(chpnt)-chpnt->voffsh,chp,
					chpnt->crate,chpnt->drate);
  }
}

/*------------------------------------*/

chan *makchan(attrib *apnt, comp *comp1, comp *comp2, int nocondens)
                 
/* make a new channel at a compartment and link
   it to the compartment's conn list. */

{
    int stype,stop,ctype;
    chan *chpnt;
    chan *addchan(comp *pnt, chan *chpnt);
    chantype *chtyp;
    attrib *npnt;
    mattrib *chapnt;
    chattrib *chapntk;
    double chn,chunitary;
    double voffm,voffh,vrev,maxcond,rextern;
    double arate,brate,crate,drate,erate,frate;
    double dtaua,dtaub,dtauc,dtaud,dtaue,dtauf;
    double taua,taub,tauc,taud,taue,tauf;
    double doffsm, doffsh, dvrev, caperm, cakd, cahc;
    double alph,bet,k1,k2,d1,d2,ca,mg;

#ifdef DEBUG
 if (debug & NCMAK && debugz & 16)
   ncfprintf (stderr,"makchan: type %d %d\n",apnt->ctype,apnt->stype);
#endif
    			/* first, figure out how much space to allocate */
  ctype = apnt->ctype;
  if (apnt->stype == NULLVAL) stype = 0;	/* channel sub-type */
  else			      stype = apnt->stype;
  chtyp = getchantype(ctype,stype);

			/* next, figure out a value for maxcond */
			/* this may affect unitary cond or N */

  chapnt = (chattrib *)apnt;
  chunitary = chapnt->unitary;
  chn = chapnt->n;

  if ((maxcond=chapnt->maxcond)==NULLVAL) {  /* if maxcond not defined */

   if (comp2 && ((synap*)comp2)->maxcond != NULLVAL)
       maxcond = ((synap*)comp2)->maxcond;   /* define maxcond from synapse */
   else {

		/* first check if maxcond can be defined from chan N */

     if (chn!=NULLVAL) {
	   double n, unitary;

	n = chn;
        if (chunitary==NULLVAL) {   
           if ((unitary=chtyp->unitary)<=0) unitary = dscu;
	   else unitary *= qcond(chtyp); 
	}
	else unitary = chunitary;
        if (unitary <= 0) n = 0;
        maxcond = n * unitary;	/* maxcond computed from N before chnoise */
     }
    else {
		/* next check if maxcond can be defined from chnoise N */

    for (stop=0,npnt=apnt->attpnt; npnt && !stop; npnt=npnt->attpnt) { 
        nattrib *napnt;

      napnt = (nattrib *)npnt;
      switch (npnt->ctype) {
	  double n, unitary;

        case CCHNOISE:
        if ((n=napnt->n)==NULLVAL) n = chn;	/* use prev.  defined N */
        if (n!=NULLVAL) {
           if ((unitary=napnt->unitary)==NULLVAL) {   
	      if (chunitary!=NULLVAL) unitary = chunitary; 
	      else {		/* use default value, chan type's first */
	        if ((unitary=chtyp->unitary)<=0) unitary = dscu;
	        else unitary *= qcond(chtyp); 
	      }
              if (unitary <= 0) n = 0;
	   }
           maxcond = n * unitary;
        }
	else chn = n;			/* save NULLVAL as evidence for below */
        break;

	case VESNOISE: break;
	case NUNIT: break;

        default: stop=1; break;		/* stop on other channels */
      }
     }   /* for (npnt;;) */
   }  /* else try to define maxcond from chnoise */

   if (chn==NULLVAL) {        /* if maxcond still not defined, use default */
           conlst *nlpnt;
           node *nodpnt;

      switch (ctype) {
      case AMPA:
      case KAINATE:
      case CGMP:
      case GABA:
      case GLY:
      case CHRC:
      case SYN2:
      case NMDA: maxcond = dmaxsyn; break;
      case CA:   maxcond = dmaxca; break;
      case NA:   maxcond = dmaxna; break;
      case ClCa: maxcond = dmaxclca; break;
      case KCa:  maxcond = dmaxk; break;
      case K:    maxcond = dmaxk; break;
      } 
    if ((nlpnt=comp1->nodlst)) nodpnt=(node*)nlpnt->conpnt;
    
 ncfprintf (stderr,"# makchan: maxcond, N for 'chan %s type %d' at node %s not defined,\n",
		findsym(ctype),stype,prnode(nodpnt->nodenm1,nodpnt->nodenm2,
					    nodpnt->nodenm3,nodpnt->nodenm4));
   ncfprintf (stderr,"#   using default %g\n",maxcond);
    }
   }  /* ! (comp2 && ((synap*)comp2)->maxcond != NULLVAL) */

  }   /* if (maxcond==NULLVAL) */

 if (maxcond <= 0 && nozerochan) {
     return ((chan *)NULL);	/* if no conductance, stop */
 }
 else { 
   chpnt = alloc_chan(ctype,stype);
   linkchan(chpnt);
   cumchan++;

   chpnt->ctype = ctype;
   chpnt->stype = stype;
   chpnt->nocondens = nocondens;	/* save for addchan() below */

   dtaue = 1.0;
   dtauf = 0.2;		/* set default noise flicker rate */
   doffsm = 0.0; doffsh = 0.0;
   dvrev = chtyp->vrev;

   switch (ctype) {		/* find default params */

      case GLY: 
      case GABA: dtaua = dtaub = dtauc = dtaud = dsyntau;
		 if (dvrev==NULLVAL) dvrev = vcl;
		 dtauf=dsyntauf;
		 break;
      case SYN2: 
      case AMPA:
      case KAINATE:
      case CGMP:
      case NMDA: dtaua = dtaub = dtauc = dtaud = dsyntau;
		 if (dvrev==NULLVAL) dvrev = .667 * vna + .333 * vk;
		 dtauf=dsyntauf;
		 break;
      case CA:   dtaua=dtaub=dcataum; dtauc=dtaud=dcatauh; 
		 doffsm = dcaoffs; doffsh = dcaoffs;
	 	 if (dvrev==NULLVAL) dvrev = vna;	/* set in makca() */
		 dtauf=dcatauf;
		break;
      case NA:   dtaua=dtaub=dnataum; dtauc=dtaud=dnatauh; 
		 doffsm = dnaoffsm; doffsh = dnaoffsh;
		 if (dvrev==NULLVAL) dvrev = vna;
		 dtauf=dnatauf;
		break;

      case K:    doffsm = dkoffsn; doffsh = dkoffsh;
		 if (dvrev==NULLVAL) dvrev = vk;
		 dtaua=dtaub=dktaum; dtauc=dtaud=dktauh; 
		 dtauf=dktauf;
		break;

      case ClCa:  
      case KCa:  doffsm = dkoffsn; doffsh = dkoffsh;
		 dtaua=dtaub=dkcataum; dtauc=dtaud=dkcatauh; 
		 if (dvrev==NULLVAL) dvrev = vk;
		 dtauf=dktauf;
		break;

      case CHRC:  doffsm = 0; doffsh = 0;
		 if (dvrev==NULLVAL) dvrev = 0;
		 dtaua=dtaub=1; dtauc=dtaud=1.0; 
		break;
  } 

  if ((taua=chapnt->taua) == NULLVAL) taua = dtaua;
  else if (taua <=0)                  taua = 1e-6;
  arate = 1.0 / taua;
  if ((taub=chapnt->taub) == NULLVAL) taub = dtaub;
  else if (taub <=0)                  taub = 1e-6;
  brate = 1.0 / taub;
  if ((tauc=chapnt->tauc) == NULLVAL) tauc = dtauc;
  else if (tauc <=0)                  tauc = 1e-6;
  crate = 1.0 / tauc;
  if ((taud=chapnt->taud) == NULLVAL) taud = dtaud;
  else if (taud <=0)                  taud = 1e-6;
  drate = 1.0 / taud;
  if ((taue=chapnt->taue) == NULLVAL) taue = dtaue;
  else if (taue <=0)                  taue = 1e-6;
  erate = 1.0 / taue;
  if ((tauf=chapnt->tauf) == NULLVAL) tauf = dtauf;
  else if (tauf <=0)                  tauf = 1e-6;
  frate = 1.0 / tauf;


  if ((rextern=chapnt->rextern) == NULLVAL) rextern = 0;

  if ((voffm=chapnt->voffsm) == NULLVAL) voffm = doffsm;
  if ((voffh=chapnt->voffsh) == NULLVAL) voffh = doffsh;

  if ((cakd=chapnt->cakd)== NULLVAL) cakd = dcakd;
  chpnt->cakd = cakd;
  if ((cahc=chapnt->cahc)== NULLVAL) cahc = dcahc;
  chpnt->cahc = cahc;
  
  if ((vrev=chapnt->vrev) == NULLVAL) {
        //double tp;

   // tp = chtyp->ionp[PNA] + chtyp->ionp[PK] + chtyp->ionp[PCL] + caperm;
   // if (tp > 0) vrev = ghkv(chtyp,dcao,dcai);
   // else        vrev = dvrev;
    vrev = dvrev;
    chpnt->setvrev = 0;
  }
  else chpnt->setvrev = 1;		/* vrev set by user */ 

  chpnt->naratio = dnai/dnao;
  chpnt->kratio  = dki/dko;
  chpnt->ions = NULL;
  if ((caperm=chapnt->caperm)!= NULLVAL) {
   	  iontab *ions;
      if (caperm < 0) caperm = 0;
      if (caperm > 1) caperm = 1;
      chpnt->ions = init_ions(chpnt->ions);
      chpnt->ions->ionp[PCA] = caperm;
  }
  chpnt->gfractab = NULL;
  chpnt->gvrevtab = NULL;

  chpnt->rextern  = rextern;

  chpnt->chtyp  = chtyp;
  switch (ctype) {

    case NMDA:
       if ((mg=((cattrib*)chapnt)->mg) == NULLVAL) mg = dnmdamg; 
       ((ntchan*)chpnt)->mg = mg;

                /* continue on to below */
    case AMPA:
    case KAINATE:
    case CGMP:
    case GABA:
    case GLY:
    case SYN2:
       chpnt->voffsm = voffm; 
       chpnt->voffsh = voffh; 
       switch (stype) {
         default:
	   case 1:
             initsschan(chpnt,comp1,comp2,chpnt->chtyp,maxcond,vrev,
			arate,brate,crate,drate,erate,frate);
	   break;
         }
        break;

    case NA:
       chpnt->voffsm = voffm; 
       chpnt->voffsh = voffh; 
       initsschan(chpnt,comp1,comp2,chpnt->chtyp,maxcond,vrev,
		arate,brate,crate,drate,erate,frate);
       if (chpnt->chtyp->hh) sethhinf(chpnt,comp1);
       break;

    case K:
       chpnt->voffsm = voffm; 
       chpnt->voffsh = voffh; 
       initsschan(chpnt,comp1,comp2,chpnt->chtyp,maxcond,vrev,
		arate,brate,crate,drate,erate,frate);
       if (chpnt->chtyp->hh) sethhinf(chpnt,comp1);
       break;

    case ClCa:
               
	if (comp1->capnt) {				/* found calcium comp*/
	     ((kcachan *)chpnt)->initfl=1;		/* initialized OK */
	}
	else {		/* calcium not set up yet, must do later */
	     ((kcachan *)chpnt)->initfl=0;	/* not properly initialized */
	}

        initsschan(chpnt,comp1,comp2,chpnt->chtyp,maxcond,vrev,
		arate,brate,crate,drate,erate,frate);
	break;

    case KCa:
       chpnt->voffsm = voffm; 
       chpnt->voffsh = voffh; 

       chapntk = (chattrib*)chapnt;
       switch (stype) {			/* select default rate constants */
            case 0:
            case 1:
            case 4:
            case 5:
             if ((d1=chapntk->d1)==NULLVAL) d1 = dsd1; 
             if ((d2=chapntk->d2)==NULLVAL) d2 = dsd2; 
             if ((k1=chapntk->k1)==NULLVAL) k1 = dsk1; 
             if ((k2=chapntk->k2)==NULLVAL) k2 = dsk2; 
	     break;
            case 2:
            case 3:
            case 6:
             if ((d1=chapntk->d1)==NULLVAL) d1 = dbd1; 
             if ((d2=chapntk->d2)==NULLVAL) d2 = dbd2; 
             if ((k1=chapntk->k1)==NULLVAL) k1 = dbk1; 
             if ((k2=chapntk->k2)==NULLVAL) k2 = dbk2; 
	     break;
       }
       switch (stype) {			/* set default rate constants */
            case 0:
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
           default: 

	  ((kcachan *)chpnt)->d1 = d1;
	  ((kcachan *)chpnt)->d2 = d2;
	  ((kcachan *)chpnt)->k1 = k1;
 	  ((kcachan *)chpnt)->k2 = k2;
            break;

       }
       switch (stype) {		/* special things that some KCa chans need */
	case 0:						/* Kca SK, HH */
	case 2:						/* Kca BK, HH */
	  if (comp1->capnt) {				/* found calcium comp*/
	     ca = comp1->capnt->cais[0];
	     ((kcachan *)chpnt)->initfl=1;		/* initialized OK */
	  }
	  else {		/* calcium not set up yet, must do later */
	     ca = dcai;
	     ((kcachan *)chpnt)->initfl=0;	/* not properly initialized */
	  }
	  alph = akcacalc((comp1->v-vext(chpnt)-ccavoff(chpnt)-chpnt->voffsm),ca,arate,d1,k1);
	  bet  = bkcacalc((comp1->v-vext(chpnt)-ccavoff(chpnt)-chpnt->voffsm),ca,arate,d2,k2);
	  ((kcachan *)chpnt)->m = alph / (alph + bet);
          ((kcachan *)chpnt)->h = 1;
	  break;
					/* same as 2: above */
					/*  but sequential state */

	case 1:				/* SKca, Markov */
	case 3:				/* BKCa, Markov */
	case 4:				/* SKCa, Markov, Hirshberg et al 1998 */
	case 5:				/* SKCa, Sah and Clements, 1999 */
	case 6:				/* BKCa, Markov, Horrigan  et al 1999 */
	default: 
               
	  if (comp1->capnt) {				/* found calcium comp*/
	     ((kcachan *)chpnt)->initfl=1;		/* initialized OK */
	  }
	  else {		/* calcium not set up yet, must do later */
	     ((kcachan *)chpnt)->initfl=0;	/* not properly initialized */
	  }
	  break;

       } /* switch (stype) */
       initsschan(chpnt,comp1,comp2,chpnt->chtyp,maxcond,vrev,
		arate,brate,crate,drate,erate,frate);
       if (chpnt->chtyp->hh) sethhinf(chpnt,comp1);
       break;

    case CA:
       chpnt->voffsm = voffm; 
       chpnt->voffsh = voffh; 
       initsschan(chpnt,comp1,comp2,chpnt->chtyp,maxcond,vrev,
		arate,brate,crate,drate,erate,frate);
       if (chpnt->chtyp->hh) sethhinf(chpnt,comp1);
      break;

    case CHRC:
       chpnt->voffsm = voffm; 
       chpnt->voffsh = voffh; 
       initsschan(chpnt,comp1,comp2,chpnt->chtyp,maxcond,vrev,
		arate,brate,crate,drate,erate,frate);
       break;

  }
  chpnt->nchan = 0;			/* set chan noise params */
  for (stop=0,npnt=apnt->attpnt; npnt && !stop; npnt=npnt->attpnt) { 
      nattrib *napnt;

      napnt = (nattrib *)npnt;
      switch (npnt->ctype) {
	   double nchan, unitary;
	   int nrseed;

 	case CCHNOISE:			/* If "chnoise" has been set */

	if ((unitary=napnt->unitary)==NULLVAL) {	/* find unitary */
	   if (chunitary!=NULLVAL) unitary = chunitary; 
	   else {		/* compute unit from N, or default value */
	     nchan = napnt->n;		
	     if (nchan==NULLVAL) nchan=chn;  /* Check for prev defined N */
	     if ((nchan==0)||(nchan==NULLVAL)) {/* N undef, use default unit */
	       if ((unitary=chtyp->unitary)<=0) unitary = dscu;
	       else unitary *= qcond(chtyp);
	     }
	     else unitary = chpnt->maxcond / nchan;
	   }
        }
	if (unitary < 0) unitary = dscu;
	if (unitary < 0) unitary = 0;
	
	nchan = napnt->n;		/*  both N, unit could be diff value */
	if (nchan==NULLVAL) nchan=chn;  /* Check for previously defined N */
	if (nchan==NULLVAL) { 	/* if N still not defined, use default unitary*/
	   if (unitary <= 0)   nchan = 0;
	   else nchan = chpnt->maxcond / unitary;	
	}
	chpnt->nchan = nchan;			/* floating point (double) number */

	if ((tauf=napnt->tauf)!=NULLVAL) {
		if (tauf <= 0) tauf = 1e-6;
		chpnt->frate = 1.0 / tauf;
	} 
        if (chaniss(ctype,stype)) maknois(chpnt);   /* set state populations */

	nrseed = napnt->rseed;
	if (nrseed >= 0) {
	   if (nrseed == 0) 
	    nrseed = noduniq(25471, chpnt->num, 15149, chpnt->num,srseed^rseed^3131);
          chpnt->cstate = makrand(RNDSIZ,"chan",cumchan);
	  if (chpnt->cstate) initstate (nrseed,chpnt->cstate,RNDSIZ);
	}

	break;

	case VESNOISE: break;
	case NUNIT: break;

	default: stop=1; break;		/* stop on other channels */
      }
  }
					/* add chan to comp's connections */
  chpnt = addchan(comp1, (chan *)chpnt);
					/*  but possibly erase new channel   */
					/*  and return pointer to old one.   */
#ifdef DEBUG
 if (debug & NCMAK && debugz & 16)
   ncfprintf (stderr,"makchan end\n");
#endif
    return (chpnt); 

  }  /* else (maxcond > 0) */


}

/*------------------------------------*/

chan *modchan(elem *cepnt, attrib *apnt, double lcap)
                 
/* Modify a channel. */

{
    int stype,stop;
    elem *epnt;
    comp *comp1;
    chan *chpnt;
    cacomp *capnt;
    chantype *chtyp;
    attrib *npnt;
    chattrib *chapnt;
    cattrib *cacompnt;
    chattrib *ochapnt;
    double maxcond,area,m,h,caperm;
    double arate,brate,crate,drate,erate,frate;
    double taua,taub,tauc,taud,taue,tauf;
    double alph,bet,k1,k2,d1,d2,ca;
    double chn, chunitary;
  
#ifdef DEBUG
 if (debug & NCMAK && debugz & 16)
   ncfprintf (stderr,"modchan: type %d %d\n",apnt->ctype,apnt->stype);
#endif

  if (!(epnt=findelem(cepnt->modif))) {
      ncfprintf (stderr,"modchan: can't find element %d\n",cepnt->modif);
      execerror ("Missing element: "," stopping... ");
      return ((chan *)NULL);
  }

  chpnt = (chan *)epnt->lptr;	/* get pointer to chan from elem */
  if (chpnt==((chan*)NULL)) return ((chan *)NULL); 
  comp1 = chpnt->comp1;

  chtyp = chpnt->chtyp;
  stype = chpnt->stype;

  if (lcap<=0) lcap = 1e-6;
  area = comp1->cap / lcap;

  chapnt = (chattrib *)apnt;
  ochapnt = (chattrib*)epnt->attpnt;
  if (ochapnt==NULL) ncfprintf (stderr,"modchan: can't find attrib\n");
  chunitary = chapnt->unitary;
  if (chunitary != NULLVAL) ochapnt->unitary = chunitary;
  chn = chapnt->n;
  if (chapnt->density != NULLVAL)  ochapnt->maxcond = chapnt->maxcond=chapnt->density*area;
  if (chapnt->ndensity != NULLVAL) ochapnt->n = chapnt->n=int(chapnt->ndensity*area*UM2CM2);
  maxcond = chpnt->maxcond;		/* get original maxcond */
  if (chapnt->maxcond != NULLVAL) ochapnt->maxcond = maxcond = chapnt->maxcond;
  else {   /* must calc new maxcond or use old val */

     if (chn!=NULLVAL) {		/* N set in chan attrib */
           double n, unitary;

        n = chn;
        if (chunitary==NULLVAL) {   
           if ((unitary=chtyp->unitary)<=0) unitary = dscu;
           else unitary *= qcond(chtyp); 
        }
        else unitary = chunitary;
        if (unitary <= 0) n = 0;
        maxcond = n * unitary;  /* maxcond computed from N before chnoise */
     }
   else {
     for (stop=0,npnt=apnt->attpnt; npnt && !stop; npnt=npnt->attpnt) {
        nattrib *napnt;
        double n,unitary;

      napnt = (nattrib *)npnt;
      switch (npnt->ctype) {

        case CCHNOISE:	
        if ((n=napnt->n)!=NULLVAL) {		/* N set in noise attrib */
           if ((unitary=napnt->unitary)==NULLVAL) {  /* unitary not defined */
	        if (n <= 0) unitary = 0;
		else unitary = chpnt->maxcond / n + 0.5;
           }
	   maxcond = n * unitary;
        }
        else if ((unitary=napnt->unitary)!=NULLVAL) {/* unitary def but not N */
		if (unitary <= 0) n = 0; 
		else n = chpnt->maxcond / unitary + 0.5;
	        maxcond = n * unitary;
	}
	else maxcond = chpnt->maxcond;		/* else use original maxcond */
        break;

        case NUNIT:
        case VESNOISE:	break;

        default: stop=1; break;
      }
    }
   }
  }
  caperm=chapnt->caperm;
  if (chapnt->caperm != NULLVAL) { 
    caperm = chapnt->caperm;
    if (caperm < 0) caperm = 0;
    if (caperm > 1) caperm = 1;
    chpnt->ions = init_ions(chpnt->ions);
    chpnt->ions->ionp[PCA] = caperm;
    ochapnt->caperm = caperm;
  }

  if (chapnt->vrev != NULLVAL)   ochapnt->vrev   = chpnt->vrev = chapnt->vrev;
  if (chapnt->voffsm != NULLVAL) ochapnt->voffsm = chpnt->voffsm = chapnt->voffsm;
  if (chapnt->voffsh != NULLVAL) ochapnt->voffsh = chpnt->voffsh = chapnt->voffsh; 

  if (chapnt->taua != NULLVAL) ochapnt->taua = taua = chapnt->taua;
  else taua = 1.0/chpnt->arate;
  if (taua <=0) taua = 1e-6;
  arate = 1.0/taua;

  if (chapnt->taub != NULLVAL) ochapnt->taub = taub = chapnt->taub;
  else taub = 1.0/chpnt->brate;
  if (taub <=0) taub = 1e-6;
  brate = 1.0/taub;

  if (chapnt->tauc != NULLVAL) ochapnt->tauc = tauc = chapnt->tauc;
  else tauc = 1.0/chpnt->crate;
  if (tauc <=0) tauc = 1e-6;
  crate = 1.0/tauc;

  if (chapnt->taud != NULLVAL) ochapnt->taud = taud = chapnt->taud;
  else taud = 1.0/chpnt->drate;
  if (taud <=0) taud = 1e-6;
  drate = 1.0/taud;

  if (chapnt->taue != NULLVAL) ochapnt->taue = taue = chapnt->taue;
  else taue = 1.0/chpnt->erate;
  if (taue <=0) taue = 1e-6;
  erate = 1.0/taue;

  if (chapnt->tauf != NULLVAL) ochapnt->tauf = tauf = chapnt->tauf;
  else tauf = 1.0/chpnt->frate;
  if (tauf <=0) tauf = 1e-6;
  frate = 1.0/tauf;

  chpnt->maxcond  = maxcond;
  chpnt->arate = arate;
  chpnt->brate = brate;
  chpnt->crate = crate;
  chpnt->drate = drate;
  chpnt->erate = erate;
  chpnt->frate = frate;

  switch (chpnt->ctype) { 	/* can't change type of chan */

    case NMDA:
       if (((cattrib*)chapnt)->mg != NULLVAL) 
		((cattrib*)ochapnt)->mg = ((ntchan*)chpnt)->mg = ((cattrib*)chapnt)->mg;

                /* continue on to below */
    case AMPA:
    case KAINATE:
    case CGMP:
    case GABA:
    case GLY:
    case SYN2:
       switch (stype) {
         default:
	   case 1:
	  dochani(chpnt,1e-8);
	   break;
         }
        break;

    case NA:
       if (chpnt->chtyp->hh) sethhinf(chpnt,comp1);
       switch (stype) {
        case 0:			/* set up equilibrium rate constants */
	   m = ((hhchan *)chpnt)->m;
	   h = ((hhchan *)chpnt)->h;
       	   ((hhchan *)chpnt)->conduct = m*m*m*h;
	  break;

         default: break;

       }
       dochani(chpnt,1e-8);
       break;

    case K:
       if (chpnt->chtyp->hh) sethhinf(chpnt,comp1);
       switch (stype) {

        default: break;

        case 0:				/* set up equilibrium rate constants */
	   m = ((hhchan *)chpnt)->m;
       	   ((hhchan *)chpnt)->conduct = (m*m)*(m*m);
	  break;

        case 2:				/* type A deactivating potassium */
	   m = ((hhchan *)chpnt)->m;
	   h = ((hhchan *)chpnt)->h;
       	   ((hhchan *)chpnt)->conduct = m*m*m*h;
	  break;

       } /* switch (stype) */
       dochani(chpnt,1e-8);
       break;

    case ClCa:
       dochani(chpnt,1e-8);
	break;

    case KCa:
       if (chpnt->chtyp->hh) sethhinf(chpnt,comp1);
         switch (stype) {
           default:
            case 0:
            case 1:
            case 4:
            case 5:
             if ((d1=chapnt->d1)==NULLVAL) d1 = dsd1; 
             if ((d2=chapnt->d2)==NULLVAL) d2 = dsd2; 
             if ((k1=chapnt->k1)==NULLVAL) k1 = dsk1; 
             if ((k2=chapnt->k2)==NULLVAL) k2 = dsk2; 
             break;
            case 2:
            case 3:
            case 6:
             if ((d1=chapnt->d1)==NULLVAL) d1 = dbd1; 
             if ((d2=chapnt->d2)==NULLVAL) d2 = dbd2; 
             if ((k1=chapnt->k1)==NULLVAL) k1 = dbk1; 
             if ((k2=chapnt->k2)==NULLVAL) k2 = dbk2; 
             break;
         }
         switch (stype) {
           default:
            case 0:
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
	  ((kcachan *)chpnt)->d1 = d1;
	  ((kcachan *)chpnt)->d2 = d2;
	  ((kcachan *)chpnt)->k1 = k1;
 	  ((kcachan *)chpnt)->k2 = k2;
            break;

        }
        switch (stype) {

        default: break;
	case 0:				/* Ca sensitive K channel */
	case 2:				
	  arate = 1/taua;
	  arate = 1/taua;


	  if (comp1->capnt) {				/* found calcium comp*/
	     ca = comp1->capnt->cais[0];
	  }
	  else {		/* calcium not set up yet, must do later */
	     ca = dcai;
	  }
	  alph = akcacalc((comp1->v-vext(chpnt)-ccavoff(chpnt)-chpnt->voffsm),ca,arate,d1,k1);
	  bet  = bkcacalc((comp1->v-vext(chpnt)-ccavoff(chpnt)-chpnt->voffsm),ca,arate,d2,k2);
	  m = alph / (alph + bet);
	  ((kcachan *)chpnt)->m = m;
          ((kcachan *)chpnt)->h = 1;
 	  ((kcachan *)chpnt)->conduct = m;
	  break;

	case 1:				/* Ca sensitive K channel */
	case 3:				
	case 4:
	case 5:
	case 6: break;

       } /* switch (stype) */
       dochani(chpnt,1e-8);
       break;

    case CA:
       if (chpnt->chtyp->hh) sethhinf(chpnt,comp1);
       switch (stype) {
        case 0:			/* set up equilibrium rate constants */
	   m = ((hhchan *)chpnt)->m;
	   h = ((hhchan *)chpnt)->h;
       	   ((hhchan *)chpnt)->conduct = m*m*m;
	  break;

        case 2:			/* set up equilibrium rate constants */
	   m = ((hhchan *)chpnt)->m;
	   h = ((hhchan *)chpnt)->h;
       	   ((hhchan *)chpnt)->conduct = m*m*h;
	  break;

	default: break;

	} /* switch (stype) */
       dochani(chpnt,1e-8);
       break;
  }

  /* If need access to old elem noise attrib, 
	need to add code to access it here. */

  for (npnt=apnt->attpnt; npnt; npnt=npnt->attpnt) { /* set noise params */
      nattrib *napnt;
      double unitary,nchan;

     napnt = (nattrib *)npnt;
     switch (npnt->ctype) {

 	case CCHNOISE:

	nchan = chpnt->nchan;
	if (napnt->n != NULLVAL) nchan = (int)(napnt->n);	

	else {		/* try to calc nchan from unitary cond */

	   if (napnt->unitary != NULLVAL) unitary = napnt->unitary;
	   else {
	     if (nchan == 0) unitary = 0;
	     else unitary = chpnt->maxcond / nchan;
	   }
	   if (unitary <= 0) nchan = 0;
	   else nchan = chpnt->maxcond/unitary + 0.5;	
	}
        chpnt->nchan = int(nchan);

	if (napnt->rseed > 0) {
	   if (chpnt->cstate) initstate (napnt->rseed,chpnt->cstate,RNDSIZ);
        }
	break;

 	case NUNIT:
 	case VESNOISE: break;

 	case CACOMP:		/* set the internal Ca during simulation */
                 cacompnt = (cattrib*)npnt;
		 if ((capnt=chpnt->comp1->capnt)) {
			     int i,cai,cao,cashell,cashell_new,caoshell;
	            if ((cai=cacompnt->cai) != NULLVAL) {
		       capnt->cai = cai;
	               if ((cashell_new=cacompnt->cashell) != NULLVAL) {
                          cashell = capnt->cashell;
			  	/* only allow same or fewer n shell */
			  if (cashell_new > cashell) cashell_new = cashell; 
			  capnt->cashell = cashell = cashell_new;
			} /* cashell */
                        cashell = capnt->cashell;
			for (i=0; i<cashell; i++) {
			   capnt->cais[i] = cai;
			}
		   }  /* cai */
	           if ((cao=cacompnt->cao) != NULLVAL) {
		       capnt->cao = cao;
		       if (capnt->caos!=NULL) {
                          caoshell = capnt->caoshell;
			  for (i=0; i<caoshell; i++) {
			     capnt->caos[i] = cao;
			  }
		       }
		   }  /* cao */
		 }
		 break;
     }
  }
#ifdef DEBUG
 if (debug & NCMAK && debugz & 16)
   ncfprintf (stderr,"modchan end\n");
#endif
  return ((chan*)chpnt); 
}

/*------------------------------------*/

chan *addchan(comp *pnt, chan *chpnt)
               
/* Add a channel to a compartment's connections. */
/* Channels will be condensed together later after */
/*  all compartments have been condensed. */

{
  if (!pnt || !chpnt) return ((chan*)NULL);

   maklst(&pnt->clst, (conn*)chpnt);	/* add to comp's connections */
   chpnt->comp1 = pnt; 
   return chpnt;		/* return pointer to channel */
}

/*------------------------------------*/

chan *addchan_extern_ca(comp *pnt, chan *chpnt)

/* Connect a calcium channel or synapse to an external compartment */

/* 1) Add the channel to the external compartment's connections (because it will add current). */
/* 2) Set channel's compe pointer to point to the external compartment. */
/* 3) Set channel's comp1 external voltage pointer to point to the external comp. */

{

  if (!pnt || !chpnt) return ((chan*)NULL);

   maklst(&pnt->clst, (conn*)chpnt);	/* add to comp's connections */
   chpnt->compe = pnt; 
   chpnt->comp1->pvext = pnt;
   //ncfprintf (stderr,"chan %s connect to comp %d\n",findsym(chpnt->ctype),pnt->num);
   return chpnt;		/* return pointer to channel */
}

/*------------------------------------*/

chan *addchan_extern(comp *pnt, chan *chpnt)

/* Connect a channel or synapse to an external compartment */

/* 1) Add the channel to the external compartment's connections (because it will add current). */
/* 2) Set channel's compe pointer to point to the external compartment. */

{

  if (!pnt || !chpnt) return ((chan*)NULL);

   maklst(&pnt->clst, (conn*)chpnt);	/* add to comp's connections */
   chpnt->compe = pnt; 
   //ncfprintf (stderr,"chan %s connect to comp %d\n",findsym(chpnt->ctype),pnt->num);
   return chpnt;		/* return pointer to channel */
}

/*------------------------------------*/

chan *addchan_extern_ph(comp *pnt, chan *chpnt)

/* add a pH compartment to a channel's external voltage */
/* voltage in external comp is added in ccavoff() */

{
  chpnt->compp = pnt;
     return chpnt;
}

/*------------------------------------*/

comp *addcomp_extern(comp *pnte, comp *pnt)

/* Connect an internal compartment to a voltage offset from an external compartment */

/*  Set comparment's external voltage pointer to point to the external comp. */

{

  if (!pnt || !pnte) return ((comp*)NULL);

   pnt->pvext = pnte;
   //ncfprintf (stderr,"comp %s connect to comp %d\n",pnt->num,pnte->num);
   return pnt;		/* return pointer to internal comp */
}

/*------------------------------------*/

double vext(conn *chpnt) 

/* return a channel's external comp voltage for ephaptic effects */
{
    double rval;

  if (chpnt->compe!=NULL) {
     rval = chpnt->compe->v;
  }
  else rval = 0;
  return rval;
}
/*------------------------------------*/

void condense_comp_chan(comp *pnt)
               
/* Check channels in a compartment's connections. */
/* If a channel of same type already exists,     */
/*   compress it with the existing one. */
/* Can't add conductances when either old or new chan has */
/*   been saved, because they might be modified after creation */
{
   int n;
   conlst *cclpnt;
   register chan *chpnt;
   register conlst *lpnt;
   conlst *olpnt,*tlpnt;
   chan *cpnt;
   kcachan *k1pnt,*k2pnt;
   kcachan *k1pnts,*k2pnts;
   register int found, ctype, stype;

  n = 0;
  if (!pnt) return;
  for (olpnt=(conlst*)NULL,cclpnt=pnt->clst; cclpnt; ) {
    n++;
    if (!(chpnt=(chan*)cclpnt->conpnt)) {
	 ncfprintf(stderr,
	"condense_comp_chan: missing connection in comp %d conn %d\n",
				pnt->num,n); /* */
	olpnt = cclpnt;
	cclpnt= cclpnt->next;
        continue;
    }
    ctype = chpnt->ctype;
    stype = chpnt->stype;

    switch (ctype) {		/* if not a membrane channel, go to next */
	case  CA:
	case  NA:
	case  ClCa:
	case KCa:   break;	/* keep checking this channel further */
	case   K:   break;	/* keep checking this channel further */
				/*  in for() loop below */
      case CGMP:
      case GABA:
      case GLY:
      case AMPA:
      case KAINATE:
      case NMDA:
      case SYN2:
	default:   olpnt = cclpnt;
		   cclpnt= cclpnt->next;
		   continue;	/* don't condense synapses, etc. */
		   break;
    }

   for (found=0,lpnt=cclpnt->next; lpnt && !found; lpnt=lpnt->next) {

    if (!(cpnt=(chan*)lpnt->conpnt)) {
	/* ncfprintf(stderr,"condense_comp_chan2: missing connection...\n"); */
 	continue;
    }
    if (cpnt==chpnt) continue;		/* skip if same chan */

        if (cpnt->ctype==ctype) { 
	 if (cpnt->stype==stype) { 
#ifdef DEBUG
	       if (debug & NCMAK && debugz & 32)
  ncfprintf (stderr,"condense_comp_chan: comp %d chan type %s %d matches\n",
		pnt->num, findsym(ctype), stype);
#endif
	   if (cpnt->nchan && chpnt->nchan) {	/* check unitary conductance */
	          double unit1, unit2, uratio;
	       unit1 = cpnt->maxcond / cpnt->nchan;
	       unit2 = chpnt->maxcond / chpnt->nchan;
	       uratio = unit1/unit2;
	       if (uratio < .8 || uratio > 1.2) {
#ifdef DEBUG
	       if (debug & NCMAK && debugz & 32)
   ncfprintf (stderr,"condense_comp_chan: maxcond %g %g unitary %g %g not same\n",
			cpnt->maxcond,chpnt->maxcond,unit1,unit2);
#endif
		continue; 			/* not same u cond */
	       }
	   }
	  if (((ctype==CA) || (cpnt->vrev==chpnt->vrev)) &&
	    (cpnt->voffsm==chpnt->voffsm) &&
	    (cpnt->voffsh==chpnt->voffsh) &&
	    ((!cpnt->nchan) == (!chpnt->nchan))&& /* both same noise status ? */
	    (cpnt->arate==chpnt->arate) &&
	    (cpnt->brate==chpnt->brate) &&
	    (cpnt->crate==chpnt->crate) &&
	    (cpnt->drate==chpnt->drate) &&
	    (!cpnt->nocondens) &&
	    (!chpnt->nocondens) )
          found = 1;
#ifdef DEBUG
           else if (debug & NCMAK && debugz & 32) {
              ncfprintf (stderr,"condense_comp_chan: other parms not same\n");
              ncfprintf (stderr,"vrev %g voffsm %g voffsh %g nchan %g arate %g brate %g crate %g drate %g saved %hd\n",cpnt->vrev,cpnt->voffsm,cpnt->voffsh,cpnt->nchan,cpnt->arate,cpnt->brate,cpnt->crate,cpnt->drate,
cpnt->nocondens);
              ncfprintf (stderr,"vrev %g voffsm %g voffsh %g nchan %g arate %g brate %g crate %g drate %g saved %hd\n",chpnt->vrev,chpnt->voffsm,chpnt->voffsh,chpnt->nchan,chpnt->arate,chpnt->brate,chpnt->crate,chpnt->drate, chpnt->nocondens);

	   }
#endif
	 }  /* if (stype) */
	} /* if (ctype) */

     switch (chpnt->ctype) {

     case K:
       if (found) {

	   switch (chpnt->stype) {

         case 0:
         case 1:
         case 3:
         case 4:
         case 5:
         case 6:
         case 7:
         case 8:
	 default: break;

        }  /* switch (chpnt->stype) */
      }
	break;

     case ClCa:
	break;

     case KCa:
       if (found) {
	 switch (chpnt->stype) {

         case 0:
         case 2:
	      k1pnt = (kcachan*)cpnt;
	      k2pnt = (kcachan*)chpnt;
     	      if ((k1pnt->d1!=k2pnt->d1) ||
     	          (k1pnt->d2!=k2pnt->d2) ||
     	          (k1pnt->k1!=k2pnt->k1) ||
     	          (k1pnt->k2!=k2pnt->k2))
	      found = 0;
              break;
         case 1:
         case 3:
         case 5:
	      k1pnts = (kcachan*)cpnt;
	      k2pnts = (kcachan*)chpnt;
     	      if ((k1pnts->d1!=k2pnts->d1) ||
     	          (k1pnts->d2!=k2pnts->d2) ||
     	          (k1pnts->k1!=k2pnts->k1) ||
     	          (k1pnts->k2!=k2pnts->k2))
	      found = 0;
              break;

         case 4:
         case 6:  break;
	 default: break;

        }  /* switch (chpnt->stype) */
      }
	break;

     case NA:
	   switch (chpnt->stype) {

		default:
		case 0:
		case 1:
		case 2: 
		case 3: 
		case 4: 
		case 5: break;

	/*	case 10: found = 0; break;  */  /* if not condense */
	   }
	   break;
		
     case CA: 
	break;		/* cacomp ptr to chan is deleted in delchan() */
		
     case AMPA:		/* do nothing because we want these to be separate */
     case KAINATE:
     case CGMP:
     case GABA:
     case GLY:
     case CHRC:
     case NMDA:
     case SYN2:
     default:
	    found = 0;
            break;

     }  /* switch (chpnt->ctype) */

     if (found) break;

   }   /* for (found=0,lpnt;;) */

// found=0;
   if (found) {
	cpnt->maxcond += chpnt->maxcond;  /* add the conductance */
	cpnt->nchan += chpnt->nchan;      /* add number of channels */
	if (chaniss(cpnt->ctype,cpnt->stype)) maknois((chan *)cpnt);

		/* patch pointers in list */
	
	if (chpnt->last) chpnt->last->next = chpnt->next;
	if (chpnt->next) chpnt->next->last = chpnt->last;
	if (chanpnt==chpnt) chanpnt=(chan*)chpnt->next;
	if (chanend==chpnt) chanend=(chan*)chpnt->last;
	delchan (chpnt);	/* no need for new (or moved) channel conn */

	tlpnt = cclpnt;		/* erase this one */
	if (olpnt) {
	   olpnt->next = cclpnt->next;
	}
	else {
	    pnt->clst = cclpnt->next;
	}
	cclpnt=cclpnt->next;
	efree (tlpnt);		/* erase conlst pointer to chan */

   }  /* if (found) */ 

   else {
       olpnt = cclpnt;
       cclpnt = cclpnt->next;
   }
  } /* for (cclpnt;;) */

      /* Finally, set nchan to an integral number */

  for (lpnt=pnt->clst; lpnt; lpnt=lpnt->next ) {
    if (!(chpnt=(chan*)lpnt->conpnt)) {
	/* ncfprintf(stderr,"set_comp_nchan: missing connection...\n"); */
        continue;
    }
    setnchan(chpnt);
  } /* for (lpnt=pnt->clst;;) */
}

/*------------------------------------*/

void delchan(chan *chpnt)

/* Delete a channel and its substructures */

{
  if (!chpnt) return;
  if (chpnt->cstate) efree (chpnt->cstate);	/* free the random noise gen */ 
  if (chpnt->conc) efree (chpnt->conc);		/* free the state concs */
  if (chpnt->gvrevtab) efree (chpnt->gvrevtab);	/* free the ion table */
  if (chpnt->gfractab) efree (chpnt->gfractab);	/* free the ion table */
  if (chpnt->ions) efree (chpnt->ions);		/* free the ion table */
  efree (chpnt);
  cumchan--;
}

/*------------------------------------*/

cacomp *makcacomp(elem *epnt, cattrib *apnt, comp *comp1, double carea)

{
      int i, cashell, caoshell, mincoresiz, nshell, nshello;
      double cai, tcai, cao, vrev, *cshp;
      double dia, r, dr, dca, dcaoo, core, cabnd, cabnda, caobnda;
      double pkm, vmax, kex, area, cmarea;
      chantype *chtyp;
      cacomp *capnt;

    if ((capnt=(cacomp *)emalloc(sizeof(cacomp))) == (cacomp*)NULL) {
      ncfprintf (stderr,"makcacomp: no space left for ca comp %d\n",cumcacomp+1);
      return ((cacomp*)NULL);  
    }
    capnt->ctype = CACOMP;
    capnt->num = cumcacomp++;
    capnt->comp1 = comp1;		/* set ca comp pntr to voltage comp */
    comp1->capnt = capnt;		/* set voltage comp pntr to ca comp */

    cao=apnt->cao;
    cai=apnt->cai;
    tcai=apnt->tcai;
    vrev=apnt->vrev;

    chtyp = getchantype(CA,0);		/* assume Ca type 0 for cai calc */

						/* If calcium not determined */
    if (cao==NULLVAL || cai==NULLVAL) {		/*  must hunt for info */
      if (vrev!=NULLVAL) {			/* Find ca conc from vrev */
        if (cai==NULLVAL) {
		double df, nf, cavoff;
	   if (cao == NULLVAL) cao = dcao;

		/* The idea here is to apply the GHK voltage equation */
		/*  backwards to calculate a value for cai. */

	   cavoff = qcavoff(chtyp->parm) * dcaspvrev;
           df = exp ( FR * (vrev-cavoff) / ktemp); /* multiplier for Ca perm */
	   nf = 1.0 / ( 1.0 + df);
           cai = ((4.0*cao*nf  + ghkvo(chtyp,chtyp->ions,vrev)) /
		exp((vrev-cavoff)*FR/ktemp) - ghkvi(chtyp,chtyp->ions,vrev)) / (4.0*nf*df);
           /* cai = cao / exp(vrev*f2rt) - dpki; */
	   if (cai < 10e-9) cai = 10e-9;
        }
        else if (cao==NULLVAL) {
          cao = dcao;				/* always use default cao */
        }	
		/* Note: when both vrev and cai are set by user */
		/*   then vrev and cai/cao will be inconsistent */
		/*  If so, then Ca current for [Ca] will be derived from cai */
		/*   and total Ca chan current will be derived from vrev. */
		/*  This represents the case when [K]i affects total Ca chan */
		/* current and thus the channel's reversal potential. */
      }
      else {  /* if (vrev==NULLVAL) */	      /* vrev unknown, use defaults */
        if (cao==NULLVAL) {
          cao = dcao;
        }
        if (cai==NULLVAL) {
          cai = dcai;
        }
      }
    } 
    if (tcai==NULLVAL) {			/* set threshold [Ca]i */
        tcai = dtcai;
	if (tcai <=0) tcai = 1e-9;
    }
    if (cai<=0.0) {
	cai = tcai;				/* limit cai to threshold */
    }
    capnt->cao  = cao;
    capnt->cai  = cai;
    capnt->tcai = tcai;

    if ((nshell=apnt->cashell)==NULLVAL) nshell=dcashell; /* num of shells */
    if (nshell < 1) {
	nshell = 1;	/* min of 1 shell */
/*	ncfprintf 
       (stderr,"# makcacomp: num of ca shells is zero or neg, setting to %d\n",
			nshell); /* */
    }
    if ((nshello=apnt->caoshell)==NULLVAL) nshello=dcaoshell; /* num of shells outside*/

    if ((area=apnt->sarea) == NULLVAL) area = dcasarea;    /* area of shell for caflux in um2 */
    if (area > 0) {
      area *= DMUM*DMUM;	    /* area was in um2 */ 
      r = sqrt (area / (4.0*MPI)); /* radius in dm (decimeters) */
    } 
    else 			    // find area from size of neural element

    /* Try to get best est. for radius and area. */
    /*  Area comes in as cm2 (from cap measurement) */
    /* Calib in dm (decimeters) because molarity is Moles/liter = M/1000 cm3 */

				/* diameter in microns */
    if   ((epnt->ctype==SPHERE) && ((dia=((sphere *)epnt)->dia)!=NULLVAL)) {
         r = dia * 0.5 * DMUM;		/* radius in dm (decimeters) */
         area = (4.0 * MPI * r*r);	/* shell surface (dm2) */
    }

    else if ((epnt->ctype==CABLE) && carea>0) { /* if cable with non-zero area */
          area = carea * DMCM2;
          r = sqrt (area / (4.0*MPI)); /* radius in dm (decimeters) */
    }
    else {				/* find radius from area (capacitance) */
      if (carea>0) {			/* in this case area comes in as cm2 */
          area = carea * DMCM2;
          r = sqrt (area / (4.0*MPI)); /* radius in dm (decimeters) */
       }
       else {
          r = CACOMPRAD * DMUM;		/* (calib in dm) make it 5 um radius */
          area = (4.0 * MPI * r*r);	/* shell surface (dm2) */
       }
    }
    // fprintf (stderr,"area %g\n",area);

    capnt->area = area;			/* save int area for redocacomp() */
    cmarea = area / DMCM2;		/* calibrate area in cm2 */

    mincoresiz = 3;
    dr = CASHD*DMUM;			/* thickness of shell in dm (1e-6) */
    if (r < (dr * (nshell+mincoresiz))) {  /* if sphere is too small for dr */
      while (r < (dr * (nshell+mincoresiz))) { 
        dr = dr * 0.5;			/* make smaller dr */
        if (dr < 1e-10) {
 	   ncfprintf (stderr,"# nc: makca: dr reduced too much, %g\n",dr);
	   break;
        }
      } 
   /* ncfprintf (stderr,"# nc: makca: nshells = %d, thickness %g um at node %s\n",
				nshell, dr/DMUM,
			 prnode (epnt->node1a,
				 epnt->node1b,
				 epnt->node1c,
				 epnt->node1d));
    */
    }
    capnt->dr = dr/DMUM;		/* save dr for printconns() */
    dca = ddca * 1e-2;			/* Diffusion const for Ca in dm2/sec*/
    dca *= exp(log(dqdc)*(tempcel-dbasetdc)/10);/* Q10 for diffusion */

    dcaoo = ddcao * 1e-2;		/* Diffusion const for Ca in dm2/sec*/
    dcaoo *= exp(log(dqdc)*(tempcel-dbasetdc)/10);/* Q10 for diffusion */

    cashell = nshell + 1;		/* core is extra "shell" */
    capnt->cashell = cashell;		/* assign pointer to ca shells */

					/* allocate space for shells & core */

    if ((cshp=(double *)emalloc(cashell*sizeof(double))) == (double*)NULL) {
      ncfprintf(stderr,"makcacomp: no space left for %d calcium shells\n",nshell);
      return ((cacomp*)NULL);  
    }
    capnt->cais = cshp;			/* assign pointer to ca shells */

    for (i=0; i<cashell; i++) {
       capnt->cais[i]  = cai;		/* set initial Ca concentration */
    }

    caoshell = nshello + 1;		/* external space is extra "shell" */
    capnt->caoshell = caoshell;		/* assign pointer to ca shells */
					/* allocate space for external shells & space */

    if (nshello>0) {
      if ((cshp=(double *)emalloc(caoshell*sizeof(double))) == (double*)NULL) {
        ncfprintf(stderr,"makcacomp: no space left for %d external calcium shells\n",nshello);
        return ((cacomp*)NULL);  
      }
      capnt->caos = cshp;			/* assign pointer to ca shells */
    }
    else capnt->caos=NULL;		/* no outside shells, use cao instead */

    if (capnt->caos!=NULL) {
      for (i=0; i<caoshell; i++) {
        capnt->caos[i]  = cao;		/* set initial Ca concentration */
      }
    }

    if (apnt->cabuf>0) {
	 double bkd, cabf, cabt, cabti, *cab, carate;

      if ((cabt=apnt->btot) ==NULLVAL) {
	   cabt = dcabt; /* Tot cabuf conc */
      	   if ((cabti=apnt->btoti)==NULLVAL) cabti = dcabti; /* Tot cabuf conc */
      }
      else {
	 if ((cabti=apnt->btoti)==NULLVAL) cabti = cabt;   /* tot specfied, but not 1st shell */
      }
      capnt->cabt = cabt;
      capnt->cabti = cabti;
						/* set rates */
      if ((cabf=apnt->bmax)==NULLVAL) {
         carate = exp(log(dqcab)*(tempcel-dbasetca)/10.0);
	 cabf = dcabf;				/* Forw rate for cabuf */
	 if (cabf==0) cabf=1;
      }
      else carate = 1.0;
      if ((bkd=apnt->bkd) ==NULLVAL) bkd = dcabr/cabf; 	/* Kd for cabuf */

		/* if vmax rate not set, use default values with temp coeff */

      capnt->cabf = cabf * timinc * carate;
      if (bkd==0) bkd = 1;
      capnt->cabr = bkd * cabf * timinc * carate;

      if ((cab=(double *)emalloc(cashell*sizeof(double))) == (double*)NULL) {
    ncfprintf(stderr,"makcacomp: no space left for %d calcium buffer\n",nshell);
        return ((cacomp*)NULL);  
      }
      capnt->cab = cab;			/* assign pointer to ca shells */
      if (capnt->cab) {
        for (i=0; i<cashell; i++) {
          if (i==0) capnt->cab[i]  = cabti;		/* Ca buffer conc */
          else capnt->cab[i]  = cabt;			/* Ca buffer conc */
        }
      }
      cabnd = 1.0;
      cabnda = 1.0;
      caobnda = 1.0;
    }
    else {				/* no cabuf */
	capnt->cabt = 0;
	capnt->cabf = 0;
	capnt->cabr = 0;
        capnt->cab = (double *)NULL;	/* assign pointer to ca shells */
        if ((cabnd=apnt->cabnd)==NULLVAL) cabnd=dcabnd; /* ratio bnd to free */
        cabnda  = 1.0 / (cabnd+1);	   /* multiplier is one more */
        caobnda = 1.0 / (cabnd+1);   		   /* multiplier is one more */
    } 
    capnt->cabnd = cabnd;		/* save for possible redocacomp() */

    if (apnt->cicr>0) {			/* set up runtime CICR */
	 double cas, vm2, vm3, ncicr, mcicr, pcicr, kacicr, kfcicr, krcicr, k1cicr, k2cicr, c1cicr;

        if ((cas=apnt->cas)==NULLVAL) cas = CASSTART; /* CICR store starting conc */
	capnt->cas = cas; 
   	if ((vm2=apnt->vm2)==NULLVAL) vm2 = vm2CICR; /* CICR store uptake vmax */
	capnt->vm2 = vm2; 
   	if ((vm3=apnt->vm3)==NULLVAL) vm3 = vm3CICR; /* CICR store release vmax */
	capnt->vm3 = vm3; 
   	if ((ncicr=apnt->ncicr)==NULLVAL) ncicr = nCICR;  /* CICR store uptake Hill Coeff */
	capnt->ncicr = ncicr; 
   	if ((mcicr=apnt->mcicr)==NULLVAL) mcicr = mCICR;  /* CICR store release Hill Coeff */
	capnt->mcicr = mcicr; 
   	if ((pcicr=apnt->pcicr)==NULLVAL) pcicr = pCICR;  /* CICR store Hill Coeff */
	capnt->pcicr = pcicr; 
   	if ((kacicr=apnt->kacicr)==NULLVAL) kacicr = kaCICR;  /* thresh CICR store release */
	capnt->kacicr = kacicr; 
   	if ((kfcicr=apnt->kfcicr)==NULLVAL) kfcicr = kfCICR;  /* passive leak CICR store */
	capnt->kfcicr = kfcicr; 
   	if ((krcicr=apnt->krcicr)==NULLVAL) krcicr = krCICR;  /* thresh CICR store release */
	capnt->krcicr = krcicr; 
	if ((k1cicr=apnt->k1cicr)==NULLVAL) k1cicr = k1CICR;  /* thresh CICR store uptake */
	capnt->k1cicr = k1cicr;
   	if ((k2cicr=apnt->k2cicr)==NULLVAL) k2cicr = k2CICR;  /* thresh CICR store uptake */
	capnt->k2cicr = k2cicr; 
   	if ((c1cicr=apnt->c1cicr)==NULLVAL) c1cicr = c1CICR;  /* ratio ER / cytoplasmic vol */
	capnt->c1cicr = c1cicr; 
    } else {
	capnt->vm2 = capnt->vm3 = 0; 
    }

    if (apnt->ip3>0) {			/* set up runtime IP3 */

	double cas2, ip3i, bip3, vip3, b2ip3, a2ip3, v2ip3, v3ip3, v4ip3, k3ip3; 
	double d1ip3, d2ip3, d3ip3, d4ip3, mip3, hip3, a3ip3, oip3, mtypeip3; 


     if ((cas2=apnt->cas2)==NULLVAL) cas2 = CAS2START; /* IP3 Ca store starting conc */
     capnt->cas2 = cas2; 
     if ((ip3i=apnt->ip3i)==NULLVAL) ip3i = IP3ISTART; /* initial [IP3] conc */
     capnt->ip3i = ip3i; 
     if ((bip3=apnt->bip3)==NULLVAL) bip3 = betaIP3; /* IP3 frac rate */
     capnt->bip3 = bip3; 
     if ((vip3=apnt->vip3)==NULLVAL) vip3 = v1IP3;   /* init IP3 store flux */
     capnt->vip3 = vip3; 
     if ((b2ip3=apnt->b2ip3)==NULLVAL) b2ip3 = b2IP3; /* IP3 store rev const for Tau M */
     capnt->b2ip3 = b2ip3; 
     if ((a2ip3=apnt->a2ip3)==NULLVAL) a2ip3 = a2IP3; /* IP3 store forward const for Tau M */
     capnt->a2ip3 = a2ip3;

     /* constants taken from Rinzel & Li, 1994 */

     if ((a3ip3=apnt->a3ip3)==NULLVAL) a3ip3 = a3IP3; /*IP3 store const IP3 Tau H,2e5 1/ms */
     capnt->a3ip3 = a3ip3; 
     if ((oip3=apnt->oip3)==NULLVAL) oip3 = oIP3; /*Hill coeff, coop pump binding ip3 store uptake, 2 */
     capnt->oip3 = oip3;
     if ((v2ip3=apnt->v2ip3)==NULLVAL) v2ip3 = v2IP3; /* IP3 store rate for Ca release,1.11/s */
     capnt->v2ip3 = v2ip3; 
     if ((v3ip3=apnt->v3ip3)==NULLVAL) v3ip3 = v3IP3; /*IP3 store const for Ca uptake, 0.9/s */
     capnt->v3ip3 = v3ip3;
     if ((v4ip3=apnt->v4ip3)==NULLVAL) v4ip3 = v4IP3; /*IP3 store const for Ca uptake, 7.11/s */
     capnt->v4ip3 = v4ip3;
     if ((k3ip3=apnt->k3ip3)==NULLVAL) k3ip3 = k3IP3; /*IP3 store kd for Ca uptake, 0.1/s*/
     capnt->k3ip3 = k3ip3;
     if ((d1ip3=apnt->d1ip3)==NULLVAL) d1ip3 = d1IP3; /*IP3, for det Q value for h gate, 0.13e-6 M */
     capnt->d1ip3 = d1ip3;
     if ((d2ip3=apnt->d2ip3)==NULLVAL) d2ip3 = d2IP3; /*IP3, for det Q value for h gate, 1.049e-6 M */
     capnt->d2ip3 = d2ip3;
     if ((d3ip3=apnt->d3ip3)==NULLVAL) d3ip3 = d3IP3; /*IP3, for det Q value for h gate, 0.9434e-6 M */
     capnt->d3ip3 = d3ip3;
     if ((d4ip3=apnt->d4ip3)==NULLVAL) d4ip3 = d4IP3; /*IP3, for det m gate at infinity, 0.82e-6 M */
     capnt->d4ip3 = d4ip3;
     if ((mip3=apnt->mip3)==NULLVAL) mip3 = MIP3START; /*IP3 store uptake m gate fraction*/ 
     capnt->mip3 = mip3;
     if ((hip3=apnt->hip3)==NULLVAL) hip3 = HIP3START; /*IP3 store uptake h gate fraction*/
     capnt->hip3 = hip3;
     if ((mtypeip3=apnt->mtypeip3)==NULLVAL) mtypeip3 = mtypeIP3;/* 0=>set to static (def, use equil inf m value) */
     capnt->mtypeip3 = mtypeip3;
    } else {
	capnt->vip3 = capnt->bip3 = capnt->v2ip3 = capnt->v3ip3 = 0; 
    }

    switch (epnt->ctype) {

    case CA:
    case SPHERE:
    default:			 /* shell factors for diffusion eqn. */

	capnt->casf0 = timinc*F2/(area*dr);  		/* caflux */
	capnt->casfn = dca * timinc/(dr*dr)*cabnda;	/* norm shells */
	core = r - nshell * dr;
	capnt->casfc = dca * timinc *3.0/(dr*core)*cabnda; /* D*t*SA/(dr*vol) */
							  /*  for core shell */
	capnt->casfno = dcaoo * timinc/(dr*dr)*caobnda;	/* norm outside shells */
        r = ddiacao * 0.5 * DMUM;			/* radius in dm (decimeters) */
	core = r - nshello * dr;
	capnt->casfco = dcaoo * timinc *3.0/(dr*core)*caobnda; /* D*t*SA/(dr*vol) */
							  /*  for core shell */
	//fprintf(stderr,"cabnda %g casfn %g casfc %g caobnda %g casfno %g casfco %g\n",
	//		cabnda,capnt->casfn,capnt->casfc,caobnda,capnt->casfno,capnt->casfco);
	break;

    case CABLE:			/* must fix this for cylindrical shells */

	capnt->casf0 = timinc*F2/(area*dr);		/* caflux (t/M/sec/A/V*/
	capnt->casfn = dca * timinc/(dr*dr)*cabnda; 	/* norm shells */
	core = r - nshell * dr;
	capnt->casfc = dca * timinc *3.0/(dr*core)*cabnda; /* D*t*SA/(dr*vol) */
							  /*  for core shell */
	capnt->casfno = dcaoo * timinc/(dr*dr)*caobnda;   /* norm shells */
        r = ddiacao * 0.5 * DMUM;			  /* radius in dm (decimeters) */
	core = r - nshello * dr;
	capnt->casfco = dcaoo * timinc *3.0/(dr*core)*caobnda; /* D*t*SA/(dr*vol) */
	break;
    }

    if (apnt->cabuf) docacompi(capnt);

    if (apnt->pump) {
      if ((pkm=apnt->pkm)==NULLVAL) pkm = dcapkm;    /* Km for ca pump */
      if ((vmax=apnt->vmax)==NULLVAL) {
		vmax = dcavmax; /* Vmax for ca pump */
      }
      vmax *= exp(log(dqcavmax) * (tempcel - dbasetca)/10.0);
      capnt->pkm   = pkm;
      capnt->vmax = vmax * cmarea;
    }
    else capnt->vmax = 0;

    if (apnt->exch) {		/* precalculate values for Na-Ca exchanger  */
				/* From de Schutter & Smolen 1998  */
				/* after DiFrancesco & Noble 1985  */
				/* and Gabbiani,Midtgaard and  Knopfel 1994 */

      if ((kex=apnt->kex)==NULLVAL) kex = dcakex; /* Exch rate for ca exch */
      capnt->kexi  = kex  * cmarea *     dnao*dnao*dnao*1e12;
      //capnt->kexo  = kex  * cmarea * cao*dnai*dnai*dnai*1e12;  // must multiply by caos[0] 
      capnt->kexo  = kex  * cmarea *     dnai*dnai*dnai*1e12;  // multiply by caos[0] in docacomp()
      capnt->ego = EG*FR/ktemp;
      capnt->egi = (EG-1.0)*FR/ktemp;
    }
   else {
      capnt->kexi = 0;
      capnt->kexo = 0;
   } 
   vrev = r2ft * log(cao/cai);   /* Nernst equation for Ca */
   capnt->vrev = vrev;
				/* surface potential due to Ca */
   capnt->cavoff = dcaspvrev* qcavoff(chtyp->parm);
  return capnt;
}

/*------------------------------------*/

void redocacomp(comp *comp1, double narea)

{
      int cashell, mincoresiz, nshell, nshello;
      double area, arearatio, r, dr, dca, dcaoo, core, cabnd, caobnd;
      cacomp *capnt;

    /* Keep track of area, estimate new size from area increment. */
    /* Calib r,area in dm (decimeters) because molarity is Moles/liter = M/1000 cm3*/

    if (!(capnt=comp1->capnt)) {
      ncfprintf (stderr,"redocacomp: no cacomp to redo %d\n",cumcacomp);
      return;  
    }
    nshell = capnt->cashell-1;	/* use original number of ca shells */
    nshello = capnt->caoshell-1;/* use original number of ca shells */
    cabnd  = capnt->cabnd;	/* get original ca bound */
    caobnd  = 1.0;		/* get original ext ca bound */
    area = capnt->area + narea;	/* find size of compartment */
    r = sqrt (area / (4.0*MPI)); /* radius in dm (decimeters) */
    mincoresiz = 3;
    dr = CASHD*DMUM;			/* thickness of shell in dm (1e-6) */
    if (r < (dr * (nshell+mincoresiz))) {  /* if sphere is too small for dr */
      while (r < (dr * (nshell+mincoresiz))) { 
        dr = dr * 0.5;			/* make smaller dr */
        if (dr < 1e-10) {
 	   ncfprintf (stderr,"# nc: redoca: dr reduced too much, %g\n",dr);
	   break;
        }
      } 
   /* ncfprintf (stderr,"# nc: redocacomp: nshells = %d, thickness %g um at node %s\n",
				nshell, dr/DMUM,
			 prnode (epnt->node1a,
				 epnt->node1b,
				 epnt->node1c,
				 epnt->node1d));
    */
    } 
    capnt->dr = dr/DMUM;		/* save dr for printconns() */
    dca = ddca * 1e-2;			/* Diffusion const for Ca in dm2/sec*/
    dca *= exp(log(dqdc)*(tempcel-dbasetdc)/10);/* Q10 for diffusion */

    dcaoo = ddcao * 1e-2;		/* Diffusion const for Ca in dm2/sec*/
    dcaoo *= exp(log(dqdc)*(tempcel-dbasetdc)/10);/* Q10 for diffusion */

    cashell = nshell + 1;		/* core is extra "shell" */
    capnt->cashell = cashell;		/* assign number of ca shells */

    capnt->casf0 = timinc*F2/(area*dr);  	/* caflux */
    capnt->casfn = dca * timinc/(dr*dr)*cabnd;	/* norm shells */
    core = r - nshell * dr;
    capnt->casfc = dca * timinc *3.0/(dr*core)*cabnd; 	/* D*t*SA/(dr*vol) */
							/*  for core shell */
    capnt->casfno = dcaoo * timinc/(dr*dr)*caobnd;	/* norm shells */
    core = ddiacao - nshello * dr;
    capnt->casfco = dcaoo * timinc *3.0/(dr*core)*caobnd; 	/* D*t*SA/(dr*vol) */
							/*  for core shell */
    arearatio = area / capnt->area;

    capnt->vmax *= arearatio;
    capnt->kexi *= arearatio;
    capnt->kexo *= arearatio;
    capnt->area = area;
				/* using CICR, add vars here scaled with area of Ca comp */
    if (capnt->vm2>0) {                   
        capnt->vm2 *= arearatio; 
        capnt->vm3 *= arearatio;
    }
    if (capnt->bip3>0) {                  /* using IP3 release */
        capnt->bip3  *= arearatio;
        capnt->vip3  *= arearatio;
        capnt->v2ip3 *= arearatio;
        capnt->v3ip3 *= arearatio;
        capnt->v4ip3 *= arearatio;
    }
}

/*------------------------------------*/

void delcacomp (cacomp *capnt)

/* delete a cacomp and all its shells and buffers */

{
  if (!capnt) return; 
  if (capnt->cais) efree (capnt->cais);	/* erase shells */
  if (capnt->caos) efree (capnt->caos);	/* erase shells */
  if (capnt->cab)  efree (capnt->cab);  /* erase buffers */
  efree (capnt);			/* erase cacomp */
  cumcacomp--;
}

/*------------------------------------*/

chan *makca(elem *epnt, cattrib *apnt, comp *comp1, double area, int saved)
    	           
/* Make a new calcium channel, and if it doesn't have one,
   make a new calcium compartment, too.  
   Area calibrated in cm2.
*/

{
    chan *chpnt,*addchan(comp *pnt, chan *chpnt);
    cacomp *capnt;
    double cao, cai, vrev;
    double maxcond;

#ifdef DEBUG
 if (debug & NCMAK && debugz & 16)
   ncfprintf (stderr,"makca: comp %d\n",comp1->num);
#endif
		/* don't make chan if maxcond 0 */ 

  if ((maxcond=apnt->maxcond)==NULLVAL) maxcond = dmaxca;

  if (!(chpnt=makchan(apnt,comp1,(comp*)NULL,saved))) {   /* make channel */
      return ((chan *)NULL); 
  }
  else {
   if (!(capnt=comp1->capnt)) { 
     if ((capnt=makcacomp(epnt,apnt,comp1,area))==NULL) { /* make new ca comp */
       ncfprintf (stderr,"makca: no space for ca comp %d\n",comp1->num);
       return ((chan*)NULL);
     }
   }

   if (!chpnt->setvrev) {
      if (capnt->caos!=NULL) cao = capnt->caos[0];
      else                   cao = capnt->cao;
      cai = capnt->cais[0];
      vrev = ghkv(chpnt->chtyp,cao,cai);   /* GHK V equation, Cai + Ki */
      chpnt->vrev = vrev;
   }

			/*  end of ca comp stuff */
#ifdef DEBUG
  if (debug & NCMAK && debugz & 16)
    ncfprintf (stderr,"makca end.\n");
#endif

   return chpnt;
 }  /* else (chpnt) */
}

/*------------------------------------*/

int findphotrec(int num1, int num2, int num3, int num4, photrec **rpnt, const char *str)
{
   photrec *pnt;

  if (num4 < 0) {
    for (pnt=recpnt; pnt; pnt=(photrec *)pnt->next) {
      if ((pnt->recnm1 == num1) && (pnt->recnm2 == num2) && 
  	  (pnt->recnm3==num3)) {
	  *rpnt = pnt;
	  return 1;
      }
    }
  } else {
    for (pnt=recpnt; pnt; pnt=(photrec *)pnt->next) {
      if ((pnt->recnm1 == num1) && (pnt->recnm2 == num2) && 
	  (pnt->recnm3==num3) && (pnt->recnm4==num4)) {
	  *rpnt = pnt;
	  return 1;
      }
    }
  }
  if (str) {
     if (num4 != NULLVAL)
  	 ncfprintf(stderr,"\n%s: can't find recep %d %d %d %d\n",
			str,num1,num2,num3,num4);
     else if (num3 != NULLVAL)
  	 ncfprintf(stderr,"\n%s: can't find recep %d %d %d\n",str,num1,num2,num3);
     else if (num2 != NULLVAL)
  	 ncfprintf(stderr,"\n%s: can't find recep %d %d\n",str,num1,num2);
     else 
  	 ncfprintf(stderr,"\n%s: can't find recep %d\n",str,num1);
   }
  return 0; 
}

/*------------------------------------*/

int findphotrec(int num1, int num2, int num3, photrec **rpnt, const char *str)

{
  return findphotrec(num1, num2, num3, NULLVAL, rpnt, str);
}

/*------------------------------------*/

void maklst(conlst **head, conn *cpnt)

/* add an entry to a connection list. */

{
    conlst *lpnt,*tmp,*oldtmp;
    int found;

  if (!head) { 
     ncfprintf (stderr,"maklst: incorrect list pointer\n");
     return;  
  }
  found = 0;
  for (oldtmp=tmp= *head; tmp; oldtmp=tmp,tmp=tmp->next){ /* find end of list */
     if (tmp->conpnt==cpnt) {
        found = 1;			/* check to see if already there */
        break;
     }
  }
  if (!found) {
 
    if ((lpnt=(conlst *)emalloc(sizeof(conlst))) == NULL) {
       ncfprintf (stderr,"maklst: no space left for conlist\n");
       return;  
    }
    if (oldtmp) oldtmp->next = lpnt;
    else *head = lpnt;
 
    lpnt->next = (conlst *)NULL;
    lpnt->conpnt = cpnt;		/* set pointer to connection */
    cumclst++;
  }
}

/*------------------------------------*/

void addlst(conlst **head, conlst *h2)

/* add one connection list (h2) to another (head). */

{
    conlst *tmp;

  if (!head) { 
     ncfprintf (stderr,"addlst: incorrect list pointer\n");
     return;  
  }
  for (tmp=h2; tmp; tmp=tmp->next){ /* find end of list */
      maklst (head, tmp->conpnt);
  }
}

/*------------------------------------*/

void dellst(conlst **head, conn *cpnt)

/* delete an entry from a connection list. */

{
    conlst *tmp,*oldtmp;
    int found;

  if (!head) { 
     ncfprintf (stderr,"lst: incorrect list pointer\n");
     return;  
  }
  if (! *head) return;			/* empty list: nothing to do */
  else if (! cpnt) return;
  else {
    for (found=0,oldtmp=tmp= *head; tmp; oldtmp=tmp,tmp=tmp->next) { 
       if (tmp->conpnt==cpnt) {
          found = 1;
          break;
       }
    }
    if (found) {
      oldtmp->next = tmp->next;
      if (*head==tmp) *head = (conlst*)NULL;
      efree (tmp);
      cumclst--;
    }
  }    /* else */
}

/*------------------------------------*/


void makdyad(dyadlst **head, synap *spnt)

/* add an entry to a synapse's dyad list */

{
    dyadlst *lpnt,*tmp,*oldtmp;

  if (!head) { 
     ncfprintf (stderr,"makdyad: incorrect list pointer\n");
     return;  
  }
  if ((lpnt=(dyadlst *)emalloc(sizeof(dyadlst))) == (dyadlst*)NULL) {
     ncfprintf (stderr,"dyadlst: no space left for dyadlist\n");
     return;  
  }
  if (! *head) *head = lpnt;	/* save head if first one */
  else {
    for (oldtmp=tmp= *head; tmp; tmp=tmp->next)  /* find end of list */
       oldtmp = tmp;			   /*  and get pointer to last */
    if (oldtmp) oldtmp->next = lpnt;
  }
 
  lpnt->sdyad = spnt;		/* set pointer to connection */
  lpnt->sdyadnum = spnt->num;
  lpnt->next = (dyadlst*)NULL;
  cumdlst++;
}

/*------------------------------------*/

int numdyad(dyadlst *head)

/* Return the number of dyads in a dyad list. */

{
    dyadlst *lpnt;
    int ndyad;

  if (!head) { 
     ncfprintf (stderr,"numdyad: incorrect list pointer\n");
     ndyad=0;  
  }
  else {
    for (ndyad=0,lpnt=head; lpnt; lpnt=lpnt->next,ndyad++) { }
  } 
  return ndyad;
}

/*------------------------------------*/

void deldyad(dyadlst **head)

/* delete a dyad list (but not the synapses it points to.) */

{
    dyadlst *lpnt,*nlpnt;

  if (!head) { 
     ncfprintf (stderr,"deldyad: incorrect list pointer\n");
     return;  
  }
  if (! *head) return;			/* empty list: nothing to do */
  else {
    for (lpnt=*head; lpnt; lpnt=nlpnt) { 
      nlpnt = lpnt->next;
      efree (lpnt);
      cumdlst--;
    }
   *head = (dyadlst*)NULL;
  } 
}

/*------------------------------------*/
