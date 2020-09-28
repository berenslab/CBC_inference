/* ncelem.cc */

/* Routines to make and access "neural elements" */

#include "nc.h"
#include "y.tab.h"
#include "ncsub.h"
#include "ncelem.h"
#include "control.h"

// #define EHASHSIZ 509             /* prime number to make better hash */
 #define EHASHSIZ 4919             /* prime number to make better hash */
// #define EHASHSIZ 19937             /* prime number to make better hash */
// #define EHASHSIZ 99929             /* prime number to make better hash */
// #define EHASHSIZ 1999771              /* prime number to make better hash */

elem *ehashtab[EHASHSIZ] = {0};    /* element elnum hash table */
elem *ehashtabe[EHASHSIZ] = {0};   /* element elnum hash table end */
elem *enhashtab[EHASHSIZ] = {0};   /* element node hash table */
elem *enhashtabe[EHASHSIZ] = {0};  /* element node hash table end */

int ehashsiz = EHASHSIZ;

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
}
#endif

#include "ncio.h"

// #define DEBUG

#ifdef DEBUG
#include "ncdebug.h"
#endif

extern elem *elempnt;
extern elem *elemend;
extern int cumelem;
extern int cumsynap;
extern int delelmt;
extern int cumattr;
extern int cumcattr;
extern int cumnattr;
extern int interp;

extern elem *elpnt;

char *emalloc(unsigned int n);
void execerror(const char *s, const char *t);
void checkelemnode(elem *elpnt);
double *copyfiltarr (int siz, double *oarr);
void copyattr(attrib *src, attrib *dest);
void efree(void *ptr);
chattrib *chanattr(elem *elpnt, short int ctype);
void delelem (elem *epnt);
double getval(const char *s);
char *prnode (int n1, int n2, int n3, int n4);
char *prnode (node *npnt);
elem *foreachn2 (elem *epnt);
elem *setforeachn2 (int etype, int na, int nb, node **ppn);

/*---------------------------------------------------*/

/*
void xdebugf(void)

{
   ncfprintf (stderr,"ehashtab %d\n",ehashtab[5]);
}
*/

/*---------------------------------------------------*/

void einithash(void)
{
   int i;

   for (i=0; i<ehashsiz; i++) {
     ehashtab[i]  = (elem*)NULL;
     ehashtabe[i] = (elem*)NULL;
   }
}

/*---------------------------------------------------*/

void eninithash(void)
{
   int i;

   for (i=0; i<ehashsiz; i++) {
     enhashtab[i]  = (elem*)NULL;
     enhashtabe[i] = (elem*)NULL;
   }
}

/*---------------------------------------------------*/

int elemhash(int num)

/* make disordered index of element num */

{
   register int i;

   i = num;
   if (i<0) i = -i;
   return (i % ehashsiz); 
}

/*---------------------------------------------------*/

int nelemhash(nodeint nodea, nodeint nodeb)

/* make disordered index of element using (ct,cn) node numbers */

{
   register int i;

   i = 1;
//   if (nodea > 1) i *= nodea;
//   if (nodeb > 1) i *= nodeb;
   i += (nodea+1) * 10000;
   i += nodeb;
   if (i<0) i = -i;
   return (i % ehashsiz); 
}

/*---------------------------------------------------*/

elem *getepnt (nodeint nodenm1, nodeint nodenm2)

/* get starting elem pointer for (ct,cn) hash */

{
  return (enhashtab[nelemhash(nodenm1,nodenm2)]);
}

/*---------------------------------------------------*/

void einstall (elem *nelem)
              
/* Install an element in elnum hash table */
/*  Uses "direct chaining" with field "hnext". */

{
    register int i;
    register elem *epnt,*elast;
    static int einitfl=0;

   if (!nelem) return;

   if (!einitfl) {			/* initialize table once at start */
        einitfl = 1;
        einithash();
   }
   nelem->hnext = NULL;
   nelem->hlast = NULL;
   i=elemhash(nelem->elnum); 		/* initial index into ehashtab*/ 
   if (!(epnt=ehashtab[i])) {            /* install directly in table */
      ehashtab[i]  = nelem;
      ehashtabe[i] = nelem;
   }
   else {                               /* otherwise, go to end of list */
//      for (; epnt; epnt=epnt->hnext)    /* find last element */
//          elast = epnt;    
//      elast->hnext = nelem;            /* put new symbol at end of list */
//      nelem->hlast = elast;

      elast = ehashtabe[i];
      elast->hnext = nelem;            /* put new symbol at end of list */
      nelem->hlast = elast;
      ehashtabe[i] = nelem;
   }
}

/*---------------------------------------------------*/

void eninstall (elem *nelem)
              
/* Install an element in elem *node* hash table */
/*  Uses "direct chaining" with field "hnnext". 
    Hash the first 2 node numbers (ct,cn) only. 
*/

{
    register int i;
    register elem *epnt,*elast;
    static int eninitfl=0;

   if (!nelem) return;
   // if (nelem->ctype!=CABLE) return;

   if (!eninitfl) {			/* initialize table once at start */
        eninitfl = 1;
        eninithash();
   }
   nelem->hnnext = NULL;
   nelem->hnlast = NULL;
   i=nelemhash(nelem->node1a,nelem->node1b); /* initial index into ehashtab*/ 
   if ((epnt=enhashtab[i])==NULL) {          /* install directly in table */
      enhashtab[i]  = nelem;
      enhashtabe[i] = nelem;
   }
   else {                               /* otherwise, go to end of list */
//      for (; epnt; epnt=epnt->hnnext)    /* find last element */
//          elast = epnt;    
//      elast->hnnext = nelem;            /* put new symbol at end of list */
//      nelem->hnlast = elast;        

      elast = enhashtabe[i];
      elast->hnnext = nelem;            /* put new symbol at end of list */
      nelem->hnlast = elast;           
      enhashtabe[i] = nelem;
   }
}

/*---------------------------------------------------*/

elem *findelem(int num)
                      
/* Find element among list of all elements.
   Elements are placed in hash table, which provides
   faster access than possible with one sequential list.
   The element itself is not moved, therefore only the 
   pointers of the element list need rearranging.
   Pointers in other lists that point to element remain correct.
*/

{
   register elem *epnt;
   int i,found;
   static int err=0;
        
  i = elemhash(num);
  for (found=0,epnt = ehashtab[i]; epnt; epnt = epnt->hnext) {
    if (epnt->elnum==num) {
       found = 1;
       break;
    }
  }

  if (found) return epnt;
  else {
     if (!err) {
       execerror("nc: can't find element.\n","");
       ncfprintf (stderr,"findelem: can't find elem %d\n",num);
       err = 1;
     }
     return (elem *)NULL; 
  }
}

/*---------------------------------------------------*/

elem *findelemne(int num)
                      
/* Find element among list of all elements.
   Elements are placed in hash table, which provides
   faster access than possible with one sequential list.
   The element itself is not moved, therefore only the 
   pointers of the element list need rearranging.
   Pointers in other lists that point to element remain correct.

   If error, don't print out anything, just return NULL.
*/

{
   elem *epnt;
   int i,found;
        
  i = elemhash(num);
  for (found=0,epnt = ehashtab[i]; epnt; epnt = epnt->hnext) {
    if (epnt->elnum==num) {
       found = 1;
       break;
    }
  }

  if (found) return epnt;
  else  return (elem *)NULL; 
}

/*---------------------------------------------------*/

elem *findelem (int etype, int na, int nb, int nc)

/* find element that connects to node */

{  
   node *npnt;
   elem *epnt;

  for (epnt=setforeachn2(etype,na,nb,&npnt);
                epnt=foreachn2(epnt); epnt = epnt->hnnext) {
	if (epnt->node1c != nc && epnt->node2c != nc) continue;
	break;
  }
  return epnt;
}

/*---------------------------------------------------*/

int eltsiz(int etype) 

{
   int esize;

  switch (etype) {

        default:      esize = 0;		break;
 	case ELEMENT: esize = sizeof(elem);	break;
 	case SPHERE:  esize = sizeof(sphere);	break;
 	case CABLE:   esize = sizeof(cable);	break;
 	case SYNAPSE: esize = sizeof(synapse);	break;
 	case ROD:     esize = sizeof(photorec);	break;
 	case CONE:    esize = sizeof(photorec);	break;
 	case CHR:     esize = sizeof(photorec);	break;
 	case VTRANSDUCER: esize = sizeof(photorec); break;
 	case ITRANSDUCER: esize = sizeof(photorec); break;
 	case ELECTRODE: esize = sizeof(electrode); break;
 	case GJ:      esize = sizeof(gapjunc);	break;
 	case PNX:     esize = sizeof(pannex);	break;
 	case RESISTOR:esize = sizeof(resistor);	break;
 	case DIODE:   esize = sizeof(diode);	break;
 	case LOAD:    esize = sizeof(loadelem);	break;
 	case GNDCAP:
 	case CAP:     esize = sizeof(capac);	break;
 	case BATT:    esize = sizeof(batt);	break;
 	case BUF:     esize = sizeof(vbuf);	break;
 	case NBUF:    esize = sizeof(nbuf);	break;
 	case CACOMP:  esize = sizeof(elem);	break;
 	case CHAN:    esize = sizeof(elem);	break;
  }
  return esize;
}

/*------------------------------------------------*/

int elemsiz(elem *epnt)

/* Return the size of a neural element in bytes. */

{
    int esize;

#ifdef DEBUG 
  if (debug & NCELEM && debugz & 16)  ncfprintf (stderr,"elemsiz %d\n",cumelem);
#endif

  if (!epnt) return(0);
  esize = eltsiz (epnt->ctype);
  return esize;
}

/*------------------------------------------*/

elem *elemcopy (elem *epnt, char *(*dalloc)(unsigned int))

/* Copy a neural element to a new location using dalloc to
    allocate space. */

{
    int i,siz;
    char *src,*dest;

  siz = elemsiz(epnt);
  src = (char*)epnt;
  if (!(dest=dalloc(siz))) {
    ncfprintf (stderr,"elemcopy, not enough space.\n");
    return (elem*)NULL;
  }
  for (i=0; i<siz; i++) {
    dest[i] = src[i];
  }
  return ((elem*)dest);
}

/*------------------------------------------*/

void elcopy (elem *src, elem *dest)

/* Copy a neural element's node information. */

{
  if (!src || !dest) return;

  dest->ctype  = src->ctype;
  dest->newelem  = src->newelem;
  dest->nocondens = src->nocondens;
  dest->node1a = src->node1a;
  dest->node1b = src->node1b;
  dest->node1c = src->node1c;
  dest->node1d = src->node1d;
  dest->node2a = src->node2a;
  dest->node2b = src->node2b;
  dest->node2c = src->node2c;
  dest->node2d = src->node2d;
  dest->attpnt = src->attpnt;
  dest->nodp1 = src->nodp1;
  dest->nodp2 = src->nodp2;
  dest->next = src->next;
  dest->last = src->last;
  dest->hnext = src->hnext;
  dest->hlast = src->hlast;
  dest->hnnext = src->hnnext;
  dest->hnlast = src->hnlast;
  dest->lptr  = src->lptr;
  dest->jnoise = src->jnoise;
  dest->rsd   = src->rsd;
  dest->modif = src->modif;
  dest->elnum = src->elnum;
  dest->region = src->region;
  dest->elabl = src->elabl;
}

/*------------------------------------------------*/

elem *tmpelem(void)

/* Set up an element to hold
   nodes and locations temporarily. */

{
  static elem nullelem = {0};
  static elem *epnt;

  epnt = &nullelem;
  epnt->ctype = ELEMENT;	/* element is default type */
  epnt->newelem = 1;
  epnt->nocondens = 0;
  epnt->elnum  = 0;
  epnt->region = NULLVAL;
  epnt->elabl  = (char *)"";
  epnt->next   = (elem *)NULL;
  epnt->last   = (elem *)NULL;
  epnt->hnext  = (elem *)NULL;
  epnt->hlast  = (elem *)NULL;
  epnt->hnnext = (elem *)NULL;
  epnt->hnlast = (elem *)NULL;
  epnt->node1a = NULLVAL;	/* reset node number */ 
  epnt->node1b = NULLVAL;
  epnt->node1c = NULLVAL;
  epnt->node1d = NULLVAL;
  epnt->node2a = NULLVAL;
  epnt->node2b = NULLVAL;
  epnt->node2c = NULLVAL;
  epnt->node2d = NULLVAL;
  epnt->nodp1  = (node *)NULL;
  epnt->nodp2  = (node *)NULL;
  epnt->attpnt = (attrib *)NULL;
  epnt->jnoise = 0;
  epnt->rsd    = 0;
  epnt->modif  = 0;
  epnt->lptr   = (elem *)NULL;
  elpnt = (elem *)NULL;		/* reset current element pointer */
  return (epnt);
}

/*------------------------------------------------*/

elem *makelem(int etype, elem *oepnt)

/* Make a new neural element, and link it
   to the elem list. Return a pointer to the 
   new element. Before doing this, though,
   check the old element (oepnt).  If it exists
   and has been already been initialized, then
   don't make a new element, but return the old
   one. */

{
    elem *epnt;
    int esize;

#ifdef DEBUG 
  if (debug & NCELEM && debugz & 16)  ncfprintf (stderr,"elem %d\n",cumelem);
#endif

  /* if modify, make null element with no defaults */

  if (oepnt && oepnt->ctype!=ELEMENT) {	/* If old element exists, */
        checkelemnode(oepnt);		/*  set its node pointers */
	return (oepnt);			/*  and return it. */ 
  }
  esize = eltsiz(etype);

#ifndef XSTIMM
    if ((epnt=(elem *)emalloc(esize)) == NULL) {
       ncfprintf (stderr,"no space left for elem %d\n",cumelem);
       return ((elem*)NULL);  
    }
    epnt->newelem = 1;

#else		/* if we are ignoring neural elements in "stim" */

    {  static char tmpelement[sizeof(synapse)*2] = {0};

       epnt = (elem *)&tmpelement;
    }
#endif

    if (oepnt) {
      elcopy (oepnt, (elem *)epnt);	/* first, copy node info. */
      checkelemnode(epnt);		/*  set elem's node pointers */
    }					/*  and the nodes' elem pointers */

    epnt->next = (elem *)NULL;		/* next, set pointers */
    if (!elempnt) elempnt = epnt;  	/* save head if first branch */
    if (elemend) {
      elemend->next = epnt;
      epnt->last = elemend;
    }
    else  {
      epnt->last = (elem *)NULL;
    }
    elemend = epnt;

    epnt->lptr  = (elem *)NULL;
    epnt->ctype  = etype;
    epnt->elnum = ++cumelem; 			/* increment total */
    einstall (epnt);				/* install new element in hash elnum table */

  switch (etype) {

        default:
 	case ELEMENT: break;

		/* default dqrm, dqri val = 1, no effect */
 case SPHERE:	
	if (interp) {
	  drm = getval("drm");
	  dcm = getval("dcm");
	  vcl = getval("vcl");
	}
	((sphere *)epnt)->dia    = NULLVAL;
	((sphere *)epnt)->Rm     = drm/exp(log(dqrm)*(tempcel-dbasetc)/10.0); 
	((sphere *)epnt)->Cm     = dcm;
	((sphere *)epnt)->vrev   = vcl;
	((sphere *)epnt)->vrest  = vcl;
	break;

 case CABLE:
	if (interp) {
	  dri = getval("dri");
	  drm = getval("drm");
	  dcm = getval("dcm");
	  vcl = getval("vcl");
	  complam = getval("complam");
	}
	((cable *)epnt)->dia    = NULLVAL;
	((cable *)epnt)->dia2   = NULLVAL;
	((cable *)epnt)->Ri     = dri/exp(log(dqri)*(tempcel-dbasetc)/10.0); 
	((cable *)epnt)->Rm     = drm/exp(log(dqrm)*(tempcel-dbasetc)/10.0); 
	((cable *)epnt)->Cm     = dcm;
	((cable *)epnt)->vrev   = vcl;
	((cable *)epnt)->vrest  = vcl;
	((cable *)epnt)->length = NULLVAL;
	((cable *)epnt)->cplam  = complam;
	break;

 case SYNAPSE:	
	if (interp) {
	  dst   = getval("dst");
	  dsms  = getval("dsms");
	  dsrrp = getval("dsrrp");
	  dstr  = getval("dstr");
	  dskd  = getval("dskd");
	  dshc  = getval("dshc");
	  dvg   = getval("dvg");
	  dsn   = getval("dsn");
	  dckd  = getval("dckd");
	  dchc  = getval("dchc");
	  dsg   = getval("dsg");
	}
	((synapse *)epnt)->thresh = dst;	/* synaptic theshold */
	((synapse *)epnt)->timec1 = (double *)NULL;
	((synapse *)epnt)->nfilt1 = NULLVAL;
	((synapse *)epnt)->timec1h = (double *)NULL;
	((synapse *)epnt)->nfilt1h = NULLVAL;
	((synapse *)epnt)->filt1hg = NULLVAL;
	((synapse *)epnt)->filt1ho = NULLVAL;
	((synapse *)epnt)->timec2 = (double *)NULL;
	((synapse *)epnt)->nfilt2 = NULLVAL;
	((synapse *)epnt)->timec3 = (double *)NULL;
	((synapse *)epnt)->nfilt3 = NULLVAL;
	((synapse *)epnt)->tfall2 = NULLVAL;
	((synapse *)epnt)->tfall3 = NULLVAL;
	((synapse *)epnt)->maxcond= NULLVAL;
	((synapse *)epnt)->rextern= NULLVAL;
	((synapse *)epnt)->maxsrate= NULLVAL;	/* maximum sustained rate */
	((synapse *)epnt)->mrrpool = NULLVAL;	/* max readily releaseable pool */
	((synapse *)epnt)->rrpool = NULLVAL;	/* readily releaseable pool */
	((synapse *)epnt)->rrpoolg = NULLVAL;	/* release gain from rrpool */
	((synapse *)epnt)->vrev   = NULLVAL;
	((synapse *)epnt)->vsize  = NULLVAL;
	((synapse *)epnt)->trconc = NULLVAL;	/* transmitter conc */
	((synapse *)epnt)->mesgconc = NULLVAL;	/* transmitter conc */
	((synapse *)epnt)->connum = 0;		/* postsynaptic connection number */
	((synapse *)epnt)->sens   = V;		/* voltage sensitivity */
	((synapse *)epnt)->ntact  = OPEN;	/* neurotransmitter action */
	((synapse *)epnt)->secmsg = 0;		/* second messenger */
	((synapse *)epnt)->nkd    = dskd;	/* half-max saturation point */
	((synapse *)epnt)->npow   = dshc;	/* hill coeff at postsyn rec */
	((synapse *)epnt)->vgain  = dvg;	/* neurotrans gain */
	((synapse *)epnt)->ngain  = dsn;	/* neurotrans gain */
	((synapse *)epnt)->caegain = NULLVAL;	/* ca expon ves release gain */
	((synapse *)epnt)->spost  = NULLVAL;
	((synapse *)epnt)->ckd    = dckd;	/* default cycG Kd */
	((synapse *)epnt)->chc    = dchc;	/* default cycG Hill coeff */
	((synapse *)epnt)->cgain  = dsg;	/* default cycG gain */
	((synapse *)epnt)->coff   = 1.0;	/* default cascade offset */
	((synapse *)epnt)->curve  = NULLVAL;
	((synapse *)epnt)->mesg1 = NULLVAL;
	((synapse *)epnt)->mesg2 = NULLVAL;
	((synapse *)epnt)->dyadelem = 0;
	break;

 case VTRANSDUCER:
 case ITRANSDUCER:
 case CONE:
 case CHR:
 case ROD:((photorec *)epnt)->xpos   = NULLVAL;
	 ((photorec *)epnt)->ypos   = NULLVAL;
	 ((photorec *)epnt)->zpos   = NULLVAL;
	 ((photorec *)epnt)->dia    = 1.674;	/* 2.2um2 area */
	 ((photorec *)epnt)->maxcond= NULLVAL;
	 ((photorec *)epnt)->pigm   = NULLVAL; 
	 ((photorec *)epnt)->pathl  = NULLVAL;
	 ((photorec *)epnt)->attf   = 1.0;
	 ((photorec *)epnt)->timec1 = NULLVAL;
	 ((photorec *)epnt)->loopg  = NULLVAL;
	 ((photorec *)epnt)->filt   = NULLVAL;
	 ((photorec *)epnt)->qc     = NULLVAL;
	 ((photorec *)epnt)->qcond  = NULLVAL;
	 ((photorec *)epnt)->linit  = NULLVAL;
	 ((photorec *)epnt)->unitary = NULLVAL;
	 ((photorec *)epnt)->save   = 0;
	 ((photorec *)epnt)->restore= 0;
	 ((photorec *)epnt)->darknoise=NULLVAL;
	 ((photorec *)epnt)->channoise=NULLVAL;
	 ((photorec *)epnt)->photnoise=NULLVAL;
	 ((photorec *)epnt)->stimchan=NULLVAL;
	 ((photorec *)epnt)->chrseed=0;
	 ((photorec *)epnt)->dkrseed=0;
	 ((photorec *)epnt)->phrseed=0;
	 break;

 case ELECTRODE:((electrode *)epnt)->dia    = NULLVAL;
	 ((electrode *)epnt)->length = NULLVAL;
	 ((electrode *)epnt)->r = NULLVAL;
	 ((electrode *)epnt)->c = NULLVAL;
	 ((electrode *)epnt)->ccomp = NULLVAL;
	 ((electrode *)epnt)->vrest = NULLVAL;
	 ((electrode *)epnt)->vrev = NULLVAL;
	 break;

 case PNX:
	 ((pannex *)epnt)->atp_gain= NULLVAL;
	 ((pannex *)epnt)->atp     = NULLVAL;
	 ((pannex *)epnt)->adp     = NULLVAL;
	 ((pannex *)epnt)->amp     = NULLVAL;
	 ((pannex *)epnt)->h       = NULLVAL;
	 ((pannex *)epnt)->atp_decr= NULLVAL;
	 ((pannex *)epnt)->pnxfunc = NULL;
 case GJ:((gapjunc *)epnt)->area   = NULLVAL;
	 ((gapjunc *)epnt)->specres= NULLVAL;
	 ((gapjunc *)epnt)->vgain  = NULLVAL;
	 ((gapjunc *)epnt)->voff   = NULLVAL;
	 ((gapjunc *)epnt)->taun   = NULLVAL;
	 ((gapjunc *)epnt)->rvgain = NULLVAL;
	 ((gapjunc *)epnt)->rvoff  = NULLVAL;
	 ((gapjunc *)epnt)->rtaun  = NULLVAL;
	 ((gapjunc *)epnt)->gmax   = NULLVAL;
	 ((gapjunc *)epnt)->gnv    = NULLVAL;
	 ((gapjunc *)epnt)->nodeca = NULLVAL;
	 ((gapjunc *)epnt)->nodecb = NULLVAL;
	 ((gapjunc *)epnt)->nodecc = NULLVAL;
	 ((gapjunc *)epnt)->nodecd = NULLVAL;
	 ((gapjunc *)epnt)->rect   = 0;
	 ((gapjunc *)epnt)->modtyp = 0;
	 ((gapjunc *)epnt)->sign   = 0;
	 break;


 case RESISTOR:	((resistor *)epnt)->r = NULLVAL;
	 	break;

 case DIODE:	((diode *)epnt)->r = NULLVAL;
	 	break;

 case LOAD:((loadelem *)epnt)->r = NULLVAL;
	 ((loadelem *)epnt)->vrev  = 0;
	 ((loadelem *)epnt)->vrest = 0;
	 break;

 case GNDCAP:
 case CAP:	((capac *)epnt)->c  = NULLVAL;
		((capac *)epnt)->vrest = NULLVAL;
	 	break;

 case BATT:	((batt *)epnt)->v  = NULLVAL;
	 	break;

 case BUF:	((vbuf *)epnt)->delay  = NULLVAL;
 		((vbuf *)epnt)->offset = NULLVAL;
 		((vbuf *)epnt)->gain   = NULLVAL;
 		((vbuf *)epnt)->tau    = NULLVAL;
 		((vbuf *)epnt)->lphp   = NULLVAL;
	 	break;
 case NBUF:	((nbuf *)epnt)->ntoffset = NULLVAL;
 		((nbuf *)epnt)->offset = NULLVAL;
 		((nbuf *)epnt)->gain   = NULLVAL;
 		((nbuf *)epnt)->ntrans  = NULLVAL;
	 	break;

  }	/* switch (etype)  to initialize the elements */

  return (epnt); 
}

/*---------------------------------------------------*/

void delelem (elem *epnt)

/* Delete an element, including its attributes and 
    its hash table entry, and list pointers but not
    any node pointers to it.
*/

{
    int i;
    attrib *apnt, *tapnt;
    elem *last,*next;

    for (apnt=epnt->attpnt; apnt; ) {   /* free element's attributes first */
        tapnt = apnt;
        switch (apnt->ctype) {
	  case CCHNOISE:
	  case VESNOISE: cumnattr--; break;
	  case CGMP:
	  case NMDA:
	  case CA:
	  case CACOMP:   cumcattr--; break;
              default:   cumattr--;  break;
        }
        apnt = apnt->attpnt;
        efree (tapnt);
    }
    next = epnt->next;
    last = epnt->last;
    if (last) last->next = next;
    if (next) next->last = last;
    if (elempnt==epnt)			/* is elempnt the one being erased? */
        elempnt = next;			/*  yes, then incr elempnt */
    if (elemend==epnt)			/* is elemend the one being erased? */
        elemend = last;			/*  yes, then decr elemend */

    i = elemhash(epnt->elnum);					/* get elem hash index */
    if (epnt->hnext!=NULL) epnt->hnext->hlast = epnt->hlast;	/* delete hash last pointer */
    if (epnt->hlast!=NULL) epnt->hlast->hnext = epnt->hnext;	/* delete hash next pointer */
    if (ehashtab[i]==epnt)  ehashtab[i] = epnt->hnext;
    if (ehashtabe[i]==epnt) ehashtabe[i] = epnt->hlast;

    i = nelemhash(epnt->node1a,epnt->node1b); 			 /* get elem node hash index */ 
    if (epnt->hnnext!=NULL) epnt->hnnext->hnlast = epnt->hnlast; /* delete node hash last pointer */
    if (epnt->hnlast!=NULL) epnt->hnlast->hnnext = epnt->hnnext; /* delete node hash next pointer */
    if (enhashtab[i]==epnt)  enhashtab[i] = epnt->hnnext;
    if (enhashtabe[i]==epnt) enhashtabe[i] = epnt->hnlast;

    switch (epnt->ctype) {
	double *filt;

     case SYNAPSE: if (filt=((synapse*)epnt)->timec1) efree (filt);
     		   if (filt=((synapse*)epnt)->timec1h)efree (filt);
     		   if (filt=((synapse*)epnt)->timec2) efree (filt);
     		   if (filt=((synapse*)epnt)->timec3) efree (filt);
		  break;
	default:	break;
    }
    efree (epnt);
    delelmt++;
/*    cumelem--;  */			/* don't decrement elem number */
}

/*----------------------------------------*/

attrib *getattrib(int elnum)

/* Return pointer to chan attribute of element, given its element
 * number */

{
   elem *elpnt;
   attrib *apnt;

if (!(elpnt=findelem(elnum))) {
      ncfprintf (stderr,"getattrib: can't find element %d\n",elnum);
      execerror ("Missing element: "," stopping... ");
      return NULL;
 }

  if ((apnt=elpnt->attpnt)==NULL) {
ncfprintf (stderr,"getattrib: can't find channel attrib for element %d\n",elnum);
       return NULL;
  }
 return apnt;
}

/*----------------------------------------*/


elem *get_elempnt(int elnum)

/* get pointer to element, given its element number */

{
    elem *epnt;

  if (!(epnt=findelem(elnum))) {
  	ncfprintf (stderr,"get_elempnt: can't find element %d\n",elnum);
  	return NULL;
  }
  else checkelemnode(epnt);
  return epnt;
}

/*---------------------------------------------------------------------*/

elem *get_elempnt(node *npnt, int elnum)

/* get pointer to element relative to node */

{
   int i;
   elem *epnt=(elem*)NULL;
   conlst *lpnt=(conlst*)NULL;

   if (!npnt) return NULL;
   for (lpnt=npnt->elemlst,i=0; lpnt && i<elnum; i++) /* find the element */
     lpnt=lpnt->next;
   if (lpnt) epnt=(elem *)lpnt->conpnt;
   if (!lpnt || !epnt) {
     ncfprintf (stderr,"get_elempnt: can't find elem %d relative to node %s\n",elnum,prnode(npnt));
     execerror ("Missing element; ","stopping..."); /* */
     return NULL;
   }
   else return (epnt);
}

/*---------------------------------------------------------------------*/

int elemtype(int elnum)

{
   elem *epnt;

   if (!(epnt=findelem(elnum))) {
      ncfprintf (stderr,"Edist: can't find element %d\n",elnum);
      return 0;
   }
   else checkelemnode(epnt);
   return epnt->ctype;
}



