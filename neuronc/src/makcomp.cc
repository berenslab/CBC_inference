/* makcomp.cc */

#include "nc.h"
#include "y.tab.h"
#include "ncsub.h"
#include "ncelem.h"
#include "ncomp.h"
#include "control.h"

// #define DEBUG

#ifdef DEBUG
#include "ncdebug.h"
#endif

#define CHASHSIZ 9949                   /* prime number to make better hash */

static comp *chashtab[CHASHSIZ] = {0};  /* comp hash table */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

#ifdef __cplusplus
}
#endif

#include "ncio.h"

extern comp *compnt;
extern comp *compend;
extern int cumcomp;

char *emalloc(unsigned int n);
void freelst (conlst *clst);
void delntcomp(ntcomp* ntpnt);
int noduniq(nodeint val1, nodeint val2, nodeint val3, nodeint val4, int seed);
char *makrand (int rndsiz, const char *mesg, int num);
void initstate (unsigned long int seed, char *state, int n);
void efree(void *ptr);
void delcacomp(cacomp *capnt);

/*---------------------------------------------------*/

/*
void xdebugf(void)

{
   ncfprintf (stderr,"chashtab %d\n",chashtab[5]);
}
*/

/*---------------------------------------------------*/

void cinithash(void)
{
   int i;

   for (i=0; i<CHASHSIZ; i++) 
     chashtab[i] = 0;
}

/*---------------------------------------------------*/

int comphash(int num)

/* Make disordered index of comp number */

{
   return (num % CHASHSIZ); 
}

/*---------------------------------------------------*/

void cinstall (comp *newcomp, int num)
              
/* Install a compartment in hash table */
/*  Uses "direct chaining" with field "hnext". */

{
    int i;
    comp *cpnt,*clast;
    static int cinitfl=0;

   if (!cinitfl) {			/* initialize table once at start */
        cinitfl = 1;
        cinithash();
   }
   newcomp->hnext = (comp*)NULL;
   i=comphash(num); 			/* initial index into chashtab*/ 
   if ((cpnt=chashtab[i])==NULL) {      /* install directly in table */
      chashtab[i] = newcomp;
   }
   else {                               /* otherwise, go to end of list */
      for (; cpnt; cpnt=cpnt->hnext)    /* find last symbol */
          clast = cpnt;    
      clast->hnext = newcomp;            /* put new comp at end of list */
   }
}

/*---------------------------------------------------*/

comp *makcomp(elem *epnt, double rm, double cap, double vrest, double vrev)

/* make a new compartment and link it to the compartment list. */

{
    comp *cpnt;
    int ctype;
 
#ifdef DEBUG
 if (debug & NCCONVERT && debugz & 16)
	 ncfprintf (stderr,"comp %d in branch %d\n",cumcomp+1,epnt->elnum);
#endif

  if ((cpnt=(comp *)emalloc(sizeof(comp))) == NULL) {
     ncfprintf (stderr,"no space left for comp %d in branch %d\n",
				cumcomp+1,epnt->elnum);
     return ((comp*)NULL);  
  }
  cpnt->next = (comp*)NULL;		/* make new comp be end of list */
  if (!compnt) compnt = cpnt;  		/* save head if first synap */
  if (compend)
    compend->next = cpnt;
  cpnt->last = compend;			/* pointer to last compartment */
  compend = cpnt;
  cpnt->nodlst = (conlst*)NULL;		/* pointer to list of nodes for comp */
  cpnt->pvext = (comp*)NULL;		/* pointer to possible external voltage comp */
  cpnt->capnt = (cacomp*)NULL;		/* pointer to possible ca comp */
  cpnt->ntpnt = (ntcomp*)NULL;		/* pointer to possible cycg comp */
  cpnt->num 	= ++cumcomp;
  cpnt->domain = 0;

  cinstall(cpnt,cpnt->num);		/* install in hash table */

  ctype = epnt->ctype;
  cpnt->ctype 	= ctype;
  if (rm<=0)   rm = 1e30;		/* check for zero size compartment */
  if (cap<=0) cap = 1e-30;
  cpnt->rm 	= rm;
  cpnt->cap 	= cap;
  cpnt->miscfl 	= 0;			/* clear misc flags (IEXT,VEXT,VBAT,CA, etc. */
  cpnt->virt    = 0;
  cpnt->extv 	= 0.0;
  cpnt->extv_old = 0.0;
  cpnt->exti 	= 0.0;
  cpnt->exti_old = 0.0;
  cpnt->v       = 0.0;
  cpnt->vph     = 0.0;
  cpnt->vext    = 0.0;
  cpnt->oldvext = 0.0;
  cpnt->oldv    = 0.0;
  cpnt->vest    = 0.0;
  cpnt->implf   = 0.0;
  cpnt->nodc    = 0.0;
  cpnt->nodcur  = 0.0;
  cpnt->tcond   = 0.0;
  cpnt->tcondt  = 0.0;
  cpnt->vrev    = 0.0;

  if (epnt->jnoise) {
	int nrseed;
     cpnt->jnoise = epnt->jnoise;
     cpnt->jstate = (char*)NULL;
     nrseed = epnt->rsd;
     if (nrseed <= 0)
	 nrseed = noduniq(15791, cpnt->num, 27211, cpnt->num,rseed^3333);
     if (nrseed > 0) {
        cpnt->jstate = makrand (RNDSIZ,"Johnson",cpnt->num);
        if (cpnt->jstate) initstate (nrseed,cpnt->jstate,RNDSIZ);
     }
  }
  switch (ctype) {		/* neural elements that control vrest, vrev */
   				/*   if they are alone in compartment */
   case CABLE:
   case SPHERE:
   case LOAD: 
   case NA: 
   case CGMP: 
   case CA: 
   case K: 
	if (vrest == NULLVAL) cpnt->v = vcl;
	else	            cpnt->v = vrest;
	cpnt->vest = cpnt->v;
	if (vrev  == NULLVAL) cpnt->vrev = vcl;
	else		    cpnt->vrev = vrev;
	break;

   default:			/* other elements can make comps, but don't */
	break;			/*    control vrest, vrev. */
  }

  cpnt->clst 	= (conlst*)NULL;
  cpnt->capnt 	= (cacomp*)NULL;
  return (cpnt); 
}

/*------------------------------------*/

void delhashcomp(comp *cpnt)

/* Delete a compartment from the hash table. */

{
    int i, found;
    comp *clast,*cpt;

    if (!cpnt) return;

	/* First, delete from hash table: */

  i = comphash(cpnt->num);
  clast = (comp*)NULL;
  for (found=0, cpt=chashtab[i]; cpt; clast=cpt, cpt=cpt->hnext) {
    if (cpt==cpnt) {
        found = 1;
        break;
    }
  }
  if (found) {
     if (clast) clast->hnext = cpt->hnext;     /* delete hash pointer */
     else chashtab[i] = cpt->hnext;	       /*  move pointer into table */
  }
}

/*---------------------------------------------------*/

void delcomp (comp *cpnt)

/* Delete a compartment. Patch pointers to the compartment list
so that the list is not broken.  Delete the comp's node list
(i.e. the list of nodes that point to the comp) by deleting each
node pointer. Don't delete node's comp pointer because the node
is already deleted.  */
{

    freelst (cpnt->clst);		/* free the connection list */
    freelst (cpnt->nodlst);		/* free the node list */
    delntcomp(cpnt->ntpnt);		/* free nt comps */
    delcacomp(cpnt->capnt);		/* free ca comp */
    if (cpnt->last) cpnt->last->next = cpnt->next;/* patch list pointers */
    else compnt = cpnt->next;
    if (cpnt->next) cpnt->next->last = cpnt->last;
    else compend = cpnt->last;
    delhashcomp(cpnt);
    efree (cpnt);			/* free the compartment */
    cumcomp--;
}

/*---------------------------------------------------*/

comp *findcomp(int num, const char *str)
                      
/* Find comp among list of all nodes.
   Comps are placed in hash table, which provides
   faster access than possible with one sequential list.
   Depends on "hnext" field in comp.
*/

{
   comp *cpnt;
   int i,found;
        
  i = comphash(num);
  for (found=0,cpnt = chashtab[i]; cpnt; cpnt = cpnt->hnext) {
    if (cpnt->num==num) {
       found = 1;
       break;
    }
  }

  if (found) return cpnt;
  else {
    if (str) ncfprintf (stderr,"%s: can't find comp %d\n",str,num);
    return (comp*)NULL; 
  }
}

/*---------------------------------------------------*/

