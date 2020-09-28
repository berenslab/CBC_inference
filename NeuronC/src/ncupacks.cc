/* Segment ncupacks in Program nc */

/* Unpacks nc and re-creates structures from parallel 
   message-passing routines.  See "ncpacks.cc".
*/

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "ncomp.h"
#include "ncelem.h"
#include "ncvirt.h"
#include "nclist.h"
#include "pvm3.h"
#include "ncio.h"

extern double timinc;

extern double pigmlen[];         /* path length through o.s. (in ncstim) */

extern recpar *rectypes[];

chan *alloc_chan(short ctype, short stype); 
void linkchan(chan *cpnt);
void maktables (double timinc);
char *makrand (int rndsiz, char *mesg, int num);
char *emalloc (int n);
void initrec();
comp *findcomp(int num, const char *str);

upk_byte  (char *val,   int n, int stride);
upk_short (short *val,  int n, int stride);
upk_long  (long  *val,  int n, int stride);
upk_int   (int *val,    int n, int stride);
upk_float (float *val,  int n, int stride);
upk_double(double *val, int n, int stride);
upk_str   (char *val);

/*------------------------------------------------*/

void upk_stconc (stconc *cpnt)

/* Unpack a "stconc" into the send buffer. */

{
  upk_double(&cpnt->cval,1,1);
  upk_double(&cpnt->cest,1,1);
  upk_double(&cpnt->dcon,1,1);
}

/*------------------------------------------------*/

chan *upk_chan ()

/* Unpack and create a "chan" from the receive buffer. */

{
   int i, one=1, zero=0;
   char temp;
   short ctype, stype;
   chan *cpnt;
   static maktabfl=0;

  upk_short(&ctype,1,1);		/* channel type */
  upk_short(&stype,1,1);		/* channel sub-type */
  cpnt = alloc_chan(ctype,stype); 
  linkchan(cpnt);

  if (!maktabfl) {
        maktabfl = 1;
        maktables(timinc);
  }
 
  cpnt->ctype = ctype;
  cpnt->stype = stype;
  upk_int(&cpnt->num,1,1);
  upk_double(&cpnt->conduct,1,1);
  upk_int(&cpnt->num1,1,1);
  upk_int(&cpnt->num2,1,1);

  upk_double(&cpnt->maxcond,1,1);
  upk_double(&cpnt->voffsm,1,1);
  upk_double(&cpnt->voffsh,1,1);
  upk_double(&cpnt->vrev,1,1);
  if (!(cpnt->chtyp=findchantype(ctype,stype)) {
    ncfprintf (stderr,"ncupacks: can't find type %d %d\n",ctype,stype);
  }
  upk_double(&cpnt->cdur,1,1);
  upk_byte(&temp,1,1);		/* 1 => cstate */
  if (temp) {
     cpnt->cstate = makrand(RNDSIZ,"chan",cumchan+1);
     upk_byte(cpnt->cstate,RNDSIZ,1);
  }
  else cpnt->cstate  = NULL;
  upk_short(&cpnt->nchan,1,1);
  upk_double(&cpnt->taua,1,1);
  upk_double(&cpnt->taub,1,1);
  switch (cpnt->ctype) {
    case NA:
      switch (cpnt->stype) {
       case 0:
         upk_double(&((hhchan *)cpnt)->m,1,1);
         upk_double(&((hhchan *)cpnt)->h,1,1);
	 break;
       case 1:
         upk_short(&((chan *)cpnt)->numstate,1,1);
         for (i=0; i<((chan *)cpnt)->numstate; i++) {
            upk_stconc(&((chan *)cpnt)->conc[i]);
         }
	 break;
       default:
         upk_double(&((hhchan *)cpnt)->m,1,1);
         upk_double(&((hhchan *)cpnt)->h,1,1);
	 break;
      }
      break;
    case K:
      switch (cpnt->stype) {
       case 0:
       case 2:
         upk_double(&((hhchan *)cpnt)->m,1,1);
         upk_double(&((hhchan *)cpnt)->h,1,1);
	 break;
       case 1:
       case 3:
       case 4:
       case 5:
       case 6:
       case 7:
       case 8:
       case 9:
      default:
         upk_short(&((chan *)cpnt)->numstate,1,1);
         for (i=0; i<((chan *)cpnt)->numstate; i++) {
            upk_stconc(&((chan *)cpnt)->conc[i]);
         }
	 break;
      }
      break;

    case ClCa:
         upk_short(&((chan *)cpnt)->numstate,1,1);
         for (i=0; i<((chan *)cpnt)->numstate; i++) {
            upk_stconc(&((chan *)cpnt)->conc[i]);
         }
	 break;

    case KCa:
      switch (cpnt->stype) {
       case 0:
       case 2:
         upk_double(&((kcachan *)cpnt)->m,1,1);
         upk_double(&((kcachan *)cpnt)->h,1,1);
         upk_double(&((kcachan *)cpnt)->d1,1,1);
         upk_double(&((kcachan *)cpnt)->d2,1,1);
         upk_double(&((kcachan *)cpnt)->k1,1,1);
         upk_double(&((kcachan *)cpnt)->k2,1,1);
	 break;
       case 1:
       case 3:
       case 4:
       case 5:
       case 6:
      default:
         upk_short(&((chan *)cpnt)->numstate,1,1);
         for (i=0; i<((chan *)cpnt)->numstate; i++) {
            upk_stconc(&((chan *)cpnt)->conc[i]);
         }
	 break;
      }
      break;
    case CA:
      switch (cpnt->stype) {
       case 0:
         upk_double(&((hhchan *)cpnt)->m,1,1);
         upk_double(&((hhchan *)cpnt)->h,1,1);
	 break;
       default:
	 break;
       }
      break;
  }
  upk_byte(&temp,1,1);			/* end of structure */
  if (temp!=VEND) {
     ncfprintf (stderr,"upk_chan: invalid format for chan ('%d')\n",temp);
     return NULL;
  }
 cumchan++;
 return cpnt;
}

/*------------------------------------------------*/

void upk_cn (conn *cpnt)

/* Unpack the conn base into the send buffer. */

{
  upk_short(&cpnt->ctype,1,1);
  upk_double(&cpnt->conduct,1,1);
  upk_int(&cpnt->num1,1,1);
  upk_int(&cpnt->num2,1,1);
}

/*------------------------------------------------*/

void upk_lpfilt (lpfilt *cpnt)

/* Unpack an lpfilt from the receive buffer. */
/*  Space has already been defined */

{
  upk_short(&cpnt->nfilt,1,1);
  upk_double(&cpnt->ftau,1,1);
  upk_double(cpnt->lfilt,NUMFILT+1,1);
  upk_double(&cpnt->tfall,1,1);
}

/*------------------------------------------------*/

void upk_ld (load *cpnt)

/* Unpack a "load" base into the send buffer. */

{
  upk_short(&cpnt->ctype,1,1);
  upk_double(&cpnt->conduct,1,1);
  upk_int(&cpnt->num1,1,1);
  upk_double(&cpnt->vrev,1,1);
}

/*------------------------------------------------*/

load *upk_load ()

/* Unpack and create a "load" from the receive buffer. */

{
   char temp;
   load *lpnt;

  if ((lpnt=(load *)emalloc(sizeof(load))) == NULL) {
    ncfprintf (stderr,"no space left for load %d\n",cumload+1);
    return (NULL);
  }
  lpnt->next = NULL;
  if (!loadpnt) loadpnt = lpnt;       /* save head if first load */
  if (loadend)
    loadend->next = lpnt;
  loadend = lpnt;
  cumload++;

  upk_ld(lpnt); 
  upk_byte(&temp,1,1);
  if (temp!=VEND) {
     ncfprintf (stderr,"upk_load: invalid format for load ('%d')\n",temp);
     return NULL;
  }
  return (lpnt);
}

/*------------------------------------------------*/

dbuf *upk_dbuf ()

/* Unpack and create a "dbuf" from the receive buffer. */

{
   dbuf *cpnt;

   int idelay;
   char temp;
   short delp;
   short *dpnt;

  if ((cpnt=(dbuf *)emalloc(sizeof(dbuf))) == NULL) {
    ncfprintf (stderr,"no space left for dbuf %d\n",cumconn+1);
    return NULL;
  }
  cpnt->next = NULL;
  if (!synpnt) synpnt = (synap *)cpnt; /* save head if first synap */
  if (synend)
     synend->next = (synap *)cpnt;
  synend = (synap *)cpnt;
  cumconn++;
 
  upk_cn(cpnt);			/* unpack the base conn */
  upk_double(&cpnt->delay,1,1);
  idelay = (int)cpnt->delay;
  if (idelay)  {
    if ((dpnt=(short int *)emalloc(idelay * sizeof(short int))) == NULL) {
      ncfprintf (stderr,"no space left for delay buffer %d\n",cumconn+1);
      return NULL;
     }
   }
   else dpnt = NULL;

  if (dpnt) {
    upk_short(dpnt,idelay,1);
    cpnt->delbuf = dpnt;
    upk_short(&delp,1,1);
    cpnt->delpnt = &dpnt[delp];
  }
  else {
    cpnt->delbuf = NULL;
    cpnt->delpnt = NULL;
  }
  upk_byte(&temp,1,1);
  if (temp!=VEND) {
     ncfprintf (stderr,"upk_dbuf: invalid format for dbuf ('%d')\n",temp);
     return NULL;
  }
  return cpnt;
}

/*------------------------------------------------*/

recep *upk_recep ()

/* Unpack a synap from the receive buffer. */

{
   char temp;
   short pigm;
   static int rectabfl=0;
   recep *cpnt;

  if ((cpnt=(recep *)emalloc(sizeof(recep))) == NULL) {
    ncfprintf (stderr,"no space left for recep %d\n",cumphotrec+1);
    return (NULL);
  }
  cpnt->next = NULL;
  if (!recpnt) recpnt = (recep *)cpnt; /* save head if first recep */
  if (recend)
     recend->next = (recep *)cpnt;
  recend = (recep *)cpnt;
  cumphotrec++;

  upk_ld(cpnt);
  upk_short(&pigm,1,1);
  
  if (!rectabfl) {
        rectabfl = 1;
        initrec();
  }
 
  cpnt->chtyp = rectypes[pigm];
  upk_double(&cpnt->xloc,1,1);
  upk_double(&cpnt->yloc,1,1);
  upk_double(&cpnt->pathl,1,1);
  upk_double(&cpnt->attf,1,1);
  upk_double(&cpnt->timec,1,1);
  upk_double(&cpnt->dnois,1,1);
  upk_short(&cpnt->filt,1,1);
  upk_short(&cpnt->pnois,1,1);
  upk_byte(&temp,1,1);		/* 1 => dstate */
  if (temp) {
     cpnt->dstate = makrand(RNDSIZ,"recep",cumphotrec+1);
     upk_byte(cpnt->dstate,RNDSIZ,1);
  }
  else cpnt->dstate  = NULL;
  upk_byte(&temp,1,1);		/* 1 => cstate */
  if (temp) {
     cpnt->pstate = makrand(RNDSIZ,"recep",cumphotrec+1);
     upk_byte(cpnt->pstate,RNDSIZ,1);
  }
  else cpnt->pstate  = NULL;
  upk_double(&cpnt->area,1,1);
  upk_double(&cpnt->aflux,1,1);
  upk_double(&cpnt->mflux,1,1);
  upk_double(&cpnt->iflux,1,1);
  upk_double(&cpnt->rhod,1,1);
  upk_double(&cpnt->gpr1,1,1);
  upk_double(&cpnt->gpr2,1,1);
  upk_double(&cpnt->gpr3,1,1);
  upk_double(&cpnt->pde,1,1);
  upk_double(&cpnt->gcyc,1,1);
  upk_double(&cpnt->cycg,1,1);
  upk_double(&cpnt->cond,1,1);
  upk_double(&cpnt->ca,1,1);
  upk_double(&cpnt->cax,1,1);
  upk_double(&cpnt->maxcond,1,1);
  upk_short(&cpnt->recnm1,1,1);
  upk_short(&cpnt->recnm2,1,1);
  upk_short(&cpnt->recnm3,1,1);
  upk_byte(&temp,1,1);
  if (temp!=VEND) {
     ncfprintf (stderr,"upk_recep: invalid format for receptor ('%d')\n",temp);
     return NULL;
  }
  return cpnt;
}

/*------------------------------------------------*/

synap *upk_synap()

/* Unpack a synap from the receive buffer. */

{
   char temp;
   synap *cpnt;

  if ((cpnt=(synap *)emalloc(sizeof(synap))) == NULL) {
    ncfprintf (stderr,"no space left for synap %d\n",cumsynap+1);
    return (NULL);
  }
  cpnt->next = NULL;
  if (!synpnt) synpnt = (synap *)cpnt; /* save head if first synap */
  if (synend)
     synend->next = (synap *)cpnt;
  synend = (synap *)cpnt;
  cumsynap++;

  upk_cn(cpnt);
  upk_double(&cpnt->vrev,1,1);
  upk_double(&cpnt->thresh,1,1);
  upk_double(&cpnt->nkd,1,1);
  upk_double(&cpnt->ngain,1,1);
  upk_double(&cpnt->ckd,1,1);
  upk_double(&cpnt->cgain,1,1);
  upk_double(&cpnt->maxcond,1,1);
  upk_double(&cpnt->cdur,1,1);
  upk_double(&cpnt->vgain,1,1);
  upk_double(&cpnt->vsize,1,1);
  upk_short(&cpnt->vsites,1,1);
  upk_short(&cpnt->nchan,1,1);
  upk_short(&cpnt->ntact,1,1);
  upk_short(&cpnt->curve,1,1);
  upk_int(&cpnt->sdyadnum,1,1);
  upk_byte(&temp,1,1);		/* 1 => cstate */
  if (temp) {
     cpnt->cstate = makrand(RNDSIZ,"recep",cumphotrec+1);
     upk_byte(cpnt->cstate,RNDSIZ,1);
  }
  else cpnt->cstate  = NULL;
  upk_byte(&temp,1,1);		/* 1 => vstate */
  if (temp) {
     cpnt->vstate = makrand(RNDSIZ,"recep",cumphotrec+1);
     upk_byte(cpnt->vstate,RNDSIZ,1);
  }
  else cpnt->vstate  = NULL;

  upk_byte(&temp,1,1);		/* 1 => filter */
  if (temp) {
       lpfilt *fpnt;

      if ((fpnt=(lpfilt *)emalloc(sizeof(lpfilt))) == NULL) {
        ncfprintf (stderr,"no space left for filter 1 %d\n",cumsynap+1);
        return (NULL);
      }
    cpnt->filt1  = fpnt;              /* filter 1 */
    upk_lpfilt(fpnt);
  }
  else cpnt->filt1  = NULL;

  upk_byte(&temp,1,1);		/* 1 => filter */
  if (temp) {
       lpfilt *fpnt;

      if ((fpnt=(lpfilt *)emalloc(sizeof(lpfilt))) == NULL) {
        ncfprintf (stderr,"no space left for filter 2 %d\n",cumsynap+1);
        return (NULL);
      }
    cpnt->filt2  = fpnt;              /* filter 2 */
    upk_lpfilt(fpnt);
  }
  else cpnt->filt2  = NULL;

  upk_byte(&temp,1,1);		/* 1 => filter */
  if (temp) {
       lpfilt *fpnt;

      if ((fpnt=(lpfilt *)emalloc(sizeof(lpfilt))) == NULL) {
        ncfprintf (stderr,"no space left for filter 3 %d\n",cumsynap+1);
        return (NULL);
      }
    cpnt->filt3  = fpnt;              /* filter 3 */
    upk_lpfilt(fpnt);
  }
  else cpnt->filt3  = NULL;

  upk_byte(&temp,1,1);
  if (temp!=VEND) {
     ncfprintf (stderr,"upk_synap: invalid format for synapse ('%d')\n",temp);
     return NULL;
  }
  return cpnt;
}

/*------------------------------------------------*/

conn *upk_conn ()

/* Unpack a conn from the receive buffer. */

{
   conn *cpnt;
   char temp;
    
  if ((cpnt=(conn *)emalloc(sizeof(conn))) == NULL) {
    ncfprintf (stderr,"no space left for conn %d\n",cumconn+1);
     return (NULL);
   }
  cumconn++;
  upk_cn(cpnt);
  upk_byte(&temp,1,1);			/* end of structure */
  if (temp!=VEND) {
     ncfprintf (stderr,"upk_chan: invalid format for conn ('%d')\n",temp);
     return NULL;
  }
  return cpnt;
}

/*------------------------------------------------*/

conlst *upk_conlst ()

/* Unpack a list of connections from the receive buffer. */

{
   conlst *head,*lpnt,*olpnt;
   char temp;
   conn *cpnt;
   int endlist;

 for (head=olpnt=NULL,endlist=0; !endlist; ) {
    upk_byte(&temp,1,1);		/* first byte of struct on conlst */
    switch (temp) {

       case VCONN:
		  cpnt = upk_conn();
		  break;

       case VCHAN:
		  cpnt = (conn *)upk_chan();
		  break;
       case VSYNAP:
		  cpnt = (conn *)upk_synap();
		  break;
       case VRECEP:
		  cpnt = (conn *)upk_recep();
		  break;
       case VLOAD:
		  cpnt = (conn *)upk_load();
		  break;
       case VBUF:
		  cpnt = (conn *)upk_dbuf();
		  break;

       case VEND: endlist = 1; break;

       default:	  ncfprintf (stderr,"upk_conlst: bad format\n");
		  endlist = 1;
		  break;
    }
   if (!cpnt) {
	ncfprintf (stderr,"upk_conlst: unable to make list\n");
	return NULL;
   }
   if ((lpnt=(conlst *)emalloc(sizeof(conlst))) == NULL) {
	ncfprintf (stderr,"upk_conlst: no space left for conlist\n");
	return NULL;
   }
   if (! head) head = lpnt;		/* set head of list */
   if (olpnt) olpnt->next = lpnt;
   lpnt->conpnt = cpnt; 
   lpnt->next = NULL;
   olpnt = lpnt;
  }
  return head;
}

/*------------------------------------------------*/

cacomp *upk_cacomp (comp *comp1)

/* Unpack a calcium compartment from the receive buffer. */
/* Must already have pointer to voltage compartment. */

{
    conlst *lpnt;
    char temp;
    cacomp *cpnt;

 if ((cpnt=(cacomp *)emalloc(sizeof(cacomp))) == NULL) {
    ncfprintf (stderr,"no space left for capnt %d\n",cumcacomp+1);
     return (NULL);
   }
 cumcacomp++;
 cpnt->next = NULL;

 upk_short(&cpnt->ctype,1,1);
 upk_short(&cpnt->stype,1,1);
 upk_int(&cpnt->num,1,1);		/* ca comp number */
 upk_int(&cpnt->vcnum,1,1);	/* number of voltage comp. */
 if (cpnt->vcnum==comp1->num) {
    cpnt->comp1 = comp1;
 }
 else cpnt->comp1 = NULL;

 upk_short(&cpnt->cashell,1,1);
 upk_double(cpnt->cais,cpnt->cashell,1);
 upk_double(&cpnt->cai,1,1);
 upk_double(&cpnt->ica,1,1);
 upk_double(&cpnt->casf0,1,1);
 upk_double(&cpnt->casfn,1,1);
 upk_double(&cpnt->casfc,1,1);
 upk_double(cpnt->caos,cpnt->caoshell,1);
 upk_double(&cpnt->cao,1,1);
 upk_double(&cpnt->casfno,1,1);
 upk_double(&cpnt->casfco,1,1);
 upk_double(&cpnt->vrev,1,1);
 upk_double(&cpnt->vmax,1,1);
 upk_double(&cpnt->pkm,1,1);
 upk_double(&cpnt->kex,1,1);
 upk_double(&cpnt->ekm,1,1);
 upk_double(&cpnt->cabnd,1,1);
 upk_double(&cpnt->ipump,1,1);
 upk_byte(&temp,1,1);			/* check for conlist */
 if (temp==VCONLST) cpnt->clst = upk_conlst();
 upk_byte(&temp,1,1);			/* end of structure */
 if (temp!=VEND) {
     ncfprintf (stderr,"upk_cacomp: invalid format for cacomp ('%d')\n",temp);
     return NULL;
  }
  return cpnt;
}

/*------------------------------------------------*/

comp *upk_comp ()

/* Unpack a voltage compartment from the receive buffer. */

{
     comp *cpnt;
     conlst *lpnt;
     char temp;
     
  if ((cpnt=(comp *)emalloc(sizeof(comp))) == NULL) {
     ncfprintf (stderr,"no space left for comp %d\n", cumcomp+1);
     return (NULL);  
  }
 cpnt->next = NULL;                    /* make new comp be end of list */
 if (!compnt) compnt = cpnt;           /* save head if first synap */
 if (compend)
    compend->next = cpnt;
 cpnt->last = compend;                 /* pointer to last compartment */
 compend = cpnt;
 cpnt->nodlst = NULL;                  /* pointer to list of nodes for comp */
 cpnt->capnt = NULL;                   /* pointer to possible ca comp */
 cumcomp++;

 upk_short(&cpnt->ctype,1,1);
 upk_short(&cpnt->miscfl,1,1);
 upk_int(&cpnt->num,1,1);
 upk_double(&cpnt->v,1,1);
 upk_double(&cpnt->oldv,1,1);
 upk_double(&cpnt->vest,1,1);
 upk_double(&cpnt->implf,1,1);
 upk_double(&cpnt->nodc,1,1);
 upk_double(&cpnt->tcond,1,1);
 upk_double(&cpnt->tcondn,1,1);
 upk_double(&cpnt->vrev,1,1);
 upk_double(&cpnt->verr,1,1);
 upk_double(&cpnt->k,1,1);
 upk_double(&cpnt->rm,1,1);
 upk_double(&cpnt->cap,1,1);
 upk_double(&cpnt->extv,1,1);
 upk_double(&cpnt->extv_old,1,1);
 upk_double(&cpnt->extvi,1,1);
 upk_double(&cpnt->exti,1,1);
 upk_short(&cpnt->g,1,1);
 upk_short(&cpnt->t,1,1);
 upk_short(&cpnt->virt,1,1);			/* virt flag 0 => not virtual */
 upk_double(&cpnt->relax,1,1);
 upk_byte(&temp,1,1);
 if (temp==VCACOMP) cpnt->capnt=upk_cacomp(cpnt); /* calcium comp goes along.*/ 
 upk_byte(&temp,1,1);
 if (temp==VCONLST) cpnt->clst=upk_conlst();
 upk_byte(&temp,1,1);			/* end of structure */
 if (temp!=VEND) {
     ncfprintf (stderr,"upk_comp: invalid format for comp ('%d')\n",temp);
     return NULL;
  }
  return cpnt;
}

/*------------------------------------------------*/

void conncomps()

/* Scan through list of compartments, checking the lists
   of connections.  For each connection, check the compartment
   numbers, find the compartments' addresses, and set the
   connection pointers.  */

{
    comp *cpnt;
    conlst *lpnt;
    conn *conpnt;

   for (cpnt=compnt; cpnt; cpnt=cpnt->next) {
      for (lpnt=cpnt->clst; lpnt; lpnt=lpnt->next) {
         conpnt = lpnt->conpnt;
         if (!conpnt) continue;
         if (!conpnt->comp1) 
             conpnt->comp1 = findnode (conpnt->num1,"conncomp1");
         switch (conpnt->ctype) {

       case AXIALRES:		/* 2 comp pointers */
       case GJ:
       case RESISTOR:
       case CAP:
       case BATT:
       case SYNAPSE:
       case BUF:
		if (!conpnt->comp2)
                  conpnt->comp2 = findnode (conpnt->num2,"conncomp2");
		break;

       case NA:			/* only one comp pointer */
       case K:
       case ClCa:
       case KCa:
       case CA:
       case ROD:
       case CONE:
       case CHR:
       case LOAD:
                  break;


        }
      }
   }
}

