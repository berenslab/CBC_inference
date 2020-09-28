/* Segment ncpacks in Program nc */

/* Packs nc structures for sending through parallel 
   message-passing routines.
*/

#include <stdio.h>
#include "nc.h"
#include "y.tab.h"
#include "ncomp.h"
#include "ncelem.h"
#include "ncvirt.h"
#include "nclist.h"
#include "pvm3.h"

pk_byte  (char *val,   int n, int stride);
pk_short (short *val,  int n, int stride);
pk_long  (long  *val,  int n, int stride);
pk_int   (int *val,    int n, int stride);
pk_float (float *val,  int n, int stride);
pk_double(double *val, int n, int stride);
pk_str   (char *val);


/*------------------------------------------------*/

void pk_stconc (stconc *cpnt)

/* Pack a "stconc" into the send buffer. */

{
  pk_double(&cpnt->cval,1,1);
  pk_double(&cpnt->cest,1,1);
  pk_double(&cpnt->dcon,1,1);
}

/*------------------------------------------------*/

void pk_chan (chan *cpnt)

/* Pack a chan into the send buffer. */

{
   int i;
   char one=1, zero=0, vchan=VCHAN, vend=VEND;

  pk_byte(&vchan,1,1);
  pk_short(&cpnt->ctype,1,1);		/* channel type */
  pk_short(&cpnt->stype,1,1);		/* channel sub-type */
  pk_int(&cpnt->num,1,1);
  pk_double(&cpnt->conduct,1,1);
  pk_int(&cpnt->num1,1,1);
  pk_int(&cpnt->num2,1,1);
  pk_double(&cpnt->maxcond,1,1);
  pk_double(&cpnt->voffsm,1,1);
  pk_double(&cpnt->voffsh,1,1);
  pk_double(&cpnt->vrev,1,1);
  /* reconstruct chtyp from na,k,ca,types[stype] */ 
  pk_double(&cpnt->cdur,1,1);
  if (cpnt->cstate) {
  	pk_byte(&one,1,1);		/* send 1 if cstate */
	pk_byte(cpnt->cstate,RNDSIZ,1);
  }
  else pk_byte(&zero,1,1);
  pk_short(&cpnt->nchan,1,1);
  pk_double(&cpnt->taua,1,1);
  pk_double(&cpnt->taub,1,1);
  switch (cpnt->ctype) {
    case NA:
      switch (cpnt->stype) {
       case 0:
         pk_double(&((hhchan *)cpnt)->m,1,1);
         pk_double(&((hhchan *)cpnt)->h,1,1);
	 break;
       case 1:
         pk_short(&((chan *)cpnt)->numstate,1,1);
         for (i=0; i<((chan *)cpnt)->numstate; i++) {
            pk_stconc(&((chan *)cpnt)->conc[i]);
         }
	 break;
       default:
         pk_double(&((hhchan *)cpnt)->m,1,1);
         pk_double(&((hhchan *)cpnt)->h,1,1);
	 break;
      }
      break;
    case K:
      switch (cpnt->stype) {
       case 0:
       case 2:
         pk_double(&((hhchan *)cpnt)->m,1,1);
         pk_double(&((hhchan *)cpnt)->h,1,1);
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
         pk_short(&((chan *)cpnt)->numstate,1,1);
         for (i=0; i<((chan *)cpnt)->numstate; i++) {
            pk_stconc(&((chan *)cpnt)->conc[i]);
         }
	 break;
      }
      break;

    case ClCa:
         pk_short(&((chan *)cpnt)->numstate,1,1);
         for (i=0; i<((chan *)cpnt)->numstate; i++) {
            pk_stconc(&((chan *)cpnt)->conc[i]);
         }
	 break;

    case KCa:
      switch (cpnt->stype) {
       case 0:
       case 2:
         pk_double(&((kcachan *)cpnt)->m,1,1);
         pk_double(&((kcachan *)cpnt)->h,1,1);
         pk_double(&((kcachan *)cpnt)->d1,1,1);
         pk_double(&((kcachan *)cpnt)->d2,1,1);
         pk_double(&((kcachan *)cpnt)->k1,1,1);
         pk_double(&((kcachan *)cpnt)->k2,1,1);
	 break;
       case 1:
       case 3:
       case 4:
       case 5:
       case 6:
      default:
         pk_short(&((chan *)cpnt)->numstate,1,1);
         for (i=0; i<((chan *)cpnt)->numstate; i++) {
            pk_stconc(&((chan *)cpnt)->conc[i]);
         }
	 break;
      }
      break;
    case CA:
      switch (cpnt->stype) {
       case 0:
         pk_double(&((hhchan *)cpnt)->m,1,1);
         pk_double(&((hhchan *)cpnt)->h,1,1);
	 break;
       default:
	 break;
       }
      break;
  }
  pk_byte(&vend,1,1);			/* end of structure */
}

/*------------------------------------------------*/

void pk_cn (conn *cpnt)

/* Pack the conn base into the send buffer. */

{
  pk_short(&cpnt->ctype,1,1);
  pk_double(&cpnt->conduct,1,1);
  pk_int(&cpnt->num1,1,1);
  pk_int(&cpnt->num2,1,1);
}
/*------------------------------------------------*/

void pk_lpfilt (lpfilt *cpnt)

/* Pack an lpfilt into the send buffer. */

{
  pk_short(&cpnt->nfilt,1,1);
  pk_double(&cpnt->ftau,1,1);
  pk_double(cpnt->lfilt,NUMFILT+1,1);
  pk_double(&cpnt->tfall,1,1);
}

/*------------------------------------------------*/

void pk_ld (load *cpnt)

/* Pack a "load" into the send buffer. */

{
  pk_short(&cpnt->ctype,1,1);
  pk_double(&cpnt->conduct,1,1);
  pk_int(&cpnt->num1,1,1);
  pk_double(&cpnt->vrev,1,1);
}

/*------------------------------------------------*/

void pk_load (load *cpnt)

/* Pack a "load" into the send buffer. */

{
   char vload=VLOAD,vend=VEND;

  pk_byte(&vload,1,1);
  pk_ld(cpnt);
  pk_byte(&vend,1,1);
}

/*------------------------------------------------*/

void pk_dbuf (dbuf *cpnt)

/* Pack a "dbuf" into the send buffer. */

{
   char vbuf=VBUF,vend=VEND;
   short dbufsiz;

  pk_byte(&vbuf,1,1);
  pk_cn(cpnt);
  pk_double(&cpnt->delay,1,1);
  pk_short(cpnt->delbuf,((int)cpnt->delay),1);
  dbufsiz = (cpnt->delbuf-cpnt->delpnt)/sizeof(short);
  pk_short(&dbufsiz,1,1);
  pk_byte(&vend,1,1);
}

/*------------------------------------------------*/

void pk_recep (recep *cpnt)

/* Pack a synap into the send buffer. */

{
   char one=1,zero=0,vrecep=VRECEP,vend=VEND;

  pk_byte(&vrecep,1,1);
  pk_ld((load*)cpnt);
  pk_short(&cpnt->chtyp->pigm,1,1);

  /* determine recep chtyp from ctype */
  pk_double(&cpnt->xloc,1,1);
  pk_double(&cpnt->yloc,1,1);
  pk_double(&cpnt->pathl,1,1);
  pk_double(&cpnt->attf,1,1);
  pk_double(&cpnt->timec,1,1);
  pk_double(&cpnt->dnois,1,1);
  pk_short(&cpnt->filt,1,1);
  pk_short(&cpnt->pnois,1,1);
  if (cpnt->dstate) {
  	pk_byte(&one,1,1);		/* send 1 if dstate */
	pk_byte(cpnt->dstate,RNDSIZ,1);
  }
  else pk_byte(&zero,1,1);
  if (cpnt->pstate) {
  	pk_byte(&one,1,1);		/* send 1 if cstate */
	pk_byte(cpnt->pstate,RNDSIZ,1);
  }
  else pk_byte(&zero,1,1);
  pk_double(&cpnt->area,1,1);
  pk_double(&cpnt->aflux,1,1);
  pk_double(&cpnt->mflux,1,1);
  pk_double(&cpnt->iflux,1,1);
  pk_double(&cpnt->rhod,1,1);
  pk_double(&cpnt->gpr1,1,1);
  pk_double(&cpnt->gpr2,1,1);
  pk_double(&cpnt->gpr3,1,1);
  pk_double(&cpnt->pde,1,1);
  pk_double(&cpnt->gcyc,1,1);
  pk_double(&cpnt->cycg,1,1);
  pk_double(&cpnt->cond,1,1);
  pk_double(&cpnt->ca,1,1);
  pk_double(&cpnt->cax,1,1);
  pk_double(&cpnt->maxcond,1,1);
  pk_short(&cpnt->recnm1,1,1);
  pk_short(&cpnt->recnm2,1,1);
  pk_short(&cpnt->recnm3,1,1);
  pk_byte(&vend,1,1);
}

/*------------------------------------------------*/

void pk_synap (synap *cpnt)

/* Pack a synap into the send buffer. */

{
   char one=1,zero=0,vsynap=VSYNAP,vend=VEND;

  pk_byte(&vsynap,1,1);
  pk_cn(cpnt);
  pk_double(&cpnt->vrev,1,1);
  pk_double(&cpnt->thresh,1,1);
  pk_double(&cpnt->nkd,1,1);
  pk_double(&cpnt->ngain,1,1);
  pk_double(&cpnt->vgain,1,1);
  pk_double(&cpnt->ckd,1,1);
  pk_double(&cpnt->chc,1,1);
  pk_double(&cpnt->cgain,1,1);
  pk_double(&cpnt->coff,1,1);
  pk_double(&cpnt->maxcond,1,1);
  pk_double(&cpnt->maxsrate,1,1);
  pk_double(&cpnt->rrpool,1,1);
  pk_double(&cpnt->maxrrpool,1,1);
  pk_double(&cpnt->vsize,1,1);
  pk_double(&cpnt->vcov,1,1);
  pk_double(&cpnt->cov,1,1);
  pk_short(&cpnt->setvrev,1,1);
  pk_short(&cpnt->vrefr,1,1);
  pk_short(&cpnt->vflag,1,1);
  pk_short(&cpnt->vsites,1,1);
  pk_short(&cpnt->sens,1,1);
  pk_short(&cpnt->ntact,1,1);
  pk_short(&cpnt->secmsg,1,1);
  pk_short(&cpnt->mesg1,1,1);
  pk_short(&cpnt->mesg2,1,1);
  pk_short(&cpnt->curve,1,1);
  pk_int(&cpnt->sdyadnum,1,1);
  if (cpnt->cstate) {
  	pk_byte(&one,1,1);		/* send 1 if cstate */
	pk_byte(cpnt->cstate,RNDSIZ,1);
  }
  else pk_byte(&zero,1,1);
  if (cpnt->vstate) {
  	pk_byte(&one,1,1);		/* send 1 if vstate */
	pk_byte(cpnt->vstate,RNDSIZ,1);
  }
  else pk_byte(&zero,1,1);
  if (cpnt->filt1) {
  	pk_byte(&one,1,1);		/* send 1 if filt1 */
	pk_lpfilt(cpnt->filt1);
  }
  else pk_byte(&zero,1,1);
  if (cpnt->filt2) {
  	pk_byte(&one,1,1);		/* send 1 if filt2 */
	pk_lpfilt(cpnt->filt2);
  }
  else pk_byte(&zero,1,1);
  if (cpnt->filt3) {
  	pk_byte(&one,1,1);		/* send 1 if filt3 */
	pk_lpfilt(cpnt->filt3);
  }
  else pk_byte(&zero,1,1);
  pk_byte(&vend,1,1);			/* end of structure */
}

/*------------------------------------------------*/

void pk_conn (conn *cpnt)

/* Pack a conn into the send buffer. */

{
   char vconn=VCONN,vend=VEND;

  pk_byte(&vconn,1,1);
  pk_cn(cpnt);
  pk_byte(&vend,1,1);			/* end of structure */
}

/*------------------------------------------------*/

void pk_conlst (conlst *cpnt, comp* comp1)

/* Pack a list of connections into the send buffer. */

{
   conlst *lpnt;
   char vconlst=VCONLST, one=1, zero=0,vend=VEND;

 if (!lpnt) {
    pk_byte(&zero,1,1);
    return;
 }
 else {
  pk_byte (&vconlst,1,1);
  for (lpnt=cpnt; lpnt; lpnt=lpnt->next) {
    if (!lpnt->conpnt) continue;
    switch (lpnt->conpnt->ctype) {  /* check that we don't send conn twice */
       case AXIALRES:
       case GJ:
       case RESISTOR:
       case CAP:
       case BATT:
       case SYNAPSE:
       case BUF:
		if (lpnt->conpnt->comp1 != comp1) continue; /* skip second */
                break;
       default: break;
    }
    switch (lpnt->conpnt->ctype) {

       case AXIALRES:
       case GJ:
       case RESISTOR:
       case CAP:
       case BATT:
		  pk_conn(lpnt->conpnt);
		  break;

       case NA:
       case K:
       case ClCa:
       case KCa:
       case CA:
		  pk_chan((chan*)lpnt->conpnt);
		  break;
       case SYNAPSE:
		  pk_synap((synap*)lpnt->conpnt);
		  break;
       case ROD:
       case CONE:
       case CHR:
		  pk_recep((recep*)lpnt->conpnt);
		  break;
       case LOAD:
		  pk_load((load*)lpnt->conpnt);
		  break;
       case BUF:
		  pk_dbuf((dbuf*)lpnt->conpnt);
		  break;

    }
  }
 pk_byte (&vend,1,1);
 }
}

/*------------------------------------------------*/

void pk_cacomp (cacomp *cpnt)

/* Pack a calcium compartment into the send buffer. */

{
    conlst *lpnt;
    char vcacomp=VCACOMP,vconlst=VCONLST,vend=VEND,zero=0;

 pk_byte(&vcacomp,1,1);		/* sending ca comp */
 pk_short(&cpnt->ctype,1,1);
 pk_short(&cpnt->stype,1,1);
 pk_int(&cpnt->num,1,1);
 pk_int(&cpnt->vcnum,1,1);	/* number of voltage comp. */
 pk_short(&cpnt->cashell,1,1);
 pk_double(cpnt->cais,cpnt->cashell,1);
 pk_double(&cpnt->cai,1,1);
 pk_double(&cpnt->ica,1,1);
 pk_double(&cpnt->casf0,1,1);
 pk_double(&cpnt->casfn,1,1);
 pk_double(&cpnt->casfc,1,1);
 pk_double(cpnt->caos,cpnt->caoshell,1);
 pk_double(&cpnt->cao,1,1);
 pk_double(&cpnt->casfno,1,1);
 pk_double(&cpnt->casfco,1,1);
 pk_double(&cpnt->vrev,1,1);
 pk_double(&cpnt->vmax,1,1);
 pk_double(&cpnt->pkm,1,1);
 pk_double(&cpnt->kex,1,1);
 pk_double(&cpnt->ekm,1,1);
 pk_double(&cpnt->cabnd,1,1);
 pk_double(&cpnt->ipump,1,1);
 if (cpnt->clst) {
	pk_byte(&vconlst,1,1);		/* Send 1 if conlst */
	pk_conlst(cpnt->clst,NULL);	/* No connections to cacomp */
  }
  else pk_byte(&zero,1,1);
 pk_byte(&vend,1,1);			/* end of structure */
}

/*------------------------------------------------*/

void pk_comp (comp *cpnt)

/* Pack a compartment into the send buffer. */

{
     conlst *lpnt;
     char zero=0,vcomp=VCOMP,vend=VEND;
     
 pk_byte(&vcomp,1,1);			/* sending comp */
 pk_short(&cpnt->ctype,1,1);
 pk_short(&cpnt->miscfl,1,1);
 pk_int(&cpnt->num,1,1);
 pk_double(&cpnt->v,1,1);
 pk_double(&cpnt->oldv,1,1);
 pk_double(&cpnt->vest,1,1);
 pk_double(&cpnt->implf,1,1);
 pk_double(&cpnt->nodc,1,1);
 pk_double(&cpnt->tcond,1,1);
 pk_double(&cpnt->tcondn,1,1);
 pk_double(&cpnt->vrev,1,1);
 pk_double(&cpnt->verr,1,1);
 pk_double(&cpnt->k,1,1);
 pk_double(&cpnt->rm,1,1);
 pk_double(&cpnt->cap,1,1);
 pk_double(&cpnt->extv,1,1);
 pk_double(&cpnt->extv_old,1,1);
 pk_double(&cpnt->extvi,1,1);
 pk_double(&cpnt->exti,1,1);
 pk_short(&cpnt->g,1,1);
 pk_short(&cpnt->t,1,1);
 pk_short(&cpnt->virt,1,1);		/* virt flag 0 => not virtual */
 pk_double(&cpnt->relax,1,1);
 if (cpnt->capnt) pk_cacomp(cpnt->capnt);/* calcium compartment goes along. */ 
 if (cpnt->clst)  pk_conlst(cpnt->clst,cpnt);
 pk_byte(&vend,1,1);			/* end of structure */

}

