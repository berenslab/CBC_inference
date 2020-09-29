/* Module ncsave in Program nc */

/* saves parameters */

#include <sys/types.h>
#include <unistd.h>
#include <string.h>

#include "ncsub.h"
#include "ncomp.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
}
#endif

#include "nc.h"
#include "y.tab.h"
#include "control.h"

#include "ncio.h"

extern comp  *compnt;
extern synap *synpnt;
extern photrec *recpnt;
extern chan *chanpnt;

extern int cumcomp;
extern int cumsynap;
extern int cumgj;
extern int cumphotrec;
extern int cumchan;
extern int cumdbuf;

extern char *infile;

void execerror(const char *s, const char *t);
char *emalloc (unsigned int n);
comp *findcomp(int num, const char *str);
int findphotrec(int num1, int num2, int num3, int num4, photrec **rpnt, const char *str);
char *print_version(double version);

static char ebuf[200];

/*---------------------------------------*/

void savemodel (const char *filnam)

/* save all state variables into a file for later restoration */

{
   int i, nf, cashell,caoshell, delay;
   double gain;
   double v, offset;
   comp *cpnt;
   photrec *rpnt;
   synap *spnt;
   chan *chpnt;
   FILE *sfil;
   cacomp *capnt;
   ntcomp *ntpnt,*npt;
   lpfilt *lpf;
   photreci *rpnti;
   phvars *vpnt;
   phvarsi *vpnti;
   
 if      (strcmp(filnam,"stdout") ==0) sfil = stdout;
 else if (strcmp(filnam,"stderr") ==0) sfil = stderr;
 else if ((sfil=fopen(filnam,"w"))==NULL) {
    sprintf (ebuf,"savemodel: can't open file '%.50s'.\n",filnam);
    execerror ("warning,",ebuf);
    return; 
 }

 fprintf (sfil,"## nc version %s pid %d running '%s'\n", 
		print_version(ncversion), getpid(),infile);
 fprintf (sfil,"## model save at %g sec\n",simtime);
 fprintf (sfil,"## comps=%d chans=%d syns=%d photorecs=%d\n",
			cumcomp,cumchan,cumsynap+cumgj+cumdbuf,cumphotrec);

 if (cumcomp >= 1) {
  fprintf (sfil,"# compartments num v relax\n");
  for (cpnt=compnt; cpnt; cpnt=cpnt->next) {            /* save compartments */
    fprintf (sfil,"c   %4d %.16g %.16g %.8g\n",cpnt->num,cpnt->v,cpnt->oldv,cpnt->relax);
    if (capnt=cpnt->capnt) {
      cashell = capnt->cashell;
      fprintf (sfil,"cai vr %.16g sh %d ",capnt->vrev,cashell);
      if (capnt->cab) {
        fprintf (sfil,"b %d",1);
        for (i=0; i<cashell; i++) {
          fprintf (sfil," %.16g %.16g",capnt->cais[i],capnt->cab[i]);
        }
      }
      else {
        fprintf (sfil,"b %d",0);
      	  for (i=0; i<cashell; i++) {
            fprintf (sfil," %.16g",capnt->cais[i]);
	  }
      }
      caoshell = capnt->caoshell;
      fprintf (sfil," cao sh %d ",caoshell);
      if (capnt->caos!=NULL)
        for (i=0; i<caoshell; i++) {
           fprintf (sfil," %.16g",capnt->caos[i]);
        }
      fprintf (sfil," %.16g %.16g",capnt->cas,capnt->cas2);	/* CICR local Ca store */
      fprintf (sfil,"\n");
    }  /* if capnt */
    if (ntpnt=cpnt->ntpnt) {
       for (npt=ntpnt; npt; npt=npt->ntpnt) {        /* save nt list */
          fprintf (sfil,"cn    %4d c %.8g\n",npt->ctype,npt->val);
       } 
    }
  } /* for (cpnt;;) */
  fprintf (sfil,"#  \n");
 } /* if (cumcomp>0) */

 if (cumchan > 0) {
  fprintf (sfil,"# kchannels ctype stype num conduct vrev nchan numstate states \n");
  for (chpnt=chanpnt; chpnt; chpnt=(chan*)chpnt->next) {         /* channels */
    fprintf (sfil,"c typ %d %d n %d g %.16g vr %.16g nch %g st %d ", 
			chpnt->ctype,chpnt->stype,
			chpnt->num,chpnt->conduct,chpnt->vrev,
			chpnt->nchan, chpnt->numstate);
    if (chpnt->chtyp->hh) {
      fprintf (sfil,"m %.16g h %.16g ",((hhchan*)chpnt)->m, ((hhchan*)chpnt)->h); /* HH */
    } 
    if (chpnt->cstate!=NULL) {
	 int rndsiz=RNDSIZ/sizeof(int);
      fprintf (sfil,"rn(%d) ",rndsiz); 			/* rng state */
      for (i=0; i<rndsiz; i++) {
        fprintf (sfil,"%d ",*((int *)&chpnt->cstate[i*sizeof(int)]));
      }
    }
    if (chpnt->nchan > 0) { 		/* if noise, print nchans in each state */
       for (i=0; i<chpnt->numstate; i++) {
	 fprintf (sfil,"%.16g %d ",chpnt->conc[i].cest,chpnt->conc[i].nchan);
       }
    } else {				/* no noise, don't print out nchan */
       for (i=0; i<chpnt->numstate; i++) {
	 fprintf (sfil,"%.16g ",chpnt->conc[i].cest);
       }
    }
    fprintf (sfil,"\n");
  }
  fprintf (sfil,"#\n");
 }

 if (cumsynap > 0 || cumgj > 0 || cumdbuf) {
  fprintf (sfil,"# synapses num conduct 1 1h 2 3\n");
  for (spnt=synpnt; spnt; spnt=(synap *)spnt->next) {      /* save synapses */

   switch (spnt->ctype) {

    case SYNAPSE:
    fprintf (sfil,"s    %4d %.16g %.8g %.8g %d %.8g ",
			spnt->num,spnt->conduct,
			spnt->vtint,spnt->vtime,spnt->vflag,
			spnt->rrpool);

    if (spnt->vstate!=NULL) {
	 int rndsiz=RNDSIZ/sizeof(int);
      fprintf (sfil,"vt(%d) ",rndsiz); 			/* ves timing rng state */
      for (i=0; i<rndsiz; i++) {
        fprintf (sfil,"%d ",*((int *)&spnt->vstate[i*sizeof(int)]));
      }
    }
    if (spnt->gstate!=NULL) {
	 int rndsiz=RNDSIZ/sizeof(int);
      fprintf (sfil,"vs(%d) ",rndsiz); 			/* ves size rng state */
      for (i=0; i<rndsiz; i++) {
        fprintf (sfil,"%d ",*((int *)&spnt->gstate[i*sizeof(int)]));
      }
    }
    fprintf (sfil,"\n");  				/* end of first line */
    if (lpf=spnt->filt1) {
       fprintf (sfil,"sfa  %4d",lpf->nfilt);
       for (i=0; i<lpf->nfilt; i++) {
            fprintf (sfil," %.8g",lpf->lfilt[i]);
       } 
       fprintf (sfil,"\n");
    }
    if (lpf=spnt->filt1h) {
       fprintf (sfil,"sfah %4d",lpf->nfilt);
       for (i=0; i<lpf->nfilt; i++) {
            fprintf (sfil," %.8g",lpf->lfilt[i]);
       } 
       fprintf (sfil,"\n");
    }
    if (lpf=spnt->filt2) {
       fprintf (sfil,"sfb  %4d",lpf->nfilt);
       for (i=0; i<lpf->nfilt; i++) {
            fprintf (sfil," %.8g",lpf->lfilt[i]);
       } 
       fprintf (sfil,"\n");
    }
    if (lpf=spnt->filt3) {
       fprintf (sfil,"sfc  %4d",lpf->nfilt);
       for (i=0; i<lpf->nfilt; i++) {
            fprintf (sfil," %.8g",lpf->lfilt[i]);
       } 
       fprintf (sfil,"\n");
    }
    break;

    case GJ:
     fprintf (sfil,"g    %4d %.16g %.8g\n", spnt->num,spnt->conduct,((gj*)spnt)->n);
    break;

    case BUF:
     gain = ((dbuf *)spnt)->gain;
     offset = ((dbuf *)spnt)->offset;
     v = ((dbuf *)spnt)->v;
     delay = int(((dbuf *)spnt)->delay);
     if ((lpf=((dbuf*)spnt)->filt)!=NULL) 
           nf = lpf->nfilt;
     else  nf = 0; 
     fprintf (sfil,"b    %4d %.8g %.8g %.8g %d %d", spnt->num, gain, offset, v, delay, nf);
     if (delay > 0) {
       fprintf (sfil," %ld", (long int)(((dbuf *)spnt)->delpnt));
       for (i=0; i<delay; i++) {
	 fprintf (sfil," %.8g",((dbuf *)spnt)->delbuf[i]);
       }
     }
     if (nf > 0) {
       for (i=0; i<nf; i++) {
            fprintf (sfil," %.8g",lpf->lfilt[i]);
       } 
     }
     fprintf (sfil,"\n");
    break;

   } /* switch (spnt->ctype) */
  }
  fprintf (sfil,"#\n");
 }

 if (cumphotrec > 0) {
  for (rpnt=recpnt; rpnt; rpnt=(photrec *)rpnt->next) {    /* save photorecs */
     if (rpnt->chtyp->stype < PINVG) {
       vpnt = &rpnt->vars;
  fprintf (sfil,"# photoreceptors typ pigm node(4) conduct rhod rhodk rstar gpr pde cycg gcyc cond ca cax ksens\n");
      fprintf (sfil,"p %d %d %-d %-d %-d %-d %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g\n",
	rpnt->chtyp->stype,
	rpnt->chtyp->pigm,
	(rpnt->recnm1<0?-1:rpnt->recnm1),
	(rpnt->recnm2<0?-1:rpnt->recnm2),
	(rpnt->recnm3<0?-1:rpnt->recnm3),
	(rpnt->recnm4<0?-1:rpnt->recnm4),
	rpnt->conduct,
	vpnt->rhod,
	vpnt->rhodk,
	vpnt->rhod2,
	vpnt->rstar,
	vpnt->gpr,
	vpnt->pde,
	vpnt->cycg,
	vpnt->gcyc,
	vpnt->cond,
	vpnt->ca,
	vpnt->cax,
	vpnt->cab,
	vpnt->ksens,
	vpnt->resp_r,
	vpnt->resp_b,
	vpnt->resp_e,
	vpnt->resp_x,
	vpnt->resp_c,
	vpnt->resp_is,
	vpnt->atten_i
	);

    } else if (rpnt->chtyp->stype==PINVG) {

       rpnti = (photreci*)rpnt;		// invergo type mouse rod photorec 
       vpnti = &rpnti->vars;
  fprintf (sfil,"# photoreceptors typ pigm node(4) conduct Rh R_Gt Gt Rn0 Rn1 Rn2 Rn3 Rn4 Rn5 RK Rn_RKpre0 Rn_RKpre1 Rn_RKpre2 Rn_RKpre3 Rn_RKpre4 Rn_RKpre5 Rn_RKpre6 Rn_RKpost0 Rn_RKpost1 Rn_RKpost2 Rn_RKpost3 Rn_RKpost4 Rn_RKpost5 Rn_RKpost6 ");
  fprintf (sfil,"Arr Rn_Arr0 Rn_Arr1 Rn_Arr2 Rn_Arr3 Rn_Arr4 Rn_Arr5 Rn_Arr6 Ops Arr_di Arr_tetra Ops_Gt Ops_G Ops_Ggtp Ggtp Rn_Gt0 Rn_Gt1 Rn_Gt2 Rn_Gt3 Rn_Gt4 Rn_Gt5 Rn_Gt6 Rn_G0 Rn_G1 Rn_G2 Rn_G3 Rn_G4 Rn_G5 Rn_G6 ");

  fprintf (sfil,"Rn_Gtp0 Rn_Gtp1 Rn_Gtp2 Rn_Gtp3 Rn_Gtp4 Rn_Gtp5 Rn_Gtp6 Gagtp Gbg PDE PDE_Gagtp PDEa_Gagtp Gagtp_PDEa_Gagtp Gagtp_aPDEa_Gagtp RGS RGS_PDEa_Gagtp RGS_Gagtp_aPDEa_Gagtp Gagdb Rect Recr_Ca Cafree Recr_Ca_RK Cabuff cGMP cond\n");

      fprintf (sfil,"p %d %d %-d %-d %-d %-d %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g ",
	rpnti->chtyp->stype,
	rpnti->chtyp->pigm,
	(rpnt->recnm1<0?-1:rpnt->recnm1),
	(rpnt->recnm2<0?-1:rpnt->recnm2),
	(rpnt->recnm3<0?-1:rpnt->recnm3),
	(rpnt->recnm4<0?-1:rpnt->recnm4),
	rpnt->conduct,
	vpnti->Rh,                      // unphosphorylated rhodopsin
	vpnti->R_Gt,                    //
	vpnti->Gt,                      //
	vpnti->Rn[0],                   // phosphorylated rhodopsin
	vpnti->Rn[1],                   // 
	vpnti->Rn[2],                   // 
	vpnti->Rn[3],                   // 
	vpnti->Rn[4],                   // 
	vpnti->Rn[5],                   // 
	vpnti->Rn[6],                   // 
	vpnti->RK,                      // rhodopsin kinase
	vpnti->Rn_RKpre[0],              // 
	vpnti->Rn_RKpre[1],
	vpnti->Rn_RKpre[2],
	vpnti->Rn_RKpre[3],
	vpnti->Rn_RKpre[4],
	vpnti->Rn_RKpre[5],
	vpnti->Rn_RKpre[6],
	vpnti->Rn_RKpost[0],
	vpnti->Rn_RKpost[1],
	vpnti->Rn_RKpost[2],
	vpnti->Rn_RKpost[3],
	vpnti->Rn_RKpost[4],
	vpnti->Rn_RKpost[5],
	vpnti->Rn_RKpost[6]
	);

      fprintf (sfil,"%-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g ",

	vpnti->Arr,                     // arrestin 
	vpnti->Rn_Arr[0],               // rhodopsin-arrestin complex
	vpnti->Rn_Arr[1],               //
	vpnti->Rn_Arr[2],               //
	vpnti->Rn_Arr[3],               //
	vpnti->Rn_Arr[4],               //
	vpnti->Rn_Arr[5],               //
	vpnti->Rn_Arr[6],               //
	vpnti->Ops,                     // ligand-free receptor 
	vpnti->Arr_di,                  // arrestin dimer
	vpnti->Arr_tetra,               // arrestin tetramer
	vpnti->Ops_Gt,                  // receptor Gt complex
	vpnti->Ops_G,                   // receptor G complex
	vpnti->Ops_Ggtp,                // receptor G complex
	vpnti->Ggtp,                    // Ggtp 
	vpnti->Rn_Gt[0],                // receptor G complex
	vpnti->Rn_Gt[1],                // 
	vpnti->Rn_Gt[2],                // 
	vpnti->Rn_Gt[3],                // 
	vpnti->Rn_Gt[4],                // 
	vpnti->Rn_Gt[5],                // 
	vpnti->Rn_Gt[6],                // 
	vpnti->Rn_G[0],                 // receptor G complex
	vpnti->Rn_G[1],                 // 
	vpnti->Rn_G[2],                 // 
	vpnti->Rn_G[3],                 // 
	vpnti->Rn_G[4],                 // 
	vpnti->Rn_G[5],                 // 
	vpnti->Rn_G[6]                  // 
	);

      fprintf (sfil,"%-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g\n",

	vpnti->Rn_Ggtp[0],              // receptor G complex
	vpnti->Rn_Ggtp[1],              // 
	vpnti->Rn_Ggtp[2],              // 
	vpnti->Rn_Ggtp[3],              // 
	vpnti->Rn_Ggtp[4],              // 
	vpnti->Rn_Ggtp[5],              // 
	vpnti->Rn_Ggtp[6],              // 
	vpnti->Gagtp,                   // 
	vpnti->Gbg,                     // 
	vpnti->PDE,                     // 
	vpnti->PDE_Gagtp,               // PDE-Gagtp complex
	vpnti->PDEa_Gagtp,              // activated PDE-Gagtp complex
	vpnti->Gagtp_PDEa_Gagtp,        // activated PDE-Gagtp complex
	vpnti->Gagtp_aPDEa_Gagtp,       // activated PDE-Gagtp complex
	vpnti->RGS,                     // 
	vpnti->RGS_PDEa_Gagtp,          // activated RGS-PDE-Gagtp complex
	vpnti->RGS_Gagtp_aPDEa_Gagtp,   // activated RGS-PDE-Gagtp complex
	vpnti->Gagdp,                   // Gagdp 
	vpnti->Rect,                    // 
	vpnti->Recr_Ca,                 // 
	vpnti->Cafree,                  // Ca free
	vpnti->Recr_Ca_RK,              // 
	vpnti->Cabuff,                  // Ca buffered
	vpnti->cGMP,                    // cyclic GMP
	vpnti->cond                     // dark conductance normalized 0 to 1
	);

    } else if (rpnt->chtyp->stype==ChR2) {
   	   photrecc *rpntc;
   	   phvarsc *vpntc;

       rpntc = (photrecc*)rpnt;		// invergo type mouse rod photorec 
       vpntc = &rpntc->vars;
      fprintf (sfil,"# photoreceptors typ pigm node(4) conduct \n");
      fprintf (sfil,"p %d %d %-d %-d %-d %-d %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g %-.9g\n",
	rpnti->chtyp->stype,
	rpnti->chtyp->pigm,
	(rpnt->recnm1<0?-1:rpnt->recnm1),
	(rpnt->recnm2<0?-1:rpnt->recnm2),
	(rpnt->recnm3<0?-1:rpnt->recnm3),
	(rpnt->recnm4<0?-1:rpnt->recnm4),
	rpnt->conduct,
	vpntc->SC1,                     // state C1 
	vpntc->SO1,                     // state O1 
	vpntc->SO2,                     // state O2 
	vpntc->SC2,                     // state C2 
	vpntc->p,                       // nonlinear intensity factor
	vpntc->cond                     // dark conductance normalized 0 to 1
	);
    }
  }
  fprintf (sfil,"#\n");
 }
 fclose(sfil);
} /* savemodel */

/*---------------------------------------*/

char *scantoken (char *cp, int n)

/* find n'th text token */

{
    int i;
    char *np;

  for (i=0; i<n; i++) {
     np = strtok (cp," ,\t\n");
     cp = NULL;
  }
  return np;
}

/*---------------------------------------*/

#define INBUFSIZ 4000

void restoremodel (const char *filnam)

/* restore all state variables from a file */

{
   int i, caoshell, compnum, n, npl, nf, ntype, snum, chnum;
   int tcomps, tchans, tsyns, tphotrecs;
   int vflag, delay;
   int cashell, cabfl;
   int ctype, stype, numstate;
   long int bpnt;
   size_t inbufsiz=INBUFSIZ;
   double gain;
   double offset, v; 
   double cais,cab;
   double val, scond, gjfrac,rrpool,vtint,vtime;
   double compv, oldv, relax;
   double nchan, chcond, vrev;
   int recnm1, recnm2, recnm3, recnm4;
   comp *cpnt,*tcpnt;
   FILE *sfil;
   static char *inbuf;
   cacomp *capnt;
   ntcomp *ntpnt,*npt;
   synap *spnt,*tspnt;
   lpfilt *lpf;
   photrec *rpnt;
   photreci *rpnti;
   chan *chpnt;
   char *cp;
   phvars *vpnt;
   phvarsi *vpnti;
   int phtype, pigm;

  if (!(sfil=fopen(filnam,"r"))) {
    sprintf (ebuf,"restoremodel: can't open file '%.50s'.\n",filnam);
    execerror ("warning,",ebuf);
    return; 
  }

  if (!(inbuf = (char *)emalloc (INBUFSIZ))) {
     ncfprintf (stderr,"Error, restoremodel: can't allocate line buffer\n");
     return;
  }

  tcomps = tchans = tsyns = tphotrecs = 0;
  for (npl=0; getline(&inbuf,&inbufsiz,sfil)>0; npl++) {		/* scan file */

   if (inbuf[0]=='#')  {				/* comps, synapses, photrecs */
     switch (inbuf[2])  {

      case ' ':	break;

      case 'c':						/* found "# compartments */
       if ((cpnt=compnt)==NULL) break;
       for (; getline(&inbuf,&inbufsiz,sfil)>0; npl++) {
          if (inbuf[0]=='#') break;  			/* stop at last comps */
          switch (inbuf[0])  {

         case 'c':					/* comp line */
           switch (inbuf[1])  {
            case ' ':					/* compartment */
            sscanf(inbuf,"c %d %lg %lg %lg",&compnum,&compv,&oldv,&relax);
	    if (cpnt->num != compnum) {			/* if not in correct order */
	       cpnt=findcomp(compnum,"restore_model");
	    }
	    if (cpnt!=NULL) {
		cpnt->v     = compv;
		cpnt->oldv  = oldv;
		cpnt->relax = relax;
		tcpnt = cpnt;				/* save ptr for calcium, nt */
	        cpnt = cpnt->next;			/* fast way of finding next comp */
		if (cpnt==NULL) cpnt=compnt;
		tcomps++;
	    }
 	    break;

            case 'a':		/* "ca" line */
	      if ((capnt=tcpnt->capnt)!=NULL) {
                n = sscanf(inbuf,"%*s vr %lg sh %d b %d",&vrev,&cashell,&cabfl); 
		capnt->vrev = vrev;
		cp=scantoken(inbuf,8);			/* find 8th token */
		if ((cashell<=capnt->cashell) && capnt->cais!=NULL) {
		  if (cabfl) {
		    for (i=0; i<cashell; i++) {
    		      if (cp==NULL) {
			sprintf (ebuf,"restoremodel: error in comp %d reading (buf) caicomp: %d.\n",compnum,i);
     	                execerror ("warning,",ebuf);
		      }
                      capnt->cais[i] = atof(cp);
		      cp=scantoken(NULL,1);
		      capnt->cab[i] = atof(cp); 
		      cp=scantoken(NULL,1);
		    }
		  } else {
		    for (i=0; i<cashell; i++) { 
    		      if (cp==NULL) {
			sprintf (ebuf,"restoremodel: error in comp %d reading caicomp: %d.\n",compnum,i);
     	                execerror ("warning,",ebuf);
		      }
                      capnt->cais[i] = atof(cp);
		      // fprintf (stderr,"%g\n",capnt->cais[i]);
		      cp=scantoken(NULL,1);
		    }
		  }
                  // n = sscanf(inbuf,"%*s sh %d",&caoshell); 
		  cp=scantoken(NULL,2);			/* skip "cao sh" */
                  caoshell = atof(cp);
		  cp=scantoken(NULL,1);			/* get caoshell */
		  if ((caoshell<=capnt->caoshell) && (capnt->caos!=NULL)) {
		    cp=scantoken(NULL,1);			/* get cao number */
		    for (i=0; i<caoshell; i++) {
    		      if (cp==NULL) {
			sprintf (ebuf,"restoremodel: error reading caocomp: %d.\n",i);
     	                execerror ("warning,",ebuf);
		      }
                      capnt->caos[i] = atof(cp);
		      cp=scantoken(NULL,1); 			/* get cao number */
		    }
		  }
                  capnt->cas = atof(cp);			/* get CICR Ca store */
		  cp=scantoken(NULL,1);
                  capnt->cas2 = atof(cp);
		  cp=scantoken(NULL,1);
		}
	      }
 	    break;

            case 'n':		/* "cn" line */
	      if ((ntpnt=tcpnt->ntpnt)!=NULL) {
	        if (sscanf(inbuf,"%*s %d c %lg\n",&ntype,&val) ==2) {
		   while (ntpnt->ctype!=ntype) {
		     ntpnt=ntpnt->ntpnt;
		     if (ntpnt==NULL) break;
		   }
		   if (ntpnt!=NULL) {
		     ntpnt->val = val;			/* set nt conc */
		   }
	        }
	      }
 	    break;
	    default: break;
	  
	      }   /* switch (inbuf[1]) */
 	    break;
	    default: break;
	   }  /* switch (inbuf[0]) */
          }  /* for all comp lines */
      break;  /* " compartment c?" */

      case 'k':						/* found "# kchannels" */
       if ((chpnt=chanpnt) == NULL) {
    	 sprintf (ebuf,"restoremodel: no channels\n");
     	 execerror ("warning,",ebuf);
       }
       for (; getline(&inbuf,&inbufsiz,sfil)>0; npl++) {
        if (inbuf[0]=='#') break;  			/* stop at last chans */
        switch (inbuf[0])  {

         case 'c':					/* chan */
	    n = sscanf (inbuf,"%*s %*s %d %d n %d g %lg vr %lg nch %lg st %d ", 
			&ctype, &stype, &chnum, &chcond, &vrev, &nchan, &numstate);
	    if (n != 7) {
    		sprintf (ebuf,"restoremodel: can't read channel.\n",filnam);
     	        execerror ("warning,",ebuf);
		return;
	    }
	    if (chpnt->num!=chnum) {
	      for (chpnt=chanpnt; chpnt; chpnt=(chan*)chpnt->next)          /* channels */
		if (chpnt->num==chnum) break;
	      if (chpnt==NULL) {
	         ncfprintf (stderr,"restoremodel, can't find chan %d\n",chnum);
     	         execerror ("warning,","can't find chan");
		return;
	      }
	    }
	    chpnt->conduct = chcond;
	    chpnt->vrev    = vrev;
	    tchans++;
	    cp = scantoken (inbuf,15);			/* find 15th token */
	    if (chpnt->chtyp->hh) {
		cp=scantoken(NULL,1);
                ((hhchan *)chpnt)->m  = atof(cp);
		cp=scantoken(NULL,2);
                ((hhchan *)chpnt)->h  = atof(cp);
		cp=scantoken(NULL,1);
	    }
	    if (chpnt->cstate!=NULL) {			/* rng state */
		int rndsiz=RNDSIZ/sizeof(int);
	      cp=scantoken(NULL,1);			/* skip over "rn" */
	      for (i=0; i<rndsiz; i++) {
		*((int *)&chpnt->cstate[i*sizeof(int)]) = atoi(cp);
		cp=scantoken(NULL,1);
	      }
	    }
	    if (chpnt->nchan > 0) {			/* channel noise */
              for (i=0; i<numstate; i++) {
		if (cp==NULL) {
		    sprintf (ebuf,"restoremodel: error reading channel state: %d.\n",i);
     	            execerror ("warning,",ebuf);
		    return;
		}
                chpnt->conc[i].cval  = atof(cp);
		cp=scantoken(NULL,1);
                chpnt->conc[i].nchan = int(atof(cp));
		cp=scantoken(NULL,1);
	      }
	    }
	    else {					/* no channel noise */
              for (i=0; i<numstate; i++) {
		if (cp==NULL) {
		    sprintf (ebuf,"restoremodel: error reading channel state: %d.\n",i);
     	            execerror ("warning,",ebuf);
		    return;
		}
                chpnt->conc[i].cest  = atof(cp);
		cp=scantoken(NULL,1);
	      }
	    }
	    chpnt = (chan *)chpnt->next;
	 break;
	default: break;

  	}  /* switch (inbuf[0]) */
       } /* for all chan lines */
       break;

      case 's':						/* found "# synapses" */
       if ((spnt=synpnt) == NULL) {
    	   sprintf (ebuf,"restoremodel: no synapses, must step/run model first.\n");
     	   execerror ("warning,",ebuf);
	  break;
       }
       for (; getline(&inbuf,&inbufsiz,sfil)>0; npl++) {
        if (inbuf[0]=='#') break;  			/* stop at last synapse */
        switch (inbuf[0])  {

         case 'b':					/* voltage buffer */
            sscanf(inbuf,"b  %d %lg %lg %lg %d %d",&snum,&gain,&offset,&v,&delay,&nf);
	    if (spnt==NULL || spnt->ctype!=BUF || spnt->num!=snum) {
               for (spnt=synpnt; spnt; spnt=(synap *)spnt->next) 
		if (spnt->ctype==BUF && spnt->num==snum) break;
	      if (spnt==NULL) 
		ncfprintf (stderr,"restoremodel, can't find buf %d\n",snum);
	    } 
	    ((dbuf *)spnt)->gain   = gain;
	    ((dbuf *)spnt)->offset = offset;
	    ((dbuf *)spnt)->v      = v;

		/* look for delay buf */

            if (delay>0) {
	      sscanf(inbuf," %ld",&bpnt); 
	      ((dbuf *)spnt)->delpnt  = (double *)bpnt;
	      ((dbuf *)spnt)->delay  = delay;
	      cp = scantoken (inbuf,9);			/* find 9th token */
	      for (i=0; i<delay; i++) { 
	        ((dbuf *)spnt)->delbuf[i] = atof(cp);
	        cp=scantoken(NULL,1);
	      }
	    }
		/* look for lp filter */

            if (nf >= 1) {
 	      if ((lpf=((dbuf*)spnt)->filt)!=NULL) {
	        if (nf > lpf->nfilt) {
    	  	  sprintf (ebuf,"restoremodel: synapse filter too big: %d.\n",nf);
     	          execerror ("warning,",ebuf);
		  return;
	        }
	  	for (i=0; i<nf; i++) {
		  if (cp==NULL) {
		       sprintf (ebuf,"restoremodel: error reading synapse filter: %d.\n",i);
     	               execerror ("warning,",ebuf);
		  }
                  lpf->lfilt[i]  = atof(cp);
		  cp=scantoken(NULL,1);
		}
	      }
	    }
	    tsyns++;
	    spnt = (synap *)spnt->next;			/* fast way of finding next synap */
	    if (spnt==NULL) spnt = synpnt;
	    break;

         case 'g':					/* gap junction */
            sscanf(inbuf,"g  %d %lg %lg",&snum,&scond,&gjfrac);
	    if (spnt==NULL || spnt->ctype!=GJ || spnt->num!=snum) {
               for (spnt=synpnt; spnt; spnt=(synap *)spnt->next) 
		if (spnt->ctype==GJ && spnt->num==snum) break;
	      if (spnt==NULL) 
		ncfprintf (stderr,"restoremodel, can't find gj %d\n",snum);
	    } 
	    spnt->conduct  = scond;
	    ((gj*)spnt)->n = gjfrac;
	    tsyns++;
	    spnt = (synap *)spnt->next;			/* fast way of finding next synap */
	    if (spnt==NULL) spnt = synpnt;
	    break;

         case 's':					/* synapse */
           switch (inbuf[1])  {

            case ' ':				/* s num conduct vting time vflag rrpool */
            sscanf(inbuf,"s %d %lg %lg %lg %d %lg",&snum,&scond,&vtint,&vtime,&vflag,&rrpool);
	    if (spnt==NULL || spnt->num!=snum) {
               for (spnt=synpnt; spnt; spnt=(synap *)spnt->next) 
		if (spnt->num==snum) break;
	      if (spnt==NULL) 
		ncfprintf (stderr,"restoremodel, can't find synapse %d\n",snum);
	    } 
	    spnt->conduct = scond;
	    spnt->vtint   = vtint;
	    spnt->vtime   = vtime;
	    spnt->vtime   = vflag;
	    spnt->rrpool  = rrpool;
	    tsyns++;
	    cp = scantoken (inbuf,8);			/* find 8th token */
	    if (spnt->vstate!=NULL) {			/* ves timing rng state */
		int rndsiz=RNDSIZ/sizeof(int);
	      cp=scantoken(NULL,1);			/* skip over "vt" */
	      for (i=0; i<rndsiz; i++) {
		*((int *)&spnt->vstate[i*sizeof(int)]) = atoi(cp);
		cp=scantoken(NULL,1);
	      }
	    }
	    if (spnt->gstate!=NULL) {			/* ves size rng state */
		int rndsiz=RNDSIZ/sizeof(int);
	      cp=scantoken(NULL,1);			/* skip over "vs" */
	      for (i=0; i<rndsiz; i++) {
		*((int *)&spnt->gstate[i*sizeof(int)]) = atoi(cp);
		cp=scantoken(NULL,1);
	      }
	    }
	    tspnt = spnt;
	    spnt = (synap *)spnt->next;			/* fast way of finding next synap */
	    if (spnt==NULL) spnt = synpnt;
 	    break;

            case 'f':		/* "sf" */
              switch (inbuf[2])  {
		case 'a':		/* "sfa" */
		 if      (inbuf[3]==' ') lpf=tspnt->filt1;
		 else if (inbuf[3]=='h') lpf=tspnt->filt1h;
		 break;
		case 'b':		/* "sfb" */
		 if (inbuf[3]==' ')      lpf=tspnt->filt2;
		 break;
		case 'c':		/* "sfc" */
		 if (inbuf[3]==' ')      lpf=tspnt->filt3;
		 break;
		default: break;
	      } /* switch (inbuf[2]) */
	      if (lpf==NULL) break;
              if (sscanf(inbuf,"%*s %d",&nf) == 1) {
	        if (nf > lpf->nfilt) {
    		  sprintf (ebuf,"restoremodel: synapse filter too big: %d.\n",nf);
     	          execerror ("warning,",ebuf);
		  return;
		}
	        cp = scantoken (inbuf,3);	/* find 3rd token */
		for (i=0; i<nf; i++) {
		  if (cp==NULL) {
		     sprintf (ebuf,"restoremodel: error reading synapse filter: %d.\n",i);
     	             execerror ("warning,",ebuf);
		  }
                  lpf->lfilt[i]  = atof(cp);
		  cp=scantoken(NULL,1);
		}
	      }
 	     break;
	     default: break;
	     } /* switch (inbuf[1]) */
 	    break;
	    default: break;
	  } /* switch (inbuf[0]) */
	} /* for getline() */

      break;  /* synapses */

      case 'p':					/* photoreceptors */
       if ((rpnt=recpnt)==NULL) {
    	   sprintf (ebuf,"restoremodel: no photoreceptors, must step/run model first.\n");
     	   execerror ("warning,",ebuf);
	  break;
       }
       for (; getline(&inbuf,&inbufsiz,sfil)>0; npl++) {
         if (inbuf[0]=='#') break;  			/* stop at last comps */
         switch (inbuf[0])  {

          case 'p':					/* photoreceptor */
           sscanf(inbuf,"%*s %d %d %d %d %d %d ", &phtype, &pigm, &recnm1, &recnm2, &recnm3, &recnm4);
	   if (rpnt->recnm1 != recnm1 &&
	       rpnt->recnm2 != recnm2 &&
	       rpnt->recnm3 != recnm3 &&
	       rpnt->recnm4 != recnm4 ) {
		   findphotrec(recnm1,recnm2,recnm3,recnm4,&rpnt,"restoremodel");
	   }
	   if (phtype < 2) {		// orig type photorecs
	    vpnt = &rpnt->vars;
            sscanf(inbuf,"%*s %*g %*g %*g %*g %*g %*g %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
		&rpnt->conduct,
		&vpnt->rhod,
		&vpnt->rhodk,
		&vpnt->rhod2,
		&vpnt->rstar,
		&vpnt->gpr,
		&vpnt->pde,
		&vpnt->cycg,
		&vpnt->gcyc,
		&vpnt->cond,
		&vpnt->ca,
		&vpnt->cax,
		&vpnt->cab,
		&vpnt->ksens,
		&vpnt->resp_r,
		&vpnt->resp_b,
		&vpnt->resp_e,
		&vpnt->resp_x,
		&vpnt->resp_c,
		&vpnt->resp_is,
		&vpnt->atten_i);


	    } else if (phtype==PINVG) {				// invergo type mouse rod photorec

		    int npp;			// number of parameters

	    rpnti = (photreci*)rpnt;
	    vpnti = &rpnti->vars;

            npp = sscanf(inbuf,"%*s %*g %*g %*g %*g %*g %*g %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
	   
	   
		&rpnti->conduct,
		&vpnti->Rh,                      // unphosphorylated rhodopsin
		&vpnti->R_Gt,                    //
		&vpnti->Gt,                      //
		&vpnti->Rn[0],                   // phosphorylated rhodopsin
		&vpnti->Rn[1],                   // 
		&vpnti->Rn[2],                   // 
		&vpnti->Rn[3],                   // 
		&vpnti->Rn[4],                   // 
		&vpnti->Rn[5],                   // 
		&vpnti->Rn[6],                   // 
		&vpnti->RK,                      // rhodopsin kinase
		&vpnti->Rn_RKpre[0],              // 
		&vpnti->Rn_RKpre[1],
		&vpnti->Rn_RKpre[2],
		&vpnti->Rn_RKpre[3],
		&vpnti->Rn_RKpre[4],
		&vpnti->Rn_RKpre[5],
		&vpnti->Rn_RKpre[6],
		&vpnti->Rn_RKpost[0],
		&vpnti->Rn_RKpost[1],
		&vpnti->Rn_RKpost[2],
		&vpnti->Rn_RKpost[3],
		&vpnti->Rn_RKpost[4],
		&vpnti->Rn_RKpost[5],
		&vpnti->Rn_RKpost[6],

		&vpnti->Arr,                     // arrestin 
		&vpnti->Rn_Arr[0],               // rhodopsin-arrestin complex
		&vpnti->Rn_Arr[1],               //
		&vpnti->Rn_Arr[2],               //
		&vpnti->Rn_Arr[3],               //
		&vpnti->Rn_Arr[4],               //
		&vpnti->Rn_Arr[5],               //
		&vpnti->Rn_Arr[6],               //
		&vpnti->Ops,                     // ligand-free receptor 
		&vpnti->Arr_di,                  // arrestin dimer
		&vpnti->Arr_tetra,               // arrestin tetramer
		&vpnti->Ops_Gt,                  // receptor Gt complex
		&vpnti->Ops_G,                   // receptor G complex
		&vpnti->Ops_Ggtp,                // receptor G complex
		&vpnti->Ggtp,                    // Ggtp 
		&vpnti->Rn_Gt[0],                // receptor G complex
		&vpnti->Rn_Gt[1],                // 
		&vpnti->Rn_Gt[2],                // 
		&vpnti->Rn_Gt[3],                // 
		&vpnti->Rn_Gt[4],                // 
		&vpnti->Rn_Gt[5],                // 
		&vpnti->Rn_Gt[6],                // 
		&vpnti->Rn_G[0],                 // receptor G complex
		&vpnti->Rn_G[1],                 // 
		&vpnti->Rn_G[2],                 // 
		&vpnti->Rn_G[3],                 // 
		&vpnti->Rn_G[4],                 // 
		&vpnti->Rn_G[5],                 // 
		&vpnti->Rn_G[6],                // 

		&vpnti->Rn_Ggtp[0],              // receptor G complex
		&vpnti->Rn_Ggtp[1],              // 
		&vpnti->Rn_Ggtp[2],              // 
		&vpnti->Rn_Ggtp[3],              // 
		&vpnti->Rn_Ggtp[4],              // 
		&vpnti->Rn_Ggtp[5],              // 
		&vpnti->Rn_Ggtp[6],              // 
		&vpnti->Gagtp,                   // 
		&vpnti->Gbg,                     // 
		&vpnti->PDE,                     // 
		&vpnti->PDE_Gagtp,               // PDE-Gagtp complex
		&vpnti->PDEa_Gagtp,              // activated PDE-Gagtp complex
		&vpnti->Gagtp_PDEa_Gagtp,        // activated PDE-Gagtp complex
		&vpnti->Gagtp_aPDEa_Gagtp,       // activated PDE-Gagtp complex
		&vpnti->RGS,                     // 
		&vpnti->RGS_PDEa_Gagtp,          // activated RGS-PDE-Gagtp complex
		&vpnti->RGS_Gagtp_aPDEa_Gagtp,   // activated RGS-PDE-Gagtp complex
		&vpnti->Gagdp,                   // Gagdp 
		&vpnti->Rect,                    // 
		&vpnti->Recr_Ca,                 // 
		&vpnti->Cafree,                  // Ca free
		&vpnti->Recr_Ca_RK,              // 
		&vpnti->Cabuff,                  // Ca buffered
		&vpnti->cGMP,                    // cyclic GMP
		&vpnti->cond);                   // dark conductance normalized 0 to 1
		if (npp < 80)  ncfprintf (stderr,"restoremodel, error reading photoreceptor params %d\n",npp);

	    } else if (phtype==ChR2) {				// channel rhodopsin from Williams et al (2013)

		    int npp;			// number of parameters
		    photrecc *rpntc;
		    phvarsc *vpntc;

	    rpntc = (photrecc*)rpnt;
	    vpntc = &rpntc->vars;

            npp = sscanf(inbuf,"%*s %*g %*g %*g %*g %*g %*g %lg %lg %lg %lg %lg %lg %lg",
		&rpnt->conduct,
		&vpntc->SC1,                     // state C1 
		&vpntc->SO1,                     // state O1 
		&vpntc->SO2,                     // state O2 
		&vpntc->SC2,                     // state C2 
		&vpntc->p,                       // nonlinear intensity factor
		&vpntc->cond);                   // dark conductance normalized 0 to 1
	    }

  	    tphotrecs++;
          break;  /* one photoreceptor line */
	  default: break;

         }  /* switch (inbuf[0]); */
       }  /* for all photoreceptor lines */
       break;
       default: break;
     }  /* switch (inbuf[2])  looking for # comps, # syns, # photrecs, # kchans, etc */
   } /* if (inbuf[0])=='#') */
  }
  // if (tcomps!=cumcomp || 
  //     tchans!=cumchan || 
  //      tsyns!=cumsynap+cumgj+cumdbuf || 
  //  tphotrecs!=cumphotrec) {
  //    fprintf (stderr,"# restore model: comps=%d chans=%d syns=%d photrecs=%d\n",
  //		tcomps, tchans, tsyns, tphotrecs);
  // }
 fclose(sfil);
}

