/* module plot_funcs.cc */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"
#include "ncio.h"

extern int chanelem[NCELTYPES][NCHANS];	    /* channel element enames (soma only) */

double spike_plot(double nodenum, double time);
double freq_color (int nplot, double xval, double yval);

#define CSIZ 50

/*------------------------ plot -----------------------------*/

double na_inact_plot(double ct, double cn, double n, double xtime)

/* Function to be called by plot(), required to have 2 arguments. */

{
    elem *epnt;
    node *npnt;
    conlst *clst;
    chan *chpnt;

  if ((npnt=nde(ct,cn,n))==NULL) return 0;
  clst = npnt->comptr->clst;
  chpnt = findchan(clst,(int)NA,(int)2);
  return (record_chan(chpnt,G,7) + record_chan(chpnt,G,8) + record_chan(chpnt,G,9)); 
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double nmda_act_plot(double ct, double cn, double n, double xtime)

/* Return sum of open and mg-flicker states from an NMDA type 2 channel. */
/* Function to be called by plot(), required to have 2 arguments. */

{
    double cond;
    elem *epnt;
    node *npnt;
    conlst *clst;
    chan *chpnt;

  if ((npnt=nde(ct,cn,n))==NULL) return 0;
  clst = npnt->comptr->clst;
  chpnt = findchan(clst,(int)NMDA,(int)2);

  cond = ((synapse *)epnt)->maxcond *  (record_chan(chpnt,G,5) + record_chan(chpnt,G,7) + 
		  record_chan(chpnt,G,9) + record_chan(chpnt,G,12)); 
  return (cond);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double nmda_act_plot(double chan_elem, double xtime)

/* Return sum of open and mg-flicker states from an NMDA type 2 channel. */
/* Function to be called by plot(), required to have 2 arguments. */

{
    elem *epnt;
    double cond;

  epnt = get_elempnt((int)chan_elem);
  cond = ((synapse *)epnt)->maxcond *  (record_chan(epnt,G,5) + record_chan(epnt,G,7) +
                   record_chan(epnt,G,9) + record_chan(epnt,G,12));
  return (cond);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double nmda_act_plot(elem *epnt, double xtime)

/* Return sum of open and mg-flicker states from an NMDA type 2 channel. */
/* Function to be called by plot(), required to have 2 arguments. */

{
    double cond;

  if (epnt==NULL) return 0;
  cond = ((synapse *)epnt)->maxcond *  (record_chan(epnt,G,5) + record_chan(epnt,G,7) +
                   record_chan(epnt,G,9) + record_chan(epnt,G,12));
  return (cond);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int set_plotnum(int ct)

{
   int plotnum;

 if(ct==dsgc || ct==gca || ct==gcb || ct==gcaoff || ct==gcboff) /* set where to put the plot */
   plotnum = 20;
 else if(ct==aii)
   plotnum = 30;
 else if(ct==sbac)
   plotnum = 32;

 else if(ct==a17)
   plotnum = 40;
 else if(ct==am)
   plotnum = 42;
 else if(ct==amh)
   plotnum = 44;

 else if(ct==dbp1)
   plotnum = 60;
 else if(ct==dbp2)
   plotnum = 61;
 else if(ct==hbp1)
   plotnum = 62;
 else if(ct==hbp2)
   plotnum = 63;
 else if(ct==rbp)
   plotnum = 64;

 else if (ct==ha || ct ==hb)
   plotnum = 80;
 else if (ct==xcone || ct ==xrod)
   plotnum = 90;
 else plotnum = 10;
 return plotnum;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_v_nod(node *npnt, double vmin, double vmax,int pcolor, const char *label, int plotnum, double psize)

{
    int c, ct, cn, n;
    char *name;
    char nbuf[CSIZ];

  if (npnt==NULL) return;

  ct = npnt->nodenm1;
  cn = npnt->nodenm2;
  n  = npnt->nodenm3;

  if (script==0) {
    if (plotnum<0) {
      plotnum = set_plotnum(ct);
    }
    
    if (pcolor==-1)	/* if user doesn't specify a color */
      c=plotnum;  	/* use graph num to set color */
    else
      c=pcolor;

    if (label==NULL || streq(label,"")) {
      name = nbuf;
      if (n==soma) sprintf(name,"V%s_%d_soma",cname[ct],cn);
      else         sprintf(name,"V%s_%d_%d",cname[ct],cn,n);
    }
    else name = (char *)label;

    if (psize < 0) {
      plot (V,npnt, vmax, vmin); plot_param(name, c, plotnum);
    } else {
      plot (V,npnt, vmax, vmin); plot_param(name, c, plotnum, psize);
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_v_nod(int ct,int cn,int n, double vmin, double vmax,int pcolor, const char *label, int plotnum, double psize)

{
   plot_v_nod(nde(ct,cn,n), vmin, vmax, pcolor, label, plotnum, psize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_vm_nod(node *npnt, double vmin, double vmax,int pcolor, const char *label, int plotnum, double psize)

{
    int c, ct, cn, n;
    char *name;
    char nbuf[CSIZ];

  if (npnt==NULL) return;

  ct = npnt->nodenm1;
  cn = npnt->nodenm2;
  n  = npnt->nodenm3;

  if (script==0) {
    if (plotnum<0) {
      plotnum = set_plotnum(ct);
    }
    
    if (pcolor==-1)	/* if user doesn't specify a color */
      c=plotnum;  	/* use graph num to set color */
    else
      c=pcolor;

    if (label==NULL || streq(label,"")) {
      name = nbuf;
      if (n==soma) sprintf(name,"VM%s_%d_soma",cname[ct],cn);
      else         sprintf(name,"VM%s_%d_%d",cname[ct],cn,n);
    }
    else name = (char *)label;

    if (psize < 0) {
      plot (VM,npnt, vmax, vmin); plot_param(name, c, plotnum);
    } else {
      plot (VM,npnt, vmax, vmin); plot_param(name, c, plotnum, psize);
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_vm_nod(int ct,int cn,int n, double vmin, double vmax,int pcolor, const char *label, int plotnum, double psize)

{
   plot_vm_nod(nde(ct,cn,n), vmin, vmax, pcolor, label, plotnum, psize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_v_lnod(int ct,int cn,int labl, double vmin, double vmax,int pcolor, const char *label, int plotnum, double psize)

/* plot V from node which has label in 8th column in morpology file */

{
    int c;
    char *name;
    char nbuf[CSIZ];

  if (labl < 0) return;
  if (nde(ct,cn,dendn_node(ct,labl))==NULL) return;

  if (script==0) {
    if (plotnum<0) {
      plotnum = set_plotnum(ct);
    }
    
    if (pcolor==-1)	/* if user doesn't specify a color */
      c=plotnum;  	/* use graph num to set color */
    else
      c=pcolor;

    name = nbuf;
    if (label==NULL || streq(label,"")) {
      if (labl==soma) sprintf(name,"V%s_%d_soma",cname[ct],cn);
      else            sprintf(name,"V%s_%d_%d",cname[ct],cn,labl);
    }
    else  if (labl==soma) sprintf(name,"V%s_%d_soma",label,cn);
    else                  sprintf(name,"V%s_%d_%d",label,cn,labl);

    if (psize < 0) {
      plot (V,nde(ct,cn,dendn_node(ct,labl)), vmax, vmin); plot_param(name, c, plotnum);
    } else {
      plot (V,nde(ct,cn,dendn_node(ct,labl)), vmax, vmin); plot_param(name, c, plotnum, psize);
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_vm_lnod(int ct,int cn,int labl, double vmin, double vmax,int pcolor, const char *label, int plotnum, double psize)

/* plot VM from node which has label in 8th column in morpology file */

{
    int c;
    char *name;
    char nbuf[CSIZ];

  if (labl < 0) return;
  if (nde(ct,cn,dendn_node(ct,labl))==NULL) return;

  if (script==0) {
    if (plotnum<0) {
      plotnum = set_plotnum(ct);
    }
    
    if (pcolor==-1)	/* if user doesn't specify a color */
      c=plotnum;  	/* use graph num to set color */
    else
      c=pcolor;

    name = nbuf;
    if (label==NULL || streq(label,"")) {
      if (labl==soma) sprintf(name,"VM%s_%d_soma",cname[ct],cn);
      else            sprintf(name,"VM%s_%d_%d",cname[ct],cn,labl);
    }
    else  if (labl==soma) sprintf(name,"VM%s_%d_soma",label,cn);
    else                  sprintf(name,"VM%s_%d_%d",label,cn,labl);

    if (psize < 0) {
      plot (VM,nde(ct,cn,dendn_node(ct,labl)), vmax, vmin); plot_param(name, c, plotnum);
    } else {
      plot (VM,nde(ct,cn,dendn_node(ct,labl)), vmax, vmin); plot_param(name, c, plotnum, psize);
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_l_nod(int ct,int cn,int n, double lmin, double lmax,int pcolor, const char *label, int plotnum, double psize)

{
    int c;
    char *name;
    char nbuf[CSIZ];

  if (n <= -1) return;
  if (nde(ct,cn,n)==NULL) return;

  if (script==0) {
    if (plotnum<0) {
      plotnum = set_plotnum(ct);
    }
    
    if (pcolor==-1)	/* if user doesn't specify a color */
      c=plotnum;  	/* use graph num to set color */
    else
      c=pcolor;

    if (label==NULL || streq(label,"")) {
      name = nbuf;
      if (n==soma) sprintf(name,"L%s_%d_soma",cname[ct],cn);
      else         sprintf(name,"L%s_%d_%d",cname[ct],cn,n);
    }
    else name = (char *)label;

    if (psize < 0) {
      plot (L,nde(ct,cn,n), lmax, lmin); plot_param(name, c, plotnum);
    } else {
      plot (L,nde(ct,cn,n), lmax, lmin); plot_param(name, c, plotnum, psize);
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_ca_nod(int ct, int cn, int n, int shell, double maxca, int pcolor, const char *label, int plotnum, double psize)
{
    int c;
    double fmin;
    char *name;
    char nbuf[CSIZ];
    node *npnt;

  if (n == -1) return;
  if ((npnt=nde(ct,cn,n))==NULL) return;
//  if (get_nfield(npnt,CACOMP)) {

    if (plotnum<0) {
      plotnum = set_plotnum(ct) + 2;
    }
    
    if (pcolor==-1)	/* if user doesn't specify a color */
      c=plotnum;	/* use graph num to set color */
    else
      c=pcolor;

    if (label==NULL || streq(label,"")) {
      name = nbuf;
      if (n==soma) sprintf(name,"Ca_%s_%d_soma",cname[ct],cn);
      else         sprintf(name,"Ca_%s_%d_%d",cname[ct],cn,n);
    }
    else name = (char *)label;

    if (psize < 0) {
      plot (CA, shell, npnt, maxca, fmin=0); plot_param(name,c,plotnum);
    } else {
      plot (CA, shell, npnt, maxca, fmin=0); plot_param(name,c,plotnum,psize);
    }
//  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_ca_nod(node *npnt, int shell, double maxca, int pcolor, const char *label, int plotnum, double psize)
{
    int c, ct, cn, n;
    double fmin;
    char *name;
    char nbuf[CSIZ];

  if (npnt==NULL) return;
  if (get_nfield(npnt,CACOMP)) {

    ct = npnt->nodenm1;
    cn = npnt->nodenm2;
    n  = npnt->nodenm3;

    if (plotnum<0) {
      plotnum = set_plotnum(ct) + 2;
    }
    
    if (pcolor==-1)	/* if user doesn't specify a color */
      c=plotnum;	/* use graph num to set color */
    else
      c=pcolor;

    if (label==NULL || streq(label,"")) {
      name = nbuf;
      if (n==soma) sprintf(name,"Ca_%s_%d_soma",cname[ct],cn);
      else         sprintf(name,"Ca_%s_%d_%d",cname[ct],cn,n);
    }
    else name = (char *)label;

    if (psize < 0) {
      plot (CA, shell, npnt, maxca, fmin=0); plot_param(name,c,plotnum);
    } else {
      plot (CA, shell, npnt, maxca, fmin=0); plot_param(name,c,plotnum,psize);
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_ca_nod(int ct, int cn, int n, double maxca, int pcolor, const char *label, int plotnum, double psize)
{
    plot_ca_nod(ct, cn, n, 1, maxca, pcolor, label, plotnum, psize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_ca_syn(synapse *s, int sh, double maxca, int pcolor, const char *label, int plotnum, double psize)
{
    node *npnt;

    if (s==NULL) return;
    if ((npnt=s->nodp1)==NULL) return;
    plot_ca_nod(npnt, sh, maxca, pcolor, label, plotnum, psize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_ca_syn(synapse *s, double maxca, int pcolor, const char *label, int plotnum, double psize)

{
    node *npnt;

    plot_ca_syn(s, 1, maxca, pcolor, label, plotnum, psize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_cabufb_nod(int ct, int cn, int n, int sh, double maxca, int pcolor, const char *label, int plotnum, double psize)
{
    int c;
    double fmin;
    char *name;
    char nbuf[CSIZ];
    node *npnt;

  if (n <= -1) return;
  if ((npnt=nde(ct,cn,n))==NULL) return;
  if (get_nfield(npnt,CACOMP)) {

    if (plotnum<0) {
      plotnum = set_plotnum(ct) + 2;
    }
    
    if (pcolor==-1)	/* if user doesn't specify a color */
      c=plotnum;	/* use graph num to set color */
    else
      c=pcolor;

    if (label==NULL || streq(label,"")) {
      name = nbuf;
      if (n==soma) sprintf(name,"Cabufb_%s_%d_soma",cname[ct],cn);
      else         sprintf(name,"Cabufb_%s_%d_%d",cname[ct],cn,n);
    }
    else name = (char *)label;

    if (psize < 0) {
      plot (CABUFB, sh, npnt, maxca, fmin=0); plot_param(name,c,plotnum);
    } else {
      plot (CABUFB, sh, npnt, maxca, fmin=0); plot_param(name,c,plotnum,psize);
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_cabufb_nod(int ct, int cn, int n, double maxca, int pcolor, const char *label, int plotnum, double psize)
{
plot_cabufb_nod(ct, cn, n, 1, maxca, pcolor, label, plotnum, psize);
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_cas_nod(int ct, int cn, int n, int sh, double maxca, int pcolor, const char *label, int plotnum, double psize)
{
    int c;
    double fmin;
    char *name;
    char nbuf[CSIZ];
    node *npnt;

  if (n <= -1) return;
  if ((npnt=nde(ct,cn,n))==NULL) return;
  if (get_nfield(npnt,CACOMP)) {

    if (plotnum<0) {
      plotnum = set_plotnum(ct) + 2;
    }
    
    if (pcolor==-1)	/* if user doesn't specify a color */
      c=plotnum;	/* use graph num to set color */
    else
      c=pcolor;

    if (label==NULL || streq(label,"")) {
      name = nbuf;
      if (n==soma) sprintf(name,"Cas_%s_%d_soma",cname[ct],cn);
      else         sprintf(name,"Cas_%s_%d_%d",cname[ct],cn,n);
    }
    else name = (char *)label;

    if (psize < 0) {
      plot (CAS, sh, npnt, maxca, fmin=0); plot_param(name,c,plotnum);
    } else {
      plot (CAS, sh, npnt, maxca, fmin=0); plot_param(name,c,plotnum,psize);
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_cas_nod(int ct, int cn, int n, double maxca, int pcolor, const char *label, int plotnum, double psize)
{
plot_cas_nod(ct, cn, n, 1, maxca, pcolor, label, plotnum, psize);
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_ph_nod(node *npnt, double maxnt, double minnt, int pcolor, const char *label, int plotnum, double psize)

{
    int c;
    int ct, cn, n;
    char *name;
    char nbuf[CSIZ];

  if (npnt==NULL) return;
  
  ct = npnt->nodenm1;
  cn = npnt->nodenm2;
  n  = npnt->nodenm3;

  if (plotnum<0) {
      plotnum = set_plotnum(ct) + 2;
  }
  if (pcolor==-1)	/* if user doesn't specify a color */
      c=plotnum;	/* use graph num to set color */
  else
      c=pcolor;

  if (label==NULL || streq(label,"")) {
      name = nbuf;
      if (n==soma) sprintf(name,"pH_%s_%d_soma",cname[ct],cn);
      else         sprintf(name,"pH_%s_%d_%d",cname[ct],cn,n);
  }
  else name = (char *)label;

  if (psize < 0) {
    plot (PH,npnt, maxnt, minnt); plot_param(name, c, plotnum);
  } else {
    plot (PH,npnt, maxnt, minnt); plot_param(name, c, plotnum, psize);
  } 
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_ph_nod(int ct, int cn, int n, double maxnt, double minnt, int pcolor, 
				const char *label, int plotnum, double psize)

{
     node *npnt;

  if ((npnt=nde(ct,cn,n))==NULL) return;
  plot_ph_nod(npnt, maxnt, minnt, pcolor, label, plotnum, psize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_spike_rate(int ct, int cn, int n, int pcolor, const char *label, int plotnum, double psize)
{
    int c;
    double fmax, fmin;
    char *name;
    char nbuf[CSIZ];
    node *npnt;

  if (n == -1) return;
  if ((npnt=nde(ct,cn,n))==NULL) return;
  
   if (plotnum < 0) {
      plotnum = set_plotnum(ct) + 5;
   };
    
   if (pcolor==-1) { c=plotnum; }
   else { c=pcolor; };

    if (label==NULL || streq(label,"")) {
      name = nbuf;
      if (n==soma) sprintf(name,"Fi[%s_%d_soma]",cname[ct],cn);
      else         sprintf(name,"Fi[%s_%d_%d]",cname[ct],cn,n);
    }
    else name = (char *)label;

   plot_func((double(*)(double, double))spike_plot,n,fmax=500,fmin=0);
		plot_param(name, c, plotnum, psize);
   plot_vpen(freq_color);
   graph_csiz (0.01);
   graph_cchar ('o');
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_na_inact(int ct, int cn, int n, int pcolor, const char *label, int plotnum, int psize)

/* Plot Na chan inactivated states. The plot is accomplished by passing element number
   of a dummy channel, which is used to find the correct channel to record at run time. */
 
{
    int c, chan_elem;
    double fmax, fmin;
    char *name;
    char nbuf[CSIZ];
    node *npnt;
    elem *epnt;


  if (n == -1) return;
  if ((npnt=nde(ct,cn,n))==NULL) return;
  
  if (plotnum < 0) {
      plotnum = set_plotnum(ct) + 6;
  }
    
  if (pcolor==-1) { c=plotnum; }
  else { c=pcolor; };

   if (label==NULL || streq(label,"")) {
     name = nbuf;
     if (n==soma) sprintf(name,"Na[%s_%d_soma]",cname[ct],cn);
     else         sprintf(name,"Na[%s_%d_%d]",cname[ct],cn,n);
   }
   else name = (char *)label;
   
   plot_func(na_inact_plot,(double)ct,(double)cn,(double)n,fmax=1,fmin=0); plot_param(name, c, plotnum, psize);
} 

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate(synapse *s, int prate, double rmin, double rmax, int pves, double fmin, double fmax, int pcond, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize)

{
    int t, ct, cn, nod, col;
    elem *cpnt;
    char name[CSIZ], vesnam[CSIZ], condnam[CSIZ];
    attrib *apnt;
    synapse *tc;

    /* look for cGMP channel on postsynaptic side */
 
    if (s==NULL) return;
    t = s->elnum;
    tc = NULL;
    ct = s->node1a;
    cn = s->node1b;
    nod = s->node1c;

    if (plotnum==-1) plotnum = set_plotnum(ct) - 1;

    if (pcolor==-1)		/* if user doesn't specify a color */
      col=plotnum;	  	/* use graph num to set color */
    else
      col=pcolor;

    for (cpnt=elempnt; cpnt=foreach(cpnt,CHAN,s->node2a,s->node2b,-1,-1,NULL,NULL,NULL,NULL);
		       cpnt=cpnt->next) {
        apnt = get_chanattr(cpnt);
	 if (apnt->ctype==CGMP) {
            tc = s;
            break;
        }
    }
    //printf("%d %d %d\n", s->node1a, s->node1b, s->node1c);

    if (!t) return;
    if (plname==NULL || streq(plname,"")) {
        sprintf(name,   "rate %s_%d_%d",cname[ct],cn,nod);
        sprintf(vesnam, "ves  %s_%d_%d",cname[ct],cn,nod);
        sprintf(condnam,"G%s_%d_%d",cname[ct],cn,nod);
    } else {
        sprintf(name,   "rate %s",plname);
        sprintf(vesnam, "ves  %s",plname);
        sprintf(condnam,"cond %s",plname);
    }
    if (prate && pves) sprintf(vesnam,"");

    if (prate) { plot(FA9,s->elnum,rmax,rmin); plot_param(name,col,plotnum,plsize); }
    if (pves)  { plot(FB4,s->elnum,fmax,fmin); plot_param(vesnam,col,plotnum,plsize); }
    if (pcond) {
       if (tc==NULL) {
	  plot(G,s->elnum,   cmax,cmin); plot_param(condnam,col,plotnum,plsize);
       }
       else {
          plot(G,tc->elnum,cmax,cmin); plot_param(condnam,col,plotnum,plsize); 
      }
    }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synratep(synapse *s, int prate, double rmin, double rmax, int pves, double fmin, double fmax, int pcond, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize)

/* print name and number of postsynaptic cell */

{
    int t, ct, cn, nod, col, ct2, cn2, nod2;
    elem *cpnt;
    char name[CSIZ], vesnam[CSIZ], condnam[CSIZ];
    attrib *apnt;
    synapse *tc;

    /* look for cGMP channel on postsynaptic side */
 
    if (s==NULL) return;
    t = s->elnum;
    tc = NULL;
    ct = s->node1a;
    cn = s->node1b;
    nod = s->node1c;
    ct2 = s->node2a;
    cn2 = s->node2b;
    nod2 = s->node2c;

    if (plotnum==-1) plotnum = set_plotnum(ct) - 1;

    if (pcolor==-1)		/* if user doesn't specify a color */
      col=plotnum;	  	/* use graph num to set color */
    else
      col=pcolor;

    for (cpnt=elempnt; cpnt=foreach(cpnt,CHAN,s->node2a,s->node2b,-1,-1,NULL,NULL,NULL,NULL);
		       cpnt=cpnt->next) {
        apnt = get_chanattr(cpnt);
	 if (apnt->ctype==CGMP) {
            tc = s;
            break;
        }
    }
    //printf("%d %d %d\n", s->node1a, s->node1b, s->node1c);

    if (!t) return;
    if (plname==NULL || streq(plname,"")) {
        sprintf(name,   "rate %s_%d %s_%d",cname[ct],nod,cname[ct2],nod2);
        sprintf(vesnam, "ves %s_%d %s_%d",cname[ct],nod,cname[ct2],nod2);
        sprintf(condnam,"G%s_%d %s_%d",   cname[ct],nod,cname[ct2],nod2);
    } else {
        sprintf(name,   "rate %s",plname);
        sprintf(vesnam, "ves  %s",plname);
        sprintf(condnam,"cond %s",plname);
    }
    if (prate && pves) sprintf(vesnam,"");

    if (prate) { plot(FA9,s->elnum,rmax,rmin); plot_param(name,col,plotnum,plsize); }
    if (pves)  { plot(FB4,s->elnum,fmax,fmin); plot_param(vesnam,col,plotnum,plsize); }
    if (pcond) {
       if (tc==NULL) {
	  plot(G,s->elnum,   cmax,cmin); plot_param(condnam,col,plotnum,plsize);
       }
       else {
          plot(G,tc->elnum,cmax,cmin); plot_param(condnam,col,plotnum,plsize); 
      }
    }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate(synapse *s, double rmin, double rmax, int pcolor, int plotnum, const char *plname, double plsize)

{
     int prate,pves,pcond;
     double fmin, fmax;

   plot_synrate(s, prate=1, rmin, rmax, pves=0, fmin=0, fmax=0, pcond=0, 0, 0, pcolor, plotnum, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synves(synapse *s, double fmin, double fmax, int pcolor, int plotnum, const char *plname, double plsize)

{
     int prate,pves,pcond;
     double rmin, rmax;
   plot_synrate(s, prate=0, rmin=0, rmax=0, pves=1, fmin, fmax, pcond=0, 0, 0, pcolor, plotnum, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_syncond(synapse *s, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize)

{
     int prate,pves,pcond;
   plot_synrate(s, prate=0, 0,0, pves=0, 0, 0, pcond=1, cmin, cmax, pcolor, plotnum, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_syncondp(synapse *s, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize)

{
     int prate,pves,pcond;
   plot_synratep(s, prate=0, 0,0, pves=0, 0, 0, pcond=1, cmin, cmax, pcolor, plotnum, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate(int ct, int cn, int nod, int prate, double rmin, double rmax, int pves, double fmin, double fmax, int pcond, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize)

/* Display synaptic output rate for a cell */
/* Use for different color from each cell */

{
    int col,t;
    char name[CSIZ], vesnam[CSIZ], condnam[CSIZ];
    elem *cpnt, *epnt, *tc;
    attrib *apnt;
    synapse *s;

  if (script==0) {
      if (plotnum==-1) plotnum = set_plotnum(ct) - 1;

    if (pcolor==-1)		/* if user doesn't specify a color */
      col=plotnum;	  	/* use graph num to set color */
    else
      col=pcolor;

    s = NULL;
    t = 0;
    for (epnt=elempnt; epnt=foreach(epnt,SYNAPSE,ct,cn,nod,-1,NULL,NULL,NULL,NULL);
		       epnt=epnt->next) {
       s = (synapse *)epnt;
       if (s->node1a==ct && s->node1b==cn) {
            t = s->elnum;
            break;
       }
    }
    if (s==NULL) return;

    //print element s->node1a, element s->node1b, element s->node1c, 
    //      element s->node2a, element s->node2b, element s->node2c;

    /* look for cGMP channel on postsynaptic side */

    tc = NULL;
    for (cpnt=elempnt; cpnt=foreach(cpnt,CHAN,s->node2a,s->node2b,-1,-1,NULL,NULL,NULL,NULL);
		       cpnt=cpnt->next) {
        apnt = get_chanattr(cpnt);
	 if (apnt->ctype==CGMP) {
            tc = epnt;
            break;
        }
    }
    //printf("%d %d %d\n", s->node1a, s->node1b, s->node1c);

    if (!t) return;
    if (plname==NULL || streq(plname,"")) {
        if (nod >= 0) {
	   sprintf(name,   "rate %s_%d_%d",cname[ct],cn,nod);
           sprintf(vesnam, "ves  %s_%d_%d",cname[ct],cn,nod);
           sprintf(condnam,"G%s_%d_%d",cname[ct],cn,nod);
	} else {
           sprintf(name,   "rate %s_%d",cname[ct],cn);
           sprintf(vesnam, "ves  %s_%d",cname[ct],cn);
	   sprintf(condnam,"G%s_%d",cname[ct],cn);
	}
    } else {
        sprintf(name,   "rate %s",plname);
        sprintf(vesnam, "ves  %s",plname);
        sprintf(condnam,"cond %s",plname);
    }
    if (prate && pves) sprintf(vesnam,"");

    if (prate) { plot(FA9,s->elnum,rmax,rmin); plot_param(name,col,plotnum,plsize); }
    if (pves)  { plot(FB4,s->elnum,fmax,fmin); plot_param(vesnam,col,plotnum,plsize); }
    if (pcond) {
       if (!tc) {
	  plot(G,s->elnum,   cmax,cmin); plot_param(condnam,col,plotnum,plsize);
       }
       else {
          plot(G,epnt->elnum,cmax,cmin); plot_param(condnam,col,plotnum,plsize); 
      }
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate(int ct, int cn, int prate, double rmin, double rmax, int pves, double fmin, double fmax, int pcond, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize)

{
   int nod;

  plot_synrate(ct, cn, nod=-1, prate, rmin, rmax, pves, fmin, fmax, pcond, cmin, cmax, pcolor, plotnum, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate(int ct, int cn, int nod, double rmin, double rmax, int pcolor, int plotnum, const char *plname, double plsize)

{
     int prate,pves,pcond;
     double fmin, fmax;

   plot_synrate(ct, cn, nod, prate=1, rmin, rmax, pves=0, fmin=0, fmax=0, pcond=0, 0, 0, pcolor, plotnum, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate(int ct, int cn, int prate, double rmin, double rmax, int pves, double fmin, double fmax, int pcolor, int plotnum, const char *plname, double plsize)
{
   int pcond;

  plot_synrate(ct, cn, prate, rmin, rmax, pves, fmin, fmax, pcond=0, 0, 0, pcolor, plotnum, plname, plsize);
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate(int ct, int cn, double rmin, double rmax, int pcolor, int plotnum, const char *plname, double plsize)

{
     int prate,pves,pcond;
     double fmin, fmax;

   plot_synrate(ct, cn, prate=1, rmin, rmax, pves=0, fmin=0, fmax=0, pcond=0, 0, 0, pcolor, plotnum, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synves(int ct, int cn, int nod, double fmin, double fmax, int pcolor, int plotnum, const char *plname, double plsize)

{
     int prate,pves,pcond;
     double rmin, rmax;
   plot_synrate(ct, cn, nod, prate=0, rmin=0, rmax=0, pves=1, fmin, fmax, pcond=0, 0, 0, pcolor, plotnum, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synves(int ct, int cn, double fmin, double fmax, int pcolor, int plotnum, const char *plname, double plsize)

{
     int prate,pves,pcond;
     double rmin, rmax;
   plot_synrate(ct, cn, prate=0, rmin=0, rmax=0, pves=1, fmin, fmax, pcond=0, 0, 0, pcolor, plotnum, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_syncond(int ct, int cn, int nod, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize)

{
     int prate,pves,pcond;
   plot_synrate(ct, cn, nod, prate=0, 0,0, pves=0, 0, 0, pcond=1, cmin, cmax, pcolor, plotnum, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_syncond(int ct, int cn, double cmin, double cmax, int pcolor, int plotnum, const char *plname, double plsize)

{
     int prate,pves,pcond;
   plot_synrate(ct, cn, prate=0, 0,0, pves=0, 0, 0, pcond=1, cmin, cmax, pcolor, plotnum, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, double rmin, double rmax, double fmax, int pcolor, const char *plname, double plsize)

/* Display synaptic output rate for a cell */

{
    int c1,c2,t;
    int plotnum;
    double fmin,fmaxx;
    char name[CSIZ], vesnam[CSIZ], condnam[CSIZ];
    elem *cpnt, *epnt, *s, *tc;
    attrib *apnt;

  if (script==0) {
      plotnum = set_plotnum(ct) - 1;

    if (pcolor==-1)		/* if user doesn't specify a color */
      c2=plotnum;	  	/* use graph num to set color */
    else
      c2=pcolor;
    c1 = magenta;

    s = NULL;
    t = 0;
    for (epnt=elempnt; epnt=foreach(epnt,SYNAPSE,ct,cn,-1,-1,NULL,NULL,NULL,NULL);
		       epnt=epnt->next) {
       s = epnt;
       if (s->node1a==ct && s->node1b==cn) {
            t = s->elnum;
            break;
       }
    }
    if (s==NULL) return;

    //print element s->node1a, element s->node1b, element s->node1c, 
    //      element s->node2a, element s->node2b, element s->node2c;

    /* look for cGMP channel on postsynaptic side */

    tc = NULL;
    for (cpnt=elempnt; cpnt=foreach(cpnt,SYNAPSE,s->node2a,s->node2b,-1,-1,NULL,NULL,NULL,NULL);
		       cpnt=cpnt->next) {
       if ((apnt=get_chanattr(cpnt))!=NULL){
	 if (apnt->ctype==CGMP) {
            tc = epnt;
            break;
         }
       }
    }
    //printf("%d %d %d\n", s->node1a, s->node1b, s->node1c);

    if (!t) return;
    if (plname) sprintf(name,"rate %s",plname);
    else        sprintf(name,"rate %s_%d",cname[ct],cn);
    if (plname) sprintf(vesnam,"ves  %s",plname);
    else        sprintf(vesnam,"ves  %s_%d",cname[ct],cn);
    if (plname) sprintf(condnam,"cond %s",plname);
    else        sprintf(condnam,"cond %s_%d",cname[ct],cn);
    plot(FA9,s->elnum,rmax,rmin); plot_param(name,c1,plotnum,plsize);
    plot(FB1,s->elnum,fmax,fmin=0); plot_param(vesnam,green,plotnum,plsize);
    plot(FC1,s->elnum,fmaxx=2,fmin=0); plot_param(vesnam,cyan,plotnum,plsize);
    plot(FC9,s->elnum,fmaxx=2,fmin=0); plot_param("2ndmsg",red,plotnum,plsize);
//    plot(G,0,s->elnum,fmax=2e-11,fmin=0); plot_param("G(0)",red,plotnum);

    if (!tc) {
	plot(G,s->elnum,fmax=1e-10,fmin=0); plot_param(condnam,blue,plotnum);
    }
    else {
	plot(G,tc->elnum,fmax=1e-10,fmin=0); plot_param(condnam,blue,plotnum); 
    }

    //plot(FC0,s->elnum,fmax=2,rmin); plot_param(name,2,plotnum); 
    //plot(FC2,s->elnum,fmax=2,rmin); plot_param(name,4,plotnum);
    //plot(FC9,s->elnum,fmax=2,rmin); plot_param(name,4,plotnum);
    //plot(G,s->elnum,fmax=3e-10,rmin); plot_param(condnam,6,plotnum);
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, double rmin, double rmax, int pcolor, double plsize)
{
    const char *plname = NULL;
    double fmax=1e-4;

   plot_synrate_out(ct, cn, rmin, rmax, fmax, pcolor, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, double rmin, double rmax,int pcolor)

{
  plot_synrate_out(ct, cn, rmin, rmax, pcolor, 1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, double rmin, double rmax,int pcolor, const char *plname)

{
    double fmax=1e-4;

  plot_synrate_out(ct, cn, rmin, rmax, fmax, pcolor, plname, 1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, double rmin, double rmax, double fmax, 
						int pcolor, double plsize)
{
    const char *plname = NULL;

   plot_synrate_out(ct, cn, rmin, rmax, fmax, pcolor, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, double rmin, double rmax, double fmax, int pcolor)

{
  plot_synrate_out(ct, cn, rmin, rmax, fmax, pcolor, 1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, double rmin, double rmax, double fmax,
						int pcolor, const char *plname)

{
  plot_synrate_out(ct, cn, rmin, rmax, fmax, pcolor, plname, 1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int ct2, int cn2, double rmin, double rmax,
			double fmax, int pcolor, int prate, const char *plname, double plsize)

/* Display synaptic output to a specific cell */

{
    int c1,c2,t;
    int plotnum;
    double fmin;
    char name[CSIZ], vesnam[CSIZ], condnam[CSIZ];
    elem *cpnt, *epnt, *s, *tc;
    attrib *apnt;

  if (script==0) {
      plotnum = set_plotnum(ct) - 1;

    if (pcolor==-1)		/* if user doesn't specify a color */
      c2=plotnum;	  	/* use graph num to set color */
    else
      c2=pcolor;
    c1 = magenta;

    s = NULL;
    t = 0;
    for (epnt=elempnt; epnt=foreach(epnt,SYNAPSE,ct,cn,-1,-1,NULL,NULL,NULL,NULL);
		       epnt=epnt->next) {
       s = epnt;
       if ((s->node2a==ct2 || ct2<0) && (s->node2b==cn2 || cn2<0)) {
            t = s->elnum;
            break;
       }
       else s=NULL;
    }
    if (s==NULL) return;

    //fprintf(stderr,"node %d %d %d %d\n", s->node1a, s->node1b, s->node1c, 
    //					 s->node2a, s->node2b, s->node2c);

    /* look for cGMP channel on postsynaptic side */

    if (!t) return;
    if (plname) sprintf(name,"rate %s",plname);
    else        sprintf(name,"rate %s_%d_%d_%d",cname[ct],cn,ct2,cn2);
    if (plname) sprintf(vesnam,"ves  %s",plname);
    else        sprintf(vesnam,"ves  %s_%d",cname[ct],cn);
    if (plname) sprintf(condnam,"cond %s",plname);
    else        sprintf(condnam,"cond %s_%d",cname[ct],cn);
    if (prate) { plot(FA9,s->elnum,rmax,rmin); plot_param(name,c2,plotnum,plsize); }
    plot(FB4,s->elnum,fmax,fmin=0); plot_param(vesnam,c2+1,plotnum,plsize);
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int ct2, int cn2, double rmin, double rmax,
			int pcolor, int prate, double plsize)
{
    const char *plname = NULL;
    double fmax=2;

   plot_synrate_out(ct, cn, ct2, cn2, rmin, rmax, fmax, pcolor, prate, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int ct2, int cn2, 
		double rmin, double rmax, int pcolor, int prate)
{
   plot_synrate_out(ct, cn, ct2, cn2, rmin, rmax, pcolor,  prate,  1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int ct2, int cn2, 
		double rmin, double rmax, int pcolor, int prate, const char *plname)
{
     double fmax=2;

   plot_synrate_out(ct, cn, ct2, cn2, rmin, rmax, fmax, pcolor,  prate,  plname, 1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int ct2, int cn2, double rmin, double rmax, 
			double fmax, int pcolor, int prate, double plsize)
{
    const char *plname = NULL;

   plot_synrate_out(ct, cn, ct2, cn2, rmin, rmax, fmax, pcolor, prate, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int ct2, int cn2, 
		double rmin, double rmax, double fmax, int pcolor, int prate)
{

   plot_synrate_out(ct, cn, ct2, cn2, rmin, rmax, fmax, pcolor,  prate,  1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int ct2, int cn2, 
		double rmin, double rmax, double fmax, int pcolor, int prate, const char *plname)
{

   plot_synrate_out(ct, cn, ct2, cn2, rmin, rmax, fmax, pcolor,  prate,  plname, 1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, double rmin, double rmax,
			double fmax, int pcolor, int prate, int pves, const char *plname, double plsize)

/* Display synaptic output to a specific cell */

{
    int c1,c2,t;
    int plotnum;
    double fmin;
    char name[CSIZ], vesnam[CSIZ], condnam[CSIZ];
    elem *cpnt, *epnt, *s, *tc;
    attrib *apnt;

  if (script==0) {
      plotnum = set_plotnum(ct) - 1;

    if (pcolor==-1)		/* if user doesn't specify a color */
      c2=plotnum;	  	/* use graph num to set color */
    else
      c2=pcolor;
    c1 = magenta;

    s = NULL;
    t = 0;
    for (epnt=elempnt; epnt=foreach(epnt,SYNAPSE,ct,cn,nod,-1,NULL,NULL,NULL,NULL);
		       epnt=epnt->next) {
       s = epnt;
       if ((s->node2a==ct2 || ct2<0) && (s->node2b==cn2 || cn2<0)) {
            t = s->elnum;
            break;
       }
       else s=NULL;
    }
    if (s==NULL) return;

    //fprintf(stderr,"plot_synrate_out: node %d %d %d %d %d %d\n", s->node1a, s->node1b, s->node1c, 
    // 					 s->node2a, s->node2b, s->node2c);

    /* look for cGMP channel on postsynaptic side */

    if (!t) return;
    //sprintf(name,"rate %s_%d_%d_%d",cname[ct],cn,ct2,cn2);
    if (plname) sprintf(name,"R %s",plname);
    else        sprintf(name,"R%s%d",cname[ct],cn);
    if (plname) sprintf(vesnam,"v %s",plname);
    else        sprintf(vesnam,"v%s%d",cname[ct],cn);
    if (plname) sprintf(condnam,"cond %s",plname);
    else        sprintf(condnam,"cond %s_%d",cname[ct],cn);
    if (prate) { plot(FA9,s->elnum,rmax,rmin); plot_param(name,c2,plotnum,plsize); }
    if (pves) { plot(FB4,s->elnum,fmax,fmin=0); plot_param(vesnam,c2+1,plotnum,plsize); }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, double rmin, double rmax,
			int pcolor, int prate, int pves, double plsize)
{
  const char *plname = NULL;
  double fmax=2;

  plot_synrate_out(ct, cn, nod, ct2, cn2, rmin, rmax, 
		  fmax, pcolor, prate, pves, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, 
		double rmin, double rmax,int pcolor, int prate)
{
    int pves;

   plot_synrate_out(ct, cn, nod, ct2, cn2, rmin, rmax, pcolor,  prate, pves=0, 1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, 
		double rmin, double rmax,int pcolor, int prate, const char *plname)
{
    int pves;
    double fmax=2;

   plot_synrate_out(ct, cn, nod, ct2, cn2, rmin, rmax, 
		   fmax=2, pcolor,  prate, pves=0, plname, 1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, double rmin, double rmax,
			double fmax, int pcolor, int prate, int pves, double plsize)
{
  const char *plname = NULL;

  plot_synrate_out(ct, cn, nod, ct2, cn2, rmin, rmax, 
		  fmax, pcolor, prate, pves, plname, plsize);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, 
		double rmin, double rmax, double fmax, int pcolor, int prate)
{
    int pves;
   plot_synrate_out(ct, cn, nod, ct2, cn2, rmin, rmax, fmax, pcolor,  prate, pves=0, 1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_synrate_out(int ct, int cn, int nod, int ct2, int cn2, double rmin, 
		double rmax, double fmax, int pcolor, int prate, const char *plname)
{
    int pves;
   plot_synrate_out(ct, cn, nod, ct2, cn2, rmin, rmax, fmax, pcolor,  prate, pves=0, plname, 1.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_i_soma(int ct, int cn) 

{
    int pen;
    double fmax, fmin;
    double plg, offtr, offb;

  if (script==0) {
    plg   =  200e-12;		/* gain of trace */
    offb  = .0;			/* trace offset base */
    offtr = .8;			/* position of trace within graph */

    plot(I,nde(ct,cn,soma),fmax=(1-offtr)*plg+offb, fmin=(0-offtr)*plg+offb);
			plot_param("Igc", pen=2, 2.0);
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_i_nod(int ct,int cn,int n, double vmin, double vmax, int pcolor,
					const char *label,int plotnum,double psize)

{
    int c;
    char *name;
    char nbuf[CSIZ];
    node *npnt;

  if ((npnt=nde(ct,cn,n))==NULL) return;

  if (script==0) {
    if (plotnum<0) {
      plotnum = set_plotnum(ct) + 1;
    }
    
    if (pcolor==-1)	/* if user doesn't specify a color */
      c=plotnum;  	/* use graph num to set color */
    else
      c=pcolor;

    if (label==NULL || streq(label,"")) {
      name = nbuf;
      if (n==soma) sprintf(name,"I%s_%d_soma",cname[ct],cn);
      else         sprintf(name,"I%s_%d_%d",cname[ct],cn,n);
    }
    else name = (char *)label;

    if (psize < 0) {
      plot(I,nde(ct,cn,n), vmax, vmin); plot_param(name,c,plotnum);
    } else {
      plot(I,nde(ct,cn,n), vmax, vmin); plot_param(name,c,plotnum,psize);
    }
  }
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_currents(int ct, int plnum, double plgain) 

{
    int pen;
    double plg, offtr, offb, pmax, pmin;

  if (script==0) {
    plg = plgain;
    offtr = 0.5;
    offb  = 0;
    pmax = (1-offtr)*plg+offb;
    pmin = (0-offtr)*plg+offb;
    if (!notinit(chanelem[ct][_NA]) && celdens[ct][0][_NA][SOMA]>0)   
		plot(G,I,chanelem[ct][_NA], pmax, pmin); plot_param("Ina", pen=1, plnum=3);
    if (!notinit(chanelem[ct][_NA5]) && celdens[ct][0][_NA5][SOMA]>0) 
		plot(G,I,chanelem[ct][_NA5], pmax, pmin); plot_param("Ina5", pen=1, plnum=3);
    if (!notinit(chanelem[ct][_NA6]) && celdens[ct][0][_NA6][SOMA]>0)   
		plot(G,I,chanelem[ct][C_NA6], pmax, pmin); plot_param("Ina6", pen=1, plnum=3);
    if (!notinit(chanelem[ct][_NA8]) && celdens[ct][0][_NA8][SOMA]>0)   
		plot(G,I,chanelem[ct][_NA8], pmax, pmin); plot_param("Ina8", pen=1, plnum=3);
    //plot(G,I,chanelem[ct][_NA], pmax*10, pmin*100); plot_param("Ina2", pen=1, plnum=3);

    if (!notinit(chanelem[ct][_KDR]) && celdens[ct][0][_KDR][SOMA]>0)   
		plot(G,I,chanelem[ct][_KDR], pmax, pmin); plot_param("Ikdr", pen=2, plnum=3);
    // if (!notinit(chanelem[ct][_KDR]))  
    //		plot(G,I,chanelem[ct][_KDR], pmax, pmin); plot_param("Ikslo", pen=5, plnum=3);
    if (!notinit(chanelem[ct][_KA]) && celdens[ct][0][_KA][SOMA]>0)   
		plot(G,I,chanelem[ct][_KA], pmax, pmin); plot_param("Ika", pen=14, plnum=3);
    if (!notinit(chanelem[ct][_KH]) && celdens[ct][0][_KH][SOMA]>0)   
		plot(G,I,chanelem[ct][_KH], pmax, pmin); plot_param("Ih", pen=4, plnum=3);
    if (!notinit(chanelem[ct][_CA]) && celdens[ct][0][_CA][SOMA]>0)   
		plot(G,I,chanelem[ct][_CA], pmax, pmin); plot_param("Ica", pen=3, plnum=3);
    if (!notinit(chanelem[ct][_CA5]) && celdens[ct][0][_CA5][SOMA]>0)   
		plot(G,I,chanelem[ct][_CA5], pmax, pmin); plot_param("Ica5", pen=3, plnum=3);
    if (!notinit(chanelem[ct][_SKCA1]) && celdens[ct][0][_SKCA1][SOMA]>0)   
		plot(G,I,chanelem[ct][_SKCA1], pmax, pmin); plot_param("Isk1", pen=5, plnum=3);
    if (!notinit(chanelem[ct][_SKCA2]) && celdens[ct][0][_SKCA2][SOMA]>0)   
		plot(G,I,chanelem[ct][_SKCA2], pmax, pmin); plot_param("Isk2", pen=8, plnum=3);
    if (!notinit(chanelem[ct][_BKCA]) && celdens[ct][0][_BKCA][SOMA]>0)   
		plot(G,I,chanelem[ct][_BKCA], pmax, pmin); plot_param("IBk", pen=7, plnum=3);
    plot(I,nde(ct,1,soma), pmax, pmin); plot_param("Iinj", pen=12, plnum=3);
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double chan_plot(double ct, double cn, double n, double ctype, double stype, double parm, double xtime)

/* Function to be called by plot() at plot time, required to have 7 arguments. */

{
    node *npnt;
    conlst *clst;
    chan *chpnt;

  if ((npnt=nde(ct,cn,n))==NULL) return 0;
  clst = npnt->comptr->clst;
  if ((chpnt=findchan(clst,(int)ctype,(int)stype))==NULL) {
      if (ninfo >= 3)
       ncfprintf (stderr,"plot_chan: can't find %d %d chan at %d %d %d\n",(int)ctype,(int)stype,(int)ct,(int)cn,(int)n);
	   return 0;
  }
  return (record_chan(chpnt,(int)parm,0)); 

}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double chan_plot(double ct, double cn, double n, double ctype, double stype, double parm, double pval, double xtime)

/* Function to be called by plot() at plot time, required to have 7 arguments. */
/*  Plot one state (pval) of a channel */

{
    node *npnt;
    conlst *clst;
    chan *chpnt;

  if ((npnt=nde(ct,cn,n))==NULL) return 0;
  clst = npnt->comptr->clst;
  if ((chpnt=findchan(clst,(int)ctype,(int)stype))==NULL) {
      if (ninfo >= 3)
       ncfprintf (stderr,"plot_chan: can't find %d %d chan at %d %d %d\n",(int)ctype,(int)stype,(int)ct,(int)cn,(int)n);
	   return 0;
  }
  return (record_chan(chpnt,(int)parm,(int)pval)); 

}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

// void plot_chan(int ct, int cn, int n, int ctype, int stype, int param, double mult, double pmax, double pmin) 

/* Plot channel current from arbitrary node.  Problem is that channel may be generated
   by an "attrib" so it has no element number.  So we make a dummy channel, and set its
   "chanelem" parameter. This causes the runtime function "record_chan()" to look for the
   compartment pointed to by the node, and find the channel with the same type and subtype 
   as the dummy channel. Save "mult" as a multiplier for the channel current at runtime.
*/

// {
//    int chan_elem;
//    node *npnt;
//    elem *epnt;
//    conlst *lpnt;
//    chan *cpnt;
// 
//    npnt = nde (ct,cn,n);
//    if (!npnt) 
//       {ncfprintf (stderr,"plot_chan_current: can't find node %d\n",n);
//    }
// 
//    epnt = at (ndn(ct,cn,n), CHAN);	/* make dummy channel */
//    make_chan(epnt,ctype,stype);		/* make dummy chan for recording */
//    chset(epnt);				/* set "chanelem" param of chan element */
//    chan_elem = epnt->elnum; 
//    plot_func(chan_plot,chan_elem,(double)param,mult,pmax,pmin);
// }

const char *chan_name (int ctype, int stype) 

{
     int n;

   switch (ctype) {
	   case AMPA:
		      switch (stype) {
			default:
			case 1: n = C_AMPA1; break;
			case 2: n = C_AMPA2; break;
			case 3: n = C_AMPA3; break;
			case 4: n = C_AMPA4; break;
			case 5: n = C_AMPA5; break;
		     } break;
	   case NMDA:  
		      switch (stype) {
			default:
			case 1: n = C_NMDA1; break;
			case 2: n = C_NMDA2; break;
		     } break;
	   case GABA: 
		      switch (stype) {
			default:
			case 1: n = C_GABA1; break;
			case 2: n = C_GABA2; break;
			case 3: n = C_GABA3; break;
			case 4: n = C_GABA4; break;
		     } break;
	   case GLY:   
		      switch (stype) {
			default:
			case 1: n = C_GLY; break;
		     } break;
	   case CGMP:
		      switch (stype) {
			default:
			case 1: n = C_CGMP1; break;
			case 2: n = C_CGMP2; break;
			case 3: n = C_CGMP3; break;
			case 4: n = C_CGMP4; break;
			case 5: n = C_CGMP5; break;
			case 6: n = C_CGMP6; break;
			case 7: n = C_CGMP7; break;
			case 8: n = C_CGMP8; break;
			case 9: n = C_CGMP9; break;
			case 10: n = C_CGMP10; break;
			case 11: n = C_CGMP11; break;
		     } break;
	    default: n = C_AMPA1; break;
      }
   if (n > C_CGMP11) n = C_CGMP11;
   return chname[n];
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_chan_states(int ct, int cn, int n, int ctype, int stype, int param, 
				double pmax, double pmin, const char *label, int plotnum, double psize) 

/* Plot channel current from arbitrary node.  */

{
      int i, numstate, pcolor;
      chantype *chpnt;
      char nbuf[CSIZ];
      char *name;

  if (n == -1) return;
  if (nde(ct,cn,n)==NULL) return;

  if (script==0) {
    if (plotnum<0) {
      plotnum = set_plotnum(ct) + 1;
    }
    if (chpnt=getchantype(ctype,stype)) numstate = chpnt->numstate;
    else numstate = 0;

    for (i=1; i<=numstate; i++) {
   	plot_func(chan_plot,(double)ct,(double)cn,(double)n,
			    (double)ctype,(double)stype,(double)param,(double)i,pmax,pmin);
        if (label==NULL || streq(label,"")) {
          name = nbuf;
          sprintf(name,"%s_%d_%s_%d_%d",chan_name(ctype,stype),i,cname[ct],cn,n);
        }
        else name = (char *)label;
        if (psize < 0) {
            plot_param(name, pcolor=i, plotnum);
         } else {
            plot_param(name, pcolor=i, plotnum, psize);
         }
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_chan(int ct, int cn, int n, int ctype, int stype, int param, int pval, double pmax, double pmin) 

/* Plot channel param from arbitrary node.  */

{
   plot_func(chan_plot,(double)ct,(double)cn,(double)n,(double)ctype,(double)stype,(double)param,
		   		  (double)pval,pmax,pmin);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_chan(int ct, int cn, int n, int ctype, int stype, int param, double pmax, double pmin) 

/* Plot channel param from arbitrary node.  */

{
  plot_chan(ct, cn, n, ctype, stype, param, 0, pmax, pmin);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_chan_current(int ct, int cn, int n, int ctype, int stype, double mult, double pmax, double pmin)
{
    plot_chan(ct, cn, n, ctype, stype, I, pmax, pmin);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_chan_current(int ct, int cn, int n, int ctype, int stype, double pmax, double pmin)
{
    double mult;

    plot_chan(ct, cn, n, ctype, stype, I, pmax, pmin);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_chan_cond(int ct, int cn, int n, int ctype, int stype, double mult, double pmax, double pmin)
{
    plot_chan(ct, cn, n, ctype, stype, G, pmax, pmin);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void plot_chan_cond(int ct, int cn, int n, int ctype, int stype, double pmax, double pmin)
{
    plot_chan(ct, cn, n, ctype, stype, G, pmax, pmin);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

struct synapse_list {
	int index;
	synapse **list;
};

#define SLISTLEN 50
#define SLISTSIZ 2000
synapse_list synap_rec_list[SLISTLEN] = {0};

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  */

synapse_list *make_synapse_list(int synlist)
{
   if (synlist>=SLISTLEN) {
	   fprintf (stderr,"# make_synapse_list: synapse list number too big %d\n",synlist);
	   synlist = SLISTLEN - 1;
   }
   if (synap_rec_list[synlist].list==NULL) {
       synap_rec_list[synlist].list = (synapse **)emalloc(SLISTSIZ*sizeof(int));
       synap_rec_list[synlist].index = 0;
   }
   if (synap_rec_list[synlist].list==NULL) {
	   ncfprintf (stderr," make_synapse_list: not enough memory for list %d\n",synlist);
	   return NULL;
   }
   return (&synap_rec_list[synlist]);	/* return list of synapses */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int synapse_add (int synlist, int ct, int cn, int nod, int ct2, int cn2)
{

/* Make array of synapse elements for recording at runtime. */
/* Synlist defines which list */

/*  Originally from expt_dsgc_sbac_bar.cc
 
   usage:  synapse_add (1, sbac, 1, 1215, sbac, 2);          // add synapse for rate plot
*/ 

   int nsyn;
   elem *epnt;
   synapse_list *s;

   if ((s=make_synapse_list(synlist))==NULL) return 0;

   for (nsyn=0,epnt=elempnt; epnt=foreach(epnt,SYNAPSE,ct,cn,nod,ct2,cn2,-1,0); epnt=epnt->next) {

     // fprintf (stderr,"synapse_add 1: list %d, %d %d %d, %d %d %d, %d %d %d\n",synlist,ct,cn,nod,
     //             epnt->node1a,epnt->node1b,epnt->node1c,
     //             epnt->node2a,epnt->node2b,epnt->node2c); /* */
            // fprintf (stderr,"synapse_add 2: %d %d %d\n",ct,cn,nod);

      if (cn2 < 0 || epnt->node2b==cn2) {
         s->list[s->index] = (synapse *)epnt;		// save pointer to synapse element
         if (++s->index>=SLISTSIZ) s->index--;
         nsyn++;
      } else continue;
   }
   return nsyn;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int synapse_add (int synlist, int ct, int cn, int nod, int ct2, int cn2, double disti, double disto)
{

/* Make array of synapse elements for recording at runtime, 
 * between disti and disto from soma of ct2,cn2. */

/* Synlist defines which list */

/*  Originally from expt_dsgc_sbac_bar.cc
 
   usage:  synapse_add (1, sbac, 1, 1215, sbac, 2);          // add synapse for rate plot
*/ 

   int nsyn;
   elem *epnt;
   synapse_list *s;
   double dist;

   if ((s=make_synapse_list(synlist))==NULL) return 0;

   for (nsyn=0,epnt=elempnt; epnt=foreach(epnt,SYNAPSE,ct,cn,nod,ct2,cn2,-1,0); epnt=epnt->next) {

     //fprintf (stderr,"synapse_add 1: list %d, %d %d %d, %d %d %d, %d %d %d\n",synlist,ct,cn,nod,
     //             epnt->node1a,epnt->node1b,epnt->node1c,
     //             epnt->node2a,epnt->node2b,epnt->node2c); /* */
            // fprintf (stderr,"synapse_add 2: %d %d %d\n",ct,cn,nod);

      if (cn2 < 0 || epnt->node2b==cn2) {

	 dist = dist2d(ndn(epnt->node2a,epnt->node2b,soma),epnt->nodp2); // find distance postsyn to soma
	 if (disti <= dist && dist <= disto) {		// if within distance limits 
             s->list[s->index] = (synapse *)epnt;	// save pointer to synapse element
             if (++s->index>=SLISTSIZ) s->index--;
	     nsyn++;
         } 
      } else continue;
   }
   return nsyn; 
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int synapse_add (int synlist, int ct, int cn, int nod, int ct2, int cn2, 
			double disti, double disto, double xloc, double yloc)
{

/* Make array of synapse elements for recording at runtime, 
 * between disti and disto from xloc, yloc. */

/* Synlist defines which list */

/*  Originally from expt_dsgc_sbac_bar.cc
 
   usage:  synapse_add (1, sbac, 1, 1215, sbac, 2);          // add synapse for rate plot
*/ 

   int nsyn;
   double dist,xr,yr;
   elem *epnt;
   synapse_list *s;

   if ((s=make_synapse_list(synlist))==NULL) return 0;

   for (nsyn=0,epnt=elempnt; epnt=foreach(epnt,SYNAPSE,ct,cn,nod,ct2,cn2,-1,0); epnt=epnt->next) {

     //fprintf (stderr,"synapse_add 1: list %d, %d %d %d, %d %d %d, %d %d %d\n",synlist,ct,cn,nod,
     //             epnt->node1a,epnt->node1b,epnt->node1c,
     //             epnt->node2a,epnt->node2b,epnt->node2c); /* */
            // fprintf (stderr,"synapse_add 2: %d %d %d\n",ct,cn,nod);

      if (cn2 < 0 || epnt->node2b==cn2) {
	 xr = xloc - epnt->nodp2->xloc;
	 yr = yloc - epnt->nodp2->yloc;
	 dist = sqrt(xr*xr + yr*yr); 		// find 2D distance from postsyn to (xloc,yloc)
	 if (disti <= dist && dist <= disto) {	// if within distance limits 
             s->list[s->index] = (synapse *)epnt;	// save pointer to synapse element
             if (++s->index>=SLISTSIZ) s->index--;
	     nsyn++;
         }
      } else continue;
   }
   return nsyn; 
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int synapse_add (int synlist, int ct, int cn, int nod, int ct2, int cn2, double vrev)
{

/* Make array of synapse elements for recording at runtime. */
/* Synlist defines which list */

/*  Originally from expt_dsgc_sbac_bar.cc
 
   usage:  synapse_add (1, sbac, 1, 1215, sbac, 2);          // add synapse for rate plot
*/ 

   int nsyn;
   elem *epnt;
   synapse_list *s;

   if ((s=make_synapse_list(synlist))==NULL) return 0;

   for (nsyn=0,epnt=elempnt; epnt=foreach(epnt,SYNAPSE,ct,cn,nod,ct2,cn2,-1,0); epnt=epnt->next) {

     //fprintf (stderr,"synapse_add 1: list %d, %d %d %d, %d %d %d, %d %d %d\n",synlist,ct,cn,nod,
     //             epnt->node1a,epnt->node1b,epnt->node1c,
     //             epnt->node2a,epnt->node2b,epnt->node2c); /* */
            // fprintf (stderr,"synapse_add 2: %d %d %d\n",ct,cn,nod);

      if (((synapse *)epnt)->vrev==vrev && ((cn2 < 0) || (epnt->node2b==cn2))) {

          s->list[s->index] = (synapse *)epnt;	// save pointer to synapse element
          if (++s->index>=SLISTSIZ) s->index--;
	  nsyn++;
      }
   }
   return nsyn; 
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int synapse_add (int synlist, int ct, int cn, int nod, int ct2, int cn2, int connum)
{

/* Make array of synapse elements for recording at runtime. */
/* Synlist defines which list */

/*  Originally from expt_dsgc_sbac_bar.cc
 
   usage:  synapse_add (1, sbac, 1, 1215, sbac, 2);          // add synapse for rate plot
*/ 

   int nsyn;
   elem *epnt;
   synapse_list *s;

   if ((s=make_synapse_list(synlist))==NULL) return 0;

   for (nsyn=0,epnt=elempnt; epnt=foreach(epnt,SYNAPSE,ct,cn,nod,ct2,cn2,-1,0); epnt=epnt->next) {

     //fprintf (stderr,"synapse_add 1: list %d, %d %d %d, %d %d %d, %d %d %d\n",synlist,ct,cn,nod,
     //             epnt->node1a,epnt->node1b,epnt->node1c,
     //             epnt->node2a,epnt->node2b,epnt->node2c); /* */
            // fprintf (stderr,"synapse_add 2: %d %d %d\n",ct,cn,nod);

      if (((synapse *)epnt)->connum==connum && ((cn2 < 0) || (epnt->node2b==cn2))) { // if correct postsyn conn number

          s->list[s->index] = (synapse *)epnt;	// save pointer to synapse element
          if (++s->index>=SLISTSIZ) s->index--;
	  nsyn++;
      }
   }
   return nsyn; 
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int synapse_add (int synlist, int ct, int cn, int nod, int ct2, int cn2, int connum, double disti, double disto)
{

/* Make array of synapse elements for recording at runtime, 
 * between disti and disto from soma of ct2,cn2. */

/* Synlist defines which list */

/*  Originally from expt_dsgc_sbac_bar.cc
 
   usage:  synapse_add (1, sbac, 1, 1215, sbac, 2);          // add synapse for rate plot
*/ 

   int nsyn;
   elem *epnt;
   double dist;
   synapse_list *s;

   if ((s=make_synapse_list(synlist))==NULL) return 0;

   for (nsyn=0,epnt=elempnt; epnt=foreach(epnt,SYNAPSE,ct,cn,nod,ct2,cn2,-1,0); epnt=epnt->next) {

     //fprintf (stderr,"synapse_add 1: list %d, %d %d %d, %d %d %d, %d %d %d\n",synlist,ct,cn,nod,
     //             epnt->node1a,epnt->node1b,epnt->node1c,
     //             epnt->node2a,epnt->node2b,epnt->node2c); /* */
            // fprintf (stderr,"synapse_add 2: %d %d %d\n",ct,cn,nod);

      if (((synapse *)epnt)->connum==connum && ((cn2 < 0) || (epnt->node2b==cn2))) {  // if correct postsyn conn number

	 dist = dist2d(ndn(epnt->node2a,epnt->node2b,soma),epnt->nodp2); // find distance postsyn to soma
	 if (disti <= dist && dist <= disto) {		// if within distance limits 

             s->list[s->index] = (synapse *)epnt;	// save pointer to synapse element
             if (++s->index>=SLISTSIZ) s->index--;
	     nsyn++;
         }
      }
   }
   return nsyn; 
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double rsyn_avg (double sl, double time)
{
   int i,n;
   int synlist;
   double avgrate=0.0;
   synapse_list *s;

   synlist = int(sl);
   s = &synap_rec_list[synlist];
   for (i=n=0; i<s->index; i++) {
         if (s->list[i]!=NULL) {
               avgrate += record_synapse(s->list[i],FA9);
	       n++;
	 }
   }
   if (n==0) return 0;
   return (avgrate/n);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double isyn_avg (double sl, double time)
{
   int i,n;
   int synlist;
   double totcur=0.0;
   synapse_list *s;

   synlist = int(sl);
   s = &synap_rec_list[synlist];
   for (i=n=0; i<s->index; i++) {
         if (s->list[i]!=NULL) {
              totcur += record_chan(s->list[i],I,0);
	      n++;
	 }
   }
   if (n==0) return 0;
   else      return (totcur/n);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double isyn_tot (double sl, double time)
{
   int i;
   int synlist;
   double totcur=0.0;
   synapse_list *s;

   synlist = int(sl);
   s = &synap_rec_list[synlist];
   for (i=0; i<s->index; i++) {
        if (s->list[i]!=NULL)
             totcur += record_chan(s->list[i],I,0);
   }
   return (totcur);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double gsyn_avg (double sl, double time)
{
   int i,n;
   int synlist;
   double totcond=0.0;
   synapse_list *s;

   synlist = int(sl);
   s = &synap_rec_list[synlist];
   for (i=n=0; i<s->index; i++) {
         if (s->list[i]!=NULL) {
                 totcond += record_chan(s->list[i],G,0);
                 n++;
   // fprintf (stderr,"gsyn_avg: %d %g cbp %d %d\n",
   //		    synlist,record_chan(s->list[i],G,0),s->list[i]->node1b,s->list[i]->node1c);
         }
  }
  if (n==0) return 0;
  else      return (totcond/n);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double gsyn_tot (double sl, double time)
{
   int i;
   int synlist;
   double totcond=0.0;
   int *synap_list;
   elem *elpnt;
   synapse *sepnt;
   synap *spnt;
   cacomp *capnt;
   double cais, vmax;
   int nshells;
   synapse_list *s;

   synlist = int(sl);
   s = &synap_rec_list[synlist];
   for (i=0; i<s->index; i++) {
         if (s->list[i]!=NULL)
               totcond += record_chan(s->list[i],G,0);
/* 
// print out list for debugging
//
	         //if (simtime > 0.05 && simtime < 0.051) {
	         if (simtime > 0.01 && simtime < 0.012) {
                  if (elpnt=findelem(synap_list[i])) {
                    sepnt = ((synapse *)elpnt);
                    spnt = (synap*)sepnt->lptr;
		    if ((capnt=spnt->comp1->capnt)!=NULL) {
                         cais = capnt->cais[0];
			 vmax = capnt->vmax;
			 nshells = capnt->cashell;
		    } else { cais = 0; vmax = -1; nshells = -1;}

		    fprintf (stderr,"sl %d prenod %d %d %4d postnod %d %d %4d cond %-9.4g v %-9.4g ca %-9.4g vmax %-8.4g cap %-9.4g %2d \n",
				synlist, sepnt->node1a, sepnt->node1b, sepnt->node1c,
					     sepnt->node2a, sepnt->node2b, sepnt->node2c,
                                             spnt->resp1->conduct, spnt->comp1->v, cais, vmax, spnt->comp1->cap,nshells);
                  }
                 }
//
//
/* */

   }
   return (totcond);
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double gnmda_syn_tot (double sl, double time)
{
   int i;
   int synlist;
   double totcond=0.0;
   synapse_list *s;

   synlist = int(sl);
   s = &synap_rec_list[synlist];
   for (i=0; i<s->index; i++) {
         if (s->list[i]!=NULL)
                totcond += nmda_act_plot(s->list[i],time);
    }
    return (totcond);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double casyn_avg (double sl, double time)
{
   int i;
   int synlist;
   double totca=0.0;
   synapse_list *s;

   synlist = int(sl);
   s = &synap_rec_list[synlist];
   for (i=0; i<s->index; i++) {
         if (s->list[i]!=NULL)
                 totca += record_synapse(s->list[i],CA);
//    fprintf (stderr,"casyn_avg: %g cbp %d %d\n",record_synapse(s->list[i],CA),s->list[i]->node1b,s->list[i]->node1c);
  }
  return (totca/s->index);
}

