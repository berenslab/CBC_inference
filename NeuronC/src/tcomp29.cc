/* tcomp29: plot of voltage vs distance along cable */


#include <stdio.h>
#include "ncfuncs.h"
#include "ncinit.h"

main(int argc, char **argv)
{
   double totlen,rnodes,nsegs,seglen,segdia;
   double stimdur,stime,start,dur;
   double max,min;
   cable *c;
   double t;
   int i;

 setptr("rnodes",      &rnodes); 
 ncinit(argc,argv);

 rnodes = 40;

 timinc = 5e-5;
 ploti = timinc;
 endexp = .01;
 implicit = 1;

 drm = 5000;
 totlen = 800;

 setvar();

 nsegs = rnodes - 1;
 seglen = totlen / nsegs;
 segdia = 2;

 for (i=1; i<=nsegs; i++) {
   c = make_cable(nd(i),nd(i+1)); c->dia=segdia; c->length=seglen;
 };

 stimdur = .0005;
 vclamp (ndn(1),-0.02, start=0, dur=stimdur);

 /* plot 7,6,5,4,3,2,1; */

 graph_x (max=nsegs*seglen,min=0);		/* volts */
 graph_y (max=-.02,min=-.07);			/* psp (volts) (expon)  */
 graph_pen (4);
 graph_init();					/* draw axes */
 
 stime = 5e-5;

 for (t=0; t<endexp; t += stime) {
    graph_restart();
    step (stime);
    if (t+1e-6 >= stimdur) graph_pen (1);
    for (i=1; i<=nsegs; i++) {
      graph (i*seglen, v(ndn(i)));
    };
 
 };  /* for (t;;) */

 ncexit();
}

