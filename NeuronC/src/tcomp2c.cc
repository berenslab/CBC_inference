/* tcomp2:  test of cable with segments and external compartment */
/*    Use to test effect of external stimulation */

/* Can set params before setvar() from command line */

#include <stdio.h>
#include "ncfuncs.h"
#include "ncinit.h"
#include "colors.h"

extern comp *compnt;

main(int argc, char **argv)
{
   double seg, start, dur, max, min;
   double predur,startime;
   comp *ext_comp,*pnt;
   cable *c;
   sphere *s;

 ncinit(argc,argv);

 timinc = 1e-6; 
 endexp = 0.05;
 complam=.01;
 relax = .15;
 crit=1e-6;
 ploti = 1e-4;
 //implicit=1;
 //euler=1;
 startime=0.01;
 // plsep = 1;

 plmax = -0.02;
 plmin = -0.07;

 // seg = 353;
 seg = 20;			// Change this to 5 to see effect of small stim region
 dri = 200;
 drm = 50000;
 dcm = 1e-6;

 setvar();

 c = make_cable (nd(1), nd(2)); c->length=seg; c->dia=0.5;
 // s = make_sphere (nd(1), 1);
 s = make_sphere (nd(3), 1);

 setxmin = 0;                    // set plot to start at 0
 predur = 0.01;
 simtime = 0 - predur;

 step (predur);

 ext_comp = ndn(3)->comptr;
 addcomp_extern(ext_comp,ndn(1)->comptr);
 
// for (pnt=compnt; pnt; pnt=pnt->next) {	/* use this loop to set many vext sites */
//    if (pnt != ext_comp)
//       addcomp_extern(ext_comp,pnt);
// }

 vclamp(nd(3), 0,     start=simtime, dur=startime);
 vclamp(nd(3), 0.005, start=startime, dur=1);
 //vclamp(nd(3), 0.005, start=startime, dur=startime);
 // vclamp(nd(3), 0,     start=startime*2, dur=1);

 plot(V,nd(1),plmax,plmin);
 plot(V,nd(2),plmax,plmin);
 plot(V,nd(3),0.025,-0.025);
 graph_pen (RED,GREEN,BLUE);
 // plot(I,nd(3),max=5e-10,min=0);
 run();
 ncexit();
}

