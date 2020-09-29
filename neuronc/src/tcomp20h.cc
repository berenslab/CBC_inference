/* tcomp20:  test of voltage clamp current recording from cones */

/* calibrate at inten 100e3 for 1e-3 dur to give
   8e-12 A peak response */


#include "ncfuncs.h"
#include "ndef.h"
#include <stdio.h>

double conerest, conerm;

#include "con2.cc"

void setptrs(void) {}

main(int argc, char **argv)
{
   double bkgr,flash,inten,v;
   double start, dur;
   double stimstart;
   double lscale,max,min;
   int n;

 ncinit(argc,argv);

 setptr("bkgr", &bkgr);
 setptr("flash", &flash);

 setvar();
 endexp = .5;
 ploti = endexp / 200;
 ploti = .001;

 drm = 5000;
 conerest = -.07;
 conerm   = 5000;

 if (notinit(bkgr))   bkgr = 10000;

 n=mconep (0,0,1,1,bkgr);
 n=mconep (0,0,2,19,bkgr);
 stim_cone (ndn(1),inten=bkgr,start=0,dur=1);
 stim_cone (ndn(2),inten=bkgr,start=0,dur=1);

 if (bkgr > 1) if (notinit(flash)) flash= bkgr * 1.0;
 lscale = max(bkgr,flash);

 stimstart = 0.1;

 stim_cone (ndn(1),inten=flash,start=stimstart,dur=.001);
 vclamp    (ndn(1),v=-.07,start=.008,dur=1.5);
 stim_cone (ndn(2),inten=flash,start=stimstart,dur=.001);
 vclamp    (ndn(2),v=-.07,start=.008,dur=1.5);

 plot (I, ndn(2),   max=0,    min=-1e-10);
 plot (I, ndn(1),   max=0,    min=-1e-10);
 //plot (V, ndn(2),   max=.040, min=-.080); 
 //plot (V, ndn(2), max=.04,  min=-.08);
 plot (L, ndn(2),   max=lscale*20,  min=0);
 
 run();
 ncexit();
}


