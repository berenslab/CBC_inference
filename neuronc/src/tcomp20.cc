/* tcomp20:  test of voltage clamp current recording
   from a cone with simultaneous voltage recording */

/* calibrate at inten 100e3 for 1e-3 dur to give
   8e-12 A peak response */


#include "ncfuncs.h"
#include <stdio.h>

double conerest, conerm;

#include "con2.cc"

void setptrs(void) {}

main(int argc, char **argv)
{
   double flash,inten,v;
   double start, dur;
   double stimstart;
   double max,min;
   int n;

 ncinit(argc,argv);

 plmax = -.01;
 plmin = -.06;
 endexp = .4;
 ploti = endexp / 200;
 ploti = .001;

 drm = 5000;
 conerest = -.07;
 conerm   = 5000;
 
 n=mcone (0,0,1);
 n=mcone (0,0,2);

 flash= 100e3;

 stimstart = 0.1;

 stim_cone (ndn(1),inten=flash,start=stimstart,dur=.001);
 vclamp    (ndn(1),v=-.07,start=.008,dur=1.5);
 stim_cone (ndn(2),inten=flash,start=stimstart,dur=.001);

 plot (I, ndn(1),   max=0,    min=-2e-10);
 plot (V, ndn(2),   max=.040, min=-.080); 
 //plot (V, ndn(2), max=.04,  min=-.08);
 plot (L, ndn(2),   max=1e7,  min=0);
 
 run();
 ncexit();
}


