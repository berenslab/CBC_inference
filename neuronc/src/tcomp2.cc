/* tcomp2:  test of cable with segments */

#include <stdio.h>
#include "ncfuncs.h"
#include "ncinit.h"

main(int argc, char **argv)
{
   double seg, start, dur, max, min;
   cable *c;

 ncinit(argc,argv);

 timinc = 1e-4; 
 endexp = 0.04;
 complam=.1;
 relax = .15;
 crit=1e-6;
 ploti = 1e-4;
 implicit=1;

 plmax = 0;
 plmin = -.07;

 seg = 353;
 drm = 5000;

 setvar();

 c = make_cable (nd(1), nd(2)); c->length=seg; c->dia=1;
 c = make_cable (nd(2), nd(3)); c->length=seg; c->dia=1;
 c = make_cable (nd(3), nd(4)); c->length=seg; c->dia=1;
 c = make_cable (nd(4), nd(5)); c->length=seg; c->dia=1;
 c = make_cable (nd(5), nd(6)); c->length=seg; c->dia=1;
 c = make_cable (nd(6), nd(7)); c->length=seg; c->dia=1;
 c = make_cable (nd(7), nd(8)); c->length=seg; c->dia=1;
 c = make_cable (nd(8), nd(9)); c->length=seg; c->dia=1;
 c = make_cable (nd(9), nd(10)); c->length=seg; c->dia=1;

 vclamp(nd(1), -0.01, start=0.01, dur=1);

 plot(V,nd(1));
 plot(I,nd(1),max=5e-10,min=0);
 plot(V,nd(2));
 plot(V,nd(3));
 plot(V,nd(4));
 plot(V,nd(5));
 plot(V,nd(6));
 plot(V,nd(7));
 run();
 ncexit();
}

