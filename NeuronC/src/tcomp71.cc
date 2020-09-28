/* tcomp66: test of AMPA synapse */


#include <stdio.h>
#include "ncfuncs.h"
#include "ncinit.h"

main(int argc, char **argv)
{
    int xxx;
    double dia, max, min;
    double start,dur,stimdur,stimtime1,stimtime2;
    sphere *d;
    synapse *s;
    elem *e;
    chattrib *a;
    nattrib *n;

 ncinit(argc,argv);

 timinc = 1e-5;
 ploti  = 1e-4;
 endexp = 0.2;
 rseed = 104;
 tempcel = 37;
 plsep = 1;

 d=make_sphere(nd(1), dia=1); d->vrest=-.07; d->Rm = 1e10;

 e = make_chan (ndn(1),CA,0);
 a = (chattrib *)e->attpnt;
 a->maxcond = 1e-10;

 stimtime1 = 0.01;
 stimdur  = 0.02;

 vclamp (ndn(1), -0.0, start=stimtime1, dur=stimdur);
 vclamp (ndn(1), -0.07, start=stimtime1+stimdur, dur=1);

 stimtime2 = 0.15;

 //plot(V,ndn(1), max=-0.0, min=-0.07);

 plot(I, e->elnum,     max=1e-12, min=0);
 plot(I,ndn(1),   max=1e-12, min=0);


 run();
 ncexit();
}

