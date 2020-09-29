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
    chattrib *a;
    nattrib *n;

 ncinit(argc,argv);

 timinc = 1e-5;
 ploti  = 1e-4;
 endexp = 0.2;
 rseed = 104;
 tempcel = 37;

 d=make_sphere(nd(1), dia=10); d->vrest=-.07;

 d=make_sphere(nd(2), dia=10); d->Rm=2000;

 s=make_synapse(ndn(1),ndn(2)); s->ntact=OPEN; s->maxcond=1e-9; 
	s->trconc=1e-3; s->curve=EXPON; s->ngain=3;
	s->nfilt2=1; s->timec2=make_filt(2); s->tfall2 = 20;
	xxx = s->elnum;

    a = make_chan ((elem *)s,AMPA,2); a->taua=1; a->taub=1;
    n = make_chnoise((elem *)s); n->n=1000; 
/*
 conn 1 to 2 synapse open maxcond 1e-9 
		resp ampa type 2 taum=1 trconc=1e-3
		chnoise=1 N=1000
		timec2 2 tfall2 = 20
		expon 3
                ename xxx;
*/


 stimtime1 = 0.01;
 stimdur  = 0.002;

 vclamp (ndn(1), -0.03, start=stimtime1, dur=stimdur);
 vclamp (ndn(1), -0.07, start=stimtime1+stimdur, dur=1);

 stimtime2 = 0.15;

 vclamp (ndn(1), -0.03, start=stimtime2, dur=stimdur);
 vclamp (ndn(1), -0.07, start=stimtime2+stimdur, dur=1);

 plot(V,ndn(1), max=-0.02, min=-0.07);
 plot(V,ndn(2), max=-0.02, min=-0.07);

 plot(FB1,   xxx, max=10, min=0);
 plot(G, 0,  xxx, max=1e-9, min=0);
 plot(G, 1,  xxx, max=2, min=0);
 plot(G, 2,  xxx, max=2, min=0);
 plot(G, 3,  xxx, max=2, min=0);
 plot(G, 4,  xxx, max=2, min=0);
 plot(G, 5,  xxx, max=2, min=0);
 plot(G, 6,  xxx, max=2, min=0);
 plot(G, 7,  xxx, max=2, min=0);
 
 run();
 ncexit();
}

