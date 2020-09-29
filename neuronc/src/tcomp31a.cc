/* tcomp31a: Test of KB (slowly inactivating K current) */


#include <stdio.h>
#include "ncfuncs.h"
#include "ncinit.h"

main(int argc, char **argv)

{
   int i,n,nchan;
   double start, dur, max, min;
   double prevolt, vpulse, prestimdur, stimdur, poststimdur;
   double dia;
   sphere *s;
   elem *e;
   chattrib *a;

 ncinit(argc,argv);

 timinc = 1e-4;
 crit   = 1e-8;
 endexp = .02;
 ploti  = 1e-4;
 
 tempcel= 10;

 drm=1e6;

 nchan = 5;

 s = make_sphere(nd(1), dia=1); s->vrev= -0.07;
 e = at(ndn(1), CHAN); a=make_chan (e,NA,0); a->vrev=0.04; a->maxcond=1e-9;

 s = make_sphere(nd(2), dia=1); s->vrev= -0.07;
 e = at(ndn(2), CHAN); a=make_chan (e,NA,0); a->vrev=0.04; a->maxcond=1e-9;

 s = make_sphere(nd(3), dia=1); s->vrev= -0.07;
 e = at(ndn(3), CHAN); a=make_chan (e,K,0); a->vrev=-0.08; a->maxcond=1e-9;

 s = make_sphere(nd(4), dia=1); s->vrev= -0.07;
 e = at(ndn(4), CHAN); a=make_chan (e,K,1); a->vrev=-0.08; a->maxcond=1e-9;

 s = make_sphere(nd(5), dia=1); s->vrev= -0.07;
 //c = at(ndn(3), CHAN); a=make_chan (e,K,3); a->vrev=-0.08; a->maxcond=1e-9;
 
 s = make_sphere(nd(5), dia=1); s->vrev= -0.07;
 //c = at(ndn(6), CHAN); a=make_chan (e,K,3); a->vrev=-0.08; a->maxcond=1e-9;

 plot (I, ndn(1), max=1e-11, min=-1e-11);
 //plot (I, ndn(2), max=1e-11, min=-1e-11);
 plot (I, ndn(3), max=1e-10, min=-1e-10);
 //plot (I, ndn(4), max=1e-10, min=-1e-10);
 //plot (I, ndn(5), max=1e-10, min=-1e-10);
 //plot (I, ndn(6), max=1e-10, min=-1e-10);

 prevolt = -.08;
 vpulse = -.07;
 prestimdur = .005;
 stimdur = .02;
 poststimdur = .05;
 
 for (i=0; i<10; i++) {
   simtime = 0;
   for (n=1; n<=nchan; n++) 
      vclamp (ndn(n),prevolt,start=simtime,dur=prestimdur);
   step (prestimdur); 
   for (n=1; n<=nchan; n++)
      vclamp (ndn(n),vpulse,start=simtime,dur=stimdur);
   step (stimdur); 
   for (n=1; n<=nchan; n++)
      vclamp (ndn(n),prevolt,start=simtime,dur=poststimdur);
   step(poststimdur); 
   vpulse  += .01;
 }
 run();
 ncexit();
}


