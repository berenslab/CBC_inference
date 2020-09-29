/* ph Buffering program */

/* R.G Smith 2013-05 */

#include <stdio.h>
#include <math.h>
#include "ncfuncs.h"
#include "ncinit.h"
#include "gprim.h"

// HCs release atp continuously in the dark

/*  
                  f1
1. atp  + atpase <-> atpx 
                  r1

                  f2
2.          atpx <-> adp + pi + atpase
                  r2

                  f3
3.           atp <-> atpo + h          pka 6.95
                  r3

 -  -  -  -  -  -  -  -  -  -  -  -  

                  f4
4. adp  + atpase <-> adpx 
                  r4

                  f5
5.          adpx <-> amp + pi + atpase
                  r5

                  f6
6.          adp  <-> adpo + h          pka 6.88
                  r6

 -  -  -  -  -  -  -  -  -  -  -  -  

                  f7
7. amp  + atpase <-> ampx 
                  r7

                  f8
8.          ampx <-> ad + pi + atpase
                  r8

 -  -  -  -  -  -  -  -  -  -  -  -  

                  f9
9.            pi <-> pio + h          pka 7.2
                  r9

 -  -  -  -  -  -  -  -  -  -  -  -  

Maybe add bicarbonate buffer with carbonic anhydrase

 -  -  -  -  -  -  -  -  -  -  -  -  

                 dd1
10.            atp -> 0 

                 dd2
11.            adp -> 0 

                 dd3
12.            amp -> 0 

                 dd4
13.             ad -> 0 

                 dd5
14.             h -> 0 

                 dd6
15.             pi -> 0 

at equlibrium: f1 * atp = r1 * adp * pi
           so: f1 = r1 * adp * pi / atp
          and: f1 / r1 = adp * pi / atp = Ka

Ka  = f1/r1
pka = -log10(f1/r1)
exp10(-pka) = f1/r1
r1 = f1/exp10(-pka);

If we assume equilibrium for dissociation of h:

Ka1 = atpo * h / atp = exp10(-6.95)
 h = exp10(-pka1) * atp / atpo

Ka2 = adpo * h / adp = exp10(-6.88)
 h = exp10(-pka2) * adp / adpo

Ka3 = pio * h / pi = exp10(-7.20)
 h = exp10(-pka3) * pi / pio

Ref for pka of atph, adph, pih2:
http://www.life.illinois.edu/crofts/bioph354/atp_hydrolysis.html

*/

/*-----------------------------------------------------*/

void drawlabel (double ypos, int pen, const char *text)
{
  gpen (pen);
  gmove (0.02, ypos);
  gdraw (0.10, ypos);
  gpen (7);
  gmove (0.11, ypos-0.008);
  gtext (text);
};

/*-----------------------------------------------------*/

 main(int argc, char **argv)

{
int i,x;
double t;

double atp, atpx, atpo, adp, adpx, adpo, amp, ampx, ad;        // state variables
double h, pi, pio;                           // state variables

double atpn, atpxn, atpon;		// new estimates
double adpn, adpxn, adpon;		// new estimates
double ampn, ampxn, adn;		// new estimates
double hn, pin, pion;			// new estimates

double dark_stim, light_stim;
double pk1, pk2, pk3;

double f1, r1, f2, r2, f3, r3, f4, r4;
double f5, r5, f6, r6, f7, r7, f8, r8, f9, r9;
double dd1, dd2, dd3, dd4, dd5, dd6;
double atpase, atpasen;

double gtime;
double tstart, tstep, tend;
double hstim;
double iter;
double ymax;
double r, atpz;

double rx1, rx2, rx3, rx4;

ncinit (argc,argv);

endexp = 10;

// command line variables
//
setptrn("tstart",  &tstart);
setptrn("tstep",   &tstep);
setptrn("tend",    &tend);

setptrn("pk1",    &pk1);
setptrn("pk2",    &pk2);
setptrn("pk3",    &pk3);

setptrn("h",     &h);
setptrn("hstim", &hstim);

setptrn("f1",    &f1);
setptrn("f2",    &f2);
setptrn("f3",    &f3);
setptrn("f4",    &f4);
setptrn("f5",    &f5);
setptrn("f6",    &f6);
setptrn("f7",    &f7);
setptrn("f8",    &f8);
setptrn("f9",    &f9);

setptrn("r1",    &r1);
setptrn("r2",    &r2);
setptr ("r3",    &r3);		// overrides pk1 if set on command line
setptrn("r4",    &r4);
setptrn("r5",    &r5);
setptr ("r6",    &r6);		// overrides pk2 if set on command line
setptrn("r7",    &r7);
setptrn("r8",    &r8);
setptr ("r9",    &r9);		// overrides pk3 if set on command line

setptrn("atpase",&atpase);

setptrn("dd1", &dd1);
setptrn("dd2", &dd2);
setptrn("dd3", &dd3);
setptrn("dd4", &dd4);
setptrn("dd5", &dd5);
setptrn("dd6", &dd6);

setptrn("dark_stim",   &dark_stim);
setptrn("light_stim",  &light_stim);

// state variables

atp  = 0.1;
atpx = 0;   // atp bound to atpase 
atpo = 0;   // atp minus h
adp  = 0.1;
adpx = 0;   // adp bound to adpase
adpo = 0;   // adp minus h
amp  = 0.1;
ampx = 0;
ad   = 0;
atpase = 0.5;

h = 1e-7;
pi = 0;
pio = 0;       // pi minus h

dark_stim  = 0.5;
light_stim = 0.0;

pk1 = 6.95;
pk2 = 6.88;
pk3 = 7.20;

//pk1 = 6.8;
//pk2 = 7.0;
//pk3 = 7.5;

 f1 = 10;
 r1 = 0.1;

 f2 = 10;
 r2 = 0.1;

 f3 = 1;

 f4 = 10;
 r4 = 0.1;

 f5 = 10;
 r5 = 0.1;

 f6 = 1;

 f7 = 10;
 r7 = 0.1;

 f8 = 10.0;
 r8 = 0.1;

 f9 = 1;

// decay rates

dd1 = 1;
dd2 = 1;
dd3 = 1;
dd4 = 1;
dd5 = 10;
dd6 = 1;

tstart = -5;
tend = endexp;
tstep = 0.5e-6;

hstim = 0e-3;

x = setvar();    // set variables from command line

// set reverse rates from pkas

if (notinit(r3)) r3 = f3/exp10(-pk1);
if (notinit(r6)) r6 = f6/exp10(-pk2);
if (notinit(r9)) r9 = f9/exp10(-pk3);

gtime = 0;
iter = 1;

ymax = 0.5;

// initialize the graph output

graph_x (endexp, 0);
graph_y (ymax, 0);
graph_y (1, 0);
graph_y (ymax, 0);
graph_y (ymax, 0);
graph_y (ymax, 0);
graph_y (2, 0);
graph_y (10, 0);
//graph_y (35,0);
graph_init();
 
//graph_pen (1, 2, 3, 4, 5, 6, 7);
graph_pen (1, 2, 3, 4, 5, 6);

for (t=tstart; t<tend; t+= tstep) { 

  if (t<=3 || t > 6) { atp += dark_stim * tstep; }
  else               { atp += light_stim * tstep; };

  if (t>=2 && t < 2.001) { h += hstim * tstep; }
  if (t>=5 && t < 5.001) { h += hstim * tstep; }

//  h += hstim * tstep;		// source of protons

  for (i=0; i<iter; i++) {

    atpn = atp + (r1 * atpx 
               - f1 * atp * atpase
               + r3 * atpo * h
               - f3 * atp
               - dd1 * atp ) * tstep;

    rx1 = adp * pi * atpase; 
    atpxn = atpx + (f1 * atp * atpase
		  - r1 * atpx
		  + r2 * rx1
		  - f2 * atpx ) * tstep;

    atpon = atpo + (f3 * atp
                  - r3 * atpo * h) * tstep;

// -  -  -  -  -  -  -  -  -  -  -  -  

    adpn = adp + (f2 * atpx 
                - r2 * rx1 
		+ r4 * adpx 
                - f4 * adp * atpase
                + r6 * adpo * h
                - f6 * adp 
                - dd2 * adp) * tstep;

    rx2 = amp * pi * atpase; 
    adpxn = adpx + (f4 * adp * atpase
		  - r4 * adpx
		  + r5 * rx2
		  - f5 * adpx ) * tstep;

    adpon = adpo + (f6 * adp
                  - r6 * adpo * h) * tstep;

// -  -  -  -  -  -  -  -  -  -  -  -  

    ampn = amp + (f5 * adpx 
                - r5 * rx2 
                + r7 * ampx 
                - f7 * amp * atpase 
                - dd3 * amp) * tstep;
 
    rx3 = ad * pi * atpase; 
    ampxn = ampx + (f7 * amp * atpase
                  - r7 * ampx
                  + r8 * rx3
		  - f8 * ampx ) * tstep;

    adn = ad + (f8 * ampx
	      - r8 * rx3      
	      - dd4 * ad) * tstep;
// -  -  -  -  -  -  -  -  -  -  -  -  

    atpasen = atpase + (r1 * atpx
		      - f1 * atp * atpase
                      + f2 * atpx
		      - r2 * rx1
		      + r4 * adpx
		      - f4 * adp * atpase
		      + f5 * adpx
		      - r5 * rx2
		      + r7 * ampx
		      - f7 * amp * atpase
		      + f8 * ampx
		      - r8 * rx3 ) * tstep;

// -  -  -  -  -  -  -  -  -  -  -  -  

    rx4 = pio * h;
    pin = pi + (f2 * atpx
              - r2 * rx1
              + f5 * adpx
              - r5 * rx2
              + f8 * ampx
              - r8 * rx3 
              + r9 * rx4 
              - f9 * pi ) * tstep;
  
    pion = pio + (f9 * pi
                - r9 * rx4) * tstep;
 
// -  -  -  -  -  -  -  -  -  -  -  -  

    hn  = h  + (f3 * atp 
              - r3 * atpo * h
              + f6 * adp
              - r6 * adpo * h 
              + f9 * pi
              - r9 * rx4 ) * tstep;

    //hn -=  (-log10(h)-dd5) * h  * tstep;
    hn -=  dd5 * h  * tstep;
    pin -= dd6 * pi * tstep;

    // if (adp<0)  adp  = 1e-40;
    // if (amp<0)  amp  = 1e-40;
    // if (h<0)  h  = 1e-40;
    // if (pi<0) pi = 1e-40;

  };  // for (i=0;;)

  atp   = atpn;
  atpx  = atpxn;
  atpo  = atpon;
  adp   = adpn;
  adpx  = adpxn;
  adpo  = adpon;
  amp   = ampn;
  ampx  = ampxn;
  ad    = adn;
  h  = hn;
  pi = pin;
  pio = pion;
  atpase = atpasen;

   if ((atpz=atp)<=0) atpz = 1e-20;
   r = (adp+amp) / atpz;

// printf ("%6.4g %6.4g %6.4g %6.4g %6.4g\n", t, atp, b, hb, h);
// printf ("%6.4g %6.4g %6.4g %6.4g %6.4g %6.4g\n", t, atp, atpase, atpx, hb, r);
// if (t>=0) print "h", h;

   if (t >= 0) {
      if (gtime >= 0.01) {
        // graph (t, atp, adp, amp, -log10(h), r);
        graph (t, atp, atpase, atpx, adp, amp, pi, -log10(h));
        gtime = 0;
      };
      gtime += tstep;
   };
};


drawlabel (0.95,  1, "atp");
drawlabel (0.92,  2, "atpase");
drawlabel (0.89,  3, "atpx");
drawlabel (0.86,  4, "adp");
drawlabel (0.83,  5, "amp");
drawlabel (0.80,  6, "pi");
drawlabel (0.77,  7, "ph");

}
