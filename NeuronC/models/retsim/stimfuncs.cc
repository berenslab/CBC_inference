/* Module stimfuncs.cc */

/* Functions to generate stimuli */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ncfuncs.h"
#include "retsim.h"
#include "ncio.h"
#include "retsim_var.h"

#define FUG   1	/* Stimulus directions */
#define PET   2
#define BOTH  3


extern double velocity;
extern double sinten;
extern int somaclamp;

int showstim = FUG;		/* when to disp stim */
int stimhalf=0;			/* ->1 for bars stimulate half cell */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/*  use "simwait()" instead of this 

 void wait(double secs)

{
    int i;
    double t,x;

#define MSEC 0.001
#define MSECWAIT 700000

simwait(secs);

//  for (t=0; t<secs; t+= MSEC) {
//      for (i=0; i<MSECWAIT; i++) x = i+t;
//  }
  x += i;
}

*/

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void stimfuncs_init(void)

{
    static int init=0;

  if (init) return;
  init = 1;
  setptrn("showstim",&showstim);
  setptrn("stimhalf",&stimhalf);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

//#define RNDUP 0.0000001
#define RNDUP 0.000000
#define BARLENGTH 1000.0

double movebar(double starttime, double xcent, double ycent, double r1, double r2, 
				double bwidth, double blength, double theta, double velocity, double sinten)

/* Move bar from one X position to another across the neural circuit */
/*  Center of bar is (xcent,ycent), and r1,r2 are radial locations for moving it. */
/*  Orientation of bar is set by "theta". */

/*  Returns the time for the end of the stimulus --  */
/*  Useful for setting the time coordinate on the plot ("endexp"). */

/* Dependent on "velocity", and "sinten" */

/* Note that "t" as defined in this proc is not the simulation time -- */
/*  it is merely a local variable used to set the starting time  */
/*  for the stimulus. */

{
     int i;
     double r, t, rstep, xbar, ybar;
     double coso, sino, timestep, timestep2;
     double orient, inten, start, dur, wavel;

  rstep = 5;			/* move stimulus in increments of 5 um */

  timestep = rstep / velocity;
  sino = -sindeg(theta);
  coso = cosdeg(theta);
  if (r1 <= r2)
  {
    for (t=starttime,r=r1; r<=r2; r+= rstep, t+=timestep)
     {
       xbar = r * coso + xcent;
       ybar = r * sino + ycent;
  
       //fprintf(stderr, "#Adding stim for bar at (%g,%g),theta=%g t=%g dur %g\n",xbar,ybar,theta,t,timestep);
       stim_bar(bwidth,blength, xbar,ybar, orient=theta, inten=sinten, start=t+RNDUP, dur=timestep);
     }
  }
  else if (r1 > r2)
  {
    for (t=starttime,r=r1; r >= r2; r-= rstep, t+=timestep)
    {
       xbar = r * coso + xcent;
       ybar = r * sino + ycent;

       //fprintf(stderr, "#Adding stim for bar at (%g,%g),theta=%g t=%g dur %g\n",xbar,ybar,theta,t,timestep);
       stim_bar(bwidth,blength, xbar,ybar, orient=theta, inten=sinten, start=t+RNDUP, dur=timestep);
    }
  }
  if (disp==1) return(0);
  return t;
}

/*-----------------------------------------------------------------------------*/

double movebar(double starttime, double xcent, double ycent, double r1, double r2, 
				double bwidth, double theta, double velocity, double sinten)

{
   return movebar(starttime, xcent, ycent, r1, r2, bwidth, BARLENGTH, theta, velocity, sinten);
}

/*-----------------------------------------------------------------------------*/

double twospot(double starttime, double xcent, double ycent, double r1, double r2, 
				double dia, double theta, double sinten, double dur, double timestep)

/*  Display two spots, one after another, to simulate motion*/
/*  Center of spot is (xcent,ycent), and r1,r2 are radial locations for moving it. */
/*  Orientation of spot is set by "theta". */

/*  Returns the time for the end of the stimulus --  */
/*  Useful for setting the time coordinate on the plot ("endexp"). */

/* Dependent on "velocity", and "sinten" */

/* Note that "t" as defined in this proc is not the simulation time -- */
/*  it is merely a local variable used to set the starting time  */
/*  for the stimulus. */

{
     double r, t, rstep, xloc, yloc;
     double coso, sino;
     double orient, inten, start;

  rstep = abs(r1-r2);			/* move stimulus from r1 to r2*/
  sino = -sindeg(theta);
  coso = cosdeg(theta);
  if(timestep < 0)
  {
	  double tmp;
	  tmp = r1;
	  r1=r2;
	  r2=tmp;
	  timestep = -timestep;
  }
   // fprintf (stderr,"timestep %g\n",timestep);

       t = starttime;
       r = r1;
       xloc = r * coso + xcent;
       yloc = r * sino + ycent;
  
       // fprintf(stderr, "#1:Adding stim for spot at (%g,%g),theta=%g t=%g dur %g\n",xloc,yloc,theta,t,dur);

       stim_spot (dia, xloc, yloc, sinten, start=t, dur);

  if (rstep>0) {
       t += timestep;
       r = r2;
       xloc = r * coso + xcent;
       yloc = r * sino + ycent;
  
       // fprintf(stderr, "#1:Adding stim for spot at (%g,%g),theta=%g t=%g dur %g\n",xloc,yloc,theta,t,dur);

       stim_spot (dia, xloc, yloc, sinten, start=t, dur);
  }

  if (disp==1) return(0);
  return t+dur;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double stepspot(double starttime,double x1,double x2,double y,double bdia, 
						double velocity, double sinten)

/* Make 2 spots separated in time moving from one X position to another across
   the neural circuit */

/*  Returns the time for the end of the stimulus --  */
/*  Useful for setting the time coordinate on the plot ("endexp"). */

/* Dependent on "velocity", and "sinten" */

/* Note that "t" as defined in this proc is not the simulation time -- */
/*  it is merely a local variable used to set the starting time  */
/*  for the stimulus. */
{
  double x, xstep, start2, start1;
  double spotdur, wait;
  double dscale, inten, start, dur;

  stimfuncs_init();
  xstep = abs (x1-x2);
  wait = xstep / velocity;
  spotdur = bdia /velocity;

  start1=starttime;
  stim_spot(bdia, x1,y, inten=sinten, start=start1, dur=spotdur);
  start2=starttime+wait;

  stim_spot(bdia, x2,y, inten=sinten, start=start2, dur=spotdur);
//  if (disp) {  //default showstim=FUG
//    display_stim (start2,dscale=5);
//    return (0);
//  }

  return starttime+2*spotdur+wait;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double spot_sine(double dia, double x, double y, double tfreq, double minten_add, double minten_mult, double contrast, double starttime, double sdur)

{
     double t, tincr, tend, sinten;

   tincr = 1e-3;
   tend = starttime + sdur;
   for (t=starttime; t<tend; t+=tincr) { 
      sinten = sin(2*PI*t*tfreq);
      stim_spot(dia, x, y, minten_add + minten_mult*sinten*contrast, t, tincr);
      // stim_spot(dia, x, y, minten, t, tincr);
   }
   return starttime+sdur;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double spot_sine(double dia, double x, double y, double tfreq, double minten_mult, double contrast, double starttime, double sdur)

{
    double minten_add;

  return spot_sine(dia, x, y, tfreq, minten_add=0, minten_mult, contrast, starttime, sdur);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double spot_chirp(double dia, double x, double y, double fstart, double fincr, 
		  double minten_add, double minten_mult, double contrast, double starttime, double sdur)

/* make sine-modulated spot with frequency increasing from parameter "fstart" at a rate of "fincr" Hz/sec */

{
     double t, tincr, sinten, f, ft;

   // tincr = 1e-3;
   tincr = 1.0001e-3;
   if (sdur <= 0) sdur = 0.001;
   for (f=fstart, ft=t=0; t<sdur; t+=tincr, ft+=tincr) { 
      sinten = sin(2*PI*t*f);
      f += fincr * ft;			   // increment frequency by fincr 
      stim_spot(dia, x, y, minten_add + minten_mult*sinten*contrast, t+starttime, tincr);
   }
   return starttime+sdur;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double spot_chirp(double dia, double x, double y, double fstart, double fincr, 
		  double minten_mult, double contrast, double starttime, double sdur)

{
     double minten_add;

  return spot_chirp(dia, x, y, fstart, fincr, minten_add=0, minten_mult, contrast, starttime, sdur);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double spot_vcontrast(double dia, double x, double y, double tfreq, double minten_add, double minten_mult, double c1, double c2, double starttime, double sdur)

/* make sine-modulated spot with contrast increasing from parameter c1 to c2 */

{
     double t, tincr, tend, sinten, c, cincr;
     double nsteps;

   // tincr = 1e-3;
   tincr = 1.0001e-3;
   if (sdur <= 0) sdur = 0.001;
   nsteps = sdur / tincr;
   tend = starttime + sdur;
   cincr = (c2 - c1) / nsteps;
   for (c=c1, t=0; t<tend; t+=tincr) {
      sinten = sin(2*PI*t*tfreq);
      stim_spot(dia, x, y, minten_add + minten_mult*sinten*c, t+starttime, tincr);
      c += cincr;				   // increment contrast 
   }
   return starttime+sdur;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double spot_vcontrast(double dia, double x, double y, double tfreq, double minten_mult, double c1, double c2, double starttime, double sdur)

{
     double minten_add;

   return spot_vcontrast(dia, x, y, tfreq, minten_add=0, minten_mult, c1, c2, starttime, sdur);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double moveannulus(double starttime,double xcent, double ycent, double r1,double r2,
			double anndia, double velocity, double sinten)

/* Move annulus from one radius to another across the neural circuit */

/*  Returns the time for the end of the stimulus --  */
/*  Useful for setting the time coordinate on the plot ("endexp"). */

/* Dependent on "velocity", and "sinten" */

/* Note that "t" as defined in this proc is not the simulation time -- */
/*  it is merely a local variable used to set the starting time  */
/*  for the stimulus. */
{
  int cellnr;
  double r, t, xstep, tstep;
  double odia, idia, awidth;
  double inten,start,dur,dscale;

  stimfuncs_init();

  xstep = 5;			/* move stimulus in increments of 5 um */
  tstep = xstep / velocity;

  awidth = anndia;		/* annulus width */

  if (r1 < r2) {
    for (t=starttime,r=r1; r<=r2; r+= xstep, t+=tstep) {
     	 odia = r*2;
         idia = odia - awidth*2;
         if (idia < 0) idia = 0;
         stim_spot(odia, xcent, ycent, inten= sinten, start=t+RNDUP, dur=tstep);
         stim_spot(idia, xcent, ycent, inten=-sinten, start=t+RNDUP, dur=tstep);

//    	 if (showstim==FUG && disp) { //default showstim=FUG
//            display_stim (starttime+0.001, dscale=5);
//            return(0);
//	 };
    };
  }
  else {
    for (t=starttime,r=r1; r >= r2; r-= xstep, t+=tstep) {
      odia = r*2;
      idia = odia - awidth*2;
      if (idia < 0) idia = 0;

      stim_spot(odia, xcent, ycent, inten= sinten, start=t+RNDUP, dur=tstep);
      stim_spot(idia, xcent, ycent, inten=-sinten, start=t+RNDUP, dur=tstep);

//      if (showstim==PET && disp) {
//        display_stim (t,dscale=5);
//        return (0);
//      }
    }
  }
  return t;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double spot_flicker(double dia, double x, double y, double tfreq, int flicker_seed, double minten_add, double minten_mult, double contrast, double starttime, double sdur)

{
     double t, tincr, tend, sinten;

   // tincr = 1e-3;
   if (tfreq<=0) tfreq = 1e-3;
   tincr = 1/tfreq;
   tend = starttime + sdur;
   initrand(1,flicker_seed);
   for (t=starttime; t<tend; t+=tincr) { 
      sinten = rrand(1) * 2.0 - 1.0;
      stim_spot(dia, x, y, minten_add + minten_mult*sinten*contrast, t, tincr);
   }
   return starttime+sdur;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

double spot_flicker(double dia, double x, double y, double tfreq, int flicker_seed, double minten_mult, double contrast, double starttime, double sdur)

{ 
	double minten_add;

  return spot_flicker(dia, x, y, tfreq, flicker_seed, minten_add=0, minten_mult, contrast, starttime, sdur);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void movesineann (double x,double y,int direction,double ann_gaussenv,double centdia,
		double phase,double speriod, double stfreq,double sinten,
		double contrast, int makenv, int sq, double starttime,double sdur)

/* Presents annulus with a moving sine wave grating */
/*  Returns the time for the end of the stimulus --  */

{
  int ct,cn,sdrift;

  if (direction==FUG) sdrift=1;
  else if (direction==PET) sdrift = -1;
  else if (direction==0) sdrift = 0;
  else {
    fprintf(stderr,"# Func movesineann error: direction not valid.\n");
    fprintf(stderr,"# Please enter direction of 0 or 1.\n");
     return;
  }

//  if (!notinit(somaclamp)) {	/* user can specify V to clamp soma */
//    vclamp (ndn(ct=sbac,cn=1,soma), somaclamp, starttime,sdur);
//  }

  /* Make a central gray area with diameter "centdia" */
  if (centdia>0)
    stim_spot(centdia, x, y, sinten, starttime, sdur);

  /* sineann stimulus with x envelope determined by ann_gaussenv */

  stim_sineann (speriod,phase,x,y,stfreq,sdrift,sinten,contrast,
		ann_gaussenv, makenv, sq, starttime,sdur);  /* */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void movesineann (double x,double y,int direction,double ann_gaussenv,double centdia,
		double phase,double speriod, double stfreq,double sinten,
		double contrast, int sq, double starttime,double sdur)

{
   int makenv;

  movesineann (x,y,direction,ann_gaussenv,centdia, phase,speriod, stfreq,sinten,
				contrast, makenv=1, sq, starttime,sdur);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void movewindmill (double x, double y, int direction, double ann_gaussenv,double centdia,
			double phase, double speriod, double stfreq, double sinten_add, 
			double sinten_mult, double scontr, int sq, double starttime, double sdur)

/* Presents moving windmill grating */
/*  Returns the time for the end of the stimulus --  */

{
    int ct,cn;
    double sdrift;

  if (direction==FUG) sdrift=1;
  else if (direction==PET) sdrift = -1;
  else if (direction==0) sdrift = 0;
  else {
    fprintf(stderr,"# Func movewindmill error: direction not valid.\n");
    fprintf(stderr,"# Please enter direction of 0 or 1.\n");
     return;
  }

//  if (!notinit(somaclamp)) {	/* user can specify V to clamp soma */
//    vclamp (ndn(ct=sbac,cn=1,soma), somaclamp, starttime,sdur);
//  }

  /* Make a central gray area with diameter "centdia" */
  if (centdia>0)
    stim_spot(centdia, x, y, sinten_add, starttime, sdur);

  /* windmill stimulus with x envelope determined by ann_gaussenv */

  stim_windmill (speriod,phase,x,y,stfreq, sdrift,sinten_add, sinten_mult,scontr,
		ann_gaussenv, sq, starttime,sdur);  /* */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void movegrating (double x,double y, double orient, double sphase, double speriod, double stfreq, int drift, 
		double sinten_add, double sinten_mult, double contrast, int sq, double start, double dur)

/* Presents moving sine/square grating */
/* sq=0 -> sine; sq=1 -> square */

{
  double wavel, mask, scale, xcent, ycent;

 stim_sine (speriod, sphase, orient, x, y, xcent=0, ycent=0, stfreq, drift, scale=1,
		 sinten_add, sinten_mult, contrast, sq, start, dur, wavel=0, mask=0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void movegrating (double x,double y, double orient, double sphase, double speriod, double stfreq, int drift, 
		double sinten_mult, double contrast, int sq, double start, double dur)

/* Presents moving sine/square grating */
/* sq=0 -> sine; sq=1 -> square */

{
  double wavel, mask, scale, xcent, ycent;
  double sinten_add;
  int stimchan;

 stim_sine (speriod, sphase, orient, x, y, xcent=0, ycent=0, stfreq, drift, scale=1,
		 sinten_mult, contrast, sq, start, dur, wavel=0, mask=0, stimchan=0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void square_wave_i (node *nd, double freq, double i, double start, double dur) 

/* Generate a square wave using current pulses */

{
    int n,ncyc;
    double halfcyc,period;

  if (freq>0) {
     period = 1.0/freq;
     halfcyc = 0.5*period;
     ncyc = dur/period;
     for (n=0; n<ncyc; n++) {
        cclamp (nd, i,start+n*period,halfcyc);
        cclamp (nd,-i,start+n*period+halfcyc,halfcyc);
     }
  } else {
   ncfprintf (stderr,"square_wave_i: zero freq\n");
   return;
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void sine_wave_i (node *nd, double freq, double i, double start, double dur, double tstep) 

/* Generate a sine wave current */

{
    int n,ncyc;
    double t,period,p2,tinc,theta;

  if (freq>0) {
     period = 1.0/freq;
     p2 = 2*PI;
     tinc = p2*tstep/period;
     for (t=theta=0; t<dur; t+=tstep, theta+=tinc) {
	   while (theta>p2) theta-=p2;
           cclamp (nd, i*sin(theta),start+t,tstep);
     }
  } else {
   ncfprintf (stderr,"sine_wave_i: zero freq\n");
   return;
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void ramp_v (node *nd, double vstart, double vstop, double start, double dur, double tstep) 

/* Generate a ramp voltage */

{
    int i, nsteps;
    double t, v, sign, vpulse, vstep;

  if (vstart <= vstop) sign = 1;
  else                 sign = -1;

  if (tstep<=0) tstep = 1e-4;
  nsteps = dur / tstep;
  vstep = (vstop - vstart) / nsteps;
  for (t=0, vpulse=vstart; (vpulse*sign)<=(vstop*sign+1e-6); t+=tstep,vpulse += vstep) {
        vclamp (nd, vpulse, start+t,tstep);
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void ramp_c (node *nd, double cstart, double cstop, double start, double dur, double tstep) 

/* Generate a ramp voltage */

{
    int i, nsteps;
    double t, v, sign, cpulse, cstep;

  if (cstart <= cstop) sign = 1;
  else                 sign = -1;

  if (tstep<=0) tstep = 1e-4;
  nsteps = dur / tstep;
  cstep = (cstop - cstart) / nsteps;
  for (t=0, cpulse=cstart; (cpulse*sign)<=(cstop*sign+1e-6); t+=tstep,cpulse += cstep) {
        cclamp (nd, cpulse, start+t,tstep);
  }
}


