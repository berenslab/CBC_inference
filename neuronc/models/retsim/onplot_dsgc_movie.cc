
/* onplot_dsgc_movie */

#include <stdio.h>
 
#ifndef NCFUNCS
 #include "ncfuncs.h"
#endif

#include "retsim.h"
#include "retsim_var.h"

double Vmax;
double Vmin;
double Vmaxg;
double Vming;
double frame_int;

int show_stim_volts;
int show_poisson_rate;
int show_actual_release;

/*-----------------------------------------------------*/

void defparams_dsgc_movie(void)

{
  setptr("Vmax",       &Vmax);
  setptr("Vmin",       &Vmin);
  setptr("Vmaxg",      &Vmaxg);
  setptr("Vming",      &Vming);
  setptr("frame_int",  &frame_int);

  setptr("show_stim_volts",     &show_stim_volts);
  setptr("show_poisson_rate",   &show_poisson_rate);
  setptr("show_actual_release", &show_actual_release);
}

/*-----------------------------------------------------*/

void onplot_dsgc_movie_init(void)

{
  if (notinit(Vmin))           Vmin      = -0.08;     // for color map
  if (notinit(Vmax))           Vmax      =  0.00;     // for color map
  if (notinit(Vming))          Vming     = -0.075;   // for plot of GC voltage
  if (notinit(Vmaxg))          Vmaxg     =  0.00;    // for plot of GC voltage
  if (notinit(frame_int))      frame_int =  0.0002;
 
  if (notinit(show_stim_volts))     show_stim_volts = 0;
  if (notinit(show_poisson_rate))   show_poisson_rate = 0;
  if (notinit(show_actual_release)) show_actual_release = 0;
}

/*-----------------------------------------------------*/

void draw_inputs (void)

{
     int only, color, cmap;
     double vmin, vmax, dscale;

	/* use cmap=3 (red) for red - blue, cmap=7 for gray scale */
	/* use cmap=4 (red) for excitatory, cmap=6 (blue) for inhib */

  	/* show voltage in presynaptic terminals */

    if (show_stim_volts) {
       display (SPHERE, MATCHING, ndt(dbp1,-1,soma), only=1, color=VCOLOR,
		  vmax=-0.04, vmin=-0.045,cmap=4,dscale=0.5);
       display (SPHERE, MATCHING, ndt(ams,-1,soma), only=1, color=VCOLOR,
		  vmax=-0.04,vmin=-0.045, cmap=6, dscale=0.5);
    }

     /* show poisson release rate for inputs */

     if (show_poisson_rate) {
       display (SPHERE, MATCHING, ndt(dbp1,-1,axtrm), only=1, color=SRCOLOR,
  		  vmax=400, vmin=0,cmap=4,dscale=5);
       display (SPHERE, MATCHING, ndt(ams,-1,axtrm), only=1, color=SRCOLOR,
  		  vmax=400,vmin=0,cmap=6,dscale=5);
     }

	/* show actual vesicle release of inputs */

     if (show_actual_release) {
       display (SPHERE, MATCHING, ndt(dbp1,-1,axtrm), only=1, color=SGCOLOR,
  		  vmax=5e-11, vmin=0, cmap=4,dscale=5);
       display (SPHERE, MATCHING, ndt(ams,-1,axtrm), only=1, color=SGCOLOR,
  		  vmax=5e-11, vmin=0, cmap=6,dscale=5);
     }
}

#include "onplot_movie.cc"

