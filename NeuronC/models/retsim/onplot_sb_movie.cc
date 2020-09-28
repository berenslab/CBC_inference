
/* onplot_sb_movie */


double V_max;
double V_min;
double frame_int;

int show_stim_volts;
int show_poisson_rate;
int show_actual_release;


/*-----------------------------------------------------*/

void defparams_sb_movie(void)

{
  setptr("V_max",     	&V_max);
  setptr("V_min",     	&V_min);
  setptr("frame_int",   &frame_int);

  setptr("show_stim_volts",     &show_stim_volts);
  setptr("show_poisson_rate",   &show_poisson_rate);
  setptr("show_actual_release", &show_actual_release);
}

/*-----------------------------------------------------*/

void onplot_sb_movie_init(void)

{
  if (notinit(frame_int))	frame_int = 0.001;
  if (notinit(V_max))   	V_max     = -0.025;
  if (notinit(V_min))   	V_min     = -0.04;

  if (notinit(show_stim_volts))     show_stim_volts = 0;
  if (notinit(show_poisson_rate))   show_poisson_rate = 0;
  if (notinit(show_actual_release)) show_actual_release = 0;
}

/*-----------------------------------------------------*/

void draw_inputs (void)

{
    int only, color, cmap;
    double vmin, vmax, dscale;
    static int runyet=0;

   if (!runyet) {
        runyet = 1;
        onplot_sb_movie_init();
    }
	/* use cmap=3 (red) for red - blue, cmap=7 for gray scale */
	/* use cmap=4 (red) for excitatory, cmap=6 (blue) for inhib */

  	/* show voltage in presynaptic terminals */

    if (show_stim_volts) {
       display (SPHERE, MATCHING, ndt(cbp,-1,soma), only=1, color=VCOLOR,
                  vmax=-0.04,vmin=-0.045,cmap=7,dscale=0.75);
    }

     /* show poisson release rate for inputs */

    if (show_poisson_rate) {
       display (SPHERE, MATCHING, ndt(cbp,-1,axtrm), only=1, color=SRCOLOR,
                  vmax=200,vmin=0,cmap=3,dscale=5);
    }

	/* show actual vesicle release of inputs */

    if (show_actual_release) {
       display (SPHERE, MATCHING, ndt(cbp,-1,axtrm), only=1, color=SGCOLOR,
                  vmax=5e-11,vmin=0,cmap=4,dscale=5);
    }
}

#include "onplot_movie.cc";

