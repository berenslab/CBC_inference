/* stimfuncs.h */

/* Functions to generate stimuli */

void simwait(double secs);

double movebar(double starttime, double xcent, double ycent, double r1, double r2, 
				double bwidth, double blength, double theta, double velocity, double sinten);

double movebar(double starttime, double xcent, double ycent, double r1, double r2, 
				double bwidth, double theta, double velocity, double sinten);

double stepspot(double starttime,double x1,double x2,double y,double bdia, 
						double velocity, double sinten);

double moveannulus(double starttime,double xcent, double ycent, double r1,double r2,
				double anndia, double velocity, double sinten);

void movesineann (double x,double y,int direction,double ann_gaussenv,double centdia,
		double phase,double speriod, double stfreq,double sinten,
		double contrast,int sq, double starttime,double sdur);

void movesineann (double x,double y,int direction,double ann_gaussenv,double centdia,
		double phase,double speriod, double stfreq,double sinten,
		double contrast, int makenv, int sq, double starttime,double sdur);

void movewindmill (double x, double y, int direction, double ann_gaussenv,double centdia,
			double phase, double speriod, double stfreq, double sinten, double scontr,
			int sq, double starttime, double sdur);

void movegrating (double x,double y, double orient, double sphase, double speriod, double stfreq, int drift,
		                double sinten, double contrast, int sq, double start, double dur);

double twospot(double starttime, double xcent, double ycent, double r1, double r2, 
				double dia, double theta, double sinten, double dur, double timestep);

double spot_sine(double dia, double x, double y, double tfreq, 
				double minten_add, double minten_mult, double contrast, double starttime, double sdur);

double spot_sine(double dia, double x, double y, double tfreq, 
				double minten_mult, double contrast, double starttime, double sdur);

double spot_chirp(double dia, double x, double y, double fstart, double fincr,
		                double minten_add, double minten_mult, double contrast, double starttime, double sdur);

double spot_chirp(double dia, double x, double y, double fstart, double fincr,
		                double minten_mult, double contrast, double starttime, double sdur);

double spot_vcontrast(double dia, double x, double y, double tfreq, double minten_add, 
				double minten_mult, double c1, double c2, double starttime, double sdur);

double spot_vcontrast(double dia, double x, double y, double tfreq, 
				double minten_mult, double c1, double c2, double starttime, double sdur);

double spot_flicker(double dia, double x, double y, double tfreq, int flicker_seed, double minten_add, double minten_mult, double contrast, double starttime, double sdur);

double spot_flicker(double dia, double x, double y, double tfreq, int flicker_seed, double minten_mult, double contrast, double starttime, double sdur);

void square_wave_i (node *nd, double freq, double i, double start, double dur);

void sine_wave_i (node *nd, double freq, double i, double start, double dur, double tstep);

void ramp_v (node *nd, double vstart, double vstop, double start, double dur, double tstep);

void ramp_c (node *nd, double cstart, double cstop, double start, double dur, double tstep);

