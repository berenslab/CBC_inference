/* module onplot.m */

/* Sets up a spike frequency plot that can be accessed by 
    calling the "spikplot()" procedure: 

 plot_freq = 1;
 if (plot_freq) plot S spikrate max 500 min 0 Char 'o' size .01 vpen spikplot;

*/

if (notinit(plot_freq)) plotfreq=0;

proc onplot () {                /* This procedure runs at plot time */
                                /*  at intervals specified by "ploti". */
  if (plot_freq) spikplot();
  totcur = calc_totcur();
};

if (plot_freq) {

/* Process intracellular voltage to find spike times, */
/* inter-spike intervals, and instantaneous spike */
/* frequencies. */

spikrate = 0;
spikint = 0;
spiktim = 0;
spikthresh = -.04;
spikyet = 0;

vh = 0;
oldvh = 0;
oldspiktim = 0;

func freq_color (nplot, xval, yval)
{
   if (yval > 50) retval = 7
   else if (yval > 20)  retval = 5
   else                 retval = 6;
   if (yval < 1) retval = -1;
   return (retval);
};

proc spikplot ()

/* calculate instantaneous frequency from inter-spike interval */

{
  vh = (V[soma] > spikthresh);

  if (vh && !oldvh) {           /* spike here */
     spiktim = time;
     if (spikyet) {
        spikint = spiktim - oldspiktim + 1e-12;
        spikrate = 1 / spikint;
     }
     else 
	 	spikyet = 1;
     oldspiktim = spiktim;
  }
  else spikrate = 0;
  oldvh = vh;
};

};
