
/* module plot_freq.m */

/* Sets up a spike frequency plot that can be accessed by calling 
   the "spikplot()" procedure:

   plot_freq = 1;
   if (plot_freq) plot spikplot max 500 min 0 Char 'o' size .01 vpen freq_color;
*/
 
if (plot_freq) {

/* Process intracellular voltage to find spike times, */
/* inter-spike intervals, and instantaneous spike */
/* frequencies. */

spikyet = 0;
oldvh = 0;
oldspiktim = 0;
spikthresh = -.04;

func freq_color (nplot, xval, yval)
{
   if (yval > 50) retval = 15
   else if (yval > 20)  retval = 14
   else                 retval = 12;
   if (yval < 1) retval = -1;
   return (retval);
};

func spikplot (plotnum,xtime)


/* calculate instantaneous frequency from inter-spike interval */

{
    local spikrate, spikint, spiktim;
    local vh;

  vh = (V[gc][soma] > spikthresh);

  if (vh && !oldvh) {           /* spike here */
     spiktim = xtime;
     if (spikyet) {
        spikint = spiktim - oldspiktim + 1e-12;
        spikrate = 1 / spikint;
     }
     else {
        spikyet = 1;
        spikrate = 0;
     };
     oldspiktim = spiktim;
  }
  else spikrate = 0;
  oldvh = vh;
  return spikrate;
 };
};

