
				/* exp((1000+EXPTHR) / EXPCONST) == 1000 */
				/* exp((1000+60)     / 154)      == 1000 */

/* set up exponential curve for synaptic transfer function.
   starts out at (0,0), goes to (1000,1000). */

{
maxsyn = 2000;
expthr = 60.0;
expcon = 154.0;

 for (i=0; i<maxsyn; i += 10) {
   val = exp((i+expthr)/expcon) - 0.5;
   if (val < 0) val = 0;
   print i, val;
 };
};
