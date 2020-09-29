/* tcomp24:  transfer function for synapse */

crit = 1e-8;
endexp = 10;
ploti=.0001;
drm=5000;
exponf = 2;
th= -.04;
dgkd = 1000;

at 1 sphere dia 1;
conn 1 to 2 synapse close expon=exponf cgain=1.001  maxcond 1e-9 thresh=th;
at 2 sphere dia 10;
conn 1 to 3 synapse close expon=exponf cgain=1.05  maxcond 1e-9 thresh=th;
at 3 sphere dia 10;
conn 1 to 4 synapse close expon=exponf cgain=1.5  maxcond 1e-9 thresh=th;
at 4 sphere dia 10;
conn 1 to 5 synapse close expon=exponf cgain=3  maxcond 1e-9 thresh=th;
at 5 sphere dia 10;
conn 1 to 6 synapse open  expon=exponf          maxcond 1e-9 thresh=th;
at 6 sphere dia 10;


/* plot V[1] max .04 min -.07;
plot V[2] max .04 min -.07;
*/

stimlen = .015;

graph X max .0 min -.05;                        /* volts */
graph Y max -.00 min -.07;                        /* psp (volts) (expon)  */
graph Y max -.00 min -.07;                        /* psp (volts) (linear) */
graph Y max -.00 min -.07;                        /* psp (volts) (linear) */
graph Y max -.00 min -.07;                        /* psp (volts) (linear) */
graph Y max -.00 min -.07;                        /* psp (volts) (linear) */
graph Y max -.00 min -.07;                        /* psp (volts) (linear) */

graph init;                                     /* draw axes */

/* stim node 2 vclamp -.07 start 0 dur 1;       /* */

for (i=0,vc=-.045; vc<-.00; vc+=.0001,i++) {
   stim node 1 vclamp vc start i * stimlen dur stimlen; /* */
   step stimlen;                                /* wait for equilbration */
   graph (V[1], V[2], V[3], V[4], V[5], V[6]); /* graph result */

};

