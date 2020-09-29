

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
}
#endif

FILE *ncstdout = stdout;
FILE *ncstderr = stderr;

main(void)
{
   double rate,dur;
   double poisdev(double rate, double dur),poidev(double xm);
   int i;

rate = 1e4; dur = 1e-2;

for (i=0; i<100; i++) {
   fprintf (ncstdout,"poisdev %g pois2 %g\n",poisdev(rate,dur),poidev(rate*dur));
}

}


double poidev (double xm)
{
   static double sq,alxm,g,oldm= -1.0;
   double em,t,y;
   double drand48(void), exp(double);

  if (xm < 1000) {
     if (xm != oldm) {
  	oldm = xm;
	g = exp(-xm);
     }
     em = -1;
     t = 1.0;
     do {
       em += 1.0;
       t *= drand48();
     } while (t > g);
  }
return em;
}

double poisdev(double rate, double dur)
                   

/* return the number of photons received in 
   time "dur", at intensity "rate".
   Normally, "rate" for photons is calibrated in terms of
   quanta / um2 / sec.
*/

{
  static double elaptim,randnum;
  double drand48(void),log(double);
  double count;

#ifdef DEBUG
  if (debug & 1) fprintf (ncstderr,"poisdev\n");
#endif

  if (rate < 1e-30) return (0);

  for (elaptim=0,count=0; elaptim<dur; )
   {
    elaptim -= log(drand48()) / rate;               /* */
    if (elaptim < dur) count++;
   }
 return (count);
}

/*------------------------------------*/

