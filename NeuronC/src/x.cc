
#include <stdio.h>
#include <math.h>

double set_int_val (double val)

{
	    double logval, tval;

      logval = floor(log10(val));
      tval = int (val / exp(logval*M_LN10) + 0.5);
      if      (tval > 0.8 && tval <= 1.2) tval = 1.0;
      else if (tval > 1.2 && tval <= 3.1) tval = 2.0;
      else if (tval > 3.1 && tval <= 8.0) tval = 5.0;
      else tval = 10.0;
      val = tval * exp(logval*M_LN10);
return val;
}


main(int argc, char **argv)

{
  double i;

  for (i=0.001; i<100; i*= 1.03) {
    printf("%g\n", set_int_val(i));
  }
}

