
/* CC file to replace celsegca.n */

#include "retsim.h"

#include "math.h"

double absval(double x){

  if (x < 0) return x * -1;
  return x;

}

double camfunc(double v, int f){

  double x, y;
  double val;

  if (f == ALPHA){
    y = -.09 * (v -63);
    x = exp (y) - 1;

    if (absval(x) > 1e-5)
      val = y / x;
    else 
      val = 1.0;

    val *= 3714.0;
  } else if (f == BETA){
    val = exp((v - 28) / -25.0);
    val *= 2;
  }

  return val*5;

}

double cahfunc(double v, int f){

  double x, y;
  double val;

  if (f == ALPHA){
    val = exp((v + 20) / -9.0);
    val *= .36;
  } else if (f == BETA){
    y = 0.05 * (v + 20);
    x = exp (y) - 1;

    if (absval(x) > 1e-5)
      val = x / y;
    else
      val = 1.0;

    val *= 52.0 / 60.0;
  }
  return val;

}

double calcca5m(double v, int f){
  return camfunc(v, f);
}

double calcca5h(double v, int f){
  return cahfunc(v, f);
}
