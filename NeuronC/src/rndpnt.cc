#include <stdio.h>
#include "drand.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "gr.h"
#include "stdplt.h"


#ifdef __cplusplus
}
#endif

int rseed=2341;
int cumrand=0;

main(void)
{
  double x,y,drand();
  int i;

  setrand (12345);

  stdplt = stdout;
  frame ("zzz");
  cpen (1);
  for (i=0; i<1e8; i++) {

    x = drand();
    y = drand();

//printf ("%g %g\n",x,y);
   move (x,y);
   draw (x,y);
  }  


}

