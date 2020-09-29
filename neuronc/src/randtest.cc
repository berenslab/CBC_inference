
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
int getpid(void);

#ifdef __cplusplus
}
#endif

#include "ncsub.h"
#include "ncomp.h"

int rseed=3332211;
int cumrand=0;

int zrand(int ngen);
void initrand (int ngen, int rseed);
char *setstate(char *arg_state);

double rrand(int ngen);
unsigned long irand(void);

/*****************************************************/

main ()

{
  int i,imax;
  int x1,x2;
  int pfl, pmax;

  pfl = 0;
  pmax = 10000000;
  imax = 2147483647;
  
  initrand (1,12345);
  x1 = zrand(1); 
  
  for (i=0; i<imax; i++) {
    x2 = zrand(1) & 0x7FFFFFFF;
  /*   printf ("%d\n",x2);  /* */
    if (x2==x1) printf ("%d\n",i); 
    if (++pfl == pmax) {
      printf ("."); 
      fflush (stdout);
      pfl = 0;
    }
  }
}

/*-----------------------------------------------------------*/

int zrand(int ngen)
{
  randstate *rpnt;
  int random(void);
  randstate *findrand(int ngen);
  void restorstate();

    return (irand());

  rpnt = findrand(ngen);
  if (rpnt) {
    setstate (rpnt->rstate);
    return (irand());
    restorstate();
  }
  else return 0;
}

