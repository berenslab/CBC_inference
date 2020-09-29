/* Module nctest in program nc */

/* Sample procedure to call interpreter from outside program */

/* To use this module, use the code in this (nctest) module as
    a template to run "nc" from your C program. 
*/

#include "nc.h"
#include "ncinit.h"
#include "ncval.h"
#include "y.tab.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif

#ifdef __cplusplus
}
#endif

int ncmain(int argc, char **argv);
double *retvar(char *name, int arg1, int arg2, int arg3, int narg);

extern stype varval[VVSIZE];    /* temp values for setting variables */
extern int varset;               /* index for setting variables */

#define NARGS 2

int main(int argc, char **argv)

{
  static const char *comnd = {"nc"};		/* command name */
  static const char *file  = {"nctest.n"};	/* first command-line argument */
  static char *ncargv[NARGS];			/* dummy argument list for "nc" */
  static double *val;

  ncargv[0] = (char *)comnd;
  ncargv[1] = (char *)file;

  varset = 0;
  varval[varset].name = (char *)"var1";		/* set a variable */
  varval[varset].val  = 55.67;
  varval[varset].type = NUMBER;
  varset++;
  if (varset >= VVSIZE) varset = VVSIZE-1;

  ncmain (NARGS,ncargv);		/* run the simulator */
  val = retvar((char*)"var1",0,0,0,0);		/* get a pointer to a variable */
  printf ("value of var1 is %g\n",*val);
  val = retvar((char*)"var2",0,0,0,1);		/* get a pointer to an array */
  printf ("value of var2 is %g\n",*val);
}

