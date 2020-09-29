
/* functions for C++ NC script language */

#ifndef LARGENUM
#include <stdio.h>
#include "nc.h"
#include "ncsub.h"
#endif

/*-----------------------------------------------------------*/

int notinit(double var)

{
   if (var==LARGENUM) return 1;
   else               return 0;
} 

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int notinit(float var)

{
   if (var==LARGENODE) return 1;
   else                return 0;
} 

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int notinit(int var)

{
   if (var==MININT) return 1;
   else             return 0;
} 

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int notinit(const char *var)

{
   if (var==NULL)   return 1;
   else             return 0;
} 

/*-----------------------------------------------------------*/
