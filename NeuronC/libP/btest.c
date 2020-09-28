/* testb: test for byteswapped machine */

#include <stdio.h>

int btest(void) 

/* returns 1 if byte-swapped from intel machines,
   0 otherwise.
*/

{
   static short intval = 1;
   static char *charval=(char *)&intval + 1;
   int t=0;

/*    printf ("%d\n",*charval); */
   return (*charval);
}

