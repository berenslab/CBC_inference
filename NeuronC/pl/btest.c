/* testb: test for byteswapped machine */

#include <stdio.h>

int btest() 

/* returns 0 if byte-swapped from intel machines,
   1 otherwise.
*/

{
   static short intval = 1;
   static char *charval=(char *)&intval + 1;
   int t=0;

/*   printf ("%d\n",*charval); /* */
   return (*charval);
}

/********************************************************/

void read_bs (char *pnum, char *p, int size)

/* Converts number in hardware-independent format to binary number: 
   Format: high byte first. Works for double, float, int, or short int.

   pnum = pointer to binary number.
   p    = pointer to hardware-independent number.
   size = size of number.
 */

{
  int i;

  if (!btest()) {   		/* if from byte-swapped Intel machine */
    for (pnum+=size-1,i=0; i<size; i++) {
      *pnum-- = *p++; 		/* reverse order of bytes */
    }
  }
  else {			/* non byte-swapped */
    switch (size) {

      case 1:   *p             = *pnum; break;
      case 2:   *(short int*)p = *(short int*)pnum; break;
      case 4:   *(int*)p       = *(int*)pnum; break;
      case 8:   *(double*)p    = *(double*)pnum; break;
    }
  }
}

/********************************************************/

void write_bs (char *pnum, char *p, int size)

/* Converts binary number to system-independent format: 
   Format: high byte first. Works for double, float, int, or short int.

   pnum = pointer to binary number.
   p    = pointer to hardware-independent number.
   size = size of number.
 */

{
  int i;

  if (!btest()) {   		/* if from byte-swapped Intel machine */
    for (pnum+=size-1,i=0; i<size; i++) {
      *p++ = *pnum--;		/* reverse order of bytes */
    }
  }
  else {			/* non byte-swapped, just write out */
    switch (size) {

      case 1:   *p             = *pnum; break;
      case 2:   *(short int*)p = *(short int*)pnum; break;
      case 4:   *(int*)p       = *(int*)pnum; break;
      case 8:   *(double*)p    = *(double*)pnum; break;
    }
  }
}
