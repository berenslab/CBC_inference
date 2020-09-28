/* getws: get a 16 bit small int from stdin */

#include <stdio.h>

int btest(void);

int getwsx(stream)
   FILE *stream;

{
   static int byteswap=0;
   static int bytetest=0;
   static short int val=0;
   static char *hival = 0;

   if (!bytetest) {
      bytetest = 1;
      byteswap = btest();
      if (byteswap) hival = (char *)&val;
      else          hival = (char *)&val+1;
   }
   val    = getc(stream);
   *hival = getc(stream); 
   return (val);
}


