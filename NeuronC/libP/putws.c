/* putws: put a 16 bit small int to stdout */

#include <stdio.h>

int btest();

void putwsx(int val, FILE *stream)
{
   static int byteswap=0;
   static int bytetest=0;
   static char *hival = 0;
   static short int sval;

   sval = val;
   if (!bytetest) {
      bytetest = 1;
      byteswap = btest();
      if (byteswap) hival = (char *)&sval;
      else          hival = (char *)&sval+1;
   }
   putc(sval,stream);
   putc(*hival,stream); 
}


