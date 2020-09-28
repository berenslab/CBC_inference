/* Ccstrip - strip all control characters except \n and \t */

#include <stdio.h>

main()
{
 int c;

 while ((c = getchar()) != EOF)
  if ((c&=0x7f,c >= ' ' && c < 0177) || c == '\t' || c == '\n')
   putchar(c);
  exit (0);
}

