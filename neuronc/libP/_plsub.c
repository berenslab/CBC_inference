/* Routines to go with "libP" library */


#include <stdio.h>
#include "stdplt.h"

/****************************************/

void cpen(int num)
{
   static int oldpen = -1;

 if (num == oldpen) return;
 else oldpen = num;
 putc ('c' | _CMD,stdplt);
 putwsx (num,stdplt);
 fflush(stdplt);
}

/****************************************/

void hidinit(void)
{
 putc ('i' | _CMD,stdplt);
 fflush(stdplt);
}

/****************************************/

void hidstart(void)
{
 putc ('h' | _CMD,stdplt);
 fflush(stdplt);
}


/****************************************/

void hidstop(void)
{
 putc ('j' | _CMD,stdplt);
 fflush(stdplt);
}

/****************************************/

void cfill( int val)
{
 putc ('o' | _CMD,stdplt);
 putc (val,stdplt);
 fflush(stdplt);
}

/****************************************/

void ctext(void)
{
 putc ('a' | _CMD,stdplt);
 fflush(stdplt);
}

/****************************************/

void cgraphics(void)
{
 putc ('g' | _CMD,stdplt);
 fflush(stdplt);
}

/****************************************/

void cpage(void)
{
 putc ('n' | _CMD,stdplt);
 fflush(stdplt);
}

/****************************************/
