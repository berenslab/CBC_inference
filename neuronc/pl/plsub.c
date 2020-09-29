/* Routines to extend "libP" library */


#include <stdio.h>
#include "stdplt.h"

#include "../h/mdef.h"

#define PLTCMD 0200

/****************************************/

void cpen(int num)

begin
   static int oldpen = 0;

 if (num == oldpen) return;
 else oldpen = num;
 putc ('c' + PLTCMD,stdout);
 putc (num,stdout);
 fflush(stdout);
end

/****************************************/

void cdash(int num)

begin
 putc ('z' + PLTCMD,stdout);
 putc (num,stdout);
 fflush(stdout);
end

/****************************************/

void hidinit(void)

begin
 putc ('i' + PLTCMD,stdout);
 fflush(stdout);
end

/****************************************/

void hidstart(void)

begin
 putc ('h' + PLTCMD,stdout);
 fflush(stdout);
end


/****************************************/

void hidstop(void)

begin
 putc ('j' + PLTCMD,stdout);
 fflush(stdout);
end

/****************************************/

void ctext(void)

begin
 putc ('a' + PLTCMD,stdout);
 putc ('p' + PLTCMD,stdout);
 fflush(stdout);
end

/****************************************/

void cgraphics(void)

begin
 putc ('g' + PLTCMD,stdout);
 putc ('p' + PLTCMD,stdout);
 fflush(stdout);
end

/****************************************/
