#include "stdplt.h"

void purge(void)
{
	putc(_CMD|'p', stdplt);
	fflush(stdplt);
  }
