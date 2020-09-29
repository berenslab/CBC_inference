#include "stdplt.h"

void erasea(void)
{
	putc(_CMD|'b', stdplt);	/* next plot abuts this one */
	fflush(stdplt);
	_dotp->_fflags |= _TXTMODE;
	_domove();
  }
