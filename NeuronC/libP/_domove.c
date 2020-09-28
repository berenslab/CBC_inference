#include "stdplt.h"

void _domove(void)
/*
 * Actually move pen to current position if necessary.
 */
  {
	register FRAME *dotp;
	register unsigned xc, yc;

	dotp = _dotp;
	xc = (long) dotp->_xcat;
	yc = (long) dotp->_ycat;
//	if (_x != xc || _y != yc || (dotp->_fflags&_TXTMODE))
	  {	dotp->_fflags &= ~_TXTMODE;
		putc(_CMD|'m', stdplt);
		_x = xc;
		_y = yc;
		putwsx(xc, stdplt);
		putwsx(yc, stdplt);
	  }
  }
