#include "stdplt.h"

void _move(double x, double y)
              
/*
 * Move pen to (x,y) in current coordinate frame.
 */
  {
	register FRAME *dotp;
	register float *ctp;

	dotp = _dotp;
	ctp = (float *)(dotp->_tcat);
/*        if (x < 0.) x = 0.; */
/*        if (y < 0.) y = 0.; */
	dotp->_x = x;
	dotp->_y = y;
	dotp->_xcat = _T(ctp, x, y) + dotp->_xocat;
	dotp->_ycat = _T(ctp, x, y) + dotp->_yocat;
	dotp->_dphase = 0.;
	dotp->_fflags &= ~_TXTMODE;
  }
