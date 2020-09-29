#include "stdplt.h"

extern int _clip(float *x1, float *y1, float *x2, float *y2, float *warr);

void _drcirc(double rad, int fill)
             
           

/*
 * Draw circle with center at current location in current coordinate frame.
 */
  {
	register FRAME *dotp;
	register FILE *rstdplt;
	static char cmd[3] = { 0 };
	char *cmdp;		/* ptr to current draw/move command */
	float xnew, ynew;	/* future dotp->_x, _y */
	float xdraw, ydraw;	/* point to draw to in "." coords */
	float xcatnew, ycatnew;	  /* future dotp->_xcat, _ycat */
	float xcatdraw, ycatdraw; /* point to draw to in "" coords */
	float drawrad;

	dotp = _dotp;
	rstdplt = stdplt;
	cmdp = cmd;

	
	if (dotp->_fat != _fat)
	  {	putc(_CMD|'f', rstdplt);
		putwsx(_fat = dotp->_fat, rstdplt);
	  }
	  {	
		_domove();
		drawrad = rad * dotp->_sxcat;
		putc(_CMD | 'k', rstdplt);	/* make circle */
		putwsx((short) drawrad, rstdplt);
		putc((char)fill, rstdplt);	/* =1 -> filled circle */
	  }
  }
