#include "stdplt.h"

double cwidth(double width, const char *fname)
               
/*
 * Set character width for textf calls in current frame.
 * The frame named by fname provides the coordinate system
 * used to evaluate width.  The width in the current coordinate
 * system is returned.
 */
  {	register FRAME *dotp, *evalp;
	char *evalfname;
	double x;

	dotp = _dotp;
	evalfname = (char *)fname;
	evalp = _find(&evalfname);
	if (*evalfname) {
		_err("cwidth: %s not found\n", _cfname((char *)fname));
	}
	dotp->_fflags |= _DOTXT;
	dotp->_cw = width*evalp->_sxcat;
	x = 1 / dotp->_sxcat;			/* these 2 lines added to */
	return (x * dotp->_cw);			/*  prevent bug in Optimizer */
/*	return(dotp->_cw/dotp->_sxcat); */
  }
