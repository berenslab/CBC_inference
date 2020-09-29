#include "stdplt.h"

void fat(int extra)
            
/*
 * Set line thickness to 1 + 'extra' line widths.
 */
  {	register FRAME *dotp;

	if ((dotp = _dotp) == &_rootx) {
		_err("fat: not valid in root\n");
	}
	dotp->_fat = extra;
  }
