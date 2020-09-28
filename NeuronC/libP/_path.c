#include "stdplt.h"
#include "string.h"

char * _path(register struct _frame *fp, register char *buf)
                     
                     
/*
 * Place full pathname of frame fp in buf.
 */
  {
	if (fp == &_rootx)
		*buf = '\0';
	else
	  {	_path(fp->_parent, buf);
		if (fp->_parent != &_rootx)
			strcat(buf, "/");
	  }
	return((char *)(strcat(buf,fp->_fname)));
  }
