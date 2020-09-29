#include "stdplt.h"

char *strcat(char *, const char *);

char *_cfname(char *fname)
              
/*
 * return ptr to complete framename of `fname' evaluated at _dotp.
 */
  {
	register FRAME *fp;
	char buf[256];

	_path((fp = _find(&fname)), buf);
	if (fp != &_rootx)
		strcat(buf, "/");
	return(strcat(buf, fname));
  }
