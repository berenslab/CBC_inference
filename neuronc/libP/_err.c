/* _err.c */

#include <stdarg.h>
#include "stdplt.h"

void _err(const char *fmt, ...)

/* Simulate printf on stderr, followed by abort(). */

{
	va_list argp;

	va_start(argp,fmt);
	fflush(stdplt);
/*	fputs("\n",stderr);		/* force Tektronix to ascii */

/*	_doprnt(fmt, &args, stderr); */
        vfprintf(stderr,fmt,argp);	/*  instead of doprnt */

	/* _cleanup(); */
	abort();			/* cause core dump */
	exit(0177);
}


void _warn(const char *fmt, ...)

/*
 * Simulate printf on stderr.
 */
  {
	va_list argp;

	va_start(argp,fmt);
	fflush(stdplt);
/*	fputs("\n",stderr);		/* force Tektronix to ascii */
/*	_doprnt(fmt, &args, stderr); */	
	vfprintf(stderr,fmt,argp); 	/*  instead of doprnt */
	/* _cleanup(); */
}
