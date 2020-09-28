/* #include "stdplt.h" */
#include "stdpltx.h"

/*
 * Standard Plot Library Data.
 */
FRAME _rootx =
  {	_DOTXT,				/* _fflags */
	&_rootx, 0, 0, "/",		/* _parent, _sibling, _child, _fname */
	0., 0.,				/* _x, _y */
	0., 0.,				/* _xo, _yo */
	0., 0.,				/* _xcat, _ycat */
	0., 0.,				/* _xocat, _yocat */
	STDSIZ, 0.,0., STDSIZ,		/* _t[2][2] */
	STDSIZ, 0.,0., STDSIZ,		/* _tcat[2][2] */
	STDSIZ,				/* _sx */
	STDSIZ,				/* _sxcat */
	STDSIZ,				/* _sy */
	STDSIZ,				/* _sycat */
	0.,				/* _rot */
	0.,				/* _rotcat */
	{0., 0., 0., 0.},		/* _dinfo[8] */
	0.,				/* _dval */
	0,				/* _fat */
	0,				/* _font */
	252.,				/* _cw (.0154*STDSIZ) */
	0.,				/* _crot */
	{0., 0., 0., 0.},		/* _wind[4] */
  };

FRAME *_dotp = &_rootx;
FILE *stdplt;
unsigned _x = 0, _y = 0;
char _fat = 0;
float _rwind[] = {0.,65532.,0.,65532.};   /* 4.*STDSIZ */
int _frlen = sizeof (FRAME);		/* for debugging */
