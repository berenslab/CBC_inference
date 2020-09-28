#include "stdplt.h"

extern int _clip(float*,float*,float*,float*,float*);

void _drtri(double x1, double y1, double x2, double y2, double x3, double y3, 
	 int fill)
/*
 * Draw triangle in current coordinate frame.
 */

  {
	register FRAME *dotp;
	register float *ctp;
        float tx1, ty1;
        float tx2, ty2;
        float tx3, ty3;

        dotp = _dotp;
        putc(_CMD| 'v',stdplt);
        putc(fill, stdplt);
        ctp = (float *)(dotp->_tcat);
        tx1 = _T(ctp, x1, y1) + dotp->_xocat;
        ty1 = _T(ctp, x1, y1) + dotp->_yocat;
        ctp = (float *)(dotp->_tcat);
        tx2 = _T(ctp, x2, y2) + dotp->_xocat;
        ty2 = _T(ctp, x2, y2) + dotp->_yocat;
        ctp = (float *)(dotp->_tcat);
        tx3 = _T(ctp, x3, y3) + dotp->_xocat;
        ty3 = _T(ctp, x3, y3) + dotp->_yocat;


        _clip(&tx1, &ty1, &tx2, &ty2, _rwind);
        _clip(&tx2, &ty2, &tx3, &ty3, _rwind);
        _clip(&tx3, &ty3, &tx1, &ty1, _rwind);

        putc(_CMD| 'v',stdplt);
        putwsx((short)tx1,stdplt);
        putwsx((short)ty1,stdplt);
        putc(_CMD| 'v',stdplt);
        putwsx((short)tx2,stdplt);
        putwsx((short)ty2,stdplt);
        putc(_CMD| 'v',stdplt);
        putwsx((short)tx3,stdplt);
        putwsx((short)ty3,stdplt);
        _move (x1,y1);
  }


