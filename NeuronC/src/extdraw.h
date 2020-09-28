#ifndef _EXTDRAW_H_
#define _EXTDRAW_H_

#include "defaultdraw.h"

/* This header declares function pointers called for non-video mode
 *  drawing methods in gprim.cc. The default implementation of
 *  these functions writes text to standard output, and is defined
 *  in defaultdraw.cc
 */

void (*ext_gtext)(char *) = default_gtext;
void (*ext_gpage)(void) = default_gpage;
void (*ext_gmove)(double x, double y) = default_gmove;
void (*ext_gdraw)(double x, double y) = default_gdraw;
void (*ext_grmove)(double x, double y) = default_grmove;
void (*ext_grdraw)(double x, double y) = default_grdraw;
void (*ext_gpen)(int x) = default_gpen;
void (*ext_gcrotate) (double x) = default_gcrotate;
void (*ext_grotate) (double x) = default_grotate;
void (*ext_gorigin) (double x, double y) = default_gorigin;
void (*ext_gframe) (char *s) = default_gframe;
void (*ext_grmframe) (char *s) = default_grmframe;
void (*ext_gcwidth) (double x) = default_gcwidth;
void (*ext_gsize) (double x) = default_gsize;
void (*ext_gdash) (int x) = default_gdash;
void (*ext_gcirc) (double rad, int fill) = default_gcirc;
void (*ext_grect) (double x1, double y1, double x2, double y2, double x3, double y3,
		double x4, double y4, int fill)  = default_grect;
void (*ext_gtri) (double x1, double y1, double x2, double y2, double x3, double y3,
		int fill)  = default_gtri;
void (*ext_gwindow) (double x1, double x2, double y1, double y2)  = default_gwindow; 
void (*ext_gpurge) (void) = default_gpurge;
void (*ext_gctext) (void) = default_gctext;
void (*ext_gcgraphics) (void) = default_gcgraphics;
void (*ext_ghinit) (void) = default_ghinit;
void (*ext_ghstart) (void) = default_ghstart;
void (*ext_ghend) (void) = default_ghend;

#endif //_EXTDRAW_H_
