/* module gprim.h in program nc */

/* graphics primitives for video and text modes */

void gtext(const char *s);
void glabel(const char *s, double x, double y);
void glabel(const char *s, int pen, double x, double y);
void gpage(void);
void gmove (double x, double y);
void gdraw (double x, double y);
void grmove (double x, double y);
void grdraw (double x, double y);
void gpen (int x);
void gcrotate (double x);
void grotate (double x);
void gorigin(double x, double y);
void gframe (const char *s);
void grmframe (const char *s);
void gcwidth (double x);
void gsize (double x);
void gdash (int x);
void gcirc (double rad, int fill);
void grect (double x1, double y1, double x2, double y2, double x3, double y3,
		double x4, double y4, int fill);
void gtri (double x1, double y1, double x2, double y2, double x3, double y3,
		int fill);
void gwindow (double x1, double x2, double y1, double y2); 
void gpurge (void);
void gctext (void);
void gcgraphics (void);
void ghinit (void);
void ghstart (void);
void ghend (void);
