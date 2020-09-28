extern void default_gtext(char *s);
extern void default_gpage(void);
extern void default_gmove(double x, double y);
extern void default_gdraw(double x, double y);
extern void default_grmove(double x, double y);
extern void default_grdraw(double x, double y);
extern void default_gpen(int x);
extern void default_gcrotate(double x);
extern void default_grotate(double x);
extern void default_gorigin(double x, double y);
extern void default_gframe(char *s);
extern void default_grmframe(char *s);
extern void default_gcwidth(double x);
extern void default_gsize(double x);
extern void default_gdash(int x);
extern void default_gcirc(double rad, int fill);
extern void default_grect(double x1, double y1, double x2, double y2, double x3, double y3,
		double x4, double y4, int fill);
extern void default_gtri(double x1, double y1, double x2, double y2, double x3, double y3,
		int fill);
extern void default_gwindow(double x1, double x2, double y1, double y2);
extern void default_gpurge(void);
extern void default_gctext(void);
extern void default_gcgraphics(void);
extern void default_ghinit(void);
extern void default_ghstart(void);
extern void default_ghend(void);
