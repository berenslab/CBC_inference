/*  Routines in graphics library libP.a */

void cgraphics(void);
void ctext(void);
void cpage(void);
void _move(double x, double y);
void _draw(double x, double y);
void _drcirc(double rad, int fill);
void _drrect(double x1, double y1, double x2, double y2, double x3, double y3,
		double x4, double y4, int fill);
void _drtri(double x1, double y1, double x2, double y2, double x3, double y3,
		int fill);
void _window(double x1, double x2, double y1, double y2, struct _frame *pl);
void _rotate(double x);
void crotate(double x, char *fname);
void _scale(double x, double y);
struct _frame *frame(const char *f);
void rmframe(char *f);
char *gfrname();
void _origin(double x, double y);
double cwidth(double width, const char *fname);
void purge(void);
int textf(const char *fmt, ...);
void cpen (int color);
void dash(int pattern);
void hidinit(void);
void hidstart(void);
void hidstop(void);
