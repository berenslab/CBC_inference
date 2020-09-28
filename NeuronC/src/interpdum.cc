/* module interp.cc in program interp */

/* Dummy routines to allow the interpreter from "nc" */
/*  to link OK */

#ifdef __cplusplus
extern "C" {
#endif

#include "stdio.h"

#ifdef __cplusplus
}
#endif

#include "ncplot.h"
#include "nc.h"

plotfr plotnod[PLOTNODSIZ];  /* node numbers, zero sample voltage */
int numplots=0;
char *stfile=0;

FILE *stimout;

int runyet=0;

void initrec() {}		/* in ncmak.cc */ 
void ncleanup() {}
void sortstimfile() {}            /* in ncmain.cc */
void efield() {}            	/* in modcode.cc */
void nfield() {}            	/* in modcode.cc */
void efrac() {}            	/* in modcode.cc */
void edist() {}            	/* in modcode.cc */
void e2dist() {}            	/* in modcode.cc */
void ndist() {}            	/* in modcode.cc */
void vplot() {}            	/* in modcode.cc */
void vplotm() {}            	/* in modcode.cc */
void incrpl() {}            	/* in modcode.cc */
void modrun() {}            	/* in modcode.cc */
void eramod() {}            	/* in modcode.cc */
void xstim() {}            	/* in modcode.cc */
void xrecparm() {}            	/* in modcode.cc */
void xrecept() {}            	/* in modcode.cc */
void noise() {}            	/* in modcode.cc */
void xsynapse() {}            	/* in modcode.cc */
void conn2sl() {}            	/* in modcode.cc */
void conn2dl() {}            	/* in modcode.cc */
void conn2d() {}            	/* in modcode.cc */
void conn2s() {}            	/* in modcode.cc */
void conn1() {}            	/* in modcode.cc */
void conn1m() {}            	/* in modcode.cc */
void conn1l() {}            	/* in modcode.cc */
void membtyp() {}            	/* in modcode.cc */
void xcachan() {}            	/* in modcode.cc */
void xchan() {}            	/* in modcode.cc */
void xsphere() {}            	/* in modcode.cc */
void xcable() {}            	/* in modcode.cc */
void xvbufd() {}            	/* in modcode.cc */
void xvbuf() {}            	/* in modcode.cc */
void xgbatt() {}            	/* in modcode.cc */
void rbatt() {}            	/* in modcode.cc */
void xgcap() {}            	/* in modcode.cc */
void rcap() {}            	/* in modcode.cc */
void xresistor() {}            	/* in modcode.cc */
void rload() {}            	/* in modcode.cc */
void xgj() {}            	/* in modcode.cc */
void xmod() {}            	/* in modcode.cc */
void initrgen() {}		/* in modcode.cc */
void erelem() {}		/* in modcode.cc */
void elimit() {}		/* in modcode.cc */
void e3dist() {}		/* in modcode.cc */
void n2dist() {}		/* in modcode.cc */
void n3dist() {}		/* in modcode.cc */

datum minf(datum&,datum&) {} 	  /* in chanfunc.cc */
datum mtau(datum&,datum&) {}	  /* in chanfunc.cc */
datum hinf(datum&,datum&) {}	  /* in chanfunc.cc */
datum htau(datum&,datum&) {}	  /* in chanfunc.cc */
datum ninf(datum&,datum&) {}  /* in chanfunc.cc */
datum ntau(datum&,datum&) {}  /* in chanfunc.cc */
datum dinf(datum&,datum&) {}  /* in chanfunc.cc */
datum dtau(datum&,datum&) {}  /* in chanfunc.cc */
datum cinf(datum&,datum&) {}  /* in chanfunc.cc */
datum ctau(datum&,datum&) {}  /* in chanfunc.cc */
datum chtau(datum&,datum&) {} /* in chanfunc.cc */
datum chinf(datum&,datum&) {} /* in chanfunc.cc */
datum ctaua(datum&,datum&) {}	/* in chanfunc.cc */
datum ctaub(datum&,datum&) {}	/* in chanfunc.cc */
datum ctauc(datum&,datum&) {}	/* in chanfunc.cc */
datum ctaud(datum&,datum&) {}	/* in chanfunc.cc */
datum ctaue(datum&,datum&) {}	/* in chanfunc.cc */
datum ctauf(datum&,datum&) {}	/* in chanfunc.cc */

void xrecord() {}            	/* in modcode.cc */
void xgausnn() {}            	/* in modcode.cc */
void foreacode() {}            	/* in code.cc */
void dispnod() {}            	/* in modcode.cc */

