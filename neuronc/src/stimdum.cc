
/* Program nc */

/* Dummy declarations for stim when compiled separately */

/* in "initchan.cc" and chanxx.cc */

/* Normally stim doesn't need to know about neural elements.  
   However in a complex circuit, it is possible that the arrangement
   of photoreceptors and their stimuli need to be determined from
   information only available from some already-defined neural
   elements.  In this case, comment out the #define XSTIMM below,
   and recompile "stim".
*/

#ifdef XSTIMM		/* define XSTIMM when stim should ignore n elems */
 #define XSTIM		/* define XSTIM for stim */
#endif                    

extern "C" {
#include <stdio.h>
}

#include "nc.h"
#include "ncsub.h"
#include "ncelem.h"
#include "ncomp.h"

#ifdef XSTIM

void initchan(void) {}

void calcratena1m(double v) {}
void calcratena1h(double v) {}

void knrate(double v, int typ, double taun) {}
void khrate(double v, int typ, double tauh) {}

void carate(double v, int typ, double tauc) {}

void maktables(double timestep) {}

double akcacalc(double v, double ca, double tau, double d1, double k1) 
	{ return 0.0;}
double bkcacalc(double v, double ca, double tau, double d2, double k2) 
	{ return 0.0;}

datum    minf(datum &v, datum &elname){}
datum    mtau(datum &v, datum &elname){}
datum    hinf(datum &v, datum &elname){}
datum    htau(datum &v, datum &elname){}
datum    (*ninf) (datum &v, datum &elname);
datum    (*ntau) (datum &v, datum &elname);
datum    (*dinf) (datum &v, datum &elname);
datum    (*dtau) (datum &v, datum &elname);
datum    (*cinf) (datum &v, datum &elname);
datum    (*ctau) (datum &v, datum &elname);
datum    (*chinf) (datum &v, datum &elname);
datum    (*chtau) (datum &v, datum &elname);
datum    ctaua(datum &v, datum &elname){}
datum    ctaub(datum &v, datum &elname){}
datum    ctauc(datum &v, datum &elname){}
datum    ctaud(datum &v, datum &elname){}
datum    ctaue(datum &v, datum &elname){}
datum    ctauf(datum &v, datum &elname){}
datum    (copychan) (datum &sc, datum &ss, datum &dc, datum &ds,
                                        datum &nstates){}
datum    (setchan_ntrans)(datum &sc, datum &ss, datum &state, datum &val){}
datum    (setchan_cond) (datum &sc, datum &ss, datum &state, datum &val){}
datum    (setchan_trans)(datum &sc,datum &ss, datum &state, datum &trans,
                                        datum &val){}
datum    (setchan_trate)(datum &sc,datum &ss, datum &state, datum &trans,
                                        datum &val){}
datum    (setchan_mul) (datum &sc, datum &ss, datum &state, datum &trans,
                                        datum &val){}
datum    (setchan_rateo)(datum &sc,datum &ss, datum &state, datum &trans,
                                        datum &val){}

double gjalpha (double v, gj *gjpnt) {}
double gjbeta  (double v, gj *gjpnt) {}

/* in "chanfunc.cc" */

chantype *getchantype(int ctype, int cnum) {return (chantype *)NULL;}
double calcchaninf(double v, chanparm *chp, double a, double b) {return 0;}
double calcchantau(double v, double tau, chanparm *chp, double a, 
					double b) {return 0;}
double ccavoff(chan *ch){return 0;};
double compcavoff(comp *pnt){return 0;};

/* in "ncomp.cc" */

void dochani(chan *chpnt, double critc) {}
 
/* in "ncplot.cc" */

void mplot(double yval, double xval, int totplots, int plotpos, int pflag) {}
void disperr(void) {}
void initpl(int charind) {}
void plotinit(int totplots) {}
void plotrst(int totplots) {}
void plotpen(int pen, int nplot) {}
void plotchar(int lchar, int mode, int nplot) {}
void plotcsiz(double size, int nplot) {}
void plotname(const char *plname, int nplot) {}
void plotn(int pl, int nplot) {}
void plotsize(double pl, int nplot) {}
void plotval(double pl, int nplot) {}
void plotarr(double *p, int maxindex) {}

const char *rootframe = {"/root"};

/* in "ncsub.cc" */

int nocond=0;                      /* no condense, in ncsub */

void resetelem(void) {}

double record (int cnod1, int cnod2, int cnod3, int cnod4, int pmod)
	 {return 0.0;}

/* in ncconvert */

void condense(void)
{};

/* in "ncnode.cc" */

#ifdef XSTIMM
node *findnode (nodeint num1, nodeint num2, nodeint num3,
		nodeint num4, char *str)
{return ((node *)NULL); }
#endif

/* in "ncrot.cc" */

void setrot(double xrot, double yrot, double zrot,
                double xcent, double ycent, double zcent,
                double rxcent, double rycent, double rzcent, double scal)
{}

/* in "ncdisp.cc" */

void ncdisp(int n1a, int n1b, int n1c, int elemtype, int exceptype,
        int na, int nb, int nc, int color, double dscale, int hide, int excl)
{}

void ncdispc(int n1a, int n1b, int n1c,
        int n2a, int n2b, int n2c, int elemtype, int exceptype,
        int na, int nb, int nc, int color, double dscale, int hide)
{}

void ncdispn(int n1a, int n1b, int n1c,
        int n2a, int n2b, int n2c, int elemtype, int exceptype,
        int na, int nb, int nc, int color, double dscale, int hide, int excl)
{}

void ncdispe (int elemnum, int color, double dscale, int hide)
{}

void drcalib (double x, double y, double length, double size, int color)
{}

void dcomp(int n1a, int n1b, int n1c, int n1d, int elemtype, int exceptype,
	int na, int nb, int nc, int nd,
	double zrange1, double zrange2,
        int color, Symbol *vpen, double (*vpenn)(int elnum, int color),
        double vmax, double vmin,
        double dscale, int hide, int excl)
{}

void makexp(void) 
{}

double qcond(chantype *chp) 
{return 0.0;}

void docacompi(cacomp *capnt)
{}

double qcavoff(chanparm *chp)
{return 0.0;}

int setcmap(datum d1) 
{return 0;}

double reccable(int elemnum, double fracdist)
{ return 0;}

void initray (double xrot, double yrot, double zrot,
        double xcent, double ycent, double zcent,
        double rxcent, double rycent, double rzcent, double scal)
{}

void set_icons(void)
{}

void dispstim(double stime,double dscale, int cmap, double flux_max, double flux_min)
{}

void initcomp(void)
{}

void ncdrnod (int n1a, int n1b, int n1c, int n1d, int elemtype, int exceptype,
	int na, int nb, int nc, int nd,
	double zrange1, double zrange2,
	int color, Symbol *vpen, double (*vpenn)(int elnum, int color),
	double dscale, int drlabel)
{}

void ncdisp(int n1a, int n1b, int n1c, int n1d, int elemtype, int exceptype,
	int na, int nb, int nc, int nd,
	double zrange1, double zrange2,
	int color, Symbol *vpen, double (*vpenn)(int elnum, int color),
	double vmax, double vmin,
	double dscale, int hide, int excl, int cmap)
{}

void ncdispc(int n1a, int n1b, int n1c, int n1d, int n2a, int n2b, int n2c,
	int n2d, int elemtype, int exceptype, int na, int nb, int nc, int nd,
	double zrange1, double zrange2,
	int color, Symbol *vpen, double (*vpenn)(int elnum, int color),
	double vmax, double vmin, double dscale, int hide, int cmap)
{}

void ncdispn(int n1a, int n1b, int n1c, int n1d, int n2a,
	int n2b, int n2c, int n2d, int elemtype, int exceptype,
	int na, int nb, int nc, int nd,
	double zrange1, double zrange2, int color,
	Symbol *vpen, double (*vpenn)(int elnum, int color),
	double dscale, double vmax, double vmin,
	int hide, int excl, int cmap)
{}

void ncdispe (int elemnum, double zrange1, double zrange2,
	int color, Symbol *vpen, double (*vpenn)(int elnum, int color),
	double vmax, double vmin,
	double dscale, int hide, int cmap)
{}

void plotfunc(double (*func)(double, double))
{}

void runsim(double time, int run)
{}

double getcurv(synap *spnt, double val)
{}

double setbind(synap *spnt, double chligand)
{}

double setcond(synap *spnt, double bound)
{return 0;}

double setbind2 (synap *spnt, double active)
{return 0;}

double rsens(photrec *rpnt, double wavel, int filt)
{return 0;}

void runrec(photrec *rpnt, int dtimestep)
{}

chattrib *make_chan (elem *epnt, int ctype, int stype)
{}

elem *foreachn2 (elem *epnt)
{}

elem *setforeachn2 (int etype, int na, int nb, node **ppn)
{}

#endif // XSTIM


