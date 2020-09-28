
/* functions for C++ NC script language */

#ifndef NCSUBH
#include "ncsub.h"
#endif

#ifndef NCOMPH 
#include "ncomp.h"
#endif

#include "ncelem.h"
#include "control.h"
#include "nc.h"
#include "y.tab.h"

#define NCFUNCS 1

// For compile problems with Mac OSX 10.9, leave out scheduler defs.
//  They are only used in morphfuncs.cc, expt_morph_props.cc, and expt_dsgc_calib.cc
//
// #include "scheduler.h"
// extern scheduler sched;

extern elem *elempnt;
extern node *nodepnt;
extern double *rcolors;

extern const char *infile;
extern const char *progname;
void ncinit(int argc, char **argv);
void ncexit(void);
void setexptvar(void);

void setptr (const char *name, int    *iptr);
void setptr (const char *name, float  *fptr);
void setptr (const char *name, double *dptr);
void setptr (const char *name, const char **sptr);

void setptrn (const char *name, int    *iptr);
void setptrn (const char *name, float  *fptr);
void setptrn (const char *name, double *dptr);
void setptrn (const char *name, const char **sptr);

node *nd(int na, int nb, int nc, int nd);
node *nd(int na);
node *nd(int na, int nb);
node *nd(int na, int nb, int nc);
node *ndn(int na, int nb, int nc, int nd);
node *ndn(int na);
node *ndn(int na, int nb);
node *ndn(int na, int nb, int nc);
node *nde(int na, int nb, int nc, int nd);
node *nde(int na);
node *nde(int na, int nb);
node *nde(int na, int nb, int nc);
node *ndt(int na, int nb, int nc, int nd);
node *ndt(int na);
node *ndt(int na, int nb);
node *ndt(int na, int nb, int nc);
node *ndt(void);
node *ndt_free(node *nd);
node *loc(node *npnt,double xloc,double yloc,double zloc);
node *loc(node *npnt,double xloc,double yloc,double zloc,int region,int dendn);
void label(node *npnt, int color);
void label(node *npnt, int color, const char *text);
void label(synapse *spnt, int color);
void label(synapse *spnt, int color, const char *text);

node *findnode(nodeint node1a, nodeint node1b, nodeint node1c, nodeint node1d, const char *s);
node *findnode(nodeint node1a, nodeint node1b, nodeint node1c, const char *s);
node *findnode(nodeint node1a, nodeint node1b, nodeint node1c);
node *findnode(nodeint node1a, nodeint node1b);
node *findnode(nodeint node1a);

elem *at (node *node1, int elemtype);
elem *at (int node1a, int node1b, int node1c, int node1d, int elemtype);
elem *at (int node1, int elemtype);
elem *at (int node1a, int node1b, int elemtype);
elem *at (int node1a, int node1b, int node1c, int elemtype);
void at(elem *epnt, node *pnode, double distfrac, node *nnode);
elem *modify (int elnum);
elem *makelem(int etype, elem *oepnt);
elem *findelem(int elnum);
elem *findelem(int etype, int na, int nb, int nc);
elem *getepnt(nodeint nodenm1, nodeint nodenm2);

elem *conn (node *n1, node *n2, int elemtype);
elem *conn (int node1a, int node1b, int node1c, int node1d,
                 int node2a, int node2b, int node2c, int node2d, int elemtype);
elem *conn (int node1, int node2, int elemtype);
elem *conn (int node1a, int node1b, int node2a, int node2b, int elemtype);
elem *conn (int node1a, int node1b, int node1c,
                 int node2a, int node2b, int node2c, int elemtype);

elem *get_elempnt(int elnum);
elem *get_elempnt(node *npnt, int elnum);
int   elemtype(int elnum);
double get_length(elem *epnt);

node *get_elemnode1 (elem *epnt);
node *get_elemnode2 (elem *epnt);

double get_nfield(node *npnt, int field, int elnum);
double get_nfield(node *npnt, int field);

double get_efield(elem *epnt, int field);
double get_efield(int  elnum, int field);
char *get_elabel(elem  *epnt, int field);
attrib *get_chanattr(elem *epnt);

double dist2d (node *n1, node *n2);
double dist3d (node *n1, node *n2);
double distzd (node *n1, node *n2);

double endist2d(elem *epnt, node *npnt);
double endist3d(elem *epnt, node *npnt);
double endistzd(elem *epnt, node *npnt);
double enfrac(elem *epnt, node *npnt);

double endist2d(int elnum, node *npnt);
double endist3d(int elnum, node *npnt);
double endistzd(int elnum, node *npnt);

electrode *make_electrode (node *node1, node *node2);
electrode *make_electrode (node *node1, node *node2, double r, double c);
electrode *make_electrode (node *node1, node *node2, double r, double c, double ccomp);
loadelem *make_shunt (node *node1, double r);
loadelem *make_shunt (node *node1, double r, double vrev);

cable *make_cable (node *node1, node *node2);
cable *make_cable (node *node1, node *node2, double dia);
cable *make_cable (node *node1, node *node2, double dia, double dia2);
photorec *make_cone (node *node1);
photorec *make_cone (node *node1, int pigm);
photorec *make_cone (node *node1, double dia);
photorec *make_cone (node *node1, int pigm, double dia);
photorec *make_cone (node *node1, double xpos, double ypos);
photorec *make_cone (node *node1, double xpos, double ypos, int stimchan);
photorec *make_rod  (node *node1);
photorec *make_rod  (node *node1, double dia);
photorec *make_rod  (node *node1, int pigm, double dia);
photorec *make_rod  (node *node1, double xpos, double ypos);
photorec *make_rod  (node *node1, double xpos, double ypos, int stimchan);
photorec *make_photorec(int ptype, node *node, double xpos, double ypos);
photorec *make_photorec(int ptype, node *node, double xpos, double ypos, int stimchan);

photorec *make_transducer(node *node, double xpos, double ypos);
photorec *make_transducer(node *node);
photorec *make_itransducer(node *node, double xpos, double ypos);
photorec *make_itransducer(node *node);

sphere *make_sphere (node *node1, double dia);
sphere *make_sphere (node *node1, double dia, double Rm);
sphere *make_sphere (node *node1, double dia, double Rm, double Cm);
synapse *make_synapse (node *node1, node *node2);
synapse *make_synapse (node *node1, node *node2, int connum);
gapjunc *make_gj (node *node1, node *node2, double maxcond);
vbuf *make_vbuf (node *node1, node *node2, double offset, double gain, double tau, double delay, int stype);
vbuf *make_vbuf (node *node1, node *node2, double offset, double gain, double tau, double delay);
nbuf *make_nbuf (node *node1, node *node2, double offset, double gain, double ntoffset, int ntrans);

double *makfiltarr (int newsiz, int oldsiz, double *oarr, double val);
chattrib *chanattr(elem *elpnt, short int ctype);
int get_chan_nstate (int elnum);
double setq10s(void);

void ename (elem *epnt, int *num);
void ename (elem *epnt, double *num);

chan *addchan_extern_ca(comp *pnt, chan *chpnt);
chan *addchan_extern(comp *pnt, chan *chpnt);
chan *addchan_extern_ph(comp *pnt, chan *chpnt);
comp *addcomp_extern(comp *epnt, comp *pnt);
chattrib *make_chan (elem *epnt, int ctype, int stype);
elem *make_chan (node *npnt, int ctype, int stype);
double get_chan_trconc (chattrib *cpnt);
int setnt (comp *cpnt, int type, double val);
double getnt (comp *cpnt, int type);
int setntfromcomp (comp *cpnt1, comp *cpnt2, double offset, double gain, int type);

nattrib *make_chnoise  (elem *epnt);
nattrib *make_vesnoise (elem *epnt);

void chset(elem *epnt);
void set_chancalc (int ctype, int stype, int np, double (*func)(double, int));

double *make_filt(double val);
double *make_filt(double val1,double val2);
double *make_filt(double val1,double val2,double val3);
double *make_filt(double val1,double val2,double val3,double val4);
double *make_filt(double val1,double val2,double val3,double val4,double val5);

int gausnn (double **xarr, double **yarr, int ncells, double density, int rsd, double reg,
                double xcent, double ycent, int ginfo);
int gausnn (double **xarr, double **yarr, double density, int rsd, double reg,
                double xsize, double ysize, double xcent, double ycent, int ginfo);
int gausnn (double **xarr, double **yarr, int ncells, int rsd, double reg,
                double xsize, double ysize, double xcent, double ycent, int ginfo);
int gausnn (double **xarr, double **yarr, double nnd, double reg, int rsd, 
                double xsize, double ysize, double xcent, double ycent, int ginfo);

int gausnn (double **xarr, double **yarr, int ncells, double density, int rsd, double reg,
                double xcent, double ycent, int first_center, int ginfo);
int gausnn (double **xarr, double **yarr, double density, int rsd, double reg,
                double xsize, double ysize, double xcent, double ycent, int first_center, int ginfo);
int gausnn (double **xarr, double **yarr, int ncells, int rsd, double reg,
                double xsize, double ysize, double xcent, double ycent, int first_center, int ginfo);
int gausnn (double **xarr, double **yarr, double nnd, double reg, int rsd, 
                double xsize, double ysize, double xcent, double ycent, int first_center, int ginfo);

void stim_node (node *npnt, double inten, double start, double dur, double wavel);
void stim_node (node *npnt, double inten, double start, double dur);
void stim_cone (node *npnt, double inten, double start, double dur, double wavel);
void stim_cone (node *npnt, double inten, double start, double dur);
void stim_rod  (node *npnt, double inten, double start, double dur, double wavel);
void stim_rod  (node *npnt, double inten, double start, double dur);

void stim_bar (double width, double length, double xloc, double yloc,
               double xcent, double ycent, double scale, double orient,
               double inten, double start, double dur, double wavel, double mask, int stimchan);

void stim_bar (double width, double length, double xloc, double yloc, double orient,
          double inten, double start, double dur, double wavel, double mask, int stimchan);

void stim_bar (double width, double length, double xloc, double yloc, double orient,
          double inten, double start, double dur, double mask, int stimchan);

void stim_bar (double width, double length, double xloc, double yloc, double orient,
          double inten, double start, double dur, double mask);

void stim_bar (double width, double length, double xloc, double yloc, double orient,
          double inten, double start, double dur);

void stim_spot (double dia, double xloc, double yloc, double xcent, double ycent,
                double scale, double inten, double start, double dur, double wavel, double mask, int stimchan);

void stim_spot (double dia, double xloc, double yloc, double xcent, double ycent,
                double scale, double inten, double start, double dur, double wavel, double mask);

void stim_spot (double dia, double xloc, double yloc, double inten,
                        double start, double dur, double wavel, double mask, int stimchan);

void stim_spot (double dia, double xloc, double yloc, double inten,
                        double start, double dur, double wavel, double mask);

void stim_spot (double dia, double xloc, double yloc, double inten, double start, double dur, double mask);

void stim_spot (double dia, double xloc, double yloc, double inten, double start, double dur);

void stim_annulus (double idia, double odia, double xloc, double yloc, 
		   double inten, double start, double dur, double wavel, double mask);

void stim_annulus (double idia, double odia, double xloc, double yloc, 
		   double inten, double start, double dur);

void stim_grating (int type, double speriod, double sphase, double orient,
		double xloc, double yloc, double tfreq, double drift,
		double inten, double contrast, double wavel,
		double xenv, double yenv, int sq, double mask, int stimchan,
		double start, double dur);

void stim_sine(double speriod, double sphase, double orient,
		double xloc, double yloc, double xcent, double ycent,
		double tfreq, int drift, double scale,
		double inten_add, double inten_mult, double contrast, int sq, 
		double start, double dur, double wavel, double mask, int stimchan);

void stim_sine(double speriod, double sphase, double orient,
		double xloc, double yloc, double xcent, double ycent,
		double tfreq, int drift, double scale,
		double inten_mult, double contrast, int sq, double start, double dur,
		double wavel, double mask, int stimchan);

void stim_sine(double speriod, double sphase, double orient, double
		xloc, double yloc, double tfreq, int drift,
		double inten_add, double inten_mult, double contrast, int sq, double start, double dur);
		
void stim_sine(double speriod, double sphase, double orient, double
		xloc, double yloc, double tfreq, int drift,
		double inten_mult, double contrast, int sq, double start, double dur);
		
void stim_gabor(double speriod, double sphase, double orient,
		double xloc, double yloc, double xcent, double ycent,
		double tfreq, double drift, double scale,
		double inten_add, double inten_mult, double contrast, 
		double xenv, double yenv, int sq,
		double start, double dur, double wavel, double mask, int stimchan);

void stim_gabor(double speriod, double sphase, double orient,
		double xloc, double yloc, double xcent, double ycent,
		double tfreq, double drift, double scale,
		double inten_mult, double contrast, double xenv, double yenv, int sq,
		double start, double dur, double wavel, double mask, int stimchan);

void stim_gabor(double speriod, double sphase, double orient, double xloc, double yloc,
		double tfreq, double drift, double inten_add, double inten_mult, double contrast,
		double xenv, double yenv, int sq, double start, double dur);
		
void stim_gabor(double speriod, double sphase, double orient, double xloc, double yloc,
		double tfreq, double drift, double inten_mult, double contrast,
		double xenv, double yenv, int sq, double start, double dur);
		
void stim_sineann(double speriod, double sphase, double xloc, double yloc,
		double xcent, double ycent, double tfreq, double drift, double scale,
		double inten_add, double inten_mult, double contrast, double xenv, int makenv, int sq,
		double start, double dur, double wavel, double mask, int stimchan);

void stim_sineann(double speriod, double sphase, double xloc, double yloc,
		double xcent, double ycent, double tfreq, double drift, double scale,
		double inten_mult, double contrast, double xenv, int makenv, int sq,
		double start, double dur, double wavel, double mask, int stimchan);

void stim_sineann(double speriod, double sphase, double xloc, double yloc,
		double xcent, double ycent, double tfreq, double drift, double scale,
		double inten_mult, double contrast, double xenv, int sq,
		double start, double dur, double wavel, double mask, int stimchan);

void stim_sineann(double speriod, double sphase, double xloc, double yloc,
		double tfreq, double drift, double inten, double contrast,
		double xenv, double start, double dur, int stimchan);

void stim_sineann(double speriod, double sphase, double xloc, double yloc,
		double tfreq, double drift, double inten_add, double inten_mult, 
		double contrast, double xenv, int sq, double start, double dur);

void stim_sineann(double speriod, double sphase, double xloc, double yloc,
		double tfreq, double drift, double inten_add, double inten_mult, 
		double contrast, double xenv, int makenv, int sq, double start, double dur);

void stim_sineann(double speriod, double sphase, double xloc, double yloc,
		double tfreq, double drift, double inten_mult, double contrast,
		double renv, int sq, double start, double dur);

void stim_sineann(double speriod, double sphase, double xloc, double yloc,
		double tfreq, double drift, double inten_mult, double contrast,
		double renv, int makenv, int sq, double start, double dur);

void stim_windmill(double speriod, double sphase, double xloc, double yloc,
		double xcent, double ycent, double tfreq, double drift, 
		double scale, double inten_add, double inten_mult, double contrast, 
		double renv, int makenv, int sq, double start, double dur, double wavel, double mask, int stimchan);

void stim_windmill(double speriod, double sphase, double xloc, double yloc,
		double xcent, double ycent, double tfreq, double drift, 
		double scale, double inten_mult, double contrast, double xenv, 
		int sq, double start, double dur, double wavel, double mask, int stimchan);

void stim_windmill(double speriod, double sphase, double xloc, double yloc,
		double tfreq, double drift, double inten, double contrast,
		double renv, double start, double dur);

void stim_windmill(double speriod, double sphase, double xloc, double yloc,
		double tfreq, double drift, double inten_add, double inten_mult, double contrast,
		double renv, int sq, double start, double dur);

void stim_windmill(double speriod, double sphase, double xloc, double yloc,
		double tfreq, double drift, double inten_mult, double contrast,
		double renv, int sq, double start, double dur);

void stim_checkerboard(double width, double height, int xn, int yn,
		double orient, double xloc, double yloc,
		double xcent, double ycent, double scale,
		double tfreq, double inten, double contrast,
		double start, double dur, double **stim_rndarr, int *stim_nfr, int chkrseed);

void stim_checkerboard(double width, double height, int xn, int yn, 
		double orient, double xloc, double yloc,
		double tfreq, double inten, double contrast,
		double start, double dur,double **stim_rndarr, int *stim_nfr, int chkrseed);

void stim_checkerboard(double width, double height, int xn, int yn, 
		double orient, double xloc, double yloc,
		double tfreq, double inten, double contrast,
		double start, double dur, int chkrseed);

void stim_checkerboard(double width, double height, int xn, int yn, 
		double orient, double xloc, double yloc,
		double tfreq, double inten, double contrast,
		double start, double dur);

void sector_mask(double width, double orient, double val, double stimtime, double dur, int stimchan);
void sector_mask(double width, double orient, double val, double stimtime, double dur);

void stim_file (const char *filename);
void comp_file (const char *filename);
void stim_backgr (double backgr, double wavel, double mask, int stimchan, double start);
void stim_backgr (double backgr, double wavel, double mask, double start);
void stim_backgr (double backgr, int stimchan, double start);
void stim_backgr (double backgr, double start);
void stim_backgr (double backgr, int stimchan);
void stim_backgr (double backgr);
void stim_blur (double blurrad, double scatter_ampl, double scatter_rad,
		                        double scatter_pow, double sscale);
void stim_blur (double blurrad, double blur_ampl, double scatter_ampl, double scatter_rad,
		                        double scatter_pow, double sscale);

void vclamp (node *npnt, double inten, double start, double dur);
void cclamp (node *npnt, double inten, double start, double dur);
void puff (node *npnt, int puffmsg, double inten, double start, double dur);

int findphotrec (int nodenm1, int nodenm2, int nodenm3, photrec **rpnt, const char *str);
int findphotrec (int nodenm1, int nodenm2, int nodenm3, int nodenm4, photrec **rpnt, const char *str);

void plot (int param, int ename);
void plot (int param, int cagparm, int ename);
void plot (int param, int cagparm, int ename, double max, double min);
void plot (int param, int ename, double max, double min);
void plot (int param, int ename, double frac);
void plot (int param, node *npnt);
void plot (int param, node *npnt, double max, double min);
void plot (int param, int cagparm, node *npnt);
void plot (int param, int cagparm, node *npnt, double max, double min);
void plot_func (double(*func)(double val, double t), double plval, double max, double min);
void plot_func (double(*func)(double val, double t), double plval, double max, double min, int pmod, int pmod2);
void plot_func (double(*func)(double val, double val2, double t), double plval, double plval2, double max, double min);
void plot_func (double(*func)(double val, double val2, double t), double plval, double plval2, double max, double min, int pmod, int pmod2);
void plot_func (double(*func)(double val, double val2, double val3, double t), double plval, double plval2, double plval3, double max, double min);
void plot_func (double(*func)(double val, double val2, double val3, double t), double plval, double plval2, double plval3, double max, double min, int pmod, int pmod2);

void plot_func (double(*func)(double val, double val2, double val3, double val4, double t), double plval, double plval2, double plval3, double plval4, double max, double min);
void plot_func (double(*func)(double val, double val2, double val3, double val4, double t), double plval, double plval2, double plval3, double plval4, double max, double min, int pmod, int pmod2);

void plot_func (double(*func)(double val, double val2, double val3, double val4, double val5, double t), double plval, double plval2, double plval3, double plval4, double plval5, double max, double min);
void plot_func (double(*func)(double val, double val2, double val3, double val4, double val5, double t), double plval, double plval2, double plval3, double plval4, double plval5, double max, double min, int pmod, int pmod2);

void plot_func (double(*func)(double val, double val2, double val3, double val4, double val5, double val6, double t), double plval, double plval2, double plval3, double plval4, double plval5, double plval6, double max, double min);
void plot_func (double(*func)(double val, double val2, double val3, double val4, double val5, double val6, double t), double plval, double plval2, double plval3, double plval4, double plval5, double plval6, double max, double min, int pmod);

void plot_func (double(*func)(double val, double val2, double val3, double val4, double val5, double val6, double val7, double t), double plval, double plval2, double plval3, double plval4, double plval5, double plval6, double plval7, double max, double min);

void plot_func (double(*func)(double val, double val2, double val3, double val4, double val5, double val6, double val7, double t), double plval, double plval2, double plval3, double plval4, double plval5, double plval6, double plval7, double max, double min, int pmod);

void plot_funcg (double(*func)(double val, double t), double plval, double max, double min);
void plot_funci (double(*func)(double val, double t), double plval, double max, double min);

void plot_var (double *var, double plval, double max, double min);
void plot_var (double *var, double plval, double max, double min, int pmod);

void plot_arr (double *var, int maxindex);

void plot_param (const char *name, int plnum, double plsize);
void plot_param (const char *name, int pen, int plnum, double plsize);
void plot_param (const char *name, int pen, int plnum, double plsize, int pmod);
void plot_param (const char *name, int pen, int plnum);
void plot_param (const char *name, int pen, int plnum, double plsize, double plval);
void plot_pen (int pen);
void plot_vpen (double (*vpen)(int, double, double));

void plotinit(int numplots);
void plotpen (int val, int i);
void plotchar (int val,int lines, int i);
void plotvpenc  (double (vpen)(int,double,double));
void plotname (const char *plname);
void plotsize (double plotsiz);
void plotmode (int pmod);
void plotfilt (int nfilt, double *timec, int nplot);
void plotfilt (int nfilt, double *timec);
void plotfilt (int nfilt, double *timec, double initval);
void plotfilt (int nfilt, double *timec, double initval, int plotnum);
void plotbessfilt(double cutoff);

void set_synapse_dr (void (*synapse_drpnt)(synapse *,int,double,double,double,double,double,int));
// void set_synapse_dr (void (*synapse_drpnt)(int,double,double,double,double,double,int));
void set_phot_dr (void (*phot_drpnt)(int,int,int,double,double,double,int));
void set_elec_dr (void (*elec_drpnt)(int,double,double,double,double,int,int));
void set_gapjunc_dr (void (*gapjunc_drpnt)(int,double,double,double,int));
void dr_node (node *npnt, double dscale);

void graph(double x, double y1);
void graph(double x, double y1, double y2);
void graph(double x, double y1, double y2, double y3);
void graph(double x, double y1, double y2, double y3, double y4);
void graph(double x, double y1, double y2, double y3, double y4, double y5);
void graph(double x, double y1, double y2, double y3, double y4, double y5, double y6);
void graph(double x, double y1, double y2, double y3, double y4, double y5, double y6, double y7);
void graph(double x, double y1, double y2, double y3, double y4, double y5, double y6, double y7, double y8);
void graph(double x, double *y, int n);

void graph_x (double xmax, double xmin);
void graph_y (double ymax, double ymin);
void graph_init (void);
void graph_restart(void);

void graph_pen (int *val, int n);
void graph_pen (int pen1);
void graph_pen (int pen1, int pen2);
void graph_pen (int pen1, int pen2, int pen3);
void graph_pen (int pen1, int pen2, int pen3, int pen4);
void graph_pen (int pen1, int pen2, int pen3, int pen4, int pen5);
void graph_pen (int pen1, int pen2, int pen3, int pen4, int pen5, int pen6);
void graph_pen (int pen1, int pen2, int pen3, int pen4, int pen5, int pen6, int pen7);
void graph_pen (int pen1, int pen2, int pen3, int pen4, int pen5, int pen6, int pen7, int pen8);

void graph_char (char val);
void graph_cchar (char val);
void graph_csiz (double val);
void graph_set (const char *name, int plnum, double plsize, double plval);
void graph_set (const char *name, int plnum, double plsize);
void graph_filt (double *val, int n);
void graph_vpen (double (*vpen)(int,double,double));

double v(node *npnt);
double v(int node1a, int node1b, int node1c, int node1d);
double v(int node1a, int node1b, int node1c);
double v(int node1a, int node1b);
double v(int node1a);
double v(int elemnum, double fracdist);

double i(node *npnt);
double i(int node1a, int node1b, int node1c, int node1d);
double i(int node1a, int node1b, int node1c);
double i(int node1a, int node1b);
double i(int node1a);
double l(node *npnt);

double record_nt (node *npnt, int nt);
double record_ca (node *npnt, int caval);
double record_ca (node *npnt, int caval, int pval);
double record_chan (int elnum, int nf);
double record_chan (int elnum, int nf, int pval);
double record_chan (elem *epnt, int nf, int pval);
double record_chan (chan *chpnt, int nf, int pval);
double record_synapse(int snum, int nf);
double record_synapse(synapse *spnt, int nf);
double record (int cnod1, int cnod2, int cnod3, int cnod4, int pmod, int pval);

elem *foreach (elem *epnt, int etype);
elem *foreach (elem *epnt, int etype, double radius, 
					double (*distfunc)(elem *e1, node *n2), 
					node *npntw);

elem *foreach (elem *epnt, int etype, int na1, int na2, int na3, int na4,
				      int nb1, int nb2, int nb3, int nb4, int xxx);
elem *foreach (elem *epnt, int etype, int na1, int na2, int na3, 
				      int nb1, int nb2, int nb3, int xxx);
elem *foreach (elem *epnt, int etype, int na1, int na2, int nb1, int nb2, int xxx);

elem *foreach (elem *epnt, int etype, int na, int nb, int nc, int nd);
elem *foreach (elem *epnt, int etype, int na, int nb, int nc);
elem *foreach (elem *epnt, int etype, int na, int nb);

elem *foreach (elem *epnt, int etype, int na, int nb, int nc, int nd,
                                      int *pa, int *pb, int *pc, int *pd);
elem *foreach (elem *epnt, int etype, int na, int nb, int nc,
                                      int *pa, int *pb, int *pc);
elem *foreach (elem *epnt, int etype, int na, int nb,
                                      int *pa, int *pb);

node *getnpnt (nodeint nodenm1, nodeint nodenm2);
elem *foreachn (elem *epnt, int etype, int na, int nb, node **ppn);
elem *foreachn2 (elem *epnt);
elem *setforeachn2 (int etype, int na, int nb, node **ppn);

elem *foreach (elem *epnt, int etype, int na, int nb, int nc, int nd,
                                      int *pa, int *pb, int *pc, int *pd,
                                        double radius, double (*distfunc)(node *n1, node *n2), 
					node *npntw);

elem *foreachn (elem *epnt, int etype, int na, int nb, 
                                        double radius, double (*distfunc)(node *n1, node *n2), 
					node *npntw, node **ppn);
elem *setforeachn3 (int etype, int na, int nb, 
                     double radius, double (*distfunc)(node *n1, node *n2), 
					node *npntw, node **ppn);
elem *foreachn3 (elem *epnt);

elem *foreach (elem *epnt, int etype, int na, int nb, int nc, int nd,
                                      int *pa, int *pb, int *pc, int *pd,
                                        double radius, double (*distfunc)(elem *e1, node *n2), 
					node *npntw);

node *foreach (node *npnt, int na, int nb, int nc, int nd,
                           int *pa, int *pb, int *pc, int *pd);

node *foreach (node *npnt, int na, int nb, int nc, int *pa, int *pb, int *pc);
node *foreach (node *npnt, int na,  int  nb, int *pa, int *pb);

node *foreach (node *npnt, int na, int nb, int nc, int nd,
                           int *pa, int *pb, int *pc, int *pd,
                           double radius, double (*distfunc)(node *n1, node *n2), node *npntw);

node *foreach (node *npnt, int na, int nb, int nc, 
                           int *pa, int *pb, int *pc, 
                           double radius, double (*distfunc)(node *n1, node *n2), node *npntw);

node *foreach (node *npnt, int na, int nb, int nc, int nd);
node *foreach (node *npnt, int na, int nb, int nc);
node *foreach (node *npnt, int na, int nb);

node *foreachn (node *npnt, int na, int nb);
node *foreachn (node *npnt, int na, int nb, int *pc);

chan *findchan (conlst *cpnt, int ctype, int stype);

void set_disp_rot( double xrot,  double yrot,  double zrot,
             double dxcent, double dycent, double dzcent,
             double rxcent, double rycent, double rzcent, double dsize);
void display_rot( double xrot,  double yrot,  double zrot);
void display_center( double dxcent,  double dycent,  double dzcent);

void display_size (double size);

void display               (int disptype, node *nd1, int dcolor);
void display (int elemtype, int disptype, node *nd1, int dcolor);

void display               (int disptype, node *nd1, int dcolor, double dscale);
void display (int elemtype, int disptype, node *nd1, int dcolor, double dscale);

void display               (int disptype, node *nd1, int exceptype, node *nexc, int dcolor, double dscale);
void display (int elemtype, int disptype, node *nd1, int exceptype, node *next, int dcolor, double dscale);

void display               (int disptype, node *nd1, double (*vpenn)(int elnum, int dcolor));
void display (int elemtype, int disptype, node *nd1, double (*vpenn)(int elnum, int dcolor));

void display               (int disptype, node *nd1, int dcolor, double vmax, double vmin);
void display (int elemtype, int disptype, node *nd1, int dcolor, double vmax, double vmin);

void display               (int disptype, node *nd1, int excl, int dcolor,
		                        double vmax, double vmin, int cmap, double dscale);

void display (int elemtype, int disptype, node *nd1, int excl, int dcolor,  
					double vmax, double vmin, int cmap, double dscale);

void display                (int disptype, node *nd1, int dcolor, double(*vpenn)(int elnum, int color), 
					int cmap, double vmax, double vmin,
             				double dscale, int hide, int excl, double stime);

void display (int elemtype, int disptype, node *nd1, int dcolor, double(*vpenn)(int elnum, int color),
             				int cmap, double vmax, double vmin,
             				double dscale, int hide, int excl, double stime);

void display               (int disptype, node *nd1, node *nd2, node *nexc, int exceptype,
             				int dcolor, double(*vpenn)(int, int),
             				int cmap, double vmax, double vmin,
             				double dscale, int hide, int excl, double stime);
void display (int elemtype, int disptype, node *nd1, node *nd2, node *nexc, int exceptype,
             				int dcolor, double(*vpenn)(int, int),
             				int cmap, double vmax, double vmin,
             				double dscale, int hide, int excl, double stime);
void display (int elemtype, int elnum, int dcolor, double(*vpenn)(int, int),
	     				int cmap, double vmax, double vmin,
	     				double dscale, int hide);
void display (int elemtype, int elnum, int dcolor, double dscale);

void display (elem *epnt, int dcolor, double dscale);

void display_z (double z2, double z1);

void display_stim (double sttime, double attime, double dscale, int cmap, double max_flux, double min_flux);
void display_stim (double sttime, double attime, double dscale, double max_flux, double min_flux);
void display_stim (double sttime, double attime, double dscale);

void display_stim (double attime, double dscale, double max_flux, double min_flux);
void display_stim (double attime, double dscale);

void display_page(void);

void set_phot_dr (void (*phot_drpnt)(int,int,int,double,double,double,int));

void drcalib    (double xcalib,double ycalib,double cline,double dsize,int dcolor);
void disp_calib (double xcalib,double ycalib,double cline,double dsize,int dcolor);
void display_calib (double cline,int dcolor);

int setcmap (int *parr, int cmapsize);
int setcmap (int cmap);

void elimit (int param, double emax, double emin);
void elimit (int elnum);
void elimit (elem *elnum);

node *erasenode(node *npnt);
elem *eraseelem(elem *epnt);
node *erase(node *npnt);
elem *erase(elem *epnt);

int setvar(void);
int setvars(int argc, char **argv);
double getval(const char *var);

int notinit(double var);
int notinit(float var);
int notinit(int   var);
int notinit(const char *var);
double *fread(const char *str, int *nlong, int *nwid);
char *emalloc(unsigned int n);
void efree(void *p);

void savemodel (const char *filnam);
void restoremodel (const char *filnam);

void ncmdlin(int argc, char **argv);
void step(double time);
void run(void);
// void exit(int n);
void initsim(void);

void __runonplot(void);
void setonplot(void (*plpnt)());
void __runonexit(void);
void set_run_on_exit(void (*runonexpnt)());
void set_run_on_step(void (*runonstpnt)());

double drand(void);
void setrand(int val);
unsigned long irand(void);
double rrand(int ngen);
void initrand (int ngen, int nrseed);
double gasdev(void);
double rgasdev(int ngen);
int  binomdev(double p, int n);
int  poisdev(double xm);
double gamdev(double a);

#include "digfilt.h"

#include "lm_funcs.h"

char *xsystem(const char *str);
double system_speed(void);
char *print_version(double version);
double elap_time(void);
void simwait(double nsecs);
int streq (const char *str1, const char *str2);
double set_int_val (double val);

