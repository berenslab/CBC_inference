# include "graph.h"
# include "stdplt.h"
double lincon[9]  = {			/* logarithmic axis tics */
		.301029957,
		.1760912591,
		.1249387366,
		.096910013,
		.079181246,
		.0669467896,
		.0579919469,
		.0511525224,
		.0457574905
		  };

double chsize[8] = {
		.0,
		.010,
		.0156,
		.02,
		.025,
		.030,
		.038,
		.047
		};

graph(char *strp, int args)
{
	struct axtype *p;
	static char gstring[1000];
	double LOG(double val);

	argp = (union ptr *)&args;
	zinit("%$ %G",strp,gstring);
	yyparse();

	pair();

	Do_min(_xdata,X);
	Do_min(_ydata,Y);

	_square((11. + ngraphs)*charw,1. -4.*charw,4.*charh,1. - charh);
	frame("Zok");

	Do_nice(&_ax[X],X,0);
	Do_scale(X);
	Do_nice(&_ax[Y],Y,0);
	Do_scale(Y);

	base();

	p = _ax;
	_axis(p++);
	_axis(p);

	Do_disp();
	_grid(p-1,p);

	frame("..");
	rmframe("Zok");
  }
