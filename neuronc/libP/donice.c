# include "graph.h"
# include "stdplt.h"

void Do_nice(struct axtype *p, int xory, int Sscale)
{
	double ticinc,ticdist;
	double _nicetic(double *min, double *max, double ticinc, int *mtic, int logflg);
	double ydist(double dy, char *fname);
	double xdist(double dx, char *fname);
	float length;
	int logflg,ticflg;
	double LOG(double val);

	ticflg = p->flags & TICS;			/* tic distance style */
	ticdist = .009;

	switch (ticflg) {

	case TIC:			/* standard tic dist */ 
	case TIC1: ticdist = .009;	/* standard tic dist, no subtics */
		  break;

	case TIC2: 			/* small tic dist */
	case TIC3: ticdist = .00675;	/* small tic dist, no subtics */
		  break;
	}

	logflg = 0;
	if(p->flags & LOGAXIS)
	  {
		min[xory] = LOG(min[xory]);
		max[xory] = LOG(max[xory]);
		logflg++;
	  }
	if(xory == X)
	  {
		do
		  {
			if(Sscale)
				length = 1.;
			else
				length = (max[X] - min[X]);
			ticinc = 8. * xdist(ticdist,"/") * length;
			ticinc = _nicetic(&min[X],&max[X],ticinc,&p->mtics,logflg);
		  }while(ticinc <= 8. * xdist(ticdist,"/") * length);
	  }
	if(xory == Y)
	  {
		do
		  {
			if(Sscale)
				length = 1.;
			else
				length = (max[Y] - min[Y]);

		     ticinc = 8. * ydist(ticdist,"/") * length;
		     ticinc = _nicetic(&min[Y],&max[Y],ticinc,&p->mtics,logflg);
		  }while(ticinc <= 8. * ydist(ticdist,"/") * length);
	  }

	p->len = max[xory] - min[xory];
	if(p->mtics == 9)
	  {
		p->ticinc = ticinc;
	  }
	else {
		if ((ticflg == TIC1) || (ticflg == TIC3)) p->mtics = 1;
		p->ticinc = ticinc / p->mtics;
	}
	p->type = xory;
  }

void Do_scale(int xory)
{
	float sx,sy,ox,oy;

	sx = sy = 1.;
	ox = oy = 0.;
	if(xory == X)
	  {
		sx = 1./(max[X] - min[X]);
		ox = -min[X];
	  }
	else
	  {
		sy = 1./(max[Y] - min[Y]);
		oy = -min[Y];
	  }
	scale(sx,sy);
	origin(ox,oy);
  }
