# include "graph.h"
# include "stdplt.h"

void _display(struct _data *xp, struct _data *yp, int dashflg)
{
	int i;
	float xmin,xmax,ymin,ymax;
	float dummy;
	double xdata,ydata,xdisp,ydisp;
	double temp;
	double LOG(double val);
	int flag;
	union ptr *x,*y;
	int n,xflag,yflag,lard;
	char symb;
	static int penflg=0;
	extern double chsize[];
	static int csiz=3;

	x = (union ptr *)&(xp->dptr);
	y = (union ptr *)&(yp->dptr);
	xflag = (xp->flgs&DATATYPE) | (_ax[X].flags&LOGAXIS);
	yflag = (yp->flgs&DATATYPE) | (_ax[Y].flags&LOGAXIS);
	n = 0;
	if(xp->npts > 0)
		n = xp->npts;
	else
		n = yp->npts;
	if(n <= 0)
		return;
	if(xp->symble)
		symb = xp->symble;
	else
		symb = yp->symble;
	if(xp->blubber > yp->blubber)
		lard = xp->blubber;
	else
		lard = yp->blubber;
	if     (xp->flgs & DTXT)
		csiz = xp->flgs & DTXT; 
	else if(yp->flgs & DTXT)
		csiz = yp->flgs & DTXT; 
	if     (xp->flgs & PEN)
		penflg = (int)xp->flgs;
	else if(yp->flgs & PEN)
		penflg = (int)yp->flgs;
	cwidth (chsize[csiz],"/");
	penflg &= PEN;
        switch (penflg) { 
	case PEN0: break;			/* dashes, no pen changes */
	case PEN1: 				/* pen changes and dashes */
	case PEN2: cpen(dashflg+1);  break;	/* pen changes without dashes */
	}

	if (penflg == PEN2) dashflg = 0;
else {
	if(xp->ddash)
		dashflg = (long int)xp->ddash;
	else if(yp->ddash)
		dashflg = (long int)yp->ddash;
	dash(dashflg);
}
	if((min[X] < max[X]) || (min[Y] < max[Y]))
	  {
		if(min[X] < max[X])
		  {
			xmin = min[X];
			xmax = max[X];
		  }
		else
		  {
			xmin = -INFINITY;
			xmax = INFINITY;
		  }
		if(min[Y] < max[Y])
		  {
			ymin = min[Y];
			ymax = max[Y];
		  }
		else
		  {
			ymin = -INFINITY;
			ymax = INFINITY;
		  }
		window(xmin,xmax,ymin,ymax);
	  }
	if(xflag == NODATA)
	  {
		xdisp = (max[X] - min[X])/n;
		xdata = min[X] - xdisp;
	  }
	if(yflag == NODATA)
	  {
		ydisp = (max[Y] - min[Y])/n;
		ydata = min[Y] - ydisp;
	  }
	fat(lard);
	for(i=0; i<n; i++)
	  {
		switch(xflag&DATATYPE)
		  {
			case NODATA:
				xdata += xdisp;
				break;
			case DOUBLE:
		/*		xdata = *x->dp++; */
				xdata = *x->dp;
				x->dp++;
				break;
			case FLOAT:
		/*		xdata = *x->fp++; */
				xdata = *x->fp;
				x->fp++;
				break;
			case INTEGER:
				xdata = *x->ip++;
				break;
			case LONG:
				xdata = *x->lp++;
				break;
			case SHORT:
				xdata = *x->hp++;
				break;
			case UNSIGNED:
				xdata = *x->up++;
				break;
		  }
		switch(yflag&DATATYPE)
		  {
			case NODATA:
				ydata += ydisp;
				break;
			case DOUBLE:
		/*		ydata = *y->dp++; */
				ydata = *y->dp;
				y->dp++;
				break;
			case FLOAT:
		/*		ydata = *y->fp++; */
				ydata = *y->fp;
				y->fp++;
				break;
			case INTEGER:
				ydata = *y->ip++;
				break;
			case LONG:
				ydata = *y->lp++;
				break;
			case SHORT:
				ydata = *y->hp++;
				break;
			case UNSIGNED:
				ydata = *y->up++;
				break;
		  }
		if(xflag&LOGAXIS)
			xdata = LOG(xdata);
		if(yflag&LOGAXIS)
			ydata = LOG(ydata);
		plot(xdata,ydata,i);
		if(symb)
		  {
			cmove(-.5,-.25);
			fat(0);
			textf("%c",symb);
			purge(); 
			fat(lard);
			cmove(.5,.25);
		  }
	  }
	fat(0);
	window(-INFINITY,INFINITY,-INFINITY,INFINITY);
	dash(0);
  }
