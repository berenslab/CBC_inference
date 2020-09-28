# include "graph.h"
# include "stdplt.h"
extern double lincon[];
extern double chsize[];

extern void fat(int);
extern void _adjtextf(int,char,char*,int);

void _axis(struct axtype *p)
{
	double disp;
	double x,y;
	float xtic,ytic;
	double max;
	double total;
	int i;
	double fabs(double);
	double xdist(double dx, char *fname);
	double ydist(double dy, char *fname);
	double pow(double, double);
	double floor(double);
	int ticflg;
	int labflg;
	static int csiz = 6;

/*	ticflg = p->flags&TICS; */
	ticflg = TIC;		/* length of tic marks are only one style */
	labflg = p->flags&LABS;
	fat(0);
	if(p->type == X)
	  {
		if (p->dflags & XTXT)
		   csiz = p->dflags & XTXT;
		cwidth(chsize[csiz],"/"); 
		xtic = 0.;
		ytic = ydist(.01,"/");
		if(ticflg == TIC)
			ytic *= -1;
		if(ticflg == TIC2)
			ytic *= -2;
		vector(p->x,p->y,p->x + p->len,p->y);
		fat(0);
		move(p->x + (p->len / 2),p->y);
		if(p->flags&NAME)
			cmove(0.,3.);
		else
			cmove(0.,-3.);
		if(p->lab)
			_adjtextf(0,'c',"%s",(int)p->lab);
		total = p->x;
		max = p->x + p->len;
	  }
	else
	  {
		if (p->dflags & YTXT)
		   csiz = (p->dflags & YTXT) >> 3;
		cwidth(chsize[csiz],"/"); 
		xtic = xdist(.01,"/");
		if(ticflg == TIC)
			xtic *= -1;
		if(ticflg == TIC2)
			xtic *= -2;
		ytic = 0.;
		vector(p->x,p->y,p->x,p->y + p->len);
		fat(0);
		move(p->x,p->y + (p->len / 2));
		if(p->flags&NAME)
			cmove(10.,0.);
		else
			cmove(-10.,0.);
		if(p->lab)
			_adjtextf(0,'v',"%s",(int)p->lab);
		total = p->y;
		max = p->y + p->len;
	  }
	fat(0);
	disp = p->ticinc;
	x = p->x;
	y = p->y;
	for(i = 0;total < max + disp/100;i++)
	  {
		if(p->mtics == 9) /* Log axis */
			disp = p->ticinc * lincon[(i%p->mtics)];
		move(x,y);
		total += disp;
		if(i%p->mtics)
		  {
			if(!(ticflg == TIC3))
			  {
				if(ticflg == TIC2)
					rmove(-xtic/4.,-ytic/4.);
				rdraw(xtic/2.,ytic/2.);
			  }
			if(p->type == X)
			  {
				x += disp;
			  }
			else
			  {
				y += disp;
			  }
		  }
		else
		  {
			if(!(ticflg == TIC3))
			  {
				if(ticflg == TIC2)
					rmove(-xtic/2.,-ytic/2.);
				rdraw(xtic,ytic);
			  }
			if(p->type == X)
			  {
				move(x,y);
				/* Make zero 0, not .00000000014... */
				if(fabs(x) < disp/2)
					x = 0.;
				x += disp;
				if(labflg == LAB2)
					continue;
				if(labflg == LAB1)
					cmove(0.,1.5);
				else
					cmove(0.,-1.5);
				if(p->flags & LOGAXIS)
					_adjtextf(1,'c',"%d",(int)floor(.5 + x - disp));
				else
					_adjtextf(0,'c',"%g",(int)(x - disp));
				/* x already incremented by disp */
			  }
			else
			  {
				move(x,y);
				/* Make zero 0, not .00000000014... */
				if(fabs(y) < disp/2)
					y = 0.;
				y += disp;
				if(labflg == LAB2)
					continue;
				if(labflg == LAB1)
				  {
					cmove(1.4,.2);
					if(p->flags & LOGAXIS)
						_adjtextf(1,'l',"%d",(int)floor(.5 + y - disp));
					else
						_adjtextf(0,'l',"%g",(int)(y - disp));
				  }
				else
				  {
					cmove(-1.,.2);
					rmove(xtic,0.);
					if(p->flags & LOGAXIS) {
					        cmove(-3.,0.);
						_adjtextf(1,'r',"%d",(int)floor(.5 + y - disp));
					        cmove(3.,0.);
					}
					else
						_adjtextf(0,'r',"%g",(int)(y - disp));
				  }
				/* y already incremented by disp */
			  }
		  }
	  }
  }
