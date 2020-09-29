# include <stdlib.h>
# include <string.h>
# include "graph.h"
# include "stdplt.h"

void _adjtextf(int log, char adj, char *fmt, int arg1, int arg2, int arg3, int arg4)
        
           		/* adjustment code */
            
/*
 * do a textf vertically centered on the current y-position, and
 * either left, center, or right adjusted on the current x-position.
 */
  {/*	FILE sbuf; */				/* for sys V */
 /* 	struct _iobuf sbuf; */			/* for vers 7 */
	char str[STRLEN];
	float xold, yold, len;

/*	sbuf._flag = _IOWRT + _IOLBF; */	/* for sys V */
/*	sbuf._flag = _IOWRT+_IOSTRG; */		/* for vers. 7 */
/*	sbuf._ptr = (unsigned char *)str;
	sbuf._cnt = STRLEN;
	sbuf._base = (unsigned char *)(-1);
	sbuf._file = -1;	/* this doesn't seem to work in sys V Venix */

	sprintf(str,fmt, arg1,arg2,arg3,arg4); 	/* use sprintf to calc length*/
/*	_doprnt(fmt, &args, &sbuf); */		/*  instead of doprnt */
	len = charw*(strlen(str));
/*	len = charw*(STRLEN - sbuf._cnt);
	putc('\0', &sbuf);  */
	fflush(stdplt);
	xold = _Xpos;
	yold = _Ypos;

	switch(adj)
	  {
	  case 'l':
		rmove(0., -.5*charh); break;
	  case 'c':
		rmove(-.5*len, -.5*charh); break;
	  case 'r':
		rmove(-len, -.5*charh); break;
	case 'v':
		cmove(0.,.5*len/charw);
		_vertextf(str);
		move(xold, yold);
		return;
	  default:
		_err("adjtextf: bad adjustment %c\n", adj);
	  }
	if(log)
	  {
		cmove(0.,-.3);
		textf("10");
		cmove(2.5,.3);
	  }
	textf("%s", str);
	move(xold, yold);
  }

void _vertextf(char *p)
{
	while(*p) {
		textf("%c",*p++);
		cmove(0.,-1.);
	}
  }


double _nicetic(double *min, double *max, double ticinc, int *mtic, int logflg)
{
	double newmin, newmax;
	int n;
	double val, temp;
	double log(double), pow(double, double), floor(double), ceil(double),xdist(double dx, char *fname),ydist(double dy, char *fname);
	double  fabs(double);

/* Determines tic spacing */
	n = floor((log(ticinc))/(LOG10));
	temp = ticinc/(pow(10.,(double)n));

	if(temp <=1.) val = 1.;
	else if (temp <=2.) val = 2.;
	else if (temp <=5.) val = 5.;
	else if (temp <10.) val = 10.;

	if(logflg && (n < 0))
	  {
		if(temp > 1.)
			val = 10;
		*mtic = 9;
	  }
	else
		*mtic = val;


/* Adjusts min and max to match tic inc */

	ticinc = val * pow(10.,(double)(n));
	newmin = floor((*min / ticinc) + EPSILON) * ticinc;
	newmax = ceil((*max / ticinc) -  EPSILON) * ticinc;
	*min = newmin; *max = newmax;
	return(ticinc);
}

void _square(float xl, float xh, float yl, float yh)
{
	origin(xl,yl);
	if(xh <= xl)
	  {
		printf("Warning: scale for graph");
		printf(" radically small in X direction\n");
		xh = xl + 1;
	  }
	if(yh <= yl)
	  {
		printf("Warning: scale for graph");
		printf(" radically small in Y direction\n");
		yh = yl + 1;
	  }
	scale(xh-xl,yh-yl);
}

void _error(char *str)
{
	printf("%s\n",str);
}

void _feror(char *str)
{
	printf("%s\n",str);
	exit(0);
}

double LOG(double val)
{
	double log(double);

	return(log(val) / LOG10);
}
