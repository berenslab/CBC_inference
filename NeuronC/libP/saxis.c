# include "graph.h"
# include "stdplt.h"

void Do_min();
void Do_nice();
void Do_scale();
void where();
void _axis();
void _clip();

void _saxis(char *strp, int doscale)
{
	struct axtype *p;
	char astring[100];
	float penx,peny;
	extern int xyflg;

	zinit("%$ %A",strp,astring);
	yyparse();
	if(xyflg == X)
		Do_min(_xdata,X);
	else
		Do_min(_ydata,Y);
	Do_nice(&_ax[xyflg],xyflg,doscale);
	if(!doscale)
		Do_scale(xyflg);
	p = &_ax[xyflg];
	if(xyflg == X)
	  {
		p->x = min[X];
		where(&penx,&peny,".");
		p->y = peny;
	  }
	else
	  {
		p->y = min[Y];
		where(&penx,&peny,".");
		p->x = penx;
	  }
	_axis(p);
	move(penx,peny);
  }

void axis(char *strp, int args)
{
	argp = (union ptr *)&args;
	_saxis(strp,1);
}
void saxis(char *strp, int args)
{
	argp = (union ptr *)&args;
	_saxis(strp,0);
}

