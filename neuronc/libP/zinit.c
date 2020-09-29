#include <string.h>
#include "graph.h"
#include "stdplt.h"

void zinit(char *source1, char *source2, char *dest)
{
	int i;
	struct axtype *p;

	ychar = 0;
	strcpy(dest,source1);
	strbeg = source2;
	strcat(dest,source2);
	scanp = dest;
	min[X] = min[Y] = INFINITY;
	max[Y] = max[X] = -INFINITY;
	for(i = 0;i<10;i++)
		_xdata[i] = _ydata[i] = 0;
	nary[X] =  nary[Y] = ngraphs = 0;
	p = _ax;
	p->flags = p->x = p->y = p->len = p->ntics = p-> mtics = p->ticinc =
	p->type = 0; p->lab = (char *)0;
	p++;
	p->flags = p->x = p->y = p->len = p->ntics = p-> mtics = p->ticinc =
	p->type = 0; p->lab = (char *)0;
	first = 0;
  }
