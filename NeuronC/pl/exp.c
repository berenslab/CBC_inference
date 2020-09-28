#include <math.h>
#include "stdplt.h"

float x[600];
float y[600];
main()
{
int	i;
	frame("Frame");
	for(i=0;i<600;i++){
		x[i]=i;
		x[i]/= 18.;
		y[i]= cos(x[i])*exp(-0.1*x[i]);
	}
	graph("%X %*l Time %*n 600 %Y %*l Frequency %*G 3",x,y);
}
