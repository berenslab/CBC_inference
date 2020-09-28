/*
 * floor and ceil-- greatest integer <= arg
 * (resp least >=)
 */

double  modf(double, double *);

double
floor(double d)
{
	double fract;

	if (d<0.0) {
		d = -d;
		fract = modf(d, &d);
		if (fract != 0.0)
			d += 1;
		d = -d;
	} else
		modf(d, &d);
	return(d);
}

double
ceil(double d)
{
	return(-floor(-d));
}


