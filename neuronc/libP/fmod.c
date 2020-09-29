


double fmod(double da, double db)
{
	register int n;

	n = da/db;
	return(da - n*db);
}
