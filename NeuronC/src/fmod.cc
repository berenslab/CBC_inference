double
fmod(double da, double db)
{
	register n;

	n = da/db;
	return(da - n*db);
}
