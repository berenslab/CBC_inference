/* mdef.h */

#define abs(x)		((x) < 0 ? -(x) : (x))
#define max(x, y)	(((x) < (y)) ? (y) : (x))
#define min(x, y)	(((x) < (y)) ? (x) : (y))
#define limit(v,M,m)	((M>m)?max(min((v),(M)),(m)):max(min((v),(m)),(M)))

