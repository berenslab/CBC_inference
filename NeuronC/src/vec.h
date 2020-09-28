/* vec.h */

struct vec {
	double x;
	double y;
	double z;
	};

vec    vcross(vec a, vec b);
double vlen (vec a);
double vdist (vec a, vec b);
double vdot(vec a, vec b);
vec    vadd(vec a, vec b);
vec    vsub(vec a, vec b);
double norm_line (vec a, vec b, vec p, vec q, vec *lp1, vec *lp2);
void vprint(vec a);

