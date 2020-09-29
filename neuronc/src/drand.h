/* drand.h */


double drand(void);
double rrand(int ngen);
void initstate(unsigned long int val, char *rstate, int size);
void drand_setstate(char *rstate);
void restorstate(void);
void setrand(int val);

