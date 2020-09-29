/* file stim.h in program stim.c */

#define NODENUM 200

typedef struct  RECNOD {
	nodeint recnm1;
	nodeint recnm2;
	nodeint recnm3;
	nodeint recnm4;
	double xpos;
	double ypos;
	double actual;
	double curval;
	double wavel;
	double stim;
	double backgnd;
	struct RECNOD *next;
	struct RECNOD *last;
	} recnod;

/* #define CONVSIZE 256  	/* for 16 bit machines */
#define CONVSIZE 512		/* for large machines */

#define BACKARR 0
#define STIMARR 1

