/* code.c */


extern "C" {

//#ifdef CPML
//#include <cpml.h>
//#endif
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "stdplt.h"
}

#include "nc.h"
#include "y.tab.h"
#include "ndef.h"
#include "control.h"

#include "lmmin.h"
#include "lm_eval.h"

#ifdef XSTIMM		/* define XSTIMM when stim should ignore neur elems */
  #define XSTIM		/* define XSTIM for stim */
#endif 

#define isnum(c)  (isdigit(c)||(c)=='+'||(c)=='-'||(c)=='.'||(c)=='e'||(c)=='E')

/* #define pushm(d) stack[stackp++] = (d)	/*  */
#define popm()   stack[--stackp]	/* function pop still needed */

/* #define CODEBUG			/*  */
/* #define popm() pop()			/*  */
#define pushm(d) push(d)		/*  */

#define MAXSIZE 134217728

#ifndef MAXSIZE
# define MAXSIZE 8191				/* max size of double array */
# define NSTACK 400
# define NPROG  2000
#else						/* larger machine */
# if MAXSIZE < 65536
#  define NSTACK 1000
#  define NPROG  4096
# else
#  define NSTACK 1048576
#  define NPROG  1048576
# endif
#endif


extern char *progname;
extern char *runfile;
extern int interp;

static datum ostack[NSTACK];    /* the stack */
static int stacksize = NSTACK;
datum *stack=ostack;            /* original stack stack */
int stackp=0;                   /* next free spot on stack */

Inst	prog[NPROG];	/* the machine */
Inst	*progp=prog;	/* next free spot for code generation */
Inst	*pc;		/* program counter during execution */
Fpa	*fpa;		/* pointer to function with arguments */
Inst	funcs[NPROG];	/* procedures and functions */
Inst	*funcbase = funcs; /* start of curent subprogram */
int	returning=0;	/* 1 if return stmt seen */
int	continuing=0;	/* 1 if continue stmt seen */
int	breaking=0;	/* 1 if break stmt seen */
int	stopping=0;	/* 1 if return, continue, or break stmt seen */

Inst *savprogp;		/* for savecode(), restorecode(), for use in fread */
int savstackp;
int savreturning;
int savcontinuing;
int savbreaking;
int savstopping;

extern double fread_val;/* interpreted value of expression for fread */



extern int	iflag;		/* stop if ^C from keyboard [or exit() ] */
extern int freadindex;		/* index for fread string in fread_interp()*/

typedef struct Frame {	/* proc/func call stack frame */
	Symbol	*sp;	/* symbol table entry */
	Inst	*retpc;	/* where to resume after return */
	int	argn;	/* n-th argument on stack */
	int	nargs;	/* number of arguments */
	int	nlocal;	/* number of local variables */
} Frame;

#define NFRAME 2000
Frame	pframe[NFRAME];
Frame	*fp = pframe;		/* frame pointer */
Frame *savfp;

static char evnon[]     = {"attempt to evaluate non-variable"};
static char varuninit[] = {"variable not initialized"};
static char uninit[]    = {"uninit "};
static char diftyp[]    = {"variables are different type"};
static char divzero[]   = {"divide by zero"};
static char nvalidtyp[] = {"invalid type for operation"};

extern char *rootframe;		/* root graphics frame */

#ifdef __cplusplus
extern "C" {
#endif

// double pow(double a, double b);
double atof(const char *);

/*
char *strcat (char *s1, const char *s2);
char *strtok(char *s, const char *delim);
char *strstr(const char *haystack, const char *needle);
char *strncpy(char *dest, const char *src, size_t n);
*/
FILE *fopen (const char *path, const char *mode);
#include "gr.h"

#ifdef __cplusplus
}
#endif

#include "ncio.h"
#include "gprim.h"

ncarray *garray(Symbol *sp);
double *garrayval(Symbol *sp);
char *emalloc(unsigned int n);
double *darr2(Symbol *sp, int narg);
void execerror (const char *s, const char *t);
void constpush(void);
void execute (Inst *pc);
void gtext(const char *s);
void glabel(const char *s, int pen, double x, double y);
void fft (double *real, double *imag, int size, int type);
void mrqmin(double x[], double y[], double sig[], int ndata, 
 double a[], double ia[], int ma, double **covar, double **alpha, double *chisq, 
 void (*funcs) (double, double [], double *, double [], int), double *alamda);

int checknum(datum d);
int checkstr(datum d);
void erarr (ncarray *arrp);
void efree(void *ptr);
void plotpen (int val, int i);
Symbol *setarrayconst(int nvals);
void dup(void);
void dupcopy(void);
void rmsym(Symbol *s);
double getvarval(const char *name);
double fread_interp(char *str);   /* intepret expressions by fread */
datum ncpow(datum &x, datum &y);
void plotarr(double *var, int maxindex);

static const char *nullstr={""};

/*-----------------------------------------------------------*/

#ifdef CODEBUG			/*  */
int codebug = 0;

/* see setdebug call in xxor(): a way to turn on tracing in interp */

void setdebug (int val)
{
  codebug = val;
}
#endif

/*-----------------------------------------------------------*/

void xdebugf(void)
{
  /* ncfprintf (stdout,"anything you want...\n"); */
}

/*-----------------------------------------------------------*/

void initcode(void) {
	int i;

	progp = prog;
	funcbase = funcs;
	stackp = 0;
	fp = pframe;
	returning = 0;
	continuing = 0;
	breaking = 0;
	stopping = 0;
        for (i=0; i<NFRAME; i++) {
             pframe[i].sp    = (Symbol *)NULL;
             pframe[i].retpc = (Inst *)NULL;
             pframe[i].argn  = 0;
             pframe[i].nargs = 0;
             pframe[i].nlocal = 0;
	}
}

void resetcode(void) {
	progp = prog;
	stackp = 0;
	fp = pframe;
	returning = 0;
	continuing = 0;
	breaking = 0;
	stopping = 0;
}

void savecode(void) {
	savprogp = progp;
	savstackp = stackp;
	savfp = fp;
	savreturning = returning;
	savcontinuing = continuing;
	savbreaking = breaking;
	savstopping = stopping;
}

void restorecode(void) {
	progp = savprogp;
	stackp = savstackp;
	fp = savfp;
	returning = savreturning;
	continuing = savcontinuing;
	breaking = savbreaking;
	stopping = savstopping;
}

/*------------------------------------------*/

void print_datum(const char *str,int d)

/* print a datum for debugging */

{
    datum *dp;

   dp = &stack[d];
   ncfprintf (stderr,"%5s %3d %3d %3d %x %g\n",
		str,dp->type,dp->atype,dp->vtype,dp->arrp,dp->val);
}

/*------------------------------------------*/

void checkstack(int sp)

/* Check stack and expand if too small. */
/*  Must not be called in case where pointer to stack is saved */
/*  and then used, because expanded stack is not in same location. */

{
   datum *nsp,*osp,*new_stack;
   int i,new_stacksize;

   if (sp >= stacksize) {		/* if stack too small */
      new_stacksize = stacksize * 2;    /*  allocate new stack 2x bigger */
      if ((new_stack=(datum*)emalloc(new_stacksize*sizeof(datum))) == NULL) {
          execerror ("stack too deep", (char *)0);
      }
      osp = stack;
      nsp = new_stack;
      for (i=0; i<stacksize; i++) {
         *nsp++ = *osp++;
      }
      if (stacksize > NSTACK) efree (stack);  /* don't free orig array */
      stack = new_stack;
      stacksize = new_stacksize;
   }
}



/*------------------------------------------*/

void push(datum d)
{
#ifdef CODEBUG
  if (codebug) ncfprintf (stderr,"push\n");
#endif
        checkstack(stackp);
	stack[stackp++] = d;
#ifdef CODEBUG
  if (codebug) print_datum("push",stackp-1);
#endif
}


void expop(void)
{
	if (stackp == 0)
		execerror ("stack underflow", (char *)0);
	--stackp;
	return;
}

datum pop(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"pop\n");
#endif
	if (stackp == 0)
		execerror ("stack underflow", (char *)0);
#ifdef CODEBUG
	  if (codebug) print_datum("pop",stackp-1);
	  if (codebug) print_datum("pop",stackp-2);
#endif
	return stack[--stackp];
}

/*------------------------------------------*/

static datum dsave;

void popsav(void)

/* temporarily save a value from the stack */

{
	if (stackp == 0)
		execerror ("stack underflow", (char *)0);
	dsave = stack[--stackp];
}

void pushsav()
/* restore the stack from previously saved value */
{
	checkstack(stackp);
	stack[stackp++] = dsave;
}

/*------------------------------------------*/

void stackmove(void)

/* move the top value on the stack down n levels */

{
     datum dsave;
     int i,n;
     int stp;

   n = (long int) (*pc++);
   if (n<0) n = -n;
   stp = stackp; 
   if ((stp - n) < 0) {
	execerror ("stackmove: stack underflow", (char *)0);
	return;
   }
   dsave = stack[--stp];		/* pop value off stack */
   for (i=0; i<n; i++,stp--) {
     stack[stp] = stack[stp-1];
   }
   stack[stp] = dsave;		/* push value on bottom */
}

/*------------------------------------------*/

int checknum(datum d)
{
  int val;

  switch(d.vtype) {
	case NUMBER:
		val = 1;
		break;
	case LITCHAR:
		val = 0;
		break;
	default:
	case STRING:
		execerror ("wrong type: must have number",0);
		val = 0;
		break;
  }
  return val;
}

/*------------------------------------------*/

int checkstr(datum d)
{
  switch(d.vtype) {
	case NUMBER:
	case LITCHAR:
	default:
		execerror ("wrong type: must have string",0);
		return 0;
		break;
	case STRING:
		return 1;
		break;
  }
}
/*------------------------------------------*/

int checkchar(datum d)

{
    int ch;

  switch(d.vtype) {
	  case STRING:
	    ch = d.str[0];
	    break;
	  case LITCHAR:
	    ch = int(d.val);
	    break;
	  case NUMBER:
	    ch = int(d.val);
	    break;
	}
   return ch;
}

/*------------------------------------------*/

void constpush(void)
{
	datum d={0};

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"constpush\n");
#endif

	d.type = CONST; 
	d.vtype = ((Symbol *)*pc)->type;	/* since literal, vtype=type */
	switch (d.vtype) {
	  case NUMBER:
	  case LITCHAR:
		d.val = ((Symbol *)*pc)->val;
		break;
	  case STRING:
		d.str = ((Symbol *)*pc)->name;
		break;
	}

	pc++;
	pushm (d);
}

/*------------------------------------------*/

void varpush(void)

/* Get variable from code line.  Push 2 datum, d1 is type, d2 is
 * number of dimensions */

{
	datum d1={0},d2={0};
        int type,nvals;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"varpush\n");
#endif
	type = (long int) (*pc++);
	d1.type = type;
	switch (type) {
	  case UNDEF:
	  case CONST:
	  case VAR:
		d1.sym = (Symbol *)(*pc++);
		d1.atype = GLOBAL;
		d1.vtype = d1.sym->vtype;    /* set NUMBER, STRING, etc.*/
		pushm(d1);
		break;
	  case ARG:
		d1.argnum = (long int)(*pc++);
		d1.atype = ARG;
		pushm(d1);
	        break;
	  case LOCALVAR:
		d1.argnum = (long int)(*pc++);
		d1.atype = LOCAL;
		pushm(d1);
	        break;
	  case ARGARRVAL:
		d1.argnum = (long int)(*pc++);
		d1.atype = ARG;
	        d2.argnum = (long int)(*pc++); /* don't convert to double */
	        d2.type = VAR;
	        pushm(d2);		  /* number of dimensions */
		pushm(d1);
	        break;
	  case LOCALARRVAL:
		d1.argnum = (long int)(*pc++);
		d1.atype = LOCAL;
		d2.argnum = (long int)(*pc++); /* don't convert to double */
		d2.type = 0;
		pushm(d2);		  /* number of dimensions */
		d1.vtype = 0;	 	  /* set NUMBER, STRING, etc.*/
		pushm(d1); 
		break;
	  case LOCALARR:
		d1.argnum = (long int)(*pc++);
		d1.atype = LOCAL;
		d2.argnum = (long int)(*pc++); /* don't convert to double */
		d2.type = 0;
		//pushm(d2);		  /* number of dimensions */
		d1.vtype = 0;	 	  /* set NUMBER, STRING, etc.*/
		pushm(d1); 
		break;
	  case ARRAYVAL:
		d1.sym = (Symbol *)(*pc++); /* type of variable */
		d1.atype = GLOBAL;
		d2.argnum = (long int)(*pc++); /* don't convert to double */
		d2.type = VAR;
		pushm(d2);		    /* number of dimensions */
		d1.vtype = d1.sym->vtype;   /* set NUMBER, STRING, etc.*/
		pushm(d1); 
		break;
	  case ARRAY:
		d1.sym = (Symbol *)(*pc++); /* type of variable */
		d1.atype = GLOBAL;
		d2.argnum = (long int)(*pc++); /* don't convert to double */
		d2.type = VAR;
		//pushm(d2);		    /* number of dimensions */
		d1.vtype = d1.sym->vtype;   /* set NUMBER, STRING, etc.*/
		pushm(d1); 
		break;
	  case ARRAYCONST:
		nvals = (long int)(*pc++);
		d1.sym = setarrayconst(nvals);
		// d1.type = ARRAY;
		d1.atype = LOCAL;
		d1.vtype = d1.sym->vtype;    /* set NUMBER, STRING, etc.*/
		pushm(d1);
		break;
	default: 
		ncfprintf (stderr,"type %d\n",d1.sym->type);
		execerror ("Assignment to non-variable", d1.sym->name);
		break;
	  }
}

/*------------------------------------------*/

void breakcode(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"breakcode\n");
#endif
	breaking = stopping = 1;
}

/*------------------------------------------*/

void contcode(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"contcode\n");
#endif
	continuing = stopping  = 1;
}

/*------------------------------------------*/

int forval(datum d)
{
   int val;
   const char *null= {""};

   switch (d.vtype) {
     case NUMBER:
     case LITCHAR:
	val = (int)(d.val);
	break;
     case STRING:
	if (d.str==NULL) d.str = null;
	val = strlen(d.str);
	break;
   }
   return val;
}

/*------------------------------------------*/

void break_fixup(int before, int after)

/* reset stack and arg counters when continuing in "for" or "while" loop */
{
    int i,nargs;
    datum d={0};

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"break_fixup\n");
#endif

	nargs = after - before;
	//ncfprintf (stderr,"nargs %d\n",nargs);
	fp->nargs -= nargs;
	fp->nlocal -= nargs;
	for (i=0; i<nargs; i++) {       /* decrement stack for local vars */
	    d = popm();
	    if (d.atype==LOCAL && d.type==LOCALARR) erarr(d.arrp);
	}
	fp->argn = stackp - 1;          /* last local var */
}

/*------------------------------------------*/

void forcode(void)
{
	datum d;
	Inst *savepc = pc;
	int savestackp = stackp;
	int val;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"forcode\n");
#endif
	execute (savepc+4);	/* init */
	d = popm();
	execute(savepc+(long int)*savepc); /* condition */
	d = popm();
	while ((val=forval(d))&&!stopping) {
		execute(savepc+(long int)*(savepc+2)); /* body */
		if (stopping) {
//			stackp = savestackp;	
			if (returning)
				break;
			stopping = 0;
			if (breaking) {
				break_fixup(savestackp,stackp);
				breaking = 0;
				break;
			}
			continuing = 0;
		}
		execute(savepc+(long int)*(savepc+1)); /* end loop expr*/
		d = popm();
		execute(savepc+(long int)*savepc); /* condition */
		d = popm();
	}
	if (!returning)
		pc = savepc+(long int)*(savepc+3);	/* next stmt */
}

/*------------------------------------------*/

void whilecode(void)
{
	datum d;
	Inst *savepc = pc;
	int savestackp = stackp;
	int val;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"whilecode\n");
#endif
	execute (savepc+2);	/* condition */
	d = popm();
	while ((val=forval(d))&&!stopping) {
		execute(savepc+(long int)*savepc);	/* body */
		if (stopping) {
			if (returning)
				break;
			stopping = 0;
			if (breaking) {
				break_fixup(savestackp,stackp);
				breaking = 0;
				break;
			}
			continuing = 0;
		}
		execute(savepc+2);	/* condition */
		d = popm();
	}
	if (!returning)
		pc = savepc+(long int)*(savepc+1);	/* next stmt */
}

/*------------------------------------------*/

void ifcode(void)
{
	datum d;
	Inst *savepc = pc;	/* then part */
	int savestackp = stackp;
	int val;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"ifcode\n");
#endif
	savestackp = stackp;
	execute(savepc+3);	/* condition */
	d = popm();
	if ((val=forval(d)))
		execute(savepc+(long int)*savepc);
	else if (*(savepc+1)) /* else part? */
		execute(savepc+(long int)*(savepc+1));
	if (!returning)
		pc = savepc+(long int)*(savepc+2); /* next stmt */
}

/*------------------------------------------*/

void copycode (Inst *src, Inst *dest, int n)

{
  int i;

  for (i=0; i<n; i++) {
    *dest++ = *src++;
  }
}

/*------------------------------------------*/

void define(Symbol *sp, int narg, Inst *p, int argcount)

/* put func/proc in symbol table */

{
     int n;
     char *s;
     Symbol *ssp;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"define\n");
#endif
   if (argcount>0) {
    if (argcount==1) s=(char *)"";
    else             s=(char *)"s";
    ncfprintf (stderr,
	"# nc: define proc: %d local variable%s.\n", argcount,s);
    ncfprintf (stderr,
	"# proc '%s': Can't define a proc or func when\n",sp->name);
    ncfprintf (stderr,"#  local variables are already defined.\n");
   }
   if (narg== -1) {		/* if substituting function/procedure name */
     ssp = (Symbol *)p; 	/* symbol table entry */
     sp->defn = ssp->defn;	/* start of code */
     sp->vtype = ssp->vtype;	/* number of arguments declared */
   }
   else if (narg== -2) {	/* if substituting bltin func name */
     ssp = (Symbol *)p; 	/* symbol table entry */
     sp->defn = (Inst*)ssp->ptr; /* start of code */
     sp->vtype = 0;		/* ignore number of arguments  */
   }
   else {
     sp->defn = funcbase;	/* start of code */
     sp->vtype = narg;	/* number of arguments declared */
     n = progp - p;
     copycode(p,funcbase,n);	/* copy func to separate func space */
     funcbase += n;		/* next code starts here */
     progp = p;
   }
}

/*------------------------------------------*/

void local(void)	/* add # of local variables to arg num in stack frame */

{
    int i,d,s,nargs,narr,targn;
    datum d1;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"local\n");
#endif
	d1.val = LARGENUM;		/* set to uninitialized value */
	d1.type = 0;			/* reset type */
	d1.atype = LOCAL;		/* reset value type */
	d1.vtype = 0;			/* reset value type */
	nargs = (long int)*pc++;
	narr  = (long int)*pc++;		/* number of local arrays */
	fp->nargs += nargs;
	fp->nlocal += nargs;
	/* ncfprintf (stderr,"local: fpnargs %d nlocal %d nargs %d narr %d\n",
				fp->nargs,fp->nlocal,nargs,narr); */
	for (i=0; i<(nargs-narr); i++) 	/* increment stack for local vars */
	    pushm(d1);
	fp->argn = stackp - 1;		/* last local var */
	for (i=0; i<narr; i++) {	/* recover local array pointers */
	  s = i-(nargs-1);
          targn = stack[fp->argn+s].type;
	  d = targn - fp->nargs;
	  if (s!=d) {			/* don't move if no extra args */
             stack[fp->argn+d].arrp  = stack[fp->argn+s].arrp;  /* move to correct place */
             stack[fp->argn+s].arrp  = (ncarray *)NULL; 
             stack[fp->argn+d].atype = LOCAL; 
             stack[fp->argn+s].type  = 0; 
             stack[fp->argn+d].vtype = stack[fp->argn+s].vtype; 
             stack[fp->argn+s].vtype = 0;
	  }
          stack[fp->argn+d].type  = LOCALARR; 
	}	
}

/*------------------------------------------*/

void locend(void)

{
    int i,nargs;
    datum d={0};

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"localend\n");
#endif

	nargs = (long int)*pc++;
	fp->nargs -= nargs;
	fp->nlocal -= nargs;
	for (i=0; i<nargs; i++) { 	/* decrement stack for local vars */
	    d = popm();
	/* ncfprintf (stderr,"locend: erase %d %d \n",
				d.type,d.atype);/* */
	    if (d.atype==LOCAL && d.type==LOCALARR) {
	 /* ncfprintf (stderr,"locend: erase %d %d ndim %d dim1 %d\n",
			d.type,d.atype,d.arrp->ndim,d.arrp->dim[0]);/* */
		erarr(d.arrp);
	    }
	}
	fp->argn = stackp - 1;		/* last local var */
}

/*------------------------------------------*/

void call(void)		/* call a function */
{
	Symbol *sp;
					/* for function */
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"call\n");
#endif
 	sp = (Symbol *)pc[0]; /* symbol table entry */
	if (fp++ >= &pframe[NFRAME-1])
		execerror(sp->name, ": call nested too deeply");
	fp->sp = sp;
	fp->nargs = (long int)pc[1];	/* # of args passed from call */
	fp->nlocal = 0;			/* reset # of local variables */
	fp->retpc = pc + 2;
	fp->argn = stackp - 1;	/* last argument (including local variables) */
	if (sp->vtype>=0 && sp->vtype != fp->nargs) 
		execerror(": wrong number of arguments in proc",sp->name);
	if (sp->type == UNDEF) {
		execerror(sp->name, ": undefined function");
		returning = stopping = 1;
	}
	else {
		execute(sp->defn);
		returning = stopping = 0;
	}
}

/*------------------------------------------*/

void callproc(Symbol *procp, int npar, double par1, double par2)
{
    datum d1,d2;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"callproc\n");
#endif
   d1.val = par1; 		/* simulate a call from interpreter */
   d1.type = VAR; 
   d1.vtype = NUMBER; 
   d2.val = par2; 
   d2.type = VAR; 
   d2.vtype = NUMBER; 
   if (npar>0) pushm (d1);
   if (npar>1) pushm (d2);
   if (fp++ >= &pframe[NFRAME-1])
		execerror(procp->name, ": call nested too deeply");
   fp->sp = procp;
   fp->nargs = npar; 
   fp->nlocal = 0;			/* reset # of local variables */
   fp->retpc = pc;
   fp->argn = stackp - 1;	/* last argument (including local variables) */
   if (procp->vtype != fp->nargs) 
		execerror(": wrong number of arguments in proc", procp->name);
   execute(procp->defn);
   returning = stopping = 0;
}

/*------------------------------------------*/

double callfunc(Symbol *funcp, int npar, double par1, double par2, 
					 double par3, double par4)
{
    datum d1,d2,d3,d4;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"callfunc\n");
#endif
   d1.val = par1; 		/* simulate a call from interpreter */
   d1.type = VAR; 
   d1.vtype = NUMBER; 
   d2.val = par2; 
   d2.type = VAR; 
   d2.vtype = NUMBER; 
   d3.val = par3; 
   d3.type = VAR; 
   d3.vtype = NUMBER; 
   d4.val = par4; 
   d4.type = VAR; 
   d4.vtype = NUMBER; 
   if (npar>0) pushm (d1);
   if (npar>1) pushm (d2);
   if (npar>2) pushm (d3);
   if (npar>3) pushm (d4);
   if (fp++ >= &pframe[NFRAME-1])
		execerror(funcp->name, ": call nested too deeply");
   fp->sp = funcp;
   fp->nargs = npar; 
   fp->nlocal = 0;			/* reset # of local variables */
   fp->retpc = pc;
   fp->argn = stackp - 1;	/* last argument (including local variables) */
   if (funcp->vtype != fp->nargs) 
		execerror(": wrong number of arguments in func", funcp->name);
   execute(funcp->defn);
   returning = stopping = 0;
   d1 = popm(); 
   return (d1.val);
}

/*------------------------------------------*/

double callfuncp(Symbol *funcp, int npar, double *par)
{
    datum d1,d2;
    int i;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"callfuncp\n");
#endif

   for (i=0; i<npar; i++) {
     d1.val = par[i]; 		/* simulate a call from interpreter */
     d1.type = VAR; 
     d1.vtype = NUMBER; 
     pushm (d1);
   }
   if (fp++ >= &pframe[NFRAME-1])
		execerror(funcp->name, ": call nested too deeply");
   fp->sp = funcp;
   fp->nargs = npar; 
   fp->nlocal = 0;			/* reset # of local variables */
   fp->retpc = pc;
   fp->argn = stackp - 1;	/* last argument (including local variables) */
   //if (funcp->vtype != fp->nargs) 
   //		execerror(": wrong number of arguments in func", funcp->name);
   execute(funcp->defn);
   returning = stopping = 0;
   d2 = popm(); 
   return (d2.val);
}

/*------------------------------------------*/

double callfuncxp(Symbol *funcp, double xval, ncarray *par, int model_num)
	/* call function for do_lmfit() below */
{
    datum d1,d2;
    int i;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"callfuncxp\n");
#endif

   d1.val = xval; 		/* simulate a call from interpreter */
   d1.type = VAR; 
   d1.vtype = NUMBER; 
   pushm (d1);
   d1.arrp = par; 		/* simulate a call from interpreter */
   d1.type = ARRAY; 
   d1.vtype = NUMBER; 
   pushm (d1);
   if (fp++ >= &pframe[NFRAME-1])
		execerror(funcp->name, ": call nested too deeply");
   fp->sp = funcp;
   fp->nargs = 2; 
   fp->nlocal = 0;			/* reset # of local variables */
   fp->retpc = pc;
   fp->argn = stackp - 1;	/* last argument (including local variables) */
   //if (funcp->vtype != fp->nargs) 
   //		execerror(": wrong number of arguments in func", funcp->name);
   execute(funcp->defn);
   returning = stopping = 0;
   d2 = popm(); 
   return (d2.val);
}

/*------------------------------------------*/

double callfuncxp2d(Symbol *funcp, double xval, double yval, ncarray *par, int model_num)
	/* call function for do_lmfit() below */
{
    datum d1,d2;
    int i;

#ifdef CODEBUG
  ncfprintf (stderr,"callfuncxp2d\n");
#endif

   d1.val = xval; 		/* simulate a call from interpreter */
   d1.type = VAR; 
   d1.vtype = NUMBER; 
   pushm (d1);
   d1.val = yval; 		
   d1.type = VAR; 
   d1.vtype = NUMBER; 
   pushm (d1);
   d1.arrp = par; 	
   d1.type = ARRAY; 
   d1.vtype = NUMBER; 
   pushm (d1);
   if (fp++ >= &pframe[NFRAME-1])
		execerror(funcp->name, ": call nested too deeply");
   fp->sp = funcp;
   fp->nargs = 3; 
   fp->nlocal = 0;			/* reset # of local variables */
   fp->retpc = pc;
   fp->argn = stackp - 1;	/* last argument (including local variables) */
   //if (funcp->vtype != fp->nargs) 
   //		execerror(": wrong number of arguments in func", funcp->name);
   execute(funcp->defn);
   returning = stopping = 0;
   d2 = popm(); 
   return (d2.val);
}

/*------------------------------------------*/

void callfunc8 (Symbol *funcp, int npar, double par1, double par2,
			double par3, double par4, double par5,
			double par6, double par7, double par8)

/* simulate a call with 8 arguments and no return */

{
    datum d[8];
    int i;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"callfunc\n");
#endif
   d[0].val = par1; 		/* simulate a call from interpreter */
   d[1].val = par2;
   d[2].val = par3;
   d[3].val = par4;
   d[4].val = par5;
   d[5].val = par6;
   d[6].val = par7;
   d[7].val = par8;
   npar = limit(npar,0,8);
   for (i=0; i<npar; i++) {
     d[i].vtype = NUMBER; 
     d[i].type = VAR;
     pushm(d[i]);
   }
   if (fp++ >= &pframe[NFRAME-1])
		execerror(funcp->name, ": call nested too deeply");
   fp->sp = funcp;
   fp->nargs = npar; 
   fp->nlocal = 0;			/* reset # of local variables */
   fp->retpc = pc;
   fp->argn = stackp - 1;	/* last argument (including local variables) */
   if (funcp->vtype != fp->nargs) 
		execerror(": wrong number of arguments in func",funcp->name);
   execute(funcp->defn);
   returning = stopping = 0;
}

/*------------------------------------------*/

void ret(void)		/* common return from func or proc */
{
	int i,npops;
	datum d={0};
#ifdef CODEBUG
    if (codebug) { ncfprintf (stderr,"ret\n");
		   ncfprintf (stderr,"nargs %x %d %d %x\n",fp,fp->nargs,fp->nlocal,fp->retpc);
    }
#endif
	 
	npops = fp->nargs;
	for (i=0; i<npops; i++) {
	    d = popm();  /* pop arguments */
	    if (d.atype==LOCAL && d.type==LOCALARR) {
		 /* ncfprintf (stderr,"locend: erase %d %d ndim %d dim1 %d\n",
			d.type,d.atype,d.arrp->ndim,d.arrp->dim[0]); /* */
		erarr(d.arrp);
	    }
	}
	pc = (Inst *)fp->retpc;
	--fp;
	returning = stopping = 1;
}

/*------------------------------------------*/

void funcret(void)	/* return from a function */
{
	datum d={0},d2;
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"funcret\n");
#endif
	if (fp->sp->type == PROCEDURE)
		execerror(fp->sp->name, "(proc) returns value");
	d = popm();	/* preserve function return value */
  	if (d.type==ARRAY) {
	  pushm(d);
	  dupcopy();	/* copy array */	
	  d = popm();	/* preserve function return value */
	  d2 = popm();	/* remove orig (possibly local) val */
	}
	ret();
/*  ncfprintf (stderr,"funcret: copy %d %d ndim %d dim1 %d\n",
		d.type,d.atype,d.arrp->ndim,d.arrp->dim[0]); /* */
	pushm(d);
}

/*------------------------------------------*/

void procret(void)	/* return from a procedure */
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"procret\n");
#endif
	if (fp->sp->type == FUNCTION)
		execerror(fp->sp->name,"func returns no value");
	ret();
}

/*------------------------------------------*/

void xexit(void)

/* exit program completely from any stack location. */

{
  ncfprintf (stderr,"Exiting...\n");
  stopping = 1;
  breaking = 1;
  iflag = 1;		/* exit more gently than exit() */
  //exit(1);
}

/*------------------------------------------*/

void getarg(int narg, int ptype, pdatum *p)	

  /* return pointer to argument */
  /* ptype = type of arg from inside procedure */
  /* type  = type of arg passed from outside procedure */

{
	static short int type, atype, vtype;
	static Symbol sp;
	char tbuf[40];
	static ncarray **arrp;

	if (narg > fp->nargs)
	    if (fp->sp) { 
	       ncfprintf (stderr,"nc: getarg: %d %d\n", narg, fp->nargs);
	       execerror(fp->sp->name, ": not passed enough arguments");
            }
	    else
	       execerror((char*)NULL,": not passed enough arguments");

	type =     stack[fp->argn + narg - fp->nargs].type;
	vtype =    stack[fp->argn + narg - fp->nargs].vtype;
	p->type  = &stack[fp->argn + narg - fp->nargs].type;
	p->vtype = &stack[fp->argn + narg - fp->nargs].vtype;
 	arrp =     &stack[fp->argn + narg - fp->nargs].arrp;
	if (ptype==LOCALARRVAL && type==LOCALARR) type = ptype;
	if (ptype==ARGARRVAL && type==ARRAY) type = ptype;
        if (ptype==ARGARRVAL)
	  if (type==VAR) {
	    execerror(fp->sp->name, ": incorrect use of argument"); 
	}
        switch (type) {

	 case LOCALARRVAL: 
 	  sp.arrp = *arrp;
	  sprintf(tbuf,"local array val %d\n",narg);
	  sp.name = tbuf;
	  sp.vtype = vtype;
	  p->val = garrayval (&sp);
          sp.type = VAR;
	  p->type = &sp.type;
	  atype = LOCAL;
	  p->atype = &atype;
	  p->vtype = &vtype;
          break;
	
	 case LOCALARR: 
	  sprintf(tbuf,"local array var %d\n",narg);
	  if (*arrp==NULL || **(double ***)arrp==NULL) {
	    execerror("array space not allocated for",tbuf);
	    return;
	  }
	  p->arrp = arrp;
	  type = ARRAY; 
	  p->type = &type;
	  atype = LOCAL;
	  if (ptype==ARG) atype = ARG;
	  p->atype = &atype;
          break;
	
	 case ARGARRVAL: 
	  sp.arrp = *arrp;
	  sprintf(tbuf,"array arg %d\n",narg);
	  sp.name = tbuf;
	  sp.vtype = vtype;
	  p->val = garrayval (&sp);
          sp.type = VAR;
	  p->type = &sp.type;
          sp.atype = ARG;
	  p->atype = &sp.atype;
	  p->vtype = &vtype;
          break;

	 case ARRAYVAL: 
	  sp.arrp = *arrp;
	  sprintf(tbuf,"array arg %d\n",narg);
	  sp.name = tbuf;
	  p->val = garrayval (&sp);
	  type = VAR; 
	  p->type = &type;
	  atype = GLOBAL;
	  p->atype = &atype;
          break;

	 case ARRAY: 
	  sprintf(tbuf,"array arg %d\n",narg);
	  if (*arrp==NULL || **(double ***)arrp==NULL) {
	    execerror("array space not allocated for",tbuf);
	    return;
	  }
	  p->arrp = arrp;
	  type = ARRAY; 
	  p->type = &type;
	  atype = GLOBAL;
	  if (ptype==ARG) atype = ARG;
	  p->atype = &atype;
	  break;

	default:
          sp.type = VAR;
	  p->atype = &sp.atype;
          sp.atype = ARG;
	  p->type = &sp.type;
	  p->val  = &stack[fp->argn + narg - fp->nargs].val;
           break;
	}
	p->vtype = &stack[fp->argn + narg - fp->nargs].vtype;
}

/*------------------------------------------*/

void bltin(void)

/* call a built-in function, with as many as five arguments */
/*  Arguments may be either NUMBER or STRING. */

{
	int i,narg;
	datum d[5]={0}, d2;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"bltin\n");
#endif
	errno=0;
	narg = (long int)*pc++;			/* get number of arguments */
	for (i=0; i<narg; i++) {
	   d2 = popm();
	   if (d2.val==LARGENUM) 
		execerror ("bltin func,",varuninit);
	   d[narg-i-1] = d2;
	}
	fpa = (Fpa *)pc++;

	switch (narg) {
	case 0: d[0] = (*fpa)();
		break;
	case 1: d[0] = (*fpa)(&d[0]);
		break;
	case 2: d[0] = (*fpa)(&d[0],&d[1]);
		break;
	default:if (narg>5) {
		    narg = 5;
		    execerror ("bltin func, ","more than 5 arguments.");
		}
		d[0] = (*fpa)(&d[0],&d[1],&d[2],&d[3],&d[4]);
		break;
	}
	pushm(d[0]);
}

/*------------------------------------------*/

Symbol *getvar(pdatum *p)

/* Get pointer to variable */
                      
{
	static datum d;
	int narg;
	static ncarray *arrp;
	static short int type;
	static short int atype;
	static Symbol locvar;

	d = popm();			/* get pointer to variable */
	if (d.type==0) {
	   execerror(evnon, 0);
           return ((Symbol *)0); 
	}
	narg = d.argnum;
	switch (d.type) {	/* if var is arg to a function */
	  case LOCALVAR:
	  case LOCALARR:
	  case LOCALARRVAL:
	  case ARGARRVAL:
	  case ARG:
		getarg(narg,d.type,p);	
		switch (*p->vtype) {
		 default:
		 case NUMBER:
	           if (*p->val==LARGENUM) { /* return error if uninitialized */
		      static char ebuf[100];
		      int argn;

		   type = ARG;
		   p->type = &type;
		   if ((argn=narg-(fp->nargs-fp->nlocal)) <=0)
		     sprintf (ebuf,"proc %.40s: #%d argument",fp->sp->name,narg);
		   else if (fp->sp)
			sprintf (ebuf,"proc %.40s: #%d local",fp->sp->name,argn);
		   else
			sprintf (ebuf,"proc: #%d local",argn);
           	   locvar.name=ebuf;	/* fake a proper symbol and name */
           	   locvar.val= *p->val;
		   return &locvar;
		  }
		  break;
		 case STRING: 
		  break;
	      }
	      break;

	  case UNDEF:
	  case CONST:
	  case VAR:
	  case ARRAY:
	     default:		/* else if var is regular variable */

	  p->type  = &d.sym->type; 
	  p->atype = &d.sym->atype; 
	  p->vtype = &d.sym->vtype; 
	  switch (d.type) {
	  case UNDEF:
	  case CONST:
	  case VAR:   p->val = &d.sym->val;
			atype = GLOBAL;
			p->atype = &atype;
		      break;
	  case ARRAYCONST:  
			arrp = garray(d.sym); 
			p->arrp = &arrp;
			type = ARRAYCONST;
			p->type = &type;
			atype = GLOBAL;
			p->atype = &atype;
			p->sym = &d.sym;
			break;
	  case ARRAY:  
			arrp = garray(d.sym); 
			p->arrp = &arrp;
			type = ARRAY;
			p->type = &type;
			atype = GLOBAL;
			p->atype = &atype;
			p->sym = &d.sym;
			break;
	  case ARRAYVAL:  
		        p->val = garrayval(d.sym); 
			type = VAR;
			p->type = &type;
			atype = GLOBAL;
			p->atype = &atype;
		      break;
	  default: execerror(evnon, 0);
		      break;
	  }
	  break;
	}

	if (p->val) {
	 if (*p->type != ARRAY)
	   switch (*p->vtype) {	/* return sym of variable if uninitialized */
	   default:
	   case NUMBER: 
	     if (*p->val==LARGENUM) return (d.sym);
	     break;
	   case STRING: 	
	     if (*p->val==LARGENUM) return (d.sym);
	     if (p->val==NULL) return (d.sym);
	     break;
	   }
	}
	else {
	  p->val = &d.val;
	  return ((Symbol *)d.sym); 
	}
      return ((Symbol *)0); 
}

/*------------------------------------------*/

void evalvar(void)		/* evaluate variable on stack */
{
	datum d1={0};
	pdatum p={0};
	Symbol *s;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"evalvar\n");
#endif
	if ((s=getvar(&p))) {
	   if (strcmp(s->name,"UNINIT")!=0) { /* ignore special uninit const */
             execerror (s->name, varuninit);
	     d1.val=LARGENUM;
	     d1.type  = VAR;	     /* make reasonable-looking variable */
	     d1.vtype = NUMBER;
	     pushm (d1);
	     return;
	   }
	}
//#ifdef CODEBUG
//    if (codebug) ncfprintf (stderr,"var %s\n",s->name);
//#endif
	d1.type  = *p.type;
	d1.atype = *p.atype;	     /* set value type: GLOBAL, ARG, LOCAL, etc.*/
	d1.vtype = *p.vtype;	     /* set value type: NUMBER, STRING, etc. */

 /*      if (d1.type==CONST) d1.vtype = NUMBER;	/* convert CONST to NUMBER */
	if (d1.type == ARRAY || d1.type==ARRAYCONST) {
		d1.arrp = *p.arrp;
	} 
	else
 	 switch (d1.vtype) {
	   case NUMBER:
	   case LITCHAR:
		d1.val  = *p.val;
		break;
	   case STRING:
		if (*p.val==LARGENUM) {	/* preserve value for "uninit" */
		  d1.val = *p.val;
		}
		else d1.str  = *p.str;
	 	break;
	 }
	pushm(d1);
}

/*------------------------------------------*/

void varnum(void)	/* evaluate whether variable on stack is number */
{
	datum d1={0};
	pdatum p={0};
	Symbol *s;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"varnum\n");
#endif
	if (!(s=getvar(&p))) {
	   if (*p.vtype==NUMBER)
		d1.val = 1;
	   else d1.val = 0;  
	}
	else    d1.val = 0;  
	d1.type  = CONST;
	d1.vtype = NUMBER;
	pushm(d1);
}

/*------------------------------------------*/

void varstr(void)	/* evaluate whether variable on stack is string */
{
	datum d1={0};
	pdatum p={0};
	Symbol *s;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"varstr\n");
#endif
	if (!(s=getvar(&p))) {
	   if (*p.vtype==STRING)
		d1.val = 1;
	   else d1.val = 0;  
	}
	else    d1.val = 0;  
	d1.type  = CONST;
	d1.vtype = NUMBER;
	pushm(d1);
}

/*------------------------------------------*/

void varchr(void)	/* evaluate whether variable on stack is char */
{
	datum d1={0};
	pdatum p={0};
	Symbol *s;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"varchar\n");
#endif
	if (!(s=getvar(&p))) {
	   if (*p.vtype==LITCHAR)
		d1.val = 1;
	   else d1.val = 0;  
	}
	else    d1.val = 0;  
	d1.type  = CONST;
	d1.vtype = NUMBER;
	pushm(d1);
}

/*------------------------------------------*/

void notinitx(void)	/* evaluate whether variable on stack is initialized */
{
	datum d1={0};
	pdatum p={0};
	Symbol *s;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"notinit\n");
#endif
	if ((s=getvar(&p))) d1.val = 1;
	else              d1.val = 0;  
	d1.type  = CONST;
	d1.vtype = NUMBER;
	pushm(d1);
}

/*------------------------------------------*/

#define MAXDIM 20

void popmul(datum *d)

/* pop a var (unevaluated) from stack */

{
   int i,npops;

  i = 0;
  d[i++] = popm();			/* var */
  switch (d[0].type) {
    case ARRAYVAL: 
    case LOCALARRVAL: 
    case ARGARRVAL: 
    d[i++] = popm();			/* number of dimensions */
    npops = i + d[1].argnum;		/* pop just dimensions */
    if (npops >= MAXDIM) execerror("popmul: too many dimensions",""); 
    for (; i<npops; i++) {
      d[i] = popm();
    }
    break;
    default:
    break;
  }
}

/*----------------------------------------------------*/

void pushmul(datum *d)

/* push a var (unevaluated) onto stack (backwards) */
{
   int i,npops;

  switch (d[0].type) {
    case ARRAYVAL: 
    case LOCALARRVAL: 
    case ARGARRVAL: 
    npops = d[1].argnum + 2;		/* pop var, ndim, dimensions */
    if (npops >= MAXDIM) execerror("pushmul: too many dimensions",""); 
    for (i=npops-1; i>=0; i--) {
      pushm(d[i]);
    }
    break;
    default:
      pushm(d[0]);
    break;
  }
}

/*----------------------------------------------------*/

void popone(void)
{
  	datum d1[MAXDIM], d2[MAXDIM];

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"popone\n");
#endif
  	popmul(d1);
  	popmul(d2);
  	pushmul(d1);
}

/*------------------------------------------*/

int getarr(datum *d, double **parr, int *siz)

/* Find the beginning of an array and its size. */

{
     int ndim,arraysiz;
     int i;
     ncarray *apnt;

if (d->type == ARRAY || d->type == ARRAYCONST) {
   apnt = d->arrp;
  *parr = apnt->arr;
  ndim = apnt->ndim;                           /* get # of dimensions */
  arraysiz = 1;
  for (i=0; i<ndim; i++) {
    arraysiz *= apnt->dim[i];                  /* mult by next lower dim */
  }
  *siz = arraysiz;
  return (1);
 }
 else return (0);
}

/*------------------------------------------*/

void getsymarr(ncarray *apnt, double **parr, int *siz)

/* Find the beginning of an array and its size. */

{
     int ndim,arraysiz;
     int i;

  *parr = apnt->arr;
  ndim = apnt->ndim;                           /* get # of dimensions */
  arraysiz = 1;
  for (i=0; i<ndim; i++) {
    arraysiz *= apnt->dim[i];                  /* mult by next lower dim */
  }
  *siz = arraysiz;
}

/*------------------------------------------*/

void xsizeof (void)

{
   datum d;
   double *p;
   int size;

   evalvar();
   d = popm();
   if (getarr (&d,&p,&size))
     d.val = (double)size;
    else 
     d.val = 1;
   d.type=CONST;
   d.vtype=NUMBER;
   pushm(d);
}

/*------------------------------------------*/

datum copyarr(datum *src, double **p, int *siz)

/* Copy array base for result of arithmetic array operation onto
   the top of the stack. Return pointer to array and its size. */

{
   datum d={0};
   Symbol sp={0};
   ncarray *apnt;
   int i, arraysiz,ndim;
   static char tnam[] = "temparr";

   apnt = src->arrp;
   ndim = apnt->ndim;                           /* get # of dimensions */
   arraysiz = 1;
   d.type = VAR;
   for (i=0; i<ndim; i++) {
     d.val = apnt->dim[i];
     arraysiz *= int(d.val);
     pushm (d);
   }
   *siz = arraysiz;
   sp.name = tnam; 
   sp.type = ARRAY; 
   darr2(&sp,ndim);
   *p = sp.arrp->arr;
   d.arrp = sp.arrp; 
   d.type = sp.type;
   d.vtype = src->vtype;
   return (d);
}

/*----------------------------------------------------*/

double f_add (double *val1, double *val2, int vtype);
double f_mul (double *val1, double *val2, int vtype);


void dobinaryop(void(faa)(int size, double *p1,double *p2,double *p3,int vtype,
			double(f_logic)(double*,double*,int),int),
		void(fav)(int size, double *p1, datum *d2, double *p3, 
			double(f_logic)(double*,double*,int),int),
		void(fva)(int size, datum *d1, double *p2, double *p3, 
			double(f_logic)(double*,double*,int),int),
		void(fvv)(datum *d1,datum *d2,double(f_logic)(double*,double*,int),int),
		double(f_logic)(double*,double*,int),
		int changetype
		)
{
#define XBUFSIZ 64
	int vtype;
	datum d1,d2,d3,d4;
	ncarray *parr1=NULL, *parr2=NULL;
	char *buf;
	
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"dobinaryop\n");
#endif
	d2 = popm();
	d1 = popm();
	//if ((d1.vtype==STRING || d2.vtype==STRING) && (d1.vtype != d2.vtype)) {
	//		execerror (diftyp,0);
	//} 
	if (d1.type == ARRAY && d2.type==ARRAY) {
		double *p1,*p2,*p3;
		int siz1,siz2,siz3;

	  if (d1.atype==TEMP) parr1 = d1.arrp;	
	  if (d2.atype==TEMP) parr2 = d2.arrp;	
          getarr(&d1, &p1, &siz1);
          getarr(&d2, &p2, &siz2);
	  if (siz1 != siz2) {
 	    execerror("arrays not same size", (char *)0);
	    pushm(d1);
	    return;
	  }
	  vtype = d1.vtype;
	  if (changetype==0) {
	     if (d1.vtype==STRING && d2.vtype==NUMBER) {
	        vtype = STRNUM;
	        d3 = copyarr(&d1, &p3, &siz3);
	     } else if (d1.vtype==NUMBER && d2.vtype==STRING) {
	        d3 = copyarr(&d2, &p3, &siz3);
	        vtype = NUMSTR;
	     }
	     else d3 = copyarr(&d1, &p3, &siz3);
	  }  else d3 = copyarr(&d1, &p3, &siz3);
	  pushm(d3);
          faa(siz1,p1,p2,p3,vtype,f_logic,changetype);     	  
	} 
	else if (d1.type == ARRAY && d2.type!=ARRAY) {
		double *p1,*p3;
		int siz1,siz3;
	
	  if (d1.atype==TEMP) parr1 = d1.arrp;	
          getarr(&d1, &p1, &siz1);
          d3 = copyarr(&d1, &p3, &siz3);
	  pushm(d3);
	  d4 = d2; 
	  if (d1.vtype==STRING && d2.vtype==NUMBER) {
	     buf = emalloc(XBUFSIZ);
             sprintf(buf,"%g",d2.val);
             d4.str = buf;
             d4.vtype=STRING;
          } 
	  else if (d1.vtype==STRING && d2.vtype==LITCHAR) {
	     buf = emalloc(XBUFSIZ);
	     buf[0] = (char)d2.val;
	     buf[1] = 0;
             d4.str = buf;
             d4.vtype=STRING;
          } 
	  else if (d1.vtype==NUMBER && d2.vtype==STRING) {
             d4.val = atof(d2.str);
             d4.vtype=NUMBER;
	  }
          fav(siz1,p1,&d4,p3,f_logic,changetype);
	} 
	else if (d1.type != ARRAY && d2.type==ARRAY) {
		double *p2,*p3;
		int siz2,siz3;
	
	  if (d2.atype==TEMP) parr2 = d2.arrp;	
          getarr(&d2, &p2, &siz2);
          d3 = copyarr(&d2, &p3, &siz3);
	  pushm(d3);
	  d4 = d1; 
	  if (d1.vtype==STRING && d2.vtype==NUMBER) {
             d4.val = atof(d1.str);
             d4.vtype=NUMBER;
          } 
	 // else if (d1.vtype==STRING && d2.vtype==LITCHAR) {
	 //    buf[0] = (char)d2.val;
	 //    buf[1] = 0;
         //    d2.str = buf;
         //    d2.vtype=STRING;
         // } 
	  else if (d1.vtype==NUMBER && d2.vtype==STRING) {
	     buf = emalloc(XBUFSIZ);
             sprintf(buf,"%g",d1.val);
             d4.str = buf;
             d4.vtype=STRING;
	  }
          fva(siz2,&d4,p2,p3,f_logic,changetype);

	} else {

          fvv(&d1,&d2,f_logic,changetype);     	  
	  d1.type = VAR;
	  pushm(d1);
        }
	if (parr1) {
		erarr(parr1);	/* if orig arrays were temp */
		 /* ncfprintf (stderr,"erasing array 1\n"); /* */
	}
	if (parr2) {
		erarr(parr2);
		/* ncfprintf (stderr,"erasing array 2\n"); /* */
	}
	d1 = popm();
	d1.atype = TEMP;
	if (changetype) d1.vtype = changetype;
	pushm(d1);
}

/*----------------------------------------------------*/

void f_aa (int size, double *p1, double *p2, double *p3, int vtype,
		double (f)(double*,double*,int), int changetype)
/* aa = array, array */

{
    int i;
    char **c1, **c2, **c3;
    double val;
    char *buf;

   if (vtype==STRING || vtype==NUMSTR || vtype==STRNUM) {
	c1 = (char**)p1;
	c2 = (char**)p2;
	c3 = (char**)p3;
  }
  for (i=0; i<size; i++) {
    if (*p1==LARGENUM || *p2==LARGENUM) {
      execerror("array",varuninit);
    }
    switch (vtype) {
    case NUMBER:
    case LITCHAR:
         *(p3++) = f(p1++, p2++, vtype);
	break;
    case STRING:
	 if (changetype==0) {
	    if (f==f_mul) {
               *(p3++) = 0;
               //*(p3++) = *(double*)(c1++);
	    } else {
               val = f((double*)(c1++), (double*)(c2++), vtype);
               *(c3++) = *(char**)&val;
	    }
	 } else {
               //*(p3++) = *(double*)(c1++);
               //*(p3++) = 0;
               *(p3++) = f((double*)(c1++), (double*)(c2++), vtype);
	 }
        break;
    case NUMSTR:		// convert array 2 to number
	 if (f==f_add) {
            val = atof(*c2++);
            *(p3++) = f(p1++, &val, NUMBER);
	 } else {
            val = f(p1++, (double*)(c2++), NUMBER);
            //val = f((double*)(c2++), p1++, STRING);
	    if (val==0) *(c3++) = (char *)nullstr;
            else *(c3++) = *(char**)&val;
	 }
        break;
    case STRNUM:		// convert array 2 to string
	 if (f==f_add) {
	     buf = emalloc(XBUFSIZ);
             sprintf(buf,"%g",*p2++);
             val = f((double*)(c1++), (double*)(&buf), STRING);
             *(c3++) = *(char**)&val;
	 } else {
             val = f((double*)(c1++), p2++, STRING);
             *(c3++) = *(char**)&val;
	 }
        break;
    }
  }
}
#undef XBUFSIZ

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void f_av (int size, double *p1, datum *d2, double *p3, 
		double(f)(double*,double*,int),int changetype)

/* av = array, variable */

{
    int i;
    double *val;
    char **c1, **c3;
    double retval;

  if (d2->vtype==STRING) {
	c1 = (char**)p1;
	c3 = (char**)p3;
  }
    
  val = &d2->val;
  if (*val==LARGENUM) {
	execerror ("f_av", varuninit);
  }
  for (i=0; i<size; i++) {
    if (*p1==LARGENUM) {
	execerror ("array", varuninit);
    }
    switch (d2->vtype) {
    case NUMBER:
    case LITCHAR:
        *(p3++) = f(p1++, val, d2->vtype);
	break;
    case STRING:
        retval = f((double*)(c1++), val, d2->vtype);
	if (changetype==NUMBER) 
           *(p3++) = retval;
	else
           *(c3++) = *(char**)&retval;
        break;
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void f_va (int size, datum *d1, double *p2, double *p3, 
		double(f)(double*,double*,int), int changetype)
/* va = variable, array */

{
    int i;
    double *val;
    char **c2, **c3;
    double retval;

  if (d1->vtype==STRING) {
	c2 = (char**)p2;
	c3 = (char**)p3;
  }
    
  val = &d1->val;
  if (*val==LARGENUM) {
	execerror ("f_va", varuninit);
  }
  for (i=0; i<size; i++) {
    if (*p2==LARGENUM) {
	execerror ("array", varuninit);
    }
    switch (d1->vtype) {
    case NUMBER:
    case LITCHAR:
        *(p3++) =  f(val, p2++, d1->vtype);
	break;
    case STRING:
        retval = f(val, (double*)(c2++), d1->vtype);
	if (changetype==NUMBER) 
          *(p3++) = retval;
	else
          *(c3++) = *(char**)&retval;
        break;
    }
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void f_vv (datum *d1, datum *d2, double(f)(double*,double*,int), int changetype)

/* vv = variable, variable */

#define XBUFSIZ 256

{
   char buf[XBUFSIZ];

  if (d1->val==LARGENUM || d2->val==LARGENUM) {
	execerror ("f_vv", varuninit);
  }
  if (d1->vtype==STRING && d2->vtype==NUMBER) {
     sprintf(buf,"%g",d2->val);
     d2->str = buf;
     d1->val = f(&d1->val, &d2->val, d1->vtype);
  } else {
     d1->val = f(&d1->val, &d2->val, d1->vtype);
  }
}
#undef XBUFSIZ

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

static char *null=(char *)"";

double f_add (double *val1, double *val2, int vtype)

{
    int len1, len2;
    char *strp;
    double val;

   switch (vtype) {
   case NUMBER:
   case LITCHAR:
	val = (*val1 + *val2);
	break;
   case STRING:
	if (*(char**)val1==NULL) val1 = (double*) &null;
	if (*(char**)val2==NULL) val2 = (double*) &null;
	len1 = strlen(*(char**)val1);
	len2 = strlen(*(char**)val2);
	strp = emalloc(len1+len2+1);
	strcpy (strp,*(char**)val1);
	strcat (strp,*(char**)val2);
	val =  *(double*)&strp;
	break;
   } 
   return val;
}

/*-----------------------------------------------------*/

void binaryop (double (f)(double*, double*, int))

{
  dobinaryop(f_aa, f_av, f_va, f_vv, f, 0);
}

/*-----------------------------------------------------*/

void logic_binaryop (double (f)(double*, double*, int))

{
  dobinaryop(f_aa, f_av, f_va, f_vv, f, NUMBER);
}

/*------------------------------------------*/

void add(void)
{

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"add\n");
#endif

 binaryop(f_add);
}

/*----------------------------------------------------*/

double f_sub (double *val1, double *val2, int vtype)

{
    double val;

   switch (vtype) {
   case NUMBER:
   case LITCHAR:
	val = (*val1 - *val2);
	break;
   case STRING:
	if (*(char**)val1==NULL) val1 = (double*) &null;
	if (*(char**)val2==NULL) val2 = (double*) &null;
	val =  *(double*)val1;
	break;
   } 
   return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void sub(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"sub\n");
#endif

  binaryop(f_sub);
}

/*----------------------------------------------------*/

double f_mul (double *val1, double *val2, int vtype)

{
    double val;

   switch (vtype) {
   case NUMBER:
	val = (*val1 * *val2);
	break;
   case LITCHAR:
	val = *val1;
	break;
   case STRING:
	val =  *val1 * *val2;
	break;
   } 
   return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void mul(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"mul\n");
#endif

  binaryop(f_mul);
}

/*----------------------------------------------------*/

void mod_aa (int size, double *p1, double *p2, double *p3, int vtype,
				double(f)(double*,double*,int),int changetype)

{
    int i;
    int denom;

  for (i=0; i<size; i++) {
    if (*p1==LARGENUM || *p2==LARGENUM) 
	execerror ("stopping in mod, array ", varuninit); 
    if ((denom=int(*(p2++)))==0) execerror("stopping in mod, ",divzero);
    *(p3++) = int(*(p1++)) % denom;
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void mod_av (int size, double *p1, datum *d2, double *p3,
				double(f)(double*,double*,int),int changetype)

{
    int i;
    int denom;

  if (d2->val==LARGENUM) 
	execerror ("stopping in mod, ", varuninit); 
  if (d2->val == 0.0)
     execerror("stopping in mod, ",divzero);
  denom = int(d2->val);
  for (i=0; i<size; i++) {
    if (*p1==LARGENUM) 
	execerror ("stopping in mod, array ", varuninit); 
   *(p3++) = int(*(p1++)) % denom; 
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void mod_va (int size, datum *d1, double *p2, double *p3,
				double(f)(double*,double*,int),int changetype)

{
    int i;
    int denom;

  if (d1->val==LARGENUM) 
	execerror ("stopping in mod, ", varuninit); 
  for (i=0; i<size; i++) {
    if (*p2==LARGENUM) 
	execerror ("stopping in mod, array ", varuninit); 
    if ((denom=int(*(p2++)))==0) execerror("stopping in mod, ",divzero);
    *(p3++) = int(d1->val) % denom; 
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void mod_vv (datum *d1, datum *d2, double(f)(double*,double*,int),int changetype)

{
   if (((int)d2->val) == 0.0) {
     execerror("stopping in mod, ",divzero);
     d1->val = LARGENUM;
     return;
   }
   if (d1->val==LARGENUM || d2->val==LARGENUM) {
	execerror ("mod argument", varuninit);
   }
   switch (d1->vtype) {
   case NUMBER:
   case LITCHAR:
	d1->val = int(d1->val) % int(d2->val);
	break;
   case STRING:
	break;
   }
}

/*------------------------------------------*/

void mod(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"mod\n");
#endif

  dobinaryop(mod_aa, mod_av, mod_va, mod_vv, NULL, 0);
}

/*----------------------------------------------------*/

void div_aa (int size, double *p1, double *p2, double *p3, int vtype,
				double(f)(double*,double*,int),int changetype)

{
    int i;
    double denom;

  for (i=0; i<size; i++) {
    if (*p1==LARGENUM || *p2==LARGENUM) 
	execerror ("stopping in div, array", varuninit); 
    if ((denom= *(p2++))==0.0) execerror("stopping in div, ",divzero);
    *(p3++) = *(p1++) / denom;
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void div_av (int size, double *p1, datum *d2, double *p3,
				double(f)(double*,double*,int),int changetype)

{
    int i;
    double val;
    
  if (d2->val==LARGENUM) 
	execerror ("stopping in div, ", varuninit); 
  if (d2->val == 0.0)
	execerror("division by zero", (char *)0);
  val = d2->val;
  for (i=0; i<size; i++) {
    if (*p1==LARGENUM) 
	execerror ("stopping in div, array ", varuninit); 
    *(p3++) = *(p1++) / val;
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void div_va (int size, datum *d1, double *p2, double *p3,
				double(f)(double*,double*,int),int changetype)

{
    int i;
    double denom,val;

  val = d1->val;
  if (val==LARGENUM) 
	execerror ("stopping in div, ", varuninit); 
  for (i=0; i<size; i++) {
    if (*p2==LARGENUM) execerror ("stopping in div, array ", varuninit); 
    if ((denom= *(p2++))==0.0) execerror("stopping in div, ",divzero);
    *(p3++) =  val / denom;
  }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void div_vv (datum *d1, datum *d2, double(f)(double*,double*,int),int changetype)

{
   if (d2->val == 0.0) {
	execerror("division by zero", (char *)0);
	d1->val = LARGENUM;
	return;
   }
   if (d1->val==LARGENUM || d2->val==LARGENUM) {
	execerror ("stopping in div, ", varuninit);
   }
   switch (d1->vtype) {
   case NUMBER:
   case LITCHAR:
	d1->val /= d2->val;
	break;
   case STRING:
	break;
   }
}

/*------------------------------------------*/
 
void xdiv(void)
{

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"div\n");
#endif

  dobinaryop(div_aa, div_av, div_va, div_vv, NULL, 0);
}

/*----------------------------------------------------*/

double logic_and (double *val1, double *val2, int vtype)
{
    static char *null=(char *)"";
    double val;

   switch (vtype) {
   case NUMBER:
   case LITCHAR:
	val = (int(*val1) && int(*val2));
	break;
   case STRING:
	if (*(char**)val1==NULL) val1 = (double*) &null;
	if (*(char**)val2==NULL) val2 = (double*) &null;
	val = (strlen(*(char**)val1)>0 && strlen(*(char**)val2)>0);
	break;
   } 
   return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void xand(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"xand\n");
#endif

  logic_binaryop (logic_and);
}

/*----------------------------------------------------*/

double logic_or (double *val1, double *val2, int vtype)
{
   static char *null=(char *)"";
   double val;

   switch (vtype) {
   case NUMBER:
   case LITCHAR:
	val = (int(*val1) || int(*val2));
	break;
   case STRING:
	if (*(char**)val1==NULL) val1 = (double*) &null;
	if (*(char**)val2==NULL) val2 = (double*) &null;
	val = (strlen(*(char**)val1)>0 || strlen(*(char**)val2)>0);
	break;
   } 
   return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void orx(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"orx\n");
#endif

  logic_binaryop (logic_or);
}

/*----------------------------------------------------*/

double logic_xor (double *val1, double *val2, int vtype)
{
    static char *null=(char *)"";
    double val;

   switch (vtype) {
   case NUMBER:
   case LITCHAR:
	val = (int(*val1) ^ int(*val2));
	break;
   case STRING:
	if (*(char**)val1==NULL) val1 = (double*) &null;
	if (*(char**)val2==NULL) val2 = (double*) &null;
	val = (strlen(*(char**)val1)>0 ^ strlen(*(char**)val2)>0);
	break;
   } 
   return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void xxor(void)
{
#ifdef CODEBUG
    // setdebug(1);
    if (codebug) ncfprintf (stderr,"xor\n");
#endif

  logic_binaryop (logic_xor);
}

/*----------------------------------------------------*/

double logic_bitand (double *val1, double *val2, int vtype)
{
    double val;

   switch (vtype) {
   case NUMBER:
   case LITCHAR:
	val = (int(*val1) & int(*val2));
	break;
   case STRING:
	val = 0;
	break;
   } 
   return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void bitand_x(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"bitand_x\n");
#endif

  logic_binaryop (logic_bitand);
}

/*----------------------------------------------------*/

double logic_bitor (double *val1, double *val2, int vtype)
{
    double val;

   switch (vtype) {
   case NUMBER:
   case LITCHAR:
	val = (int(*val1) | int(*val2));
	break;
   case STRING:
	val = 0;
	break;
   } 
   return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void bitor_x(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"bitor_x\n");
#endif

  logic_binaryop (logic_bitor);
}

/*----------------------------------------------------*/

void domonoop  (void(*fa) (int size, double *p1, double *p2, int type),
                void(*fv)(double *p, int type)
                )
{
        datum d={0},d1;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"domonoop\n");
#endif
 	d = popm();
        if (d.type == ARRAY) {
                double *p1,*p2;
                int siz1,siz2;
          d.vtype = NUMBER;
          getarr(&d, &p1, &siz1);
          d1 = copyarr(&d, &p2, &siz2);
	  pushm(d1);
          if (fa) fa(siz1,p1,p2,d.vtype);
        }
        else {
          if (fv) fv(&d.val,d.vtype);
          d.type = VAR;
          d.vtype = NUMBER;
          pushm(d);
        }
}

/*----------------------------------------------------*/


void dup(void)

/* duplicate top of stack */
/*  for type=ARRAYVAL, need to duplicate dimensions for array access as well */

{
   datum dd[MAXDIM];

  popmul(dd);
  pushmul(dd);
  pushmul(dd);
}

/*------------------------------------------*/

void dupcopy(void)

/* Duplicate top of stack and copy array */
/* Only used on evaluated expressions */

{
   datum d1,d2;
   int i,size;
   double *p1, *p2;

  d1 = popm();
  pushm(d1);
  if (d1.type==ARRAY) {
    p1 = d1.arrp->arr;
    d2 = copyarr(&d1,&p2,&size);
    pushm(d2);
    for (i=0; i<size; i++) {
      *(p2+i) = *(p1+i);
    } 
  }
  else pushm(d1);
}

/*----------------------------------------------------*/

void push1(void)

/* push the value 1 onto stack */

{
    datum d1;
  d1.type = VAR;
  d1.vtype = NUMBER;
  d1.atype = GLOBAL;
  d1.val   = 1;
  pushm(d1); 
}

/*------------------------------------------*/

void postop(void(*op)())

{
  datum d[MAXDIM]; 
			/* what's left after: */
  dup();		/* var var */
  popmul(d);		/* var      (var popped) */
  evalvar();		/* val      (var popped */
  dupcopy();		/* val valx (var popped) */
  popone();		/* valx     (var popped) */
  pushmul(d);		/* valx var */	
  dup();		/* valx var var */
  evalvar();		/* valx var val */
  push1();		/* valx var val 1 */
  op();			/* valx var val+1 */
  assign();		/* valx var=val+1 */
  d[0]=popm();		/* valx */
}

/*------------------------------------------*/

void postinc(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"postinc\n");
#endif

  postop(add);
}

void postdec(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"postdec\n");
#endif
  postop(sub);
}

/*------------------------------------------*/

void preop(void(*op)())

{
			/* what's left after: */
  dup();		/* var var */
  evalvar();		/* var val */
  push1();		/* var val 1 */
  op();			/* var val+1 */
  assign();		/* var=val+1 */
}

void preinc(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"preinc\n");
#endif
	preop(add);
}

void predec(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"predec\n");
#endif
	preop(sub);
}

void varexp(double **p1, double *p2, short int *p3)

/* Check for valid variable and
return variable pointer and expression value */

{
	datum d1;
	pdatum p={0};
	Symbol *s;

        d1 = popm();                    /* get expression to be assigned */
	if ((s=getvar(&p))) 		/* get variable to set value */
           execerror (s->name, varuninit);

        if (*p.vtype == 0 ||            /* check to make sure same type */
            *p.vtype == d1.vtype) {
          *p.vtype = d1.vtype;          /* assign value's type */
	  *p3 = *p.vtype;
          *p1 = p.val;			/* assign var address */
          switch (d1.vtype) {
            case NUMBER:
            case LITCHAR:
		*p2 = d1.val;		/* assign numerical value */
                break;
            case STRING:		/* don't assign string value here. */
                break;
          }             
	}
	else {
/*              ncfprintf (stderr,"d1 %d  d2 %d\n",d1.vtype,d2.vtype); */
/*                execerror (diftyp,d->sym->name); */
                execerror (diftyp,0);

	}
}

void opeq(void(*op) ())
{
   datum d[MAXDIM];

  popmul(d);	/* pop expression */
  dup();
  evalvar();
  pushmul(d);
  op();
  assign();
}

void addeq(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"addeq\n");
#endif
   opeq(add);
}

void subeq(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"subeq\n");
#endif
   opeq(sub);
}

void muleq(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"muleq\n");
#endif
   opeq(mul);
}

void diveq(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"diveq\n");
#endif
   opeq(xdiv);
}

void andeq(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"andeq\n");
#endif
   opeq(bitand_x);
}

void oreq(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"oreq\n");
#endif
   opeq(bitor_x);
}

/*----------------------------------------------------*/

double logic_gt (double *val1, double *val2, int vtype)
{
   double val;
   static char *null=(char *)"";

	switch (vtype) {
	  case NUMBER:
	  case LITCHAR:
		val = (double)(*val1 > *val2);
		break;
	  case STRING:
		if (*(char**)val1==NULL) val1 = (double*) &null;
		if (*(char**)val2==NULL) val2 = (double*) &null;
		val = (strcmp(*(char**)val1,*(char**)val2) > 0);
		break;
	}
	return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void gt(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"gt\n");
#endif

  logic_binaryop (logic_gt);
}

/*----------------------------------------------------*/

double logic_lt (double *val1, double *val2, int vtype)
{
   double val;
   static char *null=(char *)"";

	switch (vtype) {
	  case NUMBER:
	  case LITCHAR:
		val = (double)(*val1 < *val2);
		break;
	  case STRING:
		if (*(char**)val1==NULL) val1 = (double*) &null;
		if (*(char**)val2==NULL) val2 = (double*) &null;
		val = (strcmp(*(char**)val1, *(char**)val2) < 0);
		break;
	}
	return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void lt(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"lt\n");
#endif

  logic_binaryop (logic_lt);
}

/*----------------------------------------------------*/

double logic_ge (double *val1, double *val2, int vtype)
{
    double val;
    static char *null=(char *)"";

	switch (vtype) {
	  case NUMBER:
	  case LITCHAR:
		val = (double)(*val1 >= *val2);
		break;
	  case STRING:
		if (*(char**)val1==NULL) val1 = (double*) &null;
		if (*(char**)val2==NULL) val2 = (double*) &null;
		val = (strcmp(*(char**)val1,*(char**)val2) >= 0);
		break;
	}
	return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void ge(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"ge\n");
#endif

  logic_binaryop (logic_ge);
}

/*----------------------------------------------------*/

double logic_le (double *val1, double *val2, int vtype)
{
    double val;
    static char *null=(char *)"";

	switch (vtype) {
	  case NUMBER:
	  case LITCHAR:
		val = (double)(*val1 <= *val2);
		break;
	  case STRING:
		if (*(char**)val1==NULL) val1 = (double*) &null;
		if (*(char**)val2==NULL) val2 = (double*) &null;
		val = (strcmp(*(char**)val1,*(char**)val2) <= 0);
		break;
	}
	return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void le(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"le\n");
#endif

  logic_binaryop (logic_le);
}

/*----------------------------------------------------*/

double logic_eq (double *val1, double *val2, int vtype)
{
    double val;
    static char *null=(char *)"";

	switch (vtype) {
	  case NUMBER:
	  case LITCHAR:
		val = (double)(*val1 == *val2);
		break;
	  case STRING:
		if (*(char**)val1==NULL) val1 = (double *) &null;
		if (*(char**)val2==NULL) val2 = (double *) &null;
		val = (strcmp(*((char**)val1), *(char**)val2) == 0);
		break;
	}
	return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void eq(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"eq\n");
#endif

  logic_binaryop (logic_eq);
}

/*---------------------------------------------------------*/

double logic_ne (double *val1, double *val2, int vtype)
{
    double val;
    static char *null=(char *)"";

	switch (vtype) {
	  case NUMBER:
	  case LITCHAR:
		val = (double)(*val1 != *val2);
		break;
	  case STRING:
		if (*(char**)val1==NULL) val1 = (double *) &null;
		if (*(char**)val2==NULL) val2 = (double *) &null;
		val = (strcmp(*(char**)val1, *(char**)val2) != 0);
		break;
	}
	return val;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

void ne(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"ne\n");
#endif

  logic_binaryop (logic_ne);
}

/*---------------------------------------------------------*/

void not_a(int size, double *p1, double *p2, int vtype)

{
    int i,cmp;

  for (i=0; i<size; i++,p1++,p2++) {
	switch (vtype) {
	  case NUMBER:
	  case LITCHAR:
                *p2 = (double)(*p1 == 0.0);
		break;
	  case STRING:
		cmp = strlen(*(char**)p1);
		*p2 = (double)(cmp==0);
		break;
	}
  }
}

/*---------------------------------------------------------*/

void not_v(double *p, int vtype)

{
    int cmp;
    static char *null=(char *)"";

	switch (vtype) {
	  case NUMBER:
	  case LITCHAR:
                *p = (double)(*p == 0.0);
		break;
	  case STRING:
		if (*(char**)p==NULL) cmp = 0;
		else cmp = strlen(*(char**)p);
		*p = (double)(cmp==0);
		break;
	}
}

/*---------------------------------------------------------*/

void xnot(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"xnot\n");
#endif
      domonoop(not_a,not_v);
}

/*---------------------------------------------------------*/

void neg_a(int size, double *p1, double *p2, int vtype)

{
    int i;

  if (vtype!=STRING)
  for (i=0; i<size; i++,p1++,p2++) {
     *p2 = -*p1;
  }
}

/*---------------------------------------------------------*/

void neg_v(double *p, int vtype)

{
	switch (vtype) {
	  case NUMBER:
	  case LITCHAR:
                *p = -*p;
		break;
	  case STRING:
		break;
	}
}

/*---------------------------------------------------------*/

void negate(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"negate\n");
#endif
      domonoop(neg_a, neg_v);
}

/*---------------------------------------------------------*/

void power(void)
{
	datum d1,d2;
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"power\n");
#endif
	d2 = popm();
	d1 = popm();
	d1 = ncpow(d1, d2);
	pushm(d1);
}

/*---------------------------------------------------------*/

void assign(void)
{
	datum d1,d3,*d;
	static pdatum p={0};
	double *p1,*p2;
	int i,ndim,siz1,siz2;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"assign\n");
#endif
	d1 = popm();			/* get expression to be assigned */
	d = &stack[stackp-1];		/* get address of variable   */
	if (d->type==0) {
	      execerror("assigment to nonvariable","");
	}
	getvar(&p);			/* get variable to set value */

	if (d1.vtype ==0 || *p.vtype == 0 || /* check to make sure same type */
	    *p.vtype == d1.vtype) {
	if (d1.type==ARRAY || d1.type==ARRAYCONST) {
          getarr(&d1, &p1, &siz1);
	  if (*p.type == ARRAY) {		/* if two arrays */

            getsymarr(*p.arrp, &p2, &siz2);
	    if (siz1 != siz2) {		/* make new array */
		double *zzz;  /* not used here */

	      if (*p.atype==ARG) {
		execerror("assign: arg array not same size", (char *)0);
		pushm(d1);
		return;
	      }	
	      if (d1.type!=ARRAYCONST || siz2 < siz1) { 
		  erarr(*p.arrp); 		/* erase original array */
	         *p.type = d1.type;
	      /* *p.atype = d1.atype; */   /* don't change LOCAL to GLOBAL */
	         *p.vtype = d1.vtype;
	         d3 = copyarr(&d1,&zzz,&siz1);
	         *p.arrp = d3.arrp;	/* copy to symbol's array pointer */
	         p2 = d3.arrp->arr;	/* get pointer to values */
	         (*p.sym)->arrp = d3.arrp;
	      }
            }  
	    if (d1.vtype==STRING) {
	        for (i=0; i<siz1; i++) {
	          *(((char**)p2)+i) = *(((char**)p1)+i);
	        } 
	    } else {
	        for (i=0; i<siz1; i++) {
	          *(p2+i) = *(p1+i);
	        } 
	    }
	    *p.vtype = d1.vtype;
	  }
	  else if (*p.type==UNDEF) { 	/* make new array */
             double *zzz; /* not used here */

	    *p.type = ARRAY;
	    *p.atype = d1.atype;
	    *p.vtype = d1.vtype;
	    d3 = copyarr(&d1,&zzz,&siz1);
	    *p.arrp = d3.arrp;		/* copy to symbol's array pointer */
	    p2 = d3.arrp->arr;		/* get pointer to values */
	    if (d1.vtype==STRING) {
	        for (i=0; i<siz1; i++) {
	          *(((char**)p2)+i) = *(((char**)p1)+i);
	        } 
	    } else {
	        for (i=0; i<siz1; i++) {
	          *(p2+i) = *(p1+i);
	        } 
	    }
	  }
	  else { 	/* assignment of wrong type */
		execerror (diftyp,"");
		pushm(d1);
		return;
	  }
	  if (d1.atype==TEMP) {
		erarr (d1.arrp);
		/* ncfprintf (stderr,"erasing array\n"); /* */
	  }
	 
	}	/* end of (d1.type == ARRAY) */
	else {
	  if (*p.type == ARRAY) {		/* if different type */
		execerror (diftyp,"");
		pushm(d1);
	   	return; 
	  }
	  *p.vtype = d1.vtype;		/* assign value's type */
	  switch (d1.vtype) {
	  case NUMBER:
	  case LITCHAR:
		*p.val = d1.val;		/* assign value */
		break;
	  case STRING:
		if (d1.val==LARGENUM) {
		   *p.val = d1.val;
		}
		else {
		   if ((*p.str=(char *)emalloc(strlen(d1.str)+1))==NULL) {  /* make space for new copy */
		         execerror("No space for string copy","");	    /*  extra zero at end */
		   } else
		   strcpy(*p.str, (const char *)d1.str);      /* copy string value */
		}

		break;
	  }		
	  if (*p.type==UNDEF) *p.type=VAR;	/* make sure it's a variable*/
	 }

	}  /* if (typ==) */

	else {
         if (d->sym) {
	    if (d->argnum >0 && d->argnum < 100) {
			char tbuf[40];
		switch (d->type) {
		 case ARG:
		 case LOCALVAR:
		  sprintf (tbuf,"local var %d",d->argnum);
		  break;
		 case LOCALARR:
		  sprintf (tbuf,"local array var %d",d->argnum);
		  break;
	         }
	         execerror("Assignment of incorrect type for",tbuf);
	    } 
	    else {

	     // included this code instead of execerror() below to allow assignment of diff type
/*		
 	     switch (d1.vtype) {
	       case NUMBER:
	  	    *p.val = d1.val;
		    *p.vtype=NUMBER;
		    break;
	       case LITCHAR:
	  	    *p.val = d1.val;
		    *p.vtype=LITCHAR;
		    break;
	       case STRING:
	 	    if (d1.val==LARGENUM) {
		       *p.val = d1.val;
		    }
		    else *p.str = (char *)d1.str;
	            *p.vtype=STRING;
		    break;
	     }		
	     if (*p.type==UNDEF) *p.type=VAR;
*/

	      execerror("Assignment of incorrect type for",d->sym->name);
	    }
	  }
	}
	pushm(d1);
}

/*---------------------------------------------------------*/

void print(void)	/* pop top value from stack, print it */
{
	datum d;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"print\n");
#endif
	d = popm();
	if (freadindex>=0) {			/* if inside fread, scanning tokens */
	       	fread_val = d.val;		/* interp_fread returns fread_val */
		return;
	}

#ifndef XSTIMM
	if (d.type == ARRAY) {
             int j, n, r;
             int ndim, ncols, nrows;
	     int arraysiz;
	     double *pnt;

          getarr(&d, &pnt, &arraysiz);
	  //ncfprintf (stdout,"\n");
	  ndim = d.arrp->ndim;
	  ncols = d.arrp->dim[ndim-1];
	  if (ndim>=2) nrows = d.arrp->dim[ndim-2];
	  else         nrows = 1;
	  for (n=0; n<arraysiz; ) {
	   for (r=0; r<nrows; r++) {
	   for (j=0; j<ncols; j++,n++) {
	    if ((d.vtype==NUMBER || d.vtype==0) && *(pnt+n)==LARGENUM) ncfprintf (stdout,uninit);
	    else
	     switch (d.vtype) {
	       case 0:
	       case NUMBER:
		    ncfprintf (stdout,"%9.7g ",*(pnt+n));
		    break;
	       case LITCHAR:
		    ncfprintf (stdout,"%c ",(((char*)pnt)+n));
		    break;
	       case STRING:
		    if (**(((char**)pnt)+n)!=0)
		      ncfprintf(stdout,"%s ",*(((char **)pnt)+n));
		    break;
	      }
	    }
	    ncfprintf (stdout,"\n");
	   }
	   //ncfprintf (stdout,"\n");
	  }

	}
	else
	if (d.val==LARGENUM) ncfprintf (stdout,uninit);
	else
	 switch (d.vtype) {
	   case 0:
	   case NUMBER:
		 ncfprintf (stdout,"\t%.8g\n",d.val);
		 break;
	   case LITCHAR:
		 ncfprintf (stdout,"\t%c\n",(char)d.val);
		 break;
	   case STRING:
		 ncfprintf (stdout,"\t%s\n",d.str);
		 break;
	 }
#endif
}

/*---------------------------------------------------------*/

void prexpr(void)	/* print value of expression: no tab */
{
	datum d;
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"prexpr\n");
#endif
	d = popm();
#ifndef XSTIMM
	if (d.type == ARRAY) {
             int j, n, r;
             int ndim, ncols, nrows;
	     int arraysiz;
	     double *pnt;

          getarr(&d, &pnt, &arraysiz);
	  //ncfprintf (stdout,"\n");
	  ndim = d.arrp->ndim;
	  ncols = d.arrp->dim[ndim-1];
	  if (ndim>=2) nrows = d.arrp->dim[ndim-2];
	  else         nrows = 1;
	  for (n=0; n<arraysiz; ) {
	   for (r=0; r<nrows; r++) {
	   for (j=0; j<ncols; j++,n++) {
	    if ((d.vtype==NUMBER || d.vtype==0) && *(pnt+n)==LARGENUM) ncfprintf (stdout,uninit);
	    else
	     switch (d.vtype) {
	       case 0:
	       case NUMBER:
		    ncfprintf (stdout,"%9.7g ",*(pnt+n));
		    break;
	       case LITCHAR:
		    ncfprintf (stdout,"%c ",(((char)*pnt)+n));
		    break;
	       case STRING:
		    if (**(((char**)pnt)+n)!=0)
		      ncfprintf(stdout,"%s ",*(((char **)pnt)+n));
		    break;
	      }
	    }
	    ncfprintf (stdout,"\n");
	   }
	   //ncfprintf (stdout,"\n");
	  }

	}
	else
	if (d.val==LARGENUM) ncfprintf (stdout,uninit);
	else
	 switch (d.vtype) {
	  case 0:
	  case NUMBER:
		ncfprintf(stdout,"%.8g ", d.val);
		break;
	  case LITCHAR:
		ncfprintf(stdout,"%c ", (char)d.val);
		break;
	  case STRING:
		ncfprintf(stdout,"%s ", d.str);
		break;
	}
#endif
}

/*------------------------------------------------------*/

void crlf(void)
{
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"crlf\n");
#endif
#ifndef XSTIMM
	ncfprintf(stdout,"\n");
#endif
}

/*---------------------------------------------------------*/

void pfopen (void)		/* fopen (filename, mode) statement */

{

#define XBUFSIZ 256
	static datum  d1;
	const char *filnam,*filmode;
        FILE *ftemp;
        // char xbuf[XBUFSIZ];

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"pfopen\n");
#endif
        if (!checkstr(d1=popm())) {		/* get mode off stack */
	   return;
	}
	filmode = d1.str;
        if (!checkstr(d1=popm())) {		/* get filename off stack */
	   return; 
        }
        filnam = d1.str;			/* file to be read */
        if      (strcmp(filnam,"stdin") ==0) ftemp = stdin;
        else if (strcmp(filnam,"stdout")==0) ftemp = stdout;
        else if (strcmp(filnam,"stderr")==0) ftemp = stderr;
        else if ((ftemp=fopen(filnam,filmode)) == NULL) {
	/*
	   sprintf (xbuf,"fopen: can't open file '%.50s'.\n",filnam);
	   execerror ("warning,",xbuf);
	   return;
	*/
        }
	d1.str = (char *)ftemp;
	d1.vtype=NUMBER;
	push(d1);
#undef XBUFSIZ
}

/*---------------------------------------------------------*/

void pfclose (void)		/* fclose (filedesc) statement */

{
	static datum  d1;
        FILE *fildesc;
	int n;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"pfclose\n");
#endif
        if (!checknum(d1=popm())) {		/* get filedesc off stack */
	   return;
	}
	fildesc = (FILE *)d1.str;
        n = fclose(fildesc);
}

/*---------------------------------------------------------*/

#define MAXARGS 17			/* 1 format string, 16 args */

void pscanf(void) 			/* scanf statement */

{
#define STRLENGTH 100
	static datum  d1;
	static pdatum p[MAXARGS]={0};
	static double *v[MAXARGS];
	static char s[MAXARGS][STRLENGTH];
	int i,scase,indx,narg,rval;
	const char *fmt, *strp;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"pscanf\n");
#endif
	narg = (long int) *pc++;
	if (narg>MAXARGS) 
	   execerror("scanf: more than 17 arguments","");
	for (i=0; i<narg; i++) {
	   indx = MAXARGS-1-i;
	   getvar(&p[indx]);			/* pop var off stack */ 
	}					/*  save pointer to value */
	for (i=0; i<narg; i++) {
            p[i].val = p[MAXARGS-narg+i].val;
            p[i].vtype = p[MAXARGS-narg+i].vtype;
	}
	for (; i<MAXARGS; i++) {		/* zero the remaining ptrs */
            p[i].val = (double*)NULL;
            p[i].vtype = (short int*)NULL;
	}
	d1 = popm();				/* get format string */
	if (d1.vtype!=STRING)
	   execerror("scanf: string required for format","");

	fmt = d1.str;
	for (strp=fmt,i=0; i<narg; i++) {	/* check format string */
	   if (strp = strstr(strp,"%")) {
	     if (*(++strp)=='s') *p[i].vtype = STRING;
	     else {
		  char ch;
		while (ch=*++strp,isnum(ch)) ;	/* check for "%20s" */
	        if (*(strp)=='s') *p[i].vtype = STRING;
	     }
	   }
	}
	for (scase=i=0; i<narg; i++) {
	   switch  (*p[i].vtype) {
	    case 0: *p[i].vtype = NUMBER;
		     break;
	    case NUMBER:
	    case LITCHAR:
		break;
	    case STRING:
	        scase += 1 << i;
		break;
	   }
	}
	rval = 0;
	for (i=0; i<MAXARGS; i++) v[i] = p[i].val;

	/* ncfprintf (stderr,"scase %d\n",scase);   /* */
	switch (scase) {
	case 0:
	  rval=scanf(fmt,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 1:
	  rval=scanf(fmt,s[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 2:
	  rval=scanf(fmt,v[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 3:
	  rval=scanf(fmt,s[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 4:
	  rval=scanf(fmt,v[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 5:
	  rval=scanf(fmt,s[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 6:
	  rval=scanf(fmt,v[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 7:
	  rval=scanf(fmt,s[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	default:
          execerror("scanf: string variables must be in first 3\n","");
	  break;
	}
					/* assign the new values */
	for (i=0; i<narg; i++) {
	   int len;

	   switch (*p[i].vtype) {
	    case 0:
	    case NUMBER:
	    case LITCHAR:
		break;
	    case STRING:
		if ((len=strlen(s[i]))>0) {
		  *((char **)v[i]) = emalloc (len+1);
		  strcpy(*((char **)v[i]),s[i]);
	          *p[i].vtype = STRING;
		}
		break;
	   }
	}
   d1.val = rval;
   d1.vtype=NUMBER;
   push(d1);
}

/*---------------------------------------------------------*/

void pfscanf(void) 		/* fscanf (filedesc,"",args) statement */

{
#define STRLENGTH 100
	static datum  d1;
	static pdatum p[MAXARGS]={0};
	static double *v[MAXARGS];
	static char s[MAXARGS][STRLENGTH];
	int i,scase,indx,narg,rval,oflg;
	const char *fmt,*strp;
        FILE *ftemp;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"pfscanf\n");
#endif
	narg = (long int) *pc++;
	if (narg>MAXARGS) 
	   execerror("fscanf: more than 17 arguments","");
	for (i=0; i<narg; i++) {
	   indx = MAXARGS-1-i;
	   getvar(&p[indx]);			/* pop var off stack */ 
	}					/*  save pointer to value */
	for (i=0; i<narg; i++) {
            p[i].val = p[MAXARGS-narg+i].val;
            p[i].vtype = p[MAXARGS-narg+i].vtype;
	}
	for (; i<MAXARGS; i++) {		/* zero the remaining ptrs */
            p[i].val   = (double*)NULL;
            p[i].vtype = (short int*)NULL;
	}
	d1 = popm();				/* get format string */
	if (d1.vtype!=STRING)
	   execerror("scanf: string required for format","");
	fmt = d1.str;

	oflg = 0;
        d1 = popm();			/* get filedesc off stack */
	if (d1.vtype==STRING) {		/* must be file name, not descr */
	    const char *filnam;
	  filnam = d1.str;
          if      (strcmp(filnam,"stdin")==0) ftemp = stdin;
          else if ((oflg=1,ftemp=fopen(filnam,"r")) == NULL) {
		char xbuf[100];
	    sprintf (xbuf,"fscanf: can't open file '%.50s'.\n",filnam);
	    execerror ("warning,",xbuf);
	    return;
          }
	}
        else ftemp = (FILE *)((long int)d1.val);  /* existing file descriptor */
	if (!ftemp) {
           execerror ("fscanf","can't open file");
           return;
        }

	for (strp=fmt,i=0; i<narg; i++) {	/* check format string */
	   if (strp = strstr(strp,"%")) {
	     if (*(++strp)=='s') *p[i].vtype = STRING;
	     else {
		  char ch;
		while (ch=*++strp,isnum(ch)) ;	/* check for "%20s" */
	        if (*(strp)=='s') *p[i].vtype = STRING;
	     }
	   }
	}
	for (scase=i=0; i<narg; i++) {
	   switch   (*p[i].vtype) {
	    case 0:  *p[i].vtype = NUMBER;
		     break;
	    case NUMBER:
	    case LITCHAR:
		break;
	    case STRING:
	        scase += 1 << i;
		break;
	   }
	}
	rval = 0;
	for (i=0; i<MAXARGS; i++) v[i] = p[i].val;

	/* ncfprintf (stderr,"scase %d\n",scase);   /* */
	switch (scase) {
	case 0:
	rval=fscanf(ftemp,fmt,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 1:
	rval=fscanf(ftemp,fmt,s[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 2:
	rval=fscanf(ftemp,fmt,v[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 3:
	rval=fscanf(ftemp,fmt,s[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 4:
	rval=fscanf(ftemp,fmt,v[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 5:
	rval=fscanf(ftemp,fmt,s[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 6:
	rval=fscanf(ftemp,fmt,v[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 7:
	rval=fscanf(ftemp,fmt,s[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	default:
          execerror("fscanf: string variables must be in first 3\n","");
	  break;
	}
	if (oflg) fclose (ftemp);
					/* assign the new values */
	for (i=0; i<narg; i++) {
	   int len;

	   switch (*p[i].vtype) {
	    case 0:
	    case NUMBER:
	    case LITCHAR:
		break;
	    case STRING:
		if ((len=strlen(s[i]))>0) {
		  *((char **)v[i]) = emalloc (len+1);
		  strcpy(*((char **)v[i]),s[i]);
	          *p[i].vtype = STRING;
		}
		break;
	   }
	}
   d1.val = rval;
   d1.vtype=NUMBER;
   push(d1);
}

/*---------------------------------------------------------*/

void pfgets(void) 		/* fgets (var,n,filedesc) statement */

{
#define XBUFSIZ 1024
	static datum  d1;
	static pdatum p={0};
	int len,rval,oflg;
        FILE *ftemp;
        char *xbuf;
        size_t n;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"pfgets\n");
#endif
	oflg = 0;
        d1 = popm();			/* get filedesc off stack */
	if (d1.vtype==STRING) {		/* must be file name, not descr */
	    const char *filnam;
	  filnam = d1.str;
          if      (strcmp(filnam,"stdin")==0) ftemp = stdin;
          else if ((oflg=1,ftemp=fopen(filnam,"r")) == NULL) {
		char ebuf[100];
	    sprintf (ebuf,"fgets: can't open file '%.50s'.\n",filnam);
	    execerror ("warning,",ebuf);
	    return;
          }
	}
        else ftemp = (FILE *)((long int)d1.val);  /* existing file descriptor */
	if (!ftemp) {
           execerror ("ncfprintf","can't open file");
           return;
        }
        if (!checknum(d1=popm())) {		/* get N off stack */
	   return; 
        }
	n = int(d1.val);
	if (n > XBUFSIZ) n = XBUFSIZ;
	getvar(&p);				/* pop var off stack */ 
	xbuf = (char *)emalloc(XBUFSIZ*sizeof(char));
	if (getline (&xbuf,&n,ftemp)>0) {
	  rval = 1;
	  len = strlen(xbuf);
	  if ((*p.str != NULL) && (*p.vtype==STRING)) {
	     if (len > strlen(*p.str)) {
		efree ((char *)*p.str);
	  	*p.str = emalloc (len+1);
	     }
	  }
	  else *p.str = emalloc (len+1);
	  strcpy(*p.str,xbuf);
	  *p.vtype = STRING;
	}
	else rval = 0;
	if (oflg) fclose(ftemp);
	d1.val = rval;
	d1.vtype=NUMBER;
	push(d1);
#undef XBUFSIZ
}

/*---------------------------------------------------------*/

void pfgetc(void) 		/* fgetc (var,filedesc) statement */

{
#define XBUFSIZ 256
	static datum  d1;
	static pdatum p={0};
	int n,len,rval,oflg;
        FILE *ftemp;
        char xbuf[XBUFSIZ];

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"pfgetc\n");
#endif
	oflg = 0;
        d1 = popm();			/* get filedesc off stack */
	if (d1.vtype==STRING) {		/* is file name, not descr */
	    const char *filnam;
	  filnam = d1.str;
          if      (strcmp(filnam,"stdin")==0) ftemp = stdin;
          else if ((oflg=1,ftemp=fopen(filnam,"r")) == NULL) {
		char ebuf[100];
	    sprintf (ebuf,"fgetc: can't open file '%.50s'.\n",filnam);
	    execerror ("warning,",ebuf);
	    return;
          }
	}
        else ftemp = (FILE *)((long int)d1.val);  /* existing file descriptor */
	if (!ftemp) {
           execerror ("ncfprintf","can't open file");
           return;
        }
	getvar(&p);				/* pop var off stack */ 
	if ((xbuf[0]=fgetc (ftemp)) >= 0) {
	  rval = 1;
	  *p.val = xbuf[0];
	  *p.vtype = NUMBER;
	}
	else rval = 0;
	if (oflg) fclose(ftemp);
	d1.val = rval;
	d1.vtype=NUMBER;
	push(d1);
#undef XBUFSIZ
}

/*---------------------------------------------------------*/

void pfputc(void) 		/* fputc (var,filedesc) statement */

{
	static datum  d1;
	static pdatum p={0};
	int n,val,rval,oflg;
        FILE *ftemp;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"pfputc\n");
#endif
	oflg = 0;
        d1 = popm();			/* get filedesc off stack */
	if (d1.vtype==STRING) {		/* is file name, not descr */
	    const char *filnam;
	  filnam = d1.str;
          if      (strcmp(filnam,"stdin")==0) ftemp = stdin;
          else if ((oflg=1,ftemp=fopen(filnam,"r")) == NULL) {
		char ebuf[100];
	    sprintf (ebuf,"fputc: can't open file '%.50s'.\n",filnam);
	    execerror ("warning,",ebuf);
	    return;
          }
	}
        else ftemp = (FILE *)((long int)d1.val);  /* existing file descriptor */
	if (!ftemp) {
           execerror ("ncfprintf","can't open file");
           return;
        }
	getvar(&p);				/* pop var off stack */ 
	val = (int)*p.val;
	fputc(val,ftemp);
	if (oflg) fclose(ftemp);
	d1.val = val;
	d1.vtype=LITCHAR;
	push(d1);
}

/*---------------------------------------------------------*/

void psscanf(void) 		/* fscanf (filname,"") statement */

{
#define STRLENGTH 100
	static datum  d1;
	static pdatum p[MAXARGS]={0};
	static double *v[MAXARGS];
	static char s[MAXARGS][STRLENGTH];
	int i,scase,indx,narg,rval;
	const char *fmt, *strp, *sval;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"psscanf\n");
#endif
	narg = (long int) *pc++;
	if (narg>MAXARGS) 
	   execerror("sscanf: more than 17 arguments","");
	for (i=0; i<narg; i++) {
	   indx = MAXARGS-1-i;
	   getvar(&p[indx]);			/* pop var off stack */ 
	}					/*  save pointer to value */
	for (i=0; i<narg; i++) {
            p[i].val = p[MAXARGS-narg+i].val;
            p[i].vtype = p[MAXARGS-narg+i].vtype;
	}
	for (; i<MAXARGS; i++) {		/* zero the remaining ptrs */
            p[i].val   = (double*)NULL;
            p[i].vtype = (short int*)NULL;
	}
	checkstr(d1=popm());				/* get format string */
	fmt = d1.str;
	checkstr(d1=popm());				/* get buffer string */
	sval = d1.str;

	for (strp=fmt,i=0; i<narg; i++) {	/* check format string */
	   if (strp = strstr(strp,"%")) {
	     if (*(++strp)=='s') *p[i].vtype = STRING;
	     else {
		  char ch;
		while (ch=*++strp,isnum(ch)) ;	/* check for "%20s" */
	        if (*(strp)=='s') *p[i].vtype = STRING;
	     }
	   }
	}
	for (scase=i=0; i<narg; i++) {
	   switch (*p[i].vtype) {
	    case 0:  *p[i].vtype = NUMBER;
		     break;
	    case NUMBER:
	    case LITCHAR:
		break;
	    case STRING:
	        scase += 1 << i;
		break;
	   }
	}
	rval = 0;
	for (i=0; i<MAXARGS; i++) v[i] = p[i].val;

	/* ncfprintf (stderr,"scase %d\n",scase);   /* */
	switch (scase) {
	case 0:
	rval=sscanf(sval,fmt,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 1:
	rval=sscanf(sval,fmt,s[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 2:
	rval=sscanf(sval,fmt,v[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 3:
	rval=sscanf(sval,fmt,s[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 4:
	rval=sscanf(sval,fmt,v[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 5:
	rval=sscanf(sval,fmt,s[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 6:
	rval=sscanf(sval,fmt,v[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 7:
	rval=sscanf(sval,fmt,s[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	default:
          execerror("sscanf: string variables must be in first 3\n","");
	  break;
	}
					/* assign the new values */
	for (i=0; i<narg; i++) {
	   int len;

	   switch (*p[i].vtype) {
	    case 0:
	    case NUMBER:
	    case LITCHAR:
		break;
	    case STRING:
		if ((len=strlen(s[i]))>0) {
		  *((char **)v[i]) = emalloc (len+1);
		  strcpy(*((char **)v[i]),s[i]);
	          *p[i].vtype = STRING;
		}
		break;
	   }
	}
   d1.val = rval;
   d1.vtype=NUMBER;
   push(d1);
}

/*---------------------------------------------------------*/

void pprintf(void) 			/* printf statement */
{
	datum  d1,d[MAXARGS];
	double v[MAXARGS];
	char  *s[MAXARGS];
	int i,scase,indx,narg;
	const char *fmt;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"pprintf\n");
#endif
	narg = (long int) *pc++;
	if (narg>MAXARGS) 
	   execerror("printf: more than 17 arguments","");
	for (i=0; i<narg; i++)
	   d[MAXARGS-1-i] = popm();
	d1 = popm();
	if (d1.vtype!=STRING)
	   execerror("printf: string required for format","");
	for (scase=i=0; i<narg; i++) {
	   indx = MAXARGS-narg+i;
	   switch (d[indx].vtype) {
	    case 0:
	    case NUMBER:
	    case LITCHAR:
	   	if (d[indx].val==LARGENUM) {	/* print "uninit" */
		     v[i] = -LARGENUM;		/*  but must be number */
		}
		else v[i] = d[indx].val;	
		break;
	    case STRING:
	   	if (d[indx].val==LARGENUM) {	/* print "uninit" */
	   	 s[i] = uninit;
		}else {
	   	 s[i] = (char *)d[indx].str;	
		}
	        scase += 1 << i;
		break;
	   }
	}
	fmt = d1.str;
	/* ncfprintf (stderr,"scase %d\n",scase);   /* */
	switch (scase) {
	case 0:
	  ncfprintf(stdout,fmt,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 1:
	  ncfprintf(stdout,fmt,s[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 2:
	  ncfprintf(stdout,fmt,v[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 3:
	  ncfprintf(stdout,fmt,s[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 4:
	  ncfprintf(stdout,fmt,v[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 5:
	  ncfprintf(stdout,fmt,s[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 6:
	  ncfprintf(stdout,fmt,v[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 7:
	  ncfprintf(stdout,fmt,s[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	default:
	  execerror("printf: string variables must be in first 3\n","");
	  break;
	}
	fflush(stdout);
}

/*---------------------------------------------------------*/

void psprintf(void) 			/* sprintf statement */
{
	datum  d1,d[MAXARGS];
	pdatum  p={0};
	double v[MAXARGS];
	char  *s[MAXARGS];
	int i,scase,indx,narg,len;
        static char *cbuf;
	const char *fmt;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"psprintf\n");
#endif
	narg = (long int) *pc++;
	if (narg>MAXARGS) 
	   execerror("sprintf: more than 17 arguments","");
	for (i=0; i<narg; i++)
	   d[MAXARGS-1-i] = popm();
	d1 = popm();
	if (d1.vtype!=STRING)
	   execerror("printf: string required for format","");
        getvar(&p);				/* get pointer to output str */
        if (*p.type == UNDEF) *p.type = VAR;	/* set output type to var */
        *p.vtype = STRING;			/* set output vtype to string */
	for (scase=i=0; i<narg; i++) {
	   indx = MAXARGS-narg+i;
	   switch (d[indx].vtype) {
	    case 0:
	    case NUMBER:
	    case LITCHAR:
	   	if (d[indx].val==LARGENUM) {	/* print "uninit" */
		     v[i] = -LARGENUM;		/*  but must be number */
		}
		else v[i] = d[indx].val;	
		break;
	    case STRING:
	   	if (d[indx].val==LARGENUM) {	/* print "uninit" */
	   	 s[i] = uninit;
		}else {
	   	 s[i] = (char *)d[indx].str;	
		}
	        scase += 1 << i;
		break;
	   }
	}
	fmt = d1.str;
        cbuf = emalloc(strlen(fmt)+5000);
        cbuf[0]=0;
	/* ncfprintf (stderr,"scase %d\n",scase);   /* */
	switch (scase) {
	case 0:
	  sprintf(cbuf,fmt,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 1:
	  sprintf(cbuf,fmt,s[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 2:
	  sprintf(cbuf,fmt,v[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 3:
	  sprintf(cbuf,fmt,s[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 4:
	  sprintf(cbuf,fmt,v[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 5:
	  sprintf(cbuf,fmt,s[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 6:
	  sprintf(cbuf,fmt,v[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 7:
	  sprintf(cbuf,fmt,s[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	default:
	  execerror("sprintf: string variables must be in first 3\n","");
	  break;
	}
	if ((len=strlen(cbuf))>0) {		/* copy the new string to var */
	  *p.str = emalloc (len+1);
	  strcpy(*p.str,cbuf);
	}
	efree(cbuf);
}

/*---------------------------------------------------------*/

void pfprintf(void) 			/* ncfprintf (stderr,"") statement */
{
	datum  d1,d[MAXARGS];
	double v[MAXARGS];
	char  *s[MAXARGS];
	int i,scase,indx,narg,oflg;
	FILE *ftemp;
	const char *fmt;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"pfprintf\n");
#endif
	narg = (long int) *pc++;
	if (narg>MAXARGS) 
	   execerror("fprintf: more than 17 arguments","");
	for (i=0; i<narg; i++)
	   d[MAXARGS-1-i] = popm();
	if (!checkstr(d1=popm())) {	/* get format string off stack */
	    return;
	}
	fmt = d1.str;
	if (d1.vtype!=STRING)
	   execerror("fprintf: string required for format","");
	for (scase=i=0; i<narg; i++) {
	   indx = MAXARGS-narg+i;
	   switch (d[indx].vtype) {
	    case 0:
	    case NUMBER:
	    case LITCHAR:
	   	if (d[indx].val==LARGENUM) {	/* print "uninit" */
		     v[i] = -LARGENUM;		/*  but must be number */
		}
		else v[i] = d[indx].val;	
		break;
	    case STRING:
	   	if (d[indx].val==LARGENUM) {	/* print "uninit" */
	   	 s[i] = uninit;
		}else {
	   	 s[i] = (char *)d[indx].str;	
		}
	        scase += 1 << i;
		break;
	   }
	}
	oflg = 0;			/* file open flag */
        d1 = popm();			/* get filedesc off stack */
	if (d1.vtype==STRING) {		/* must be file name, not descr */
	    const char *filnam;
	  filnam = d1.str;
          if      (strcmp(filnam,"stdout")==0)  ftemp = stdout;
          else if (strcmp(filnam,"stderr")==0)  ftemp = stderr;
          else if ((oflg=1,ftemp=fopen(filnam,"w")) == NULL) {
		char ebuf[100];
	    sprintf (ebuf,"fprintf: can't open file '%.50s'.\n",filnam);
	    execerror ("warning,",ebuf);
	    return;
          }
	}
        else ftemp = (FILE *)((long int)d1.val);  /* existing file descriptor */
	if (!ftemp) {
           execerror ("fprintf","can't open file");
           return;
        }

	/* ncfprintf (stderr,"scase %d\n",scase);   /* */
	switch (scase) {
	case 0:
	ncfprintf(ftemp,fmt,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 1:
	ncfprintf(ftemp,fmt,s[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 2:
	ncfprintf(ftemp,fmt,v[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 3:
	ncfprintf(ftemp,fmt,s[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 4:
	  ncfprintf(ftemp,fmt,v[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 5:
	  ncfprintf(ftemp,fmt,s[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 6:
	  ncfprintf(ftemp,fmt,v[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 7:
	  ncfprintf(ftemp,fmt,s[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
        default:
	  execerror("fprintf: string variables must be in first 3\n","");
          break;
	}
	fflush(ftemp);
	if (oflg) fclose(ftemp);
}

/*---------------------------------------------------------*/

void txtf(void) 			/* textf statement */
{
	datum  d1,d[MAXARGS];
	double v[MAXARGS];
	char  *s[MAXARGS];
	const char *fmt;
	char *cbuf;
	int i,scase,indx,narg;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"txtf\n");
#endif
	narg = (long int) *pc++;
	if (narg>MAXARGS) 
	   execerror("textf: more than 17 arguments","");
	for (i=0; i<narg; i++)
	   d[MAXARGS-1-i] = popm();
	d1 = popm();
	if (d1.vtype!=STRING)
	   execerror("textf: string required for format","");
	fmt = d1.str;
	for (scase=i=0; i<narg; i++) {
	   indx = MAXARGS-narg+i;
	   switch (d[indx].vtype) {
	    case 0:
	    case NUMBER:
	    case LITCHAR:
	   	if (d[indx].val==LARGENUM) {	/* print "uninit" */
		     v[i] = -LARGENUM;		/*  but must be number */
		}
		else v[i] = d[indx].val;	
		break;
	    case STRING:
	   	if (d[indx].val==LARGENUM) {	/* print "uninit" */
	   	 s[i] = uninit;
		}else {
	   	 s[i] = (char *)d[indx].str;	
		}
	        scase += 1 << i;
		break;
	   }
	}
        cbuf = emalloc(strlen(fmt)+5000);
        cbuf[0]=0;

	/* ncfprintf (stderr,"scase %d\n",scase);   /* */
	switch (scase) {
	case 0:
	sprintf(cbuf,fmt,v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 1:
	sprintf(cbuf,fmt,s[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 2:
	sprintf(cbuf,fmt,v[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 3:
	sprintf(cbuf,fmt,s[0],s[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 4:
	  sprintf(cbuf,fmt,v[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 5:
	  sprintf(cbuf,fmt,s[0],v[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 6:
	  sprintf(cbuf,fmt,v[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	case 7:
	  sprintf(cbuf,fmt,s[0],s[1],s[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],v[10],v[11],v[12],v[13],v[14],v[15]);
	  break;
	default:
	  execerror("txtf: string variables must be in first 3\n","");
	  break;
	}
    gtext(cbuf);
    efree(cbuf);
}

/*---------------------------------------------------------*/
#define MAXARGGS 5

void gglabel(void) 			/* glabel statement */
{
	datum  d1,d[MAXARGGS];
	double v[MAXARGGS];
	char  *s[MAXARGGS];
	const char *fmt;
	char *cbuf;
	int i,indx,narg;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"gglabel\n");
#endif
	narg = (long int) *pc++;
	if (narg>4) 
	   execerror("glabel: more than 4 arguments","");
	for (i=0; i<narg; i++)
	   d[MAXARGGS-1-i] = popm();
	d1 = popm();
	if (d1.vtype!=STRING)
	   execerror("glabel: string required for format","");
	fmt = d1.str;
	for (i=0; i<narg; i++) {
	   indx = MAXARGGS-narg+i;
	   switch (d[indx].vtype) {
	    case 0:
	    case NUMBER:
	    case LITCHAR:
	   	if (d[indx].val==LARGENUM) {	/* print "uninit" */
		     v[i] = -LARGENUM;		/*  but must be number */
		}
		else v[i] = d[indx].val;	
		break;
	    case STRING:
	   	if (d[indx].val==LARGENUM) {	/* print "uninit" */
	   	 s[i] = uninit;
		}else {
	   	 s[i] = (char *)d[indx].str;	
		}
		break;
	   }
	}
    // fprintf (stderr,"%s %g %g %g %g\n",d1.str,v[0],v[1], v[2]);
    glabel(d1.str,v[0], v[1], v[2]);
}
#undef MAXARGGS

/*---------------------------------------------------------*/

int garrsiz(Symbol *apnt)
                
/*  find absolute size of array */

{
   int i,ndim,size;
   ncarray *ap;

   ap = apnt->arrp;			/* pointer to array */
   ndim = ap->ndim;
   for (size=1,i=0; i<ndim; i++) {
	size *= ap->dim[i];		/* multiply dims for size */
   }
   return size;
}

/*---------------------------------------------------------*/

void dlocarray(void)	/* allocate space for local array */
			/* stack var points to array */
{
	int argn,ndim;
	double *dpnt;
	Symbol sp;
	char tnam[40];
	datum d={0};

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"dlocarray\n");
#endif
	argn = (long int) *pc++;
	ndim = (long int) *pc++;                     /* number of args */
	sprintf (tnam,"local array arg %d\n",argn);
	sp.name = tnam; 
	sp.type = ARRAY; 
	dpnt = darr2(&sp,ndim);
	d.arrp  = sp.arrp;
	d.type = argn;			/* save arg number temporarily */
	d.atype = LOCAL; 
	d.vtype  = 0;
	pushm (d);
}

/*---------------------------------------------------------*/

void darray(void)	/* allocate space for array */

/* Allocate memory space for array and
 its descriptor block, and set descriptor block.
	
	Desciptor block:

	0	pointer to start of array
	1	number of dimensions in array
	2	dimension 1
	3	dimension 2
		etc.
*/

#define MAXDIM  20				/* maximum # of dimensions */

{
	int ndim;
	Symbol *sp;
	double *dpnt;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"darray\n");
#endif
	sp = (Symbol *)(*pc++);			/* pointer to array symbol */
	ndim = (long int) *pc++;			/* number of args */
	dpnt = darr2(sp,ndim);
}

/*---------------------------------------------------------*/

void getrange(int *dim1, double *low, double *high, double *incr)

/* get array initializer from range */

{
    datum d1,d2,d3;
    double range,lincr;
    int dimsav;

  dimsav = *dim1;
  switch (dimsav) {
   case 2:
    lincr = 1;
    break;
   case 3:
     d1=popm();
     lincr = d1.val;
     break;
  } 
  d2=popm();
  *high = d2.val;
  d3=popm();
  *low = d3.val;
  if (lincr==0) lincr = 1;
  range = *high - *low;
  if (range<0 && lincr>0) lincr = -lincr;
  *dim1 = int(abs(range / lincr) + 1 + 0.1* lincr/range);
  *incr = lincr;
  //ncfprintf (stderr,"%g %g %g %g\n",*low,*incr,*high, .1*lincr/range);
  pushm(d3);
  pushm(d2);
  if (dimsav==3) pushm(d1);
}

/*---------------------------------------------------------*/

double *darr2(Symbol *sp, int ndim)		/* define array */
{
	int i,dim1,darr[MAXDIM],x;
	datum d={0};
	long int size, dsize;
	ncarray *apnt;
	Symbol *arrsym;

	if (sp->type==VAR) {
		execerror(sp->name,"array name already defined");
		return((double*)NULL);
	}

        if (sp->name) arrsym = lookup(sp->name);
	if (arrsym!=NULL) {	/* if array has already been defined, erase it */
 	    if (arrsym->type==ARRAY && arrsym->val!=LARGENUM)
	       erarr(arrsym->arrp);
	}

	if (ndim==0) {	    /* if array size was defined by initialization */
	  double high, incr, low; 

	    ndim = 1;
	    dim1 = (long int) *pc++;	/* size of init list */
	    if (dim1<0) {
	      dim1 = -dim1;
	      getrange (&dim1, &low, &high, &incr);
	    }
	    d.val = dim1;	
	    pushm(d);
	}
	size = (ndim + 1) * sizeof(int) + sizeof (char *);
					/* allocate space for limits */
	apnt = sp->arrp = (ncarray *)emalloc((unsigned int)size);
	if (!apnt) {
		execerror("no space for array ",sp->name);
		return((double*)NULL);
	}
	apnt->ndim = ndim;
	size = 1;
	for (i=0; i<ndim; i++) {
	  d = popm();				/* get dimensions from stack */
	  x = (int)(d.val);			/* this line fixes bug */
	  darr[(MAXDIM-1)-i] = x;		/* save backwards in temp arr */
	  size *= x;				/* */
	}
	for (i=0; i<ndim; i++) {
	  apnt->dim[i] = darr[MAXDIM-ndim+i];	/* put first dim first */
	}
	if (size == 0) {
		// execerror("array dim is zero",sp->name);
		// return((double*)NULL);
	}
	if (size < 0) {
		execerror("array dim is negative",sp->name);
		exit(0);
		return((double*)NULL);
	}
	else if (size > MAXSIZE) {
		execerror("array too large:",sp->name);
		exit(0);
		return((double*)NULL);
	}
	dsize = size * sizeof(double);		
/*	ncfprintf(stderr,"size %u\n",dsize); 	/* */

					/* allocate space for array */
	apnt->arr = (double *)emalloc((unsigned int)dsize);	
	if (! apnt->arr) {
		execerror("run out of space for array ",sp->name);
		return ((double*)NULL);
	}
	for (i=0; i<size; i++) 
	  apnt->arr[i] = LARGENUM;	/* make array uninitialized */
	sp->type = ARRAY;
	//sp->vtype = NUMBER;
	return (apnt->arr);
}

/*------------------------------------------------------*/

ncarray *garray(Symbol *sp)	/* return a pointer to whole array */
               
{
	if (sp->arrp==NULL || *(double **)sp->arrp==NULL) {
	   execerror("array space not allocated for",sp->name);
	   return ((ncarray*)NULL);
	}
	return (sp->arrp);	/* return pointer to array */
}

/*------------------------------------------------------*/

double
*garrayval(Symbol *sp)	/* return a pointer to an array value */
               
{
	datum d;
	int i,narg,darr[MAXDIM],ndim;
	int *index,offset,x;
	ncarray *apnt;
	double *parr;

	d = popm();				/* get number of indexes */
	narg = d.argnum;
	index = darr+MAXDIM;
	for (i=0; i<narg; i++) {
	  d = popm();				/* get indexes from stack */
	  x = (int)(d.val);
	  *(--index) = x;			/* save backwards in temp arr*/
	}
	if (sp->arrp==NULL || *(double **)sp->arrp==NULL) {
	   execerror("array space not allocated for",sp->name);
	   return ((double *)NULL);
	}
	parr = *(double **)sp->arrp;	/* get pointer to array */

	apnt = sp->arrp;
	ndim = apnt->ndim;			/* get # of dimensions */
	if (narg != ndim) {
	   ncfprintf (stderr,"%s has ndim %d, tried %d\n",sp->name,ndim,narg);
	   execerror("array indices different than dimensions of",
					sp->name);
	   return ((double *)NULL);
	}
	index = darr+MAXDIM-narg;		/* pointer to first index */
	for (i=0; i<ndim; i++) {
	       int indx;
	  indx = *index++;
	  if (indx >= apnt->dim[i]) {
        ncfprintf 
	(stderr,"illegal index size for array '%s', dim %d, index=%d, max=%d\n",
		sp->name,i+1,indx,apnt->dim[i]-1);
            execerror("illegal index size for",sp->name);
	  }
	}
	offset = 0;
	index = darr+MAXDIM-narg;		/* pointer to first index */
	for (i=1; i<ndim; i++) {
	  offset += *index++;			/* add the index */
	  offset *= apnt->dim[i];		/* mult by next lower dim */
	}
	offset += *index;			/* add the final index */
	switch (sp->vtype) {
	   case NUMBER:
	   case LITCHAR:
  	       parr += offset;				/* calculate array pointer */
	       break;
	   case STRING:
  	       parr = (double *)(((char **)parr) + offset); /* calculate array pointer */
	};

	return (parr);
}

/*------------------------------------------------------*/

void xdims (void)	/* return dimensions of an array */

{
     int i,ndim,*pdim;
     pdatum p={0};
     datum d1={0},d2={0};
     double *dpnt;
     Symbol newarr;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"dims\n");
#endif

   getvar(&p);				/* get existing array from stack */
   if (*p.type==ARRAY) {
     newarr.name = (char *)"_dimarr"; 
     newarr.type = ARRAY; 
     ndim = (*p.arrp)->ndim;
     pdim = (*p.arrp)->dim;
     d2.val = ndim;			/* size of new array = ndim of arr1 */
     pushm(d2);
     dpnt=darr2(&newarr,1);		/* new array has 1 dim */
     for (i=0; i<ndim; i++) {
        *(dpnt+i) = *(pdim+i);		/* copy values into array */
     }
     d1.type = ARRAY; 
     d1.vtype = NUMBER;
     d1.arrp = newarr.arrp;
     pushm(d1);
   }
}

/*------------------------------------------------------*/

Symbol *setarrayconst(int nvals)

/* Get array values from stack and put into array */

{
   int i, rangefl;
   datum d1,d2;
   static Symbol newarr={0};
   ncarray *apnt;

  newarr.name = (char *)"setarrayconst";
  newarr.type = ARRAYCONST;
  if (nvals < 0) {		/* if range initializer specified */
      double low,incr,high;
    nvals = -nvals;
    rangefl = nvals;
    getrange (&nvals, &low, &high, &incr); /* get range values */
  }
  else rangefl = 0;
  d2.val = nvals;		/* size of new array (1 dim) */
  pushm(d2);
  darr2(&newarr,1);		/* new array has 1 dim */
  apnt = newarr.arrp;
  if (rangefl) {
	double val, low, high, incr;
        int sinit;
	datum dn;

        nvals = rangefl;
	sinit = nvals;
	getrange (&nvals, &low, &high, &incr); /* get range values */
	if (sinit==2) { dn=popm(); dn=popm(); }
	else          { dn=popm(); dn=popm(); dn=popm(); }
	for (val=low,i=0; i<nvals; i++, val+=incr) {	/* initialize array */
           apnt->arr[i] = val; 	/* set array val */
	}
  }
  else {
   for (i=nvals-1; i>=0; i--) {	/* initialize array */
    d1=popm();
    newarr.vtype = d1.vtype; /* set array type (number, etc.) */
    switch (d1.vtype) {
        case NUMBER:
        case LITCHAR:
          apnt->arr[i] = d1.val; 	/* set array val */
          break;
        case STRING:
          ((char **)apnt->arr)[i] = (char *)d1.str; 	/* set array val */
        break;
     }
  }
 }
 d2.type = VAR;
 d2.val = 0;
 //pushm (d2);		/* number of indexes (for garray()) */
 return (&newarr);
}

/*------------------------------------------------------*/

void initarr(void)	/* initialize array */

{
	Symbol *sp;
	datum d1,*d;
	int arg,i,ndim,narg,ninit,size,rangefl;
        ncarray *apnt;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"initarr\n");
#endif
	sp = (Symbol *)(*pc++);			/* pointer to array symbol */
	ninit = (long int) *pc++;		/* number of initialize vals */
	if (ninit<0) {				/* range { expr : expr : expr } */
	  ninit = -ninit;
	  rangefl = 1;
	}
	else rangefl = 0;
	narg = (long int) sp;
	if (narg==0) {
	   d = &stack[stackp-1-ninit];		/* get local arr from stack */
	   apnt = d->arrp;
	   arg = 1;
	}
        else if (narg > 0 && narg < 1000) {  /* if var is arg to a function */
             d = &stack[fp->argn + narg - fp->nargs];
	     apnt = d->arrp;
	     arg = 1;
	}
        else {
	     apnt = sp->arrp;
	     arg = 0;
	}
        ndim = apnt->ndim;                      /* get # of dimensions */

        for (size=1,i=0; i<ndim; i++) {		/* calculate array size */
          size *= apnt->dim[i];
	}
	if (rangefl) {				/* read low, incr, high */
	     double low,incr,high,val;		/*  for range */
             int lsize,sinit;
	     datum dn;
	
	  if (arg) d->vtype = NUMBER;		/* set array type */
	  else    sp->vtype = NUMBER; 
	  sinit = ninit;
	  getrange (&ninit, &low, &high, &incr); /* get range values */
	  if (sinit==2) { dn=popm(); dn=popm(); }
	  else          { dn=popm(); dn=popm(); dn=popm(); }
	  lsize = size;
	  if (ninit < lsize) lsize = ninit;
	  for (val=low,i=0; i<lsize; i++, val+=incr) {	/* initialize array */
	     apnt->arr[i] = val; 			/* set array val */
	  }
	} else {
	 for (i=ninit-1; i>=0; i--) {	/* initialize array */
	   d1=popm();
	   if (arg) d->vtype = d1.vtype;	/* set array type (number, etc.) */
	   else    sp->vtype = d1.vtype; /* set array type (number, etc.) */
	   if (i<size) 
	    switch (d1.vtype) {
	     case NUMBER:
	     case LITCHAR:
	       apnt->arr[i] = d1.val; 	/* set array val */
	       break;
	     case STRING:
	       ((char **)apnt->arr)[i] = (char *)d1.str; 	/* set array val */
	       break;
	    }
	 }
	}
	if (ninit==1) {		/* if only 1 value, set whole array to it */
	  for (i=1; i<size; i++) {
	   switch (d1.vtype) {
	    case NUMBER:
	    case LITCHAR:
	      apnt->arr[i] = d1.val; 	/* set all of array to 'val' */
	      break;
	    case STRING:
	      ((char **)apnt->arr)[i] = (char *)d1.str; 	/* set array val */
	       ninit = size;		/* don't erase array below */
	      break;
	   }
	  }
	}
	if (d1.vtype==STRING)
	  for (i=size-1; i>=ninit; i--) {	/* initialize string array */
              ((char **)apnt->arr)[i] = (char*)NULL;  /* init array val */
	  }
}

/*------------------------------------------------------*/

void dofft(void) /* fourier transform on an array */
		 /* array must be:   arr[length][2]
		    where arr is the array name;
		    and length is the length of the array.
		    arr[0] is transformed into the amplitude fft, and
	 	    arr[1] is transformed into the phase fft.
		*/
{
	int size1;
	double *realpt, *imagpt;
	Symbol *sp,*param;
	ncarray *arrp;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"dofft\n");
#endif
	param = (Symbol *)(*pc++);		/* pointer to fft type symbol */
	sp = (Symbol *)(*pc++);			/* pointer to array symbol */
	arrp = sp->arrp;		/* get pointer to address of array */
	if (!arrp) {
	   execerror("fft: array space not allocated for",sp->name);
	   return;
	}
	realpt = *(double **)arrp;
	size1 = arrp->dim[1];				/* get second dim */
	imagpt = realpt + size1;
	fft (realpt,imagpt,size1,param->type);
}

/*------------------------------------------------------*/

void do_lmfit(void) 
	
/* Levenberg-Marquardt curve fitting on an array. */
/* Adapted from lmfit-2.4 from sourceforge.net */

/* lmfit(func,data[2],coeff[1])
  data array must be: arr[2][length], arr[0] is x, and arr[1] is y.
  coeff array must be: coeff[], which is starting value of coefficients.  
*/

{
	int m_dat, n_p;
	double *datax, *datay, *coeff;
	Symbol *sp,*func;
	ncarray *datap,*coeffp;
	lm_control_type control;
	lm_data_type data;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"do_lmfit\n");
#endif
   func = (Symbol *)(*pc++);		/* pointer to nc function to be fit */
   sp = (Symbol *)(*pc++);		/* pointer to array symbol */
   datap = sp->arrp;			/* get pointer to address of array */
   if (!datap) {
      execerror("do_lmfit: array space not allocated for",sp->name);
      return;
   }
   sp = (Symbol *)(*pc++);		/* pointer to array symbol */
   coeffp = sp->arrp;			/* get pointer to address of array */
					/* coefficients can be anything but */
   					/*  {0,0,...0} */
   if (!coeffp) {
      execerror("do_lmfit: array space not allocated for",sp->name);
      return;
   }
   datax = datap->arr;
   m_dat = datap->dim[1];		/* get second dim */
   datay = datax + m_dat;

   //fprintf (stderr,"%g %g %g %g\n",datay[0],datay[1],datay[2],datay[3]);

   coeff = coeffp->arr;		/* get coefficients */
   n_p = coeffp->dim[0];		/* number of coeff */

   lm_initialize_control(&control);

   data.user_funcp = func;		/* the user-defined fit function */
   data.user_t = datax;
   data.user_y = datay;
   data.user_n_p = n_p;

   					/* perform the fit */
   //fprintf (stderr,"A %d %d\n",n_p,m_dat);
   lm_minimize (m_dat,n_p,coeff,NULL,lm_evaluate_interp,lm_print_interp,&data,&control);
   if (info>=1) ncfprintf(stdout,"status: %s after %d evaluations\n",
		                 lm_shortmsg[control.info], control.nfev);
}

/*------------------------------------------------------*/

void do_lmfit2d(void) 
	
/* Levenberg-Marquardt curve fitting on an array. */
/* Adapted from lmfit-2.4 from sourceforge.net */

/* lmfit(func,data[2],coeff[1])
  data array must be: arr[2][length], arr[0] is x, then y, and arr[1] is z.
  coeff array must be: coeff[], which is starting value of coefficients.  
*/

{
	int m_dat, n_p;
	double *datax, *datay, *coeff;
	Symbol *sp,*func;
	ncarray *datap,*coeffp;
	lm_control_type control;
	lm_data_type data;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"do_lmfit\n");
#endif
   func = (Symbol *)(*pc++);		/* pointer to nc function to be fit */
   sp = (Symbol *)(*pc++);		/* pointer to array symbol */
   datap = sp->arrp;			/* get pointer to address of array */
   if (!datap) {
      execerror("do_lmfit: array space not allocated for",sp->name);
      return;
   }
   sp = (Symbol *)(*pc++);		/* pointer to array symbol */
   coeffp = sp->arrp;			/* get pointer to address of array */
					/* coefficients can be anything but */
   					/*  {0,0,...0} */
   if (!coeffp) {
      execerror("do_lmfit: array space not allocated for",sp->name);
      return;
   }
   datax = datap->arr;
   m_dat = datap->dim[1];		/* get second dim */
   datay = datax + m_dat;

   //fprintf (stderr,"%g %g %g %g\n",datay[0],datay[1],datay[2],datay[3]);

   coeff = coeffp->arr;		/* get coefficients */
   n_p = coeffp->dim[0];		/* number of coeff */

   lm_initialize_control(&control);

   data.user_funcp = func;		/* the user-defined fit function */
   data.user_t = datax;
   data.user_y = datay;
   data.user_n_p = n_p;

   					/* perform the fit */
   //fprintf (stderr,"A %d %d\n",n_p,m_dat);
   lm_minimize (m_dat,n_p,coeff,NULL,lm_evaluate_2d_interp,lm_print_2d_interp,&data,&control);
   if (info>=1) ncfprintf(stdout,"status: %s after %d evaluations\n",
		                 lm_shortmsg[control.info], control.nfev);
}

/*------------------------------------------------------*/

void erarr (ncarray *arrp)

{
  if (arrp) {
    arrp->ndim=0;
    efree((char *)((double *)arrp->arr));	/* erase array */
    efree((char *)arrp);		/* erase array descriptor */
  }
}

/*------------------------------------------------------*/

void erasarr(void)			/* erase array */

{
	Symbol *sp, *arrsym;
	
#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"erasarr\n");
#endif
	sp = (Symbol *)(*pc++);
	rmsym (sp);
}

/*------------------------------------------------------*/

double *retvar (char *arrname, int arg1, int arg2, int arg3, int narg)

/* Allow access to variable or array from outside simulator.
    Given name and indices of array, return pointer to element in array.
    If number of indices is 0, assume it is a normal variable.
*/

{
   Symbol *arrsym; 
   double *arrpnt;
   datum d0={0},d1={0},d2={0},d3={0};

  if (arrname) arrsym = lookup(arrname);
  if (!arrname || !arrsym) {
         ncfprintf (stderr,"retarr: invalid name for array\n");
         return ((double*)NULL);
  }
  switch (narg) {

    case 0:
          if (arrsym->type==VAR)
	      arrpnt = &arrsym->val;
	  else arrpnt = (double*)NULL;
	  break;
    case 1:
	  d1.val = arg1;
	  d0.sym = (Symbol *)(long int)narg;
	  pushm (d1);
	  pushm (d0);
	  arrpnt = garrayval(arrsym);
	  break;
    case 2:
	  d1.val = arg1;
	  d2.val = arg2;
	  d0.sym = (Symbol *)(long int)narg;
	  pushm (d1);
	  pushm (d2);
	  pushm (d0);
	  arrpnt = garrayval(arrsym);
	  break;
    case 3:
	  d1.val = arg1;
	  d2.val = arg2;
	  d3.val = arg3;
	  d0.sym = (Symbol *)(long int)narg;
	  pushm (d1);
	  pushm (d2);
	  pushm (d3);
	  pushm (d0);
	  arrpnt = garrayval(arrsym);
	  break;
    default:
	  arrpnt = (double*)NULL;
	  break;
  }
 return (arrpnt);
}

/*------------------------------------------------------*/

void varread(void)	/* read into variable */
{
	datum d1={0},d2={0};
	pdatum p={0};
	char vbuf[120];
	int len;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"varread\n");
#endif
	d2 = stack[stackp-1];
	getvar(&p);
	switch (*p.vtype) {
	  case 0:
		switch (fscanf(stdin, "%lf", p.val)) {
		case EOF:
			d1.val = *p.val = 0.0;
			break;
		case 0:
			vbuf[0]=0;
			switch (fscanf(stdin, "%s",vbuf)) {
			case EOF:
				d1.val = *p.val = 0;
				break;
			case 0:
				execerror("non-string read into", d2.sym->name);
				break;
			default:
				d1.val = 1.0;
				if ((len=strlen(vbuf))>0) {
				  *p.str = emalloc (len+1);
				  strcpy(*p.str,vbuf);
				}
				break;
			}
        		*p.vtype = STRING;
			break;
		default:
			d1.val = 1.0;
        		*p.vtype = NUMBER;
			break;
		}
		break;
	  case NUMBER:
	  default:
		switch (fscanf(stdin, "%lf", p.val)) {
		case EOF:
			d1.val = *p.val = 0.0;
			break;
		case 0:
			execerror("non-number read into", d2.sym->name);
			break;
		default:
			d1.val = 1.0;
			break;
		}
        	*p.vtype = NUMBER;
		break;

	  case STRING:
		vbuf[0]=0;
		switch (fscanf(stdin, "%s",vbuf)) {
		case EOF:
			d1.val = *p.val = 0;
			break;
		case 0:
			execerror("non-string read into", d2.sym->name);
			break;
		default:
			d1.val = 1.0;
			if ((len=strlen(vbuf))>0) {
			  *p.str = emalloc (len+1);
			  strcpy(*p.str,vbuf);
			}
			break;
		}
        	*p.vtype = STRING;
		break;
	}
	/* *p.type = VAR; */
	d1.vtype = NUMBER;
	d1.type = VAR;
	pushm(d1);
}


/*------------------------------------------------------*/

#define TBUFSIZ 1024
// #define TBUFSIZ 32

struct pa { 
    int    size;
    int    windx;
    int    rindx;
    struct pa *next; 
    double arr[TBUFSIZ]; 
 };

/* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

pa *makpa (pa *last)

{ 
   pa *bpnt;

   if (!(bpnt=(pa*)emalloc(sizeof(pa)))) {
     ncfprintf(stderr,"nc: makpa, can't allocate.\n");
   }
   bpnt->next = (pa *)NULL;
   bpnt->size = TBUFSIZ;
   bpnt->rindx = 0;
   bpnt->windx = 0;
   if (last) last->next = bpnt;
   return (bpnt); 
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

void delpa (pa *bpnt)

{
  if (!bpnt) return;
  delpa (bpnt->next); 
  efree (bpnt);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

pa *storpa(pa *bpnt, double val)

/* Store a numeric value in a pa buffer. */
/*  Allocate a new buffer if the old one is full. */
/* The first time, bpnt is set to NULL, and a new buffer */
/*  is created. On succeeding calls, a new value is stored */
/*  and indx is incremented. */

{
     int makfl;

   if      (!bpnt) 		       makfl = 2;  /* the first time */
   else if (bpnt->windx >= bpnt->size) makfl = 1;  /* overflow */
   else 			       makfl = 0;  /* normal storage */
   if (makfl) {
      bpnt = makpa(bpnt);
   }
   if (makfl<2) bpnt->arr[bpnt->windx++] = val;
   return bpnt;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

pa *storpas(pa *bpnt, char *val)

/* Store a string value in a pa buffer. */
/*  Allocate a new buffer if the old one is full. */
/* The first time, bpnt is set to NULL, and a new buffer */
/*  is created. On succeeding calls, a new value is stored */
/*  and indx is incremented. */

{
     int makfl;

   if      (!bpnt) 		       makfl = 2;  /* the first time */
   else if (bpnt->windx >= bpnt->size) makfl = 1;  /* overflow */
   else 			       makfl = 0;  /* normal storage */
   if (makfl) {
      bpnt = makpa(bpnt);
   }
   if (makfl<2) ((char**)(bpnt->arr))[bpnt->windx++] = val;
   return bpnt;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

double getpa(pa **bpnt)

{
   pa *b;

   b = *bpnt;
   if (!b) {
	 return 0;
   }
   if (b->rindx >= b->size) {
      if (!(b=b->next)) {
	 *bpnt = NULL;
	 return 0;
      }
   } 
   *bpnt = b;
   return (b->arr[b->rindx++]);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

char *getpas(pa **bpnt)

{
   pa *b;
   char *val;

   b = *bpnt;
   if (!b) {
	 return 0;
   }
   if (b->rindx >= b->size) {
      if (!(b=b->next)) {
	 *bpnt = NULL;
	 return 0;
      }
   } 
   *bpnt = b;
   return (((char **)(b->arr))[b->rindx++]);
}

/*------------------------------------------------------*/

#define STRSIZ 15000

pa *do_fread(const char *filnam, int *nlongp, int *nwidp, int vtyp)

{
    int j,nlong,nwid,strl,linenum;
    double val;
    FILE *ftemp;
    char *cp,*np,*strp;
    pa *bst;
    pa *bpnt;
    Symbol *s;
    char *cbuf;
    char *str, *strtmp;
    char xbuf[150];
    static char ebuf[100];
    size_t strsiz = STRSIZ;
    static const char delim[] = {" ,:\t\n\r"};

 if (strcmp(filnam,"stdin")==0) ftemp = stdin;
 else if ((ftemp=fopen(filnam,"r")) == NULL) {
    	 // sprintf (ebuf,"fread: can't open file '%.50s'.\n",filnam);
         // execerror ("warning,",ebuf);
   return NULL;
 }
 //str = emalloc(STRSIZ);			/* make buffer to hold a line */
 str = NULL;
 strsiz = 0;
 for (linenum=1; (((strl=getline(&str,&strsiz,ftemp))>0) && *str== '#'); linenum++); /* get first line */
 if (strl<=0) {
    fclose (ftemp);
    sprintf (xbuf,"fread: file '%.50s' is zero length.\n",filnam);
    execerror ("warning,",xbuf);
    return NULL;
 }
 cp=strstr(str,"#");			/* look for # after values */
 if (cp) {
    *cp = '\0';				/* terminate string at # */
 }
				/* must use buffered read for stdin */

 bst = storpa((pa*)NULL,0);		/* set up temporary buffer */
 bpnt = bst;

 cbuf = emalloc(STRSIZ); 
 strncpy(cbuf,str,STRSIZ-1);

 /* Since we're scanning the first line twice, strtok() needs cbuf to be
    a different address than str, because it remembers the address of 
    the array it is given. Otherwise cbuf would not be needed.
  */

 /* Count how many text tokens on line. */
 /* This can be numbers, variables, or expressions without spaces */

 for (cp=np=cbuf,nwid=0; cp; nwid++) {	/* how many tokens on line */
    cp = strtok(np,delim);		/* get pointer to next token */
    np = (char *)NULL;
 }
 nwid--;
 efree(cbuf);

 for (cp=np=str,j=0; j<nwid; j++) {			/* scan first line */
   cp = strtok(np,delim);
   if (cp==NULL) {
      sprintf (xbuf,"Error in format of file '%.50s' while reading array.\n",filnam);
      execerror ("warning,",xbuf);
   }
   np = (char *)NULL;	

   if (vtyp==NUMBER) {
    if (fread_expr) {		/* allows file read by fread to contain expressions (no spaces) */
      val = fread_interp(cp); 
      bpnt = storpa(bpnt,val);
    } else {
      if (isalpha(*cp) || (*cp=='_')) {

       val = getvarval(cp);			/* look in "command line" variable list */
       if (val!=LARGENUM) bpnt = storpa(bpnt,val);
       else {
         s = lookup(cp);			/* look in orig symbol table */ 
         if (s) bpnt = storpa(bpnt,s->val);
         else {
	  sprintf (xbuf,
          "# Missing def for var '%.20s' in file '%.50s' while reading array, first line.\n", 
							cp,filnam);
           execerror ("warning,",xbuf);
           bpnt = storpa(bpnt,0);
	 }
       }
      } 
      else   {  double tval;
         tval = atof(cp);
         bpnt = storpa(bpnt,tval);
      }
    }
   } else {   /* (vtyp!=NUMBER) */
      strtmp = emalloc(strlen(cp)); 
      strcpy (strtmp,cp);
      bpnt = storpas(bpnt,strtmp);
   }
 }   /* for j */
 efree(str);
     
						/* read in lines from file */
 for (nlong=1, linenum++; ; nlong++,linenum++) {
   str = NULL;
   strsiz = 0;
   while ((strl=getline(&str,&strsiz,ftemp)>0) && *str=='#') linenum++;/* get next line */
   if (strl<=0) break;
   np = str;
   cp=strstr(str,"#");			/* look for # after values */
   if (cp) {
     *cp = '\0';				/* terminate string at # */
   }

 for (j=0; j<nwid; j++) {			/* scan lines */
   cp = strtok(np,delim);
   if (cp==NULL) {
       nlong--; 
       if (j==0) break;				/* allow empty lines */
       else { 
	  ncfprintf (stderr,"# %s: error in format of file '%.50s' while reading array, line %d.\n", 
			  	progname,filnam,linenum);
          break;
       }
      //sprintf (xbuf,"Error in format of file '%.50s' while reading array, line %d.\n", filnam,linenum);
      //execerror ("warning,",xbuf);
      //	    break;
	    /* execerror ("warning,",xbuf); */
   }
   np = (char *)NULL;	

   if (vtyp==NUMBER) {
    if (fread_expr) {	       /* allows file read by fread to contain expressions (no embedded spaces) */
     val = fread_interp(cp); 
     bpnt = storpa(bpnt,val);
    } else {
     if (isalpha(*cp) || (*cp=='_')) {
       val = getvarval(cp);		/* look in "command line" variable list */
       if (val!=LARGENUM) bpnt = storpa(bpnt,val);
       else {
         s = lookup(cp); 
         if (s) bpnt = storpa(bpnt,s->val);
         else {
	  sprintf (xbuf,
          "# Missing def for var '%.20s' in file '%.50s' while reading array, line %d.\n",cp,filnam,linenum);
           execerror ("warning,",xbuf);
           bpnt = storpa(bpnt,0);
	 }
       }
      } 
      else   {  double tval;
            tval = atof(cp);
            bpnt = storpa(bpnt,tval);
      }
    }
   } else {   /* (vtyp!=NUMBER) */
      strtmp = emalloc(strlen(cp)); 
      strcpy (strtmp,cp);
      bpnt = storpas(bpnt,strtmp);
   }
 }   /* for j */
 efree(str);
     
 }  /* for (nlong=1;;) */
 *nlongp = nlong;
 *nwidp  = nwid;
  fclose (ftemp);
  return bst;
}

#undef STRSIZ

/*------------------------------------------------------*/

void xfread(void)

/* define and read array from file */

/* Can be 1 or 2 dimensional; 
   Number of lines defines first dimension;
   Number of elements on first line (first row) defines second dimension;
   Variables as well as numbers are allowed;
   but these variables must already exist in symbol table;
   Comment lines start with "#".
*/

{
    const char *filnam;
    Symbol *var,*s;
    int i,j,narg,nwid,nlong;
    datum d1={0},d2={0};
    pdatum p1={0},p2={0};
    char *cp,*np,*strp;
    double *dpnt;
    char xbuf[150];
    pa *bst;
    pa *bpnt;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"xfread\n");
#endif
    var = (Symbol *)(*pc++);		/* array to be created */
    narg = (long int) *pc++;		/* number of var's (dimensions) */
    if (narg==2) getvar(&p2);		/* get pointer to second dim */
    getvar(&p1);			/* get pointer to first dim */

    if (!checkstr(d1=popm())) {
	return; 
    }
    filnam = d1.str;		/* file to be read */
    bst = do_fread(filnam,&nlong,&nwid,NUMBER);
    if (bst==NULL) {
       sprintf (xbuf,"# Can't read file '%.50s'\n", filnam);
       execerror ("stopping,",xbuf);
    }
    d1.val = nlong;
    d2.val = nwid;
    d1.type = d2.type = CONST;
    d1.vtype = d2.vtype = NUMBER;

    if (nwid > 1 && nlong > 1) {
	pushm(d1);
	pushm(d2);
	dpnt=darr2(var,2);
    }
    else {
	if (nlong > 1) pushm(d1);
	else 	       pushm(d2);
	dpnt=darr2(var,1);
    } 
    if (dpnt==NULL) {
       sprintf (xbuf,"# Can't allocate array '%.20s' in file '%.50s'\n", var->name,filnam);
       execerror ("stopping,",xbuf);
    } 
    else {
     bpnt = bst;				/* reset buffer pointer for read */
     for (i=0; i<nlong; i++) {
        int indx;
       indx = i*nwid; 
       for (j=0; j<nwid; j++) {
          *(dpnt+indx+j) = getpa(&bpnt);	/* copy values into array */
       }
     } 
    delpa(bst);				/* delete temporary buffer */
    *p1.val = nlong;
    *p1.vtype = var->vtype = NUMBER;
    if(narg==2) {
	 *p2.val = nwid;
         *p2.vtype = NUMBER;
    }
  }
}

/*------------------------------------------------------*/

void xfreads(void)

/* define and read string array from file */

/* Can be 1 or 2 dimensional; 
   Number of lines defines first dimension;
   Number of elements on first line (first row) defines second dimension;
   Reads strings only, no variable interpretation
   Comment lines start with "#".
*/

{
    const char *filnam;
    Symbol *var,*s;
    int i,j,narg,nwid,nlong;
    datum d1={0},d2={0};
    pdatum p1={0},p2={0};
    char *cp,*np,*strp;
    char **dpnt;
    char xbuf[150];
    pa *bst;
    pa *bpnt;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"xfreads\n");
#endif
    var = (Symbol *)(*pc++);		/* array to be created */
    narg = (long int) *pc++;		/* number of var's (dimensions) */
    if (narg==2) getvar(&p2);		/* get pointer to second dim */
    getvar(&p1);			/* get pointer to first dim */

    if (!checkstr(d1=popm())) {
	return; 
    }
    filnam = d1.str;		/* file to be read */
    bst = do_fread(filnam,&nlong,&nwid,STRING);
    if (bst==NULL) {
       sprintf (xbuf,"# Can't read file '%.50s'\n", filnam);
       execerror ("stopping,",xbuf);
    }
    d1.val = nlong;
    d2.val = nwid;
    d1.type = d2.type = CONST;
    d1.vtype = d2.vtype = NUMBER;
    var->vtype = STRING;

    if (nwid > 1 && nlong > 1) {
	pushm(d1);
	pushm(d2);
	dpnt=(char **)darr2(var,2);
    }
    else {
	if (nlong > 1) pushm(d1);
	else 	       pushm(d2);
	dpnt=(char **)darr2(var,1);
    } 
    if (dpnt==NULL) {
       sprintf (xbuf,"# Can't allocate array '%.20s' in file '%.50s'\n", var->name,filnam);
       execerror ("stopping,",xbuf);
    } 
    else {
     bpnt = bst;				/* reset buffer pointer for read */
     for (i=0; i<nlong; i++) {
        int indx;
       indx = i*nwid; 
       for (j=0; j<nwid; j++) {
          *(dpnt+indx+j) = getpas(&bpnt);	/* copy values into array */
       }
     } 
    delpa(bst);					/* delete temporary buffer */
    *p1.val = nlong;
    *p1.vtype = NUMBER;
    if(narg==2) {
	 *p2.val = nwid;
         *p2.vtype = NUMBER;
    }
  }
}

/*------------------------------------------------------*/

double *fread (const char *filnam, int *nlongp, int *nwidp)

{
  int i, j, nlong, nwid;
  pa *bpnt,*bst;
  double *dpnt;

  if ((bst=do_fread(filnam,nlongp,nwidp,NUMBER))==NULL) {
     return NULL;
  }
  nlong = *nlongp;
  nwid  = *nwidp;
  dpnt = (double *)emalloc (nlong*nwid*sizeof(double));
  bpnt = bst;                                /* reset buffer pointer for read */
  for (i=0; i<nlong; i++) {
     int indx;
    indx = i*nwid;
    for (j=0; j<nwid; j++) {
       *(dpnt+indx+j) = getpa(&bpnt);        /* copy values into array */
    }
  }
 delpa(bst);                         /* delete temporary buffer */
 return (dpnt);
}

/*------------------------------------------------------*/

double *freads (const char *filnam, int *nlongp, int *nwidp)

{
  int i, j, nlong, nwid;
  pa *bpnt,*bst;
  double *dpnt;

  if ((bst=do_fread(filnam,nlongp,nwidp,STRING))==NULL) {
     return NULL;
  }
  nlong = *nlongp;
  nwid  = *nwidp;
  dpnt = (double *)emalloc (nlong*nwid*sizeof(double));
  bpnt = bst;                                /* reset buffer pointer for read */
  for (i=0; i<nlong; i++) {
     int indx;
    indx = i*nwid;
    for (j=0; j<nwid; j++) {
       *(dpnt+indx+j) = getpa(&bpnt);        /* copy values into array */
    }
  }
 delpa(bst);                         /* delete temporary buffer */
 return (dpnt);
}

/*------------------------------------------------------*/

void xfwrite(void)

/* write array to file.
 
   Array can be 1 or 2 dimensional; 
   Size of first dimension defines number of lines;
   Size of second dimension defines number of columns;
*/

{
    const char *filnam;
    FILE *ftemp;
    Symbol *var;
    double *vpnt;
    int i,j;
    int dims,dim1,dim2,**apnt,*dpnt,vtyp;
    datum d;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"xfwrite\n");
#endif
  if (!checkstr(d=popm())) {
	 pc+=2;
         return;
  }
  filnam = d.str;			/* file to be written */
  var = (Symbol *)(*pc++);		/* array */
  apnt = (int **)var->arrp;		/* pointer to pointer to array */
  dpnt = (int *)(apnt + 1);		/* pointer to number of dims */
  dims = *dpnt;				/* number of dimensions */
  dim1 = *(dpnt+1);			/* size of dimensions */
  dim2 = *(dpnt+2);
  vpnt = (double *)*apnt;
  vtyp = var->vtype;

  if (strcmp(filnam,"stdout")==0) 
	ftemp = stdout;
  else if ((ftemp=fopen(filnam,"w")) == NULL) {
    ncfprintf (stderr,"fwrite: can't open file %s\n",filnam);
    return;
  }

  switch (dims) {

  case 2:
   for (i=0; i<dim1; i++) {
     for (j=0; j<dim2; j++,vpnt++) {
       if (j==(dim2-1)) {
	  if (vtyp==NUMBER) ncfprintf (ftemp,"%g",*vpnt);
	  else              ncfprintf (ftemp,"%s",(char*)(long)*vpnt);
       } else {
	  if (vtyp==NUMBER) ncfprintf (ftemp,"%g ",*vpnt);
	  else              ncfprintf (ftemp,"%s ",(char*)(long)*vpnt);
       }
     }
     ncfprintf (ftemp,"\n");
   }
   break;

   case 1:
     for (i=0; i<dim1; i++) {
	  if (vtyp==NUMBER) ncfprintf (ftemp,"%g\n",*vpnt++);
	  else              ncfprintf (ftemp,"%s",(char*)(long)*vpnt++);
     }
    break;

    default:
    break;

  }  /* switch dims */

  fclose(ftemp);
}

/*------------------------------------------------------*/

void xunlink(void)

/* erase file.  */

{
    const char *filnam;
    datum d;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"xunlink\n");
#endif
  if (!checkstr(d=popm())) {
	 pc+=2;
         return;
  }
  filnam = d.str;			/* file to be written */

  if (unlink(filnam) != 0) {
    ncfprintf (stderr,"unlink: can't remove file %s\n",filnam);
    return;
  }
}

/*------------------------------------------------------*/

void getflds(void)

/* Read a line from a stream and put text tokens into array */

{
      int i, j, nwid,strl;
      char *cp,*np;
      static char *str=NULL;
      static datum d1,d2;
      FILE *fd;
      size_t strsiz = 0;
      static const char delim[] = {" ,:\t\n\r"};
      static char **arr = NULL;
      static ncarray *arrp = NULL;

   nwid = 1000;					/* set max fields = 100) */
   if (arrp!=NULL) {				/* erase contents of previous Fld array */
       if (arrp->arr!=NULL) {
           for (i=0; i<arrp->dim[0]; i++) {
               if (arr[i]!=NULL) efree(arr[i]);
           }
       }
   }
   else { 					/* create new array */
       arrp = (ncarray *)emalloc(sizeof(ncarray));
       arr = (char **)emalloc(nwid*sizeof(char *)); 
       arrp->arr = (double *)arr;
   }
   arrp->ndim = 1;
   arrp->dim[0] = 0;

   if (!checknum(d1=popm())) {                  /* get filedesc off stack */
         d2.type=UNDEF;
         push (d2);
         return;
   }
   fd = (FILE *)d1.str;
   str = NULL;					/* getline will make buffer to hold line */
   strsiz = 0;
   if ((strl=getline(&str,&strsiz,fd))<0) {	/* get line of text */
       fclose (fd);
       arrp->dim[0] = 0;
       d2.type = ARRAY;
       d2.vtype = STRING;
       d2.arrp = arrp;
       push (d2);
       if (str!=NULL) efree(str); 
       return;
   }
   np = str;
   arr[0] = (char *)emalloc(strlen(str)+1); 
   strncpy (arr[0], str, strlen(str)-1); 	/* don't copy \n at end of line */
   for (cp=np=str,j=1; j<nwid; j++) {  
       if ((cp=strtok(np,delim))==NULL) break;
       arr[j] = (char *)emalloc(strlen(cp)+1); 
       strcpy (arr[j], cp); 
       np = NULL;
   }
   arrp->dim[0] = j;
   d2.type = ARRAY;
   d2.vtype = STRING;
   d2.arrp = arrp; 
   push(d2);
   if (str!=NULL) efree(str); 
}

/*------------------------------------------------------*/

void gplot(void)

/* low-level plotting commands */

{
   int narg,fill;
   Symbol *param;
   datum d1,d2;
   double x,y;
   double x1,y1,x2,y2,x3,y3,x4,y4;
   const char *str=NULL;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"gplot\n");
#endif
 param  = (Symbol *)*pc++;
 narg = (long int) *pc++;
 if (param->type==GRECT) {
     d1 = popm(); fill = (int)(d1.val);
     d1 = popm(); y4 = d1.val;
     d1 = popm(); x4 = d1.val;
     d1 = popm(); y3 = d1.val;
     d1 = popm(); x3 = d1.val;
     d1 = popm(); y2 = d1.val;
     d1 = popm(); x2 = d1.val;
     d1 = popm(); y1 = d1.val;
     d1 = popm(); x1 = d1.val;
  }
 else if (param->type==GWINDOW) {
     d1 = popm(); y2 = d1.val;
     d1 = popm(); x2 = d1.val;
     d1 = popm(); y1 = d1.val;
     d1 = popm(); x1 = d1.val;
  }
 else {
   if (narg > 1) {
     d2 = popm();
     y = d2.val;
   }
   else y = 0;
   if (narg>0) d1 = popm();
   switch (d1.vtype) {
     default:
     case 0:
     case NUMBER:
      x   = d1.val; break;
     case STRING:
      str = d1.str; break;
   }
 }
#ifndef XSTIMM
 switch (param->type) {


   case GMOVE:  gmove (x,y);
		break;

   case GDRAW:  gdraw (x,y);
		break;

   case GRMOVE: grmove (x,y);
		break;

   case GRDRAW: grdraw (x,y);
		break;

   case GPEN:   gpen ((int) x);
		break;

#ifndef XSTIM
   case GVPEN:  plotpen (int(x), int(y-1));
		break;
#endif

   case GROT:   grotate (x);
		break;

   case GCROT:  gcrotate (x);
		break;

   case GORIG:  gorigin (x,y);
		break;

   case GFRAME:	if (str!=NULL) gframe (str);
		break;

   case GCWID:	gcwidth (x);
		break;

   case GSIZ:	gsize (x);
		break;

   case GDASH:	gdash ((int)x);
		break;

   case GCIRC:	gcirc (x,(int)y);		/* radius, fill */
		break;

   case GWINDOW:gwindow (x1,x2,y1,y2);
		break;

   case GRECT:	grect (x1,y1,x2,y2,x3,y3,x4,y4,fill);
		break;

   case GPURGE:	gpurge ();
		break;
	}
#endif
}

/*------------------------------------------------------*/

#define GRPH

#ifdef GRPH

#include "ncplot.h"

extern plotfr plotnod[PLOTNODSIZ];  /* holds nodes to be recorded */
extern int plotnum;            /* number of plots to be displayed */
extern int vidmode;		/* =1 -> graphics, =0 -> text */

void mplot (double y, double x, int n, int i, int pflag);
void varcopy();
void varcopyu();
void plotinit(int plotnum);
void plotrst(int plotnum);
void plotpen (int val, int i);
void plotchar (int val,int lines, int i);
void plotcsiz (double val, int plotnum);
void plotfilt (int nfilt, double *timec, int plotnum);
void plotvpen  (Symbol *vis, int plotnum);
void plotname (const char *plname, int plotnum);
void plotn (int pl, int plotnum);
void plotsize (double plotsiz, int plotnum);
void plotval (double plotvl, int plotnum);

/*------------------------------------------------------*/

void graphparams (int type, int narg, int plotn);

void grph(void)

/* make a graph */
{
  datum d1,d2;
  Symbol *param;
  int pl;
  static int ch,narg,i,n,npops;
  static double xval,yval[PLOTNODSIZ];
  double getval(const char *str);

 param  = (Symbol *)*pc++;
 narg = (long int) *pc++;			/* (x,y1,y2,...,yn) */
 if (param == 0) {				/* get arguments for graph */
	n = narg - 1;
	if (n >= PLOTNODSIZ) n=PLOTNODSIZ-1;
/*	if (n > plotnum) n = plotnum; /* */
	for (i=n-1; i>=0; i-- ) {
	   d2 = popm();
	   yval[i] = d2.val;
	}
	checknum(d1 = popm());
	xval = d1.val;
#ifndef XSTIMM
	for (i=0; i<n; i++) {
	   mplot (yval[i],xval,n,i,1);
	}
#endif
 }
 else {
 vidmode = (int)getval("vidmode");
 switch (param->type) {
   case X:			/* bug: can't use for indiv separate graphs */
   case Y:    npops = 2;
	      break;
   case SIZE: npops = 1;
	      break;
   case INIT:
   case FILT:
   case RESTART:
	      npops=0;
	      break;
   case PEN:
   case VPEN:
   case PLNAME:
   case PLNUM:
   case PLSIZE:
   case CHAR:
   case CCHAR:
	      npops = 0;
	      break;
 }

 if (npops==2) checknum(d2 = popm());
 if (npops>=1) checknum(d1 = popm());

#ifndef XSTIMM
 switch (param->type) {
   case X:			 /* bug: can't use for indiv separate graphs */
	if (!strcmp(gfrname(),rootframe))	 /* if root fr */
	      plotnum = -1;		   /* number of plots for this graph */
	pl = plotnum+1;
	if (d1.val < d2.val) { double t; t = d1.val; d1.val=d2.val; d2.val=t; }
	plotnod[pl].pxmax  = d1.val;
	plotnod[pl].pxmin  = d2.val;
	plotnod[pl].xrange = d1.val - d2.val;
	break;

   case Y:
	plotnum++;
	if (strcmp(gfrname(),rootframe))	 /* if diff than root frame */
	strncpy(plotnod[plotnum].plframe,gfrname(),PNAMSIZ-1); /* cp frame nm */
	if (d1.val < d2.val) { double t; t = d1.val; d1.val=d2.val; d2.val=t; }
	plotnod[plotnum].pymax  = d1.val;
	plotnod[plotnum].pymin  = d2.val;
	plotnod[plotnum].yrange = d1.val - d2.val;
	plotnod[plotnum].pmod  = GREC;
	break;

  case INIT:
	if (interp) varcopy();
	plotinit(plotnum+1);
	break;

  case RESTART:
	if (interp) varcopy();
	plotrst(plotnum);
	break;
  case PEN:
  case CHAR:
  case CCHAR:
  case SIZE:
  case PLNAME:
  case PLNUM:
  case PLSIZE:
  case PLVAL:
  case FILT:
  case VPEN:
       graphparams (param->type, narg, plotnum+1);
  }
 } /* else */
}

/*------------------------------------------------------------*/

void plparams (void) 
{
    int narg;
    Symbol *param;

 param  = (Symbol *)*pc++;
 narg = (long int) *pc++;			/* (x,y1,y2,...,yn) */
 graphparams (param->type, narg, plotnum);
}

/*------------------------------------------------------------*/

void graphparams (int type, int narg, int plotnum)
{
    static int ch, i;
    datum d1;
   #define TTIMEC 10
    double ttimec[TTIMEC];

  switch (type) { 
  case PEN:
	if (narg>0) {			/* If on line by itself, */
	  for (i=narg; i>0; i-- ) {	/* set plot colors before graph init */
	   checknum(d1 = popm());	/* Can also be used w/graph*/
	   plotpen ((int)d1.val,i-1);
	  }
	}
	else {				/* If after "max min" then set color */
	   checknum(d1 = popm());	/*  for this plot only */
	   plotpen ((int)d1.val,plotnum);
	}
	break;

  case CHAR:
	if (narg>0) {			/* If on line by itself, */
	  for (i=narg; i>0; i-- ) {	/*  set plot chars before graph init */
	   ch = checkchar(d1 = popm());	/*(Can also be used w/graphing)*/
	   plotchar (ch,LINES,i-1);
	  }
	}
	else {				/* If after "max min" then set char */
	   ch = checkchar(d1 = popm());	/*  for this plot only */
	   plotchar (ch,LINES,plotnum);
	}
	break;

  case CCHAR:
	if (narg>0) {			/* If on line by itself, */
	  for (i=narg; i>0; i-- ) {	/*  set plot chars before graph init */
	   ch = checkchar(d1 = popm());	/*(Can also be used w/graphing)*/
	   plotchar (ch,NOLINES,i-1);
	  }
	}
	else {				/* If after "max min" then set char */
	   ch = checkchar(d1 = popm());	/*  for this plot only */
	   plotchar (ch,NOLINES,plotnum);
	}
	break;
  case SIZE:
        checkstr(d1=popm());
	plotcsiz (d1.val,plotnum);
	break;
  case PLNAME:
        checkstr(d1=popm());
	plotname(d1.str,plotnum);
	break;
  case PLNUM:
        checknum(d1=popm());
	plotn((int)d1.val,plotnum);
	break;
  case PLSIZE:
        checknum(d1=popm());
	plotsize(d1.val,plotnum);
	break;
  case PLVAL:
        checknum(d1=popm());
	plotval(d1.val,plotnum);
	break;
  case PLARR:
        d1=popm();
	if (d1.type==ARRAY) {
	  int maxindex;
          if (d1.arrp!=NULL) {
	     if (d1.arrp->ndim==2)
	       maxindex = d1.arrp->dim[0];
	  }
	  plotarr(d1.arrp->arr,maxindex);
	}
	break;
  case FILT:
	if (interp) varcopyu();
        if (narg >= TTIMEC) narg = TTIMEC - 1;
        for (i=0; i<TTIMEC; i++) {
            ttimec[i] = 0;
        }
        for (i=narg-1; i>=0; i--) {
            checknum(d1 = popm());
            ttimec[i] = d1.val;
        }
#ifndef XSTIM
	plotfilt (narg,ttimec,plotnum);
#endif		/* ! XSTIM */
	break;

  case VPEN: {
	  Symbol *vpen;
	vpen = (Symbol *)*(pc-1);
	if (interp) varcopyu();
	
#ifndef XSTIM
	plotvpen (vpen,plotnum);
#endif		/* ! XSTIM */
	}
	break;

  }  /* switch */
#endif		/* ! XSTIMM */

}

#endif /* GRPH */

/*------------------------------------------------------*/

Inst *code(Inst f)	/* install one instruction or operand */
	       
{
	Inst *oprogp = progp;

	if (progp >= &prog[NPROG])
		execerror("program too big", (char *)0);
	*progp++ = f;
	return oprogp;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

//Inst *code(int val)	/* install one instruction or operand */
//{
//  return code ((Inst)(long int)val);
//}

/*------------------------------------------------------*/

void execute(Inst *p)
{
	Inst func;

	for (pc = p; *pc != STOPC && !stopping && !iflag; ) {
/*	ncfprintf (stderr,"%d %d\n",(int)*pc,(int)expop);  /* debug pc */
                func = *pc++;   /* inside func, pc must point to "arg"  */
                ((Inst )(*func))();
/*              (*(*pc++))();   /* old way, doesn't always work */
				/*  because pc is incr after func call */
	}

}

void edit(void)
{
    char sysbuf[40];
    datum d1;

#ifdef CODEBUG
    if (codebug) ncfprintf (stderr,"edit\n");
#endif
        if (!checkstr(d1=popm())) return;
	sprintf(sysbuf,"vi %s", d1.str);
	system (sysbuf);
}

/*------------------------------------*/

//#define CBUFSIZ 2048
//#define CBUFSIZ 4096
//#define CBUFSIZ 8192
//#define CBUFSIZ 16384
#define CBUFSIZ 32768

void xsystem(void) 

/* Run a system command.  Return its stdout and */
/*  push it on the stack.  If it is numeric, convert */
/*  it into a number and return it as a numeric value. */

{
  int c,i;
  datum d1={0},d2={0};
  Symbol *param;
  FILE *cstdout;
  char *cbuf,*strp;
  char *comndline;

 if (!(cbuf=emalloc(CBUFSIZ))) {
    execerror ("can't allocate char buf", (char *)0);
    d2.val = 0;
    push (d2);
    return;	
 }
 param = (Symbol*)*pc++;
 if (param) {		/* string inside back quotes */
   comndline=(char *)param;
 } else {		/* get expression from stack */
   d1 = popm(); 
#ifndef XSTIMM
   if (d1.type == ARRAY) {
     int i,arraysiz;
     double *pnt;

   getarr(&d1, &pnt, &arraysiz);
   for (i=0; i<arraysiz; i++) {
    if (d1.vtype!=STRING && *(pnt+i)==LARGENUM) 
     comndline = uninit;
    else
    switch (d1.vtype) {
      case 0:
      case NUMBER:
	    break;
      case LITCHAR:
	    cbuf[i] = (char)*(pnt+i);
	    break;
      case STRING:
	     comndline = (char *)(pnt+i);
	    break;
     } 
   }  
   if (d1.vtype==LITCHAR) comndline = cbuf;

 }  /* if (type==ARRAY) */

 else
 switch (d1.vtype) {
   case 0:
   case NUMBER:
	 break;
   case LITCHAR:
	 *comndline = (char)(long int)d1.val;
	 break;
   case STRING:
	 comndline = (char *)d1.str;
	 break;
  }
#endif
 } /* get expression from stack */

#ifndef XSTIMM
 if (cstdout=popen(comndline,"r")) { 
    for (i=0; i<20; i++) cbuf[i] = 0;
    for (i=0;(c = getc(cstdout)) != EOF && i<CBUFSIZ-1; i++) {
      cbuf[i] = c; 		/* copy output of pipe into buffer */
    }
    cbuf[i++] = (char)NULL; 
    pclose(cstdout);

    strp = emalloc(i);		/* save space by shortening string */
    strncpy (strp,cbuf,i);
    d2.str = strp;
    efree(cbuf);
  }
  else d2.str = NULL;
#endif
 d2.vtype = STRING;
 d2.type = CONST; 
 push (d2);
}

/*------------------------------------*/

char *xsystem (const char *str)

{
    char c;
    int i;
    char *cbuf;
    char *strp=NULL;
    FILE *cstdout;

 if (!(cbuf=emalloc(CBUFSIZ))) {
    execerror ("can't allocate char buf", (char *)0);
    return (NULL);	
 }
#ifndef XSTIMM
 if (cstdout=popen(str,"r")) { 
   for (i=0; i<10; i++) cbuf[i] = 0;
   for (i=0;(c = getc(cstdout)) != EOF && i<CBUFSIZ-1; i++) {
     cbuf[i] = c; 		/* copy output of pipe into buffer */
   }
   if ((i>0) && ((cbuf[i-1]=='\n') || (cbuf[i-1]=='\r'))) i--; 
   cbuf[i++] = (char)NULL; 

   strp = emalloc(i);		/* save space by shortening string */
   strncpy (strp,cbuf,i);
   efree(cbuf);
   pclose(cstdout);
  }
#endif
  return (strp);
}

/*------------------------------------*/

void set_argv_array(int argc, char **argv)

/* make argv array for interpreter */

{
     double argcc;
     ncarray *argvarr;
     Symbol *s,*v;

   if ((argcc = argc - 2) < 0) {
      execerror ("can't allocate argv array", "argc invalid");
   };
   s = install ("argc",VAR,argcc);
   s->atype = GLOBAL;
   s->vtype = NUMBER;

   if ((argvarr = (ncarray *)emalloc(sizeof(ncarray))) == NULL) {
      execerror ("can't allocate argv array", (char *)0);
      return;
   }
   if ((v = install ("argv",VAR,0)) == NULL) {
      execerror ("can't allocate argv array", (char *)0);
      return;
   }
   argvarr->arr = (double *)(argv+2); 		/* remove first 2 argv values "nc -c" */
   argvarr->ndim = 1;
   argvarr->dim[0] = argcc;
   v->arrp = argvarr;
   v->type = ARRAY;
   v->atype = GLOBAL;
   v->vtype = STRING;
}
