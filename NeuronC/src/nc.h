/* nc.h */

struct Symbol;
struct ncarray;

struct datum {				/* interpreter stack type */
	short type;			/* UNDEF,CONST,VAR,ARRAY,ARRVAL etc. */
	short atype;			/* GLOB, LOCAL, ARG, etc. */
	short vtype;			/* NUMBER, LITCHAR, STRING, etc. */
	union {
		double val; 		/* numeric value on stack */
		int   argnum;		/* arg number index for stack var */
		ncarray *arrp;		/* pointer to array */
		const char *str; 	/* string value on stack */ 
		Symbol *sym;		/* symbol table entry */
	};
};

struct pdatum {				/* interpreter stack type */
	short *type;			/* UNDEF,CONST,VAR,ARRAY,ARRVAL etc. */
	short *atype;			/* GLOB, LOCAL, ARG, etc. */
	short *vtype;			/* NUMBER, LITCHAR, STRING, etc. */
	Symbol **sym;
	union {
		double *val; 		/* numeric value on stack */
		int   *argnum;		/* arg number index for stack var */
		ncarray **arrp;		/* pointer to array */
		char  **str; 		/* string value on stack */ 
	};
};

typedef void (*Inst)();
typedef datum (*Fpa)(...);
#define STOPC	((Inst) 0) 

struct 	Symbol { 			/* symbol table entry */
	char *name;
	short type;			/* UNDEF,CONST,VAR,ARRAY,ARRAYVAL etc.*/
	short atype;			/* GLOB, LOCAL, ARG, etc. */
	short vtype;			/* NUMBER, LITCHAR, STRING, etc. */
	short ntype;			/* VINT, VFLOAT, VDOUBLE, etc */
	union {
		double	val;		/* VAR */
		ncarray  *arrp;		/* pointer to array */
		Fpa 	ptr;		/* BLTIN */
		Inst 	*defn;		/* FUNCTION, PROCEDURE */
		int	argnum;		/* ARG number for formal def */
		const char *str;		/* STRING */
	};
	union {
	    int    *iptr;
	    char   *cptr;
	    char   **sptr;
	    float  *fptr;
	    double *dptr;
	};
	Symbol	*next;			/* pointer to next Symbol */
	Symbol	*last;			/* pointer to last Symbol */
};

struct 	ncarray { 			/* array type, alloc on heap */
		double	*arr;		/* pointer to values */
		int 	ndim;		/* number of dimensions */
		int 	dim[1];		/* size of first dimension */
};

typedef struct SYMTYPE {
	const char *name;
	int type;
	union {
	    double val;
	    int    *iptr;
	    char   *cptr;
	    char   **sptr;
	    float  *fptr;
	    double *dptr;
	};
}stype;

Symbol *install(const char *, int, double), *lookup(const char *), 
	*lookupnt(const char *, int), *lookupnta(const char *, int, int starg);

/* #define LARGENUM 1e38 */		/* undef for float */
#define LARGENUM 1e307			/* undef for double */

#define LARGENODE 1e38			/* undef for node vals: float */
#define LARGNOD (LARGENODE * 0.999)

#define PCOND 1                         /* print compartments (conductances) */
#define PCOMP 2                         /* print compartments (as spheres) */

#define DISP  1                         /* allow "display" statement to run */
#define DCOMP 2                         /* display compartments */
#define DCONN 4                         /* display connections */
#define DNODE 8                         /* display nodes */
#define DSTIM 16                        /* display stimulus */
#define DMOVIE 32                       /* display movie */

int yylex(void);            /* interp */
void yyerror(const char *s);      /* report compile-time error */

extern "C" {
	int system(const char *str);
}

void expop();
datum pop(void);
void push(datum d);
void popsav(void);
void pushsav(void);
void stackmove(void);
void evalvar(void), add(void), sub(void);
void notinitx(void), varnum(void), varstr(void), varchr(void); 
void mul(void), xdiv(void), mod(void);
void negate(void), power(void);
void preinc(void), predec(void), postinc(void), postdec(void);
void popone(void);
void addeq(void),subeq(void),muleq(void),diveq(void),andeq(void),oreq(void);
void pfprintf(void), pprintf(void), psprintf(void);
void pscanf(void), pfscanf(void), psscanf(void), getflds(void);
void pfopen(void), pfclose(void), pfgets(void), pfgetc(void), pfputc(void);
void crlf(void), xsystem(void), xdims(void);
void darray(void), initarr(void), erasarr(void), erelem(void), eranode(void); 
void xcable(void), xgj(void), xsynapse(void), xphotrec(void); 
void xrecparm(void), xdiode(void);
void rload(void), xresistor(void), xsphere(void), xelec(void), xvbuf(void);
void xvbufd(void), xvbufg(void), xvbufo(void);
void xnbuf(void), xnbufo(void), xnbuft(void), xnbufg(void), xnbufn(void);
void conn1(void), conn1l(void), conn2s(void), conn2sl(void);
void conn2d(void), conn2dl(void);
void xnode(), relative(), modrun(void), xstim(void), membtyp(void), xset(void);
void vplot(void), vplotm(void);
void rbatt(void), xgbatt(void), xgcap(void), rcap(void), xrecord(void);
void grph(void), gplot(void), plparams(void), gglabel(void); 
void txtf(void), xfread(void), xfreads(void), xfwrite(void), xunlink(void);
void xchan(void), xcachan(void);
void noise(void), xmod(void), eedist(void);
void e3dist(void), e2dist(void), ezdist(void), efrac(void);
void n3dist(void), n2dist(void), nzdist(void);
void xgausnn(void), incrpl(void), elimit(void), initrgen(void);
void dispnod(void), dispndn(), disprot();
void efield(void), cfield(void), nfield(void);
void conn1m(void), dofft(void), do_lmfit(void), do_lmfit2d(void), gettype(void); 
void defnonly(const char *s);
void define(Symbol *sp, int narg, Inst *p, int starg);
void erasarg(int starg);

extern Inst *progp, *progbase, prog[], *code(Inst f);
void assign(void), bltin(void), varpush(void), constpush(void);
void print(void), local(void), locend(void), dlocarray(void);
void prexpr(void), bitand_x(void), bitor_x(void), varread(void);
void gt(void), lt(void), eq(void), ge(void), le(void), ne(void);
void xand(void), orx(void), xnot(void), xxor(void);
void ifcode(void), whilecode(void), call(void), arg(), argassign(), xexit();
void funcret(void), procret(void), xsizeof(void);
void forcode(void),breakcode(void),contcode(void);
void edit(), eramod(void), pushfil(void), foreacode(void), xdebugf(void);
void savemod(void), restoremod(void);

