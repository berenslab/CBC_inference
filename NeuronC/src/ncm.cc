/* segment ncm in program nc */

/* Main routines for interpreter */

extern "C" {
#include <stdio.h>
#include <ctype.h>
#include <signal.h>
#include <setjmp.h>
#include <string.h>

//jmp_buf ncbegin;

//void exit(int code);
//char *strchr(const char *ptr, int n); /* `index()' in some systems */
}

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <typeinfo>
using namespace std;

#include "nc.h"
#include "y.tab.h"
#include "script.h"
#include "ncio.h"

/* typedef int __sighandler_t; 	/* */
/* typedef void(*)(int) __sighandler_t;    /* */

char	*progname;
int	lineno = 0;	/* line number for error reports */
int	charpos = 0;	/* char position for error reports */
int	indef=0;	/* =1 means inside proc or func definition */
int	inlocal=0;	/* =1 means inside local variable definition */
int	savlocal=0;	/* =1 means inside local variable definition */
int	argcount=0;	/* count of arguments and local variables in block */
int	arrcount=0;	/* count of local arrays in block */
int	starg=0;	/* starting argcount for local vars in { } block */
int	formal=0;	/* =1 means looking for formal parameters to func */ 
char	*runfile;       /* interactive file name */
istream *tfin;		/* temp input stream pointer */
int	numfiles = 0;
int	ncerror = 0;	/* set by warning() for plotinit() */
int	iflag=0;	/* interrupt flag =1 -> program has received ^C */
int	set_tty = 0;	/* =1 -> input file is a tty, needs echo, etc. */

extern const char *einfile;	/* interactive file name */
extern istream *fin;	/* temp input stream pointer */
extern int cppfl;

extern int vidmode;
extern int *istat;
extern Inst *pc;

#define FSTACKSIZ 10
Script	*fstack[FSTACKSIZ]={NULL};

char *emalloc(unsigned int len);

void execerror(const char *s, const char *t);
int backslash(int c);
int follow(int expect, int ifyes, int ifno);
void warning (const char *s, const char *t);
void initcode();
int yyparse();
void execute(Inst *pc);
void disperr();
void resetcode();
void savecode();
void restorecode();
void ungetcx(int ch, istream *fil);
void ungetcc(int ch, istream *fil);
int checkstr(datum d);
datum pop(void);
int isa_tty();
double callproc(Symbol *procp, int npar, double par1, double par2);
double getvarval(const char *name);

#define LINBUFSIZ 240
#define MAXLINES 10

/*-------------------------------------------------------*/

static char linebuf[MAXLINES][LINBUFSIZ] = {0};
static char *cbuf = linebuf[0];
static int lineind[MAXLINES]={0};
static int charpi=lineind[0];
static int charpo=0;
static int charind=0;
static int linenum=0;
static int linelim=0;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int putcc(char ch, istream *fil)

/* Put char to stdout if input is stdin. */

{
  if (fil== &cin) { cout << ch; return ch; }
  else return ch;
}	

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void insertc(int ch)

/* Expand the buffer by one char and insert char. */

{
    int i;

  if (++charpi >= LINBUFSIZ) charpi = LINBUFSIZ-1;
  for (i=charpi; i>charind; i--) cbuf[i] = cbuf[i-1];
  cbuf[charind++] = ch;
  if (charind > charpi) charind = charpi;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void deletec()

/* Delete char and reduce the buffer by one. */

{
    int i;

  for (i=charind; i<charpi; i++) {
    cbuf[i] = cbuf[i+1];
  }
  cbuf[i] = ' ';
  if (--charpi < 0) charpi = 0;
  if (charind > charpi) charind = charpi;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void refreshline(istream *fil)

{
    int i, ch;

  putcc('\r',fil);
  for (i=0; i<charpi; i++) putcc(' ',fil);
  putcc('\r',fil);
  for (i=0; i<charpi; i++)
    if ((ch=cbuf[i]) >= 0x20) putcc(ch,fil);
  putcc('\r',fil);
  for (i=0; i<charind; i++) {
    if ((ch=cbuf[i]) >= 0x20) putcc(ch,fil);
  };
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void eraseline(istream *fil)

{
    int i, ch;

  putcc('\r',fil);
  for (i=0; i<charpi; i++) {
     putcc(' ',fil);
  };
  putcc('\r',fil);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int setlinebuf(int num)

{
    static int oldlinenum = 0;

  lineind[oldlinenum] = charpi;		/* save old index */
  if (num >= MAXLINES) num = 0;		/* make circular */
  if (num < 0) num = MAXLINES-1;
  cbuf = linebuf[num];
  charpi = lineind[num];
  if (charind>charpi) charind = charpi;
  charpo = 0;
  oldlinenum = num;
  return (num);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int getcc (istream *fil)

/* A simple line editor to handle cursor motions. */
/*  This procedure replaces getc(), and requires no erase processing in stty. */
/*  See init_tty() below. */

{
  int i,stop,ch,cho,chx;

 if (!set_tty) return fil->get();
 else {
  if (charpo==charpi) {
   for (stop=0; !stop; ) {
     ch = fil->get();
     if (ch != EOF) ch  &= ~0x80;
     //ncfprintf (stderr,"'%x'",ch);
     if (ch >= 0x20 && ch < 0x7f) {
	insertc(ch);
	refreshline(fil);
     }
     else {	/* control chars */
       switch (ch) {
         case 0x04: ch = EOF;
		    break;
         case 0x1b: 
		    if((ch=fil->get())=='[') {
		      switch (ch=fil->get()) {
		       case 'A': eraseline(fil);			/* up */
				 linenum = setlinebuf(--linenum);
				 refreshline(fil);
				 break;
		       case 'B': eraseline(fil);			/* down */
				 linenum = setlinebuf(++linenum);
				 refreshline(fil);
				 break;
		       case 'C': if (charind<charpi) 			/* right */
				     charind++;  
				  else insertc(' ');
       		       		  refreshline(fil);
				 break;
		       case 'D': if (charind>0) charind--;		/* left */
       		       		  refreshline(fil);
				  break;
			default: break;
		      }
		    }
		    break;

         case '\r': 
         case '\n': cbuf[charpi++] = ch;
  		    if (charpi >= LINBUFSIZ) charpi = LINBUFSIZ-1;
		    putcc('\n',fil);
		    break;
         case 0x7f: 
         case '\b': if (charind>0) {
		       charind--; 
		       eraseline(fil);
		       deletec();
       		       refreshline(fil);
		    } 
		    break;
         default:   break; 
       }
     }
     switch (ch) {

        default: break; 
      case EOF:  return(ch); break;
      case '\n':
      case '\r': stop = 1; break;
     }
   }
  }
  if (charpo < charpi) {
    cho = cbuf[charpo++];
    if (cho==EOF || charpo>=charpi) {
	setlinebuf (++linenum);
	charpo = charind = charpi = 0;
    }
    //ncfprintf (stderr,"'%x'",cho);
    return cho;
  }
  else return EOF;
 }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void ungetcc(int ch, istream *fil)

/* put back one char in buffer */

{
    int i;

  if (!set_tty) fil->putback(ch);
  else {
    if (charpo==0) {			/* if no space for char */
      for (i=charpi; i>=0; i--) {	/* move buffer one char */
        cbuf[i+1] = cbuf[i];
      }
      charpi++;
      if (charpi >= LINBUFSIZ) charpi = LINBUFSIZ-1;
      charpo = 1;
    }
    cbuf[--charpo] = ch;
  }
};

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

void grabline(istream *fil)

{
   int ch,stop;

   for (stop=0; !stop; ) {
     ch = getcc(fil);
     switch (ch) {

        default: break; 
      case '\n':
      case '\r': stop = 1; break;
     }
   }
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#define FREADSTRLEN 100
char freadstr[FREADSTRLEN]={0};
int freadindex= -1;

int getcx(istream *fil)
             
/* get input char and count but ignore newlines */

{
  int ch;

 if (freadindex>=0) {
   ch = freadstr[freadindex];
   if (++freadindex>=FREADSTRLEN) freadindex--;
 }
 else { 
  ch = getcc(fil);
  charpos++;
  if (ch == '\r') {
    charpos = 0;
    ch = ' ';
  }
  if (ch == '\n') {
    lineno++;
    charpos = 0;
    ch = ' ';
  }
  if (ch=='#' && charpos==1) {
    grabline(fil);
    lineno++;
    charpos = 0;
    ch = ' ';
  }
  if (ch=='/') {				/* '//' is a single line comment  */
    int nextch;				/* Rob & Mark, June '96 	  */
    nextch=getcc(fil);
    if (nextch=='/') {
      grabline(fil);
      lineno++;
      charpos = 0;
      ch = ' ';
    }
    else {
      ungetcx(nextch,fil);
    }
  }   
 }
return (ch);
#undef LINBUFSIZ
}

int getcy(istream *fil)
             
/* get input char and count,  leave newlines */

{
  int ch;

 ch = getcc(fil);
 charpos++;
 if (ch == '\n') {
   lineno++;
   charpos = 0;
 }
return (ch);
}

/*-----------------------------------------------------*/

void ungetcx(int ch, istream *fil)
{
  if (freadindex>=0) {
   if (freadindex>0) freadindex--;
  } else {
    if(charpos>0) charpos--;
    ungetcc(ch,fil);
  }
}

/*-----------------------------------------------------*/

void pushfil (void)
{
    const char *ifil;
    datum d1;

  if (!checkstr(d1=pop())) return;
  ifil = d1.str;
  /*if (cppfl) {
           FILE *ftemp;

     if (ifil) {
         if ((ftemp=freopen (ifil,"r",stdin))==NULL) {
           ncfprintf (stderr,"Error: nc: cpp can't open input file %s\n",ifil);
         }
     }
     if ((ftemp=popen("cpp -P","r"))==NULL) {
         ncfprintf (stderr,"Error: nc: cpp can't open input file\n");
     }
     else tfin = ftemp;
  }*/
  /*else */
  if ((tfin = new ifstream(ifil)) == NULL)
  {
      ncfprintf(stderr,"%s: can't open %s\n", progname, ifil);
      return; 
  }
  if (numfiles < FSTACKSIZ)
  {
      if (fstack[numfiles] == NULL) fstack[numfiles] = new Script;
      fstack[numfiles]->fp = fin;
      fstack[numfiles]->name = einfile;
      fstack[numfiles]->lineno = lineno;
       
      fin = tfin;
      einfile = ifil;
      lineno = 1;
      charpos = 0;
      numfiles++;
      return;
    }
  else
    {
      ncfprintf (stderr,"%s: too many reentrant files\n",progname);
      return;
    }
}

istream *popfil (void)
{
   if (numfiles > 0)
    {
       if (typeid(ifstream *) == typeid(fin)) 
	   ((ifstream *) fin)->close();		/* close old file */
       numfiles--;
       einfile = fstack[numfiles]->name;
       lineno  = fstack[numfiles]->lineno;
       return (fin = fstack[numfiles]->fp); 
    }
   else return ((istream *) NULL);
}

int yylex(void)         /* interp */
{
#define SBUFLEN 4090
	char sbuf[SBUFLEN];
	int c, ch, retval;

	for (c=getcx(fin); c == ' ' || c == '\t'; c=getcx(fin))
		;
	while (c == EOF) {
    	   if (popfil()!=NULL) {
    		for (c=getcx(fin); c == ' ' || c == '\t'; c=getcx(fin))
			;
	   } 
	   else {
		   retval=0;
	 	   break;
	   }
	}
       do {
	if (c == '/') {                   /* comment */ 
  	  if ((ch = getcx(fin)) != '*') {
        	ungetcx(ch,fin);
   	  }
    	  else {                          /* inside comment here */
   	   for (ch=getcx(fin); ch != EOF; ch=getcx(fin))
     	     if (ch == '*')
        	if ((ch = getcx(fin)) == '/') {
	           while ((c=getcx(fin)) == ' ' || c == '\t')
		     ;
	           while (c == EOF)
	             if (popfil())
		       {
			 while ((c=getcx(fin)) == ' ' || c == '\t')
			   ;
		       }
	             else { retval=EOF; break; }	

		    break;       /* out of for (ch=getcx ; ; ) */
        	}		
		else ungetcx(ch,fin);
  	    }    /* else */
 	  }   /* if */
	} while (c == '/' && (ungetcx(ch=getcx(fin),fin),ch=='*'));

	if (c=='.'||isdigit(c)) {           /* number */
	int base=10;
	int h;
    	double d;
    	Symbol *s;
		static char *num,numbuf[120];

        	ungetcx(c, fin);	/* look for num with exponent: */
		num = numbuf;
		c= *num++ =getcx(fin);
		if (isdigit(c))		{c= *num++ =getcx(fin);}
		if (c=='x' || c=='X')	{c= *num++ =getcx(fin); 
		   base = 16;
		   while (isalnum(c))   {c= *num++ =getcx(fin);}
		} else {
		  while (isdigit(c))	{c= *num++ =getcx(fin);}
		  if (c=='.')		{c= *num++ =getcx(fin);}
		  while (isdigit(c))	{c= *num++ =getcx(fin);}
		  if (c=='e'||c=='E')	{c= *num++ =getcx(fin);
   		     if (c=='+'||c=='-')	{c= *num++ =getcx(fin);}
   		     while (isdigit(c))	{c= *num++ =getcx(fin);}
		  }
		}
		ungetcx(c,fin);
		*(--num) = 0;
		switch (base) {
		  case 10: sscanf(numbuf,"%lf",&d); break;
    		  case 16: sscanf(numbuf,"%x",&h); d = h; break;
		}
        	if ((s=lookup(numbuf)) == 0) s = install(numbuf, NUMBER, d);
        	yylval.sym = s;
        	retval= NUMBER;
	}
	else if (isalpha(c)||c=='_') {
        	Symbol *s,*tsym;
		double val;

        	char *p = sbuf;
        	do {
        		if (p >= sbuf + sizeof(sbuf) - 1) {
            			*p = '\0';
                		execerror("name too long", sbuf);
            		}
            		*p++ = c;
        	} while ((c=getcx(fin)) != EOF && (isalnum(c) || c=='_'));
        	ungetcx(c,fin);
        	*p = '\0';

		if ((s=lookup(sbuf)) == 0) {
		  val = getvarval(sbuf);    /* look in list of "command line" vars or predef vars */
       		  if (val!=LARGENUM) {
		     s = install(sbuf,VAR,val);
		     s->vtype = NUMBER;
		  }
		  else s = install(sbuf,UNDEF,LARGENUM);
		  //	s = install(sbuf,UNDEF,LARGENUM);
		}
        	if (formal) {             /* if formal param def */
      		      Symbol *tsym;
  		  if ((tsym=lookupnt(sbuf,ARG)) != 0)
               	       execerror("duplicate formal arg:", sbuf);
          	  s = install(sbuf,ARG,LARGENUM);
        	}
		else if (inlocal) {       /* if local variable */
  		      Symbol *tsym;
		    if (s->type!=DIM) {
  		      if ((tsym=lookupnta(sbuf,LOCALVAR,inlocal)) != 0)
             		execerror("duplicate local variable:", sbuf);
          	      s = install(sbuf,LOCALVAR,LARGENUM);
		    }
      		}
					/* symbol found inside a func */
		if (indef)  { 
					/* check if it is local first */
            	    if      ((tsym=lookupnt(sbuf,ARG)) != 0)      s = tsym; 
		} 
		if ((tsym=lookupnt(sbuf,LOCALVAR)) != 0) s = tsym;
		if ((tsym=lookupnt(sbuf,LOCALARR)) != 0) s = tsym;
		
        	yylval.sym = s;
        	retval= (s->type == UNDEF ? VAR : s->type);
	}
	else if (c == '\'') { /* quoted char */
		c = getcx(fin);		/* get literal char */
		c = backslash(c);
        	//yylval.sym = install("", LITCHAR, (double)(int) c);
		sbuf[0]=c;
		sbuf[1]=0;
        	yylval.sym = install(sbuf, LITCHAR, (double)(int) c);
		c = getcx(fin);		/* move through end quote */
        	retval=LITCHAR;
	}
	else if (c == '"') { /* quoted string */
		char *p;
		for (p = sbuf; (c=getcx(fin)) != '"'; p++) {
	        	if (c == '\n' || c == EOF)
				execerror("missing quote", "");
	        	if (p >= sbuf + sizeof(sbuf) - 1) {
				*p = '\0';
				execerror ("string too long",sbuf);
	        	}
	        	*p = backslash(c);
        	}
        	*p = 0;
		yylval.sym = install (sbuf,STRING,0.0);
		retval=STRING;
	}
	else if (c == '`') { /* command line */
		char *p;
		for (p = sbuf; (c=getcy(fin)) != '\n' && c != '`'; p++) {
			if (c == EOF) break;
			if (p >= sbuf + sizeof(sbuf) - 1) {
		   		*p = '\0';
				execerror ("command line too long",sbuf);
		   	}
			*p = backslash(c);
		}
		*p = 0;
		p = emalloc(strlen(sbuf)+1);
		yylval.sym = (Symbol *)p;
		strcpy(p, sbuf);
		retval=COMNDLINE;
	}
	else switch (c) {
	case '+':       retval=follow('=', ADDEQ, follow ('+', INCROP, '+'));
			 break;
        case '-':       retval=follow('=', SUBEQ, follow ('-', DECROP,
			 			  follow ('>', PFIELD, '-')));
			 break;
		case '*':       retval=follow('=', MULEQ, '*'); break;
        case '/':       retval=follow('=', DIVEQ, '/'); break;
        case '^':       retval=follow('^', XOR, '^'); break;
        case '>':       retval=follow('=', GE, GT); break;
        case '<':       retval=follow('=', LE, LT); break;
        case '=':       retval=follow('=', EQ, '='); break;
        case '!':       retval=follow('=', NE, NOT); break;
        case '|':       retval=follow('=', OREQ,  follow('|', OR, BITOR)); 
			 break;
        case '&':       retval=follow('=', ANDEQ, follow('&', AND, BITAND)); 
			 break;
        case '\n':      lineno++; retval='\n'; break;
        default:        retval=c; break;
	}

   /* ncfprintf (stderr,"yylex token: %d \n", retval); /* */
   return retval;
}


int backslash(int c)    /* get next char with \'s interpreted */
              
{
        static char transtab[] = "b\bf\fn\nr\rt\t";

        if (c != '\\')
                return c;
        c = getcx(fin);
        if (islower(c) && strchr(transtab, c))
                return strchr(transtab, c)[1];
        return c;
}
int follow(int expect, int ifyes, int ifno)     /* look ahead for >=, etc. */
{
        int c;

        c = getcx(fin);
        if (c == expect)
                return ifyes;
        ungetcx(c, fin);
        return ifno;
}

void defnonly(const char *s)     /* warn if illegal definition */
{
        if (!indef)
                execerror(s, "used outside definition");
}

void yyerror(const char *s)      /* report compile-time error */
{
        warning (s, (char *)0);
}

void execerror(const char *s, const char *t) /* recover from run-time error */
{

        warning(s, t);
	charpos = 0;
}

void fpecatch(int i)      /* floating point exceptions */
{
        execerror ("floating point exception", (char *)0);
}

int isa_tty()

/* check to see if the input is from a terminal */
/* the "tty -s" command returns 0 if standard input is a terminal */ 

{
  return (fin== &cin && system("tty -s")==0);
}

void onexit(int val)

{
  if (set_tty) system ("stty sane"); 
  //longjmp(ncbegin, 0); 			/* return from "run_interp()" below */
  exit(val);
}

void init_tty(void)

{
  system ("stty raw");
  system ("stty isig brkint");
  system ("stty -echo -echoctl");
  system ("stty opost");
  set_tty = 1;
}


/*------------------------------------*/

static Symbol *run_on_ex_pnt = (Symbol *)NULL;
static void (*run_on_exit_pnt)()= NULL;

void __runonexit(void)

/* call either interpreted or compiled exit routine */

{
	   double d=0;

  if (run_on_ex_pnt) { callproc (run_on_ex_pnt,0,d,d); }
  if (run_on_exit_pnt) run_on_exit_pnt();
}

/*------------------------------------*/

void set_run_on_exit(void (*runonexpnt)())

{
  run_on_exit_pnt = runonexpnt;
}


/*------------------------------------*/

void set_run_on_ex(void)

/* look for interpreter defined exit routine */

{
    Symbol *f;

  if (f=lookup("run_on_exit")) run_on_ex_pnt = f;
  else                         run_on_ex_pnt = (Symbol*)NULL;
}

/*------------------------------------*/


void onintr(int i)
{
	(void) signal (SIGINT,onintr);   /* "onintr" at ^C */
        //if (++iflag >= 1) onexit(0);
	__runonexit();			/* possibly run user routine before exit */
	//exit(0);
	onexit(0);
}

void run_interp(void)   /* execute until EOF */
{
	//if (setjmp(ncbegin) > 0) return;	/* set jump environment */
	lineno = 1;
	charpos = 0;
	initcode();
        if (isa_tty()) init_tty();
	istat = (int*)signal(SIGINT, SIG_IGN);  /* save original status */
        if (istat != (int*)SIG_IGN) 
	   (void) signal (SIGINT,onintr);   /* "onintr" at ^C */
	  istat = (int *)signal(SIGTERM, SIG_IGN);  /* save original status */
	if (istat != (int*)SIG_IGN)
	   (void) signal (SIGTERM,onintr);      /* "onintr" at TERM */
	  istat = (int *)signal(SIGPIPE, SIG_IGN);  /* save original status */
	if (istat != (int*)SIG_IGN)
	   (void) signal (SIGPIPE,onintr);      /* "onintr" at TERM */
        signal(SIGFPE, fpecatch);
        for (resetcode(); yyparse() && !iflag; resetcode()) {
                execute(prog);
        }
        iflag = 0;
}

/*------------------------------------*/

// char freadstr[FREADSTRLEN]=0;
// int freadindex= -1;
double fread_val=0;

double fread_interp(char *str) /* intepret expressions by fread */

{
    int len;
    Inst *progpf;

  savecode();		/* save the program, stack, and frame pointers */
  len = strlen(str);
  if (len < FREADSTRLEN) {
    strcpy (freadstr,str);
    strcpy (&freadstr[len],"; "); /* add semicolon at end for print */
    freadindex = 0;  	/* set getcx() to read freadstr */
  }
  //fprintf (stderr,"1prog %x progp %x\n",prog,progp);
  progpf = progp;	/* save start of expression */
  yyparse();		/* parse the token from file */
  execute(progpf);	/* run expression, set fread_val in print "expr;" */
  
  //fprintf (stderr,"2prog %x progp %x\n",prog,progp);
  freadindex = -1;     /* switch back to normal getcx() */
  restorecode();       /* restore pointers */
  return fread_val;
}

/*------------------------------------*/

void disperror ()
{
        ncerror = 1;
	disperr();				/* display error message */
}

/*------------------------------------*/

void warning(const char *s, const char *t)   /* print warning message */
                    
{
        int i;

	disperror();
        if (strcmp(einfile,"stdin") == 0) {
          for (i=0; i<charpos-1; i++)           /* print arrow under error */
            ncfprintf (stderr," ");
          ncfprintf (stderr,"^\n");
        }
        if (s)
	   ncfprintf(stderr, "%s: %s", progname, s); /* print lineno with error */
        if (t)
                ncfprintf(stderr, " %s", t);
        if (einfile)
                ncfprintf(stderr," in file '%s'", einfile);
        ncfprintf (stderr, " near line %d char %d.\n", lineno,charpos);
 /*       while (c != ';' && c != EOF)
                c = getcx(fin);  /* flush rest of input line */
}


