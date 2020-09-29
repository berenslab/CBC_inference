/* symbol.c */

#include "nc.h"
#include "y.tab.h"
#include "ncsetvar.h"

#define HASHSIZ 199                     /* prime number to make better hash */

static Symbol *hashtab[HASHSIZ] = {0};  /* symbol hash table */
static int hashmade = 0;

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <string.h>

#ifdef __cplusplus
}
#endif

#include "ncio.h"

char *emalloc(unsigned int n);
void efree(void *ptr);
void erarr (ncarray *arrp);

/*---------------------------------------------------*/

void inithash(void)
{
   int i;

   for (i=0; i<HASHSIZ; i++) 
     hashtab[i] = 0;
   hashmade = 1;
}

/*---------------------------------------------------*/

int hashsym(const char *s)
             

/* make disordered index of symbol name */

{
   int i;

   if (!s) return (0);
   for (i=0; *s; s++) i += *s;
   return (i%HASHSIZ); 
}

/*---------------------------------------------------*/

void hinstall (Symbol *sp)
              

/* Install a symbol in hash table */
/*  Uses "direct chaining" in a doubly-linked list. */
/*  Store latest entries at beginning of list */

{
    int i;
    Symbol *spnt;

   sp->next = NULL;
   sp->last = NULL;
   i = hashsym(sp->name);               /* initial index into hashtab */ 
   if (!(spnt=hashtab[i])) {            /* install directly in table */
      hashtab[i] = sp;
   }
   else {                               /* otherwise, put at start hashtab entry */
      sp->next = hashtab[i];
      hashtab[i]->last = sp;
      hashtab[i] = sp;
   }
}

/*---------------------------------------------------*/

Symbol *lookup(const char *s)      /* Find s in symbol table, given its name */
			     /*  Ignore STRING entries */                
{
        Symbol *sp;
        int i;
        
        i = hashsym(s);
        for (sp = hashtab[i]; sp != (Symbol *)0; sp = sp->next) {
                if (sp->name) {
		        if (sp->type==STRING) continue;
			if (strcmp(sp->name, s) == 0) return sp;
		}
        }
        return (Symbol *)NULL;    /* 0 ==> not found */
}

/*---------------------------------------------------*/

Symbol *lookups(const char *s)      /* Find s in symbol table, given its name */
{
        Symbol *sp;
        int i;
        
        i = hashsym(s);
        for (sp = hashtab[i]; sp != (Symbol *)0; sp = sp->next) {
                if (sp->name) {
			if (strcmp(sp->name, s) == 0) return sp;
		}
        }
        return (Symbol *)NULL;    /* 0 ==> not found */
}

/*---------------------------------------------------*/

Symbol *lookupnt(const char *s, int t)/* find s in symbol table, given name and type */
			        /*  Ignore STRING entries */                
	      
{
        Symbol *sp;
        int i;
        
        i = hashsym(s);
        for (sp = hashtab[i]; sp != (Symbol *)0; sp = sp->next) {
            if (sp->name) {
		if (sp->type==STRING) continue;
                if ((strcmp(sp->name, s) == 0) && (sp->type == t))
                        return sp;
	    }
        }
        return (Symbol *)NULL;    /* 0 ==> not found */
}

/*---------------------------------------------------*/

Symbol *lookupnta(const char *s, int t, int starg) 
	/* find s in symbol table, given its name, type, and 
           beginning argnum for its { } block. */
	      
{
        Symbol *sp;
        int i;
        
        i = hashsym(s);
        for (sp = hashtab[i]; sp != (Symbol *)0; sp = sp->next) {
                if ((strcmp(sp->name, s) == 0) && 
		    (sp->type == t) &&
		    (sp->argnum>=starg))
                        return sp;
        }
        return (Symbol *)NULL;    /* 0 ==> not found */
}

/*---------------------------------------------------*/

double getval(const char *s)        /* Find symbol in symbol table, given its name */
			/*  and return its value */
                
{
        Symbol *sp;
        int i,found=0;
	double retval;

        i = hashsym(s);
        for (sp = hashtab[i]; sp != (Symbol *)0; sp = sp->next) {
                if (sp->name)
                  if (strcmp(sp->name, s) == 0) {
                        retval= sp->val;
			found = 1;
			break;
		  }
        }
	if (!found) {
	  ncfprintf (stderr,"getval: can't find '%s'\n",s);
	  retval = 0;
	}
	// fprintf (stderr,"%22s %10.4g \n",sp->name,retval);
        return retval;       /* not found ==> error, large num */
}

/*---------------------------------------------------*/

datum printsym()

{
	datum d = {0};
        Symbol *sp;
        int i,j;

    for (i=j=0; i<HASHSIZ; i++) {
        for (sp = hashtab[i]; sp != (Symbol *)0; sp = sp->next,j++) {
		if (sp->type==BLTIN)
                  ncfprintf (stdout,"%s\n",sp->name);
		else
                  ncfprintf (stdout,"%s %g\n",sp->name,sp->val);
	}
    }
    d.val = j;
    d.vtype = NUMBER;
    return d;
}

/*---------------------------------------------------*/

char *findsym(int num)              /* find a symbol's name, given its type */
                
{
        Symbol *sp;
        int i;

        for (i=0; i<HASHSIZ; i++) 
          for (sp = hashtab[i]; sp != (Symbol *)0; sp = sp->next)
                if (sp->type == num)
                        return sp->name;
        return 0;       /* 0 ==> not found */
}

/*---------------------------------------------------*/

Symbol *findarg(int num)            /* find an argument given its number */
                

{
        Symbol *sp;
        int i;

        for (i=0; i<HASHSIZ; i++) 
          for (sp = hashtab[i]; sp != (Symbol *)0; sp = sp->next)
                if (sp->argnum == num)
                        return sp;
        return 0;       /* 0 ==> not found */
}

/*---------------------------------------------------*/

void erasarg(int argnum)           /* erase all arguments after func def */

{                               /*  erase all local variables, too. */
        Symbol *sp,*oldsp,*next;
        int i;

        for (i=0; i<HASHSIZ; i++) 
          for (oldsp = sp = hashtab[i]; sp != (Symbol *)0; sp = next)
            {
                next = sp->next;
/* if (sp->type==LOCALVAR) ncfprintf (stderr,"%s %d ",sp->name,sp->argnum); */
/* if (sp->type==ARG)ncfprintf(stderr,"%s %d %d * ",sp->name,sp->argnum,argnum);*/

                if ((sp->type == ARG || sp->type==LOCALVAR ||
			sp->type == LOCALARR) && sp->argnum >= argnum)
                 {
                   if (sp == hashtab[i]) {
			if (next) next->last = NULL; 
			hashtab[i] = next;
		   }
                   else {
			oldsp->next = next;
			if (next) next->last = oldsp;
		   }
                   if (sp->name) efree(sp->name);
                   efree(sp);
                 }
                else
                   oldsp = sp; 
           }
 /* ncfprintf (stderr,"\n"); */
}

/*---------------------------------------------------*/

Symbol *install(const char *s, int t, double d)        /* install s in symbol table */
                
{
        Symbol *sp;

	if (!(sp=lookupnt(s,t))) {
          sp = (Symbol *) emalloc(sizeof(Symbol));
          sp->name = emalloc(strlen(s)+1); /* +1 for '\0' */
          strcpy(sp->name, s);
          sp->type = t;
          sp->vtype = 0;
          sp->ntype = 0;
          sp->val = d;
          hinstall (sp);
        }
/*      ncfprintf (stderr,"'%s' %ld\n",s,sp);     /* */
        return sp;
}

/*---------------------------------------------------*/

void rmsym(Symbol *spnt)

/* remove a symbol */

{
    int i;

 if (!hashmade) inithash();
 if (spnt) { 
     i = hashsym(spnt->name);               /* initial index into hashtab */ 
     if (hashtab[i] == spnt) {
        hashtab[i] = spnt->next;
     } 
     if (spnt->last) spnt->last->next = spnt->next;
     if (spnt->next) spnt->next->last = spnt->last;
     if (spnt->type==ARRAY) erarr (spnt->arrp);
     spnt->type = UNDEF;
     //efree (spnt->name);// don't erase name, just leave for new def
     //efree (spnt);      // don't remove symbol because 
			  // pre-compiled code points to it
 }
}


/*---------------------------------------------------*/

void erasym(const char *s)

/* erase a symbol in table given its name */

{
    Symbol *spnt;

  if (spnt=lookup(s)) rmsym(spnt);   /* Find s in symbol table, given its name */
};

/*---------------------------------------------------*/

void erasymtab(void)

/* erase symbol table */

{
   Symbol *next,*spnt;
   int i;

  if (!hashmade) inithash();
  for (i=0; i<HASHSIZ; i++) {
     for (spnt=hashtab[i]; spnt; ) {
       next = spnt->next;
       efree (spnt->name);
       if (spnt->type==ARRAY) erarr (spnt->arrp);
       efree (spnt);
       spnt = next;
     } 
     hashtab[i] = 0;
  }
}


/*---------------------------------------------------*/
/* Set symbol table and variable from symbol table pointer */

void setval(Symbol *sym, double val, int vtyp)         

/* set value of symbol, given its pointer. */

{ 
    if (sym) {
        if (vtyp!=STRING) sym->val = val;
	if (vtyp) sym->vtype = (short int)vtyp;	/* if value type is set */
	if (sym->dptr!=NULL) {
	   switch (sym->ntype) {       
		case VINT:    *sym->iptr = val; break;
		case VDOUBLE: 
		default:      *sym->dptr = val; break;
	   }
	}
	// fprintf (stderr,"%20s %10.4g %x\n",sym->name,val,sym->dptr);
    }
}

/*---------------------------------------------------*/

void setval(Symbol *sym, int val, int vtyp)         

/* set value of symbol, given its pointer. */

{ 
    if (sym) {
        if (vtyp!=STRING) sym->val = val;
	if (vtyp) sym->vtype = (short int)vtyp;	/* if value type is set */
	if (sym->iptr!=NULL) *sym->iptr = val;
	// fprintf (stderr,"%20s %10.4g %x\n",sym->name,val,sym->iptr);
    }
}


/*---------------------------------------------------*/

void setval(Symbol *sym, char *cptr, int vtyp)         

/* set value of symbol, given its pointer. */

{ 
    if (sym) {
        if (vtyp==STRING) sym->str = cptr;
	if (vtyp) sym->vtype = (short int)vtyp;	/* if value type is set */
	if (sym->sptr!=NULL) *sym->sptr = cptr; 
	// fprintf (stderr,"%20s %10 %x\n",sym->name,cptr,sym->sptr);
    }
}



/*---------------------------------------------------*/

/* Set pointer to variable for later access */

void setval(Symbol *sym, double *val, int vtyp)         

/* set value of symbol, given its pointer. */

{ 
    if (sym) {
        if (vtyp!=STRING) sym->val = *val;
	if (vtyp) sym->vtype = (short int)vtyp;	/* if value type is set */
	sym->dptr = val;
	sym->ntype = VDOUBLE;
    }
}

/*---------------------------------------------------*/

void setval(Symbol *sym, int *val, int vtyp)         

/* set value of symbol, given its pointer. */

{ 
    if (sym) {
        if (vtyp!=STRING) sym->val = *val;
	if (vtyp) sym->vtype = (short int)vtyp;	/* if value type is set */
	sym->iptr = val;
	sym->ntype = VINT;
    }
}


/*---------------------------------------------------*/

void setval(Symbol *sym, char **cptr, int vtyp)         

/* set value of symbol, given its pointer. */

{ 
    if (sym) {
        if (vtyp==STRING) sym->str = *cptr;
	if (vtyp) sym->vtype = (short int)vtyp;	/* if value type is set */
	sym->sptr = cptr;
	sym->ntype = VSTRING;
    }
}


