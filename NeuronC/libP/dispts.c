# include "graph.h"
# include "stdplt.h"

void zinit(char*,char*,char*);
void yyparse(void);
void pair(void);
void Do_disp(void);

void _dispts(char *strp, int doscale)
{
	char astring[100];

	zinit("%$ %D",strp,astring);
	yyparse();
	pair();
	Do_disp();
  }

void dispts(char *strp, int args)
{
	argp = (union ptr *)&args;
	_dispts(strp,1);
  }
