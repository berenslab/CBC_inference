
# line 2 "graph.y"
# include "graph.h"
# include "stdio.h"
struct axtype *_axptr;
struct _data *_datptr;
int starflg;
int yyline;
int xyflg;
int oxyflg;
int i,j;
char *charp,*chrp;
char *sp;
int frag;
double atof(const char *);
#define yyclearin yychar = -1
#define yyerrok yyerrflag = 0
extern int yychar;
extern int yyerrflag;
#ifndef YYMAXDEPTH
#define YYMAXDEPTH 150
#endif
#ifndef YYSTYPE
#define YYSTYPE int
#endif
YYSTYPE yylval, yyval;
typedef int yytabelem;
# define YYERRCODE 256

# line 369 "graph.y"



yylex(void)
{
	extern int yylval;
	starflg = 0;
	while((*scanp != '%') && (*scanp != '\0'))
	  {
		ychar++;
		scanp++;
	  }
	if(*scanp++ == '\0')
		return(0);
	ychar++;
	*(scanp - 1) = 0;	/* To make previous string null terminated */
	if(*scanp == '*')
	  {
		(starflg)++;
		scanp++;
		ychar++;
	  }
	if(starflg)
		yylval = (int)(scanp + 1); /* pointer to string after token */
	else
		yylval = (int)argp;
	ychar++;
	return((int)*scanp++);
  }
yyerror(char *s)
{
	char *str;

	ychar -= 5;
	str = strbeg;
	fprintf(stderr,"\n%s\n:  ",s);
	fprintf(stderr,"error near underbar\n");
	while(str - strbeg < ychar)
		putc(*str++,stderr);
	fprintf(stderr,"_%s\n",scanp);
  }
yytabelem yyexca[] ={
-1, 1,
	0, -1,
	-2, 0,
	};
# define YYNPROD 72
# define YYLAST 211
yytabelem yyact[]={

    50,     3,    59,    77,    56,    65,    29,    60,    11,    13,
    42,    62,    19,     5,     9,    73,     6,    46,    39,     4,
    31,    41,    28,    23,    26,    32,    33,    30,    44,    24,
     7,    27,    48,    47,    53,    49,    55,    21,    58,    54,
    12,    14,    37,    61,    64,    63,    51,    66,    16,    45,
    38,    57,    52,    40,    59,    67,    56,    65,    18,    60,
    15,     8,     2,    62,    25,     1,     0,     0,    59,    43,
    56,     0,     0,    60,    11,    13,    42,    62,     0,    79,
     0,     0,     0,     0,    39,    36,    53,    41,    55,     0,
    58,    54,     0,     0,    35,    61,    64,    63,    68,    75,
    53,     0,    55,    57,    58,    54,    12,    14,    37,    61,
    64,    63,     0,     0,    71,    76,    38,    57,    50,    40,
    59,     0,    56,    70,    74,    60,     0,     0,     0,    62,
    34,     0,     0,     0,     0,     0,     0,     0,     0,    59,
     0,    56,    65,     0,    60,     0,     0,     0,    62,    75,
    48,    47,    53,    49,    55,     0,    58,    54,     0,    69,
    78,    61,    64,    63,    51,    76,     0,     0,     0,    57,
    52,    53,     0,    55,    74,    58,    54,    10,     0,     0,
    61,    64,    63,    17,    20,    22,     0,     0,    57,     0,
     0,     0,     0,    17,     0,     0,    20,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,    72,     0,     0,
    78 };
yytabelem yypact[]={

 -1000, -1000,   -35,   -52,   -80,   -80,   -80,   -80,   -80, -1000,
 -1000, -1000, -1000, -1000, -1000,   -80, -1000, -1000,   -80, -1000,
 -1000, -1000, -1000, -1000,   -66, -1000,     0, -1000,    52,   -14,
 -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,
 -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,
 -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,
 -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,
 -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,    71 };
yytabelem yypgo[]={

     0,    65,    62,    61,    60,    58,    37,    14,   177,    29,
    27,    20,    25,    26,    48,    24,   130,    94,    85,    69,
     3,    12,    22,     6,    15 };
yytabelem yyr1[]={

     0,     1,     1,     1,     1,     2,     3,     3,     7,     9,
     9,    10,    10,    10,     8,     8,     8,     8,     4,     4,
    14,    15,    15,    11,    11,    11,    11,    11,    11,    11,
    11,    11,    11,    12,    12,    12,     5,     5,    21,    22,
    22,    22,    22,    22,    22,    13,    13,    13,    13,    13,
    13,    19,    16,    16,    16,    16,    16,    16,    16,    16,
    17,    17,    20,    18,     6,    23,    23,    24,    24,    24,
    24,    24 };
yytabelem yyr2[]={

     0,     9,     9,     9,     9,     1,     5,     3,     5,     5,
     1,     3,     3,     3,     3,     3,     3,     3,     5,     3,
     5,     5,     1,     3,     3,     3,     3,     3,     3,     3,
     3,     3,     3,     3,     3,     3,     5,     3,     5,     5,
     5,     5,     5,     5,     1,     3,     3,     3,     3,     3,
     3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
     3,     3,     3,     3,     9,     5,     1,     3,     3,     3,
     3,     3 };
yytabelem yychk[]={

 -1000,    -1,    -2,    36,    71,    65,    68,    82,    -3,    -7,
    -8,    88,   120,    89,   121,    -4,   -14,    -8,    -5,   -21,
    -8,    -6,    -8,    -7,    -9,   -14,   -15,   -21,   -22,   -23,
   -10,   -11,   -12,   -13,   -16,   -17,   -18,   108,   116,    84,
   119,    87,    76,   -19,   -20,   115,    83,    99,    98,   101,
    66,   112,   118,   100,   105,   102,    70,   117,   104,    68,
    73,   109,    77,   111,   110,    71,   -11,   -13,   -19,   -16,
   -17,   -18,    -8,   -24,   -17,   -19,   -18,   -20,   -16,   -23 };
yytabelem yydef[]={

     5,    -2,     0,     0,     0,     0,     0,     0,     1,     7,
    10,    14,    15,    16,    17,     2,    19,    22,     3,    37,
    44,     4,    66,     6,     8,    18,    20,    36,    38,     0,
     9,    11,    12,    13,    23,    24,    25,    26,    27,    28,
    29,    30,    31,    32,    33,    34,    35,    45,    46,    47,
    48,    49,    50,    52,    53,    54,    55,    56,    57,    58,
    59,    60,    61,    63,    51,    62,    21,    39,    40,    41,
    42,    43,    66,    65,    67,    68,    69,    70,    71,    64 };
/* @(#)yaccpar */

/*
** Skeleton parser driver for yacc output
*/

/*
** yacc user known macros and defines
*/
#define YYERROR		goto yyerrlab
#define YYACCEPT	return(0)
#define YYABORT		return(1)
#define YYBACKUP( newtoken, newvalue )\
{\
	if ( yychar >= 0 || ( yyr2[ yytmp ] >> 1 ) != 1 )\
	{\
		yyerror( "syntax error - cannot backup" );\
		goto yyerrlab;\
	}\
	yychar = newtoken;\
	yystate = *yyps;\
	yylval = newvalue;\
	goto yynewstate;\
}
#define YYRECOVERING()	(!!yyerrflag)
#ifndef YYDEBUG
/*#	define YYDEBUG	1	/* make debugging available */
#endif

/*
** user known globals
*/
int yydebug;			/* set to 1 to get debugging */

/*
** driver internal defines
*/
#define YYFLAG		(-1000)

/*
** global variables used by the parser
*/
YYSTYPE yyv[ YYMAXDEPTH ];	/* value stack */
int yys[ YYMAXDEPTH ];		/* state stack */

YYSTYPE *yypv;			/* top of value stack */
int *yyps;			/* top of state stack */

int yystate;			/* current state */
int yytmp;			/* extra var (lasts between blocks) */

int yynerrs;			/* number of errors */
int yyerrflag;			/* error recovery flag */
int yychar;			/* current input token number */



/*
** yyparse - return 0 if worked, 1 if syntax error not recovered from
*/
int
yyparse()
{
	register YYSTYPE *yypvt;	/* top of value stack for $vars */

	/*
	** Initialize externals - yyparse may be called more than once
	*/
	yypv = &yyv[-1];
	yyps = &yys[-1];
	yystate = 0;
	yytmp = 0;
	yynerrs = 0;
	yyerrflag = 0;
	yychar = -1;

	goto yystack;
	{
		register YYSTYPE *yy_pv;	/* top of value stack */
		register int *yy_ps;		/* top of state stack */
		register int yy_state;		/* current state */
		register int  yy_n;		/* internal state number info */

		/*
		** get globals into registers.
		** branch to here only if YYBACKUP was called.
		*/
	yynewstate:
		yy_pv = yypv;
		yy_ps = yyps;
		yy_state = yystate;
		goto yy_newstate;

		/*
		** get globals into registers.
		** either we just started, or we just finished a reduction
		*/
	yystack:
		yy_pv = yypv;
		yy_ps = yyps;
		yy_state = yystate;

		/*
		** top of for (;;) loop while no reductions done
		*/
	yy_stack:
		/*
		** put a state and value onto the stacks
		*/
#if YYDEBUG
		/*
		** if debugging, look up token value in list of value vs.
		** name pairs.  0 and negative (-1) are special values.
		** Note: linear search is used since time is not a real
		** consideration while debugging.
		*/
		if ( yydebug )
		{
			register int yy_i;

			printf( "State %d, token ", yy_state );
			if ( yychar == 0 )
				printf( "end-of-file\n" );
			else if ( yychar < 0 )
				printf( "-none-\n" );
			else
			{
				for ( yy_i = 0; yytoks[yy_i].t_val >= 0;
					yy_i++ )
				{
					if ( yytoks[yy_i].t_val == yychar )
						break;
				}
				printf( "%s\n", yytoks[yy_i].t_name );
			}
		}
#endif /* YYDEBUG */
		if ( ++yy_ps >= &yys[ YYMAXDEPTH ] )	/* room on stack? */
		{
			yyerror( "yacc stack overflow" );
			YYABORT;
		}
		*yy_ps = yy_state;
		*++yy_pv = yyval;

		/*
		** we have a new state - find out what to do
		*/
	yy_newstate:
		if ( ( yy_n = yypact[ yy_state ] ) <= YYFLAG )
			goto yydefault;		/* simple state */
#if YYDEBUG
		/*
		** if debugging, need to mark whether new token grabbed
		*/
		yytmp = yychar < 0;
#endif
		if ( ( yychar < 0 ) && ( ( yychar = yylex() ) < 0 ) )
			yychar = 0;		/* reached EOF */
#if YYDEBUG
		if ( yydebug && yytmp )
		{
			register int yy_i;

			printf( "Received token " );
			if ( yychar == 0 )
				printf( "end-of-file\n" );
			else if ( yychar < 0 )
				printf( "-none-\n" );
			else
			{
				for ( yy_i = 0; yytoks[yy_i].t_val >= 0;
					yy_i++ )
				{
					if ( yytoks[yy_i].t_val == yychar )
						break;
				}
				printf( "%s\n", yytoks[yy_i].t_name );
			}
		}
#endif /* YYDEBUG */
		if ( ( ( yy_n += yychar ) < 0 ) || ( yy_n >= YYLAST ) )
			goto yydefault;
		if ( yychk[ yy_n = yyact[ yy_n ] ] == yychar )	/*valid shift*/
		{
			yychar = -1;
			yyval = yylval;
			yy_state = yy_n;
			if ( yyerrflag > 0 )
				yyerrflag--;
			goto yy_stack;
		}

	yydefault:
		if ( ( yy_n = yydef[ yy_state ] ) == -2 )
		{
#if YYDEBUG
			yytmp = yychar < 0;
#endif
			if ( ( yychar < 0 ) && ( ( yychar = yylex() ) < 0 ) )
				yychar = 0;		/* reached EOF */
#if YYDEBUG
			if ( yydebug && yytmp )
			{
				register int yy_i;

				printf( "Received token " );
				if ( yychar == 0 )
					printf( "end-of-file\n" );
				else if ( yychar < 0 )
					printf( "-none-\n" );
				else
				{
					for ( yy_i = 0;
						yytoks[yy_i].t_val >= 0;
						yy_i++ )
					{
						if ( yytoks[yy_i].t_val
							== yychar )
						{
							break;
						}
					}
					printf( "%s\n", yytoks[yy_i].t_name );
				}
			}
#endif /* YYDEBUG */
			/*
			** look through exception table
			*/
			{
				register int *yyxi = yyexca;

				while ( ( *yyxi != -1 ) ||
					( yyxi[1] != yy_state ) )
				{
					yyxi += 2;
				}
				while ( ( *(yyxi += 2) >= 0 ) &&
					( *yyxi != yychar ) )
					;
				if ( ( yy_n = yyxi[1] ) < 0 )
					YYACCEPT;
			}
		}

		/*
		** check for syntax error
		*/
		if ( yy_n == 0 )	/* have an error */
		{
			/* no worry about speed here! */
			switch ( yyerrflag )
			{
			case 0:		/* new error */
				yyerror( "syntax error" );
				goto skip_init;
			yyerrlab:
				/*
				** get globals into registers.
				** we have a user generated syntax type error
				*/
				yy_pv = yypv;
				yy_ps = yyps;
				yy_state = yystate;
				yynerrs++;
			skip_init:
			case 1:
			case 2:		/* incompletely recovered error */
					/* try again... */
				yyerrflag = 3;
				/*
				** find state where "error" is a legal
				** shift action
				*/
				while ( yy_ps >= yys )
				{
					yy_n = yypact[ *yy_ps ] + YYERRCODE;
					if ( yy_n >= 0 && yy_n < YYLAST &&
						yychk[yyact[yy_n]] == YYERRCODE)					{
						/*
						** simulate shift of "error"
						*/
						yy_state = yyact[ yy_n ];
						goto yy_stack;
					}
					/*
					** current state has no shift on
					** "error", pop stack
					*/
#if YYDEBUG
#	define _POP_ "Error recovery pops state %d, uncovers state %d\n"
					if ( yydebug )
						printf( _POP_, *yy_ps,
							yy_ps[-1] );
#	undef _POP_
#endif
					yy_ps--;
					yy_pv--;
				}
				/*
				** there is no state on stack with "error" as
				** a valid shift.  give up.
				*/
				YYABORT;
			case 3:		/* no shift yet; eat a token */
#if YYDEBUG
				/*
				** if debugging, look up token in list of
				** pairs.  0 and negative shouldn't occur,
				** but since timing doesn't matter when
				** debugging, it doesn't hurt to leave the
				** tests here.
				*/
				if ( yydebug )
				{
					register int yy_i;

					printf( "Error recovery discards " );
					if ( yychar == 0 )
						printf( "token end-of-file\n" );
					else if ( yychar < 0 )
						printf( "token -none-\n" );
					else
					{
						for ( yy_i = 0;
							yytoks[yy_i].t_val >= 0;
							yy_i++ )
						{
							if ( yytoks[yy_i].t_val
								== yychar )
							{
								break;
							}
						}
						printf( "token %s\n",
							yytoks[yy_i].t_name );
					}
				}
#endif /* YYDEBUG */
				if ( yychar == 0 )	/* reached EOF. quit */
					YYABORT;
				yychar = -1;
				goto yy_newstate;
			}
		}/* end if ( yy_n == 0 ) */
		/*
		** reduction by production yy_n
		** put stack tops, etc. so things right after switch
		*/
#if YYDEBUG
		/*
		** if debugging, print the string that is the user's
		** specification of the reduction which is just about
		** to be done.
		*/
		if ( yydebug )
			printf( "Reduce by (%d) \"%s\"\n",
				yy_n, yyreds[ yy_n ] );
#endif
		yytmp = yy_n;			/* value to switch over */
		yypvt = yy_pv;			/* $vars top of value stack */
		/*
		** Look in goto table for next state
		** Sorry about using yy_state here as temporary
		** register variable, but why not, if it works...
		** If yyr2[ yy_n ] doesn't have the low order bit
		** set, then there is no action to be done for
		** this reduction.  So, no saving & unsaving of
		** registers done.  The only difference between the
		** code just after the if and the body of the if is
		** the goto yy_stack in the body.  This way the test
		** can be made before the choice of what to do is needed.
		*/
		{
			/* length of production doubled with extra bit */
			register int yy_len = yyr2[ yy_n ];

			if ( !( yy_len & 01 ) )
			{
				yy_len >>= 1;
				yyval = ( yy_pv -= yy_len )[1];	/* $$ = $1 */
				yy_state = yypgo[ yy_n = yyr1[ yy_n ] ] +
					*( yy_ps -= yy_len ) + 1;
				if ( yy_state >= YYLAST ||
					yychk[ yy_state =
					yyact[ yy_state ] ] != -yy_n )
				{
					yy_state = yyact[ yypgo[ yy_n ] ];
				}
				goto yy_stack;
			}
			yy_len >>= 1;
			yyval = ( yy_pv -= yy_len )[1];	/* $$ = $1 */
			yy_state = yypgo[ yy_n = yyr1[ yy_n ] ] +
				*( yy_ps -= yy_len ) + 1;
			if ( yy_state >= YYLAST ||
				yychk[ yy_state = yyact[ yy_state ] ] != -yy_n )
			{
				yy_state = yyact[ yypgo[ yy_n ] ];
			}
		}
					/* save until reenter driver code */
		yystate = yy_state;
		yyps = yy_ps;
		yypv = yy_pv;
	}
	/*
	** code supplied by user is placed in this switch
	*/
	switch( yytmp )
	{
		
case 5:
# line 22 "graph.y"
{ i = j = 0;} break;
case 14:
# line 35 "graph.y"
 {	
		nary[X]++;
		xyflg = X;
		_datptr = _xdata[i++] = (struct _data *)calloc(1,sizeof(struct _data));
		if(_datptr == 0)
			_feror("No space available: calloc");
		_datptr->flgs |= MINMAX;
		if(starflg)
			_datptr->flgs |= NODATA;
		else
			_datptr->dptr = (argp++)->dp;
		_axptr = &_ax[xyflg];
		} break;
case 15:
# line 48 "graph.y"
 {	
		nary[X]++;
		xyflg = X;
		_datptr = _xdata[i++] = (struct _data *)calloc(1,sizeof(struct _data));
		if(_datptr == 0)
			_feror("No space available: calloc");
		if(starflg)
			_datptr->flgs |= NODATA;
		else
			_datptr->dptr = (argp++)->dp;
		_axptr = &_ax[xyflg];
		} break;
case 16:
# line 60 "graph.y"
 {
		nary[Y]++;
		xyflg = Y;
		_datptr = _ydata[j++] = (struct _data *)calloc(1,sizeof(struct _data));
		if(_datptr == 0)
			_feror("No space available: calloc");
		_datptr->flgs |= MINMAX;
		if(starflg)
			_datptr->flgs |= NODATA;
		else
			_datptr->dptr = (argp++)->dp;
		_axptr = &_ax[xyflg];
		} break;
case 17:
# line 73 "graph.y"
 {
		nary[Y]++;
		xyflg = Y;
		_datptr = _ydata[j++] = (struct _data *)calloc(1,sizeof(struct _data));
		if(_datptr == 0)
			_feror("No space available: calloc");
		if(starflg)
			_datptr->flgs |= NODATA;
		else
			_datptr->dptr = (argp++)->dp;
		_axptr = &_ax[xyflg];
		} break;
case 21:
# line 91 "graph.y"
{
		if(first)
		  {
			if(xyflg != oxyflg)
			  {
				yyerror("Saxis: Non matching xy descriptor");
				exit(1);
			  }
		  }
		first = 1;
		oxyflg = xyflg;
		} break;
case 22:
# line 103 "graph.y"
 {
		if(first)
		  {
			if(xyflg != oxyflg)
			  {
				yyerror("Saxis: Non matching xy descriptor");
				exit(1);
			  }
		  }
		first = 1;
		oxyflg = xyflg;
	  } break;
case 26:
# line 120 "graph.y"
{	chrp = (char *)yypvt[-0];
		if(starflg)
		  {
			if(xyflg == X)
				charp = _axptr->lab = xlabbuf;
			else
				charp = _axptr->lab = ylabbuf;
			while((*chrp == ' ') || (*chrp == '\t'))
				chrp++;
			while((*chrp != '%') && (*chrp != '\0'))
				*charp++ = *chrp++;
			*charp = '\0';
		  }
		else
			_axptr->lab = (argp++)->cp;} break;
case 27:
# line 135 "graph.y"
{	
		if(starflg)
			frag = atoi(yypvt[-0]);
		else
			frag = *(argp++)->ip;
		if((frag >= 0) && (frag <= 3))
			_axptr->flags |= frag;
		else
			printf("Bad Tic (%t) type %d\n",frag); } break;
case 28:
# line 144 "graph.y"
{	
		if(starflg)
			frag = atoi(yypvt[-0]);
		else
			frag = *(argp++)->ip;
		if((frag >= 0) && (frag <= 3))
			_axptr->flags |= frag << 2;
		else
			printf("Bad Tic (%t) type %d\n",frag); } break;
case 29:
# line 153 "graph.y"
{	
		if(starflg)
			frag = atoi(yypvt[-0]);
		else
			frag = *(argp++)->ip;
		if((frag >= 0) && (frag <= 7))
			_axptr->dflags |= frag;
		else
			printf("Bad size (%s) type %d\n",frag); } break;
case 30:
# line 162 "graph.y"
{	
		if(starflg)
			frag = atoi(yypvt[-0]);
		else
			frag = *(argp++)->ip;
		if((frag >= 0) && (frag <= 7))
			_axptr->dflags |= frag << 3;
		else
			printf("Bad size (%S) type %d\n",frag); } break;
case 31:
# line 171 "graph.y"
{	chrp = (char *)yypvt[-0];
		_axptr->flags |= NAME;
		if(starflg)
		  {
			if(xyflg == X)
				charp = _axptr->lab = xlabbuf;
			else
				charp = _axptr->lab = ylabbuf;
			while((*chrp == ' ') || (*chrp == '\t'))
				chrp++;
			while((*chrp != '%') && (*chrp != '\0'))
				*charp++ = *chrp++;
			*charp = '\0';
		  }
		else
			_axptr->lab = (argp++)->cp;} break;
case 34:
# line 190 "graph.y"
 {
		if(starflg)
		  {
			if(xyflg == X)
				_axptr->y = atof(yypvt[-0]);	
			else
				_axptr->x = atof(yypvt[-0]);
		  }
		else
		  {
			if(xyflg == X)
				_axptr->y = *((double *)argp); /*argp->dv;*/
			else
				_axptr->x = *((double *)argp); /* argp->dv; */
			argp += 4;
		  }
		_axptr->flags |= SHIFT;
	  } break;
case 35:
# line 208 "graph.y"
 {
		if(starflg)
		  {
			if(xyflg == X)
				_axptr->y = atof(yypvt[-0]);	
			else
				_axptr->x = atof(yypvt[-0]);
		  }
		else
		  {
			if(xyflg == X)
				_axptr->y = *((double *)argp); /*argp->dv; */
			else
				_axptr->x = *((double *)argp); /*argp->dv; */
			argp += 4;
		  }
		_axptr->flags |= SHIFT1;
	  } break;
case 45:
# line 243 "graph.y"
 {
		if(starflg)
		  {
			for(sp = (char *)yypvt[-0];*sp != '%';sp++)
				if(*sp == ' ')
					continue;
				else
					break;
			_datptr->symble = *sp;
		  }
		else
			_datptr->symble = *(argp++)->cp;
		} break;
case 46:
# line 256 "graph.y"
 {
		if(starflg)
			_datptr->blubber = atoi(yypvt[-0]);
		else
			_datptr->blubber = *(argp++)->ip;
		} break;
case 47:
# line 262 "graph.y"
{	chrp = (char *)yypvt[-0];
		if(starflg)
		  {
			charp = _datptr->labl = (char *)calloc(1,30);
			if(charp == 0)
				_feror("No space available: calloc");
			while((*chrp == ' ') || (*chrp == '\t'))
				chrp++;
			while((*chrp != '%') && (*chrp != '\0'))
				*charp++ = *chrp++;
			*charp = '\0';
		  }
		else
			_datptr->labl = (argp++)->cp;} break;
case 48:
# line 276 "graph.y"
 {
		if(starflg)
		  {
		_datptr->ddash = (float *)calloc(1,16);
		if(_datptr->ddash == 0)
			_feror("No space available: calloc");
		if(sscanf(yypvt[-0],"%f%f%f%f",_datptr->ddash,_datptr->ddash+1,
		  _datptr->ddash + 2,_datptr->ddash + 3) != 4)
		  {
		   fprintf(stderr,"Not enough arguments for 'B' command\n");	
		   fprintf(stderr,"%f %f %f %f\n",*_datptr->ddash
			,*(_datptr->ddash+1),
			*(_datptr->ddash + 2),*(_datptr->ddash + 3));
			_datptr->ddash = 0;
		  }
		  }
		else
			_datptr->ddash = (argp++)->fp;
		} break;
case 49:
# line 295 "graph.y"
{					/* pen changes */
		if(starflg)
			frag = atoi(yypvt[-0]);
		else
			frag = *(argp++)->ip;
		if((frag >= 0) && (frag <= 3))
			_datptr->flgs |= frag << 4; /* shift 3 left 4 for 060 */
	    } break;
case 50:
# line 303 "graph.y"
{					/* char size for point label */
		if(starflg)
			frag = atoi(yypvt[-0]);
		else
			frag = *(argp++)->ip;
		if((frag >= 0) && (frag <= 7))
			_datptr->flgs |= frag; 
	    } break;
case 51:
# line 311 "graph.y"
{
		if(starflg)
			_datptr->npts = atoi(yypvt[-0]);
		else
			_datptr->npts = *(argp++)->ip; } break;
case 52:
# line 318 "graph.y"
{_datptr->flgs |= INTEGER;} break;
case 53:
# line 319 "graph.y"
{_datptr->flgs |= INTEGER;} break;
case 54:
# line 320 "graph.y"
{_datptr->flgs |= FLOAT;} break;
case 55:
# line 321 "graph.y"
{_datptr->flgs |= DOUBLE;} break;
case 56:
# line 322 "graph.y"
{_datptr->flgs |= UNSIGNED;} break;
case 57:
# line 323 "graph.y"
{_datptr->flgs |= SHORT;} break;
case 58:
# line 324 "graph.y"
{_datptr->flgs |= LONG;} break;
case 59:
# line 325 "graph.y"
{_datptr->flgs |= LONG;} break;
case 60:
# line 327 "graph.y"
{
		_axptr->flags |= MINSET;
		if(starflg)
			min[xyflg] = atof(yypvt[-0]);
		else
		  {
			min[xyflg] = *(double *)argp;  /* argp->dv; */
			argp += 4;
		  }
	     } break;
case 61:
# line 338 "graph.y"
{	_axptr->flags |= MAXSET;
		if(starflg)
			max[xyflg] = atof(yypvt[-0]);
		else
		  {
			max[xyflg] = *(double *)argp;  /* argp->dv; */
			argp += 4;
		  }
	   } break;
case 62:
# line 347 "graph.y"
{
		if(starflg)
			frag = atoi(yypvt[-0]);
		else
			frag = *(argp++)->ip;
		if((frag >= 0) && (frag <= 3))
			_axptr->flags |= frag << 5;
		else
			printf("Bad Gridtype (%%G) type %o\n",frag);
	   } break;
case 63:
# line 357 "graph.y"
{
		_axptr->flags |= LOGAXIS;
	    } break;
	}
	goto yystack;		/* reset registers in driver code */
}
