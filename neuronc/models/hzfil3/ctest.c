/* Program colortest */
/* Makes a color chip chart for */
/*  color PostScript printer. */

/* Prints PostScript program on standard output. */

/* Latest Mod	Aug 92 R.G.Smith */

#include <stdio.h>

FILE *pictin;

double start;
int startfl;

double atof();
double sqrt();

/****************************************/

main(argc,argv)
   int argc;
   char **argv;

{
   char *cptr;
   int i;
   FILE *temp,*freopen();
	 
 pictin = stdin;
 if (argc==1)			/* if user needs help */
   run();
 else
 do					/* if there are any switches */
  {
   argc--; argv++;
   cptr = *argv;
   if (argc)
    {
     if (*cptr == '-')
      {
       cptr++;
       switch (*cptr)
        {

          case 'i':
		signal (2,1);
		break;
     
      	  case 'b': 
		argv++; argc--;
		start = atof(*argv);
		startfl = 1;
		break;

      	  case 'e': 
		break;

	  default:
		fprintf (stderr,"colort: unknown switch '%s'\n",*argv);
		exit();

        }  /* switch */
      }	   /* if */
     else
      {
       if((pictin=fopen(cptr,"r"))==NULL)
         {
           fprintf(stderr,"oplist: cannot open %s\n",cptr);
	   fflush (stderr);
           continue;
         }
       run();
       if (argc <= 1) break;
     }
    }
   else run();
  }
 while (argc > 0);
}

/****************************************/

#define printbox(x,y,r,g,b) printb((x),(y),(double)(r),(double)(g),(double)(b))

run()

{
   int i,j,k;
   int col,row,xloc,yloc,yo;
   double r,g,b,bstart;

printf ("%%!PS-Adobe1.1\n");
printf ("%%%%Creator: R.G. Smith.\n");
printf ("%%%%Title: Color Test\n");
printf ("%%%%CreationDate: Mon Aug 10 16:23:29 1992\n");
printf ("%%%%Pages: 1\n");
printf ("%%%%BoundingBox: 0 0 3000 3000\n");
printf ("%%%%EndComments\n");
printf ("%%%%EndProlog\n");
printf ("%%%%Page: 0 1\n\n");

printf ("newpath\n");
printf ("save\n");
printf ("%% This is for text at the top of page..\n");
printf ("/CommentFont\n");
printf ("  /Helvetica findfont 10 scalefont def\n");
printf ("/Date (Mon Aug 10 16:23:29 1992\n");
printf (") def\n");
printf ("%% This is where you can scale, rotate and translate the image\n");
printf ("%%     if you feel the need to manipulate the image further...\n");
printf ("0 rotate\n");
printf ("50 108 translate\n");
printf ("1 1 scale\n");
printf ("1 setlinewidth\n");
printf ("stroke\n");
printf ("0 690 moveto\n");
printf ("CommentFont setfont\n");
printf ("Date show\n");
printf ("%% Color data begins here...\n\n");

printf ("/mt {moveto} def\n");
printf ("/lt {lineto} def\n");
printf ("/str 20 string def\n");
printf ("/pr {str cvs show } def\n");
printf ("/Helvetica findfont 7 scalefont setfont \n");
printf ("\n");
printf ("/cbox {\n");
printf ("gsave\n");
printf ("3 copy\n");
printf ("gsave\n");
printf ("  35 0 rmoveto\n");
printf ("  40 0 rlineto\n");
printf ("  0 40 rlineto\n");
printf ("  -40 0 rlineto\n");
printf ("  0 -40 rlineto\n");
printf ("  closepath\n");
printf ("  setrgbcolor\n");
printf ("  fill\n");
printf ("grestore\n");
printf ("0 setgray\n");
printf ("0 20 rmoveto \n");
printf ("3 1 roll\n");
printf ("2 1 roll\n");
printf ("pr ( )show pr ( )show  pr\n");
printf ("stroke\n");
printf ("grestore\n");
printf ("} def\n");
printf ("\n");

printf ("0.5 setlinewidth\n");

printf ("0 -50 translate\n\n");

printf ("\n%% first, the subtractive primaries:\n\n");

printbox (0,  680,1,0,1);
printbox (90, 680,1,1,0);
printbox (180,680,0,1,1);

printf ("\n%% next, the additive primaries:\n\n");

printbox (0,  630,1,0,0);
printbox (90, 630,0,1,0);
printbox (180,630,0,0,1);

printf ("\n");

bstart = 0;
for (b=bstart,yo=0; yo<600; b+=0.2,yo+=300) {
 for (xloc=col=0; col<6; col++,xloc+=90) {
  g = col * .2;
  for (yloc=0,row=0; row<6; row++,yloc+=50) {
     r = row * .2; 
     printbox (xloc,yloc+yo,r,g,b);
  }
 }
}

printf ("stroke showpage\n");
printf ("\n");
printf ("20 50 translate\n\n");

bstart = .4;
for (b=bstart,yo=0; yo<600; b+=0.2,yo+=300) {
 for (xloc=col=0; col<6; col++,xloc+=90) {
  g = col * .2;
  for (yloc=0,row=0; row<6; row++,yloc+=50) {
     r = row * .2; 
     printbox (xloc,yloc+yo,r,g,b);
  }
 }
}

printf ("stroke showpage\n");
printf ("\n");
printf ("20 50 translate\n\n");

bstart = .8;
for (b=bstart,yo=0; yo<600; b+=0.2,yo+=300) {
 for (xloc=col=0; col<6; col++,xloc+=90) {
  g = col * .2;
  for (yloc=0,row=0; row<6; row++,yloc+=50) {
     r = row * .2; 
     printbox (xloc,yloc+yo,r,g,b);
  }
 }
}

printf ("stroke showpage\n");

printf ("restore\n");
printf ("clear\n");

}

/****************************************/

printb(xloc,yloc,r,g,b)
   int xloc,yloc;
   double r,g,b;

{
  printf ("%-3d %-3d mt ",xloc,yloc);
  printf ("%3g %3g %3g cbox\n",r,g,b);
}

