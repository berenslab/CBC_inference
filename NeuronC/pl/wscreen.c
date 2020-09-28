
/* scrtest.c */

/* Runs PC monochrome video screen */

#include <stdio.h>
#include "../h/mdef.h"

int xmcur,ymcur;		/* loc of text cursor on monochrome screen */
int mloc;			/* loc of text cursor in monochr window */

#define SCRCOLS 80		/* size of screen */
#define SCRROWS 25
#define LINELEN 160		/* length of chars + attrib for one line */
#define MBASE   0xb000		/* segment for mono screen */

#define CR	0x0d
#define LF	0x0a
#define TAB	0x09
#define SPACE	0x20

#define WHITE 	7

FILE *instream;
char inbuf[BUFSIZE];

/****************************************/

main(argc,argv)
   int argc;
   char **argv;

begin
   FAST TEXT *cptr;
   int i;
   FILE *temp,*freopen();
	 
 instream = stdin;
 if (argc==1)			/* if user needs help */
    run();
 else
 do					/* if there are any switches */
  begin
   argc--; argv++;
   cptr = *argv;
   if (argc)
    begin
     if (*cptr == '-')
      begin
       cptr++;
       switch (*cptr)
        begin

          case 'i':
		signal (2,1);
		break;

          case 'e':
		winit();
		exit();
		break;

        end  /* switch */
      end	   /* if */
     else
      begin
       if((instream=fopen(cptr,"r"))==NULL)
         begin
           fprintf(stderr,"Wscreen: cannot open %s\n",cptr);
	   fflush (stderr);
           continue;
         end
       run();
       if (argc <= 1) break;
     end
    end
   else run();
  end
 while (argc > 0);
end


/****************************************/

run()

/* get user input from input stream */

{
 int c;

 winit();
 while ((c = getc(instream)) != EOF)
   putcm(c);
 exit (0);
}

/*-----------------------------------------*/

winit()
{
    int i,endloc;

  xmcur = ymcur = 0;
  mloc = 0;
  endloc = SCRROWS * LINELEN;
  for (i=0; i<endloc; i+=2)		/* erase the screen to spaces */
	memwrit (i,MBASE,SPACE);  

/*  for (i=LINELEN; i<endloc; i+=2)
	memwrit (i,MBASE,'b'); 
*/
}

/*-----------------------------------------*/
 
xymcur(x,y)
  int x,y;
{
  while (x >= SCRCOLS)
    {
      x -= SCRCOLS;
      y++;
    }
  while (y >= SCRROWS) 
     {
        mscroll();
        y--;
     }
  mloc = x + SCRCOLS * y; 
  xmcur = x;
  ymcur = y;
}

/*-----------------------------------------*/

putcm (ch)
   int ch;

{
	int xtab;

  switch (ch)
   {
 
   case CR:	xymcur (0,ymcur);
		break;
   case LF:	xymcur (0,ymcur+1);
		break;
   case TAB:	xtab = ((xmcur+1) / 8 + 1) * 8;
		xymcur (xtab,ymcur); 
		break;
   default:	memwrit (mloc<<1,MBASE,ch);	  /* even addrs = char */
		memwrit ((mloc<<1)+1,MBASE,WHITE);  /* odd  addrs = attribs */
		xymcur(xmcur+1,ymcur);
		break;
  }
}


/*-----------------------------------------*/

mscroll()

{
    int i,len;

  len = (SCRROWS-1) * LINELEN;
  blkmov(MBASE, LINELEN, MBASE, 0, 1, 1, 1, len);
  for (i=0; i<LINELEN; i+=2)
    memwrit (len+i,MBASE,SPACE);		/* erase bottom line */
    
}

