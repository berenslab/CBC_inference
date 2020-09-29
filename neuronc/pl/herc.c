/* Driver for Hercules card */

/* Draws line segments in color on Hercules Hi-Res board 
*/

/*  Latest mod  22-Feb-87	R.G.Smith       */

/* colors */

#define BLACK   0
#define BLUE    1
#define GREEN   2
#define CYAN    3
#define RED     4
#define MAGENTA 5
#define BROWN   6
#define WHITE   7
#define GRAY    8
#define LTBLUE  9
#define LTGREEN 10
#define LTCYAN  11
#define LTRED   12
#define LTMAG   13
#define YELLOW  14
#define BRTWHT  15


#define XHI 479
#define YHI 347

#define STDSCL (STDSIZ/348) 		/* scale from 348 to STDSIZ */

#define abs(a)  ((a) < 0 ? -(a) : (a))
#define outb(a,b) io_outb(a,b)
#define inb(a)    io_inb(a)

static int x,y;			/* point values for line drawing */
static int scrcolor;

char inb();

/********************************************/

/*   Hercules hi-res graphics board driver */

/*   Simple driver to start Hi-Res (720 * 348) mode */

#define VIDBUF 0xb8000		/* start of video buffer seg */
#define INDEX  0x3b4		/* index register */
#define CNTRL  0x3b8		/* control register */
#define DSTAT  0x3ba		/* display status register */
#define SOFTSW 0x3bf		/* softw switch port,0=text,1=page1,3=both g */
#define VERTRET 0x80		/* vertical retrace bit in display status reg */

#define VIDBUFSIZE 2048 	/* size of one page in graphics screen */

#define SCRNON 8
#define GRPH   0x82
#define TEXT   0x20

static char gtable[12] = {		/* contents for 6845 registers */
			0x35,0x2d,0x2e,0x07,
			0x5b,0x02,0x57,0x57,
			0x02,0x03,0x00,0x00 };

static char ttable[12] = {0x61,0x50,0x52,0x0f,
	 		0x19,0x06,0x19,0x19,
	 		0x02,0x0d,0x0b,0x0c };

/********************************************/

hercioctl(dev,cmd,addr)
int dev;
int cmd;
int *addr;

/* all Herc control routines are ioctl calls.

   cmd

   0	set Herc back to text mode
   1	set Herc to graphics mode
   2    draw vector from (addr[0],addr[1]) to (addr[2],addr[3])
	 with color addr[4]
   3    erase screen to black

*/   

{
   static int *lindat;
   static int x1,y1,x2,y2;

   switch (cmd) {

   case 0:	tmodeh();		/* return to text mode */
		break;

   case 1:	gmodeh();		/* set Herc to graphics mode */
		break;

   case 2:	lindat = addr;		/* make line from (x1,y1) to (x2,y2) */
		x1 = fuword(lindat++);
		y1 = fuword(lindat++);
		x2 = fuword(lindat++);
		y2 = fuword(lindat++);
 		y1 = YHI - y1;
 		y2 = YHI - y2;
		scrcolor = fuword(lindat++);
		linsegh(x1,y1,x2,y2);
		break;

   case 3:	clrch();
		break;
  }
}
		
/**************************************************/


linsegh (x1,y1,x2,y2)
   int x1,y1,x2,y2;
   
/* Bresenham's Algorithm;
   see Newman and Sproul, 2nd Ed.  for explanation.

 */

{
        static int deltax,deltay,dx,dy;
        static int incr1,incr2,d,xend,yend;
        static int xinc,yinc;

 deltax = x2 - x1;
 deltay = y2 - y1;
 dx = abs(deltax); dy = abs(deltay);

 if (dy <= dx)
   {
     if (deltax < 0)
       {
         yinc = ((deltay < 0) ? 1 : -1); 
         x = x2; y = y2; xend = x1;
       }
     else   /* deltax >= 0 */
       {
         yinc = ((deltay < 0) ? -1 : 1);
         x = x1; y = y1; xend = x2;
       }
     pointh();
     d = (dy << 1) - dx;
     incr1 = dy << 1;
     incr2 = (dy - dx) << 1;
     while (x < xend)
       {
         x++;
         if (d < 0) d += incr1;
         else     { y += yinc; d += incr2; }
         pointh();
       }
    } 
  else   /* dy > dx */
    {
      if (deltay < 0)
        {
          xinc = ((deltax < 0) ? 1 : -1);
          x = x2; y = y2; yend = y1;
        }
      else
        {
          xinc = ((deltax < 0) ? -1 : 1);
          x = x1; y = y1; yend = y2;
        }
      pointh();
      d = (dx << 1) - dy;
      incr1 = dx << 1;
      incr2 = (dx - dy) << 1;
      while (y < yend)
       {
         y++;
         if (d < 0) d += incr1;
         else     { x += xinc; d += incr2; }
         pointh();
       }
   }  /* else */

} /* drawto */

/********************************************/

gmodeh()

{
   static int i;

/*  if ((iofd = open("/dev/iomem",O_RDWR)) < 0) {
 *   fprintf (stderr,"Can't open /dev/iomem\n");
 *   return;
 * }
 * q = phys(0,0,0x5c0);		/* set 0xb8000 as extra segment page */
 

  outb(SOFTSW,3);		/* allow graphics mode */

  outb(CNTRL,GRPH);		/* set to graphics mode with screen off */

  for(i=0; i<12; i++)		/* set params for 720 x 348 graphics mode */
   outib(INDEX,i,gtable[i]);

  outb(CNTRL,GRPH+SCRNON);	/* screen on, page 0 */
}

/********************************************/

tmodeh()

{
   static int i;


/*  if ((iofd = open("/dev/iomem",O_RDWR)) < 0) {
 *   fprintf (stderr,"Can't open /dev/iomem\n");
 *   return;
 * }
 */

  outb(CNTRL,TEXT);		/* change back to text mode with screen off */

  for(i=0; i<12; i++)		/* set params for text mode */
   outib(INDEX,i,ttable[i]);

  outb(CNTRL,TEXT+SCRNON);	/* screen on, page 0 */

}

/**************************************************/

clrch()

/* erase screen on Hercules board */

{
  register int i,j;
  long int loc;
  static int erase[16] = 0;

 loc = ((long)VIDBUF);
 for (i=0,j=0; i<VIDBUFSIZE; i++,j+=16)
  kpcopy (erase,loc+j,16);			/* send color byte */

}

/**************************************************/

pointh ()

{
  register int bit,ty;
  int val;
  long int loc;
  static char bittab[] = {1,2,4,8,0x10,0x20,0x40,0x80};

  
  loc = ((long)VIDBUF) + ((y&3) << 13) + 90 * (y >> 2) + (x>>3);
  bit = 7-(x&7);
  pkcopy (loc,&val,1);				/* latch byte */
  if (scrcolor & 0x80)	      val ^=  bittab[bit];
  else if (scrcolor & 0x0f)   val |=  bittab[bit];
  else			      val &= ~bittab[bit];
  kpcopy (&val,loc,1);

}

/**************************************************/

/*
outb(port,val)
   int port;
   int val;

{
   if (iofd > 0) {
    lseek (iofd,(long)port,0);
    write (iofd,&val,1); 
  }
}
*/

/**************************************************/

/*
char inb(port)

{
  static char val;

 if (iofd > 0) {
   lseek (iofd,(long)port,0);
   read(iofd,&val,1);
   return val;
 }
 else return ((char)0);
}
*/

/**************************************************/

outib(port,index,val)
  int port;
  int index,val;

{
  outb(port,index);
  outb(port+1,val);
}


