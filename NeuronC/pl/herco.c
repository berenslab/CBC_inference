
/*  Segment HERC.C */

/* System-dependent graphics and text I/O;

   For use with Hercules type graphics cards where text
   can't be displayed at the same time as graphics.

*/

/*  Latest mod  10-June-86	R.G.Smith       */

#include <stdio.h>
#include <fcntl.h>

/********************************************/

/*   Hercules hi-res graphics board driver */

/*   Simple driver to start Hi-Res (720 * 348) mode */

#define VIDBUF 0xb800		/* start of video buffer seg */
#define INDEX  0x3b4		/* index register */
#define CNTRL  0x3b8		/* control register */
#define DSTAT  0x3ba		/* display status register */
#define SOFTSW 0x3bf		/* softw switch port,0=text,1=page1,3=both g */
#define VERTRET 0x80		/* vertical retrace bit in display status reg */

#define VIDBUFSIZE 0x4000	/* size of one page in graphics screen */

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

int iofd;

/********************************************/

gmodeh()

{
   static int i,q;

  if ((iofd = open("/dev/iomem",O_RDWR)) < 0) {
    fprintf (stderr,"Can't open /dev/iomem\n");
    return;
  }

  q = phys(0,0,0x5c0);		/* set 0xb8000 as extra segment page */

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


  if ((iofd = open("/dev/iomem",O_RDWR)) < 0) {
    fprintf (stderr,"Can't open /dev/iomem\n");
    return;
  }
  outb(CNTRL,TEXT);		/* change back to text mode with screen off */

  for(i=0; i<12; i++)		/* set params for text mode */
   outib(INDEX,i,ttable[i]);

  outb(CNTRL,TEXT+SCRNON);	/* screen on, page 0 */

  close (iofd);
}

/**************************************************/

vidclr()

/* erase screen on video board */

{
  register int i,j;

 for (i=0,j=0; i<VIDBUFSIZE; i++,j+=2)
   putesw(0,j);
     
}

/**************************************************/

writdot (loc,bit,color)
  int loc,bit,color;

{
  int val;
  static char bittab[] = {1,2,4,8,0x10,0x20,0x40,0x80};

  val = getesb(loc);
  if (color & 0x80)	   val ^=  bittab[bit];
  else if (color & 0x0f)   val |=  bittab[bit];
  else			   val &= ~bittab[bit];
  putesb(val,loc);

}

/**************************************************/


outb(port,val)
   int port;
   int val;

{
   if (iofd > 0) {
    lseek (iofd,(long)port,0);
    write (iofd,&val,1); 
  }
}

/**************************************************/

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

/**************************************************/

outib(port,index,val)
  int port;
  int index,val;

{
  outb(port,index);
  outb(port+1,val);
}

/********************************************/

