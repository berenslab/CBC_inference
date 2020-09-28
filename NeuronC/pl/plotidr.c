/* Segment PLOTDR */

/* Low level driver routines for HP 7475 plotter.  */

/* Latest Mod 	10-Mar-82		R.G.Smith	*/

#include <stdio.h>
#include <fcntl.h>
#include <termio.h>
#include "../h/mdef.h"


#define HANDSH '\05'		/* ENQ control char */
#define ESC    '\033'		/* ESC control char */
#define LF     '\012'
#define CLSBRK '\175'		/* Terminator char for plotter (nop) */
 

/* Plotter variables */

static TEXT plotbuf [BUFSIZE];
static int bufcnt = 0;
static int fdp = 0;
static char comnd[80];
static BOOL plset = 0;

static struct termio ttysto=0,ttystn=0;

/********************************************/

int plopen(dev)
    TEXT *dev;

/* Open plotter port. Default device is /dev/com2. */

begin
  extern int fdp;

 if ((dev == NULL) or (*dev == NULL)) dev = "/dev/tty01";
 if ((fdp = open (dev,O_RDWR|O_SYNC)) == -1)
   begin
    chmod(dev,0666);
    if ((fdp = open (dev,O_RDWR|O_SYNC)) == -1) begin
      fprintf (stderr, "Plotter port does not open\r\n");
      return (F);
    end
   end 

 if (not plset) ioctl (fdp,TCGETA,&ttysto);		/* get old speed */
 ioctl (fdp,TCGETA,&ttystn);
 ttystn.c_cflag = B9600 | CS8 | CREAD | HUPCL | CLOCAL;	/* set 9600 baud */
 ttystn.c_lflag &= ~(ECHO | ECHOK);
 ttystn.c_iflag |= IXON | IXOFF;
 
 ioctl (fdp,TCSETA,&ttystn);				/* set new speed */

 plset = T;
 return (T);
end

/****************************************/

plclose()

begin
  extern int fdp;

 if (plset) ioctl (fdp,TCSETAW,&ttysto);		/* set old speed */
 plset = F;
 close (fdp);
end

/****************************************/

plout (str,n)
   TEXT *str;
   int n; 

/* Send a string to plotter */

begin
  extern int fdp;
  
  write (fdp,str,n);
end

/****************************************/

char plin ()

/* Get a char sent from plotter */

begin
   char ch;
   extern int fdp;

  read (fdp, &ch, 1);
  return (ch);
end

/****************************************/

plon()

/* Send command to turn on plotter */

begin
  plout ("\033.Y",3);
end

/****************************************/

ploff()

begin
  plout ("\033.Z",3);
end

/****************************************/

pflush()

/* Flush the plotter buffer. */

/* Send the handshake char to the plotter.
   This tells the plotter to send the return handshake
string back to the computer.  Because UNIX saves all
chars input in the char queue, we must ignore the extra
chars we receive when the terminal is not a login device.
In is case, the first chars we receive are some numbers
and then a '\r' (see explanation in main() above).  
   The handshake string is set to be "g\r" (2 chars), which
is easily discriminated from numerical chars.
   After the handshake string is received, send the buffer up
to the buffer index, and then clear the buffer.
 
*/


begin
end

/****************************************/

phand()

/* Routine to initialize plotter handshake mode.
Set plotter output mode, and handshake char.
After ENQ char '\05' is sent to plotter, it
responds with a 'g' when space for another block
of commands is available in its memory.
See HP 7221A manual for details. */

begin
/* plout ("\033.M;;10;13:",11);	/* Send output mode command seq. */
 				/* no delay or trigger; echo term = <LF>. */
				/* Output terminate char is <CR>. */

 plout ("\033.I81;;17:",11); 	/* Send handshake mode command */
 			/* Use Xon-Xoff (17-19) software handshake mode */
 plout ("\033.N;19:",7);	/* set Xoff trigger character */

 plout (";",1);		/* Terminator for commands (NOP) */
end

/****************************************/


