/* Segment PLOTDR */

/* Low level driver routines for HP 7221A plotter.  */

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
 ttystn.c_lflag &= ~(ICANON | ISIG);		/* sg_flags |= RAW; */
 ttystn.c_cc[VMIN] = 1;			/* return after VMIN chars */
 ttystn.c_cc[VTIME] = 0;			/* no timeouts */
 ttystn.c_cflag = B2400 | CS8 | CREAD | HUPCL | CLOCAL;	/* set 2400 baud */
 ttystn.c_lflag &= ~(ECHO | ECHOK);
 
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
   int i;

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
  ioctl (fdp,TCFLSH,2);			/* flush tty input and output queues */
  plout ("\033.(",3);
end

/****************************************/

ploff()

begin
    int i;

  plout ("\033.)",3);
  for (i=0; i<64; i++)
     plout ("\0\0\0\0\0\0\0\0)",8);
  ioctl (fdp,TCFLSH,2);			/* flush tty input and output queues */
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
   int i;
   char ch;

  ch = HANDSH;
  plout (&ch,1);		/* Send handshake char */
  for (i=bufcnt; i<BUFSIZE; i++) /* Fill rest of buffer with terminator */
    plotbuf [i] = CLSBRK;
  ch = plin();
  if (ch <= '9')
   while ((ch = plin()) <= '9');  /* Skip numbers (the buffer size report) */
/* printf ("%x ",ch);		/* debug */
  ch = plin();			/* Wait for return to OK sending buf */
/* printf ("%x ",ch);		/* debug */
  plout (plotbuf, bufcnt);	/* Send entire buffer; last char is NOP */
  if (bufcnt < BUFSIZE-1)
    plout ("\175",1);		/* Send NOP terminator if buffer not full */
  bufcnt = 0;
end

/****************************************/

pbuf (ch)
   char ch;

/* Routine to put char. in output buffer
for plotter.  When buffer is full, it is
dumped to the standard output. One char is
left at end of buffer for a terminator 
char for last command in buffer.  Use this
routine for all graphics commands to plotter.
*/

begin
 if (bufcnt <  (BUFSIZE-1)) plotbuf [bufcnt++] = ch;
 if (bufcnt >= (BUFSIZE-1))   	/* Is buffer full?      */
   pflush ();			/* Yes, dump it to file */
end

/****************************************/

pbuf2 (ch1,ch2)
    char ch1,ch2;

/* Send two chars to buffer. */

begin
  pbuf (ch1);
  pbuf (ch2);
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
 plout ("\033.M:",11);	/* Send output mode command seq. */
 			/* no delay or trigger; no echo term char. */
			/* Output terminate char is <CR> (default). */

 plout ("\033.H512;5;103:",13);  /* Send handshake mode command */
 			/* 512 = Block size for handshake query */
 			/* ENQ = handshake enable char (05). */
 			/* 'g' is handshake return for OK to send (103d) */
			/* Handshake is ended with out term char <CR> */
 plout ("}",1);		/* Terminator for commands (NOP) */
end

/****************************************/


