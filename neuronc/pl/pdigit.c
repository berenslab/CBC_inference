/* Program PDIGIT */

/* Digitize program for HP 7221A plotter.  */

/* Latest Mod 	10-Mar-82		R.G.Smith	*/

#include <stdio.h>
#include <sgtty.h>
#include "../h/mdef.h"


#define HANDSH '\05'		/* ENQ control char */
#define ESC    '\033'		/* ESC control char */
#define LF     '\012'
#define CLSBRK '\0175'		/* Terminator char for plotter (nop) */
 

/* Plotter variables */

static int fdp = 0;

#include <signal.h>			/* used by main for user interrupts */
#include <setjmp.h>
jmp_buf  sjbuf; 
void (*istat)();			/* look in "UNIX programming" */

extern int iflag = 0;			/* interrupt flag used by mdisplay */

double fmod();

/****************************************/

onintr()

begin
  extern BOOL iflag;

  if (iflag++ > 1)
    begin
     plout("\033.J",3);			/* abort plotter control seq */
     plout("\033.K",3);			/* abort plotter graphics */
     plout ("v\100",2);				/* put away pen */
     ploff();
     pflush();
     exit(); 
    end
  signal (SIGINT,onintr);		/* go to "onintr" at ^C */
end

/****************************************/

main (argc,argv)
    int argc;
    TEXT *argv [];

/* Routine to digitize points using the HP 7221A plotter.
   To run, turn on plotter, and connect to a serial port,
   preferably tty2.  Run "pdigit" and the plotter will 
   initialize, and the "enter" light on the plotter will
   light up, not flashing.  Whenever you press the "enter"
   button, pdigit prints the coordinates of the point on the
   standard output.
*/

/*
Format for calling this routine:

	pdigit       		(use when plotter is connected to tty)
	pdigit -l		(don't use )
	pdigit /dev/ttyx	(use to connect to different device file) 
	pdigit -l /dev/ttyx

	where /dev/ttyx is any tty device file

   There are two safe modes to use this program:

	1) Connect tty to plotter, connect plotter to /dev/tty2
line, login on tty2, and run "pdigit".

	2) Set entry for tty2 in /etc/ttys for no login. 
Connect plotter to /dev/tty2 line, login on any other terminal,
and run "pdigit /dev/ttyp", where ttyp is the plotter tty line.
 
*/

begin
   int ch;
   TEXT *str,*pldev;			/* Plotter device file name */
   BOOL lmode;
   TEXT point [20];
   int x,y;
   extern void (*istat)(),onintr();	/* look in "UNIX programming" */
 
  lmode = F;				/* Default plotter line isn't login */
  pldev = NULL;			/* Default tty line for plotter */
  if (argc >= 2)
    begin
     str = argv[1];
     if (argc == 3) pldev = argv[2]; 
     if (*str == '-')			/* If first arg is flag ... */
      begin
       if ((*(++str) == 'l') or (*str == 'L'))
        lmode = T;
       else
        begin
         fprintf (stderr, "Unknown flag: %s\r\n",argv[1]);
         return;
        end
      end
     else if (argc == 2) pldev = argv[1]; 
    end

  if (not plopen(pldev)) return;	/* Open plotter device file */
  plon();				/* Turn on plotter */
  phand();				/* Initialize handshake mode */
  preset();
  plimit(520,380,15720,10380);		/* Cal plotter by thousandths of inch */
  pgrid(15200,10000);
  if (lmode) plout ("\033.B}",4);	/* Tell plotter to send buffer space */
					/* This connects plotter to fdp file */
  pflush();				/* send cmds from pset2 above */ 
  
  iflag = F;
  istat = signal (SIGINT, SIG_IGN);	/* save original status */
  if (istat != SIG_IGN)
     signal (SIGINT,onintr);		/* go to "onintr" at ^C */

  while (iflag == 0)			/* Read chars one at a time;  */
   begin
    plout ("\033.D}",4);		/* Send digitize command to plotter */
    str = point;
    while ((*str++ = plin()) != '\n');	/* get plotter's response */
    *str = NULL;
    sscanf (point,"%d,%d",&x,&y);
    ploff();
    fprintf (stdout,"%d %d\n",x,y);
    fprintf (stderr,"\r");
    plon();
  end
  pbuf ('}');				/* Send nop to terminate last cmd */
  pflush ();				/* Flush out chars left in buffer */
  ploff();				/* Turn off plotter */
  plclose();				/* Reset plotter file */
end





/****************************************/

psend(ch)
   char ch;

/* char send routine needed by plotsub */

begin
  pbuf(ch);
end



