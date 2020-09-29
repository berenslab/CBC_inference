/* Program HPLOT */

/* Driver routines for HP 7221A plotter.  */

/* Latest Mod 	11-Mar-82		R.G.Smith	*/

#include <stdio.h>
#include <sgtty.h>
#include "../h/mdef.h"


#define HANDSH '\05'		/* ENQ control char */
#define ESC    '\033'		/* ESC control char */
#define LF     '\012'
#define CLSBRK '\175'		/* Terminator char for plotter (nop) */
 
#include <signal.h>
void (*istat)();		/* Used by interrupt routine setup below */
				/* see "UNIX Programming, second edition" */

/****************************************/

/* Routine to plot using the HP 7221A plotter.
Turn on the plotter, and initialize the handshake
mode.  Send the standard input to the plotter device
file. 

Format for calling this routine:

	hplot			(use if plotter is connected to /dev/tty2)
	hplot -l		(use from other port if /dev/tty2 is logged in)
	hplot /dev/ttyx		(use if plotter is connected to another port)
	hplot -l /dev/ttyx	(use to run plotter from another tty line)

	where /dev/ttyx is any tty device file

To run:

   Normally, connect plotter to /dev/tty2 and don't log in /dev/tty2;
to run plotter, use  

	"hplot"

   If both terminal and plotter are connected to /dev/tty2, log in
/dev/tty2 and use:

	"hplot"  
 
   If both terminal and plotter are connected to another port, use:

	"hplot /dev/tty" or "hplot /dev/ttyx" (x = 3,4...)

*********************************************************

Explanation of "-l" flag:  (normally don't use this option)

   The problem here is that if the terminal line to the plotter
(and possibly a terminal) is set in "/etc/ttys" to be a login device,
UNIX does not connect the line to this process until a <CR> or <LF>
come in.  Thus we must send an instruction to the plotter to make it
give an extra line of input to the computer, before we can expect to
be able to do handshaking (see pflush() and phand() below). 
   The -l option sends this extra instruction after the plotter has
been turned on, and the output mode has been set, so that the chars
sent will be followed by a '\r' = 015 (13d)
   If the plotter port is not logged in, or if it is logged in and
the terminal is also connected to the same port, the "-l" option
should not be used.
*/

int ch;
TEXT *cptr,*pldev;			/* Plotter device file name */
BOOL lmode;


main(argc,argv)
int argc; char **argv;

begin
    
  lmode = F;				/* Default plotter line isn't login */
  pldev = 0;				/* Default tty line for plotter */

	do	/* loop over arg list */
	   {
		argc--; argv++;
		cptr= *argv;
		if(argc)
		 {
			if (*cptr == '-')
		         {
				cptr++;
				switch(*cptr)
				 {
					case 'L':
			  		case 'l':	     /* special mode */
						lmode = T;
						break;
					case 'i':	     /* no interrupt */
				/*		signal (2,1); */
						break;

				  }		/* switch */
			  }			/* if *cptr */
			else
			 {			/* *cptr not= '-' */
				pldev = cptr;
			 }
		 }				/* if (argc) */
  	} while(argc);
  run();
}

/**************************************/

onintr()

begin
 plout("\033.J",3);			/* abort plotter control seq */
 plout("\033.K",3);			/* abort plotter graphics */
 plout ("v\100",2);				/* put away pen */
 ploff();
 plclose();
 exit();
end

/**************************************/

run()

begin
 istat = signal (SIGINT, SIG_IGN);	/* save orig staus */
 if (istat != SIG_IGN)
 	 signal (SIGINT, onintr);	/* set up interrupt routine */

  if (not plopen(pldev)) return;	/* Open plotter device file */
  plon();				/* Turn on plotter */
  phand();				/* Initialize handshake mode */
  if (lmode) plout ("\033.B",3);	/* Tell plotter to send buffer space */
					/* This connects plotter to fdp file */
  while ((ch = getc (stdin)) != EOF)	/* Read chars one at a time;  */
   begin
    if (ch == 'z') pflush();		/* flush pen move (end of cur file) */
    else pbuf (ch);			/* Send char to plotter */ 
   end
  pbuf ('}');				/* Send nop to terminate last cmd */
  pflush ();				/* Flush out chars left in buffer */
  ploff();				/* Turn off plotter */
  plclose();				/* Reset plotter file */
end



