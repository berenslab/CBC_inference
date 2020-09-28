/* textmod */

#include <stdio.h>
#include "../h/vdef.h"

#define TEXT  1
#define VIDEO 2
#define ERASE 4

int scrrev=0;
int scrntyp=SCRNTYP;
int scrcolor;
int operation=TEXT;

FILE *fopen();

main(argc,argv)
int argc; char **argv;
   {
	double atof();
	register char *cptr;

	do	/* loop over arg list */
	   {
		argc--; argv++;
		cptr= *argv;
		if(argc && (*cptr == '-'))
		   {
			cptr++;
			switch(*cptr)
			  {
			  case 'b': 
				scrntyp = MONOCHR;
				break;
			  case 'e': 
				scrntyp = ENHANCED;
				break;
			  case 'g': 
				scrntyp = VGA;
				break;
			  case 'h': 
				scrntyp = HERCULES;
				break;
			  case 't':
				scrntyp = TX4014;
				break;
			  case 'c':
				operation = ERASE;
				break;
			  case 'v':
				operation = VIDEO;
				break;
			  default:
				break;
			  }
			continue;
		   }
		if (operation==TEXT) txmod();
		else if (operation==VIDEO) vidmod();
		else if (operation==ERASE) txeras();

	   } while(argc);
	exit();
   }


txmod()

{
  switch (scrntyp) {
 
  case ENHANCED:   tmodee();
		  break;

  case HERCULES:  tmodeh();
		 break;

  case VGA: 	tmodev();
		 break;

/*  case MONOCHR:  openpl();
		 close(); 
		 system("erase");
		 break;  		/* */
  }
}

vidmod()

{
  switch (scrntyp) {
 
  case ENHANCED:   gmodee();
		  break;

  case HERCULES:  gmodeh();
		 break;

  case VGA:  	gmodev();
		 break;

/*  case MONOCHR:  openpl();
		 closepl();
		break;		/* */
  }
}


txeras()

{
  switch (scrntyp) {
 
  case ENHANCED:   gmodee();
		   clrce();
		  break;

  case HERCULES:  gmodeh();
		  clrch();
		 break;

  case VGA:  	gmodev();
		  clrcv();
		 break;

/*  case MONOCHR:  erasee();
		break;		/* */
  }
}


setorg()
{}

tdot()
{}

