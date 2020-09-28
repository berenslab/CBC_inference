
/* sample header HEADER.C */

#include <stdio.h>
#include "../h/mdef.h"

FILE *pictin;
char inbuf[BUFSIZE];

/****************************************/

main(argc,argv)
   int argc;
   char **argv;

begin
   FAST TEXT *cptr;
   int i;
   FILE *temp,*freopen();
	 
 pictin = stdin;
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

        end  /* switch */
      end	   /* if */
     else
      begin
       if((pictin=fopen(cptr,"r"))==NULL)
         begin
           fprintf(stderr,"Header: cannot open %s\n",cptr);
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

/* get user input from stdin */


begin
/* label(pictin); */
end

