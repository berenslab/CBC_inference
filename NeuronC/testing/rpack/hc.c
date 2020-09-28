/* hc - horizontal cat */

#include <stdio.h>
#include <string.h>

FILE *file1,*file2,*fopen();

#define BUFSIZE 512

/****************************************/

main(argc,argv)
   int argc;
   char **argv;

{
   char *cptr;
   FILE *fopen();
	 
 if (argc==1)			/* if user needs help */
    printf ("\nusage: hc file1 file2\n");
 else
 if (argc==3)
 do					/* if there are any switches */
  {
   argc--;
   if (argc)
      {
       argc--; argv++;
       cptr = *argv;
       if((file1=fopen(cptr,"r"))==NULL)
         {
           fprintf(stderr,"hc: cannot open %s\n",cptr);
           continue;
         }

       argc--; argv++;
       cptr = *argv;
       if((file2=fopen(cptr,"r"))==NULL)
         {
           fprintf(stderr,"hc: cannot open %s\n",cptr);
           continue;
         }

    }      /* if (argc) */
   run();
  }
 while (argc > 0);
}


run()

{
  static char linbuf1[BUFSIZE*2]=0;
  static char linbuf2[BUFSIZE]=0;
  static char *s1,*s2,*p1,*p2;
  int i;

 for ( s1=fgets(linbuf1,BUFSIZE,file1),
       s2=fgets(linbuf2,BUFSIZE,file2); 
	s1 || s2;
		   s1 = fgets(linbuf1,BUFSIZE,file1),  
  	   	   s2 = fgets(linbuf2,BUFSIZE,file2) )
  {
  p1 = strchr(linbuf1,'\n');
  if (p1) *p1 = 0;

  p2 = strchr(linbuf2,'\n');
  if (p2) *p2 = 0;

  if (linbuf2[0])  {
    strcat (linbuf1," ");
    strcat (linbuf1,linbuf2);
  }

  printf ("%s\n",linbuf1);
  
  for (i=0; i<BUFSIZE; i++) {
    linbuf1[i] = 0;
    linbuf2[i] = 0;
  }
 } /* while */
 exit (0);
}

