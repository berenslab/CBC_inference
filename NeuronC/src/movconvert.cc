/* Program movconvert */

/* integrates changes in separate frames for movie */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

double atof(const char *);

#ifdef __cplusplus
}
#endif


#define min(a,b) ((a)>(b) ? (b) : (a))
#define max(a,b) ((a)>(b) ? (a) : (b))

char *filen;
static int nw=1;        /* =1 -> 1 byte, 0-255; =2 bytes, -> 65535 */
static int nframes=1;
static char backgr=0;
static char white=255;
void run(void);

FILE *textout;
FILE *outfil;
FILE *fp;

/* -------------------------------------------------------------- */

int main(int argc, char **argv)
{
   char *cptr;
   FILE *freopen(const char *, const char *, FILE *);
         
 outfil = NULL;
 textout = stdout;
 if (argc==1)                   /* if user needs help */
   run();
 else
 do                                     /* if there are any switches */
  {
   argc--; argv++;
   cptr = *argv;
   if (argc)
    {
     if (*cptr == '-')
      {
       cptr++;
       switch (*cptr)
        {
     
          case 'b': 
                argv++; argc--;
                backgr = int(atof(*argv));
                break;

          case 'f': 
                argv++; argc--;
                filen = *argv;
                break;

          case 'n': 
                argv++; argc--;
                nframes = int(atof(*argv));
                break;

          default:
                fprintf (stderr,"gausnn: unknown switch '%s'\n",*argv);
                exit(2);

        }  /* switch */
      }    /* if */
     else
      {
       if((outfil=fopen(cptr,"w"))==NULL)
         {
           fprintf(stderr,"gausnn: cannot open %s\n",cptr);
           fflush (stderr);
           continue;
         }
       run();
       if (argc <= 1) break;
     }
    }
   else run();
  }
 while (argc > 0);
}

/* -------------------------------------------------------------- */

void run(void)
{

#define CBUFSIZ 200

    int i,j,n;
    int width, height,picsize,bits;
    static int backgr_set = 0;
    FILE *ftemp;
    static char filnam[CBUFSIZ];
    static char sbuf[CBUFSIZ];
    static char cbuf[CBUFSIZ];
    char *parr, *marr;

 parr = marr = NULL; 
 for (n=1; n<=nframes; n++) {
  sprintf (filnam,"%s%04d.ppm", filen,n);
  if (!(ftemp=fopen (filnam,"r"))) {
    fprintf (stderr,"movconvert: can't find file '%s', stopping.\n",filnam);
    exit(2);
  };
  fprintf (stderr,"%s\n", filnam);
  fgets (&sbuf[0],4,ftemp);
  if (!strstr(sbuf,"P6")) {
    fprintf (stderr,"movconvert: incorrect file type '%s'\n",sbuf);
    exit(2);
  }
  fgets (&cbuf[0],CBUFSIZ-1,ftemp);
  fscanf(ftemp,"%d",&width);
  fscanf(ftemp,"%d\n",&height);
  fscanf(ftemp,"%d\n",&bits);

  picsize = height*width*3;
  if (!parr) {
    if (!(parr = (char *)malloc( ((long int) (sizeof(char)*picsize))))) {
      fprintf (stderr,"movconvert: can't allocate buffer space %d\n",
			width*height);
      exit(2);
    }
    if (!(marr = (char *)malloc( ((long int) (sizeof(char)*picsize))))) {
      fprintf (stderr,"movconvert: can't allocate buffer space %d\n",
			width*height);
      exit(2);
    }
  }
  if (!backgr_set) {
    for (i=0; i<picsize; i++) 		/* set background */
      marr[i]   = backgr;
    backgr_set = 1;
  }
  fread (parr,picsize,nw,ftemp);	/* read in picture */
  for (i=0; i<picsize; i+=3) {		/* copy only non-background */
     if ((parr[i]   != backgr) ||
         (parr[i+1] != backgr) ||
         (parr[i+2] != backgr)) {
           marr[i]   = parr[i];
           marr[i+1] = parr[i+1];
           marr[i+2] = parr[i+2];
     }
  }
  //fprintf (stderr,"'%s'", cbuf);
  fclose (ftemp);
  ftemp = fopen (filnam,"w");
  //ftemp = fopen ("test.ppm","w");
  fprintf(ftemp,"P6\n");		/* write header */
  fprintf(ftemp,"%s",cbuf);
  fprintf(ftemp,"%d %d\n",width, height);
  fprintf(ftemp,"%d\n",bits);
  fwrite (marr,picsize,nw,ftemp);	/* write picture */
  fclose (ftemp);
 }
}

/* -------------------------------------------------------------- */
