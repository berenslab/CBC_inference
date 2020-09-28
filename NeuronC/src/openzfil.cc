
/* Open file with checking for .gz, .bz2, 
   either explicit or implicit */

/* R.G.Smith Feb 2009 */


#include <stdio.h>
#include <string.h>

#define FILENAMLEN 200

static int info=0;

/*---------------------------------------------------------*/

FILE *openfilz (char *filnam, char *infilm, const char *fm)

/* Open file for reading.  Check if file is .gz or .bz2, 
   or has a .gz or .bz2 version in same directory.  If so,
   open pipe to read the file with zcat or bzcat. Store filename
   actually used in '*infilm'.
*/

{
   static char infilz[FILENAMLEN];
   FILE *fildesc;
   const char *progname="openfilz";

     /* check if filename is .gz */

   if (info>=3) fprintf (stderr,"# opening file '%s'\n",filnam); 

   strcpy(infilm,filnam);

   if (strstr(filnam,".gz") != NULL) {
       strcpy(infilm,filnam);
       strcpy(infilz,"zcat ");
       strncat(infilz,filnam,FILENAMLEN-5);
       if ((fildesc=popen(infilz,fm))==NULL) {
           fprintf (stderr,"%s: can't open file '%s'\n",progname,filnam);
          fildesc = NULL;
       }
   }
      /* check if filename is .bz2 */

   else if (strstr(filnam,".bz2") != NULL) {
       strcpy(infilm,filnam);
       strcpy(infilz,"bzcat ");
       strncat(infilz,filnam,FILENAMLEN-6);
       if ((fildesc=popen(infilz,fm))==NULL) {
          fprintf (stderr,"%s: can't open file '%s'\n",progname,filnam);
          fildesc = NULL;
       }
   }
      /* check if filename is .xz */

   else if (strstr(filnam,".xz") != NULL) {
       strcpy(infilm,filnam);
       strcpy(infilz,"xzcat ");
       strncat(infilz,filnam,FILENAMLEN-6);
       if ((fildesc=popen(infilz,fm))==NULL) {
          fprintf (stderr,"%s: can't open file '%s'\n",progname,filnam);
          fildesc = NULL;
       }
   }
      /* check if filename is .zip */

   else if (strstr(filnam,".zip") != NULL) {
       strcpy(infilm,filnam);
       strcpy(infilz,"unzip -p ");
       strncat(infilz,filnam,FILENAMLEN-9);
       if ((fildesc=popen(infilz,fm))==NULL) {
          fprintf (stderr,"%s: can't open file '%s'\n",progname,filnam);
          fildesc = NULL;
       }
   }
      /* Try to open user-specified non .gz file, if it doesn't then check if .gz version exists */
 
   else if ((fildesc=fopen(filnam,fm)) == NULL) {
      strcpy(infilz,filnam);
      strncat(infilz,".gz",FILENAMLEN-5);
      if ((fildesc=fopen(infilz,fm))!=NULL){ /* look for .gz file */
         if (info>=3) fprintf (stderr,"# trying  file '%s'\n",infilz); 
         fclose(fildesc);
         strcpy(infilm,infilz);
         strcpy(infilz,"zcat ");          /* try to uncompress .gz file */
         strncat(infilz,filnam,FILENAMLEN-5);
         strncat(infilz,".gz",FILENAMLEN-5);
         if ((fildesc=popen(infilz,fm))==NULL) {
            fprintf (stderr,"%s: can't open file '%s'\n", progname,filnam);
            fildesc = NULL;
         }
      }

      else {   /* else check if .bz2 version exists */
        strcpy(infilz,filnam);
        strncat(infilz,".bz2",FILENAMLEN-5);
        if ((fildesc=fopen(infilz,fm))!=NULL){ /* look for .bz2 file */
           if (info>=3) fprintf (stderr,"# trying  file '%s'\n",infilz); 
           fclose(fildesc);
           strcpy(infilm,infilz);
           strcpy(infilz,"bzcat ");          /* try to uncompress .bz2 file */
           strncat(infilz,filnam,FILENAMLEN-5);
           strncat(infilz,".bz2",FILENAMLEN-5);
           if ((fildesc=popen(infilz,fm))==NULL) {
              fprintf (stderr,"%s: can't open file '%s'\n", progname,filnam);
              fildesc = NULL;
           }
        }

        else {   /* else check if .xz version exists */
        strcpy(infilz,filnam);
        strncat(infilz,".xz",FILENAMLEN-5);
        if ((fildesc=fopen(infilz,fm))!=NULL){ /* look for .xz file */
           if (info>=3) fprintf (stderr,"# trying  file '%s'\n",infilz); 
           fclose(fildesc);
           strcpy(infilm,infilz);
           strcpy(infilz,"xzcat ");          /* try to uncompress .xz file */
           strncat(infilz,filnam,FILENAMLEN-6);
           strncat(infilz,".xz",FILENAMLEN-9);
           if ((fildesc=popen(infilz,fm))==NULL) {
              fprintf (stderr,"%s: can't open file '%s'\n", progname,filnam);
              fildesc = NULL;
           }
        }

        else {  /* else check if .zip version exists */
        strcpy(infilz,filnam);
        strncat(infilz,".zip",FILENAMLEN-5);
        if ((fildesc=fopen(infilz,fm))!=NULL){ /* look for .zip file */
           if (info>=3) fprintf (stderr,"# trying  file '%s'\n",infilz); 
           fclose(fildesc);
           strcpy(infilm,infilz);
           strcpy(infilz,"unzip -p ");          /* try to uncompress .zip file */
           strncat(infilz,filnam,FILENAMLEN-5);
           strncat(infilz,".zip",FILENAMLEN-5);
           if ((fildesc=popen(infilz,fm))==NULL) {
              fprintf (stderr,"%s: can't open file '%s'\n", progname,filnam);
              fildesc = NULL;
           }
        } else {
          fprintf (stderr,"%s: can't find file '%s'\n",progname,filnam);
          fildesc=NULL;
        }
       }
      }
     }
   }    /* if (fopen(filnam,fm) == NULL) */

   //if (info>=3) fprintf (stderr,"# saving filename '%s' for rewind.\n",infilm); 

   return fildesc;
}

