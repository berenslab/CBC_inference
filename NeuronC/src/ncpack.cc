/* Segment ncpack in Program nc */

/* Packing routines for low level structures */
/* Outputs structures to 

    1) PVM tasks, 
 or
    2) to a text file.  */

#include <stdio.h>
#include "pvm3.h"

extern pkfilefl;
extern FILE *filout;

/****************************************************/

int pk_byte (char *cp, int nitem, int stride)
{
    int info;

  if (!pkfilefl) 
    info = pvm_pkbyte (cp,nitem,stride);
  else {
      int i;
   for (i=0; i<nitem; i++,cp+=stride) 
     ncfprintf (filout,"%d ",*cp);
  }
 return info;
}

int pk_short (short *np, int nitem, int stride)
{
    int info;

  if (!pkfilefl) 
    info = pvm_pkshort (np,nitem,stride);
  else {
      int i;
   for (i=0; i<nitem; i++,np+=stride) 
     ncfprintf (filout,"%d ",*np);
  }
 return info;
}

int pk_long (long *np, int nitem, int stride)
{
    int info;

  if (!pkfilefl) 
    info = pvm_pklong (np,nitem,stride);
  else {
      int i;
   for (i=0; i<nitem; i++,np+=stride) 
     ncfprintf (filout,"%d ",*np);
  }
 return info;
}

int pk_int (int *np, int nitem, int stride)
{
    int info;

  if (!pkfilefl) 
    info = pvm_pkint (np,nitem,stride);
  else {
      int i;
   for (i=0; i<nitem; i++,np+=stride) 
     ncfprintf (filout,"%d ",*np);
  }
 return info;
}

int pk_float(float *fp, int nitem, int stride)
{
    int info;

  if (!pkfilefl) 
    info = pvm_pkfloat(fp, nitem, stride);
  else {
      int i;
   for (i=0; i<nitem; i++,fp+=stride) 
     ncfprintf (filout,"%.8g ",*fp);
  }
 return info;
}

int pk_double(double *dp, int nitem, int stride)
{
    int info;

  if (!pkfilefl) 
    info = pvm_pkdouble(dp, nitem, stride);
  else {
      int i;
   for (i=0; i<nitem; i++,dp+=stride) 
     ncfprintf (filout,"%.16g ",*dp);
  }
 return info;
}

int pk_str(char *sp)
{
  int info;

  if (!pkfilefl) 
    info = pvm_pkstr(sp);
  else 
     ncfprintf (filout,"%s ",sp);
 return info;
}

