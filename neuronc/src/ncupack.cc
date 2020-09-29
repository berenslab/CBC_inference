/* Segment ncunpack in Program nc */

/* Unpacking routines for low level structures */
/* Inputs structures from 

    1) PVM tasks, 
 or
    2) a text file.  */

#include <stdio.h>
#include "pvm3.h"

extern pkfilefl;
extern FILE *filout;

/****************************************************/

int upk_byte (char *cp, int nitem, int stride)
{
    int info;

  if (!pkfilefl) 
    info = pvm_upkbyte (cp,nitem,stride);
  else {
      int i;
   for (i=0; i<nitem; i++,cp+=stride) 
     fscanf (filout,"%d ",cp);
  }
 return info;
}

int upk_short (short *np, int nitem, int stride)
{
    int info;

  if (!pkfilefl) 
    info = pvm_upkshort (np,nitem,stride);
  else {
      int i;
   for (i=0; i<nitem; i++,np+=stride) 
     fscanf (filout,"%d ",np);
  }
 return info;
}

int upk_long (long *np, int nitem, int stride)
{
    int info;

  if (!pkfilefl) 
    info = pvm_upklong (np,nitem,stride);
  else {
      int i;
   for (i=0; i<nitem; i++,np+=stride) 
     fscanf (filout,"%d ",np);
  }
 return info;
}

int upk_int (int *np, int nitem, int stride)
{
    int info;

  if (!pkfilefl) 
    info = pvm_upkint (np,nitem,stride);
  else {
      int i;
   for (i=0; i<nitem; i++,np+=stride) 
     fscanf (filout,"%d ",np);
  }
 return info;
}

int upk_float(float *fp, int nitem, int stride)
{
    int info;

  if (!pkfilefl) 
    info = pvm_upkfloat(fp, nitem, stride);
  else {
      int i;
   for (i=0; i<nitem; i++,fp+=stride) 
     fscanf (filout,"%.8g ",fp);
  }
 return info;
}

int upk_double(double *dp, int nitem, int stride)
{
    int info;

  if (!pkfilefl) 
    info = pvm_upkdouble(dp, nitem, stride);
  else {
      int i;
   for (i=0; i<nitem; i++,dp+=stride) 
     fscanf (filout,"%.16g ",dp);
  }
 return info;
}

int upk_str(char *sp)
{
  int info;

  if (!pkfilefl) 
    info = pvm_upkstr(sp);
  else 
     fscanf (filout,"%s ",sp);
 return info;
}

