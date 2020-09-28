/* segment ncv in program nc */

/* Print version routine */

#include "nc.h"
#include "y.tab.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#ifdef CPML
#include <cpml.h>
#else
#include <math.h>
#endif

#ifdef __cplusplus
}
#endif

extern int ncversion[];

/*------------------------------------*/

char *print_version (double version)

{
     int v1, v2, v3;    
     static char vbuf[40];

    v1 = int(version);
    v2 = int((version-v1) * 10 + 0.0001);
    v3 = int((((version-v1) * 10) - v2) * 100 + 0.0001);
    if (v3 > 0)
      sprintf (vbuf,"%d.%d.%d",v1,v2,v3); 
    else if (v2 > 0)
      sprintf (vbuf,"%d.%d",v1,v2); 
    else if (v1 > 0)
      sprintf (vbuf,"%d",v1); 
  return vbuf;
}

/*------------------------------------*/

datum print_ncversion (datum &x)

{
     datum d={0};

  d.str = print_version (x.val);
  d.vtype = STRING;
  return d;
}

/*------------------------------------*/

double scanversion (char *vbuf)

/* Read in a version number from a text file. */

{
    int v0,v1,v2;
    double version;

  sscanf (vbuf,"%d.%d.%d",&v0,&v1,&v2);
  version = v0 + v1*1e-1 + v2*1e-3;
  return version;
}
