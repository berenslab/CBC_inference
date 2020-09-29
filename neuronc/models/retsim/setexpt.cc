
/* module setexpt in program nc */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dlfcn.h>

#include "ncio.h"
#include "ncfuncs.h"
#include "retsim.h"
#include "retsim_var.h"

void nullfunc(void) {}

void (*defparams)(void)=nullfunc;	/* Define params for the experiment */
                                	/*  defined in the experiment file */

void (*setparams)(void)=nullfunc;	/* Set or possibly change the neuron params */
                                	/*  defined in the experiment file */

void (*setdens)(void)=nullfunc;		/* Set or possibly change the density params */
                                	/*  defined in the experiment file */

void (*addcells)(void)=nullfunc;	/* Add extra cells */
                                	/*  defined in the experiment file */

void (*addlabels)(void)=nullfunc;	/* Add node labels */
                                	/*  defined in the experiment file */

void (*runexpt)(void)=nullfunc;		/* run the experiment */

/*----------------------------------------------------------------*/

void setexpt(void)

/* Simple plugin definition for nc experiments */

{
#define XBUFSIZ 80
    int i;
    void *lib_handle;
    char *error;
    char exptlib[XBUFSIZ];
    void (*tempfunc)(void);

  setexptvar(); 			/* set "expt" command line variable */

  if (strcmp(expt,"help")==0) return;

  strncpy (exptlib,"./expt_",XBUFSIZ);
  strncat (exptlib,expt,XBUFSIZ);
// #ifdef MACOSX
//   strncat (exptlib,".dylib",XBUFSIZ);
// #else
  strncat (exptlib,".so",XBUFSIZ);
// #endif

  lib_handle = dlopen(exptlib, RTLD_LAZY);
  if (!lib_handle) { ncfprintf(stderr, "setexpt: %s\n", dlerror()); exit(1); }

  tempfunc = (void (*)(void)) dlsym(lib_handle, "_Z9defparamsv");
  if ((error = dlerror()) == NULL)  { defparams = tempfunc; }

  tempfunc = (void (*)(void)) dlsym(lib_handle, "_Z9setparamsv");
   if ((error = dlerror()) == NULL)  { setparams = tempfunc; }

  tempfunc = (void (*)(void)) dlsym(lib_handle, "_Z7setdensv");
   if ((error = dlerror()) == NULL)  { setdens = tempfunc; }

  tempfunc = (void (*)(void)) dlsym(lib_handle, "_Z8addcellsv");
   if ((error = dlerror()) == NULL)  { addcells = tempfunc; }

  tempfunc = (void (*)(void)) dlsym(lib_handle, "_Z9addlabelsv");
   if ((error = dlerror()) == NULL)  { addlabels = tempfunc; }

  runexpt   = (void (*)(void)) dlsym(lib_handle, "_Z7runexptv");
   if ((error = dlerror()) != NULL)  { ncfprintf(stderr, "%s\n", error); exit(1); }

  // dlclose(lib_handle);

#undef XBUFSIZ 
}


