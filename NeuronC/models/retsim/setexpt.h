/* setexpt.h */
/* includes for setexpt.cc */

/* functions in the experiment file, to be linked with retsim */

extern void (*defparams)(void);	/* Define params for the experiment */
                   		/*  defined in the experiment file */

extern void (*setparams)(void);	/* Set or possibly change the neuron params */
                        	/*  defined in the experiment file */

extern void (*setdens)(void);	/* Set or possibly change the density params */
                        	/*  defined in the experiment file */

extern void (*addcells)(void);	/* Add extra cells  */
                        	/*  defined in the experiment file */

extern void (*addlabels)(void);	/* Add node labels  */
                        	/*  defined in the experiment file */

extern void (*runexpt)(void);		/* run the experiment */

