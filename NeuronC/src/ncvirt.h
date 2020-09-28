/* Header file ncvirt.h for program "nc". */
/* Structures for virtual compartments. */

struct virtobj {                   	/* structure for virtual list */
        short int structype;           	/* type of structure: 1=>compartment */
        int num;			/* unique id number */
	int tid;			/* task id # (object's location) */
	virtobj *next;			/* pointer to next object */
	virtobj *last;
       };

struct updatobj {                 	/* struct for comm betw. virt comps */
	comp *addr;			/* address of compartment */
        short int structype;           	/* type of structure: 1=>compartment */
	short int field;		/* field that needs updating in comp */
	int compnum;			/* unique number for comp. */
	int tid;			/* task where structure is located */
	updatobj *next;
	updatobj *last;
       };

#define VEND    0		/* end of structure */ 
#define VCOMP   1		/* type of structure for virt: compartment */
#define VCONLST 2		/* conlst */ 
#define VCONN   3		/* connection */ 
#define VRECEP  4		/* photoreceptor */ 
#define VCACOMP 5		/* calcium comp */ 
#define VSYNAP  6		/* synapse */ 
#define VLOAD   7		/* load */ 
#define VBUF    8		/* delay buffer */ 
#define VCHAN   9		/* membrane channel */ 

#define VVOLT 4			/* volt field of struct comp */


