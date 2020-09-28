

/* GLOBAL DEFINITIONS */

/*  Latest Mod	18-Dec-81		R.G.Smith	*/

/* type defs */

#define FAST register
typedef char TEXT;
typedef short int BOOL;
typedef short int ssiz;		/* link size for incrementing pointers */

/*	macros
 */

#define abs(x)		((x) < 0 ? -(x) : (x))
#define isalpha(c)	(islower(c) || isupper(c))
#define isdigit(c)	('0' <= (c) && (c) <= '9')
#define islower(c)	('a' <= (c) && (c) <= 'z')
#define isupper(c)	('A' <= (c) && (c) <= 'Z')
#define iswhite(c)	((c) <= ' ' || 0177 <= (c)) /* ASCII only */
#define max(x, y)	(((x) < (y)) ? (y) : (x))
#define min(x, y)	(((x) < (y)) ? (x) : (y))
#define tolower(c)	(isupper(c) ? ((c) + ('a' - 'A')) : (c))
#define toupper(c)	(islower(c) ? ((c) - ('a' - 'A')) : (c))

/*	Convenience definitions
*/

#define T 1
#define F 0

/* Logical operators - not bitwise!! */

#define not !
#define and &&
#define or  ||

#define begin {
#define end }

/* Moved to 32 bit 80386 by R.G.Smith, 6/89.
   Idea is to keep all link structures identical
   with their previous counterparts, using  ints,
   and use regular 4-byte (long) ints for everything else.
*/
 
typedef struct 			/* Link structure for all links */
	begin
	  ssiz len,typ;		/* Length of link, class */
	  ssiz section;
	  ssiz area,branch;
	  ssiz xbio,ybio;
	  ssiz linknum;
	  ssiz cellnum;
	  ssiz color;
	  ssiz pts;
	end link;		/* basic link structure */

typedef struct 			/* Link structure for all links */
	begin
	  ssiz len,typ;		/* Length of link, class */
	  ssiz section;
	  ssiz area,branch;
	  ssiz xbio,ybio;
	  ssiz linknum;
	  ssiz cellnum;
	  ssiz color;
	  ssiz nxnum,nxsec;
	  ssiz lanum,lasec;
	  ssiz len2;
	end link2;		/* line link structure */


/* filter defs */
		 
#define linklen  (sizeof(link)/sizeof(ssiz))

/* defs for fpar[][] */

#define SPEC 0			/* Spec number in dpar */
#define PNUM 1			/* pen number */
#define HIDF 2			/* hide flag */
#define DASHN 3			/* dash number */
#define SYMB 4			/* symbol shape and size */
#define STRP 5			/* string pointer for filter spec */
#define FPARWID 6		/* width of fpar[][] */

#define FSIZE 50		/* Filter spec size	*/
#define DFSIZE 200		/* Display filter pointer array size */
#define FGSPCLEN 10		/* Length of flag spec */
#define SYMTABLEN 2000		/* Symbol table length */
#define FLGSIZ 200		/* Flag table length   */

#define NLABELS 50		/* number of labels in picture */

#define VARRSIZ 90		/* display paramater value array */

/* Miscellaneous constants */

#define	MAXFRAME  5000		/* maximum frame for stage */
#define	MAXSEC	  1000		/* Max section number */

#define BTYP 0
#define TTYP 2
#define FTYP 3
#define GTYP 4
#define CTYP 5
#define LTYP 6

#define FORW 1
#define REV  0


/*
 * Variables for file + section buffers
 */

#define BUFSIZE 512			/* File block size */
#define BSIZ   1024			/* Define actual block for compiler */
#define DSIZ   (BSIZ - 6)		/* Define actual data for compiler */
#define BLOCKSIZE   BSIZ * 4		/* Each section has 4K words */
#define DATASIZE    (BLOCKSIZE - 6)	/* Actual data space is less */
#define SECBLKS (2 * BLOCKSIZE / BUFSIZE) /* Number of sector blocks in sect */
#define VERSN	5000			/* "version no." - typ field in file */
#define TXTFIL 1			/* include buffer with file */
#define DATFIL 0			/* no buffer */
#define READ  0
#define WRITE 1
#define RW    2
#define NUMBUFS 3			/* number of section buffers in mont */

typedef	 ssiz sectdata [DSIZ];

typedef	struct 
	begin
	  ssiz len;			/* Length of block */
	  ssiz typ;			/* Type or version num */
	  ssiz section;			/* Section number */
	  ssiz tmpw;			/* Temp word; not used */
	  ssiz free;			/* The free data pointer */
	  ssiz seq;			/* The max sequence til now */
	 sectdata data;			/* The data array for a sect. */
	end 
	sectblock;

typedef struct
	begin
	  int	fd;			/* Process vector file */
	  int	flen;			/* file size        */
	  int	blk;			/* Block pointer    */
	  int	indx;			/* Char index       */
        end pfil;

typedef struct
	begin
	  int	fd;			/* Text file */
	  int	flen;			/* file size        */
	  int	blk;			/* Block pointer    */
	  int	indx;			/* Char index       */
	 TEXT   buf[BUFSIZE];		/* The file buffer  */
        end tfil;


#define MACLEN 80			/* Length of macro command */
#define NMACS  14			/* number of buts and func keys */
#define MESGLEN 180			/* Length of global mesg */
#define NAMSIZE 20			/* Length of name for "mfile" */


/* 	Display modes for screen 	*/

#define DISPXY	1			/* display changed x,y loc */
#define COORDIS 2			/* display only x,y loc */
#define NEWDISP 4			/* re-display same points only */
#define NEWPTS  8			/* re-compute display buffer */

#define CODX -126			/* Codes for trace pt. format */
#define CODY -127
#define CODXY -128

#define BITPSCAL (400./8.)		/* scale factors for bitpad input */
#define MICHAELCONST (BITPSCAL/2.71)

#define FERR  0				/* error return */
#define FEND  1				/* file end */
#define FNORM 2 			/* normal exit */

