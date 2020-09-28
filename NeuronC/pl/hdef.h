/* Segment HDEF.H in program MPICT */

/* Global defs for hidden line table */

/* Hidden area is represented by entries in
hidtab, which are of type "hvect".
The entries represent the hidden area
polygon which is defined by a sequence
of points in a linked list. Each point in the 
list has a pointer to the next point.

   New segments are represented by entries in
"stab", of type "svect".  New segments start
with a new point, and may have a sequence of new
points followed by an intercept for the last vector
in the segment.
   Each entry in "stab" has a pointer to the first new point,
and a pointer to the last point (the intercept) in the segment.

*/

typedef struct
	 begin
	  struct hvect *vlnk;	/* link to next vector */
	  int snum;		/* section number for delete area */
	  long int xv;		/* x,y values for point */
	  long int yv;
	 end hvect;

typedef struct
	 begin
	  struct hvect *vlnk;	/* start of segment in hidtab */
	  struct hvect *intp;	/* last point (intercept) for this seg */
	 end svect;

#define XINF 2147483647

/* #define XINF 32767 */

#define HSIZE 800
#define ISIZE 200
#define SSIZE 200

#define INTERCEPT 1111 


