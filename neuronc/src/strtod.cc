/* module strtod */
/* converts char strings to a floating point value */

#include <stdio.h>

#define isnum(c)  (isdigit(c)||(c)=='+'||(c)=='-'||(c)=='.'||(c)=='e'||(c)=='E')
#define isdigit(c)	('0' <= (c) && (c) <= '9')
#define iswhite(c)	((c) > 0 && ((c) <= ' ' || 0177 <= (c)))

double strtod(const char *s, char **p)

/* Return a double value converted from
   the string s.  Set the string pointer
   p pointing to the character that terminates the scan.
   Skip over optional whitespace at the beginning
   of the string.
*/

{
   char *str;
   double val;

/*   double atof(char *s);  */

   if (!s) return(0.0);
   str = (char *)s; 
   while (iswhite(*str)) str++; 
   if (p) *p = (char *)s;
   if (! isnum(*str)) return(0.0);		/* return if not recognized */
   if (*str== 0) return(0.0);			/* return if at end */

   sscanf (str,"%lf",&val);
 
/*   val = atof(str); 				/* */
   if ((*str=='+')||(*str=='-')) str++;		/* optional '+' or '-' */
   while (isdigit(*str)) str++; 		/* parse the mantissa */
   if (*str=='.') str++;			/* optional '.' */
   while (isdigit(*str)) str++; 		/* parse the fraction */
   if ((*str=='e')||(*str=='E')) {		/* parse the 'e' */
	str++;	
   	if ((*str=='+')||(*str=='-')) str++;	/* optional '+' or '-' */
   	while (isdigit(*str)) str++; 		/* parse the exponent */	
   }
   if (p) *p = str;
   return (val);
}

