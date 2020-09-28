
#define abs(x) ((x) < 0 ? -(x) : (x))


double
modf(double value, double *iptr)
                       

/* Returns the positive fractional part of value
   and stores the integer part indirectly through
   iptr.
*/

{
   double ipart, fpart;

   ipart = (long int) value;
   fpart = value - ipart;
   if (fpart < 0)
     {
       fpart += 1.0;
       ipart -= 1.0;
     }
   *iptr = ipart;
   return (fpart);
}

/*------------------------------------------*/

frexp(double value, int *eptr)
               
            

/* Returns the mantissa of a double value as
   a double quantity, x, of magnitude less
   than 1 and stores an integer n such that
   value = x*2**n indirectly through eptr.
*/

{
    int i;

  for (i=0; abs(value) >= 1.0; i++)
     value /= 2.0;
  *eptr = i;
  return (value); 
}
