/* Segment hdot in program vid */

hdot(x,y)
   int x,y;

/* Write dots on the Hercules compatible graphics card;
    348 * 720 resolution  */


{
 if (x > 719) x = 719;
 if (x < 0)   x = 0;
 if (y > 347) y = 347;
 if (y < 0)   y = 0;


   
/*	The formula for x,y position in graphics page: */

/*  writdot (0x2000 * (y%4) + 90 * (y/4) + x/8, 7 - (x%8) );	*/

/*      A faster way to make formula: */

 writdot (((y&3) << 13) + 90 * (y >> 2) + (x>>3), 7-(x&7)); 

}



