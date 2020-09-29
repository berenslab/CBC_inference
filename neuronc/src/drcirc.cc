drcirc (double x, double y, double rad)
{ 
   double theta, chord;

 move (x,y);
 rpmove (rad,-PI/2);
 if      (rad > .2) chord = .05;
 else if (rad > .1) chord = .1;
 else               chord = .2;
 rpmove (rad*chord/2,-PI);
 for (theta=0.0; theta < PI*2; theta += chord)
   rpdraw (rad*chord,theta);

}


