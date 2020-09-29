/* procedure to connect cones with gap junctions */

proc congj(xi,xj,xh,xk) {

source = conarr[xi][xj][0];		/* source cone */
sz     = conarr[xi][xj][1];

dest   = conarr[xh][xk][0];		/* dest cone */
dz     = conarr[xh][xk][1];

found=0;
 for (zi=0; zi<nconn; zi++) {		/* check for previous connection */
   if (gconnect[source][zi]==dest) found = 1;
   if (gconnect[source][zi]<0) break;
 };

 if (found == 0) {
    gconnect[source][zi] = dest;		/* remember connection */
    conn [source][sz] to [dest][dz] gj gjcond;
/*  snode = conarr[xi][xj][2]++;
    dnode = conarr[xh][xk][2]++;
    conn [source][sz]    to [source][snode] cable dia .3 length 3;
    conn [source][snode] to [dest][dnode]   gj gjcond; 
    conn [dest][dnode]   to [dest][dz]      cable dia .3 length 3;
*/
 };
 for (zi=0; zi<nconn; zi++) {		/* remember reverse conn */
   if (gconnect[dest][zi]==source) found = 1;
   if (gconnect[dest][zi]<0) break;
 };
 if (found == 0) 
   gconnect[dest][zi] = source;
};


