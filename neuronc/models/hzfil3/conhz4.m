/* procedure to connect horizontal cells with gap junctions */

proc conhz(xi,xj,xh,xk) {

source = hzarr[xi][xj][0];		/* source hz */
sn     = hzarr[xi][xj][1];		/* source hz node */

dest   = hzarr[xh][xk][0];		/* dest hz */
dn     = hzarr[xh][xk][1];		/* dest hz node */

found=0;
 for (zi=0; zi<nhconn; zi++) {		/* check for previous connection */
   if (hzconn[source][zi]==dest) found = 1;
   if (hzconn[source][zi]<0) break;
 };

 if (found == 0) {
/*    conn [sn][soma] to [sn][1000+dn] cable length=50 dia=1; /* make cable */
/*    conn [sn][1000+dn] to [dn][1000+sn] gj hgjcond;/* make connection */
 /*   conn [dn][soma] to [dn][1000+sn] cable length=50 dia=1; /* make cable */
    conn [sn][soma] to [dn][soma] gj hgjcond;/* make connection */
    hzconn[source][zi] = dest;		        /* remember connection */
 };

 for (zi=0; zi<nhconn; zi++) {		/* check for reverse connection */
   if (hzconn[dest][zi]==source) found = 1;
   if (hzconn[dest][zi]<0) break;
 };
 if (found == 0) 			/* remember reverse connection */
   hzconn[dest][zi] = source;
};


