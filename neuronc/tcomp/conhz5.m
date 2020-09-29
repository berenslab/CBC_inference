/* procedure to connect horizontal cells with gap junctions */

proc conhz(xi,xj,xh,xk) {

source = hzarr[hztype][xi][xj][0];		/* source hz */
sn     = hzarr[hztype][xi][xj][1];		/* source hz node */

dest   = hzarr[hztype][xh][xk][0];		/* dest hz */
dn     = hzarr[hztype][xh][xk][1];		/* dest hz node */

found=0;
 for (zi=0; zi<nhconn; zi++) {		/* check for previous connection */
   if (hzconn[hztype][source][zi]==dest) found = 1;
   if (hzconn[hztype][source][zi]<0) break;
 };

 if (found == 0) {
/*    conn [sn][soma] to [sn][1000+dn] cable length=50 dia=1; /* make cable */
/*    conn [sn][1000+dn] to [dn][1000+sn] gj hgjcond;/* make connection */
 /*   conn [dn][soma] to [dn][1000+sn] cable length=50 dia=1; /* make cable */
    conn [sn][soma] to [dn][soma] gj hgjcond;/* make connection */
    hzconn[hztype][source][zi] = dest;	        /* remember connection */
 };

 for (zi=0; zi<nhconn; zi++) {		/* check for reverse connection */
   if (hzconn[hztype][dest][zi]==source) found = 1;
   if (hzconn[hztype][dest][zi]<0) break;
 };
 if (found == 0) 			/* remember reverse connection */
   hzconn[hztype][dest][zi] = source;
};


