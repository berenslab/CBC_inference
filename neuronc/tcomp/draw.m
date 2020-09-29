/* draw the contents of a file */

fread con1, cone1;


proc dcone(xloc,yloc,scale) {

 gmove (xloc,yloc);

 for (i=0,x=1; x>= -1e6; i++) {
   x = cone1[i][0];
   y = cone1[i][1];
   if (x>=0) gdraw (x*scale+xloc,y*scale+yloc);
 };
};

dcone(.5,.5,.05);

