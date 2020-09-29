/* module pcrf.m */

/* Plots cone rf with legend */

colsep  = .1;

/*-----------------------------------------------*/

proc ptext (pen,x,y,text) {

   gpen (pen);
   gmove (x,y);
   gtext (text);
};

/*-----------------------------------------------*/

proc pleg (pen,x1,y,text) {

   gpen (pen);
   gdash (pen-1);
   gmove (x1,y);
   grdraw (colsep-.025,0.);
   gmove (x1+colsep, y);
   gtext (text);
};

/*-----------------------------------------------*/

lcol    = .55;
linsp   = .04;
toprow = .95;
dim r[12];

proc drawlgndn (file,hdr,all,alg,L1,L2,L3,L4,L5,L6) {

/* draw normalized legend */

rcol    = lcol + colsep;

for (i=0; i<12; i++) {
   r[i] = toprow - i * linsp;
};

ptext (7, lcol, r[0],  file);

ptext (7, lcol, r[2],  hdr);
ptext (7, lcol, r[3],  all);
ptext (7, rcol, r[3],  alg);
if (L1) pleg  (1, lcol, r[4],  L1);
if (L2) pleg  (2, lcol, r[5],  L2);
if (L3) pleg  (3, lcol, r[6],  L3);
if (L4) pleg  (4, lcol, r[7],  L4);
if (L5) pleg  (5, lcol, r[8],  L5);
if (L6) pleg  (6, lcol, r[9],  L6);

/*
ptext (7, 0.52, 0.02, "microns");
gcrot (90);
ptext (7, 0.08, 0.45, "Normalized Response");
*/

erase r;
};

/*-----------------------------------------------*/

proc drawlgnda (file,hdr,all,alg,L1,L2,L3,L4,L5,L6) {

/* draw absolute legend */

rcol    = lcol + colsep;

for (i=0; i<12; i++) {
   r[i] = toprow - i * linsp;
};

ptext (7, lcol, r[0],  file);

ptext (7, lcol, r[2],  hdr);
ptext (7, lcol, r[3],  all);
ptext (7, rcol, r[3],  alg);
if (L1) pleg  (1, lcol, r[4],  L1);
if (L2) pleg  (2, lcol, r[5],  L2);
if (L3) pleg  (3, lcol, r[6],  L3);
if (L4) pleg  (4, lcol, r[7],  L4);
if (L5) pleg  (5, lcol, r[8],  L5);
if (L6) pleg  (6, lcol, r[9],  L6);

ptext (7, 0.52, 0.02, "microns");
gcrot (90);
ptext (7, 0.06, 0.48, "Absolute Response  (mV)");

erase r;
};

/*-----------------------------------------------*/

proc drawlgnda2 (file,hdr,all,alg,L1,L2,L3,L4,L5,L6) {

/* draw absolute legend */
/* second label is pen 1 */

rcol    = lcol + colsep;

for (i=0; i<12; i++) {
   r[i] = toprow - i * linsp;
};

ptext (7, lcol, r[0],  file);

ptext (7, lcol, r[2],  hdr);
ptext (7, lcol, r[3],  all);
ptext (7, rcol, r[3],  alg);
if (L1) pleg  (2, lcol, r[4],  L1);
if (L2) pleg  (1, lcol, r[5],  L2);
if (L3) pleg  (3, lcol, r[6],  L3);
if (L4) pleg  (4, lcol, r[7],  L4);
if (L5) pleg  (5, lcol, r[8],  L5);
if (L6) pleg  (6, lcol, r[9],  L6);

ptext (7, 0.52, 0.02, "microns");
gcrot (90);
ptext (7, 0.06, 0.48, "Absolute Response  (mV)");

erase r;
};

/*-----------------------------------------------*/

