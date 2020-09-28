
dim ans[20];

proc hex(dec) {

for (i=0; i<20; i++) {
  ans[i] = 0;
};

x = dec;
if (x < 0) 
   x += 65536;

for (i=0; i<20; i++) {

if (x <= 0) break;
y = int(x/16);
z = x - (y * 16);
ans[i] = z;
x = y;
};

for (i=19; i>=0; i--) {
  val = ans[i];
  if (val > 0) break; 
};

base = 1;
for (i; i>=0; i--) {
  val = ans[i];
  if (val < 10)  
      printf ("%g",val);
  if (val >= 10) {
     if (val == 10) printf ("A"); 
     if (val == 11) printf ("B"); 
     if (val == 12) printf ("C"); 
     if (val == 13) printf ("D"); 
     if (val == 14) printf ("E"); 
     if (val == 15) printf ("F"); 
  };

};

 printf ("\n");
};

