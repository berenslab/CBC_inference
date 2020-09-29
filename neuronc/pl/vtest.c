

#include <stdio.h>

int scrntyp;
int scrcolor=5;
int scrrev;

main()

{
  gmodev();
  movev(1,1); drawv(16000,16000); 
  drawv(16000,1); drawv(1,1); 
  scrcolor = 3;
  movev(2000,2000); drawv (10000,2000); 
  getc(stdin);
  clrcv(); 
  getc(stdin);
  tmodev(); 
return;


}

tdot()
{}
setorg()
{}

