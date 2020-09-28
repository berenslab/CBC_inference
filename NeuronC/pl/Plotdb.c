#
/*  A program to list the intermediate plot file
 *  in a readable form.
 *
 *  Commands c,h,i,j added for plotter and hidden line routine.
 *    R.G.Smith   7-Oct-83
 *  Commands a,g added for graphics mode
 *    R.G.Smith   1-July-86
 */

#include        <stdio.h>
#define SCALE   16384.0
#define PLTCMD  0200
#define CMDMSK  0177
#define ENDRAS  -1
#define ESC     033
#define ESCMAP  '\\'
char inbuf[BUFSIZ];
extern int getwsx();

FILE *pltin;

void doplot();

int main(int argc, char **argv)

   {
        register char *cptr;
        FILE *temp, *freopen();

        pltin= stdin;
        setbuf(pltin,inbuf);
        do
           {
                argc--; argv++;
                cptr= *argv;
                if(argc)
                        if((temp=freopen(cptr,"r",pltin))==NULL)
                           {
                                printf("Plotdb: cannot open %s\n",cptr);
                                continue;
                           }
                          else pltin= temp;
                doplot();
           } while( argc );
   }

void doplot(void)
   {
        int i,size,orient,fat,dotsize,fill;
        char tfont;
        float x,y,r,xmin,xmax,ymin,ymax,angle,tsize;
        int nlines,nesc,offset,count,header,omin,omax,cmin,cmax;
        register int c;
        while((c=getc(pltin)) != EOF)
           {
                if ((c&PLTCMD) == 0 || ((c&CMDMSK) < 'a' || (c&CMDMSK) > 'z'))
                   {
                        printf("straight text:\n");
                        do
                          {     if (c&PLTCMD)
                                        putchar(ESCMAP);
                                putchar(c&CMDMSK);
                          }
                        while(((c=getc(pltin))&PLTCMD)==0);
                        if( c == EOF ) break;
                        ungetc(c,pltin);
                        continue;
                   }
                switch(c &= CMDMSK)
                  {
                  case 'e':
                  case 'b':
                        for(i=0;i<10;i++) printf("**");
                        printf("\n%c:\n",c);
                        break;
                  case 'm':
                  case 'd':
                        x=getwsx(pltin)/SCALE;
                        y=getwsx(pltin)/SCALE;
                        printf("%c: %-9.6f %-9.6f\n",c,x,y);
                        break;
                  case 'k':
                        r=getwsx(pltin)/SCALE;
                        fill= getc(pltin)&0377;
                        printf("%c: circle radius %-9.6g fill %d\n", c,r,fill);
			break;
                  case 'l': 		/* rectangle */
			{
			double x1,y1,x2,y2,x3,y3,x4,y4,temp;

                        fill= getc(pltin)&0377;
                        temp= getc(pltin)&0377;	/* "m" */
                        x1=getwsx(pltin)/SCALE;
                        y1=getwsx(pltin)/SCALE;
                        temp= getc(pltin)&0377;	/* "d" */
                        x2=getwsx(pltin)/SCALE;
                        y2=getwsx(pltin)/SCALE;
                        temp= getc(pltin)&0377;	/* "d" */
                        x3=getwsx(pltin)/SCALE;
                        y3=getwsx(pltin)/SCALE;
                        temp= getc(pltin)&0377;	/* "d" */
                        x4=getwsx(pltin)/SCALE;
                        y4=getwsx(pltin)/SCALE;
                        printf("%c: rectangle x1 %-4g y1 %-4g x2 %-4g y2 %-4g x3 %-4g y3 %-4g x4 %-4g y4 %-4g fill %d\n",
			 c,x1,y1,x2,y2,x3,y3,x4,y4,fill);
			}
			break;
                  case 'o': 		/* fill */
                        fill= getc(pltin)&0377;
                        printf("%c: fill %d\n",c,fill);
			break;
                  case 'v': 		/* triangle */
			{
			double x1,y1,x2,y2,x3,y3,temp;

                        fill= getc(pltin)&0377;
                        temp= getc(pltin)&0377;	/* "m" */
                        x1=getwsx(pltin)/SCALE;
                        y1=getwsx(pltin)/SCALE;
                        temp= getc(pltin)&0377;	/* "d" */
                        x2=getwsx(pltin)/SCALE;
                        y2=getwsx(pltin)/SCALE;
                        temp= getc(pltin)&0377;	/* "d" */
                        x3=getwsx(pltin)/SCALE;
                        y3=getwsx(pltin)/SCALE;
                        printf("%c: triangle x1 %-4g y1 %-4g x2 %-4g y2 %-4g x3 %-4g y3 %-4g fill %d\n",
			 c,x1,y1,x2,y2,x3,y3,fill);
			}
			break;
                  case 's':
                        tsize= (getc(pltin)&0377)*64./SCALE;
                        tfont= getc(pltin)&0377;
                        if (tfont == '\0') tfont = 'A';
                        angle = getwsx(pltin)*360./32768.;
                        printf("s: font=%c, width=%f, angle=%f degrees\n",
                                tfont, tsize, angle);
                        break;
                  case 't':
                        printf("t: ");
                        while((c=getc(pltin)))
                          {     if(c==ESC)
                                  {     switch(c=getc(pltin))
                                          {
                                          case 'u':
                                          case 'U':
                                          case 'd':
                                          case 'D':
                                          case 'b':
                                          case 's':
                                          case 'S':
                                          case 'f':
                                          case 'F':
                                          case 'g':
                                                putchar(ESCMAP);
                                                putchar(c);
                                                break;
                                          default:
                                                putchar(ESCMAP);
                                                putchar(c);
                                                putchar('?');
                                                break;
                                          }
                                  }
                                else    printf("%c",c);
                          }
                        printf("\n");
                        break;
                  case 'f':
                        fat=getwsx(pltin);
                        printf("f: %d\n",fat);
                        break;
                  case 'w':
                        xmin=getwsx(pltin)/SCALE;
                        xmax=getwsx(pltin)/SCALE;
                        ymin=getwsx(pltin)/SCALE;
                        ymax=getwsx(pltin)/SCALE;
                        printf("w: %-9.6f %-9.6f %-9.6f %-9.6f\n",
					xmin,xmax,ymin,ymax);
                        break;
                  case 'r':
                        dotsize=getc(pltin);
                        nlines=0;
                        nesc = 0;
                        cmax= 0;
                        cmin= 10000;
                        omax= 0;
                        omin= 10000;
                        while((header=getwsx(pltin)) != ENDRAS)
                          {
                                count= header&0377;
                                offset= (header>>8)&0377;
                                if (count == 0377)      /* escape */
                                  {     switch (offset)
                                          {
                                          case 't':
                                          case 'T':
                                                for (count=0; getc(pltin)>0; )
                                                        count++;
                                                if (count&1)
                                                        getc(pltin);
                                                break;
                                          }
                                        nesc++;
                                  }
                                else            /* raster line */
                                  {     if(cmin>count) cmin=count;
                                        if(cmax<count) cmax=count;
                                        if(omin>offset) omin=offset;
                                        if(omax<offset) omax=offset;
                                        while(count--) getwsx(pltin);
                                        nlines++;
                                  }
                           }
                        printf("r: %d pels/dot, %d lines, %d escapes,\n",
                                dotsize, nlines, nesc);
                        printf("offset(min,max)= (%d,%d), ",omin,omax);
                        printf("count(min,max)= (%d,%d)\n",cmin,cmax);
                        break;
                  case 'x':
                        printf("display X-hairs\n");
                        break;
                  case 'p':
                        printf("purge plot buffers\n");
                        break;
                  case 'n':
                        printf("n: new page\n");
                        break;
		  case 'c':
			printf ("c: change pen color to %d\n",getwsx(pltin));
                        break;
		  case 'h':
			printf ("h: remove hidden lines\n");
                        break;
		  case 'i':
			printf ("i: initialize hidden line table\n");
                        break;
		  case 'j':
			printf ("j: no hidden lines\n");
			break;
		  case 'g':
			printf ("g: graphics mode\n");
			break;
		  case 'a':
			printf ("a: text mode\n");
			break;
                  default:
                        printf("unknown command %c %o\n",c,c);
           }
   }
   }
