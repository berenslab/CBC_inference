#include <stdio.h>
#include "stdplt.h"

#define NGRAPHS 8

FILE	*fgraph;
int	i,npoints[NGRAPHS],ngraphs;
char	label[50]	="";
char	Label[50]	="";
char	min[20]	="";
char	max[20]	="";
char	pen[20]	="";
char	tic[20]	="";
char	Tic[20]	="";
char	siz[20]	="";
char	Siz[20]	="";
char	Dsiz[10]="";
char    logc[10] ="";
char    Logc[10] ="";
char	grid[50]	="%*G 3 ";
char	points[50]	="";
char	line[120];
char	cmd[500];
float	x[1000],y[1000];
FILE	*fopen();
int np4,np5,np6;
float *argvals[20];		/* array to hold arguments passed to graph() */

main(argc,argv)
char	*argv[];
{
    int i;

	fgraph=stdin;
	argc--;
	argv++;
	while(argc--){
		if(**argv=='-')
			switch(argv[0][1]){
				case 'f':
					argc--;
					if(**(++argv)=='-'){
						fgraph=stdin;
						readin();
						break;
					}
					fgraph=fopen(*argv,"r");
					readin();
					break;
				case 'l':
					argc--;
					sprintf(label,"%%*l %s ",*(++argv));
					break;
				case 'L':
					argc--;
					sprintf(Label,"%%*l %s ",*(++argv));
					break;
				case 'm':
					argc--;
					sprintf(min,"%%*m %s ",*(++argv));
					break;
				case 'M':
					argc--;
					sprintf(max,"%%*M %s ",*(++argv));
					break;
				case 'p':
					argc--;
					sprintf(pen,"%%*p %s ",*(++argv));
					break;
				case 't':
					argc--;
					sprintf(tic,"%%*t %s ",*(++argv));
					break;
				case 'T':
					argc--;
					sprintf(Tic,"%%*t %s ",*(++argv));
					break;
				case 'w':
					argc--;
					sprintf(siz,"%%*w %s ",*(++argv));
					break;
				case 'W':
					argc--;
					sprintf(Siz,"%%*W %s ",*(++argv));
					break;
				case 'v':
					argc--;
					sprintf(Dsiz,"%%*v %s ",*(++argv));
					break;
				case 'g':
					argc--;
					sprintf(grid,"%%*G %s ",*(++argv));
					break;
				case 'c':
					argc--;
					sprintf(points,"%%*c %s ",*(++argv));
					break;
				case 'C':
					argc--;
					sprintf(points,
					"%%*c %s %%*B 0 1 0 1 ",*(++argv));
					break;
				case 'o':
					sprintf(logc,"%%o ");
					break;
				case 'O':
					sprintf(Logc,"%%o ");
					break;
				default:
					break;
			}
		else
			fprintf(stderr,"error\n");
		argv++;
	}
	if(*npoints==0)readin();
	frame("abc");
	np4 = npoints[0]+npoints[1]+npoints[2]+npoints[3];
	np5 = np4+npoints[4];
        np6 = np5+npoints[5];
	i = 0;		/* arguments passed in array instead of on stack */
	argvals[i++]=x;
	argvals[i++]=y;
	argvals[i++]=x+npoints[0];
	argvals[i++]=y+npoints[0];
	argvals[i++]=x+npoints[0]+npoints[1];
	argvals[i++]=y+npoints[0]+npoints[1];
	argvals[i++]=x+npoints[0]+npoints[1]+npoints[2];
	argvals[i++]=y+npoints[0]+npoints[1]+npoints[2];
	argvals[i++]=x+npoints[0]+npoints[1]+npoints[2]+npoints[3];
	argvals[i++]=y+npoints[0]+npoints[1]+npoints[2]+npoints[3];
	argvals[i++]=x+np5;
	argvals[i++]=y+np5;
	argvals[i++]=x+np6;
	argvals[i++]=y+np6;
	argvals[i++]=x+np6+npoints[6];
	argvals[i++]=y+np6+npoints[6];
	argvals[i++]=x+np6+npoints[6]+npoints[7];
	argvals[i++]=y+np6+npoints[6]+npoints[7];
	graph(cmd,argvals);
}
readin()
{
float	*xx,*yy;
int	np,i;
char	*cp;
	xx= x;
	yy= y;
	for(i=0;i<ngraphs;i++){
		xx+=npoints[i];
		yy+=npoints[i];
	}
	np=0;
	if(fgraph==NULL)exit();
	while(fgets(line,120,fgraph)!=NULL){
			if (*line == '#') continue;
			switch(sscanf(line,"%f %f",xx++,yy++)){
				case(0):
					printf("Bad x,y value line %d\n",np+1);
					goto done;
				case(1):
					*(yy-1)= *(xx-1);
					*(xx-1)=np+1;
			}
		np++;
		}
done:
	cp=cmd;
	while(*cp++);
	*(cp-1)=' ';
	npoints[i]=np;
	ngraphs++;
	if (ngraphs >=NGRAPHS) {
		fprintf(stderr,"graph: too many plot files\n");
		ngraphs--;
		return;
	}		
	sprintf(cp,"%%X %s %s %s %s %%Y %s %s %s %s %s %s %s %s %s %s %%*n %d ",
	 label,tic,siz,logc,Label,Logc,min,max,pen,Tic,Siz,grid,Dsiz,points,np);
	*label= *Label= *min= *max= *pen= *tic= *Tic= *grid= *points=0;
	*logc= *Logc= *siz= *Siz= *Dsiz= 0;
}
