/*
 * prog to compile the symbol descriptions into the form
 * used by subroutine putsym
 *
 * The one argument to symcomp is the generic name to be
 * used in the output structures (i.e. ascii). If no arg
 * is given the 'XXX' is used.
 *
 * The input characters are designed on a (0-4) x (0-7) grid.
 * Each line of the input file defines one character.
 * Each line starts with the character being defined followed
 * by a ':'. If the next character is a 'g' then the character
 * is taken to be a descender (like g) and the special drop bit
 * is incorporated into the compiled character code. The
 * rest of the line consists of moves ('m') or draws ('d'), each
 * followed by two one digit numbers specifing the x-y co-ords
 * of the operation. No embedded blanks are allowed anywhere
 * in the line. Descenders are designed in a raised position.
 * The drop bit will cause the character to be lowered two units
 * when it is plotted. The first character defined is will be
 * the character for all undefined letters in the 96 letter
 * alphabet. Typically the blank is defined first.
 */
#include <stdio.h>
#include <stdlib.h>

#define NSYMBOL 96
#define DROPBIT 02000

int address[NSYMBOL];
char vectors[2000];
char dname[] = "XXX";

int main(int argc, char **argv)
{
	register int *symadd;
	register char *symvec;
	int i,add,count,x,y,num,sym,c;
	char vecent;
	char *name;

	argc--; 
	argv++;
	if(argc)	name= *argv;
	else		name= dname;
	for(i=0;i<NSYMBOL;i++) address[i]= 0;
	symadd= address;
	symvec= vectors;
	count=0;
	while((sym=getchar()) != EOF)
	{
		if((c=getchar()) != ':')
		{
			fprintf(stderr, "error colon trap %c %c\n",sym,c);
			exit(1);
		}
		add=count;
		while((c=getchar()) != '\n')
		{
			if(c=='g')
			{
				add |= DROPBIT;
				continue;
			}
			x=getchar() -'0';
			y=getchar() -'0';
			vecent=0;
			switch(c)
			{
			case 'm': 
				vecent= 0200;
			case 'd': 
				vecent |= (x<<4)+y;
				break;
			default: 
				fprintf(stderr, "unknown plot command %c %c \n",sym,c);
				exit(1);
			}
			symvec[count++]= vecent;
		}
		symvec[count-1] |= 010;
		symadd[sym-040]= add;
	}
	/*
	 * print the address table
	 */
	printf("int %sadd[] = {",name);
	for(i=0;i<NSYMBOL;i++)
	{
		if(i%10) printf("0%o",address[i]);
		else  printf("\n\t0%o",address[i]);
		if( i != (NSYMBOL-1) ) printf(",");
	}
	printf("\n};\n\n");
	/*
	 * print the vector table
	 */
	printf("char %svec[] = {",name);
	for(i=0;i<count;i++)
	{
		num= vectors[i];
		num &= 0377;
		if(i%16) printf("0%o",num);
		else   printf("\n\t0%o",num);
		if( i != (count-1) ) printf(",");
	}
	printf("\n};\n\n");
	/*
	 * print the structure
	 */
	printf("struct {\n\tint *saddr;\n\tchar *svec;\n");
	printf("} %s = {\n\t%sadd,\n\t%svec\n};\n",name,name,name);
	/*
	 * list the symbols compiled
	 */
	fprintf(stderr, "\nfollowing symbols compiled\n\n");
	num=0;
	for(i=0;i<NSYMBOL;i++)
	{
		c= i+040;
		if(vectors[i]== 0) continue;
		else   fprintf(stderr, "%c",c);
		num++;
		if((num%20)==0) fprintf(stderr, "\n");
	}
	fprintf(stderr, "\n\nhigh water mark %d bytes\n",count);
	exit(0);
}
