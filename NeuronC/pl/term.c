#include <stdio.h>
#include <curses.h>

#define	STDCW	192		/* world coordinates */
#define	XMAX	16384
#define	YMAX	16384
#define PLTCMD	0200		/* flag bit for plot commands */
#define CMDMSK	0177		/* mask for value of plot commands */
#define	NDIV	5

int	tsize, tangle, tfont, fat;
int	xold, yold, xorigin, yorigin;

int	xmin, xmax, ymin, ymax;	/* device coordinates */
int	xbot, ybot;
double	abot, bbot, cbot;
double	sclx, scly;
int	tty;
int	rflag;			/* print reference scale */

extern int getwsx();
double	atof();

main(argc, argv)
char **argv;
{
	register char *arg;
	char ch;

/*	tty = open("/dev/tty", 0); */
	initscr();
	cbot = 1;
	while (--argc > 0) {
		arg = *++argv;
		if (*arg++ != '-')
			break;
		switch (*arg++) {
		case 'x':
			abot = atof(arg);
			xbot = abot * 16384.;
			rflag++;
			break;
		case 'y':
			bbot = atof(arg);
			ybot = bbot * 16384.;
			rflag++;
			break;
		case 's':
			cbot = atof(arg);
			rflag++;
			break;
		case 'r':
			rflag++;
			break;
		}
	}
	xmax = COLS;
	ymax = LINES;
	sclx = cbot * 5.5 / 7.5 * xmax / 16384.;
	scly = cbot * ymax / 16384.;
	if (argc > 0) for (; ;) {
		if (freopen(*argv++, "r", stdin))
			doplot();
		else
			perror(*--argv);
		if (--argc <= 0)
			break;
/*		read(tty, &ch, 1); */
	} else
		doplot();
	exit(0);
}

line(x1, y1, x2, y2)
{
	register x, y;
	int dx, dy;
	int xinc, yinc;
	register res1;
	int res2;

	x1 = (x1 - xbot) * sclx;
	y1 = LINES - 1 - (y1 - ybot) * scly;
	x2 = (x2 - xbot) * sclx;
	y2 = LINES - 1 - (y2 - ybot) * scly;
	if (clip(&x1, &y1, &x2, &y2))
		return;
	xinc = 1;
	yinc = 1;
	x = x1;
	y = y1;
	if ((dx = x2-x) < 0) {
		xinc = -1;
		dx = -dx;
	}
	if ((dy = y2-y) < 0) {
		yinc = -1;
		dy = -dy;
	}
	res1 = 0;
	res2 = 0;
	if (dx >= dy) while (x != x2) {
		mvaddch(y, x, '+');
		if (res1 > res2) {
			res2 += dx - res1;
			res1 = 0;
			y += yinc;
		}
		res1 += dy;
		x += xinc;
	} else while (y != y2) {
		mvaddch(y, x, '+');
		if (res1 > res2) {
			res2 += dy - res1;
			res1 = 0;
			x += xinc;
		}
		res1 += dx;
		y += yinc;
	}
	mvaddch(y2, x2, '+');
}

doplot()
{
	register int i;
	register int c;
	int xnew, ynew;
	int eras = 0;
	char ch;

	tsize = STDCW;
	while ((c = getc(stdin)) != EOF) {
		if ((c&PLTCMD) == 0 || ((c&CMDMSK) < 'a' || (c&CMDMSK) > 'z')) {
			do putc(c, stderr);
				while (((c = getc(stdin)) & PLTCMD) == 0);
			if (c == EOF) break;
			ungetc(c, stdin);
			continue;
		}
		/*
		 * c is a lower case letter with PLTCMD on
		 */
		switch (c & CMDMSK) {
		case 'm':		/* move */
			xold = xorigin = getwsx(stdin);
			yold = yorigin = getwsx(stdin);
			break;
		case 'd':		/* draw */
			xnew = getwsx(stdin);
			ynew = getwsx(stdin);
			if(fat) {
				if(abs(xnew-xold) >= abs(ynew-yold))
					for(i = -(fat/2); i <= (fat+1)/2; i++)
						line(xold,yold+i,xnew,ynew+i);
				else
					for(i = -(fat/2); i <= (fat+1)/2; i++)
						line(xold+i,yold,xnew+i,ynew);
			}
			else
				line(xold,yold,xnew,ynew);
			xold = xnew;
			yold = ynew;
			break;
		case 's':	/* set angle and size for text */
			tsize = (getc(stdin)&0377)*64;
			tfont = getc(stdin)&0377;
			tangle = getwsx(stdin);
			break;
		case 't':		/* text */
			putsym();
			break;
		case 'p':		/* purge plot buffers */
		        fflush (stdout);
			break;
		case 'b':		/* break */
		case 'e':		/* erase */
			/* ignore first 'e' or 'b' */
			if (eras == 0) {
				eras = 1;
				break;
			}
			purge();
	/*		read(tty, &ch, 1); */
			break;
		case 'f':		/* fat */
			fat = getwsx(stdin);
			break;
		case 'r':		/* raster data */
			fprintf(stderr, "raster data not in yet\n");
			return;
		case 'w':		/* window */
			xmin = (unsigned)getwsx(stdin)*sclx;
			if(xmin<0) xmin = 0;
			xmax = (unsigned)getwsx(stdin)*sclx;
			if(xmax>XMAX) xmax = XMAX;
			ymin = (unsigned)getwsx(stdin)*scly;
			if(ymin<0) ymin = 0;
			ymax = (unsigned)getwsx(stdin)*scly;
			if(ymax>YMAX) ymax = YMAX;
			break;
		case 'c':		/* change pen command */
			getc(stdin);
			break;
		case 'x':		/* display X-hairs */
		case 'n':		/* no op */
		case 'h':		/* hidden line commands R.G.S */
		case 'i':
		case 'j':
			break;
		default: 		/* error */
			fprintf(stderr,"invalid plot command %c\n",c&CMDMSK);
			return;
		}
	}
	purge();
}

/* code for clip */
#define code(x,y) (x<xmin?1:(x>xmax?2:0))|(y<ymin?4:(y>ymax?8:0))

/* window the plot */
clip(x1, y1, x2, y2)
int *x1, *y1, *x2, *y2;
{
	register int c1, c2, temp;
	int swap;

	c1 = code(*x1, *y1);
	c2 = code(*x2, *y2);
	swap = 0;
	if (!(c1||c2))	/* line completely in bounds */
		return(0);
	while (c1|c2) {
		if (c1&c2)	/* line completely out of bounds */
			return(1);
		if (!c1) {
			/* interchange endpoints */
			temp = *x1; *x1 = *x2; *x2 = temp;
			temp = *y1; *y1 = *y2; *y2 = temp;
			temp = c1; c1 = c2; c2 = temp;
			swap = ~swap;
		}
		if (c1<4) {
			/* move endpoint in x */
			temp = (c1&2?xmax:xmin);
			*y1 = solve(temp, *x1, *y1, *x2, *y2);
			*x1 = temp;
		} else {
			/* move endpoint in y */
			temp = (c1&8?ymax:ymin);
			*x1 = solve(temp, *y1, *x1, *y2, *x2);
			*y1 = temp;
		}
		c1 = code(*x1, *y1);
	}
	/* put endpoints in order */
	if (swap) {
		temp = *x1; *x1 = *x2; *x2 = temp;
		temp = *y1; *y1 = *y2; *y2 = temp;
	}
	return(0);
}

/* solve linear eqn in integers */
solve(pnot, p1, q1, p2, q2)
{
	register int pmid, qmid;

	if (p1>p2) {
		pmid = p1; p1 = p2; p2 = pmid;
		qmid = q1; q1 = q2; q2 = qmid;
	}
	if (pnot<=p1)
		return(q1);
	if (pnot>=p2)
		return(q2);
	/* iterate until convergence */
	for(; ;) {
		pmid = (p1+p2) >> 1;
		qmid = (q1+q2) >> 1;
		if (pmid<pnot) {
			p1 = pmid; q1 = qmid;
		} else if (pmid>pnot) {
			p2 = pmid; q2 = qmid;
		} else return(qmid);
	}
}

purge()
{
	register int i;

	if (rflag)
		for (i = 1; i < NDIV; i++) {
			move(LINES-1, COLS*i/NDIV - 2);
			printw("%2.2f", abot + i*7.5/(cbot*NDIV*5.5));
			move(LINES*(NDIV-i)/NDIV - 1, 0);
			printw("%2.2f", bbot + i/(cbot*NDIV));
		}
	refresh();
	mvcur(-1, -1, 0, 0);
	clear();
}
