CFLAGS = -O -I/usr/lib -I. -I../libP -I../pl
#
# for Mac, use this line and uncomment 
# CFLAGS = -O -I/opt/X11/include -I/usr/lib -I.
#
HD = .
PL = ../pl
B = ../bin
CC=gcc

PDOBS = pdigit.o plotdr.o plotsub.o


#for Venix/386:
#
#all:    vid textmod hpflt ibmflt hplot iplot labels posit Plotdb
#
#  for all drivers:
#
#VOBJS = vid.o mprinta.o mprintc.o mprinte.o mprintv.o mprinth.o mprintl.o \
#	mprintx.o hfile.o hvirt.o mvid.o iopriv.o putsym.o getws.o btest.o
#
#  for just VGA and PS drivers:
#
#VOBJS = vid.o mprinta.o mprintc.o mprintv.o mvid.o iopriv.o \
#	putsym.o getws.o btest.o
#
#vid:	$(VOBJS)
#	$(CC) -s -o vid $(VOBJS) -lrtx -lX11_s -lpt -lsocket -lnet -lnsl_s -lm
#	@size vid
#	$B/access vid
#
#       On Venix, must "su root" before running
#       the following commands (equivalent to running
#	the "access vid" command).  This allows VGA and EGA
#       drivers to access the hardware:
#
#	chown root vid
#	chmod ug+s vid
#
#
# for most workstations:
# (remember to comment out unwanted driver defs in mvid.c, e.g. THERC)
#
VOBJS = vid.o mprinta.o mprintb.o mprintc.o mprintxf.o mprintx.o \
	mvid.o putsym.o getws.o btest.o

all:    vid labels posit Plotdb
	ln -f vid ../bin

vid:	$(VOBJS)
	$(CC) -o vid $(VOBJS) -L/usr/X11R6/lib -lX11 -lm
	@size vid

chars.h: asciichar greekchar 
	make symcomp
	./symcomp ascii <asciichar >chars.h
	echo "" >>chars.h
	./symcomp greek <greekchar >>chars.h

putsym.h: chars.h
	touch putsym.h

hpflt:	hpflt.o plothp.o getws.o btest.o
	$(CC) -s -o hpflt hpflt.o plothp.o getws.o btest.o -lm
	@size hpflt

ibmflt:	ibmflt.o plotibm.o getws.o btest.o
	$(CC) -s -o ibmflt ibmflt.o plotibm.o getws.o btest.o -lm
	@size ibmflt

labels: labels.o label.o mdraw.o 
	$(CC) -s -o labels labels.o label.o mdraw.o -L../libP -lP -lm
	@size labels

pltest: pltest.o
	$(CC) -s -o pltest pltest.o cwidth.o -lP -lm
	@size pltest

Tek:	Tek.o getws.o btest.o
	$(CC) -s -o Tek Tek.o getws.o btest.o -lm

Plotdb:	Plotdb.o getws.o btest.o
	$(CC) -s -o Plotdb Plotdb.o getws.o btest.o

posit:	posit.o getws.o putws.o btest.o
	$(CC) -s -o posit getws.o posit.o putws.o btest.o
	@size posit

hvirt.o: hvirt.c
	$(CC) -c -DHIDE12 $(CFLAGS) hvirt.c

textmod:textmod.o mprinte.o mprinth.o mprintv.o iopriv.o
	$(CC) -s -o textmod textmod.o mprinte.o mprinth.o mprintv.o \
	iopriv.o -lrtx -lm 
	@size textmod

wscreen:wscreen.o memio.o
	$(CC) -s -o wscreen wscreen.o memio.o -lm 
	@size wscreen

scrtest: scrtest.o vidio.o 
	$(CC) -s -o scrtest scrtest.o vidio.o -lm
	@size scrtest

hplot:	hplot.o plotdr.o
	$(CC) -s -o hplot hplot.o plotdr.o
	@size hplot 

iplot:	iplot.o plotidr.o
	$(CC) -s -o iplot iplot.o plotidr.o
	@size iplot 

pdigit:	$(PDOBS)
	$(CC) -s -o  pdigit $(PDOBS)
	@size pdigit

plabel:	plabel.o plotsub.o
	$(CC) -s -o plabel plabel.o plotsub.o
	@size plabel

graph:	graph.o ../libP/xgraph.o ../libP/y.tab.o
	$(CC) -s -o graph graph.o ../libP/xgraph.o ../libP/y.tab.o -L../libP -lP -lm

# label.o mdraw.o graph.o plsub.o: stdplt.h

graphx:	graph.o misc.o donice.o
	$(CC) -s -o graph graph.o misc.o donice.o -lP -lm

symcomp: symcomp.c
	$(CC) -s -o symcomp symcomp.c

vtest:	vtest.o mprintv.o iopriv.o 
#	$(CC) -s -o vtest vtest.o mprintm.o -ltxpc
#	$(CC) -s -o vtest vtest.o mprinte.o iopriv.o -lrtx
	$(CC) -s -o vtest vtest.o mprintv.o iopriv.o -lrtx
#	$(CC) -s -o vtest vtest.o egadriv.o iopriv.o inout.o -ltxpc
	$B/access vtest

erasee: erasee.o mprinte.o
	$(CC) -s -o erasee erasee.o mprinte.o

clean:
	rm -f vid graph labels Plotdb posit chars.h *.o

backup:
	@echo "Please insert the 'pl' backup disk; type <RET> when ready..."
	@line
	cd ..;tar cvf /dev/fd0  pl/{*.c,*.s,*.h,*.hlp,asciichar,greekchar,icon*,makefile*}
	@echo pl is backed up...

addbackup:
	@echo "Please insert the 'pl' backup disk; type <RET> when ready..."
	@line
	cd ..;tar rvf /dev/fd0  pl/{*.c,*.s,*.h,*.hlp,asciichar,greekchar,icon*,makefile*}
	@echo pl is backed up...

# dependencies for plotter routines

Tek.o vid.o hpflt.o ibmflt.o: putsym.h

vid.o textmod.o mvid.o: ./vdef.h

hplot.o pdigit.o plabel.o plotdr.o plotsub.o label.o \
labels.o: mdef.h

