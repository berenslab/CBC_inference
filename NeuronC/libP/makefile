CFLAGS= -I. -O
CC=gcc

LIBOBJS= _cat.o _cfname.o _clip.o _domove.o _draw.o _err.o _find.o _fmatch.o \
_inv2.o _move.o _origin.o _path.o _plsub.o _rotate.o _scale.o _window.o \
_xhair.o _drcirc.o _drrect.o _drtri.o axis.o base.o cfont.o crotate.o \
cwidth.o dash.o data.o display.o dispts.o dodisp.o domin.o donice.o erase.o \
erasea.o eval.o fat.o fgrid.o fmod.o frame.o fwhere.o grid.o minmax.o misc.o \
pair.o prframe.o purge.o rmframe.o saxis.o textf.o where.o xdist.o \
xhair.o ydist.o zinit.o putws.o btest.o

all:	libP.a y.tab.o xgraph.o

libP.a:	$(LIBOBJS)
	ar vr libP.a $(LIBOBJS)
	ranlib libP.a

$(LIBOBJS): stdplt.h

clean:
	rm -f $(LIBOBJS) libP.a

axis.o base.o display.o dispts.o dodisp.o domin.o donice.o fgrid.o \
grid.o minmax.o misc.o pair.o saxis.o xgraph.o y.tab.o zinit.o: graph.h

y.tab.o: graph.y
	../yaccsrc/yacc -d graph.y
	cc $(CFLAGS) -c y.tab.c

backup:
	@echo "Please insert the 'plot' backup disk; type <RET> when ready..."
	@line
	cd ..; tar cvf /dev/fd0 libP/{*.c,*.h,*.y,makefile,N*,R*}
	@echo plot is backed up.

addbackup:
	@echo "Please insert the 'plot' backup disk; type <RET> when ready..."
	@line
	cd ..; tar rvf /dev/fd0 libP/{*.c,*.h,*.y,makefile,N*,R*}
	@echo plot is backed up.
 


