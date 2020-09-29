CFLAGS = -O3
# CFLAGS = -O3 -fPIC
# CFLAGS = -g -fPIC
# CFLAGS = -g 
# FCFLAGS = -pg 
# FCFLAGS = -O3 
# CFLAGS = -pg 
# CFLAGS = -g -Wall -Wmissing-prototypes
# CFLAGS = -O5 -DCPML -ffast-math -fomit-frame-pointer -finline-functions
# CFLAGS = -O5 -ffast-math -fomit-frame-pointer -finline-functions
# CFLAGS = -O3 -non_shared  # for DEC UNIX
LDFLAGS = -s
# LDFLAGS = -g
# LDFLAGS = -g
# LDFLAGS = -pg
# LDFLAGS = -s -non_shared  # for DEC UNIX
YFLAGS = -d

NCOBJS = main.o 
NOBJS = nc.o ncmain.o ncfuncs.o ncstimfuncs.o ncsetvar.o ncm.o ncv.o ncomp.o ncsub.o \
	ncconvert.o ncmak.o ncnode.o code.o modcode.o vecsub.o init.o math.o symbol.o \
	initchan.o chanfunc.o channa1.o channa2.o channa3.o channa5.o chanif.o channa6.o \
	chank1.o chankca.o chank3.o chank4.o chanca1.o \
	chancgmp.o chansyn.o chanampa.o channmda.o changaba.o \
	ncstim.o ncplot.o emalloc.o drand.o rndev.o rnd_taus.o rnd_mt.o strtod.o \
	ncdisp.o ncsymb.o ncrot.o prcomp.o gprim.o defaultdraw.o gausnn.o ncray.o \
	ncelem.o makcomp.o ncstimin.o synfilt.o ncio.o ncsave.o scheduler.o

SOBJS = nc.o ncmain.o main.o ncfuncss.o ncsetvar.o ncm.o ncv.o codes.o init.o math.o symbol.o \
	gr.o drand.o rnd_taus.o rnd_mt.o stimsub.o modcodes.o ncmaks.o convarr.o emalloc.o \
	strtod.o gprim.o defaultdraw.o gausnn.o ncelem.o makcomp.o ncstimin.o rndev.o \
	synfilt.o ncnode.o vecsub.o ncio.o ncsave.o

GOBJS = gmain.o gausnnx.o ncio.o

# STIM = XSTIMM		# to make stim ignore neural elems for speed
STIM = XSTIM		# to make stim that knows about with neural elements 

#MCC = mcc -n4
#CC = mcc -n4 -n7 -n8
#CC = mcc -n2 -n3		# for 80387 co-processor
SHELL=/bin/csh
# CXX = g++
# CC = gcc
YACC=../yaccsrc/yacc

AR=ar
ARFLAGS=crv
LIBDIR=.
LIBNC=libnc.a

.SUFFIXES: .cc

#%.o:	%.cc
#	$(CXX) -c $(CFLAGS) $<

.cc.o:	
	$(CXX) -c $(CFLAGS) $*.cc

all:	stim nc ncv nd plotmod gausnn


libnc.a: $(NOBJS)
	$(AR) $(ARFLAGS) $(LIBDIR)/$(LIBNC) $(NOBJS)

libnc.so.1.0.1: $(NOBJS)
	gcc -shared -Wl,-soname,libnc.so.1 -o libnc.so.1.0.1 $(NOBJS)

nc:	$(NCOBJS) $(NOBJS) libnc.a
	$(CXX) $(LDFLAGS) -o nc $(NCOBJS) $(NOBJS) -L. -lnc -L../libP -lP -lm

ncdynamic:$(NCOBJS) libnc.so.1.0.1
	$(CXX) $(LDFLAGS) -o nc $(NCOBJS) -L. -lnc -L../libP -lP -lm

#	$(CXX) -s -i -z -o nc $(NOBJS) -lP -lm
#	$(CXX) ../libnc.a $(LDFLAGS) -o nc -L../libP -lP -lm
#	$(CXX) $(LDFLAGS) -o nc $(NOBJS) -L../libP -lP -lffm -lm
#	# -lffm = use libffm.a fast math library
#	$(CXX) $(LDFLAGS) -o nc $(NOBJS) -L../libP -lP -lcpml
#	# -lcpml = use Compac Portable Math library instead of -lm
#                   50% faster than -lffm
#	$(CXX) -pg       -o nc $(NOBJS) -lbmon -L../libP -lP -lm
#	$(CXX) $(LDFLAGS) -o nc $(NOBJS) -L../libP -lP -lmalloc -lm
#	$(CXX) $(LDFLAGS) -o nc $(NOBJS) -lP 
	@size nc; echo 

ncv:	nc
	ln -s nc ncv

nd:	nc
	ln -s nc nd

ncbnl:	ncbnl.o rndev.o drand.o rnd_taus.o rnd_mt.o emalloc.o 
	$(CXX) -o ncbnl ncbnl.o rndev.o drand.o rnd_taus.o rnd_mt.o emalloc.o -lm

ncbmon:	$(NOBJS)
	$(CXX) -pg -o nc $(NOBJS) -lbmon -L../libP -lP -lm
	@size nc; echo 

nctest:	$(NOBJS) nctest.o
	$(CXX) -o nctest nctest.o $(NOBJS) -L../libP -lP -lm
	@size nc; echo 

stim:	$(SOBJS)
#	$(CXX) $(LDFLAGS) -o stim $(SOBJS) -lmalloc
	$(CXX) $(LDFLAGS) -o stim $(SOBJS) -lm
#	$(CXX) $(LDFLAGS) -o stim $(SOBJS) -lffm -lm
#	# -lffm = use libffm.a fast math library
	@size stim; echo 

modcodes.o: modcode.cc
	$(CXX) -c $(CFLAGS) -D$(STIM) modcode.cc -o modcodes.o

ncmaks.o: ncmak.cc
	$(CXX) -c $(CFLAGS) -D$(STIM) ncmak.cc -o ncmaks.o

ncfuncss.o: ncfuncs.cc
	$(CXX) -c $(CFLAGS) -D$(STIM) ncfuncs.cc -o ncfuncss.o

code.o:	code.cc
	$(CXX) -c $(CFLAGS) code.cc

codes.o: code.cc
	$(CXX) -c $(CFLAGS) -D$(STIM) code.cc -o codes.o

rnd_taus.o: rnd_taus.cc
	$(CXX) -c $(CFLAGS) rnd_taus.cc 

rndev.o: rndev.cc
	$(CXX) -c $(CFLAGS) rndev.cc 

ncomp.o: ncomp.cc
	$(CXX) -c $(CFLAGS) ncomp.cc

ncmain.o: ncmain.cc
#	$(CXX) $(CFLAGS) -DDEC_ALPHA -c ncmain.cc
	$(CXX) -c $(CFLAGS) ncmain.cc

plotmod: plotmod.o ncv.o ncplot.o ncsymb.o gprim.o defaultdraw.o strtod.o \
	synfilt.o emalloc.o ncio.o
	$(CXX) $(LDFLAGS) -o plotmod plotmod.o ncv.o ncplot.o ncsymb.o \
	gprim.o defaultdraw.o strtod.o emalloc.o synfilt.o ncio.o -L../libP -lP -lm
#	gprim.o strtod.o emalloc.o synfilt.o -L../libP -lP -lcpml
	@size plotmod; echo 

gausnnx.o: gausnn.cc
	$(CXX) -c $(CFLAGS) -DNONC gausnn.cc -o gausnnx.o

gausnn: $(GOBJS) drand.o rnd_taus.o rnd_mt.o emalloc.o
	$(CXX) $(LDFLAGS) -o gausnn $(GOBJS) drand.o rnd_taus.o rnd_mt.o \
	emalloc.o -lm
#	emalloc.o -lcpml

rndpnt: rndpnt.o drand.o rnd_taus.o
	$(CXX) $(LDFLAGS) -o rndpnt rndpnt.o drand.o rnd_taus.o -L../libP -lP -lm

binomdev: binomdev.o drand.o rnd_taus.o emalloc.o
	$(CXX) $(LDFLAGS) -o binomdev binomdev.o drand.o rnd_taus.o emalloc.o -lm

gamdev: gamdev.o drand.o rnd_taus.o emalloc.o
	$(CXX) -o gamdev gamdev.o drand.o rnd_taus.o emalloc.o -lm

interp:	nc.o ncmain.o ncm.o code.o init.o math.o symbol.o ncplot.o \
	emalloc.o drand.o rnd_taus.o strtod.o ncsymb.o gprim.o defaultdraw.o interpdum.o
	$(CXX) $(CFLAGS) -o interp nc.o ncmain.o ncm.o code.o init.o \
	math.o symbol.o ncplot.o emalloc.o drand.o rnd_taus.o strtod.o ncsymb.o \
	gprim.o interpdum.o -L../libP -lP -lm

randtest: randtest.o rnd_taus.o emalloc.o
	$(CXX) -o randtest randtest.o rnd_taus.o emalloc.o -lm

linesp:	linesp.c
	$(CXX) -o linesp linesp.c -lm

linespr:linespr.c
	$(CXX) -o linespr linespr.c -lm

linespr2:linespr2.c
	$(CXX) -o linespr2 linespr2.c -lm

test:	test.o 
	$(CXX) $(CFLAGS) -o test test.o -L/usr/local/lib -lgcc -lm

xxx:	xxx.o rndev.o drand.o rnd_taus.o emalloc.o
	$(CXX) $(CFLSGS) -o xxx xxx.o rndev.o drand.o rnd_taus.o emalloc.o -lm

rndtest:rndtest.o rnd_taus.o rnd_glibc.o  rnd_mt.o emalloc.o
	$(CXX) -O3 -o rndtest rndtest.o rnd_taus.o rnd_glibc.o rnd_mt.o emalloc.o -lm

wave:	wave.o 
	$(CXX) $(LDFLAGS) -o wave wave.o -lm
#	$(CXX) $(LDFLAGS) -o wave wave.o -lcpml

tcomp2: tcomp2.o $(NOBJS)
	$(CXX) $(LDFLAGS) -o tcomp2 tcomp2.o -L. -lnc -L../libP -lP -lm
	
tcomp2f: tcomp2f.o $(NOBJS)
	$(CXX) $(LDFLAGS) -o tcomp2f tcomp2f.o -L. -lnc -L../libP -lP -lm
	
tcomp20: tcomp20.o $(NOBJS)
	$(CXX) $(LDFLAGS) -o tcomp20 tcomp20.o -L. -lnc -L../libP -lP -lm
	
tcomp29: tcomp29.o $(NOBJS)
	$(CXX) $(LDFLAGS) -o tcomp29 tcomp29.o -L. -lnc -L../libP -lP -lm
	
tcomp31a: tcomp31a.o $(NOBJS)
	$(CXX) $(LDFLAGS) -o tcomp31a tcomp31a.o -L. -lnc -L../libP -lP -lm
	
tcomp66: tcomp66.o $(NOBJS)
	$(CXX) $(LDFLAGS) -o tcomp66 tcomp66.o -L. -lnc -L../libP -lP -lm

turtlesens: turtlesens.o 
	$(CXX) $(LDFLAGS) -o turtlesens turtlesens.o -lm
	
wave.h:	wave
	./wave > wave.h

ncstim.o stimsub.o: wave.h 

y.tab.h: nc.o

nc.o:	nc.y 
	$(YACC) -d nc.y && mv y.tab.c y.tab.cc	
#	bison -y nc.y
	$(CXX) $(CFLAGS) -c y.tab.cc 
	mv y.tab.o nc.o || rm y.tab.cc

initchan.o: channa1.o channa2.o channa3.o channa5.o chanif.o channa6.o chank1.o chank3.o \
chank4.o chankca.o chanca1.o channmda.o chanampa.o chancgmp.o changaba.o \
chansyn.o 

code.o modcode.o ncdisp.o ncplot.o ncray.o ncstim.o ncsub.o stimsub.o: ndef.h

ncdisp.o ncplot.o ncray.o: colors.h

chanampa.o chanca1.o chancgmp.o chanfunc.o changaba.o chank1.o chank3.o \
chank4.o chankca.o channa1.o channa2.o channa3.o channa5.o chanif.o channa6.o channmda.o \
chansyn.o initchan.o makcomp.o modcode.o ncdisp.o ncconvert.o ncmain.o \
ncmak.o ncomp.o ncplot.o ncstim.o ncsub.o prcomp.o stimcode.o stimsub.o \
synfilt.o modcode.o: control.h

init.o plotmod.o: controlx.h

code.o ncdisp.o ncplot.o: gprim.h

code.o gprim.o init.o modcode.o ncplot.o ncsub.o ncsymb.o stimcode.o \
y.tab.o: gr.h

chanampa.o chanca1.o chancgmp.o chanfunc.o changaba.o chank1.o chank3.o \
chank4.o chankca.o channa1.o channa2.o channa3.o channa5.o chanif.o channa6.o channmda.o \
chansyn.o code.o init.o initchan.o interpdum.o makcomp.o math.o modcode.o \
ncdisp.o ncelem.o ncconvert.o ncm.o ncmain.o ncmak.o ncnode.o ncomp.o ncpacks.o \
ncplot.o ncray.o ncrot.o ncstim.o ncsub.o nctest.o ncupacks.o prcomp.o \
simann.o stimcode.o stimsub.o symbol.o synfilt.o test.o y.tab.o: nc.h

chanfunc.o chankca.o init.o makcomp.o modcode.o ncdisp.o ncelem.o ncconvert.o \
ncmak.o ncnode.o ncomp.o ncpacks.o ncplot.o ncray.o ncstim.o ncsub.o \
ncupacks.o plotmod.o prcomp.o stimcode.o stimsub.o: ncelem.h

ncmak.o ncpacks.o ncsub.o ncupacks.o plotmod.o: nclist.h

chanampa.o chanca1.o chancgmp.o chanfunc.o changaba.o chank1.o chank3.o \
chank4.o chankca.o channa1.o channa2.o channa3.o channa5.o chanif.o channa6.o channmda.o \
chansyn.o drand.o initchan.o makcomp.o modcode.o ncdisp.o ncconvert.o \
ncmak.o ncomp.o ncpacks.o ncplot.o ncrot.o ncslv.o ncstim.o ncsub.o \
ncupacks.o ncvirt.o plotmod.o prcomp.o randtest.o simann.o stimcode.o \
stimsub.o synfilt.o wave.o ncsave.o: ncomp.h

chanampa.o chanca1.o chancgmp.o chanfunc.o changaba.o chank1.o chank3.o \
chank4.o chankca.o channa1.o channa2.o channa3.o channa5.o chanif.o channa6.o channmda.o \
chansyn.o init.o initchan.o math.o modcode.o ncconvert.o ncmak.o ncplot.o \
ncstim.o ncstimin.o ncsub.o prcomp.o rndev.o: nconst.h

code.o interpdum.o modcode.o ncplot.o ncsub.o plotmod.o: ncplot.h

binom.o chanampa.o chanca1.o chancgmp.o chanfunc.o changaba.o chank1.o \
chank3.o chank4.o chankca.o channa1.o channa2.o channa3.o channa5.o chanif.o channa6.o \
channmda.o chansyn.o convarr.o drand.o init.o initchan.o makcomp.o math.o \
modcode.o ncbnl.o ncdisp.o ncelem.o ncconvert.o ncmak.o ncnode.o ncomp.o \
ncplot.o ncray.o ncrot.o ncstim.o ncstimin.o ncsub.o plotmod.o prcomp.o \
rnd_glibc.o randtest.o rndev.o simann.o stimsub.o synfilt.o wave.o: ncsub.h

binom.o chankca.o code.o convarr.o drand.o defaultdraw.o emalloc.o emalloc2.o gausnn.o \
gmain.o gprim.o init.o initchan.o makcomp.o math.o modcode.o ncdisp.o ncelem.o \
ncm.o ncmain.o ncmak.o ncnode.o ncomp.o ncplot.o ncray.o ncrot.o ncstim.o \
ncsub.o ncsymb.cc ncupacks.o ncvirt.o prcomp.o rndev.o stimsub.o symbol.o \
synfilt.o vecsub.o: ncio.h

init.o ncmain.o nctest.o plotmod.o test.o: ncval.h

ncpacks.o ncslv.o ncupacks.o ncvirt.o: ncvirt.h

code.o frame.o gprim.o gr.o init.o modcode.o ncdisp.o ncplot.o ncsymb.o \
plotmod.o rndpnt.o: stdplt.h

convarr.o modcode.o stimcode.o stimsub.o: stim.h

gprim.o: extdraw.h

wave.o: turtle.h

modcode.o vecsub.o: vec.h

ncstim.o stimsub.o: wave.h

chanampa.o chanca1.o chancgmp.o chanfunc.o changaba.o chank1.o chank3.o \
chank4.o chankca.o channa1.o channa2.o channa3.o channa5.o chanif.o channa6.o channmda.o \
chansyn.o code.o init.o initchan.o makcomp.o math.o modcode.o ncdisp.o \
ncelem.o ncconvert.o ncm.o ncmain.o ncmak.o ncnode.o ncomp.o ncpacks.o ncplot.o \
ncray.o ncstim.o ncsub.o nctest.o ncupacks.o prcomp.o stimcode.o symbol.o \
synfilt.o stimsub.o test.o: y.tab.h

x.tab.h:	y.tab.h
	-cmp -s x.tab.h y.tab.h || cp y.tab.h x.tab.h

clean:
	rm -f $(NOBJS) $(SOBJS) wave wave.o y.tab.cc y.tab.h nc libnc.a stim \
	 plotmod plotmod.o gausnn $(GOBJS) core

backup:	
	@echo "Please insert the 'nc' backup disk; type <RET> when ready..."
#	@line
	cd ..;tar cvf /dev/fd0 nc/{*.cc,*.c,*.y,*.h,*.m,fib*,acker*,pfilt,neurc*,nc.pov,makefile}
	@echo nc is backed up...

backup1:	
	@echo "Please insert the 'nc' backup disk drive f0"
	@echo "Type <RET> when ready..."
	@line
	cd ..; tar cv nc/{*.cc,*.c,*.y,*.h,*.m,fib*,acker*,pfilt,nc.pov,neurc*}
	@echo g++ nc is backed up...

databak:
	@echo "Please insert the 'nc data' backup disk,"
	@echo " and type <RET> when ready..."
	@line
	cd ..; tar cv nc/{acker*,tcomp*,*.n,*.m,*.doc,makefile}
	@echo nc is backed up...

backhz:
	@echo "Please insert the 'hz data' backup disk,"
	@echo " and type <RET> when ready..."
	@line
	tar cv makefile hz*
	@echo hz files are backed up...

copy:
	cd ..; tar cvf imodc imod/{*.cc,*.c,*.y,*.h,fib*,acker*,tcomp*,nc.pov,makefile}
	@echo nc is backed up...

