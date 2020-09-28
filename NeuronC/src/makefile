# CFLAGS = -Ofast
CFLAGS = -O3
# CFLAGS =      # for Mac OSX 
# CFLAGS = -g
# CFLAGS = -g -O3 
# FCFLAGS = -pg 
# FCFLAGS = -O3 
# CFLAGS = -pg 
# CFLAGS = -g -Wall -Wmissing-prototypes
# CFLAGS = -O5 -DCPML -ffast-math -fomit-frame-pointer -finline-functions
# CFLAGS = -O5 -ffast-math -fomit-frame-pointer -finline-functions
# CFLAGS = -O3 -non_shared  # for DEC UNIX
# LDFLAGS = -s
# LDFLAGS = -g
# LDFLAGS = -g
# LDFLAGS = -pg
# LDFLAGS = -s -non_shared  # for DEC UNIX
YFLAGS = -d

NCOBJS = main.o 
NOBJS = nc.o ncmain.o ncfuncs.o ncstimfuncs.o ncsetvar.o ncm.o ncv.o ncomp.o ncsub.o \
	ncconvert.o ncmak.o ncnode.o code.o modcode.o vecsub.o init.o math.o symbol.o \
	initchan.o chanfunc.o channa1.o channa2.o channa3.o channa5.o chanif.o channa6.o \
	chank1.o chankca.o chank3.o chank4.o chank5.o chank10.o chanclca.o chanca1.o \
	chancgmp.o chansyn.o chanampa.o channmda.o changaba.o changly.o chanchr2.o ncstim.o \
	ncplot.o emalloc.o drand.o rndev.o rnd_taus.o rnd_taus113.o rnd_tinymt64.o rnd_mt.o strtod.o \
	ncdisp.o ncsymb.o ncrot.o prcomp.o gprim.o defaultdraw.o gausnn.o ncray.o \
	ncelem.o makcomp.o ncstimin.o synfilt.o digfilt.o makfilter.o complex.o ncio.o \
	ncsave.o scheduler.o lm_funcs.o lmmin.o lm_eval.o lm_evali.o stimsub.o convarr.o \
	rodi.o chr2.o

SOBJS = nc.o ncmains.o main.o ncsetvar.o ncm.o ncv.o codes.o init.o math.o symbol.o \
	gr.o drand.o rnd_taus.o rnd_taus113.o rnd_tinymt64.o rnd_mt.o stimsub.o modcodes.o \
	ncmak.o convarr.o emalloc.o strtod.o gprim.o defaultdraw.o gausnn.o ncelem.o makcomp.o \
	ncstimin.o rndev.o synfilt.o digfilt.o makfilter.o complex.o ncnode.o vecsub.o ncio.o \
	ncsave.o lm_funcs.o lmmin.o lm_eval.o lm_evali.o rodi.o chr2.o stimdum.o

GOBJS = gmain.o gausnnx.o ncio.o

# STIM = XSTIMM		# to make stim ignore neural elems for speed
STIM = XSTIM		# to make stim that knows about neural elements 

#MCC = mcc -n4
#CC = mcc -n4 -n7 -n8
#CC = mcc -n2 -n3		# for 80387 co-processor
# CXX = g++
# CC = gcc
YACC=../yaccsrc/yacc

AR=ar
ARFLAGS=crv
LIBDIR=.
LIBNC=libnc.a
NCFIL = stim nc ncv nd plotmod gausnn

.SUFFIXES: .cc

#%.o:	%.cc
#	$(CXX) -c $(CFLAGS) $<

.cc.o:	
	$(CXX) -c $(CFLAGS) -I../libP -I../pl $*.cc

all:	$(NCFIL)
	ln -f $(NCFIL) ../bin

libnc.a: $(NOBJS)
	$(AR) $(ARFLAGS) $(LIBDIR)/$(LIBNC) $(NOBJS)
#	ranlib libnc.a

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

ncbnl:	ncbnl.o rndev.o drand.o rnd_taus.o rnd_taus113.o rnd_mt.o emalloc.o 
	$(CXX) -o ncbnl ncbnl.o rndev.o drand.o rnd_taus.o rnd_taus113.o rnd_mt.o emalloc.o -lm

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

ncmains.o: ncmain.cc
	$(CXX) -c $(CFLAGS) -D$(STIM) ncmain.cc -o ncmains.o

modcodes.o: modcode.cc
	$(CXX) -c $(CFLAGS) -D$(STIM) -I../libP -I../pl modcode.cc -o modcodes.o

stimdum.o: stimdum.cc
	$(CXX) -c $(CFLAGS) -D$(STIM) stimdum.cc -o stimdum.o

codes.o: code.cc
	$(CXX) -c $(CFLAGS) -D$(STIM) -I../libP -I../pl code.cc -o codes.o

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
	synfilt.o digfilt.o emalloc.o ncio.o notinit.o openzfil.o
	$(CXX) $(LDFLAGS) -o plotmod plotmod.o ncv.o ncplot.o ncsymb.o gprim.o \
	defaultdraw.o strtod.o emalloc.o synfilt.o digfilt.o makfilter.o complex.o \
	ncio.o notinit.o openzfil.o -L../libP -lP -lm
#	gprim.o strtod.o emalloc.o synfilt.o -L../libP -lP -lcpml
	@size plotmod; echo 

gausnnx.o: gausnn.cc
	$(CXX) -c $(CFLAGS) -DNONC gausnn.cc -o gausnnx.o

gausnn: $(GOBJS) drand.o rnd_taus.o rnd_taus113.o rnd_mt.o emalloc.o
	$(CXX) $(LDFLAGS) -o gausnn $(GOBJS) drand.o rnd_taus.o rnd_taus113.o rnd_mt.o \
	emalloc.o -lm
#	emalloc.o -lcpml

gaussfit: libnc.a gaussfit.o 							# link with libnc.a
	$(CXX) $(LDFLAGS) gaussfit.o -o gaussfit -L. -lnc -L../libP -lP -lm  

gaussfitn: gaussfitn.o lm_funcs.o lm_eval.o lmmin.o ncsetvar.o ncio.o	        # link separately with lmfuncs
	$(CXX) $(LDFLAGS) gaussfitn.o lm_funcs.o lm_eval.o lmmin.o ncsetvar.o ncio.o -o gaussfitn -lm   

modelfit: libnc.a modelfit.o 							# link with libnc.a
	$(CXX) $(CFLAGS) $(LDFLAGS) modelfit.o -o modelfit -L. -lnc -L../libP -lP -lm  

# modelfit: modelfit.o lm_funcs.o lm_eval.o lmmin.o ncsetvar.o ncio.o	        # link separately with lmfuncs
#	$(CXX) $(CFLAGS) $(LDFLAGS) modelfit.o lm_funcs.o lm_eval.o lmmin.o ncsetvar.o ncio.o -o modelfit -lm   

makgauss: makgauss.o ncsetvar.o ncio.o	   
	$(CXX) $(LDFLAGS) makgauss.o ncsetvar.o ncio.o -o makgauss -lm   

lm_test: lm_test.o lm_funcs.o lm_eval.o lmmin.o	ncio.o
	$(CXX) $(LDFLAGS) lm_test.o lm_funcs.o lm_eval.o lmmin.o ncio.o -o lm_test -lm   

lm_test_sim: lm_test_sim.o lm_funcs.o lm_eval.o lmmin.o	ncio.o
	$(CXX) $(LDFLAGS) lm_test_sim.o lm_funcs.o lm_eval.o lmmin.o ncio.o -o lm_test_sim -L. -lnc -L../libP -lP -lm   

rndpnt: rndpnt.o drand.o rnd_taus.o rnd_taus113.o rnd_tinymt64.o rnd_mt.o ncio.o emalloc.o
	$(CXX) $(LDFLAGS) -o rndpnt rndpnt.o drand.o rnd_taus.o rnd_taus113.o rnd_tinymt64.o rnd_mt.o \
		ncio.o emalloc.o -L../libP -lP -lm

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

rndtest:rndtest.o rnd_taus.o rnd_taus113.o rnd_glibc.o rnd_mt.o emalloc.o ncio.o
	$(CXX) -O3 -o rndtest rndtest.o rnd_taus.o rnd_taus113.o rnd_glibc.o rnd_mt.o emalloc.o ncio.o -lm

test_digfilt: test_digfilt.o digfilt.o rndev.o drand.o rnd_taus.o emalloc.o ncio.o
	        $(CXX) -o test_digfilt test_digfilt.o digfilt.o rndev.o drand.o rnd_taus.o emalloc.o ncio.o -lm

wave:	wave.o 
	$(CXX) $(LDFLAGS) -o wave wave.o -lm
#	$(CXX) $(LDFLAGS) -o wave wave.o -lcpml

tcomp2: tcomp2.o libnc.a
	$(CXX) $(LDFLAGS) -o tcomp2 tcomp2.o -L. -lnc -L../libP -lP -lm
	
tcomp2c: tcomp2c.o libnc.a
	$(CXX) $(LDFLAGS) -o tcomp2c tcomp2c.o -L. -lnc -L../libP -lP -lm
	
tcomp2f: tcomp2f.o libnc.a
	$(CXX) $(LDFLAGS) -o tcomp2f tcomp2f.o -L. -lnc -L../libP -lP -lm
	
tcomp20: con2.cc tcomp20.cc libnc.a
	$(CXX) $(LDFLAGS) -o tcomp20 tcomp20.cc -L. -lnc -L../libP -lP -lm
	
tcomp20h: con2.cc tcomp20h.cc libnc.a
	$(CXX) $(LDFLAGS) -o tcomp20h tcomp20h.cc -L. -lnc -L../libP -lP -lm
	
tcomp29: tcomp29.o libnc.a
	$(CXX) $(LDFLAGS) -o tcomp29 tcomp29.cc -L. -lnc -L../libP -lP -lm
	
tcomp31a: tcomp31a.o libnc.a
	$(CXX) $(LDFLAGS) -o tcomp31a tcomp31a.o -L. -lnc -L../libP -lP -lm
	
tcomp66: tcomp66.o libnc.a
	$(CXX) $(LDFLAGS) -o tcomp66 tcomp66.o -L. -lnc -L../libP -lP -lm

tcomp71: tcomp71.o libnc.a
	$(CXX) $(LDFLAGS) -o tcomp71 tcomp71.o -L. -lnc -L../libP -lP -lm

poisdist: poisdist.o libnc.a
	$(CXX) $(LDFLAGS) -o poisdist poisdist.o -L. -lnc -L../libP -lP -lm

gausinteg: gausinteg.o libnc.a
	$(CXX) $(LDFLAGS) -o gausinteg gausinteg.o -L. -lnc -L../libP -lP -lm

gausdist: gausdist.o libnc.a
	$(CXX) $(LDFLAGS) -o gausdist gausdist.o -L. -lnc -L../libP -lP -lm

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
chank4.o chank5.o chank10.o chankca.o chanclca.o chanca1.o channmda.o chanampa.o chancgmp.o changaba.o \
changly.o chansyn.o 

code.o modcode.o ncdisp.o ncomp.o ncplot.o ncray.o ncstim.o ncstimfuncs.o ncsub.o stimsub.o: ndef.h

ncdisp.o ncplot.o ncray.o: colors.h

chanampa.o chanca1.o chancgmp.o chanfunc.o changaba.o changly.o chanif.o chank1.o chank3.o \
chank4.o chank5.o chank10.o chankca.o chanclca.o channa1.o channa2.o channa3.o channa5.o \
channa6.o channmda.o \
chansyn.o code.o initchan.o lm_evali.o makcomp.o modcode.o ncconvert.o \
ncdisp.o ncelem.o ncfuncs.o ncmain.o ncmak.o ncomp.o ncplot.o \
ncsave.o ncstim.o ncstimfuncs.o ncsub.o prcomp.o stimsub.o \
synfilt.o: control.h

init.o plotmod.o: controlx.h

code.o modcode.o ncdisp.o ncfuncs.o ncplot.o: gprim.h

code.o gprim.o init.o modcode.o ncplot.o ncsub.o ncsymb.o rndpnt.o y.tab.o: ../pl/gr.h

chanampa.o chanca1.o chancgmp.o chanfunc.o changaba.o changly.o chanif.o chank1.o chank3.o \
chank4.o chank5.o chank10.o chankca.o chanclca.o channa1.o channa2.o channa3.o channa5.o \
channa6.o channmda.o \
chansyn.o code.o init.o initchan.o interpdum.o lm_eval.o lm_evali.o makcomp.o \
math.o modcode.o ncconvert.o ncdisp.o ncelem.o ncfuncs.o ncm.o ncmain.o ncmak.o \
ncnode.o ncomp.o ncpacks.o ncplot.o ncray.o ncrot.o ncsave.o ncsetvar.o ncstim.o \
ncstimfuncs.o ncsub.o nctest.o ncupacks.o ncv.o prcomp.o simann.o stimdum.o  \
stimsub.o symbol.o synfilt.o y.tab.o: nc.h

chanfunc.o chankca.o chanclca.o init.o makcomp.o modcode.o ncconvert.o ncdisp.o ncelem.o \
ncfuncs.o ncmak.o ncnode.o ncomp.o ncpacks.o ncplot.o ncray.o ncstim.o \
ncstimfuncs.o ncsub.o ncsub.newsyn.o ncupacks.o plotmod.o prcomp.o stimdum.o \
stimsub.o: ncelem.h

ncmak.o ncpacks.o ncsub.o ncupacks.o plotmod.o: nclist.h

chanampa.o chanca1.o chancgmp.o chanfunc.o changaba.o changly.o chanif.o chank1.o \
chank3.o chank4.o chank5.o chank10.o chankca.o chanclca.o channa1.o channa2.o channa3.o channa5.o \
channa6.o channmda.o chansyn.o drand.o gausnn.o initchan.o makcomp.o \
modcode.o ncconvert.o ncdisp.o ncfuncs.o ncmak.o ncomp.o ncpacks.o \
ncplot.o ncrot.o ncsave.o ncstim.o ncstimfuncs.o ncsub.o ncupacks.o \
plotmod.o prcomp.o randtest.o simann.o stimdum.o stimsub.o synfilt.o  \
wave.o: ncomp.h

chanampa.o chanca1.o chancgmp.o chanfunc.o changaba.o changly.o chanif.o chank1.o \
chank3.o chank4.o chank5.o chank10.o chankca.o chanclca.o channa1.o channa3.o channa5.o channa6.o \
channmda.o chansyn.o init.o initchan.o math.o modcode.o ncbnl.o \
ncconvert.o ncmak.o ncplot.o ncstim.o ncstimfuncs.o ncstimin.o ncsub.o \
prcomp.o rndev.o stimsub.o: nconst.h

code.o interpdum.o modcode.o ncplot.o ncstimfuncs.o ncsub.o plotmod.o: ncplot.h

binom.o chanampa.o chanca1.o chancgmp.o chanfunc.o changaba.o changly.o chanif.o \
chank1.o chank3.o chank4.o chank5.o chank10.o chankca.o chanclca.o channa1.o channa2.o channa3.o \
channa5.o channa6.o channmda.o chansyn.o convarr.o drand.o gausnn.o \
init.o initchan.o makcomp.o math.o modcode.o ncconvert.o ncdisp.o \
ncelem.o ncfuncs.o ncmak.o ncnode.o ncomp.o ncplot.o ncray.o ncrot.o \
ncsave.o ncsetvar.o ncstim.o ncstimfuncs.o ncstimin.o ncsub.o plotmod.o \
prcomp.o random.o randtest.o rnd_glibc.o rndev.o simann.o stimdum.o \
stimsub.o synfilt.o wave.o: ncsub.h

binom.o chankca.o chanclca.o code.o convarr.o defaultdraw.o drand.o emalloc.o \
gausnn.o gaussfitn.o gmain.o gprim.o init.o initchan.o \
lm_funcs.o makcomp.o makgauss.o math.o matrix.o modcode.o ncconvert.o \
ncdisp.o ncelem.o ncfuncs.o ncm.o ncmain.o ncmak.o ncnode.o \
ncomp.o ncplot.o ncray.o ncrot.o ncsave.o ncsetvar.o ncstim.o ncstimfuncs.o \
ncsub.o ncsymb.o ncupacks.o ncvirt.o plotmod.o prcomp.o rndev.o stimsub.o \
symbol.o synfilt.o vecsub.o: ncio.h

init.o ncmain.o nctest.o plotmod.o test.o: ncval.h

ncpacks.o ncslv.o ncupacks.o ncvirt.o: ncvirt.h

code.o frame.o gprim.o gr.o init.o modcode.o ncdisp.o ncfuncs.o ncplot.o \
ncstimfuncs.o ncsymb.o plotmod.o rndpnt.o: ../libP/stdplt.h

convarr.o modcode.o stimsub.o: stim.h

gprim.o: extdraw.h

wave.o: turtle.h goldfish.h

modcode.o vecsub.o: vec.h

ncstim.o stimsub.o: wave.h

code.o gaussfit.o gaussfitn.o lm_eval.o lm_evali.o lm_funcs.o makgauss.o: lm_eval.h

code.o gaussfit.o gaussfitn.o lmmin.o lm_eval.o lm_evali.o lm_funcs.o makgauss.o: lmmin.h

chanampa.o chanca1.o chancgmp.o chanfunc.o changaba.o changly.o chanif.o chank1.o \
chank3.o chank4.o chank5.o chank10.o chankca.o chanclca.o channa1.o channa2.o channa3.o channa5.o \
channa6.o channmda.o chansyn.o code.o init.o initchan.o makcomp.o \
math.o modcode.o ncconvert.o ncdisp.o ncelem.o ncfuncs.o ncm.o ncmain.o \
ncmak.o ncnode.o ncomp.o ncpacks.o ncplot.o ncray.o ncsave.o ncstim.o \
ncstimfuncs.o ncsub.o nctest.o ncupacks.o ncv.o prcomp.o stimsub.o symbol.o \
synfilt.o: y.tab.h

init.o, ncmain.o, ncsetvar.o symbol.o: ncsetvar.h 

x.tab.h:	y.tab.h
	-cmp -s x.tab.h y.tab.h || cp y.tab.h x.tab.h

clean:
	rm -f $(NOBJS) $(SOBJS) wave wave.o y.tab.cc y.tab.h nc libnc.a stim \
	 plotmod plotmod.o notinit.o openzfil.o rndpnt.o gausnn $(GOBJS) core

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

