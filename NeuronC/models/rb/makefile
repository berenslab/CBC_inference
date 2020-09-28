# CFLAGS = -O3
# CFLAGS = -O3 -fPIC
CFLAGS = -g 
# CFLAGS = -g -fPIC
# CFLAGS = -pg -fPIC 
# CFLAGS = -pg 
# CFLAGS = -g -Wall -Wmissing-prototypes
# CFLAGS = -O5 -DCPML -ffast-math -fomit-frame-pointer -finline-functions
# CFLAGS = -O5 -ffast-math -fomit-frame-pointer -finline-functions
# CFLAGS = -O3 -non_shared  # for DEC UNIX
# LDFLAGS = -s
LDFLAGS = -g
# LDFLAGS = -pg
# LDFLAGS = -s -non_shared  # for DEC UNIX

## NB: to get the expt*.so files to load, must define path.
## Either use ldconfig, or define LD_LIBRARY_PATH:
##
## bash:    LD_LIBRARY_PATH = .; export LD_LIBRARY_PATH
## tcsh:    setenv LD_LIBRARY_PATH .
##

NC_HOME = ~/nc

ifndef NC_HOME
	NCCDIR=../../src
	LIBPDIR=../../libP
else
	NCCDIR=$(NC_HOME)/src
	LIBPDIR=$(NC_HOME)/libP
endif

#NCCDIR = ../../src
#LIBPDIR = ../../libP

SHELL=/bin/csh
# CXX = g++
# CC = gcc

.SUFFIXES: .cc .so

.cc.o:	
	$(CXX) -c $(CFLAGS) -I $(NCCDIR) $*.cc
.cc.so:
	$(CXX) -shared -Wl,-soname,$*.so -o $*.so $*.cc -I $(NCCDIR) -ldl

LIBNC=libnc.a
NCLIBS = $(NCCDIR)/libnc.a $(LIBPDIR)/libP.a
NCLIBCALL  = -L $(NCCDIR) -lnc -L $(LIBPDIR) -lP

all:	rb3 

printvars:
	@echo "NeuronC Home: $(NC_HOME)"
	@echo "      NCCDIR: $(NCCDIR)"
	@echo "     LIBPDIR: $(LIBPDIR)"

rb3.o: $(NCCDIR)/ncfuncs.h

rb3.o: $(NCCDIR)/ncio.h

$(EXPTOBJS): retcolors.h

rb3:	rb3.o
	$(CXX) $(LDFLAGS) -o rb3 rb3.o $(NCLIBCALL) -lm 

clean:
	rm -f rb3 rb3.o $(RETOBJS) $(EXPTOBJS) retsim libnc.a core

