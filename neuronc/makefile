
NCFIL = nc stim plotmod gausnn 

all:  libP.a vid yacc ncfil 

yacc:	
	cd yaccsrc; make yacc

vid:	
	cd pl; make vid; ln -f vid ../bin

libP.a:
	cd libP; make libP.a 

ncfil: 
	cd src; make $(NCFIL); ln -f $(NCFIL) ../bin

clean:
	cd src; make clean
	cd yaccsrc; make clean
	cd libP; make clean
	cd pl; make clean
	cd bin; rm -f $(NCFIL)

backup:
	cd ..; tar czf nc.tgz nc

