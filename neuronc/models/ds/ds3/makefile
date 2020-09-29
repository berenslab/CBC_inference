
DSOBJS= ds3a.n ds3.n dsdefs.n dsseg.n gcdens.n amdens.n \
	dsconv2 dsconv3 \
	makefile howtorun.txt readme.txt nonlin.txt

dsdist: $(DSOBJS)
	mkdir -p dsdist
	cp $(DSOBJS) dsdist
	tar cvzf dsdist.tgz dsdist
	rm -r dsdist

dsbackup: dsdist
	rsyncc dsdist.tgz retina:/home/ftp/pub/rob

