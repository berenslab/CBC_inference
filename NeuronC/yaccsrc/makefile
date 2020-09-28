CFLAGS= -O
CC=gcc

all:	yacc

cmp:	all
	cmp yacc ${DESTDIR}/usr/bin/yacc

cp:	all
	-strip yacc
	cp yacc ${DESTDIR}/usr/bin/yacc
	rm yacc *.o

yacc: y1.o y2.o y3.o y4.o
	$(CC) $(CFLAGS) -o yacc y1.o y2.o y3.o y4.o
	size yacc

y1.o y2.o y3.o y4.o: dextern files

clean:
	rm -f *.o

backup:
	@echo "Please insert the 'yacc' backup disk; type <RET> when ready..."
	@line
	cd ..; tar cv yacc/{*.c,*.y,*.h,files,dextern,lib,*makefile}
	@echo yacc is backed up...
