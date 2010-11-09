all: straighten run

straighten: straighten.c util.c util.h debug.c debug.h Makefile
	gcc -g -O3 straighten.c util.c debug.c -lg2 -lgd -lm -lX11 -lXext -largtable2 -lgsl -lgslcblas -lrt -o straighten

.PHONY : run
run:
	./straighten ../TU3335.raw --width 2437 --height 1851

.PHONY : debug
debug:
	echo "run ../TU3335.raw --width 2437 --height 1851" > gdb.tmp
	gdb -x gdb.tmp ./straighten
	rm gdb.tmp

.PHONY : clean
clean:
	rm -f gdb.tmp straighten
