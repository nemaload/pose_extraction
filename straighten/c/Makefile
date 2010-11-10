CFILES=straighten.c util.c debug.c
GCCOPTS=-O3 -lg2 -lgd -lm -lX11 -lXext -largtable2 -lgsl -lgslcblas -lrt -o straighten
RUNOPTS=../TU3335.raw --width 2437 --height 1851 -pp

all: straighten run

straighten: debug.h Makefile $(CFILES)
	gcc $(CFILES) $(GCCOPTS)

.PHONY : run
run: straighten
	./straighten $(RUNOPTS)

.PHONY : debug
debug:
	gcc -g $(CFILES) $(GCCOPTS)
	echo "run $(RUNOPTS)" > gdb.tmp
	gdb -x gdb.tmp ./straighten
	rm gdb.tmp

.PHONY : profile
profile:
	gcc -g -pg $(CFILES) $(GCCOPTS)
	./straighten $(RUNOPTS)
	gprof straighten gmon.out | less

.PHONY : clean
clean:
	rm -f gdb.tmp straighten
