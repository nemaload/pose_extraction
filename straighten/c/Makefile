CFILES=straighten.c util.c debug.c
GCCOPTS=-lg2 -lgd -lm -lX11 -lXext -largtable2 -lgsl -lgslcblas -lrt -o straighten
RUNOPTS=../TU3335.raw --width 2437 --height 1851 -p
DEBUGOPTS=-DNO_PROGRESS_BARS -DGSL_RANGE_CHECK_OFF -g

all: straighten run

straighten: debug.h Makefile $(CFILES)
	gcc -O3 $(CFILES) $(GCCOPTS)

.PHONY : run
run: straighten
	./straighten $(RUNOPTS)

.PHONY : debug
debug:
	gcc $(DEBUGOPTS) $(CFILES) $(GCCOPTS)
	echo "run $(RUNOPTS)" > gdb.tmp
	gdb -x gdb.tmp ./straighten
	rm gdb.tmp

.PHONY : ddd
ddd:
	gcc $(DEBUGOPTS) $(CFILES) $(GCCOPTS)
	ddd --gdb ./straighten

.PHONY : profile
profile:
	gcc -O3 -g -pg $(CFILES) $(GCCOPTS)
	./straighten $(RUNOPTS)
	gprof straighten gmon.out | less
	rm gmon.out

.PHONY : clean
clean:
	rm -f gdb.tmp straighten
