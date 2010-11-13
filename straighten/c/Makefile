CFILES=straighten.c util.c debug.c graphics.c argboiler.c
GCCOPTS=-lg2 -lgd -lm -lX11 -lXext -largtable2 -lgsl -lgslcblas -lrt -o straighten
RUNOPTS=../TU3335.raw --width 2437 --height 1851 -p #--out-slice=80
DEBUGOPTS=-DNO_PROGRESS_BARS -DGSL_RANGE_CHECK_OFF -g
OPTOPTS=-mtune=native -march=native -O3

all: straighten run

straighten: debug.h util.h graphics.h argboiler.h argboiler_args.h Makefile $(CFILES)
	gcc $(OPTOPTS) $(CFILES) $(GCCOPTS)

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
	gcc $(OPTOPTS) -g -pg $(CFILES) $(GCCOPTS)
	./straighten $(RUNOPTS)
	gprof straighten gmon.out | less
	rm gmon.out

.PHONY : clean
clean:
	rm -f gdb.tmp straighten
