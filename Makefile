STRAIGHTEN_C=src/straighten/straighten.c src/common/util.c src/common/debug.c
STRAIGHTEN_H=src/common/debug.h src/common/util.h src/common/argboiler.h src/straighten/argboiler_args.h 
STRAIGHTEN_O=-o bin/straighten
GCCOPTS=-lg2 -lgd -lm -lX11 -lXext -largtable2 -lgsl -lgslcblas -lrt
RUNOPTS=images/TU3335.raw --width 2437 --height 1851 -p -n4 --out-slice=80 -s 2.0
DEBUGOPTS=-DNO_PROGRESS_BARS -DGSL_RANGE_CHECK_OFF -g
OPTOPTS=-mtune=native -march=native -O3

all: bin/straighten run

bin/straighten: Makefile $(STRAIGHTEN_C) $(STRAIGHTEN_H)
	gcc $(OPTOPTS) $(STRAIGHTEN_C) $(GCCOPTS) $(STRAIGHTEN_O)

.PHONY : noprogress
noprogress: Makefile $(STRAIGHTEN_C) $(STRAIGHTEN_H)
	gcc -DNO_PROGRESS_BARS $(OPTOPTS) $(STRAIGHTEN_C) $(GCCOPTS) $(STRAIGHTEN_O)

.PHONY : run
run: bin/straighten
	bin/straighten $(RUNOPTS)

.PHONY : debug
debug:
	gcc $(DEBUGOPTS) $(STRAIGHTEN_C) $(GCCOPTS) $(STRAIGHTEN_O)
	echo "run $(RUNOPTS)" > gdb.tmp
	gdb -x gdb.tmp bin/straighten
	rm gdb.tmp

.PHONY : ddd
ddd:
	gcc $(DEBUGOPTS) $(STRAIGHTEN_C) $(GCCOPTS) $(STRAIGHTEN_O)
	ddd --gdb bin/straighten

.PHONY : profile
profile:
	gcc $(OPTOPTS) -g -pg $(STRAIGHTEN_C) $(GCCOPTS) $(STRAIGHTEN_O)
	bin/straighten $(RUNOPTS)
	gprof bin/straighten gmon.out | less
	rm gmon.out

.PHONY : clean
clean:
	rm -f gdb.tmp bin/straighten
