STRAIGHTEN_C=src/straighten/straighten.c src/common/util.c src/common/debug.c
STRAIGHTEN_H=src/common/debug.h src/common/util.h src/common/argboiler.h src/straighten/argboiler_args.h 
STRAIGHTEN_O=-o bin/straighten
SVD_C=src/svd/svd.c src/common/util.c src/common/debug.c
SVD_H=src/common/debug.h src/common/util.h src/common/argboiler.h
SVD_O=-o bin/svd
SAMPLER_C=src/sampler/sampler.c src/common/util.c src/common/debug.c
SAMPLER_H=src/common/debug.h src/common/util.h src/common/argboiler.h
SAMPLER_O=-o bin/sampler
DIMENSIONS_C=src/dimensions/dimensions.c src/common/util.c src/common/debug.c
DIMENSIONS_H=src/common/debug.h src/common/util.h src/common/argboiler.h
DIMENSIONS_O=-o bin/dimensions
GCCOPTS=-lg2 -lgd -lm -lX11 -lXext -largtable2 -lgsl -lgslcblas -lrt -lfftw3 -lSDL /usr/lib/svdlibc/libsvd.a
DEBUGOPTS=-DNO_PROGRESS_BARS -DGSL_RANGE_CHECK_OFF -g
OPTOPTS=-mtune=native -march=native -O3

all: bin/straighten

bin/straighten: Makefile $(STRAIGHTEN_C) src/straighten/initial_guess.c src/straighten/refine_backbone.c src/straighten/restack.c $(STRAIGHTEN_H)
	gcc $(OPTOPTS) $(STRAIGHTEN_C) $(GCCOPTS) $(STRAIGHTEN_O)

bin/dimensions: Makefile $(DIMENSIONS_C) $(DIMENSIONS_H)
	gcc $(OPTOPTS) $(DIMENSIONS_C) $(GCCOPTS) $(DIMENSIONS_O)

bin/svd: Makefile $(SVD_C) $(SVD_H)
	gcc $(OPTOPTS) $(SVD_C) $(GCCOPTS) $(SVD_O)

bin/sampler: Makefile $(SAMPLER_C) $(SAMPLER_H)
	gcc $(OPTOPTS) $(SAMPLER_C) $(GCCOPTS) $(SAMPLER_O)

.PHONY : noprogress
noprogress: Makefile $(STRAIGHTEN_C) $(STRAIGHTEN_H)
	gcc -DNO_PROGRESS_BARS $(OPTOPTS) $(STRAIGHTEN_C) $(GCCOPTS) $(STRAIGHTEN_O)

.PHONY : debug_straighten
debug_straighten:
	gcc $(DEBUGOPTS) $(STRAIGHTEN_C) $(GCCOPTS) $(STRAIGHTEN_O)
	echo "run $(RUNOPTS)" > gdb.tmp
	gdb -x gdb.tmp bin/straighten
	rm gdb.tmp

.PHONY : ddd_straighten
ddd_straighten:
	gcc $(DEBUGOPTS) $(STRAIGHTEN_C) $(GCCOPTS) $(STRAIGHTEN_O)
	ddd --gdb bin/straighten

.PHONY : profile_straighten
profile_straighten:
	gcc $(OPTOPTS) -g -pg $(STRAIGHTEN_C) $(GCCOPTS) $(STRAIGHTEN_O)
	bin/straighten $(RUNOPTS)
	gprof bin/straighten gmon.out | less
	rm gmon.out

.PHONY : clean
clean:
	rm -f gdb.tmp bin/straighten
