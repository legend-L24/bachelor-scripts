CPPFLAGS = -DBENCHMARK
CFLAGS = -pipe -g -std=c99 -Wall -march=native -O2
OBJECTS = bench_rbtree bench_metric bench_bump bench_cryst

all: $(OBJECTS)

bench_rbtree: rbtree.o bench_rbtree.o
bench_metric: rng.o metric.o utils.o bench_metric.o -lm
bench_bump: rng.o rbtree.o metric.o utils.o \
	cryst_base.o cryst_extra.o cryst_read.o bench_bump.o -lm
bench_cryst: rng.o rbtree.o metric.o utils.o \
	cryst_base.o cryst_extra.o cryst_read.o bench_cryst.o -lm

clean:
	rm -f *.o

distclean: clean
	rm -f $(OBJECTS)

