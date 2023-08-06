#I change O2 to O1,if process ,should return back
CFLAGS = -pipe -g -std=c99 -Wall -march=native -O2
OBJECTS = decr_bias decr_grid decr_mcs decr_lsa decr_sac decr_sas
#OBJECTS = decr_mcs

all: $(OBJECTS)
decr_bias: rng.o rbtree.o metric.o utils.o cryst_base.o cryst_extra.o \
        cryst_read.o optim.o decr_bias.o -lm
decr_grid: rng.o rbtree.o metric.o utils.o cryst_base.o cryst_extra.o \
        cryst_read.o optim.o decr_grid.o -lm
decr_mcs: rng.o rbtree.o metric.o utils.o cryst_base.o cryst_extra.o \
	cryst_read.o optim.o decr_mcs.o -lm
decr_lsa: rng.o rbtree.o metric.o utils.o cryst_base.o cryst_extra.o \
	cryst_read.o optim.o decr_lsa.o -lm
decr_sac: rng.o rbtree.o metric.o utils.o cryst_base.o cryst_extra.o \
	cryst_read.o optim.o netio.o decr_sac.o -lm -lzmq
decr_sas: rng.o rbtree.o metric.o utils.o cryst_base.o cryst_extra.o \
	cryst_read.o optim.o netio.o decr_sas.o -lm -lpthread -lzmq

clean:
	rm -f *.o

distclean: clean
	rm -f $(OBJECTS)

