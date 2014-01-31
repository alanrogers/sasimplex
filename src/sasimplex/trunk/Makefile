# Where the executable files will be copied
destination := $(HOME)/bin

opt := -DNDEBUG -O3  -finline-functions  # For full optimization
#opt :=  -O0 -fno-inline-functions -DDEBUG     # For debugging
#prof := -pg -rdynamic                    # For profiling
prof :=
incl := -I/usr/local/include -I/opt/local/include

CC := gcc

# Flags to determine the warning messages issued by the compiler
warn := \
 -Wall \
 -Wcast-align \
 -Wcast-qual \
 -Wmissing-declarations \
 -Wmissing-prototypes \
 -Wnested-externs \
 -Wpointer-arith \
 -Wstrict-prototypes \
 -Wno-unused-parameter \
 -Wno-unused-function \
 -Wshadow \
 -Wundef \
 -Wwrite-strings

CFLAGS := -g -std=gnu99 $(warn) $(incl) $(opt) $(prof) $(osargs)

lib := -L/usr/local/lib -lgsl -lgslcblas -lpthread -lm

.c.o:
	$(CC) $(CFLAGS) $(incl) -c -o ${@F}  $<

XSASIMPLEX := xsasimplex.o sasimplex.o
xsasimplex : $(XSASIMPLEX)
	$(CC) $(CFLAGS) -o $@ $(XSASIMPLEX) $(lib)

XSIMPLEX := xsimplex.o
xsimplex : $(XSIMPLEX)
	$(CC) $(CFLAGS) -o $@ $(XSIMPLEX) $(lib)

# Make dependencies file
depend : *.c *.h
	echo '#Automatically generated dependency info' > depend
	$(CC) -MM $(incl) *.c >> depend

clean :
	rm -f *.o *~ gmon.out *.tmp core core.* vgcore.*

#include depend

.SUFFIXES:
.SUFFIXES: .c .o
.PHONY: clean
