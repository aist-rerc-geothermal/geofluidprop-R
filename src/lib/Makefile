# compile just the C-callable library
library:
	( $(MAKE) lib )

# remove object files, but keep the compiled programs and library archives
# clean:
# 	( $(MAKE) clean )

# CC		= gcc
COPTS		= -O2
#WARN        = -Wall -Wno-unknown-pragmas
#WARN        = -Wno-unknown-pragmas
INCDIR      = 
INCDIR      = -I./geofluidprop/src
CFLAGS		= $(COPTS) $(WARN) $(INCDIR)  -std=gnu99
LLIBS		= -lm

SOURCES = $(wildcard ./geofluidprop/*.c ./geofluidprop/src/util/*.c ./geofluidprop/src/model/*/*.c)
SOURCES := $(filter-out $(wildcard ./geofluidprop/src/model/*/*_prost*), $(SOURCES))
#SOURCES := $(filter-out $(wildcard src/**/test*.c src/**/**/test*.c),$(SOURCES))
#SOURCES := $(filter-out $(wildcard src/interface/*.* src/spline/*.* src/sbtl/*.*),$(SOURCES))
$(warning SOURCES = $(SOURCES))
OBJS = $(SOURCES:.c=.o)

#$(info        PWD = $(shell pwd)  )

LIB = libgeofluidprop.a

lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -rcv $@ $(OBJS)

clean:
	-rm -f $(OBJS)
#	-rm -f src/util/*.o *?~ core
	-rm -f $(LIB)

