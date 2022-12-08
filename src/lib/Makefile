# compile just the C-callable library
library:
	( $(MAKE) lib )

# remove object files, but keep the compiled programs and library archives
# clean:
# 	( $(MAKE) clean )

# CC		= gcc
#COPTS		= -O2
#WARN        = -Wall -Wno-unknown-pragmas
#WARN        = -Wno-unknown-pragmas
INCDIR      = -I./geofluidprop/src
CFLAGS		= $(COPTS) $(WARN) $(INCDIR)  -std=gnu99
#LLIBS		= -lm

SOURCES = \
	geofluidprop/src/util/utility.c\
	geofluidprop/src/util/gauss3.c\
	geofluidprop/src/model/iapws/IAPS-84.c\
	geofluidprop/src/model/iapws/IAPWS-IF97-region1.c\
	geofluidprop/src/model/iapws/IAPWS-IF97-region2.c\
	geofluidprop/src/model/iapws/IAPWS-IF97-region3.c\
	geofluidprop/src/model/iapws/IAPWS-IF97-boundary23.c\
	geofluidprop/src/model/iapws/IAPWS-IF97-sat.c\
	geofluidprop/src/model/iapws/IAPWS-IF97.c\
	geofluidprop/src/model/iapws/IAPWS-95.c\
	geofluidprop/src/model/iapws/IAPWS-Melt-11.c\
	geofluidprop/src/model/iapws/IAPWS-SAT-92.c\
	geofluidprop/src/model/iapws/IAPWS-ThermalConductivity-11.c\
	geofluidprop/src/model/iapws/IAPWS-Viscosity-85.c\
	geofluidprop/src/model/iapws/IAPWS-Viscosity-08.c\
	geofluidprop/src/model/driesner07/Driesner2007_H2ONaCl.c\
	geofluidprop/src/model/driesner07/Driesner2007_NaCl.c\
	geofluidprop/src/model/klyukinetal17/KlyukinEtAl2017.c

#SOURCES = $(wildcard ./geofluidprop/*.c ./geofluidprop/src/util/*.c ./geofluidprop/src/model/*/*.c)
#SOURCES := $(filter-out $(wildcard ./geofluidprop/src/model/*/*_prost*), $(SOURCES))
#SOURCES := $(filter-out $(wildcard src/**/test*.c src/**/**/test*.c),$(SOURCES))
#SOURCES := $(filter-out $(wildcard src/interface/*.* src/spline/*.* src/sbtl/*.*),$(SOURCES))
#$(warning SOURCES = $(SOURCES))
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

