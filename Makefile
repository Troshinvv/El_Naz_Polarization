# =======================================
#   Makefile for MpdRoot Analysis team
#     by V. Kireyeu and A. Mudrokh
#             v2017-04-18
#   Modified by E. Nazarova
#
#  This Makefile is for usage with macro 
#  compilation via the system compiler.
#  Done to be compatible with the new version of mpdroot
#  Do before starting script:
#  source /cvmfs/nica.jinr.ru/sw/os/login.sh
#  module add mpddev
#  export MPDROOT=/scratch2/nazarova/mpd
#  source $MPDROOT/config/env.sh
# =======================================

# Used compiler
CC=g++

# Makefile will proceed all files with .cc suffix
# SOURCES=$(wildcard src/*.cc)
# Output executable has the same name without suffix
# OBJECTS=$(patsubst src/%.cc,%, $(SOURCES))

SOURCES=$(wildcard *.cc)
OBJECTS=$(patsubst %.cc,%, $(SOURCES))

ROOTCONFIG := root-config

CFLAGS := $(shell $(ROOTCONFIG) --cflags)
CFLAGS += -I${FAIRROOT_ROOT}/include
CFLAGS += -I$(FAIRROOT_ROOT)/base/event 
CFLAGS += -I$(FAIRLOGGER_ROOT)/include
CFLAGS += -I$(FMT_ROOT)/include
CFLAGS += -I$(PYTHIA6_ROOT)
CFLAGS += -I$(VMC_ROOT)/include/vmc
CFLAGS += -I$(MPDROOT)/include
CFLAGS += -I$(PWD)/

LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LDFLAGS += $(shell $(ROOTCONFIG) --glibs)
LDFLAGS += -L$(FAIRROOT_ROOT)/lib
LDFLAGS += -L$(PYTHIA6_ROOT)/lib
LDFLAGS += -L$(MPDROOT)/lib
LDFLAGS += -L$(PWD)/ 

LDFLAGS += -lMpdMCStack -lPassive -lBase -lMpdField -lUniGenFormat -lMpdDst -lMpdMcDst -lLHETrack
LDFLAGS += -lMpdBase -lMpdGeneralGenerator -lMpdGen -lMpdMiniEvent
LDFLAGS += -lKalman -ltpc -lTof -lEtof -lMinuit -lpythia6 -lEGPythia6 -lMathMore
LDFLAGS += -lZdc -lFfd -lLHETrack -lSts -lMpdGen -lMpdPid -lEG -lMpdMiniEvent -g -O0 -std=c++17

all:	$(OBJECTS)

# $(OBJECTS): % : src/%.cc
$(OBJECTS): % : %.cc
	$(CC) $< -o $@ $(CFLAGS) $(LDFLAGS)

clean:
	rm -vf $(OBJECTS)

