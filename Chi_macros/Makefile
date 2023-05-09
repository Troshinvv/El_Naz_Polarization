# =======================================
#   Makefile for MpdRoot Analysis team
#     by V. Kireyeu and A. Mudrokh
#             v2017-04-18
#
#  This Makefile uses for analysis script 
#  compilation via the system compiler.
#  Please do not forget to run env.sh
#  script from MpdRoot build directory!!! 
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
CFLAGS += -I$(FAIRLOGGER_ROOT)/include
CFLAGS += -I$(FMT_ROOT)/include
CFLAGS += -I$(PYTHIA6_ROOT)
CFLAGS += -I$(VMC_ROOT)/include/vmc
CFLAGS += -I$(MPDROOT)/include

LDFLAGS := $(shell $(ROOTCONFIG) --ldflags) 
LDFLAGS += $(shell $(ROOTCONFIG) --glibs)
LDFLAGS += -L$(FAIRROOT_ROOT)/lib
LDFLAGS += -L$(PYTHIA6_ROOT)/lib
LDFLAGS += -L$(MPDROOT)/lib

LDFLAGS += -lMpdMCStack -lPassive -lBase -lMpdField -lUniGenFormat -lMpdDst -lMpdMcDst
LDFLAGS += -lMpdBase -lMpdGeneralGenerator -lMpdGen -lMpdMiniEvent -lMpdPhysics -ltpcClusterHitFinder
LDFLAGS += -lKalman -ltpc -lTof -lEtof -lMinuit -lpythia6 -lEGPythia6 -lMathMore -ltpcDigitizer
LDFLAGS += -lZdc -lFfd -lLHETrack -lSts -lMpdPid -lEG -lMpdMiniEvent -g -O0 -std=c++17

#LDFLAGS += -lBase -lDbUtils -lEmc -lEtof -lEventDisplay -lFfd -lKalman -lLHETrack -lMpdBase
#LDFLAGS += -lMpdCentralityAll -lMpdDst -lMpdEventPlaneAll -lMinuit -lpythia6 -lEGPythia6 -lMathMore -lEG
#LDFLAGS += -lMpdField -lMpdFsiTools -lMpdGeneralGenerator -lMpdGenFactory -lMpdGen -lMpdMcDst
#LDFLAGS += -lMpdMCStack -lMpdMiniEvent -lMpdPairKKTracks -lMpdPhotons -lMpdPhysics -lMpdPid
#LDFLAGS += -lmulti -lPassive -lQA -lSts -lTof -ltpcAlignment -ltpcClusterHitFinder -ltpcDigitizer
#LDFLAGS += -ltpcFairTpc -ltpcGeometry -ltpcPid -ltpc -lUniCommon -lUniGenFormat -lZdc  -g -O0 -std=c++17

all:	$(OBJECTS)

# $(OBJECTS): % : src/%.cc
$(OBJECTS): % : %.cc
	$(CC) $< -o $@ $(CFLAGS) $(LDFLAGS) 

clean:
	rm -vf $(OBJECTS)

