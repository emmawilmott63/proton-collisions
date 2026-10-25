PYTHIA = $(HOME)/Desktop/projects/pythia8318
FASTJET = $(HOME)/Desktop/projects/fastjet-install

CXX = c++
CXXFLAGS = -O2 -std=c++17 $(shell root-config --cflags) -I$(PYTHIA)/include -I$(FASTJET)/include
LDFLAGS = $(shell root-config --libs) -L$(PYTHIA)/lib -lpythia8 -L$(FASTJET)/lib -lfastjet
RPATH = -Wl,-rpath,$(PYTHIA)/lib -Wl,-rpath,$(FASTJET)/lib

# Generation programs (share generateCommon.h)
GEN = jets-charm jets-inc

# Analysis programs (share analysisUtils.h)
ANA = sanity mass buildtree features

# TMVA programs (share analysisUtils.h and tmvaCommon.h); need -lTMVA, which
# root-config --libs does not include
TMVA_PROGS = trainTMVA applyTMVA

all: $(GEN) $(ANA) $(TMVA_PROGS) reconstruct-d0

$(GEN): %: macros/%.cc macros/generateCommon.h macros/analysisUtils.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)

$(ANA): %: macros/%.cc macros/analysisUtils.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)

$(TMVA_PROGS): %: macros/%.cc macros/analysisUtils.h macros/tmvaCommon.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) -lTMVA $(RPATH)

reconstruct-d0: macros/reconstructD0.cc
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)
	
clean:
	rm -f $(GEN) $(ANA) $(TMVA_PROGS) test reconstruct-d0

.PHONY: all clean
