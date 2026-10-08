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

all: $(GEN) $(ANA)

$(GEN): %: macros/%.cc macros/generateCommon.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)

$(ANA): %: macros/%.cc macros/analysisUtils.h
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)

test: macros/test.cc
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)

clean:
	rm -f $(GEN) $(ANA) test

.PHONY: all clean
