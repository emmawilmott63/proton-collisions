PYTHIA = $(HOME)/Desktop/projects/pythia8318
FASTJET = $(HOME)/Desktop/projects/fastjet-install

CXX = c++
CXXFLAGS = -O2 -std=c++17 $(shell root-config --cflags) -I$(PYTHIA)/include -I$(FASTJET)/include
LDFLAGS = $(shell root-config --libs) -L$(PYTHIA)/lib -lpythia8 -L$(FASTJET)/lib -lfastjet
RPATH = -Wl,-rpath,$(PYTHIA)/lib -Wl,-rpath,$(FASTJET)/lib

test: macros/test.cc
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)

generateD0: macros/generateD0.cc
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)

reconstructD0: macros/reconstructD0.cc
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)

sanityChecks: macros/sanityChecks.cc
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(RPATH)

clean:
	rm -f generateD0 reconstructD0 sanityChecks test
