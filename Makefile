PYTHIA = $(HOME)/Desktop/projects/pythia8318
FASTJET = $(HOME)/Desktop/projects/fastjet-install

CXX = c++
CXXFLAGS = -O2 -std=c++17 $(shell root-config --cflags) \
           -I$(PYTHIA)/include -I$(FASTJET)/include
LDFLAGS = $(shell root-config --libs) \
          -L$(PYTHIA)/lib -lpythia8 \
          -L$(FASTJET)/lib -lfastjet

generateD0: macros/generateD0.cc
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -f generateD0
