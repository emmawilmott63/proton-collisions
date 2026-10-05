#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "TFile.h"

int main() {
    Pythia8::Pythia p;
    p.readString("HardQCD:all = on");
    p.readString("Beams:eCM = 200.");
    p.init();
    p.next();

    fastjet::JetDefinition jets(fastjet::antikt_algorithm, 0.4);
    TFile f("data/test.root", "RECREATE");

    return 0;
}
