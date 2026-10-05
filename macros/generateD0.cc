#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"

#include <vector>

int main() {

    // -----------------------------
    // Pythia setup
    // -----------------------------
    Pythia8::Pythia p;

    p.readString("HardQCD:all = on");
    p.readString("Beams:eCM = 200.");
    p.readString("PhaseSpace:pTHatMin = 5.");
    p.readString("PhaseSpace:pTHatMax = 10.");

    p.init();

    // -----------------------------
    // ROOT output
    // -----------------------------
    TFile out("data/d0.root", "RECREATE");

    TTree tree("Particles", "Final-state particles");

    // Event information
    int event;

    // Particle information
    std::vector<int> pid;
    std::vector<int> charge;
    std::vector<float> px;
    std::vector<float> py;
    std::vector<float> pz;
    std::vector<float> energy;

    // Truth information:
    // index of the D0 parent in the Pythia event record.
    // -1 means this particle was not identified as
    // a daughter of a D0 -> K pi decay.
    std::vector<int> d0Parent;

    tree.Branch("event", &event);

    tree.Branch("pid", &pid);
    tree.Branch("charge", &charge);
    tree.Branch("px", &px);
    tree.Branch("py", &py);
    tree.Branch("pz", &pz);
    tree.Branch("energy", &energy);

    tree.Branch("d0Parent", &d0Parent);

    // -----------------------------
    // Event loop
    // -----------------------------
    const int nEvents = 10000;

    for (event = 0; event < nEvents; ++event) {

        if (!p.next())
            continue;

        // Clear vectors for this event
        pid.clear();
        charge.clear();
        px.clear();
        py.clear();
        pz.clear();
        energy.clear();
        d0Parent.clear();

        // -------------------------------------------------
        // Find D0 -> K pi decays in the Pythia event record
        // -------------------------------------------------

        std::vector<int> d0DaughterIndices;
        std::vector<int> d0ParentIndices;

        for (int i = 0; i < p.event.size(); ++i) {

            int id = p.event[i].id();

            // D0 or anti-D0
            if (id != 421 && id != -421)
                continue;

            int d1 = p.event[i].daughter1();
            int d2 = p.event[i].daughter2();

            if (d1 <= 0 || d2 <= 0)
                continue;

            // We are looking for the direct two-body decay
            // D0  -> K- pi+
            // D0bar -> K+ pi-
            if (d2 != d1 + 1)
                continue;

            int id1 = p.event[d1].id();
            int id2 = p.event[d2].id();

            bool isD0ToKPi =
                (id == 421 &&
                 id1 == -321 &&
                 id2 == 211);

            bool isD0barToKPi =
                (id == -421 &&
                 id1 == 321 &&
                 id2 == -211);

            if (!isD0ToKPi && !isD0barToKPi)
                continue;

            // Store the two daughter indices.
            d0DaughterIndices.push_back(d1);
            d0DaughterIndices.push_back(d2);

            d0ParentIndices.push_back(i);
            d0ParentIndices.push_back(i);
        }

        // -------------------------------------------------
        // Store all final-state particles
        // -------------------------------------------------

        for (int i = 0; i < p.event.size(); ++i) {

            if (!p.event[i].isFinal())
                continue;

            pid.push_back(p.event[i].id());
            charge.push_back(static_cast<int>(p.event[i].charge()));

            px.push_back(static_cast<float>(p.event[i].px()));
            py.push_back(static_cast<float>(p.event[i].py()));
            pz.push_back(static_cast<float>(p.event[i].pz()));
            energy.push_back(static_cast<float>(p.event[i].e()));

            // Default: not from a reconstructed D0 -> K pi decay
            int parent = -1;

            // Check whether this particle is one of our
            // truth D0 daughters.
            for (size_t j = 0; j < d0DaughterIndices.size(); ++j) {

                if (i == d0DaughterIndices[j]) {
                    parent = d0ParentIndices[j];
                    break;
                }
            }

            d0Parent.push_back(parent);
        }

        tree.Fill();
    }

    tree.Write();
    out.Close();

    return 0;
}