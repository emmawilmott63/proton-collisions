#ifndef GENERATE_COMMON_H
#define GENERATE_COMMON_H

// Shared Pythia -> ROOT generation code.
// Used by jets-charm.cc and jets-inc.cc so the two
// samples are produced by IDENTICAL code and differ only in the
// hard-process selection.

#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

enum class SampleMode {
    Charm,      // every event has a c cbar hard process
    Inclusive   // all hard QCD processes: realistic, rare charm
};

inline int generateSample(SampleMode mode,
                          int nEvents,
                          int seed,
                          const std::string& outName) {

    const std::string modeName =
        (mode == SampleMode::Charm) ? "charm" : "inclusive";

    gSystem->mkdir("data", true);

    // -----------------------------
    // Pythia setup
    // -----------------------------
    Pythia8::Pythia p;

    if (mode == SampleMode::Charm) {
        // Charm-enriched: only g g -> c cbar and q qbar -> c cbar.
        // NOTE: this excludes flavour excitation and gluon splitting
        // (g -> c cbar in the shower), which are in HardQCD:all.
        p.readString("HardQCD:hardccbar = on");
    } else {
        p.readString("HardQCD:all = on");
    }

    p.readString("Beams:eCM = 200.");
    p.readString("PhaseSpace:pTHatMin = 5.");
    p.readString("PhaseSpace:pTHatMax = 10.");

    p.readString("Random:setSeed = on");
    p.readString("Random:seed = " + std::to_string(seed));

    // Keep the console quiet
    p.readString("Next:numberShowEvent = 0");

    if (!p.init())
        return 1;

    // -----------------------------
    // ROOT output
    // -----------------------------
    TFile out(outName.c_str(), "RECREATE");

    TTree tree("Particles", "Final-state particles");

    // Event information
    int event;

    // Final-state particle information
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

    // Truth information about ALL D0 / D0bar in the event
    // (one entry per D0, whatever it decays to).
    std::vector<int>   d0TruthId;     // +421 or -421
    std::vector<float> d0TruthPt;
    std::vector<float> d0TruthY;      // rapidity
    std::vector<int>   d0TruthToKPi;  // 1 if direct D0 -> K pi decay

    tree.Branch("event", &event);

    tree.Branch("pid", &pid);
    tree.Branch("charge", &charge);
    tree.Branch("px", &px);
    tree.Branch("py", &py);
    tree.Branch("pz", &pz);
    tree.Branch("energy", &energy);

    tree.Branch("d0Parent", &d0Parent);

    tree.Branch("d0TruthId", &d0TruthId);
    tree.Branch("d0TruthPt", &d0TruthPt);
    tree.Branch("d0TruthY", &d0TruthY);
    tree.Branch("d0TruthToKPi", &d0TruthToKPi);

    // Bookkeeping for consistency checks
    long long nGenerated       = 0;
    long long nD0All           = 0;
    long long nD0KPi           = 0;
    long long nDaughtersStored = 0;
    long long nEventsWithD0KPi = 0;

    // -----------------------------
    // Event loop
    // -----------------------------
    for (event = 0; event < nEvents; ++event) {

        if (!p.next())
            continue;

        ++nGenerated;

        // Clear vectors for this event
        pid.clear();
        charge.clear();
        px.clear();
        py.clear();
        pz.clear();
        energy.clear();
        d0Parent.clear();

        d0TruthId.clear();
        d0TruthPt.clear();
        d0TruthY.clear();
        d0TruthToKPi.clear();

        // -------------------------------------------------
        // Find D0 decays in the Pythia event record
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

            // Skip intermediate copies of the same particle
            // so each physical D0 is counted once.
            if (d1 > 0 && p.event[d1].id() == id)
                continue;

            // We are looking for the direct two-body decay
            // D0    -> K- pi+
            // D0bar -> K+ pi-
            bool isKPi = false;

            if (d1 > 0 && d2 > 0 && d2 == d1 + 1) {

                int id1 = p.event[d1].id();
                int id2 = p.event[d2].id();

                isKPi =
                    (id ==  421 && id1 == -321 && id2 ==  211) ||
                    (id == -421 && id1 ==  321 && id2 == -211);
            }

            // Truth record for every D0
            d0TruthId.push_back(id);
            d0TruthPt.push_back(static_cast<float>(p.event[i].pT()));
            d0TruthY.push_back(static_cast<float>(p.event[i].y()));
            d0TruthToKPi.push_back(isKPi ? 1 : 0);

            ++nD0All;

            if (!isKPi)
                continue;

            ++nD0KPi;

            // Store the two daughter indices.
            d0DaughterIndices.push_back(d1);
            d0DaughterIndices.push_back(d2);

            d0ParentIndices.push_back(i);
            d0ParentIndices.push_back(i);
        }

        if (!d0ParentIndices.empty())
            ++nEventsWithD0KPi;

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

            // Default: not from a D0 -> K pi decay
            int parent = -1;

            // Check whether this particle is one of our
            // truth D0 daughters.
            for (size_t j = 0; j < d0DaughterIndices.size(); ++j) {

                if (i == d0DaughterIndices[j]) {
                    parent = d0ParentIndices[j];
                    ++nDaughtersStored;
                    break;
                }
            }

            d0Parent.push_back(parent);
        }

        tree.Fill();
    }

    tree.Write();
    out.Close();

    // -----------------------------
    // Summary and consistency checks
    // -----------------------------
    const double nGen = std::max<double>(nGenerated, 1.0);

    std::cout << "\n===== Generation summary =====\n"
              << "Sample:                     " << modeName << "\n"
              << "Seed:                       " << seed << "\n"
              << "Events generated:           " << nGenerated << "\n"
              << "Cross section (sigmaGen):   " << p.info.sigmaGen()
              << " mb  (+/- " << p.info.sigmaErr() << ")\n"
              << "D0 + D0bar in events:       " << nD0All
              << "  (" << nD0All / nGen << " per event)\n"
              << "Direct D0 -> K pi decays:   " << nD0KPi
              << "  (fraction of all D0: "
              << double(nD0KPi) / std::max(nD0All, 1LL) << ")\n"
              << "Events with a D0 -> K pi:   " << nEventsWithD0KPi
              << "  (" << 100.0 * nEventsWithD0KPi / nGen << " %)\n"
              << "K/pi daughters stored:      " << nDaughtersStored
              << "  (expected " << 2 * nD0KPi << ")\n";

    if (nDaughtersStored != 2 * nD0KPi)
        std::cout << "WARNING: some D0 -> K pi daughters are not final-state "
                     "particles, so truth matching will miss them.\n";

    std::cout << "Output written to " << outName << "\n"
              << "==============================\n";

    return 0;
}

#endif
