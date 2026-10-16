
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    const std::string inputName =
        (argc > 1) ? argv[1] : "data/d0_inc_5_100.root";

    const std::string outputName =
        (argc > 2) ? argv[2] : "data/d0_reco_5_100.root";

    // --------------------------------------------------
    // Open input file and get the trees
    // --------------------------------------------------
    TFile* input = TFile::Open(inputName.c_str(), "READ");

    if (!input || input->IsZombie()) {
        std::cerr << "Error: could not open " << inputName << '\n';
        return 1;
    }

    TTree* particles = nullptr;
    TTree* runInfo = nullptr;

    input->GetObject("Particles", particles);
    input->GetObject("RunInfo", runInfo);

    if (!particles || !runInfo) {
        std::cerr << "Error: missing Particles or RunInfo tree\n";
        input->Close();
        return 1;
    }

    // --------------------------------------------------
    // Read run metadata and calculate event weights
    //
    // weight = sigmaGen_mb / nGenerated
    // --------------------------------------------------
    int seed = 0;
    long long nGenerated = 0;
    double sigmaGen = 0.0;

    runInfo->SetBranchAddress("seed", &seed);
    runInfo->SetBranchAddress("nGenerated", &nGenerated);
    runInfo->SetBranchAddress("sigmaGen_mb", &sigmaGen);

    std::map<int, double> eventWeights;

    for (Long64_t i = 0; i < runInfo->GetEntries(); ++i) {
        runInfo->GetEntry(i);

        if (nGenerated <= 0 || sigmaGen < 0.0) {
            std::cerr << "Warning: invalid run metadata for seed "
                      << seed << '\n';
            continue;
        }

        eventWeights[seed] = sigmaGen / nGenerated;

        std::cout << "Seed " << seed
                  << ": cross section = " << sigmaGen << " mb"
                  << ", generated events = " << nGenerated
                  << ", event weight = " << eventWeights[seed]
                  << " mb\n";
    }

    // --------------------------------------------------
    // Connect particle-tree branches
    // --------------------------------------------------
    int runSeed = 0;

    std::vector<int>* pid = nullptr;
    std::vector<int>* charge = nullptr;
    std::vector<float>* px = nullptr;
    std::vector<float>* py = nullptr;
    std::vector<float>* pz = nullptr;
    std::vector<float>* energy = nullptr;

    particles->SetBranchAddress("runSeed", &runSeed);
    particles->SetBranchAddress("pid", &pid);
    particles->SetBranchAddress("charge", &charge);
    particles->SetBranchAddress("px", &px);
    particles->SetBranchAddress("py", &py);
    particles->SetBranchAddress("pz", &pz);
    particles->SetBranchAddress("energy", &energy);

    // --------------------------------------------------
    // Output file and histograms
    // --------------------------------------------------
    TFile* output = TFile::Open(outputName.c_str(), "RECREATE");

    if (!output || output->IsZombie()) {
        std::cerr << "Error: could not create " << outputName << '\n';
        input->Close();
        return 1;
    }

    TH1D hMassOS(
        "hMass_OS",
        "Opposite-sign K#pi candidates;M_{K#pi} [GeV];Candidates",
        500, 1.6, 2.1
    );

    TH1D hMassSS(
        "hMass_SS",
        "Same-sign K#pi control;M_{K#pi} [GeV];Candidates",
        500, 1.6, 2.1
    );

    TH1D hMassOSWeighted(
        "hMass_OS_weighted",
        "Weighted opposite-sign K#pi;M_{K#pi} [GeV];Weighted candidates [mb]",
        500, 1.6, 2.1
    );

    TH1D hMassSSWeighted(
        "hMass_SS_weighted",
        "Weighted same-sign K#pi;M_{K#pi} [GeV];Weighted candidates [mb]",
        500, 1.6, 2.1
    );

    hMassOS.Sumw2();
    hMassSS.Sumw2();
    hMassOSWeighted.Sumw2();
    hMassSSWeighted.Sumw2();

    // --------------------------------------------------
    // Reconstruct all K-pi pairs within each event
    // --------------------------------------------------
    Long64_t eventsProcessed = 0;
    Long64_t eventsWithoutWeight = 0;
    Long64_t nOS = 0;
    Long64_t nSS = 0;

    for (Long64_t entry = 0;
         entry < particles->GetEntries();
         ++entry) {

        particles->GetEntry(entry);

        if (!pid || !charge || !px || !py || !pz || !energy) {
            std::cerr << "Error: missing particle data at entry "
                      << entry << '\n';
            continue;
        }

        auto weightIt = eventWeights.find(runSeed);

        if (weightIt == eventWeights.end()) {
            ++eventsWithoutWeight;
            continue;
        }

        const double eventWeight = weightIt->second;
        ++eventsProcessed;

        const size_t nParticles = pid->size();

        if (charge->size() != nParticles ||
            px->size() != nParticles ||
            py->size() != nParticles ||
            pz->size() != nParticles ||
            energy->size() != nParticles) {
            std::cerr << "Warning: mismatched particle-vector sizes at entry "
                      << entry << '\n';
            continue;
        }

        // Find charged kaons (|PDG ID| = 321)
        // and charged pions (|PDG ID| = 211).
        std::vector<size_t> kaons;
        std::vector<size_t> pions;

        for (size_t i = 0; i < nParticles; ++i) {
            if (std::abs(pid->at(i)) == 321 &&
                charge->at(i) != 0) {
                kaons.push_back(i);
            }

            if (std::abs(pid->at(i)) == 211 &&
                charge->at(i) != 0) {
                pions.push_back(i);
            }
        }

        // Combine each kaon with each pion in this event.
        for (size_t k : kaons) {
            for (size_t pi : pions) {
                const double e =
                    energy->at(k) + energy->at(pi);

                const double x =
                    px->at(k) + px->at(pi);

                const double y =
                    py->at(k) + py->at(pi);

                const double z =
                    pz->at(k) + pz->at(pi);

                const double massSquared =
                    e * e - x * x - y * y - z * z;

                // Protect against tiny negative values from rounding.
                if (massSquared < -1e-5) {
                    continue;
                }

                const double mass =
                    std::sqrt(std::max(0.0, massSquared));

                const bool oppositeSign =
                    charge->at(k) * charge->at(pi) < 0;

                if (oppositeSign) {
                    hMassOS.Fill(mass);
                    hMassOSWeighted.Fill(mass, eventWeight);
                    ++nOS;
                } else {
                    hMassSS.Fill(mass);
                    hMassSSWeighted.Fill(mass, eventWeight);
                    ++nSS;
                }
            }
        }
    }

    // --------------------------------------------------
    // Save output
    // --------------------------------------------------
    output->cd();

    hMassOS.Write();
    hMassSS.Write();
    hMassOSWeighted.Write();
    hMassSSWeighted.Write();

    output->Close();
    input->Close();

    std::cout << "\nReconstruction complete\n"
              << "Input: " << inputName << '\n'
              << "Output: " << outputName << '\n'
              << "Events processed: " << eventsProcessed << '\n'
              << "Events without a matching weight: "
              << eventsWithoutWeight << '\n'
              << "Opposite-sign candidates: " << nOS << '\n'
              << "Same-sign candidates: " << nSS << '\n';

    return 0;
}
