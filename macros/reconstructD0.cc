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
        delete input;
        return 1;
    }

    // Run weight = generated cross section / successful generated events.
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
            std::cerr << "Warning: invalid run metadata for seed " << seed << '\n';
            continue;
        }
        eventWeights[seed] = sigmaGen / nGenerated;
        std::cout << "Seed " << seed << ": cross section = " << sigmaGen
                  << " mb, generated events = " << nGenerated
                  << ", event weight = " << eventWeights[seed] << " mb\n";
    }

    // Particle branches.
    int runSeed = 0;
    std::vector<int>* pid = nullptr;
    std::vector<int>* charge = nullptr;
    std::vector<float>* px = nullptr;
    std::vector<float>* py = nullptr;
    std::vector<float>* pz = nullptr;
    std::vector<float>* energy = nullptr;
    std::vector<int>* d0Parent = nullptr;

    particles->SetBranchAddress("runSeed", &runSeed);
    particles->SetBranchAddress("pid", &pid);
    particles->SetBranchAddress("charge", &charge);
    particles->SetBranchAddress("px", &px);
    particles->SetBranchAddress("py", &py);
    particles->SetBranchAddress("pz", &pz);
    particles->SetBranchAddress("energy", &energy);

    // This branch is written by generateCommon.h. It stores the PYTHIA event
    // record index of the D0 parent for daughters of a direct D0 -> K pi decay;
    // all other final-state particles have -1.
    if (!particles->GetBranch("d0Parent")) {
        std::cerr << "Error: input Particles tree has no d0Parent branch.\n"
                  << "Regenerate the input ROOT file with the current generator.\n";
        input->Close();
        delete input;
        return 1;
    }
    particles->SetBranchAddress("d0Parent", &d0Parent);

    TFile* output = TFile::Open(outputName.c_str(), "RECREATE");
    if (!output || output->IsZombie()) {
        std::cerr << "Error: could not create " << outputName << '\n';
        input->Close();
        delete input;
        return 1;
    }

    constexpr int nBins = 500;
    constexpr double massMin = 1.6;
    constexpr double massMax = 2.1;

    TH1D hMassOS("hMass_OS",
        "Opposite-sign K#pi candidates;M_{K#pi} [GeV/c^{2}];Candidates",
        nBins, massMin, massMax);
    TH1D hMassSS("hMass_SS",
        "Same-sign K#pi control;M_{K#pi} [GeV/c^{2}];Candidates",
        nBins, massMin, massMax);
    TH1D hMassOSWeighted("hMass_OS_weighted",
        "Weighted opposite-sign K#pi;M_{K#pi} [GeV/c^{2}];Weighted yield [mb]",
        nBins, massMin, massMax);
    TH1D hMassSSWeighted("hMass_SS_weighted",
        "Weighted same-sign K#pi;M_{K#pi} [GeV/c^{2}];Weighted yield [mb]",
        nBins, massMin, massMax);

    // Truth-matched subset: the K and pi must share a nonnegative d0Parent index.
    TH1D hMassTruth("hMass_truth",
        "Truth-matched D^{0}#rightarrow K#pi;M_{K#pi} [GeV/c^{2}];Candidates",
        nBins, massMin, massMax);
    TH1D hMassTruthWeighted("hMass_truth_weighted",
        "Weighted truth-matched D^{0}#rightarrow K#pi;M_{K#pi} [GeV/c^{2}];Weighted yield [mb]",
        nBins, massMin, massMax);

    // Opposite-sign candidates not truth-matched to the same generated D0.
    TH1D hMassOSComb("hMass_OS_combinatorial",
        "Opposite-sign non-truth-matched candidates;M_{K#pi} [GeV/c^{2}];Candidates",
        nBins, massMin, massMax);
    TH1D hMassOSCombWeighted("hMass_OS_combinatorial_weighted",
        "Weighted opposite-sign non-truth-matched candidates;M_{K#pi} [GeV/c^{2}];Weighted yield [mb]",
        nBins, massMin, massMax);

    for (TH1D* h : {&hMassOS, &hMassSS, &hMassOSWeighted, &hMassSSWeighted,
                    &hMassTruth, &hMassTruthWeighted, &hMassOSComb,
                    &hMassOSCombWeighted}) {
        h->Sumw2();
    }

    Long64_t eventsProcessed = 0;
    Long64_t eventsWithoutWeight = 0;
    Long64_t nOS = 0, nSS = 0, nTruth = 0, nOSCombinatorial = 0;
    Long64_t nTruthParentPairs = 0;

    for (Long64_t entry = 0; entry < particles->GetEntries(); ++entry) {
        particles->GetEntry(entry);

        if (!pid || !charge || !px || !py || !pz || !energy || !d0Parent) {
            std::cerr << "Error: missing particle data at entry " << entry << '\n';
            continue;
        }

        const auto weightIt = eventWeights.find(runSeed);
        if (weightIt == eventWeights.end()) {
            ++eventsWithoutWeight;
            continue;
        }
        const double eventWeight = weightIt->second;
        ++eventsProcessed;

        const size_t nParticles = pid->size();
        if (charge->size() != nParticles || px->size() != nParticles ||
            py->size() != nParticles || pz->size() != nParticles ||
            energy->size() != nParticles || d0Parent->size() != nParticles) {
            std::cerr << "Warning: mismatched particle-vector sizes at entry "
                      << entry << '\n';
            continue;
        }

        std::vector<size_t> kaons, pions;
        for (size_t i = 0; i < nParticles; ++i) {
            if (std::abs(pid->at(i)) == 321 && charge->at(i) != 0)
                kaons.push_back(i);
            if (std::abs(pid->at(i)) == 211 && charge->at(i) != 0)
                pions.push_back(i);
        }

        for (size_t k : kaons) {
            for (size_t pi : pions) {
                const double e = energy->at(k) + energy->at(pi);
                const double x = px->at(k) + px->at(pi);
                const double y = py->at(k) + py->at(pi);
                const double z = pz->at(k) + pz->at(pi);
                const double massSquared = e*e - x*x - y*y - z*z;

                if (massSquared < -1e-5)
                    continue;
                const double mass = std::sqrt(std::max(0.0, massSquared));

                const bool oppositeSign = charge->at(k) * charge->at(pi) < 0;
                // Truth match requires both particles to be tagged as daughters
                // of the exact same generated D0 in the PYTHIA event record.
                const int parentK = d0Parent->at(k);
                const int parentPi = d0Parent->at(pi);
                const bool truthMatched =
                    parentK >= 0 && parentPi >= 0 && parentK == parentPi;

                if (truthMatched)
                    ++nTruthParentPairs;

                if (oppositeSign) {
                    hMassOS.Fill(mass);
                    hMassOSWeighted.Fill(mass, eventWeight);
                    ++nOS;

                    if (truthMatched) {
                        hMassTruth.Fill(mass);
                        hMassTruthWeighted.Fill(mass, eventWeight);
                        ++nTruth;
                    } else {
                        hMassOSComb.Fill(mass);
                        hMassOSCombWeighted.Fill(mass, eventWeight);
                        ++nOSCombinatorial;
                    }
                } else {
                    hMassSS.Fill(mass);
                    hMassSSWeighted.Fill(mass, eventWeight);
                    ++nSS;
                }
            }
        }
    }

    output->cd();
    hMassOS.Write();
    hMassSS.Write();
    hMassOSWeighted.Write();
    hMassSSWeighted.Write();
    hMassTruth.Write();
    hMassTruthWeighted.Write();
    hMassOSComb.Write();
    hMassOSCombWeighted.Write();

    output->Close();
    input->Close();
    delete output;
    delete input;

    std::cout << "\nReconstruction complete\n"
              << "Input: " << inputName << '\n'
              << "Output: " << outputName << '\n'
              << "Events processed: " << eventsProcessed << '\n'
              << "Events without a matching weight: " << eventsWithoutWeight << '\n'
              << "Opposite-sign candidates: " << nOS << '\n'
              << "Same-sign candidates: " << nSS << '\n'
              << "Truth-matched OS candidates: " << nTruth << '\n'
              << "OS non-truth-matched candidates: " << nOSCombinatorial << '\n'
              << "All-sign pairs with the same nonnegative D0 parent index: "
              << nTruthParentPairs << '\n';

    return 0;
}
