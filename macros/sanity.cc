#include "analysisUtils.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// sanity: basic look at the generated data.
//   - jet pT
//   - leading-jet pT
//   - final-state particle IDs
//
// Usage: ./sanity [inputFile] [--tag=name]
//   inputFile  default data/d0_charm.root

int main(int argc, char* argv[]) {

    const std::string inputName =
        getInputFile(argc, argv, "data/d0_charm.root");
    const std::string tag = getTag(argc, argv, inputName);

    // -----------------------------
    // Settings
    // -----------------------------
    const double jetR         = 0.4;   // anti-kt radius
    const double jetPtMin     = 2.0;   // GeV/c, minimum jet pT kept
    const double jetMaxAbsEta = 0.6;   // |eta| acceptance for jets
    const int    nPidBins     = 12;    // species shown in the ID plot

    gSystem->mkdir("data", true);
    gSystem->mkdir("plots", true);

    // -----------------------------
    // Open input
    // -----------------------------
    TFile input(inputName.c_str(), "READ");

    TTree* tree = nullptr;
    input.GetObject("Particles", tree);

    if (!tree) {
        std::cerr << "Could not find TTree 'Particles' in "
                  << inputName << "\n";
        return 1;
    }

    std::vector<int>*   pid    = nullptr;
    std::vector<float>* px     = nullptr;
    std::vector<float>* py     = nullptr;
    std::vector<float>* pz     = nullptr;
    std::vector<float>* energy = nullptr;

    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("energy", &energy);

    // -----------------------------
    // Histograms (detached from the input file)
    // -----------------------------
    TH1F* hJetPt = new TH1F("hJetPt",
        "Jet p_{T};Jet p_{T} [GeV/c];Jets", 60, 0, 30);
    TH1F* hLeadJetPt = new TH1F("hLeadJetPt",
        "Leading jet p_{T};Leading jet p_{T} [GeV/c];Events", 60, 0, 30);

    hJetPt->SetDirectory(nullptr);
    hLeadJetPt->SetDirectory(nullptr);

    std::map<int, long long> pidCounts;   // |PDG ID| -> count

    long long nEvents = 0;
    long long nParticlesTotal = 0;
    long long nJetsTotal = 0;

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);

    // -----------------------------
    // Event loop
    // -----------------------------
    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {

        tree->GetEntry(iEvent);
        ++nEvents;

        std::vector<fastjet::PseudoJet> jetInputs;

        for (size_t i = 0; i < pid->size(); ++i) {

            const int absId = std::abs(pid->at(i));
            pidCounts[absId]++;

            // Neutrinos are invisible: leave them out of the jets
            if (absId == 12 || absId == 14 || absId == 16)
                continue;

            jetInputs.emplace_back(px->at(i), py->at(i),
                                   pz->at(i), energy->at(i));
        }

        nParticlesTotal += pid->size();

        fastjet::ClusterSequence cs(jetInputs, jetDef);
        std::vector<fastjet::PseudoJet> jets =
            fastjet::sorted_by_pt(cs.inclusive_jets(jetPtMin));

        bool filledLeading = false;

        for (const auto& jet : jets) {

            if (std::abs(jet.eta()) > jetMaxAbsEta)
                continue;

            ++nJetsTotal;
            hJetPt->Fill(jet.pt());

            // jets are pT-sorted, so the first one passing the cut leads
            if (!filledLeading) {
                hLeadJetPt->Fill(jet.pt());
                filledLeading = true;
            }
        }
    }

    // -----------------------------
    // Particle-ID histogram (most common species)
    // -----------------------------
    const std::map<int, std::string> names = {
        {22,   "#gamma"},
        {211,  "#pi^{#pm}"},
        {321,  "K^{#pm}"},
        {2212, "p/#bar{p}"},
        {2112, "n/#bar{n}"},
        {11,   "e^{#pm}"},
        {13,   "#mu^{#pm}"},
        {130,  "K^{0}_{L}"},
        {310,  "K^{0}_{S}"},
        {3122, "#Lambda"},
        {3222, "#Sigma^{+}"},
        {3112, "#Sigma^{-}"},
        {3312, "#Xi^{-}"},
        {3334, "#Omega^{-}"},
        {12,   "#nu_{e}"},
        {14,   "#nu_{#mu}"},
        {16,   "#nu_{#tau}"}
    };

    std::vector<std::pair<int, long long>> sorted(pidCounts.begin(),
                                                  pidCounts.end());
    std::sort(sorted.begin(), sorted.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    const int nShown = std::min<int>(nPidBins, sorted.size());

    TH1F* hPid = new TH1F("hPid",
        "Final-state particle IDs;Particle ID (|PDG ID|);Particles",
        std::max(nShown, 1), 0, std::max(nShown, 1));
    hPid->SetDirectory(nullptr);
    hPid->SetStats(0);

    for (int b = 0; b < nShown; ++b) {
        const int id = sorted[b].first;
        hPid->SetBinContent(b + 1, sorted[b].second);

        auto it = names.find(id);
        hPid->GetXaxis()->SetBinLabel(
            b + 1, it != names.end() ? it->second.c_str()
                                     : std::to_string(id).c_str());
    }

    // -----------------------------
    // Console summary
    // -----------------------------
    const double nEv = std::max<double>(nEvents, 1.0);

    std::cout << "\n===== sanity summary =====\n"
              << "Input:                 " << inputName << "\n"
              << "Events read:           " << nEvents << "\n"
              << "Final-state particles: " << nParticlesTotal
              << "  (" << nParticlesTotal / nEv << " per event)\n"
              << "Jets (pT>" << jetPtMin << ", |eta|<" << jetMaxAbsEta << "): "
              << nJetsTotal << "  (" << nJetsTotal / nEv << " per event)\n"
              << "Top species (|PDG ID|: count):\n";

    for (int b = 0; b < nShown; ++b)
        std::cout << "  " << sorted[b].first << ": " << sorted[b].second << "\n";

    std::cout << "==========================\n";

    // -----------------------------
    // Save histograms
    // -----------------------------
    TFile output(("data/sanity" + tag + ".root").c_str(), "RECREATE");
    hJetPt->Write();
    hLeadJetPt->Write();
    hPid->Write();
    output.Close();
    input.Close();

    // -----------------------------
    // Plots
    // -----------------------------
    gStyle->SetOptStat(1110);   // name, entries, mean, std dev

    hJetPt->SetLineColor(kBlue + 1);
    hJetPt->SetFillColorAlpha(kBlue + 1, 0.25);
    hLeadJetPt->SetLineColor(kBlue + 1);
    hLeadJetPt->SetFillColorAlpha(kBlue + 1, 0.25);
    hPid->SetLineColor(kAzure + 1);
    hPid->SetFillColorAlpha(kAzure + 1, 0.25);

    // Combined canvas
    {
        TCanvas c("c", "sanity", 1800, 600);
        c.Divide(3, 1);

        c.cd(1);
        gPad->SetLogy();
        hJetPt->Draw("HIST");

        c.cd(2);
        hLeadJetPt->Draw("HIST");

        c.cd(3);
        gPad->SetLogy();
        gPad->SetBottomMargin(0.12);
        hPid->Draw("HIST");

        c.SaveAs(("plots/sanity" + tag + ".png").c_str());
    }

    // Individual PNGs
    {
        TCanvas c1("c1", "Jet pT", 800, 600);
        c1.SetLogy();
        hJetPt->Draw("HIST");
        c1.SaveAs(("plots/jet_pt" + tag + ".png").c_str());

        TCanvas c2("c2", "Leading jet pT", 800, 600);
        hLeadJetPt->Draw("HIST");
        c2.SaveAs(("plots/leading_jet_pt" + tag + ".png").c_str());

        TCanvas c3("c3", "Particle IDs", 900, 600);
        c3.SetLogy();
        c3.SetBottomMargin(0.12);
        hPid->Draw("HIST");
        c3.SaveAs(("plots/particle_ids" + tag + ".png").c_str());
    }

    delete hJetPt;
    delete hLeadJetPt;
    delete hPid;

    return 0;
}
