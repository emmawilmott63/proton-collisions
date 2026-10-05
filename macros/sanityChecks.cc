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
#include <set>
#include <string>
#include <vector>

int main() {

    // -----------------------------
    // Settings
    // -----------------------------
    const double jetR         = 0.4;   // anti-kt radius
    const double jetPtMin     = 2.0;   // GeV/c, minimum jet pT kept
    const double jetMaxAbsEta = 0.6;   // |eta| acceptance for jets
    const int    nPidBins     = 12;    // number of species shown in the PID plot

    gSystem->mkdir("data", true);
    gSystem->mkdir("plots", true);

    // -----------------------------
    // Open input
    // -----------------------------
    TFile input("data/d0.root", "READ");

    TTree* tree = nullptr;
    input.GetObject("Particles", tree);

    if (!tree) {
        std::cerr << "Could not find TTree 'Particles' in data/d0.root\n";
        return 1;
    }

    std::vector<int>*   pid      = nullptr;
    std::vector<int>*   charge   = nullptr;
    std::vector<float>* px       = nullptr;
    std::vector<float>* py       = nullptr;
    std::vector<float>* pz       = nullptr;
    std::vector<float>* energy   = nullptr;
    std::vector<int>*   d0Parent = nullptr;

    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("energy", &energy);
    tree->SetBranchAddress("d0Parent", &d0Parent);

    // -----------------------------
    // Histograms (detached from the input file)
    // -----------------------------
    auto make = [](const char* name, const char* title,
                   int nBins, double lo, double hi) {
        TH1F* h = new TH1F(name, title, nBins, lo, hi);
        h->SetDirectory(nullptr);
        return h;
    };

    TH1F* hJetPt = make("hJetPt",
        "All jets;Jet p_{T} [GeV/c];Jets", 60, 0, 30);
    TH1F* hLeadJetPt = make("hLeadJetPt",
        "Leading jet;Leading jet p_{T} [GeV/c];Events", 60, 0, 30);
    TH1F* hNJets = make("hNJets",
        "Jets per event;N_{jets};Events", 10, -0.5, 9.5);
    TH1F* hJetNConst = make("hJetNConst",
        "Jet constituents;Constituents per jet;Jets", 40, -0.5, 39.5);

    TH1F* hNParticles = make("hNParticles",
        "Final-state particles per event;N_{particles};Events", 100, 0, 400);
    TH1F* hNCharged = make("hNCharged",
        "Charged particles per event;N_{charged};Events", 80, 0, 200);
    TH1F* hPartPt = make("hPartPt",
        "Particle p_{T};p_{T} [GeV/c];Particles", 100, 0, 15);
    TH1F* hPartEta = make("hPartEta",
        "Particle pseudorapidity;#eta;Particles", 120, -6, 6);
    TH1F* hD0PerEvent = make("hD0PerEvent",
        "D^{0}#rightarrowK#pi decays per event;N_{D^{0}#rightarrowK#pi};Events",
        6, -0.5, 5.5);

    std::map<int, long long> pidCounts;   // |PDG ID| -> count

    long long nEvents = 0;
    long long nParticlesTotal = 0;
    long long nJetsTotal = 0;
    long long nD0Total = 0;

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);

    // -----------------------------
    // Event loop
    // -----------------------------
    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {

        tree->GetEntry(iEvent);
        ++nEvents;

        std::vector<fastjet::PseudoJet> jetInputs;
        int nCharged = 0;
        std::set<int> d0Set;

        for (size_t i = 0; i < pid->size(); ++i) {

            const int    absId = std::abs(pid->at(i));
            const double pt    = std::hypot(px->at(i), py->at(i));

            pidCounts[absId]++;
            hPartPt->Fill(pt);

            // eta is undefined for pt = 0 (beam-axis particles)
            if (pt > 1e-6)
                hPartEta->Fill(std::asinh(pz->at(i) / pt));

            if (charge->at(i) != 0)
                ++nCharged;

            if (d0Parent->at(i) >= 0)
                d0Set.insert(d0Parent->at(i));

            // Neutrinos are invisible: leave them out of the jets
            if (absId == 12 || absId == 14 || absId == 16)
                continue;

            jetInputs.emplace_back(px->at(i), py->at(i),
                                   pz->at(i), energy->at(i));
        }

        nParticlesTotal += pid->size();
        hNParticles->Fill(pid->size());
        hNCharged->Fill(nCharged);

        hD0PerEvent->Fill(d0Set.size());
        nD0Total += d0Set.size();

        // -----------------------------
        // Jet clustering
        // -----------------------------
        fastjet::ClusterSequence cs(jetInputs, jetDef);
        std::vector<fastjet::PseudoJet> jets =
            fastjet::sorted_by_pt(cs.inclusive_jets(jetPtMin));

        int nJets = 0;
        bool filledLeading = false;

        for (const auto& jet : jets) {

            if (std::abs(jet.eta()) > jetMaxAbsEta)
                continue;

            ++nJets;
            hJetPt->Fill(jet.pt());
            hJetNConst->Fill(jet.constituents().size());

            // jets are pT-sorted, so the first one passing the cut leads
            if (!filledLeading) {
                hLeadJetPt->Fill(jet.pt());
                filledLeading = true;
            }
        }

        hNJets->Fill(nJets);
        nJetsTotal += nJets;
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

    TH1F* hPid = make("hPid",
        "Final-state particle species;|PDG ID|;Particles",
        std::max(nShown, 1), 0, std::max(nShown, 1));

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
    std::cout << "\n===== Sanity-check summary =====\n"
              << "Events read:                " << nEvents << "\n"
              << "Final-state particles:      " << nParticlesTotal
              << "  (" << double(nParticlesTotal) / std::max(nEvents, 1LL)
              << " per event)\n"
              << "Jets (pT>" << jetPtMin << ", |eta|<" << jetMaxAbsEta << "): "
              << nJetsTotal
              << "  (" << double(nJetsTotal) / std::max(nEvents, 1LL)
              << " per event)\n"
              << "Truth D0 -> K pi decays:    " << nD0Total << "\n"
              << "Top species (|PDG ID|: count):\n";

    for (int b = 0; b < nShown; ++b)
        std::cout << "  " << sorted[b].first << ": " << sorted[b].second << "\n";

    std::cout << "================================\n";

    // -----------------------------
    // Save histograms
    // -----------------------------
    TFile output("data/sanity_checks.root", "RECREATE");

    for (TH1F* h : {hJetPt, hLeadJetPt, hNJets, hJetNConst, hPid,
                    hNParticles, hNCharged, hPartPt, hPartEta, hD0PerEvent})
        h->Write();

    output.Close();
    input.Close();

    // -----------------------------
    // Plots
    // -----------------------------
    gStyle->SetOptStat(1110);   // name, entries, mean, std dev

    struct Item { TH1F* h; bool logY; int color; };
    const std::vector<Item> items = {
        {hJetPt,       true,  kBlue + 1},
        {hLeadJetPt,   false, kBlue + 1},
        {hNJets,       true,  kBlue + 1},
        {hJetNConst,   false, kBlue + 1},
        {hPid,         true,  kAzure + 1},
        {hNParticles,  false, kRed + 1},
        {hNCharged,    false, kRed + 1},
        {hPartPt,      true,  kGreen + 2},
        {hPartEta,     false, kGreen + 2}
    };

    // One summary canvas
    TCanvas summary("summary", "Sanity checks", 1800, 1500);
    summary.Divide(3, 3);

    for (size_t k = 0; k < items.size(); ++k) {
        summary.cd(k + 1);
        gPad->SetLogy(items[k].logY);
        items[k].h->SetLineColor(items[k].color);
        items[k].h->SetFillColorAlpha(items[k].color, 0.25);
        items[k].h->Draw("HIST");
    }

    summary.SaveAs("plots/sanity_checks.png");

    // Individual PNGs of the two plots you asked about, for convenience
    {
        TCanvas c1("c1", "Jet pT", 800, 600);
        c1.SetLogy();
        hJetPt->Draw("HIST");
        c1.SaveAs("plots/jet_pt.png");

        TCanvas c2("c2", "Particle IDs", 900, 600);
        c2.SetLogy();
        c2.SetBottomMargin(0.12);
        hPid->Draw("HIST");
        c2.SaveAs("plots/particle_ids.png");

        TCanvas c3("c3", "D0 per event", 800, 600);
        c3.SetLogy();
        hD0PerEvent->Draw("HIST");
        c3.SaveAs("plots/d0_per_event.png");
    }

    for (TH1F* h : {hJetPt, hLeadJetPt, hNJets, hJetNConst, hPid,
                    hNParticles, hNCharged, hPartPt, hPartEta, hD0PerEvent})
        delete h;

    return 0;
}
