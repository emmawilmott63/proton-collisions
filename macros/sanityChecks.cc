#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
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

// Usage:
//   ./sanityChecks [inputFile] [tag]
//
//   inputFile  default data/d0_charm.root
//   tag        optional label appended to output names; default is taken
//              from the file name (data/d0_inclusive.root -> "inclusive")
//              -> plots/sanity_checks_inclusive.png

int main(int argc, char* argv[]) {

    const std::string inputName = (argc > 1) ? argv[1] : "data/d0_charm.root";

    std::string tag;
    if (argc > 2) {
        tag = std::string("_") + argv[2];
    } else {
        std::string stem = inputName.substr(inputName.find_last_of('/') + 1);
        if (stem.size() > 5 && stem.substr(stem.size() - 5) == ".root")
            stem.resize(stem.size() - 5);
        if (stem.rfind("d0_", 0) == 0)
            stem = stem.substr(3);
        tag = "_" + stem;
    }

    // -----------------------------
    // Settings
    // -----------------------------
    const double jetR         = 0.4;   // anti-kt radius
    const double jetPtMin     = 2.0;   // GeV/c, minimum jet pT kept
    const double jetMaxAbsEta = 0.6;   // |eta| acceptance for jets
    const int    nPidBins     = 12;    // number of species shown in the PID plot

    // Reconstruction cuts: KEEP IN SYNC with reconstructD0.cc
    const double minDaughterPt = 1.0;
    const double minPairPt     = 3.0;
    const double maxAbsEta     = 1.0;

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

    // Truth D0 branches exist only if the file was made by the newer
    // generateD0; handle older files gracefully.
    std::vector<float>* d0TruthPt    = nullptr;
    std::vector<float>* d0TruthY     = nullptr;
    std::vector<int>*   d0TruthToKPi = nullptr;

    const bool haveTruth = tree->GetBranch("d0TruthPt") != nullptr;

    if (haveTruth) {
        tree->SetBranchAddress("d0TruthPt", &d0TruthPt);
        tree->SetBranchAddress("d0TruthY", &d0TruthY);
        tree->SetBranchAddress("d0TruthToKPi", &d0TruthToKPi);
    } else {
        std::cout << "NOTE: no d0Truth* branches found; regenerate with the "
                     "updated generateD0 to get D0 spectra.\n";
    }

    // -----------------------------
    // Histograms (detached from the input file)
    // -----------------------------
    auto make = [](const char* name, const char* title,
                   int nBins, double lo, double hi) {
        TH1F* h = new TH1F(name, title, nBins, lo, hi);
        h->SetDirectory(nullptr);
        return h;
    };

    // Jets
    TH1F* hJetPt = make("hJetPt",
        "All jets;Jet p_{T} [GeV/c];Jets", 60, 0, 30);
    TH1F* hLeadJetPt = make("hLeadJetPt",
        "Leading jet;Leading jet p_{T} [GeV/c];Events", 60, 0, 30);
    TH1F* hNJets = make("hNJets",
        "Jets per event;N_{jets};Events", 10, -0.5, 9.5);
    TH1F* hJetNConst = make("hJetNConst",
        "Jet constituents;Constituents per jet;Jets", 40, -0.5, 39.5);

    // Event-level and particle-level
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

    // Identified-hadron pT spectra (normalised per event at the end)
    TH1F* hPtPi = make("hPtPi",
        ";p_{T} [GeV/c];1/N_{ev} dN/dp_{T} [(GeV/c)^{-1}]", 100, 0, 15);
    TH1F* hPtK = make("hPtK",
        ";p_{T} [GeV/c];1/N_{ev} dN/dp_{T} [(GeV/c)^{-1}]", 100, 0, 15);
    TH1F* hPtP = make("hPtP",
        ";p_{T} [GeV/c];1/N_{ev} dN/dp_{T} [(GeV/c)^{-1}]", 100, 0, 15);

    // D0 truth spectra
    TH1F* hD0Pt = make("hD0Pt",
        "D^{0} p_{T} spectrum;p_{T} [GeV/c];1/N_{ev} dN/dp_{T} [(GeV/c)^{-1}]",
        40, 0, 20);
    TH1F* hD0KPiPt = make("hD0KPiPt", ";p_{T} [GeV/c];", 40, 0, 20);
    TH1F* hD0AccPt = make("hD0AccPt", ";p_{T} [GeV/c];", 40, 0, 20);
    TH1F* hD0Y = make("hD0Y",
        "D^{0} rapidity;y;D^{0}", 60, -3, 3);

    std::map<int, long long> pidCounts;   // |PDG ID| -> count

    // Exact mean pT: 0 = pi, 1 = K, 2 = p, 3 = D0
    double    sumPt[4] = {0, 0, 0, 0};
    long long cntPt[4] = {0, 0, 0, 0};

    long long nEvents = 0;
    long long nParticlesTotal = 0;
    long long nJetsTotal = 0;
    long long nD0Total = 0;
    long long nD0AllTruth = 0;
    long long nD0KPiTruth = 0;
    long long nD0InAcc = 0;
    long long nBadTruthPairs = 0;

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

        // D0 parent index -> indices of its stored daughters
        std::map<int, std::vector<size_t>> d0Daughters;

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

            // Identified hadron spectra
            if (absId == 211) { hPtPi->Fill(pt); sumPt[0] += pt; ++cntPt[0]; }
            if (absId == 321) { hPtK->Fill(pt);  sumPt[1] += pt; ++cntPt[1]; }
            if (absId == 2212){ hPtP->Fill(pt);  sumPt[2] += pt; ++cntPt[2]; }

            if (d0Parent->at(i) >= 0)
                d0Daughters[d0Parent->at(i)].push_back(i);

            // Neutrinos are invisible: leave them out of the jets
            if (absId == 12 || absId == 14 || absId == 16)
                continue;

            jetInputs.emplace_back(px->at(i), py->at(i),
                                   pz->at(i), energy->at(i));
        }

        nParticlesTotal += pid->size();
        hNParticles->Fill(pid->size());
        hNCharged->Fill(nCharged);

        hD0PerEvent->Fill(d0Daughters.size());
        nD0Total += d0Daughters.size();

        // -----------------------------
        // Truth D0 spectra
        // -----------------------------
        if (haveTruth) {
            for (size_t k = 0; k < d0TruthPt->size(); ++k) {

                const double pt = d0TruthPt->at(k);

                hD0Pt->Fill(pt);
                hD0Y->Fill(d0TruthY->at(k));
                sumPt[3] += pt;
                ++cntPt[3];
                ++nD0AllTruth;

                if (d0TruthToKPi->at(k) == 1) {
                    hD0KPiPt->Fill(pt);
                    ++nD0KPiTruth;
                }
            }
        }

        // -----------------------------
        // D0 -> K pi decays whose daughters pass the
        // reconstruction cuts (i.e. could be found)
        // -----------------------------
        for (const auto& entry : d0Daughters) {

            const std::vector<size_t>& d = entry.second;

            if (d.size() != 2) {
                ++nBadTruthPairs;
                continue;
            }

            const size_t a = d[0];
            const size_t b = d[1];

            const double ptA = std::hypot(px->at(a), py->at(a));
            const double ptB = std::hypot(px->at(b), py->at(b));

            if (ptA < minDaughterPt || ptB < minDaughterPt)
                continue;

            if (std::abs(std::asinh(pz->at(a) / ptA)) > maxAbsEta ||
                std::abs(std::asinh(pz->at(b) / ptB)) > maxAbsEta)
                continue;

            const double pairPt = std::hypot(px->at(a) + px->at(b),
                                             py->at(a) + py->at(b));

            if (pairPt < minPairPt)
                continue;

            hD0AccPt->Fill(pairPt);
            ++nD0InAcc;
        }

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
    // Acceptance x efficiency vs pT (before per-event scaling)
    // -----------------------------
    TH1F* hD0Eff = static_cast<TH1F*>(hD0AccPt->Clone("hD0Eff"));
    hD0Eff->SetDirectory(nullptr);
    hD0Eff->SetTitle("Fraction of D^{0}#rightarrowK#pi passing reconstruction "
                     "cuts;D^{0} p_{T} [GeV/c];Acceptance #times efficiency");
    hD0Eff->Divide(hD0AccPt, hD0KPiPt, 1.0, 1.0, "B");

    // -----------------------------
    // Normalise spectra: per event, per GeV/c
    // -----------------------------
    const double nEv = std::max<double>(nEvents, 1.0);

    for (TH1F* h : {hPtPi, hPtK, hPtP, hD0Pt, hD0KPiPt, hD0AccPt})
        h->Scale(1.0 / (nEv * h->GetBinWidth(1)));

    hPtPi->SetTitle("#pi^{#pm}");
    hPtK->SetTitle("K^{#pm}");
    hPtP->SetTitle("p, #bar{p}");

    hD0KPiPt->SetTitle("D^{0}#rightarrowK#pi (truth)");
    hD0AccPt->SetTitle("D^{0}#rightarrowK#pi passing cuts");

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
    auto meanPt = [&](int k) {
        return cntPt[k] > 0 ? sumPt[k] / cntPt[k] : 0.0;
    };

    std::cout << "\n===== Sanity-check summary =====\n"
              << "Input:                      " << inputName << "\n"
              << "Events read:                " << nEvents << "\n"
              << "Final-state particles:      " << nParticlesTotal
              << "  (" << double(nParticlesTotal) / nEv << " per event)\n"
              << "Jets (pT>" << jetPtMin << ", |eta|<" << jetMaxAbsEta << "): "
              << nJetsTotal
              << "  (" << double(nJetsTotal) / nEv << " per event)\n"
              << "Truth D0 -> K pi (matched): " << nD0Total << "\n";

    if (haveTruth) {
        std::cout << "All truth D0 + D0bar:       " << nD0AllTruth
                  << "  (" << double(nD0AllTruth) / nEv << " per event)\n"
                  << "  of which direct K pi:     " << nD0KPiTruth
                  << "  (fraction "
                  << double(nD0KPiTruth) / std::max(nD0AllTruth, 1LL) << ")\n"
                  << "  passing reco cuts:        " << nD0InAcc
                  << "  (fraction of K pi "
                  << double(nD0InAcc) / std::max(nD0KPiTruth, 1LL) << ")\n";
    }

    if (nBadTruthPairs > 0)
        std::cout << "WARNING: " << nBadTruthPairs
                  << " D0 parents did not have exactly 2 stored daughters.\n";

    std::cout << "Mean pT [GeV/c]:  pi " << meanPt(0)
              << "   K " << meanPt(1)
              << "   p " << meanPt(2);
    if (haveTruth)
        std::cout << "   D0 " << meanPt(3);
    std::cout << "\nTop species (|PDG ID|: count):\n";

    for (int b = 0; b < nShown; ++b)
        std::cout << "  " << sorted[b].first << ": " << sorted[b].second << "\n";

    std::cout << "================================\n";

    // -----------------------------
    // Save histograms
    // -----------------------------
    TFile output(("data/sanity_checks" + tag + ".root").c_str(), "RECREATE");

    std::vector<TH1F*> all = {
        hJetPt, hLeadJetPt, hNJets, hJetNConst, hPid,
        hNParticles, hNCharged, hPartPt, hPartEta, hD0PerEvent,
        hPtPi, hPtK, hPtP, hD0Pt, hD0KPiPt, hD0AccPt, hD0Eff, hD0Y
    };

    for (TH1F* h : all)
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

    // Summary canvas
    {
        TCanvas summary("summary", "Sanity checks", 1800, 1500);
        summary.Divide(3, 3);

        for (size_t k = 0; k < items.size(); ++k) {
            summary.cd(k + 1);
            gPad->SetLogy(items[k].logY);
            items[k].h->SetLineColor(items[k].color);
            items[k].h->SetFillColorAlpha(items[k].color, 0.25);
            items[k].h->Draw("HIST");
        }

        summary.SaveAs(("plots/sanity_checks" + tag + ".png").c_str());
    }

    // Individual PNGs
    {
        TCanvas c1("c1", "Jet pT", 800, 600);
        c1.SetLogy();
        hJetPt->Draw("HIST");
        c1.SaveAs(("plots/jet_pt" + tag + ".png").c_str());

        TCanvas c2("c2", "Particle IDs", 900, 600);
        c2.SetLogy();
        c2.SetBottomMargin(0.12);
        hPid->Draw("HIST");
        c2.SaveAs(("plots/particle_ids" + tag + ".png").c_str());

        TCanvas c3("c3", "D0 per event", 800, 600);
        c3.SetLogy();
        hD0PerEvent->Draw("HIST");
        c3.SaveAs(("plots/d0_per_event" + tag + ".png").c_str());
    }

    // pT spectra canvas
    {
        gStyle->SetOptStat(0);

        TCanvas spectra("spectra", "pT spectra", 1400, 1100);
        spectra.Divide(2, 2);

        // (1) Identified hadron spectra
        spectra.cd(1);
        gPad->SetLogy();
        hPtPi->SetLineColor(kBlue + 1);
        hPtK->SetLineColor(kRed + 1);
        hPtP->SetLineColor(kGreen + 2);
        hPtPi->SetMinimum(1e-6);
        hPtPi->SetTitle("Charged hadron p_{T} spectra");
        hPtPi->Draw("HIST");
        hPtK->Draw("HIST SAME");
        hPtP->Draw("HIST SAME");

        TLegend leg1(0.65, 0.65, 0.88, 0.88);
        leg1.AddEntry(hPtPi, "#pi^{#pm}", "l");
        leg1.AddEntry(hPtK, "K^{#pm}", "l");
        leg1.AddEntry(hPtP, "p, #bar{p}", "l");
        leg1.Draw();

        // (2) D0 spectra
        spectra.cd(2);
        gPad->SetLogy();
        hD0Pt->SetLineColor(kBlack);
        hD0KPiPt->SetLineColor(kMagenta + 1);
        hD0AccPt->SetLineColor(kOrange + 1);
        hD0Pt->SetMinimum(1e-7);
        hD0Pt->SetTitle("D^{0} p_{T} spectra (truth)");
        hD0Pt->Draw("HIST");
        hD0KPiPt->Draw("HIST SAME");
        hD0AccPt->Draw("HIST SAME");

        TLegend leg2(0.45, 0.65, 0.88, 0.88);
        leg2.AddEntry(hD0Pt, "All D^{0}", "l");
        leg2.AddEntry(hD0KPiPt, "D^{0}#rightarrowK#pi", "l");
        leg2.AddEntry(hD0AccPt, "D^{0}#rightarrowK#pi passing cuts", "l");
        leg2.Draw();

        // (3) Acceptance x efficiency
        spectra.cd(3);
        hD0Eff->SetLineColor(kOrange + 1);
        hD0Eff->SetMinimum(0);
        hD0Eff->SetMaximum(1.1);
        hD0Eff->Draw("E");

        // (4) D0 rapidity
        spectra.cd(4);
        hD0Y->SetLineColor(kBlack);
        hD0Y->Draw("HIST");

        spectra.SaveAs(("plots/pt_spectra" + tag + ".png").c_str());
    }

    for (TH1F* h : all)
        delete h;

    return 0;
}
