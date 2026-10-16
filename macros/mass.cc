#include "analysisUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// mass: K pi invariant-mass spectra over 0 - 3 GeV/c^2.
//   - K pi candidates (opposite signs)
//   - K pi candidates (same signs)
//   - K pi candidates (OS - SS), with the truth-matched D0 as a line
// Dashed vertical lines mark K*(892), K*(1430) and the D0.
//
// Usage: ./mass [inputFile] [--option=value ...]
//   inputFile  default data/d0_charm.root
//
// Options (defaults in brackets):
//   --tag=NAME            output tag            [from the input file name]
//   --min-daughter-pt=X   K and pi pT >= X      [0]     GeV/c
//   --max-daughter-pt=X   K and pi pT <= X      [none]  GeV/c
//   --min-pair-pt=X       pair pT >= X          [0]     GeV/c
//   --max-pair-pt=X       pair pT <= X          [none]  GeV/c
//   --max-eta=X           |eta| <= X            [none]
//   --smear-pt=X          relative pT smearing of K and pi   [0]
//   --smear-eta=X         absolute eta smearing              [0]
//   --smear-phi=X         absolute phi smearing (rad)        [0]
//   --seed=N              random seed for the smearing       [12345]
//   --pt-edges=a,b,c,...  also make OS - SS plots in bins of pair pT,
//                         e.g. --pt-edges=0,1,2,3,5,8
//
// The cuts default to NONE so the whole spectrum is visible. Tighten them
// to narrow things down, e.g.
//   ./mass --min-daughter-pt=0.6 --max-daughter-pt=2.2 --max-eta=1
//
// K and pi energies are always recomputed from (pT, eta, phi) with the PDG
// masses, so smearing a momentum changes the energy consistently.

namespace {

const double kMassK  = 0.493677;
const double kMassPi = 0.13957039;

struct Resonance {
    double      mass;    // GeV/c^2
    const char* label;
    int         color;
};

// Neutral K*(892), the K*(1430) region (K0* and K2* sit within ~10 MeV,
// both broad), and the D0.
const std::vector<Resonance> resonances = {
    {0.8955, "K*(892)",  kGray + 2},
    {1.43,   "K*(1430)", kGray + 2},
    {1.8648, "D^{0}",    kGreen + 2}
};

// Draw dashed vertical lines at the resonance masses on the current pad.
void drawResonanceLines() {

    gPad->Update();
    const double y1 = gPad->GetUymin();
    const double y2 = gPad->GetUymax();

    for (const auto& r : resonances) {

        TLine* line = new TLine(r.mass, y1, r.mass, y2);
        line->SetLineStyle(2);
        line->SetLineColor(r.color);
        line->Draw();

        TLatex* text = new TLatex(r.mass, y2 - 0.02 * (y2 - y1), r.label);
        text->SetTextAngle(90);
        text->SetTextAlign(33);
        text->SetTextSize(0.03);
        text->SetTextColor(r.color);
        text->Draw();
    }
}

std::vector<double> parseList(const std::string& s) {
    std::vector<double> v;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ','))
        if (!item.empty())
            v.push_back(std::atof(item.c_str()));
    return v;
}

// The truth D0 is tiny next to the OS - SS spectrum, so a scaled copy is
// drawn as the line (only scaled down if it would be taller than OS - SS).
double truthScale(double maxSub, double maxTrue) {
    return (maxTrue > 0 && maxSub > 0 && maxSub < maxTrue)
               ? maxSub / maxTrue : 1.0;
}

}  // namespace

int main(int argc, char* argv[]) {

    const std::string inputName =
        getInputFile(argc, argv, "data/d0_charm.root");
    const std::string tag = getTag(argc, argv, inputName);

    // -----------------------------
    // Cuts (default: none)
    // -----------------------------
    const double minDaughterPt = getNumber(argc, argv, "min-daughter-pt", 0.0);
    const double maxDaughterPt = getNumber(argc, argv, "max-daughter-pt", 1e9);
    const double minPairPt     = getNumber(argc, argv, "min-pair-pt", 0.0);
    const double maxPairPt     = getNumber(argc, argv, "max-pair-pt", 1e9);
    const double maxAbsEta     = getNumber(argc, argv, "max-eta", 100.0);

    // -----------------------------
    // Optional smearing of the K and pi momenta
    // -----------------------------
    const double smearPt  = getNumber(argc, argv, "smear-pt", 0.0);
    const double smearEta = getNumber(argc, argv, "smear-eta", 0.0);
    const double smearPhi = getNumber(argc, argv, "smear-phi", 0.0);
    const int    seed     = static_cast<int>(getNumber(argc, argv, "seed", 12345));
    TRandom3 rng(seed);

    // -----------------------------
    // Optional pair-pT bins
    // -----------------------------
    std::vector<double> edges = parseList(getOption(argc, argv, "pt-edges", ""));
    std::sort(edges.begin(), edges.end());
    const int nPtBins = edges.size() >= 2 ? static_cast<int>(edges.size()) - 1 : 0;

    // -----------------------------
    // Mass axis
    // -----------------------------
    const int    nBins   = 600;   // 5 MeV/c^2 per bin
    const double massMin = 0.0;
    const double massMax = 3.0;
    const int    nBinsPt = 150;   // 20 MeV/c^2 per bin, for the pair-pT plots

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
    std::vector<int>*   d0Parent = nullptr;

    tree->SetBranchAddress("pid", &pid);
    tree->SetBranchAddress("charge", &charge);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("d0Parent", &d0Parent);

    // -----------------------------
    // Histograms
    // -----------------------------
    TH1F hMassOS("hMassOS",
        "K#pi candidates (opposite signs);M_{K#pi} [GeV/c^{2}];Candidates",
        nBins, massMin, massMax);

    TH1F hMassSS("hMassSS",
        "K#pi candidates (same signs);M_{K#pi} [GeV/c^{2}];Candidates",
        nBins, massMin, massMax);

    TH1F hMassTrue("hMassTrue",
        "D^{0} (truth);M_{K#pi} [GeV/c^{2}];Candidates",
        nBins, massMin, massMax);

    // Detach from the input file so closing it doesn't delete them
    hMassOS.SetDirectory(nullptr);
    hMassSS.SetDirectory(nullptr);
    hMassTrue.SetDirectory(nullptr);

    // Per-pair-pT-bin histograms
    std::vector<TH1F*> binOS, binSS, binTrue;
    for (int b = 0; b < nPtBins; ++b) {
        const std::string n = std::to_string(b);
        binOS.push_back(new TH1F(("binOS" + n).c_str(), "", nBinsPt, massMin, massMax));
        binSS.push_back(new TH1F(("binSS" + n).c_str(), "", nBinsPt, massMin, massMax));
        binTrue.push_back(new TH1F(("binTrue" + n).c_str(), "", nBinsPt, massMin, massMax));
        binOS[b]->SetDirectory(nullptr);
        binSS[b]->SetDirectory(nullptr);
        binTrue[b]->SetDirectory(nullptr);
    }

    // -----------------------------
    // Event loop
    // -----------------------------
    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {

        tree->GetEntry(iEvent);

        const size_t n = pid->size();

        // Smear each K and pi ONCE per event, then pair them
        std::vector<double> tpx(n, 0.0), tpy(n, 0.0), tpz(n, 0.0), tE(n, 0.0);
        std::vector<double> tpt(n, 0.0), teta(n, 0.0);
        std::vector<char>   usable(n, 0);

        for (size_t k = 0; k < n; ++k) {

            const int absId = std::abs(pid->at(k));
            if (absId != 321 && absId != 211)
                continue;

            double pt = std::hypot(px->at(k), py->at(k));
            if (pt < 1e-6)
                continue;

            double eta = std::asinh(pz->at(k) / pt);
            double phi = std::atan2(py->at(k), px->at(k));

            if (smearPt  > 0) pt  *= (1.0 + rng.Gaus(0.0, smearPt));
            if (smearEta > 0) eta += rng.Gaus(0.0, smearEta);
            if (smearPhi > 0) phi += rng.Gaus(0.0, smearPhi);

            if (pt < 1e-3) pt = 1e-3;   // keep the track physical

            const double m = (absId == 321) ? kMassK : kMassPi;

            tpt[k]  = pt;
            teta[k] = eta;
            tpx[k]  = pt * std::cos(phi);
            tpy[k]  = pt * std::sin(phi);
            tpz[k]  = pt * std::sinh(eta);
            tE[k]   = std::sqrt(tpx[k] * tpx[k] + tpy[k] * tpy[k]
                                + tpz[k] * tpz[k] + m * m);
            usable[k] = 1;
        }

        for (size_t i = 0; i < n; ++i) {

            // K candidates
            if (std::abs(pid->at(i)) != 321 || !usable[i])
                continue;

            if (tpt[i] < minDaughterPt || tpt[i] > maxDaughterPt ||
                std::abs(teta[i]) > maxAbsEta)
                continue;

            for (size_t j = 0; j < n; ++j) {

                // pi candidates
                if (std::abs(pid->at(j)) != 211 || !usable[j])
                    continue;

                // Don't pair a particle with itself
                if (i == j)
                    continue;

                if (tpt[j] < minDaughterPt || tpt[j] > maxDaughterPt ||
                    std::abs(teta[j]) > maxAbsEta)
                    continue;

                const bool oppositeSign = charge->at(i) * charge->at(j) < 0;
                const bool sameSign     = charge->at(i) * charge->at(j) > 0;

                // Invariant mass
                const double totalPx = tpx[i] + tpx[j];
                const double totalPy = tpy[i] + tpy[j];
                const double totalPz = tpz[i] + tpz[j];
                const double totalE  = tE[i] + tE[j];

                const double pairPt = std::hypot(totalPx, totalPy);
                if (pairPt < minPairPt || pairPt > maxPairPt)
                    continue;

                const double massSquared =
                    totalE * totalE
                    - totalPx * totalPx
                    - totalPy * totalPy
                    - totalPz * totalPz;

                if (massSquared <= 0)
                    continue;

                const double mass = std::sqrt(massSquared);

                const bool isTruth =
                    d0Parent->at(i) >= 0 &&
                    d0Parent->at(i) == d0Parent->at(j);

                if (oppositeSign)
                    hMassOS.Fill(mass);

                if (sameSign)
                    hMassSS.Fill(mass);

                // Truth matching: both come from the same D0 -> K pi
                if (isTruth)
                    hMassTrue.Fill(mass);

                // Pair-pT bin
                for (int b = 0; b < nPtBins; ++b) {
                    if (pairPt >= edges[b] && pairPt < edges[b + 1]) {
                        if (oppositeSign) binOS[b]->Fill(mass);
                        if (sameSign)     binSS[b]->Fill(mass);
                        if (isTruth)      binTrue[b]->Fill(mass);
                        break;
                    }
                }
            }
        }
    }

    // -----------------------------
    // Background subtraction: OS - SS
    // -----------------------------
    TH1F* hMassSub = static_cast<TH1F*>(hMassOS.Clone("hMassSub"));
    hMassSub->SetDirectory(nullptr);
    hMassSub->SetTitle(
        "K#pi candidates (OS - SS);M_{K#pi} [GeV/c^{2}];Candidates");
    hMassSub->Add(&hMassSS, -1.0);

    // -----------------------------
    // Console summary
    // -----------------------------
    std::cout << "\n===== mass summary =====\n"
              << "Input:  " << inputName << "\n"
              << "Cuts:   daughter pT in [" << minDaughterPt << ", "
              << (maxDaughterPt < 1e8 ? std::to_string(maxDaughterPt) : "inf")
              << "], pair pT in [" << minPairPt << ", "
              << (maxPairPt < 1e8 ? std::to_string(maxPairPt) : "inf")
              << "], |eta| < " << maxAbsEta << "\n"
              << "Smear:  pT " << smearPt << ", eta " << smearEta
              << ", phi " << smearPhi << "\n"
              << "OS candidates:      " << hMassOS.GetEntries() << "\n"
              << "SS candidates:      " << hMassSS.GetEntries() << "\n"
              << "Truth-matched D0:   " << hMassTrue.GetEntries() << "\n";

    {
        // Strongest OS - SS peak between 0.8 and 1.0 (K*(892) region)
        TAxis* ax = hMassSub->GetXaxis();
        ax->SetRangeUser(0.80, 1.00);
        std::cout << "OS - SS maximum in 0.80-1.00 GeV/c^2 at "
                  << ax->GetBinCenter(hMassSub->GetMaximumBin())
                  << " GeV/c^2 (K*(892) expected near 0.89-0.90)\n";
        ax->SetRange(0, 0);   // restore the full range
    }

    std::cout << "==========================\n";

    // -----------------------------
    // Pair-pT bins: OS - SS
    // -----------------------------
    std::vector<TH1F*> binSub;
    for (int b = 0; b < nPtBins; ++b) {
        TH1F* h = static_cast<TH1F*>(binOS[b]->Clone(("binSub" + std::to_string(b)).c_str()));
        h->SetDirectory(nullptr);
        h->Add(binSS[b], -1.0);
        h->SetTitle(Form("OS - SS, pair p_{T} #in [%.3g, %.3g) GeV/c;"
                         "M_{K#pi} [GeV/c^{2}];Candidates",
                         edges[b], edges[b + 1]));
        binSub.push_back(h);
    }

    // -----------------------------
    // Save histograms
    // -----------------------------
    TFile output(("data/mass" + tag + ".root").c_str(), "RECREATE");
    hMassOS.Write();
    hMassSS.Write();
    hMassSub->Write();
    hMassTrue.Write();
    for (int b = 0; b < nPtBins; ++b) {
        binSub[b]->Write();
        binTrue[b]->SetName(("binTrue" + std::to_string(b)).c_str());
        binTrue[b]->Write();
    }
    output.Close();
    input.Close();

    // -----------------------------
    // Plots
    // -----------------------------
    gStyle->SetOptStat(0);

    hMassOS.SetLineColor(kBlue + 1);
    hMassSS.SetLineColor(kRed + 1);
    hMassSub->SetLineColor(kBlack);

    const double maxSub  = hMassSub->GetMaximum();
    const double maxTrue = hMassTrue.GetMaximum();
    const double scale   = truthScale(maxSub, maxTrue);

    TH1F* hTrueDraw = static_cast<TH1F*>(hMassTrue.Clone("hTrueDraw"));
    hTrueDraw->SetDirectory(nullptr);
    hTrueDraw->Scale(scale);
    hTrueDraw->SetLineColor(kGreen + 2);
    hTrueDraw->SetLineWidth(2);

    hMassSub->SetMaximum(1.2 * maxSub);

    auto drawSub = [&]() {
        hMassSub->Draw("HIST");
        hTrueDraw->Draw("HIST SAME");
        drawResonanceLines();

        TLegend* leg = new TLegend(0.50, 0.75, 0.88, 0.88);
        leg->AddEntry(hMassSub, "OS - SS", "l");
        leg->AddEntry(hTrueDraw,
                      Form("D^{0} truth (#times%.2g)", scale), "l");
        leg->Draw();
    };

    // Combined canvas
    {
        TCanvas c("c", "mass", 1800, 600);
        c.Divide(3, 1);

        c.cd(1);
        hMassOS.Draw("HIST");
        drawResonanceLines();

        c.cd(2);
        hMassSS.Draw("HIST");
        drawResonanceLines();

        c.cd(3);
        drawSub();

        c.SaveAs(("plots/mass" + tag + ".png").c_str());
    }

    // Individual PNGs
    {
        TCanvas c1("c1", "OS", 800, 600);
        hMassOS.Draw("HIST");
        drawResonanceLines();
        c1.SaveAs(("plots/kpi_os" + tag + ".png").c_str());

        TCanvas c2("c2", "SS", 800, 600);
        hMassSS.Draw("HIST");
        drawResonanceLines();
        c2.SaveAs(("plots/kpi_ss" + tag + ".png").c_str());

        TCanvas c3("c3", "OS - SS", 800, 600);
        drawSub();
        c3.SaveAs(("plots/kpi_os_minus_ss" + tag + ".png").c_str());
    }

    // Pair-pT-binned OS - SS
    if (nPtBins > 0) {

        const int nCols = std::min(nPtBins, 3);
        const int nRows = (nPtBins + nCols - 1) / nCols;

        TCanvas cb("cb", "mass in pair pT bins", 600 * nCols, 450 * nRows);
        cb.Divide(nCols, nRows);

        std::vector<TH1F*> truthDraw;

        for (int b = 0; b < nPtBins; ++b) {

            cb.cd(b + 1);

            binSub[b]->SetLineColor(kBlack);

            const double mSub  = binSub[b]->GetMaximum();
            const double mTrue = binTrue[b]->GetMaximum();
            const double s     = truthScale(mSub, mTrue);

            TH1F* t = static_cast<TH1F*>(binTrue[b]->Clone(
                ("binTrueDraw" + std::to_string(b)).c_str()));
            t->SetDirectory(nullptr);
            t->Scale(s);
            t->SetLineColor(kGreen + 2);
            t->SetLineWidth(2);
            truthDraw.push_back(t);

            binSub[b]->SetMaximum(1.2 * std::max(mSub, 1.0));
            binSub[b]->Draw("HIST");
            t->Draw("HIST SAME");
            drawResonanceLines();

            TLegend* leg = new TLegend(0.45, 0.78, 0.88, 0.88);
            leg->AddEntry(t, Form("D^{0} truth (#times%.2g)", s), "l");
            leg->Draw();
        }

        cb.SaveAs(("plots/mass_ptbins" + tag + ".png").c_str());

        for (TH1F* t : truthDraw) delete t;
    }

    delete hTrueDraw;
    delete hMassSub;
    for (TH1F* h : binOS)   delete h;
    for (TH1F* h : binSS)   delete h;
    for (TH1F* h : binTrue) delete h;
    for (TH1F* h : binSub)  delete h;

    return 0;
}
