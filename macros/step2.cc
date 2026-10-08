#include "analysisUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// Step 2: K pi invariant-mass spectra over 0 - 3 GeV/c^2.
//   - K pi candidates (opposite signs)
//   - K pi candidates (same signs)
//   - K pi candidates (OS - SS), with the truth-matched D0 as a line
//
// Usage: ./step2 [inputFile] [--tag=name]
//                [--min-daughter-pt=X] [--min-pair-pt=X] [--max-eta=X]
//   inputFile  default data/d0_charm.root
//
// The cuts default to NONE so the whole spectrum (K*(892), D0, ...)
// is visible. Tighten them from the command line to narrow things down,
// e.g.  ./step2 --min-daughter-pt=1 --min-pair-pt=3 --max-eta=1

int main(int argc, char* argv[]) {

    const std::string inputName =
        getInputFile(argc, argv, "data/d0_charm.root");
    const std::string tag = getTag(argc, argv, inputName);

    // -----------------------------
    // Cuts (default: none)
    // -----------------------------
    const double minDaughterPt = getNumber(argc, argv, "min-daughter-pt", 0.0);
    const double minPairPt     = getNumber(argc, argv, "min-pair-pt", 0.0);
    const double maxAbsEta     = getNumber(argc, argv, "max-eta", 100.0);

    // -----------------------------
    // Mass axis
    // -----------------------------
    const int    nBins   = 600;   // 5 MeV/c^2 per bin
    const double massMin = 0.0;
    const double massMax = 3.0;

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

    // -----------------------------
    // Event loop
    // -----------------------------
    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {

        tree->GetEntry(iEvent);

        for (size_t i = 0; i < pid->size(); ++i) {

            // K candidates
            if (std::abs(pid->at(i)) != 321)
                continue;

            const double ptK = std::hypot(px->at(i), py->at(i));
            const double etaK =
                (ptK > 1e-6) ? std::asinh(pz->at(i) / ptK) : 0.0;

            if (ptK < minDaughterPt || std::abs(etaK) > maxAbsEta)
                continue;

            for (size_t j = 0; j < pid->size(); ++j) {

                // pi candidates
                if (std::abs(pid->at(j)) != 211)
                    continue;

                // Don't pair a particle with itself
                if (i == j)
                    continue;

                const double ptPi = std::hypot(px->at(j), py->at(j));
                const double etaPi =
                    (ptPi > 1e-6) ? std::asinh(pz->at(j) / ptPi) : 0.0;

                if (ptPi < minDaughterPt || std::abs(etaPi) > maxAbsEta)
                    continue;

                const bool oppositeSign = charge->at(i) * charge->at(j) < 0;
                const bool sameSign     = charge->at(i) * charge->at(j) > 0;

                // Invariant mass
                const double totalPx = px->at(i) + px->at(j);
                const double totalPy = py->at(i) + py->at(j);
                const double totalPz = pz->at(i) + pz->at(j);
                const double totalE  = energy->at(i) + energy->at(j);

                if (std::hypot(totalPx, totalPy) < minPairPt)
                    continue;

                const double massSquared =
                    totalE * totalE
                    - totalPx * totalPx
                    - totalPy * totalPy
                    - totalPz * totalPz;

                if (massSquared <= 0)
                    continue;

                const double mass = std::sqrt(massSquared);

                if (oppositeSign)
                    hMassOS.Fill(mass);

                if (sameSign)
                    hMassSS.Fill(mass);

                // Truth matching: both come from the same D0 -> K pi
                if (d0Parent->at(i) >= 0 &&
                    d0Parent->at(i) == d0Parent->at(j)) {

                    hMassTrue.Fill(mass);
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
    std::cout << "\n===== Step 2 summary =====\n"
              << "Input:  " << inputName << "\n"
              << "Cuts:   daughter pT > " << minDaughterPt
              << ", pair pT > " << minPairPt
              << ", |eta| < " << maxAbsEta << "\n"
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
    // Save histograms
    // -----------------------------
    TFile output(("data/step2" + tag + ".root").c_str(), "RECREATE");
    hMassOS.Write();
    hMassSS.Write();
    hMassSub->Write();
    hMassTrue.Write();
    output.Close();
    input.Close();

    // -----------------------------
    // Plots
    // -----------------------------
    gStyle->SetOptStat(0);

    hMassOS.SetLineColor(kBlue + 1);
    hMassSS.SetLineColor(kRed + 1);
    hMassSub->SetLineColor(kBlack);

    // The truth D0 is tiny next to the OS - SS spectrum, so a scaled copy
    // is drawn as the line; the scale factor is shown in the legend.
    const double maxSub  = hMassSub->GetMaximum();
    const double maxTrue = hMassTrue.GetMaximum();
    const double scale =
        (maxTrue > 0 && maxSub > 0) ? 0.5 * maxSub / maxTrue : 1.0;

    TH1F* hTrueDraw = static_cast<TH1F*>(hMassTrue.Clone("hTrueDraw"));
    hTrueDraw->SetDirectory(nullptr);
    hTrueDraw->Scale(scale);
    hTrueDraw->SetLineColor(kGreen + 2);
    hTrueDraw->SetLineWidth(2);

    hMassSub->SetMaximum(1.2 * maxSub);

    auto drawSub = [&]() {
        hMassSub->Draw("HIST");
        hTrueDraw->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.50, 0.75, 0.88, 0.88);
        leg->AddEntry(hMassSub, "OS - SS", "l");
        leg->AddEntry(hTrueDraw,
                      Form("D^{0} truth (#times%.2g)", scale), "l");
        leg->Draw();
    };

    // Combined canvas
    {
        TCanvas c("c", "Step 2", 1800, 600);
        c.Divide(3, 1);

        c.cd(1);
        hMassOS.Draw("HIST");

        c.cd(2);
        hMassSS.Draw("HIST");

        c.cd(3);
        drawSub();

        c.SaveAs(("plots/step2" + tag + ".png").c_str());
    }

    // Individual PNGs
    {
        TCanvas c1("c1", "OS", 800, 600);
        hMassOS.Draw("HIST");
        c1.SaveAs(("plots/kpi_os" + tag + ".png").c_str());

        TCanvas c2("c2", "SS", 800, 600);
        hMassSS.Draw("HIST");
        c2.SaveAs(("plots/kpi_ss" + tag + ".png").c_str());

        TCanvas c3("c3", "OS - SS", 800, 600);
        drawSub();
        c3.SaveAs(("plots/kpi_os_minus_ss" + tag + ".png").c_str());
    }

    delete hTrueDraw;
    delete hMassSub;

    return 0;
}
