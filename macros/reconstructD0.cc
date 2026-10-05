#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include <vector>
#include <cmath>
#include <algorithm>

int main() {

    // -----------------------------
    // Analysis cuts (tune these!)
    // -----------------------------
    const double minDaughterPt = 1.0;   // GeV/c, each of K and pi
    const double minPairPt     = 3.0;   // GeV/c, K pi pair
    const double maxAbsEta     = 1.0;   // |eta| of each daughter

    // Make sure the output folders exist
    gSystem->mkdir("data", true);
    gSystem->mkdir("plots", true);

    // -----------------------------
    // Open input
    // -----------------------------
    TFile input("data/d0.root", "READ");

    TTree* tree = nullptr;
    input.GetObject("Particles", tree);

    if (!tree) {
        return 1;
    }

    // -----------------------------
    // Variables from the TTree
    // -----------------------------
    int event;

    std::vector<int>* pid = nullptr;
    std::vector<int>* charge = nullptr;

    std::vector<float>* px = nullptr;
    std::vector<float>* py = nullptr;
    std::vector<float>* pz = nullptr;
    std::vector<float>* energy = nullptr;

    std::vector<int>* d0Parent = nullptr;

    tree->SetBranchAddress("event", &event);

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

    // All opposite-sign K pi combinations
    TH1F hMassOS(
        "hMassOS",
        "Opposite-sign K#pi candidates;M_{K#pi} [GeV/c^{2}];Candidates",
        200,
        1.6,
        2.1
    );

    // Same-sign K pi combinations (background estimate)
    TH1F hMassSS(
        "hMassSS",
        "Same-sign K#pi candidates;M_{K#pi} [GeV/c^{2}];Candidates",
        200,
        1.6,
        2.1
    );

    // Truth-matched D0 candidates
    TH1F hMassTrue(
        "hMassTrue",
        "Truth-matched D^{0} candidates;M_{K#pi} [GeV/c^{2}];Candidates",
        200,
        1.6,
        2.1
    );

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

        // -----------------------------------
        // Loop over all particles
        // -----------------------------------

        for (size_t i = 0; i < pid->size(); ++i) {

            // K candidates
            if (std::abs(pid->at(i)) != 321)
                continue;

            const double ptK = std::hypot(px->at(i), py->at(i));
            if (ptK < minDaughterPt)
                continue;

            const double etaK = std::asinh(pz->at(i) / ptK);
            if (std::abs(etaK) > maxAbsEta)
                continue;

            for (size_t j = 0; j < pid->size(); ++j) {

                // pi candidates
                if (std::abs(pid->at(j)) != 211)
                    continue;

                // Don't pair a particle with itself
                if (i == j)
                    continue;

                const double ptPi = std::hypot(px->at(j), py->at(j));
                if (ptPi < minDaughterPt)
                    continue;

                const double etaPi = std::asinh(pz->at(j) / ptPi);
                if (std::abs(etaPi) > maxAbsEta)
                    continue;

                // -----------------------------------
                // Determine whether the pair is
                // opposite-sign or same-sign
                // -----------------------------------

                bool oppositeSign =
                    charge->at(i) * charge->at(j) < 0;

                bool sameSign =
                    charge->at(i) * charge->at(j) > 0;

                // -----------------------------------
                // Calculate invariant mass
                // -----------------------------------

                double totalPx =
                    px->at(i) + px->at(j);

                double totalPy =
                    py->at(i) + py->at(j);

                double totalPz =
                    pz->at(i) + pz->at(j);

                double totalE =
                    energy->at(i) + energy->at(j);

                // Pair pT cut
                if (std::hypot(totalPx, totalPy) < minPairPt)
                    continue;

                double massSquared =
                    totalE * totalE
                    - totalPx * totalPx
                    - totalPy * totalPy
                    - totalPz * totalPz;

                if (massSquared <= 0)
                    continue;

                double mass = std::sqrt(massSquared);

                // -----------------------------------
                // Fill OS / SS distributions
                // -----------------------------------

                if (oppositeSign)
                    hMassOS.Fill(mass);

                if (sameSign)
                    hMassSS.Fill(mass);

                // -----------------------------------
                // Truth matching (same cuts applied)
                // -----------------------------------

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
    hMassSub->SetTitle("OS - SS;M_{K#pi} [GeV/c^{2}];Candidates");
    hMassSub->Add(&hMassSS, -1.0);

    // -----------------------------
    // Save histograms
    // -----------------------------
    TFile output("data/d0_reconstructed.root", "RECREATE");

    hMassOS.Write();
    hMassSS.Write();
    hMassTrue.Write();
    hMassSub->Write();

    output.Close();
    input.Close();

    // -----------------------------
    // Plots
    // -----------------------------
    gStyle->SetOptStat(0);

    // Plot 1: OS only (same filename as before)
    {
        TCanvas canvas("canvas", "D0 invariant mass", 800, 600);
        hMassOS.SetLineColor(kBlue);
        hMassOS.Draw();
        canvas.SaveAs("plots/d0_mass_os.png");
    }

    // Plot 2: summary (OS vs SS, and OS-SS vs truth)
    {
        TCanvas canvas2("canvas2", "D0 summary", 1400, 600);
        canvas2.Divide(2, 1);

        hMassSS.SetLineColor(kRed);
        hMassTrue.SetLineColor(kGreen + 2);
        hMassSub->SetLineColor(kBlack);

        // Left: OS and SS overlaid
        canvas2.cd(1);
        hMassOS.SetMaximum(1.3 * std::max(hMassOS.GetMaximum(),
                                          hMassSS.GetMaximum()));
        hMassOS.Draw("HIST");
        hMassSS.Draw("HIST SAME");

        TLegend leg1(0.6, 0.75, 0.88, 0.88);
        leg1.AddEntry(&hMassOS, "Opposite-sign", "l");
        leg1.AddEntry(&hMassSS, "Same-sign", "l");
        leg1.Draw();

        // Right: background-subtracted vs truth-matched
        canvas2.cd(2);
        hMassSub->SetMaximum(1.3 * std::max(hMassSub->GetMaximum(),
                                            hMassTrue.GetMaximum()));
        hMassSub->Draw("HIST");
        hMassTrue.Draw("HIST SAME");

        TLegend leg2(0.6, 0.75, 0.88, 0.88);
        leg2.AddEntry(hMassSub, "OS - SS", "l");
        leg2.AddEntry(&hMassTrue, "Truth-matched D^{0}", "l");
        leg2.Draw();

        canvas2.SaveAs("plots/d0_mass_summary.png");
    }

    delete hMassSub;

    return 0;
}
