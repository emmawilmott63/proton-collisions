#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include <vector>
#include <cmath>

int main() {

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

    // Same-sign K pi combinations
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

            for (size_t j = 0; j < pid->size(); ++j) {

                // pi candidates
                if (std::abs(pid->at(j)) != 211)
                    continue;

                // Don't pair a particle with itself
                if (i == j)
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
                // Truth matching
                // -----------------------------------

                if (d0Parent->at(i) >= 0 &&
                    d0Parent->at(i) == d0Parent->at(j)) {

                    hMassTrue.Fill(mass);
                }
            }
        }
    }

    // -----------------------------
    // Save histograms
    // -----------------------------

    TFile output("data/d0_reconstructed.root", "RECREATE");

    hMassOS.Write();
    hMassSS.Write();
    hMassTrue.Write();

    output.Close();
    input.Close();

    // -----------------------------
    // Make a quick plot
    // -----------------------------

    TCanvas canvas("canvas", "D0 invariant mass", 800, 600);

    hMassOS.Draw();

    canvas.SaveAs("plots/d0_mass_os.png");

    return 0;
}