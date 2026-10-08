#include "analysisUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// features: look at the candidate tree before training.
// For every feature, overlays signal (truth D0 -> K pi pairs) and
// background (all other opposite-sign K pi pairs), each normalised to
// unit area, and prints a simple separation score.
//
// Usage: ./features [candidateFile] [--tag=name]
//   candidateFile  default data/candidates_charm.root

namespace {

struct Feature {
    const char* name;    // branch name in the Candidates tree
    const char* label;   // x-axis title
    int    nBins;
    double lo, hi;
};

}  // namespace

int main(int argc, char* argv[]) {

    const std::string inputName =
        getInputFile(argc, argv, "data/candidates_charm.root");
    const std::string tag = getTag(argc, argv, inputName);

    gSystem->mkdir("plots", true);

    const std::vector<Feature> features = {
        {"ptK",          "K p_{T} [GeV/c]",               50,  0.0, 10.0},
        {"ptPi",         "#pi p_{T} [GeV/c]",             50,  0.0, 10.0},
        {"etaK",         "K #eta",                        40, -1.1,  1.1},
        {"etaPi",        "#pi #eta",                      40, -1.1,  1.1},
        {"dEta",         "#Delta#eta (K - #pi)",          50, -2.2,  2.2},
        {"dPhi",         "#Delta#phi (K - #pi)",          50, -3.2,  3.2},
        {"dR",           "#DeltaR (K, #pi)",              50,  0.0,  4.0},
        {"openAngle",    "Opening angle [rad]",           50,  0.0,  3.2},
        {"ptBalance",    "|pT_{K} - pT_{#pi}| / (pT_{K} + pT_{#pi})", 50, 0.0, 1.0},
        {"ptHard",       "Harder daughter p_{T} [GeV/c]", 50,  0.0, 10.0},
        {"ptSoft",       "Softer daughter p_{T} [GeV/c]", 50,  0.0,  6.0},
        {"ptSum",        "Scalar p_{T} sum (K + #pi) [GeV/c]", 50, 0.0, 15.0},
        {"pairPt",       "Pair p_{T} [GeV/c]",            50,  0.0, 15.0},
        {"pairEta",      "Pair #eta",                     40, -1.1,  1.1},
        {"cosThetaStar", "cos#theta* (K, pair rest frame)", 40, -1.0, 1.0},

        // Jet association
        {"sameJet",      "K and #pi in the same jet (0/1)", 2,  0.0,  2.0},
        {"jetPt",        "Jet p_{T} [GeV/c]",             50,  0.0, 30.0},
        {"zJet",         "Pair p_{T} / jet p_{T}",        50,  0.0,  1.5},
        {"dRJet",        "#DeltaR (pair, jet axis)",      50,  0.0,  5.0},
        {"jetNConst",    "Jet constituents",              50,  0.0, 60.0},
        {"jetMass",      "Jet mass [GeV/c^{2}]",          50,  0.0,  8.0},
        {"jetEta",       "Jet #eta",                      40, -2.0,  2.0},

        // Isolation
        {"isoRel02",     "Charged isolation, R = 0.2",    50,  0.0,  5.0},
        {"isoRel03",     "Charged isolation, R = 0.3",    50,  0.0,  5.0},
        {"isoRel",       "Charged isolation, R = 0.4",    50,  0.0,  5.0},
        {"isoRel06",     "Charged isolation, R = 0.6",    50,  0.0,  5.0},
        {"isoNeutralRel","Neutral isolation, R = 0.4",    50,  0.0,  5.0},

        // Event level
        {"nCharged",     "Charged multiplicity (|#eta|<1)", 50, 0.0, 50.0},
        {"nJets",        "Number of jets (|#eta|<1)",      8,  0.0,  8.0},
        {"sphericityT",  "Transverse sphericity",         40,  0.0,  1.0},
        {"thrustT",      "Transverse thrust",             40,  0.6,  1.0},

        {"mass",         "M_{K#pi} [GeV/c^{2}] (not a training variable)", 50, 1.6, 2.1}
    };

    // -----------------------------
    // Open input
    // -----------------------------
    TFile input(inputName.c_str(), "READ");

    TTree* tree = nullptr;
    input.GetObject("Candidates", tree);

    if (!tree) {
        std::cerr << "Could not find TTree 'Candidates' in "
                  << inputName << "\n";
        return 1;
    }

    int   isSignal = 0;
    float weight   = 1.0f;
    tree->SetBranchAddress("isSignal", &isSignal);
    tree->SetBranchAddress("weight", &weight);

    std::vector<float> values(features.size(), 0.0f);

    for (size_t k = 0; k < features.size(); ++k) {
        if (tree->SetBranchAddress(features[k].name, &values[k]) < 0) {
            std::cerr << "Branch '" << features[k].name
                      << "' not found in " << inputName << "\n";
            return 1;
        }
    }

    // -----------------------------
    // Histograms
    // -----------------------------
    std::vector<TH1F*> hSig, hBkg;

    for (const auto& f : features) {

        TH1F* s = new TH1F((std::string("sig_") + f.name).c_str(),
            (std::string(";") + f.label + ";Fraction of candidates").c_str(),
            f.nBins, f.lo, f.hi);
        TH1F* b = new TH1F((std::string("bkg_") + f.name).c_str(),
            (std::string(";") + f.label + ";Fraction of candidates").c_str(),
            f.nBins, f.lo, f.hi);

        s->SetDirectory(nullptr);
        b->SetDirectory(nullptr);
        s->Sumw2();
        b->Sumw2();

        hSig.push_back(s);
        hBkg.push_back(b);
    }

    long long nSig = 0, nBkg = 0;

    // -----------------------------
    // Loop over candidates
    // -----------------------------
    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {

        tree->GetEntry(i);

        if (isSignal == 1) ++nSig;
        else               ++nBkg;

        for (size_t k = 0; k < features.size(); ++k) {
            if (isSignal == 1) hSig[k]->Fill(values[k]);
            else               hBkg[k]->Fill(values[k], weight);
        }
    }

    input.Close();

    // -----------------------------
    // Normalise (including under/overflow) and score
    // -----------------------------
    struct Score { std::string name; double sep; };
    std::vector<Score> scores;

    for (size_t k = 0; k < features.size(); ++k) {

        const double iS = hSig[k]->Integral(0, hSig[k]->GetNbinsX() + 1);
        const double iB = hBkg[k]->Integral(0, hBkg[k]->GetNbinsX() + 1);

        if (iS > 0) hSig[k]->Scale(1.0 / iS);
        if (iB > 0) hBkg[k]->Scale(1.0 / iB);

        // Total-variation distance between the two normalised
        // distributions: 0 = identical shapes, 1 = no overlap.
        double sep = 0.0;
        for (int b = 0; b <= hSig[k]->GetNbinsX() + 1; ++b)
            sep += std::abs(hSig[k]->GetBinContent(b) -
                            hBkg[k]->GetBinContent(b));
        sep *= 0.5;

        scores.push_back({features[k].name, sep});
    }

    // -----------------------------
    // Console summary
    // -----------------------------
    std::cout << "\n===== features summary =====\n"
              << "Input:            " << inputName << "\n"
              << "Signal rows:      " << nSig << "\n"
              << "Background rows:  " << nBkg << "\n"
              << "Separation (0 = identical shapes, 1 = no overlap):\n";

    std::sort(scores.begin(), scores.end(),
              [](const Score& a, const Score& b) { return a.sep > b.sep; });

    for (const auto& s : scores)
        std::cout << "  " << s.name << ": " << s.sep << "\n";

    std::cout << "==========================\n";

    // -----------------------------
    // Plot
    // -----------------------------
    gStyle->SetOptStat(0);

    const int nPads = static_cast<int>(features.size()) + 1;   // + legend
    const int nCols = 6;
    const int nRows = (nPads + nCols - 1) / nCols;

    TCanvas canvas("canvas", "Features", 400 * nCols, 330 * nRows);
    canvas.Divide(nCols, nRows);

    for (size_t k = 0; k < features.size(); ++k) {

        canvas.cd(k + 1);

        hSig[k]->SetLineColor(kBlue + 1);
        hSig[k]->SetLineWidth(2);
        hBkg[k]->SetLineColor(kRed + 1);
        hBkg[k]->SetLineWidth(2);

        hSig[k]->SetMaximum(1.3 * std::max(hSig[k]->GetMaximum(),
                                           hBkg[k]->GetMaximum()));
        hSig[k]->SetMinimum(0);

        hSig[k]->Draw("HIST");
        hBkg[k]->Draw("HIST SAME");
    }

    // Legend in the first free pad
    canvas.cd(features.size() + 1);
    TLegend leg(0.1, 0.4, 0.9, 0.6);
    leg.AddEntry(hSig[0], "Signal: truth D^{0}#rightarrowK#pi", "l");
    leg.AddEntry(hBkg[0], "Background: other OS K#pi pairs", "l");
    leg.Draw();

    canvas.SaveAs(("plots/features" + tag + ".png").c_str());

    for (TH1F* h : hSig) delete h;
    for (TH1F* h : hBkg) delete h;

    return 0;
}
