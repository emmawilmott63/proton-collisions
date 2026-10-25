#include "tmvaCommon.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/Reader.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// applyTMVA: evaluate the trained BDT on the candidate tree, scan the BDT
// cut, and check the mass shape before / after the cut.
//
// Usage: ./applyTMVA [candidatesFile] [--option=value ...]
//   candidatesFile    default data/candidates_inclusive.root
//
// Options (defaults in brackets):
//   --method=BDT|BDTG   which trained method                  [BDT]
//   --dataset=NAME      TMVA dataset / weight directory       [dataset]
//   --events=test|all   test = ODD event numbers only (never seen in
//                       training); "all" includes training events and is
//                       over-optimistic                       [test]
//   --mass-peak=X       D0 mass for the signal window         [1.86484]
//   --sig-window=X      half-width of the signal window       [0.05]
//   --mass-min=X --mass-max=X   mass histogram range          [1.6, 2.1]
//   --min-bkg=X         ignore cuts leaving less background than this
//                       (weighted yield, full-sample units)   [10]
//   --cut=X             use this BDT cut for the mass plots instead of the
//                       best-significance cut
//   --vars / --drop     MUST match the options used in trainTMVA
//   --out=PATH          histograms   [data/bdt_response<tag>.root]
//
// Yields are weighted counts in the signal window and, for --events=test,
// multiplied by 2 so they refer to the full generated sample.
// S/sqrt(S+B) uses these generator-level counts, not a cross section.
//
// Mass-sculpting check: the BDT never sees "mass", but mass is a function
// of ptK, ptPi, etaK, etaPi, dPhi (all usable inputs), so the BDT can learn
// it. The program prints the fraction of background that lies in the
// signal window before and after the cut. A large increase means the cut is
// shaping the background around the D0 mass; this also happens a little for
// physical reasons (harder pairs sit at higher mass), so treat it as a
// diagnostic, not proof.

int main(int argc, char* argv[]) {

    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    // Histograms must not attach to the input file (it is closed before
    // they are destroyed); they are written explicitly to the output file.
    TH1::AddDirectory(kFALSE);

    const std::string inputName =
        getInputFile(argc, argv, "data/candidates_inclusive.root");
    const std::string tag     = getTag(argc, argv, inputName);
    const std::string dataset = getOption(argc, argv, "dataset", "dataset");
    const std::string method  = getOption(argc, argv, "method", "BDT");
    const std::string events  = getOption(argc, argv, "events", "test");
    const std::string outName =
        getOption(argc, argv, "out", "data/bdt_response" + tag + ".root");

    const double massPeak  = getNumber(argc, argv, "mass-peak", 1.86484);
    const double sigWindow = getNumber(argc, argv, "sig-window", 0.05);
    const double massMin   = getNumber(argc, argv, "mass-min", 1.6);
    const double massMax   = getNumber(argc, argv, "mass-max", 2.1);
    const double minBkg    = getNumber(argc, argv, "min-bkg", 10.0);
    const double fixedCut  = getNumber(argc, argv, "cut", -999.0);

    if (events != "test" && events != "all") {
        std::cerr << "--events must be test or all\n";
        return 1;
    }
    if (events == "all")
        std::cerr << "Warning: --events=all includes the training sample; "
                     "performance will look better than it is.\n";
    const double scale = (events == "test") ? 2.0 : 1.0;

    const std::vector<std::string> vars = getTMVAVars(argc, argv);

    const std::string weightFile =
        dataset + "/weights/TMVAClassification_" + method + ".weights.xml";
    if (gSystem->AccessPathName(weightFile.c_str())) {
        std::cerr << "Cannot find " << weightFile
                  << " (run trainTMVA first, same --dataset/--method)\n";
        return 1;
    }

    gSystem->mkdir("data", true);
    gSystem->mkdir("plots", true);

    // -----------------------------
    // Input tree
    // -----------------------------
    TFile input(inputName.c_str(), "READ");
    if (input.IsZombie()) {
        std::cerr << "Could not open " << inputName << "\n";
        return 1;
    }
    TTree* tree = nullptr;
    input.GetObject("Candidates", tree);
    if (!tree) {
        std::cerr << "No TTree 'Candidates' in " << inputName << "\n";
        return 1;
    }

    std::vector<float> vals(vars.size(), 0.f);
    float b_mass = 0, b_weight = 1;
    int   b_event = 0, b_isSignal = 0;

    for (size_t i = 0; i < vars.size(); ++i)
        if (tree->SetBranchAddress(vars[i].c_str(), &vals[i]) < 0) {
            std::cerr << "Branch '" << vars[i] << "' not found\n";
            return 1;
        }
    tree->SetBranchAddress("mass", &b_mass);
    tree->SetBranchAddress("weight", &b_weight);
    tree->SetBranchAddress("event", &b_event);
    tree->SetBranchAddress("isSignal", &b_isSignal);

    // -----------------------------
    // TMVA reader: same variables, same order, same spectator as training
    // -----------------------------
    TMVA::Reader reader("!Color:!Silent");
    for (size_t i = 0; i < vars.size(); ++i)
        reader.AddVariable(vars[i].c_str(), &vals[i]);
    reader.AddSpectator("mass", &b_mass);
    reader.BookMVA(method.c_str(), weightFile.c_str());

    // -----------------------------
    // Score every selected candidate
    // -----------------------------
    struct Row { float score, mass, w; bool sig; };
    std::vector<Row> rows;

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (events == "test" && (b_event % 2) != 1)
            continue;
        rows.push_back({static_cast<float>(reader.EvaluateMVA(method.c_str())),
                        b_mass, b_weight, b_isSignal == 1});
    }
    std::cout << "Scored " << rows.size() << " candidates (" << events
              << " events)\n";

    // -----------------------------
    // Histograms
    // -----------------------------
    const int nScore = 400;
    const int nMass  = 100;

    TH1D hScoreSig("hScore_sig", ";BDT response;Fraction of candidates",
                   nScore, -1, 1);
    TH1D hScoreBkg("hScore_bkg", ";BDT response;Fraction of candidates",
                   nScore, -1, 1);
    TH1D hSigWin("hScoreWindow_sig", "signal in window", nScore, -1, 1);
    TH1D hBkgWin("hScoreWindow_bkg", "background in window", nScore, -1, 1);

    for (const Row& r : rows) {
        (r.sig ? hScoreSig : hScoreBkg).Fill(r.score, r.w);
        if (std::abs(r.mass - massPeak) < sigWindow)
            (r.sig ? hSigWin : hBkgWin).Fill(r.score, r.w);
    }

    // sAbove[i] / bAbove[i] = window yield with score >= cut(i)
    std::vector<double> cutValue(nScore), sAbove(nScore), bAbove(nScore);
    {
        double s = 0, b = 0;
        for (int bin = nScore; bin >= 1; --bin) {
            s += hSigWin.GetBinContent(bin);
            b += hBkgWin.GetBinContent(bin);
            const int i = bin - 1;
            cutValue[i] = hSigWin.GetXaxis()->GetBinLowEdge(bin);
            sAbove[i] = s * scale;
            bAbove[i] = b * scale;
        }
    }
    const double sTotal = sAbove[0];
    const double bTotal = bAbove[0];

    if (sTotal <= 0 || bTotal <= 0) {
        std::cerr << "No signal or no background in the mass window "
                  << massPeak << " +- " << sigWindow << "\n";
        return 1;
    }

    auto significance = [&](int i) {
        const double s = sAbove[i], b = bAbove[i];
        return (s + b > 0) ? s / std::sqrt(s + b) : 0.0;
    };

    int best = 0;
    for (int i = 0; i < nScore; ++i)
        if (bAbove[i] >= minBkg && significance(i) > significance(best))
            best = i;

    int chosen = best;
    if (fixedCut > -900) {
        chosen = 0;
        for (int i = 0; i < nScore; ++i)
            if (cutValue[i] <= fixedCut)
                chosen = i;
    }
    const double cut = cutValue[chosen];

    // -----------------------------
    // Console tables
    // -----------------------------
    std::cout << std::fixed << std::setprecision(3)
              << "\nSignal window: " << massPeak << " +- " << sigWindow
              << " GeV/c^2\n"
              << "No BDT cut:  S = " << sTotal << "  B = " << bTotal
              << "  S/B = " << sTotal / bTotal
              << "  S/sqrt(S+B) = " << significance(0) << "\n\n"
              << "  eff_S   cut     bkg_rej        S          B        S/B   S/sqrt(S+B)\n";

    for (double target : {0.95, 0.90, 0.80, 0.70, 0.50, 0.30}) {
        int k = 0;
        for (int i = 0; i < nScore; ++i)
            if (sAbove[i] / sTotal >= target)
                k = i;
        std::cout << std::setw(7) << sAbove[k] / sTotal
                  << std::setw(8) << cutValue[k]
                  << std::setw(11) << 1.0 - bAbove[k] / bTotal
                  << std::setw(10) << sAbove[k]
                  << std::setw(11) << bAbove[k]
                  << std::setw(11) << (bAbove[k] > 0 ? sAbove[k] / bAbove[k] : 0)
                  << std::setw(11) << significance(k) << "\n";
    }

    std::cout << "\nBest cut (B >= " << minBkg << "): BDT > " << cutValue[best]
              << "  S/sqrt(S+B) = " << significance(best)
              << "  (eff_S " << sAbove[best] / sTotal
              << ", bkg rejection " << 1.0 - bAbove[best] / bTotal << ")\n";
    if (fixedCut > -900)
        std::cout << "Using --cut: BDT > " << cut << "\n";

    // -----------------------------
    // Mass shape before / after the cut
    // -----------------------------
    TH1D hMassBkgAll("hMass_bkg_all", ";M_{K#pi} [GeV/c^{2}];Normalized",
                     nMass, massMin, massMax);
    TH1D hMassBkgCut("hMass_bkg_cut", ";M_{K#pi} [GeV/c^{2}];Normalized",
                     nMass, massMin, massMax);
    TH1D hMassSigAll("hMass_sig_all", ";M_{K#pi} [GeV/c^{2}];Yield",
                     nMass, massMin, massMax);
    TH1D hMassSigCut("hMass_sig_cut", ";M_{K#pi} [GeV/c^{2}];Yield",
                     nMass, massMin, massMax);

    double bkgAll = 0, bkgAllWin = 0, bkgCutSum = 0, bkgCutWin = 0;
    for (const Row& r : rows) {
        const bool inWin = std::abs(r.mass - massPeak) < sigWindow;
        const bool pass  = r.score >= cut;
        if (r.sig) {
            hMassSigAll.Fill(r.mass, r.w);
            if (pass) hMassSigCut.Fill(r.mass, r.w);
        } else {
            hMassBkgAll.Fill(r.mass, r.w);
            bkgAll += r.w;
            if (inWin) bkgAllWin += r.w;
            if (pass) {
                hMassBkgCut.Fill(r.mass, r.w);
                bkgCutSum += r.w;
                if (inWin) bkgCutWin += r.w;
            }
        }
    }

    std::cout << "\nMass-sculpting check (BDT > " << cut << "):\n"
              << "  background fraction inside the signal window: before cut "
              << std::setprecision(4) << (bkgAll > 0 ? bkgAllWin / bkgAll : 0)
              << ", after cut "
              << (bkgCutSum > 0 ? bkgCutWin / bkgCutSum : 0) << "\n";

    // -----------------------------
    // Output file and plots
    // -----------------------------
    TFile out(outName.c_str(), "RECREATE");
    hScoreSig.Write();
    hScoreBkg.Write();
    hSigWin.Write();
    hBkgWin.Write();
    hMassBkgAll.Write();
    hMassBkgCut.Write();
    hMassSigAll.Write();
    hMassSigCut.Write();

    const std::string stem = "plots/tmva_" + method + tag;

    {   // BDT response
        TH1D* s = static_cast<TH1D*>(hScoreSig.Clone("s_norm"));
        TH1D* b = static_cast<TH1D*>(hScoreBkg.Clone("b_norm"));
        if (s->Integral() > 0) s->Scale(1.0 / s->Integral());
        if (b->Integral() > 0) b->Scale(1.0 / b->Integral());
        s->SetLineColor(kBlue + 1);  s->SetLineWidth(2);
        b->SetLineColor(kRed + 1);   b->SetLineWidth(2);
        TCanvas c("c1", "", 800, 600);
        c.SetLogy();
        s->SetMaximum(2.0 * std::max(s->GetMaximum(), b->GetMaximum()));
        s->Draw("HIST");
        b->Draw("HIST SAME");
        TLegend leg(0.15, 0.75, 0.40, 0.88);
        leg.AddEntry(s, "Signal (true D^{0}#rightarrowK#pi)", "l");
        leg.AddEntry(b, "Background", "l");
        leg.Draw();
        c.SaveAs((stem + "_response.png").c_str());
    }

    {   // S/sqrt(S+B) vs cut
        std::vector<double> x, y;
        for (int i = 0; i < nScore; ++i)
            if (bAbove[i] >= minBkg) {
                x.push_back(cutValue[i]);
                y.push_back(significance(i));
            }
        if (!x.empty()) {
            TGraph g(static_cast<int>(x.size()), x.data(), y.data());
            g.SetTitle(";BDT cut;S/#sqrt{S+B} in signal window");
            g.SetLineWidth(2);
            TCanvas c("c2", "", 800, 600);
            g.Draw("AL");
            c.SaveAs((stem + "_significance.png").c_str());
        }
    }

    {   // background mass shape before / after the cut
        TH1D* a = static_cast<TH1D*>(hMassBkgAll.Clone("a_norm"));
        TH1D* b = static_cast<TH1D*>(hMassBkgCut.Clone("b_norm2"));
        if (a->Integral() > 0) a->Scale(1.0 / a->Integral());
        if (b->Integral() > 0) b->Scale(1.0 / b->Integral());
        a->SetLineColor(kBlack);    a->SetLineWidth(2);
        b->SetLineColor(kRed + 1);  b->SetLineWidth(2);
        a->SetMinimum(0);
        a->SetMaximum(1.4 * std::max(a->GetMaximum(), b->GetMaximum()));
        TCanvas c("c3", "", 800, 600);
        a->Draw("HIST");
        b->Draw("HIST SAME");
        TLegend leg(0.55, 0.75, 0.88, 0.88);
        leg.AddEntry(a, "Background, no cut", "l");
        leg.AddEntry(b, ("Background, BDT > " +
                         std::to_string(cut).substr(0, 5)).c_str(), "l");
        leg.Draw();
        c.SaveAs((stem + "_bkgmass.png").c_str());
    }

    out.Close();
    input.Close();

    std::cout << "\nHistograms: " << outName << "\nPlots:      " << stem
              << "_{response,significance,bkgmass}.png\n";
    return 0;
}
