#include "tmvaCommon.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/Reader.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
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
//   --events=test|train|all
//                       test  = ODD event numbers only (never seen in
//                               training)
//                       train = EVEN event numbers (the training sample;
//                               use it to check for overtraining by
//                               comparing the ROC integral with test)
//                       all   = everything (over-optimistic)  [test]
//   --sig-window=X      half-width of the signal window       [auto]
//   --mass-peak=X       centre of the signal window           [auto]
//   --window-nsigma=X   auto window = X * robust sigma of the signal
//                       mass (16-84% quantiles, true D0 -> K pi rows,
//                       median as centre)                     [2]
//   --min-window=X      floor for the auto window             [0.005]
//   --mass-min=X --mass-max=X   mass histogram range          [1.6, 2.1]
//   --lumi-pb=X         integrated luminosity for the yields, in pb^-1
//                       (needs xsWeight in the candidate tree) [1]
//   --min-bkg-rows=N    ignore cuts leaving fewer than N unweighted
//                       background rows in the window (avoids picking a
//                       cut in a statistical fluctuation)      [20]
//   --cut=X             use this BDT cut for the mass plots instead of the
//                       best-significance cut
//   --presel=EXPR       row selection; MUST match --presel used in trainTMVA
//                       (rows outside the training selection are not scored)
//   --vars / --drop / --drop-group
//                       MUST match the options used in trainTMVA
//   --out=PATH          histograms   [TMVA-output/results/bdt_response<tag>.root]
//
// Yields are expected counts in the signal window for the chosen
// luminosity: weight * xsWeight [mb] * L [mb^-1]. For --events=test they
// are multiplied by 2 because only half the events are scored.
// If the candidate tree has no xsWeight, yields are plain weighted counts
// of the generated sample.
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

    if (!checkOptions(argc, argv,
            {"method", "dataset", "events", "out", "tag", "mass-peak",
             "sig-window", "window-nsigma", "min-window", "mass-min",
             "mass-max", "lumi-pb", "min-bkg-rows", "cut", "presel", "vars",
             "drop", "drop-group"}))
        return 1;

    const std::string inputName =
        getInputFile(argc, argv, "data/candidates_inclusive.root");
    const std::string tag     = getTag(argc, argv, inputName);
    const std::string dataset = getOption(argc, argv, "dataset", "dataset");
    const std::string method  = getOption(argc, argv, "method", "BDT");
    const std::string events  = getOption(argc, argv, "events", "test");
    const std::string outName =
        getOption(argc, argv, "out",
                  outDir("results") + "/bdt_response" + tag + ".root");

    // Signal window: by default taken from the signal mass distribution
    // itself (median +- nSigma * robust sigma), so it follows the smearing.
    const double massPeakOpt  = getNumber(argc, argv, "mass-peak", -1.0);
    const double sigWindowOpt = getNumber(argc, argv, "sig-window", -1.0);
    const double windowNSigma = getNumber(argc, argv, "window-nsigma", 2.0);
    const double minWindow    = getNumber(argc, argv, "min-window", 0.005);
    const double massMin   = getNumber(argc, argv, "mass-min", 1.6);
    const double massMax   = getNumber(argc, argv, "mass-max", 2.1);
    const double lumiPb    = getNumber(argc, argv, "lumi-pb", 1.0);
    const double minBkgRows = getNumber(argc, argv, "min-bkg-rows", 20.0);
    const double fixedCut  = getNumber(argc, argv, "cut", -999.0);

    if (events != "test" && events != "train" && events != "all") {
        std::cerr << "--events must be test, train or all\n";
        return 1;
    }
    if (events != "test")
        std::cerr << "Warning: --events=" << events << " includes training "
                     "events; performance will look better than it is.\n";
    // test and train each hold half of the events, so double the yields
    const double scale = (events == "all") ? 1.0 : 2.0;

    const std::vector<std::string> vars = getTMVAVars(argc, argv);
    const std::string presel = getOption(argc, argv, "presel", "");

    const std::string weightFile =
        weightFilePath(dataset, method);
    if (gSystem->AccessPathName(weightFile.c_str())) {
        std::cerr << "Cannot find " << weightFile
                  << " (run trainTMVA first, same --dataset/--method; "
                     "results from before the TMVA-output/ folder existed can "
                     "be moved with: make migrate-output)\n";
        return 1;
    }

    gSystem->mkdir(outDir("results").c_str(), true);
    gSystem->mkdir(outDir("plots").c_str(), true);

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

    if (!presel.empty()) {
        gROOT->cd();   // the selected copy lives in memory
        tree = tree->CopyTree(presel.c_str());
        std::cout << "Preselection '" << presel << "': "
                  << tree->GetEntries() << " rows\n";
    }

    std::vector<float> vals(vars.size(), 0.f);
    float b_mass = 0, b_weight = 1;
    double b_xsWeight = 1.0;
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

    // 1 pb^-1 = 1e9 mb^-1 (xsWeight is in mb per generated event)
    const bool hasXs = tree->GetBranch("xsWeight") != nullptr;
    const double lumiMb = lumiPb * 1.0e9;
    if (hasXs)
        tree->SetBranchAddress("xsWeight", &b_xsWeight);
    else
        std::cerr << "Warning: no xsWeight branch; yields are generated-"
                     "sample counts, not luminosity-scaled.\n";

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
    struct Row { float score, mass; double w; float wFlat; bool sig; };
    std::vector<Row> rows;

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (events == "test" && (b_event % 2) != 1)
            continue;
        if (events == "train" && (b_event % 2) != 0)
            continue;
        double w = b_weight;
        if (hasXs)
            w *= b_xsWeight * lumiMb;
        rows.push_back({static_cast<float>(reader.EvaluateMVA(method.c_str())),
                        b_mass, w, b_weight, b_isSignal == 1});
    }
    std::cout << "Scored " << rows.size() << " candidates (" << events
              << " events)\n";

    // -----------------------------
    // Signal window from the signal mass distribution
    // -----------------------------
    double massPeak = 1.86484, sigWindow = 0.05, sigmaEff = 0.0;
    {
        std::vector<double> ms;
        for (const Row& r : rows)
            if (r.sig && std::abs(r.mass - 1.86484) < 0.3)
                ms.push_back(r.mass);
        if (ms.size() >= 20) {
            std::sort(ms.begin(), ms.end());
            const auto q = [&](double f) {
                return ms[static_cast<size_t>(f * (ms.size() - 1))];
            };
            massPeak = q(0.5);
            sigmaEff = 0.5 * (q(0.84) - q(0.16));
            sigWindow = std::max(windowNSigma * sigmaEff, minWindow);
        } else {
            std::cerr << "Warning: fewer than 20 signal rows; using default "
                         "window 1.86484 +- 0.05\n";
        }
        if (massPeakOpt > 0)  massPeak  = massPeakOpt;
        if (sigWindowOpt > 0) sigWindow = sigWindowOpt;
        std::cout << "Signal mass: median " << std::setprecision(5)
                  << massPeak << ", robust sigma " << sigmaEff
                  << " GeV/c^2 -> window +- " << sigWindow << "\n";
    }

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
    TH1D hBkgRows("hScoreWindow_bkg_rows",
                  "background rows in window (unweighted)", nScore, -1, 1);
    // flat-weight response (training weights, all masses) for the ROC integral
    TH1D hAucSig("hAuc_sig", "", nScore, -1, 1);
    TH1D hAucBkg("hAuc_bkg", "", nScore, -1, 1);
    // errors from the Monte Carlo weights (a few heavy rows dominate)
    hSigWin.Sumw2();
    hBkgWin.Sumw2();

    for (const Row& r : rows) {
        (r.sig ? hScoreSig : hScoreBkg).Fill(r.score, r.w);
        (r.sig ? hAucSig : hAucBkg).Fill(r.score, r.wFlat);
        if (std::abs(r.mass - massPeak) < sigWindow) {
            (r.sig ? hSigWin : hBkgWin).Fill(r.score, r.w);
            if (!r.sig)
                hBkgRows.Fill(r.score);
        }
    }

    // ROC integral = P(signal score > background score), flat weights
    double auc = 0.0;
    {
        const double sTot = hAucSig.Integral(), bTot = hAucBkg.Integral();
        double sAb = 0.0;
        for (int bin = nScore; bin >= 1; --bin) {
            const double sb = hAucSig.GetBinContent(bin);
            const double bb = hAucBkg.GetBinContent(bin);
            if (sTot > 0 && bTot > 0)
                auc += (bb / bTot) * (sAb + 0.5 * sb) / sTot;
            sAb += sb;
        }
    }
    std::cout << "ROC integral (test sample, flat weights): "
              << std::setprecision(4) << auc << "\n";

    // sAbove[i] / bAbove[i] = window yield with score >= cut(i)
    std::vector<double> cutValue(nScore), sAbove(nScore), bAbove(nScore);
    std::vector<double> bRowsAbove(nScore), sErr2(nScore), bErr2(nScore);
    {
        double s = 0, b = 0, nb = 0, se2 = 0, be2 = 0;
        for (int bin = nScore; bin >= 1; --bin) {
            s += hSigWin.GetBinContent(bin);
            b += hBkgWin.GetBinContent(bin);
            nb += hBkgRows.GetBinContent(bin);
            se2 += hSigWin.GetBinError(bin) * hSigWin.GetBinError(bin);
            be2 += hBkgWin.GetBinError(bin) * hBkgWin.GetBinError(bin);
            const int i = bin - 1;
            bRowsAbove[i] = nb;
            sErr2[i] = se2 * scale * scale;
            bErr2[i] = be2 * scale * scale;
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

    // Monte Carlo statistical error of S/sqrt(S+B), propagated from the
    // summed squared weights of the signal and background in the window
    auto significanceErr = [&](int i) {
        const double s = sAbove[i], b = bAbove[i], d = s + b;
        if (d <= 0)
            return 0.0;
        const double d15 = std::pow(d, 1.5);
        const double dS = (s + 2.0 * b) / (2.0 * d15);
        const double dB = -s / (2.0 * d15);
        return std::sqrt(dS * dS * sErr2[i] + dB * dB * bErr2[i]);
    };

    int best = 0;
    for (int i = 0; i < nScore; ++i)
        if (bRowsAbove[i] >= minBkgRows && significance(i) > significance(best))
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
              << "Yields: " << (hasXs ? "expected counts for L = " +
                                   std::to_string(lumiPb) + " pb^-1"
                                 : std::string("generated-sample counts"))
              << "\n"
              << "No BDT cut:  S = " << sTotal << "  B = " << bTotal
              << "  S/B = " << sTotal / bTotal
              << "  S/sqrt(S+B) = " << significance(0) << " +- "
              << significanceErr(0) << "\n\n"
              << std::setw(7) << "eff_S" << std::setw(8) << "cut"
              << std::setw(11) << "bkg_rej" << std::setw(14) << "S"
              << std::setw(16) << "B" << std::setw(11) << "S/B"
              << std::setw(13) << "S/sqrt(S+B)" << "\n";

    for (double target : {0.95, 0.90, 0.80, 0.70, 0.50, 0.30}) {
        int k = 0;
        for (int i = 0; i < nScore; ++i)
            if (sAbove[i] / sTotal >= target)
                k = i;
        std::cout << std::setw(7) << sAbove[k] / sTotal
                  << std::setw(8) << cutValue[k]
                  << std::setw(11) << 1.0 - bAbove[k] / bTotal
                  << std::setw(14) << sAbove[k]
                  << std::setw(16) << bAbove[k]
                  << std::setw(11) << (bAbove[k] > 0 ? sAbove[k] / bAbove[k] : 0)
                  << std::setw(13) << significance(k) << "\n";
    }

    std::cout << "\nBest cut (>= " << minBkgRows << " background rows): BDT > " << cutValue[best]
              << "  S/sqrt(S+B) = " << significance(best) << " +- "
              << significanceErr(best) << " (MC stat.)"
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

    hMassBkgAll.Sumw2();
    hMassBkgCut.Sumw2();
    double bkgAll = 0, bkgAllWin = 0, bkgCutSum = 0, bkgCutWin = 0;
    double bkgCutSumW2 = 0;
    long long bkgCutRows = 0;
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
                bkgCutSumW2 += r.w * r.w;
                ++bkgCutRows;
                if (inWin) bkgCutWin += r.w;
            }
        }
    }

    std::cout << "\nMass-sculpting check (BDT > " << cut << "):\n"
              << "  background fraction inside the signal window: before cut "
              << std::setprecision(4) << (bkgAll > 0 ? bkgAllWin / bkgAll : 0)
              << ", after cut "
              << (bkgCutSum > 0 ? bkgCutWin / bkgCutSum : 0) << "\n"
              << "  background rows passing: " << bkgCutRows
              << ", effective rows (sumW^2/sumW2): "
              << std::setprecision(1)
              << (bkgCutSumW2 > 0 ? bkgCutSum * bkgCutSum / bkgCutSumW2 : 0)
              << "  (a few heavy low-pThat rows can dominate)\n";

    // -----------------------------
    // Output file and plots
    // -----------------------------
    TFile out(outName.c_str(), "RECREATE");
    hScoreSig.Write();
    hScoreBkg.Write();
    hSigWin.Write();
    hBkgWin.Write();
    hBkgRows.Write();
    hMassBkgAll.Write();
    hMassBkgCut.Write();
    hMassSigAll.Write();
    hMassSigCut.Write();

    const std::string stem = outDir("plots") + "/tmva_" + method + tag;

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
        std::vector<double> x, y, ey;
        for (int i = 0; i < nScore; ++i)
            if (bRowsAbove[i] >= minBkgRows) {
                x.push_back(cutValue[i]);
                y.push_back(significance(i));
                ey.push_back(significanceErr(i));
            }
        if (!x.empty()) {
            TGraphErrors g(static_cast<int>(x.size()), x.data(), y.data(),
                           nullptr, ey.data());
            g.SetTitle(";BDT cut;S/#sqrt{S+B} in signal window");
            g.SetLineWidth(2);
            g.SetFillColorAlpha(kGray + 1, 0.4);   // band = MC statistical error
            TCanvas c("c2", "", 800, 600);
            g.Draw("A3");
            g.Draw("LX SAME");
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
        a->Draw("E");
        b->Draw("E SAME");
        TLegend leg(0.55, 0.75, 0.88, 0.88);
        leg.AddEntry(a, "Background, no cut", "l");
        leg.AddEntry(b, ("Background, BDT > " +
                         std::to_string(cut).substr(0, 5)).c_str(), "l");
        leg.Draw();
        c.SaveAs((stem + "_bkgmass.png").c_str());
    }

    // One summary line per (tag, method), for comparing configurations
    // (e.g. smearing levels): TMVA-output/results/tmva_summary.csv. Re-running the same
    // tag and method replaces its row. For an overtraining check give the
    // --events=train run its own --tag (e.g. --tag=smear0p03_train) so it
    // does not replace the test row.
    {
        std::ostringstream row;
        row << std::setprecision(6) << tag.substr(1) << "," << method << ","
            << events << "," << (hasXs ? lumiPb : 0.0) << "," << massPeak
            << "," << sigWindow << "," << sigmaEff << "," << sTotal << ","
            << bTotal << "," << significance(0) << "," << cutValue[best]
            << "," << sAbove[best] << "," << bAbove[best] << ","
            << significance(best) << "," << significanceErr(best) << ","
            << sAbove[best] / sTotal << "," << 1.0 - bAbove[best] / bTotal
            << "," << auc;
        upsertCsvRow(summaryCsvPath(),
            "tag,method,events,lumi_pb,mass_peak,window,sigma_eff,S_nocut,"
            "B_nocut,signif_nocut,best_cut,S_best,B_best,signif_best,"
            "signif_best_err,effS_best,bkgrej_best,auc",
            row.str(), 2);
    }

    out.Close();
    input.Close();

    std::cout << "\nHistograms: " << outName << "\nPlots:      " << stem
              << "_{response,significance,bkgmass}.png\n";
    return 0;
}
