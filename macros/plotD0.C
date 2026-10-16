// Plot D0 -> K pi invariant-mass histograms from reconstructD0.cc.
//
// Run from the project root:
//   root -l -b -q 'macros/plotD0.C("data/d0_reco_5_100.root")'
//
// Optional output directory:
//   root -l -b -q 'macros/plotD0.C("data/d0_reco_5_100.root", "plots")'
//
// Expected histograms: hMass_OS, hMass_SS, hMass_OS_weighted,
//                      hMass_SS_weighted
//
// Outputs:
//   mass_weighted-recon.png       weighted OS, SS, and OS-SS
//   mass_unweighted-recon.png     unweighted OS, SS, and OS-SS
//   kpi_os_recon.png              OS candidates, weighted/unweighted side by side
//   kpi_ss_recon.png              SS candidates, weighted/unweighted side by side

#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"

void plotD0(const char* inputName = "data/d0_reco_5_100.root",
            const char* outputDir = "plots")
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleSize(0.045, "XY");
    gStyle->SetLabelSize(0.04, "XY");
    gStyle->SetLegendBorderSize(0);
    gSystem->mkdir(outputDir, true);

    TFile* input = TFile::Open(inputName, "READ");
    if (!input || input->IsZombie()) {
        std::cerr << "ERROR: Could not open input file: " << inputName << std::endl;
        if (input) { input->Close(); delete input; }
        return;
    }

    auto getHistogram = [&](const char* name) -> TH1* {
        TH1* h = dynamic_cast<TH1*>(input->Get(name));
        if (!h)
            std::cerr << "ERROR: Histogram '" << name << "' was not found in "
                      << inputName << std::endl;
        return h;
    };

    TH1* hOS  = getHistogram("hMass_OS");
    TH1* hSS  = getHistogram("hMass_SS");
    TH1* hOSw = getHistogram("hMass_OS_weighted");
    TH1* hSSw = getHistogram("hMass_SS_weighted");

    if (!hOS || !hSS || !hOSw || !hSSw) {
        std::cerr << "Check histogram names in reconstructD0.cc and the input file."
                  << std::endl;
        input->Close();
        delete input;
        return;
    }

    // Make independent clones for display and subtraction, leaving input histograms untouched.
    TH1* hOSDisplay = dynamic_cast<TH1*>(hOS->Clone("hOSDisplay"));
    TH1* hSSDisplay = dynamic_cast<TH1*>(hSS->Clone("hSSDisplay"));
    TH1* hOSwDisplay = dynamic_cast<TH1*>(hOSw->Clone("hOSwDisplay"));
    TH1* hSSwDisplay = dynamic_cast<TH1*>(hSSw->Clone("hSSwDisplay"));
    TH1* hDiff = dynamic_cast<TH1*>(hOS->Clone("hOSMinusSS"));
    TH1* hDiffW = dynamic_cast<TH1*>(hOSw->Clone("hOSMinusSSWeighted"));
    hDiff->SetDirectory(nullptr);
    hDiffW->SetDirectory(nullptr);
    hOSDisplay->SetDirectory(nullptr);
    hSSDisplay->SetDirectory(nullptr);
    hOSwDisplay->SetDirectory(nullptr);
    hSSwDisplay->SetDirectory(nullptr);

    // Subtract same-sign from opposite-sign, bin by bin.
    // This is a useful control plot, but is not automatically a calibrated
    // combinatorial-background subtraction: OS/SS normalization may be needed.
    hDiff->Add(hSS, -1.0);
    hDiffW->Add(hSSw, -1.0);

    for (TH1* h : {hOSDisplay, hOSwDisplay}) {
        h->SetLineColor(kBlue + 1);
        h->SetLineWidth(2);
        h->SetMarkerColor(kBlue + 1);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.55);
    }
    for (TH1* h : {hSSDisplay, hSSwDisplay}) {
        h->SetLineColor(kRed + 1);
        h->SetLineWidth(2);
        h->SetMarkerColor(kRed + 1);
        h->SetMarkerStyle(24);
        h->SetMarkerSize(0.55);
    }
    hDiff->SetLineColor(kGreen + 2);
    hDiff->SetLineWidth(2);
    hDiffW->SetLineColor(kGreen + 2);
    hDiffW->SetLineWidth(2);

    const std::string out(outputDir);

    // 1. Weighted mass plot: OS and SS on the left; OS - SS on the right.
    TCanvas* cWeighted = new TCanvas("cMassWeightedRecon", "Weighted mass reconstruction", 1300, 550);
    cWeighted->Divide(2, 1);

    cWeighted->cd(1);
    gPad->SetLeftMargin(0.13);
    hOSwDisplay->SetTitle("Weighted K#pi candidates;M_{K#pi} (GeV/c^{2});Weighted yield (mb)");
    hOSwDisplay->Draw("HIST");
    hSSwDisplay->Draw("HIST SAME");
    TLegend* legW = new TLegend(0.58, 0.76, 0.87, 0.88);
    legW->AddEntry(hOSwDisplay, "Opposite-sign (OS)", "l");
    legW->AddEntry(hSSwDisplay, "Same-sign (SS)", "l");
    legW->Draw();

    cWeighted->cd(2);
    gPad->SetLeftMargin(0.13);
    hDiffW->SetTitle("Weighted OS - SS;M_{K#pi} (GeV/c^{2});OS - SS (mb)");
    hDiffW->Draw("HIST");
    cWeighted->SaveAs((out + "/mass_weighted-recon.png").c_str());
    cWeighted->SaveAs((out + "/mass_weighted-recon.pdf").c_str());

    // 2. Unweighted mass plot: OS and SS on the left; OS - SS on the right.
    TCanvas* cUnweighted = new TCanvas("cMassUnweightedRecon", "Unweighted mass reconstruction", 1300, 550);
    cUnweighted->Divide(2, 1);

    cUnweighted->cd(1);
    gPad->SetLeftMargin(0.13);
    hOSDisplay->SetTitle("K#pi candidates;M_{K#pi} (GeV/c^{2});Candidates");
    hOSDisplay->Draw("E");
    hSSDisplay->Draw("E SAME");
    TLegend* legU = new TLegend(0.58, 0.76, 0.87, 0.88);
    legU->AddEntry(hOSDisplay, "Opposite-sign (OS)", "lep");
    legU->AddEntry(hSSDisplay, "Same-sign (SS)", "lep");
    legU->Draw();

    cUnweighted->cd(2);
    gPad->SetLeftMargin(0.13);
    hDiff->SetTitle("Unweighted OS - SS;M_{K#pi} (GeV/c^{2});OS - SS candidates");
    hDiff->Draw("HIST");
    cUnweighted->SaveAs((out + "/mass_unweighted-recon.png").c_str());
    cUnweighted->SaveAs((out + "/mass_unweighted-recon.pdf").c_str());

    // 3. OS candidates, weighted and unweighted side by side.
    TCanvas* cOS = new TCanvas("cKpiOSRecon", "OS candidates", 1200, 520);
    cOS->Divide(2, 1);
    cOS->cd(1);
    gPad->SetLeftMargin(0.13);
    hOSDisplay->SetTitle("Opposite-sign candidates (unweighted);M_{K#pi} (GeV/c^{2});Candidates");
    hOSDisplay->Draw("HIST");
    cOS->cd(2);
    gPad->SetLeftMargin(0.13);
    hOSwDisplay->SetTitle("Opposite-sign candidates (weighted);M_{K#pi} (GeV/c^{2});Weighted yield (mb)");
    hOSwDisplay->Draw("HIST");
    cOS->SaveAs((out + "/kpi_os_recon.png").c_str());
    cOS->SaveAs((out + "/kpi_os_recon.pdf").c_str());

    // 4. SS candidates, weighted and unweighted side by side.
    TCanvas* cSS = new TCanvas("cKpiSSRecon", "SS candidates", 1200, 520);
    cSS->Divide(2, 1);
    cSS->cd(1);
    gPad->SetLeftMargin(0.13);
    hSSDisplay->SetTitle("Same-sign candidates (unweighted);M_{K#pi} (GeV/c^{2});Candidates");
    hSSDisplay->Draw("HIST");
    cSS->cd(2);
    gPad->SetLeftMargin(0.13);
    hSSwDisplay->SetTitle("Same-sign candidates (weighted);M_{K#pi} (GeV/c^{2});Weighted yield (mb)");
    hSSwDisplay->Draw("HIST");
    cSS->SaveAs((out + "/kpi_ss_recon.png").c_str());
    cSS->SaveAs((out + "/kpi_ss_recon.pdf").c_str());

    std::cout << "\nSaved plots to " << outputDir << "/:" << std::endl
              << "  mass_weighted-recon.png" << std::endl
              << "  mass_unweighted-recon.png" << std::endl
              << "  kpi_os_recon.png" << std::endl
              << "  kpi_ss_recon.png" << std::endl;
    std::cout << "Note: OS - SS is a simple bin-by-bin subtraction; "
                 "background normalization may be needed for physics interpretation."
              << std::endl;

    input->Close();
    delete input;
}
