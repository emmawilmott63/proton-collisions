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

    // Consistent colors: opposite-sign blue, same-sign red.
    for (TH1* h : {hOS, hOSw}) {
        h->SetLineColor(kBlue + 1);
        h->SetLineWidth(2);
        h->SetMarkerColor(kBlue + 1);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.55);
    }
    for (TH1* h : {hSS, hSSw}) {
        h->SetLineColor(kRed + 1);
        h->SetLineWidth(2);
        h->SetMarkerColor(kRed + 1);
        h->SetMarkerStyle(24);
        h->SetMarkerSize(0.55);
    }

    const std::string out(outputDir);

    // 1. Unweighted OS and SS comparison.
    TCanvas* c1 = new TCanvas("cMassOSSS", "Unweighted OS vs SS", 900, 700);
    c1->SetLeftMargin(0.13);
    c1->SetBottomMargin(0.12);
    hOS->SetTitle("D^{0} candidate invariant mass;M_{K#pi} (GeV/c^{2});Candidates");
    hOS->Draw("E");
    hSS->Draw("E SAME");
    TLegend* leg1 = new TLegend(0.59, 0.76, 0.87, 0.87);
    leg1->AddEntry(hOS, "Opposite-sign K#pi", "lep");
    leg1->AddEntry(hSS, "Same-sign K#pi", "lep");
    leg1->Draw();
    TLatex label1;
    label1.SetNDC();
    label1.SetTextSize(0.032);
    label1.DrawLatex(0.14, 0.94, "Unweighted candidates");
    c1->SaveAs((out + "/mass_os_vs_ss.png").c_str());
    c1->SaveAs((out + "/mass_os_vs_ss.pdf").c_str());

    // 2. Cross-section-weighted OS and SS comparison.
    TCanvas* c2 = new TCanvas("cMassWeightedOSSS", "Weighted OS vs SS", 900, 700);
    c2->SetLeftMargin(0.13);
    c2->SetBottomMargin(0.12);
    hOSw->SetTitle("Weighted D^{0} candidate invariant mass;M_{K#pi} (GeV/c^{2});Weighted yield (mb)");
    hOSw->Draw("HIST");
    hSSw->Draw("HIST SAME");
    TLegend* leg2 = new TLegend(0.59, 0.76, 0.87, 0.87);
    leg2->AddEntry(hOSw, "Opposite-sign K#pi", "l");
    leg2->AddEntry(hSSw, "Same-sign K#pi", "l");
    leg2->Draw();
    TLatex label2;
    label2.SetNDC();
    label2.SetTextSize(0.032);
    label2.DrawLatex(0.14, 0.94, "Cross-section weighted");
    c2->SaveAs((out + "/mass_weighted_os_vs_ss.png").c_str());
    c2->SaveAs((out + "/mass_weighted_os_vs_ss.pdf").c_str());

    // 3. OS distributions side by side; scales/units are intentionally separate.
    TCanvas* c3 = new TCanvas("cMassOSWeightedComparison", "OS weighted comparison", 1100, 500);
    c3->Divide(2, 1);
    c3->cd(1);
    gPad->SetLeftMargin(0.13);
    hOS->SetTitle("Opposite-sign (unweighted);M_{K#pi} (GeV/c^{2});Candidates");
    hOS->Draw("HIST");
    c3->cd(2);
    gPad->SetLeftMargin(0.13);
    hOSw->SetTitle("Opposite-sign (weighted);M_{K#pi} (GeV/c^{2});Weighted yield (mb)");
    hOSw->Draw("HIST");
    c3->SaveAs((out + "/mass_unweighted_vs_weighted.png").c_str());
    c3->SaveAs((out + "/mass_unweighted_vs_weighted.pdf").c_str());

    // 4. SS distributions side by side.
    TCanvas* c4 = new TCanvas("cMassSSWeightedComparison", "SS weighted comparison", 1100, 500);
    c4->Divide(2, 1);
    c4->cd(1);
    gPad->SetLeftMargin(0.13);
    hSS->SetTitle("Same-sign (unweighted);M_{K#pi} (GeV/c^{2});Candidates");
    hSS->Draw("HIST");
    c4->cd(2);
    gPad->SetLeftMargin(0.13);
    hSSw->SetTitle("Same-sign (weighted);M_{K#pi} (GeV/c^{2});Weighted yield (mb)");
    hSSw->Draw("HIST");
    c4->SaveAs((out + "/mass_ss_unweighted_vs_weighted.png").c_str());
    c4->SaveAs((out + "/mass_ss_unweighted_vs_weighted.pdf").c_str());

    std::cout << "\nSaved plots to " << outputDir << "/:" << std::endl
              << "  mass_os_vs_ss.png" << std::endl
              << "  mass_weighted_os_vs_ss.png" << std::endl
              << "  mass_unweighted_vs_weighted.png" << std::endl
              << "  mass_ss_unweighted_vs_weighted.png" << std::endl;

    input->Close();
    delete input;
}
