#include "analysisUtils.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// buildtree: turn the particle tree into a CANDIDATE tree for TMVA.
// One row per opposite-sign K pi pair that passes the preselection.
//
// Usage: ./buildtree [inputFile] [--option=value ...]
//   inputFile        default data/d0_charm.root
//
// Options (defaults in brackets):
//   --out=PATH            output file          [data/candidates<tag>.root]
//   --tag=NAME            output tag           [derived from the input name]
//   --min-daughter-pt=X   preselection         [0.5]   GeV/c
//   --max-eta=X           preselection         [1.0]
//   --mass-min=X          preselection window  [1.6]   GeV/c^2
//   --mass-max=X          preselection window  [2.1]   GeV/c^2
//   --smear-pt=X          relative pT smearing of K and pi   [0]
//   --smear-eta=X         absolute eta smearing              [0]
//   --smear-phi=X         absolute phi smearing (rad)        [0]
//   --bkg-keep=X          fraction of background rows kept   [1.0]
//                         (kept rows get weight 1/X; signal weight is 1)
//   --seed=N              random seed                        [12345]
//
// Branches written:
//   label / bookkeeping : event, isSignal, weight
//   NOT for training    : mass   (training on it would sculpt the peak)
//   features            : ptK, ptPi, etaK, etaPi, dEta, dPhi, dR,
//                         openAngle, ptBalance, pairPt, pairEta,
//                         cosThetaStar, jetPt, zJet, dRJet, isoRel

namespace {

const double kPi     = 3.14159265358979323846;
const double kMassK  = 0.493677;
const double kMassPi = 0.13957039;

double deltaPhi(double a, double b) {
    double d = a - b;
    while (d >  kPi) d -= 2 * kPi;
    while (d < -kPi) d += 2 * kPi;
    return d;
}

struct Track {
    double pt = 0, eta = 0, phi = 0;
    double px = 0, py = 0, pz = 0, e = 0;
};

// Build a track from (pt, eta, phi) and a mass hypothesis.
Track makeTrack(double pt, double eta, double phi, double mass) {
    Track t;
    t.pt  = pt;
    t.eta = eta;
    t.phi = phi;
    t.px  = pt * std::cos(phi);
    t.py  = pt * std::sin(phi);
    t.pz  = pt * std::sinh(eta);
    t.e   = std::sqrt(t.px * t.px + t.py * t.py + t.pz * t.pz + mass * mass);
    return t;
}

}  // namespace

int main(int argc, char* argv[]) {

    const std::string inputName =
        getInputFile(argc, argv, "data/d0_charm.root");
    const std::string tag = getTag(argc, argv, inputName);

    const std::string outName =
        getOption(argc, argv, "out", "data/candidates" + tag + ".root");

    // -----------------------------
    // Options
    // -----------------------------
    const double minDaughterPt = getNumber(argc, argv, "min-daughter-pt", 0.5);
    const double maxAbsEta     = getNumber(argc, argv, "max-eta", 1.0);
    const double massMin       = getNumber(argc, argv, "mass-min", 1.6);
    const double massMax       = getNumber(argc, argv, "mass-max", 2.1);

    const double smearPt  = getNumber(argc, argv, "smear-pt", 0.0);
    const double smearEta = getNumber(argc, argv, "smear-eta", 0.0);
    const double smearPhi = getNumber(argc, argv, "smear-phi", 0.0);

    const double bkgKeep = getNumber(argc, argv, "bkg-keep", 1.0);
    const int    seed    = static_cast<int>(getNumber(argc, argv, "seed", 12345));

    if (bkgKeep <= 0.0 || bkgKeep > 1.0) {
        std::cerr << "--bkg-keep must be in (0, 1]\n";
        return 1;
    }

    // Jets used for the jet-based features
    const double jetR     = 0.4;
    const double jetPtMin = 2.0;
    const double isoR     = 0.4;   // cone for the isolation feature

    TRandom3 rng(seed);

    gSystem->mkdir("data", true);

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
    // Output tree
    // -----------------------------
    TFile output(outName.c_str(), "RECREATE");
    TTree out("Candidates", "K pi candidates");

    int   b_event = 0, b_isSignal = 0;
    float b_weight = 1, b_mass = 0;
    float b_ptK = 0, b_ptPi = 0, b_etaK = 0, b_etaPi = 0;
    float b_dEta = 0, b_dPhi = 0, b_dR = 0;
    float b_openAngle = 0, b_ptBalance = 0;
    float b_pairPt = 0, b_pairEta = 0, b_cosThetaStar = 0;
    float b_jetPt = 0, b_zJet = 0, b_dRJet = 0, b_isoRel = 0;

    out.Branch("event", &b_event);
    out.Branch("isSignal", &b_isSignal);
    out.Branch("weight", &b_weight);
    out.Branch("mass", &b_mass);

    out.Branch("ptK", &b_ptK);
    out.Branch("ptPi", &b_ptPi);
    out.Branch("etaK", &b_etaK);
    out.Branch("etaPi", &b_etaPi);
    out.Branch("dEta", &b_dEta);
    out.Branch("dPhi", &b_dPhi);
    out.Branch("dR", &b_dR);
    out.Branch("openAngle", &b_openAngle);
    out.Branch("ptBalance", &b_ptBalance);
    out.Branch("pairPt", &b_pairPt);
    out.Branch("pairEta", &b_pairEta);
    out.Branch("cosThetaStar", &b_cosThetaStar);
    out.Branch("jetPt", &b_jetPt);
    out.Branch("zJet", &b_zJet);
    out.Branch("dRJet", &b_dRJet);
    out.Branch("isoRel", &b_isoRel);

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);

    long long nSignal = 0;
    long long nBackground = 0;
    long long nEvents = 0;

    // -----------------------------
    // Event loop
    // -----------------------------
    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {

        tree->GetEntry(iEvent);
        ++nEvents;

        const size_t n = pid->size();

        // -------------------------------------------------
        // Jets from the (unsmeared) final-state particles
        // -------------------------------------------------
        std::vector<fastjet::PseudoJet> jetInputs;
        for (size_t k = 0; k < n; ++k) {
            const int absId = std::abs(pid->at(k));
            if (absId == 12 || absId == 14 || absId == 16)
                continue;
            jetInputs.emplace_back(px->at(k), py->at(k),
                                   pz->at(k), energy->at(k));
        }

        fastjet::ClusterSequence cs(jetInputs, jetDef);
        std::vector<fastjet::PseudoJet> jets =
            fastjet::sorted_by_pt(cs.inclusive_jets(jetPtMin));

        // Transverse momentum and direction of every particle
        // (unsmeared; used for the isolation sum)
        std::vector<double> pt(n), eta(n), phi(n);
        for (size_t k = 0; k < n; ++k) {
            pt[k]  = std::hypot(px->at(k), py->at(k));
            eta[k] = (pt[k] > 1e-6) ? std::asinh(pz->at(k) / pt[k]) : 0.0;
            phi[k] = std::atan2(py->at(k), px->at(k));
        }

        // -------------------------------------------------
        // Smear each K and pi candidate ONCE per event
        // -------------------------------------------------
        std::vector<Track> trk(n);
        for (size_t k = 0; k < n; ++k) {

            const int absId = std::abs(pid->at(k));
            if (absId != 321 && absId != 211)
                continue;

            double p  = pt[k];
            double e  = eta[k];
            double ph = phi[k];

            if (smearPt  > 0) p  *= (1.0 + rng.Gaus(0.0, smearPt));
            if (smearEta > 0) e  += rng.Gaus(0.0, smearEta);
            if (smearPhi > 0) ph += rng.Gaus(0.0, smearPhi);

            if (p < 1e-3) p = 1e-3;   // keep the track physical

            trk[k] = makeTrack(p, e, ph,
                               absId == 321 ? kMassK : kMassPi);
        }

        // -------------------------------------------------
        // Pair every K with every pi
        // -------------------------------------------------
        for (size_t i = 0; i < n; ++i) {

            if (std::abs(pid->at(i)) != 321)
                continue;

            const Track& K = trk[i];

            if (K.pt < minDaughterPt || std::abs(K.eta) > maxAbsEta)
                continue;

            for (size_t j = 0; j < n; ++j) {

                if (std::abs(pid->at(j)) != 211)
                    continue;

                // opposite-sign pairs only
                if (charge->at(i) * charge->at(j) >= 0)
                    continue;

                const Track& P = trk[j];

                if (P.pt < minDaughterPt || std::abs(P.eta) > maxAbsEta)
                    continue;

                // Pair kinematics
                const double sumPx = K.px + P.px;
                const double sumPy = K.py + P.py;
                const double sumPz = K.pz + P.pz;
                const double sumE  = K.e + P.e;

                const double m2 = sumE * sumE - sumPx * sumPx
                                  - sumPy * sumPy - sumPz * sumPz;
                if (m2 <= 0)
                    continue;

                const double mass = std::sqrt(m2);
                if (mass < massMin || mass > massMax)
                    continue;

                // Truth label: both from the same D0 -> K pi decay
                const bool isSignal =
                    d0Parent->at(i) >= 0 &&
                    d0Parent->at(i) == d0Parent->at(j);

                // Optional background prescale (weights restore the
                // original proportions)
                double weight = 1.0;
                if (!isSignal) {
                    if (bkgKeep < 1.0 && rng.Uniform() > bkgKeep)
                        continue;
                    weight = 1.0 / bkgKeep;
                }

                const double pairPt  = std::hypot(sumPx, sumPy);
                const double pairPhi = std::atan2(sumPy, sumPx);
                const double pairEta =
                    (pairPt > 1e-6) ? std::asinh(sumPz / pairPt) : 0.0;

                // Angles between the daughters
                const double dEta = K.eta - P.eta;
                const double dPhi = deltaPhi(K.phi, P.phi);
                const double dR   = std::hypot(dEta, dPhi);

                const double pK  = std::sqrt(K.px * K.px + K.py * K.py + K.pz * K.pz);
                const double pPi = std::sqrt(P.px * P.px + P.py * P.py + P.pz * P.pz);
                double cosOpen =
                    (K.px * P.px + K.py * P.py + K.pz * P.pz) / (pK * pPi);
                cosOpen = std::max(-1.0, std::min(1.0, cosOpen));
                const double openAngle = std::acos(cosOpen);

                const double ptBalance =
                    std::abs(K.pt - P.pt) / (K.pt + P.pt);

                // Decay angle of the K in the pair rest frame,
                // measured relative to the pair direction
                double cosThetaStar = 0.0;
                const double pPair =
                    std::sqrt(sumPx * sumPx + sumPy * sumPy + sumPz * sumPz);
                if (pPair > 1e-9) {
                    const double nx = sumPx / pPair;
                    const double ny = sumPy / pPair;
                    const double nz = sumPz / pPair;
                    const double pPar = K.px * nx + K.py * ny + K.pz * nz;
                    const double beta = pPair / sumE;
                    const double gamma = sumE / mass;
                    const double pStarPar = gamma * (pPar - beta * K.e);
                    const double eStar = gamma * (K.e - beta * pPar);
                    const double pStar = std::sqrt(
                        std::max(eStar * eStar - kMassK * kMassK, 0.0));
                    if (pStar > 1e-9)
                        cosThetaStar = pStarPar / pStar;
                }

                // Nearest jet to the pair direction.
                // No jet in the event: jetPt = 0, zJet = 0, dRJet = 9.
                double jetPt = 0.0, zJet = 0.0, dRJet = 9.0;
                for (const auto& jet : jets) {
                    const double dr = std::hypot(
                        jet.eta() - pairEta,
                        deltaPhi(jet.phi(), pairPhi));
                    if (dr < dRJet) {
                        dRJet = dr;
                        jetPt = jet.pt();
                    }
                }
                if (jetPt > 0)
                    zJet = pairPt / jetPt;

                // Isolation: charged pT in a cone around the pair axis
                // (excluding the K and the pi), relative to the pair pT
                double isoSum = 0.0;
                for (size_t k = 0; k < n; ++k) {
                    if (k == i || k == j || charge->at(k) == 0)
                        continue;
                    const double dr = std::hypot(
                        eta[k] - pairEta, deltaPhi(phi[k], pairPhi));
                    if (dr < isoR)
                        isoSum += pt[k];
                }
                const double isoRel = (pairPt > 1e-6) ? isoSum / pairPt : 0.0;

                // Fill the row
                b_event        = static_cast<int>(iEvent);
                b_isSignal     = isSignal ? 1 : 0;
                b_weight       = static_cast<float>(weight);
                b_mass         = static_cast<float>(mass);
                b_ptK          = static_cast<float>(K.pt);
                b_ptPi         = static_cast<float>(P.pt);
                b_etaK         = static_cast<float>(K.eta);
                b_etaPi        = static_cast<float>(P.eta);
                b_dEta         = static_cast<float>(dEta);
                b_dPhi         = static_cast<float>(dPhi);
                b_dR           = static_cast<float>(dR);
                b_openAngle    = static_cast<float>(openAngle);
                b_ptBalance    = static_cast<float>(ptBalance);
                b_pairPt       = static_cast<float>(pairPt);
                b_pairEta      = static_cast<float>(pairEta);
                b_cosThetaStar = static_cast<float>(cosThetaStar);
                b_jetPt        = static_cast<float>(jetPt);
                b_zJet         = static_cast<float>(zJet);
                b_dRJet        = static_cast<float>(dRJet);
                b_isoRel       = static_cast<float>(isoRel);

                out.Fill();

                if (isSignal) ++nSignal;
                else          ++nBackground;
            }
        }
    }

    out.Write();
    output.Close();
    input.Close();

    // -----------------------------
    // Console summary
    // -----------------------------
    std::cout << "\n===== buildtree summary =====\n"
              << "Input:              " << inputName << "\n"
              << "Output:             " << outName << "\n"
              << "Events read:        " << nEvents << "\n"
              << "Preselection:       daughter pT > " << minDaughterPt
              << ", |eta| < " << maxAbsEta
              << ", " << massMin << " < mass < " << massMax << "\n"
              << "Smearing:           pT " << smearPt
              << ", eta " << smearEta << ", phi " << smearPhi << "\n"
              << "Signal rows:        " << nSignal << "\n"
              << "Background rows:    " << nBackground
              << "  (kept fraction " << bkgKeep << ")\n";

    if (nSignal > 0)
        std::cout << "Background / signal (unweighted rows): "
                  << double(nBackground) / nSignal << "\n";

    std::cout << "=============================\n";

    return 0;
}
