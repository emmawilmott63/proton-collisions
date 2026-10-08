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
//
//   daughter / pair     : ptK, ptPi, ptHard, ptSoft, ptSum, ptBalance,
//                         etaK, etaPi, dEta, dPhi, dR, openAngle,
//                         pairPt, pairEta, cosThetaStar
//
//   jet association     : sameJet    1 if K and pi are constituents of the
//                                    same jet, else 0
//                         jetPt, jetEta, jetMass, jetNConst
//                                    properties of the REFERENCE jet: the
//                                    shared jet if sameJet, otherwise the
//                                    jet nearest the pair direction
//                         zJet       pairPt / jetPt
//                         dRJet      distance from the pair direction to
//                                    the reference-jet axis
//                         (no jet in the event: jetPt = jetEta = jetMass =
//                          jetNConst = zJet = 0, dRJet = 9, sameJet = 0)
//
//   isolation           : isoRel02, isoRel03, isoRel, isoRel06
//                                    charged pT in a cone of R = 0.2, 0.3,
//                                    0.4 (isoRel), 0.6 around the pair axis,
//                                    / pair pT, K and pi excluded
//                         isoNeutralRel
//                                    same for neutral (non-neutrino)
//                                    particles, R = 0.4
//
//   event level         : nCharged, nJets, sphericityT, thrustT
//                                    computed from charged particles (and
//                                    jets) with |eta| < 1, pT > 0.2 GeV/c,
//                                    from the UNSMEARED particles
//
// Every feature branch is written as float (including the counts and the
// sameJet flag) so downstream code can read them all the same way.

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

bool isNeutrino(int absId) {
    return absId == 12 || absId == 14 || absId == 16;
}

// Transverse event shapes from a list of (px, py).
//   sphericityT : 2 * lambda_min / (lambda_min + lambda_max) of the
//                 2x2 momentum tensor.  0 = pencil-like, 1 = isotropic.
//   thrustT     : max over axes n of sum|p.n| / sum pT.
//                 1 = pencil-like, 2/pi = isotropic.
// With fewer than two particles the event counts as pencil-like
// (sphericityT = 0, thrustT = 1).
void transverseShapes(const std::vector<double>& vpx,
                      const std::vector<double>& vpy,
                      double& sphericityT, double& thrustT) {

    sphericityT = 0.0;
    thrustT     = 1.0;

    const size_t m = vpx.size();
    if (m < 2)
        return;

    double sxx = 0, syy = 0, sxy = 0, sumPt = 0;
    for (size_t k = 0; k < m; ++k) {
        sxx += vpx[k] * vpx[k];
        syy += vpy[k] * vpy[k];
        sxy += vpx[k] * vpy[k];
        sumPt += std::hypot(vpx[k], vpy[k]);
    }
    if (sumPt <= 0)
        return;

    // Eigenvalues of [[sxx, sxy], [sxy, syy]]
    const double tr   = sxx + syy;
    const double disc = std::sqrt(std::max(
        0.25 * (sxx - syy) * (sxx - syy) + sxy * sxy, 0.0));
    const double lamMin = 0.5 * tr - disc;
    if (tr > 0)
        sphericityT = 2.0 * lamMin / tr;

    // Scan the axis angle over [0, pi)
    const int nAngles = 180;
    double best = 0.0;
    for (int a = 0; a < nAngles; ++a) {
        const double ang = kPi * a / nAngles;
        const double nx = std::cos(ang);
        const double ny = std::sin(ang);
        double s = 0.0;
        for (size_t k = 0; k < m; ++k)
            s += std::abs(vpx[k] * nx + vpy[k] * ny);
        best = std::max(best, s);
    }
    thrustT = best / sumPt;
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

    // Isolation cones (isoCones[2] is the original 0.4 isolation, and the
    // radius also used for the neutral sum)
    const int    nCones          = 4;
    const double isoCones[nCones] = {0.2, 0.3, 0.4, 0.6};
    const int    kNeutralCone    = 2;

    // Acceptance for the event-level quantities
    const double evMaxAbsEta = 1.0;
    const double evMinPt     = 0.2;   // GeV/c

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

    // daughters / pair
    float b_ptK = 0, b_ptPi = 0, b_ptHard = 0, b_ptSoft = 0, b_ptSum = 0;
    float b_etaK = 0, b_etaPi = 0;
    float b_dEta = 0, b_dPhi = 0, b_dR = 0;
    float b_openAngle = 0, b_ptBalance = 0;
    float b_pairPt = 0, b_pairEta = 0, b_cosThetaStar = 0;

    // jets
    float b_sameJet = 0;
    float b_jetPt = 0, b_jetEta = 0, b_jetMass = 0, b_jetNConst = 0;
    float b_zJet = 0, b_dRJet = 0;

    // isolation
    float b_isoRel02 = 0, b_isoRel03 = 0, b_isoRel = 0, b_isoRel06 = 0;
    float b_isoNeutralRel = 0;

    // event level
    float b_nCharged = 0, b_nJets = 0, b_sphericityT = 0, b_thrustT = 0;

    out.Branch("event", &b_event);
    out.Branch("isSignal", &b_isSignal);
    out.Branch("weight", &b_weight);
    out.Branch("mass", &b_mass);

    out.Branch("ptK", &b_ptK);
    out.Branch("ptPi", &b_ptPi);
    out.Branch("ptHard", &b_ptHard);
    out.Branch("ptSoft", &b_ptSoft);
    out.Branch("ptSum", &b_ptSum);
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

    out.Branch("sameJet", &b_sameJet);
    out.Branch("jetPt", &b_jetPt);
    out.Branch("jetEta", &b_jetEta);
    out.Branch("jetMass", &b_jetMass);
    out.Branch("jetNConst", &b_jetNConst);
    out.Branch("zJet", &b_zJet);
    out.Branch("dRJet", &b_dRJet);

    out.Branch("isoRel02", &b_isoRel02);
    out.Branch("isoRel03", &b_isoRel03);
    out.Branch("isoRel", &b_isoRel);
    out.Branch("isoRel06", &b_isoRel06);
    out.Branch("isoNeutralRel", &b_isoNeutralRel);

    out.Branch("nCharged", &b_nCharged);
    out.Branch("nJets", &b_nJets);
    out.Branch("sphericityT", &b_sphericityT);
    out.Branch("thrustT", &b_thrustT);

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);

    long long nSignal = 0;
    long long nBackground = 0;
    long long nEvents = 0;
    long long nSameJetSig = 0, nSameJetBkg = 0;

    // -----------------------------
    // Event loop
    // -----------------------------
    const Long64_t nEntries = tree->GetEntries();

    for (Long64_t iEvent = 0; iEvent < nEntries; ++iEvent) {

        tree->GetEntry(iEvent);
        ++nEvents;

        const size_t n = pid->size();

        // -------------------------------------------------
        // Jets from the (unsmeared) final-state particles.
        // user_index remembers which particle each input is, so
        // we can tell which jet a given K or pi ended up in.
        // -------------------------------------------------
        std::vector<fastjet::PseudoJet> jetInputs;
        for (size_t k = 0; k < n; ++k) {
            if (isNeutrino(std::abs(pid->at(k))))
                continue;
            jetInputs.emplace_back(px->at(k), py->at(k),
                                   pz->at(k), energy->at(k));
            jetInputs.back().set_user_index(static_cast<int>(k));
        }

        fastjet::ClusterSequence cs(jetInputs, jetDef);
        std::vector<fastjet::PseudoJet> jets =
            fastjet::sorted_by_pt(cs.inclusive_jets(jetPtMin));

        // jetOf[k] = index into `jets` of the jet containing particle k,
        // or -1 if it is in no jet above jetPtMin
        std::vector<int> jetOf(n, -1);
        for (size_t jj = 0; jj < jets.size(); ++jj)
            for (const auto& c : jets[jj].constituents())
                jetOf[c.user_index()] = static_cast<int>(jj);

        // Transverse momentum and direction of every particle
        // (unsmeared; used for the isolation sums and event shapes)
        std::vector<double> pt(n), eta(n), phi(n);
        for (size_t k = 0; k < n; ++k) {
            pt[k]  = std::hypot(px->at(k), py->at(k));
            eta[k] = (pt[k] > 1e-6) ? std::asinh(pz->at(k) / pt[k]) : 0.0;
            phi[k] = std::atan2(py->at(k), px->at(k));
        }

        // -------------------------------------------------
        // Event-level quantities (once per event)
        // -------------------------------------------------
        std::vector<double> chPx, chPy;
        for (size_t k = 0; k < n; ++k) {
            if (charge->at(k) == 0)
                continue;
            if (pt[k] < evMinPt || std::abs(eta[k]) > evMaxAbsEta)
                continue;
            chPx.push_back(px->at(k));
            chPy.push_back(py->at(k));
        }

        int nJetsEv = 0;
        for (const auto& jet : jets)
            if (std::abs(jet.eta()) < evMaxAbsEta)
                ++nJetsEv;

        double sphericityT = 0.0, thrustT = 1.0;
        transverseShapes(chPx, chPy, sphericityT, thrustT);

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

                // Harder / softer daughter and scalar pT sum
                const double ptHard = std::max(K.pt, P.pt);
                const double ptSoft = std::min(K.pt, P.pt);
                const double ptSum  = K.pt + P.pt;

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

                // Same-jet flag. The jets use the unsmeared particles, so
                // this does not change with the smearing options.
                const bool sameJet =
                    jetOf[i] >= 0 && jetOf[i] == jetOf[j];

                // Reference jet: the shared jet if there is one, otherwise
                // the jet nearest the pair direction.
                int refJet = -1;
                double dRJet = 9.0;

                if (sameJet) {
                    refJet = jetOf[i];
                    dRJet = std::hypot(
                        jets[refJet].eta() - pairEta,
                        deltaPhi(jets[refJet].phi(), pairPhi));
                } else {
                    for (size_t jj = 0; jj < jets.size(); ++jj) {
                        const double dr = std::hypot(
                            jets[jj].eta() - pairEta,
                            deltaPhi(jets[jj].phi(), pairPhi));
                        if (dr < dRJet) {
                            dRJet = dr;
                            refJet = static_cast<int>(jj);
                        }
                    }
                }

                double jetPt = 0.0, jetEta = 0.0, jetMass = 0.0;
                double jetNConst = 0.0, zJet = 0.0;
                if (refJet >= 0) {
                    jetPt     = jets[refJet].pt();
                    jetEta    = jets[refJet].eta();
                    jetMass   = jets[refJet].m();
                    jetNConst = static_cast<double>(
                        jets[refJet].constituents().size());
                    if (jetPt > 0)
                        zJet = pairPt / jetPt;
                }

                // Isolation: pT in cones around the pair axis (K and pi
                // excluded), relative to the pair pT.
                //   charged: four cone sizes
                //   neutral: one cone, neutrinos left out
                double isoCharged[nCones] = {0.0, 0.0, 0.0, 0.0};
                double isoNeutral = 0.0;

                for (size_t k = 0; k < n; ++k) {
                    if (k == i || k == j)
                        continue;

                    const double dr = std::hypot(
                        eta[k] - pairEta, deltaPhi(phi[k], pairPhi));

                    if (charge->at(k) != 0) {
                        for (int c = 0; c < nCones; ++c)
                            if (dr < isoCones[c])
                                isoCharged[c] += pt[k];
                    } else if (!isNeutrino(std::abs(pid->at(k)))) {
                        if (dr < isoCones[kNeutralCone])
                            isoNeutral += pt[k];
                    }
                }

                const double invPairPt = (pairPt > 1e-6) ? 1.0 / pairPt : 0.0;

                // Fill the row
                b_event        = static_cast<int>(iEvent);
                b_isSignal     = isSignal ? 1 : 0;
                b_weight       = static_cast<float>(weight);
                b_mass         = static_cast<float>(mass);

                b_ptK          = static_cast<float>(K.pt);
                b_ptPi         = static_cast<float>(P.pt);
                b_ptHard       = static_cast<float>(ptHard);
                b_ptSoft       = static_cast<float>(ptSoft);
                b_ptSum        = static_cast<float>(ptSum);
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

                b_sameJet      = sameJet ? 1.0f : 0.0f;
                b_jetPt        = static_cast<float>(jetPt);
                b_jetEta       = static_cast<float>(jetEta);
                b_jetMass      = static_cast<float>(jetMass);
                b_jetNConst    = static_cast<float>(jetNConst);
                b_zJet         = static_cast<float>(zJet);
                b_dRJet        = static_cast<float>(dRJet);

                b_isoRel02     = static_cast<float>(isoCharged[0] * invPairPt);
                b_isoRel03     = static_cast<float>(isoCharged[1] * invPairPt);
                b_isoRel       = static_cast<float>(isoCharged[2] * invPairPt);
                b_isoRel06     = static_cast<float>(isoCharged[3] * invPairPt);
                b_isoNeutralRel = static_cast<float>(isoNeutral * invPairPt);

                b_nCharged     = static_cast<float>(chPx.size());
                b_nJets        = static_cast<float>(nJetsEv);
                b_sphericityT  = static_cast<float>(sphericityT);
                b_thrustT      = static_cast<float>(thrustT);

                out.Fill();

                if (isSignal) {
                    ++nSignal;
                    if (sameJet) ++nSameJetSig;
                } else {
                    ++nBackground;
                    if (sameJet) ++nSameJetBkg;
                }
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

    if (nSignal > 0)
        std::cout << "K and pi in the same jet:  signal "
                  << 100.0 * nSameJetSig / nSignal << " %";
    if (nBackground > 0)
        std::cout << ",  background "
                  << 100.0 * nSameJetBkg / nBackground << " %";
    std::cout << "\n";

    std::cout << "=============================\n";

    return 0;
}
