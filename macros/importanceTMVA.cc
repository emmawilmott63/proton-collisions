#include "tmvaCommon.h"

#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/Reader.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <string>
#include <vector>

// importanceTMVA: permutation importance of the trained BDT's input
// variables, plus how strongly each variable is correlated with the mass.
//
// For each variable (and each feature group) the column is shuffled among
// the events, which destroys its information while keeping its distribution;
// the importance is the drop in the ROC integral (AUC) of the trained BDT.
// Variables in a group are shuffled together (same permutation), so the
// group number measures the information the BDT would lose without that
// kind of feature, even when the variables overlap.
//
// Usage: ./importanceTMVA [candidatesFile] [--option=value ...]
//   candidatesFile   default data/candidates_inclusive.root
//
// Options (defaults in brackets):
//   --dataset=NAME   trained dataset directory               [dataset]
//   --method=NAME    trained method                          [BDT]
//   --events=test|train   which half of the events            [test]
//   --max-bkg=N      keep at most ~N background rows (all signal rows are
//                    kept; kept background rows are re-weighted)   [200000]
//   --repeats=N      shuffles per variable (gives the error)  [3]
//   --seed=N         random seed                              [12345]
//   --presel=EXPR    row selection used in training
//   --vars / --drop / --drop-group   MUST match trainTMVA
//   --tag=NAME / --out=PATH   output [data/importance<tag>.csv]
//
// Weights: flat "weight" only (the same weights as the default training).
//
// Reading the output:
//   auc_drop        baseline AUC minus AUC with this variable shuffled
//                   (larger = the BDT relies on it more)
//   err             spread of auc_drop over the repeats
//   corr_mass_bkg   Pearson correlation of the variable with the K-pi mass
//                   in background rows. A variable that is both important
//                   and strongly mass-correlated can sculpt the background
//                   mass shape (relevant for sideband fits).
// Caveat: with correlated variables, single-variable numbers underestimate
// importance (the BDT can use a correlated partner), so also read the
// group numbers, or retrain without a group (scripts/run_ablation.sh).

namespace {

struct Scored { float score; float w; bool sig; };

// P(signal score > background score) + 0.5 P(equal), weighted
double weightedAuc(std::vector<Scored>& v) {
    std::sort(v.begin(), v.end(),
              [](const Scored& a, const Scored& b) { return a.score < b.score; });
    double sTot = 0, bTot = 0;
    for (const Scored& x : v)
        (x.sig ? sTot : bTot) += x.w;
    if (sTot <= 0 || bTot <= 0)
        return 0.5;
    double bBelow = 0, area = 0;
    size_t i = 0;
    while (i < v.size()) {
        size_t j = i;
        double s = 0, b = 0;
        while (j < v.size() && v[j].score == v[i].score) {
            (v[j].sig ? s : b) += v[j].w;
            ++j;
        }
        area += s * (bBelow + 0.5 * b);
        bBelow += b;
        i = j;
    }
    return area / (sTot * bTot);
}

struct Result {
    std::string name, type, group;
    double drop = 0, err = 0, corr = 0;
    bool hasCorr = false;
};

}  // namespace

int main(int argc, char* argv[]) {

    if (!checkOptions(argc, argv,
            {"dataset", "method", "events", "max-bkg", "repeats", "seed",
             "presel", "vars", "drop", "drop-group", "tag", "out"}))
        return 1;

    const std::string inputName =
        getInputFile(argc, argv, "data/candidates_inclusive.root");
    const std::string tag     = getTag(argc, argv, inputName);
    const std::string dataset = getOption(argc, argv, "dataset", "dataset");
    const std::string method  = getOption(argc, argv, "method", "BDT");
    const std::string events  = getOption(argc, argv, "events", "test");
    const std::string presel  = getOption(argc, argv, "presel", "");
    const std::string outName =
        getOption(argc, argv, "out", "data/importance" + tag + ".csv");
    const double maxBkg = getNumber(argc, argv, "max-bkg", 200000);
    const int repeats =
        std::max(1, static_cast<int>(getNumber(argc, argv, "repeats", 3)));
    const unsigned seed =
        static_cast<unsigned>(getNumber(argc, argv, "seed", 12345));

    if (events != "test" && events != "train") {
        std::cerr << "--events must be test or train\n";
        return 1;
    }
    const std::string parity =
        (events == "test") ? "(event%2)==1" : "(event%2)==0";

    const std::vector<std::string> vars = getTMVAVars(argc, argv);
    const size_t nVars = vars.size();

    const std::string weightFile =
        dataset + "/weights/TMVAClassification_" + method + ".weights.xml";
    if (gSystem->AccessPathName(weightFile.c_str())) {
        std::cerr << "Cannot find " << weightFile
                  << " (run trainTMVA first, same --dataset/--method)\n";
        return 1;
    }

    // -----------------------------
    // Load rows (all signal, a subsample of the background)
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
        gROOT->cd();
        tree = tree->CopyTree(presel.c_str());
    }

    std::vector<float> vals(nVars, 0.f);
    float b_mass = 0, b_weight = 1;
    int b_event = 0, b_isSignal = 0;
    for (size_t i = 0; i < nVars; ++i)
        if (tree->SetBranchAddress(vars[i].c_str(), &vals[i]) < 0) {
            std::cerr << "Branch '" << vars[i] << "' not found\n";
            return 1;
        }
    tree->SetBranchAddress("mass", &b_mass);
    tree->SetBranchAddress("weight", &b_weight);
    tree->SetBranchAddress("event", &b_event);
    if (tree->SetBranchAddress("isSignal", &b_isSignal) < 0) {
        std::cerr << "Branch 'isSignal' not found\n";
        return 1;
    }

    const Long64_t nBkgAll =
        tree->GetEntries(("isSignal==0&&" + parity).c_str());
    const double keepProb =
        nBkgAll > 0 ? std::min(1.0, maxBkg / static_cast<double>(nBkgAll)) : 1.0;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    std::vector<float> X;               // row-major, nRows x nVars
    std::vector<float> mass, w;
    std::vector<char> isSig;

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (((b_event % 2) != 0) != (events == "test"))
            continue;
        const bool sig = (b_isSignal == 1);
        if (!sig && uni(rng) > keepProb)
            continue;
        X.insert(X.end(), vals.begin(), vals.end());
        mass.push_back(b_mass);
        w.push_back(sig ? b_weight
                        : static_cast<float>(b_weight / keepProb));
        isSig.push_back(sig ? 1 : 0);
    }
    const size_t nRows = mass.size();
    size_t nSigRows = 0;
    for (char s : isSig) nSigRows += (s != 0);
    std::cout << "Rows: " << nRows << " (" << nSigRows << " signal, "
              << nRows - nSigRows << " background, background keep fraction "
              << std::setprecision(3) << keepProb << ")\n";
    if (nSigRows < 100 || nSigRows == nRows) {
        std::cerr << "Not enough signal or background rows\n";
        return 1;
    }

    // -----------------------------
    // Trained BDT
    // -----------------------------
    float massF = 0;
    TMVA::Reader reader("!Color:!Silent");
    for (size_t i = 0; i < nVars; ++i)
        reader.AddVariable(vars[i].c_str(), &vals[i]);
    reader.AddSpectator("mass", &massF);
    reader.BookMVA(method.c_str(), weightFile.c_str());

    // AUC with the columns in varIdx replaced by their values from the row
    // given by perm[r] (an empty varIdx is the unshuffled baseline)
    auto evaluateAuc = [&](const std::vector<int>& varIdx,
                           const std::vector<int>& perm) {
        std::vector<Scored> sc(nRows);
        for (size_t r = 0; r < nRows; ++r) {
            for (size_t k = 0; k < nVars; ++k)
                vals[k] = X[r * nVars + k];
            for (int idx : varIdx)
                vals[idx] = X[static_cast<size_t>(perm[r]) * nVars + idx];
            massF = mass[r];
            sc[r] = {static_cast<float>(reader.EvaluateMVA(method.c_str())),
                     w[r], isSig[r] != 0};
        }
        return weightedAuc(sc);
    };

    std::vector<int> identity(nRows);
    std::iota(identity.begin(), identity.end(), 0);

    const double baseline = evaluateAuc({}, identity);
    std::cout << "Baseline AUC (" << events << " events, flat weights): "
              << std::setprecision(5) << baseline << "\n" << std::flush;

    auto permutationImportance = [&](const std::vector<int>& varIdx,
                                     double& meanDrop, double& err) {
        std::vector<double> drops;
        for (int rep = 0; rep < repeats; ++rep) {
            std::vector<int> perm = identity;
            std::shuffle(perm.begin(), perm.end(), rng);
            drops.push_back(baseline - evaluateAuc(varIdx, perm));
        }
        meanDrop = std::accumulate(drops.begin(), drops.end(), 0.0) /
                   static_cast<double>(drops.size());
        double var = 0;
        for (double d : drops)
            var += (d - meanDrop) * (d - meanDrop);
        err = drops.size() > 1
                  ? std::sqrt(var / static_cast<double>(drops.size() - 1))
                  : 0.0;
    };

    std::vector<Result> singles, groups;

    for (size_t k = 0; k < nVars; ++k) {
        Result r;
        r.name = vars[k];
        r.type = "variable";
        r.group = featureGroup(vars[k]);
        permutationImportance({static_cast<int>(k)}, r.drop, r.err);

        // correlation with the mass in background rows
        double n = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
        for (size_t i = 0; i < nRows; ++i) {
            if (isSig[i])
                continue;
            const double x = X[i * nVars + k], y = mass[i];
            n += 1; sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
        }
        const double vx = sxx - sx * sx / n, vy = syy - sy * sy / n;
        if (n > 2 && vx > 0 && vy > 0) {
            r.corr = (sxy - sx * sy / n) / std::sqrt(vx * vy);
            r.hasCorr = true;
        }
        singles.push_back(r);
        std::cout << "  " << std::setw(14) << std::left << r.name
                  << std::right << std::fixed << std::setprecision(5)
                  << r.drop << "\n" << std::flush;
    }

    std::map<std::string, std::vector<int>> groupIdx;
    for (size_t k = 0; k < nVars; ++k)
        groupIdx[featureGroup(vars[k])].push_back(static_cast<int>(k));
    for (const auto& g : groupIdx) {
        Result r;
        r.name = g.first;
        r.type = "group";
        r.group = g.first;
        permutationImportance(g.second, r.drop, r.err);
        groups.push_back(r);
        std::cout << "  group " << std::setw(10) << std::left << r.name
                  << std::right << r.drop << "\n" << std::flush;
    }

    const auto byDrop = [](const Result& a, const Result& b) {
        return a.drop > b.drop;
    };
    std::sort(singles.begin(), singles.end(), byDrop);
    std::sort(groups.begin(), groups.end(), byDrop);

    // -----------------------------
    // Report
    // -----------------------------
    std::cout << std::fixed << std::setprecision(5)
              << "\n===== permutation importance (AUC drop) =====\n"
              << "Baseline AUC " << baseline << ", " << repeats
              << " shuffle(s) each\n\nGroups:\n";
    for (const Result& r : groups)
        std::cout << "  " << std::setw(12) << std::left << r.name << std::right
                  << std::setw(10) << r.drop << " +- " << r.err << "\n";
    std::cout << "\nVariables (corr = correlation with mass in background):\n";
    for (const Result& r : singles)
        std::cout << "  " << std::setw(14) << std::left << r.name << std::right
                  << std::setw(10) << r.drop << " +- " << std::setw(8) << r.err
                  << "   corr " << std::setprecision(3) << std::setw(7)
                  << r.corr << std::setprecision(5) << "   [" << r.group
                  << "]\n";

    gSystem->mkdir("data", true);
    std::ofstream csv(outName);
    csv << "name,type,group,auc_drop,err,corr_mass_bkg,baseline_auc,events\n";
    for (const std::vector<Result>* list : {&groups, &singles})
        for (const Result& r : *list)
            csv << r.name << "," << r.type << "," << r.group << ","
                << r.drop << "," << r.err << ","
                << (r.hasCorr ? std::to_string(r.corr) : std::string()) << ","
                << baseline << "," << events << "\n";
    std::cout << "\nWritten: " << outName << "\n";

    input.Close();
    return 0;
}
