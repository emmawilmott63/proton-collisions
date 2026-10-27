#include "tmvaCommon.h"

#include "TCut.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Types.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// trainTMVA: train a boosted decision tree (TMVA) to separate signal
// (true D0 -> K pi) from background K pi pairs.
//
// Usage: ./trainTMVA [candidatesFile] [--option=value ...]
//   candidatesFile    default data/candidates_inclusive.root
//
// ---- Input: one file or two ------------------------------------------------
//   One file (default): signal = rows with isSignal = 1, background = rows
//   with isSignal = 0, both from candidatesFile (simulation).
//   Two files (data-driven background, e.g. STAR):
//     --signal-file=PATH       signal rows (isSignal = 1), e.g. embedding
//     --background-file=PATH   background rows, e.g. real-data mass sidebands
//   Both files must contain a TTree "Candidates" with the same feature
//   branches plus event, mass and weight (see docs/TMVA_REFERENCE.md).
//   If the background file has no isSignal branch, every row is a candidate.
//
// ---- Row selection (hard conditions) ---------------------------------------
//   --presel=EXPR       cut on BOTH classes, e.g. "pairPt>1.0"  (--cut is the
//                       old name). Use the same --presel in applyTMVA.
//   --sig-cut=EXPR      extra cut on the signal class only
//   --bkg-cut=EXPR      extra cut on the background class only, e.g. a mass
//                       sideband: "mass<1.80||mass>1.93"
//   (EXPR is a ROOT TTree expression; && and || are allowed.)
//
// ---- Variables -------------------------------------------------------------
//   --vars=a,b,c        use exactly these branches   (see tmvaCommon.h)
//   --drop=a,b          remove branches from the default list
//   --drop-group=g,h    remove feature groups (daughters, pair, jet,
//                       isolation, event)
//
// ---- BDT parameters (defaults in brackets) ---------------------------------
//   --method=NAME       TMVA method title; names the weights file. BDTG
//                       defaults to gradient boosting, else AdaBoost  [BDT]
//   --boost-type=T      AdaBoost | RealAdaBoost | Grad | Bagging       [AdaBoost]
//   --ntrees=N          number of trees                                [400]
//   --depth=N           max tree depth                                 [3]
//   --min-node=X        min leaf size, percent of the training events  [2.0]
//   --ncuts=N           grid points per variable when scanning splits  [20]
//   --separation=S      GiniIndex | CrossEntropy | MisClassificationError
//                       | SDivSqrtSPlusB   (not used by Grad)          [GiniIndex]
//   --ada-beta=X        AdaBoost learning rate                         [0.5]
//   --shrinkage=X       Grad learning rate                             [0.1]
//   --bagging=0|1       train each tree on a random subsample          [1]
//   --bag-fraction=X    subsample fraction                             [0.5]
//   --prune-method=M    NoPruning | ExpectedError | CostComplexity     [NoPruning]
//   --prune-strength=X  pruning strength (ignored for NoPruning)       [0]
//   --bdt-options=STR   raw TMVA options appended last (override anything),
//                       e.g. "DoBoostMonitor:SkipNormalization"
//
// ---- Weights and other -----------------------------------------------------
//   --train-weights=flat|xs
//                       flat: weight only. xs: weight * xsWeight (physical
//                       mixture of pThat bins; a few heavy rows can
//                       dominate)                                      [flat]
//   --norm-mode=M       EqualNumEvents | NumEvents | None: how class
//                       weights are balanced in training   [EqualNumEvents]
//   --dataset=NAME      model name: TMVA-output/models/<NAME>/         [dataset]
//   --out=PATH          TMVA output file [TMVA-output/training/tmva_training<tag>.root]
//
// Train / test split: candidates with EVEN event numbers train, ODD event
// numbers test. Splitting by event (not by candidate) keeps all candidates
// of one event on the same side, so event-level features cannot leak
// between the two samples. "event" only needs to be an integer that is
// unique per event within a file.
//
// "mass" is only a spectator: it is NOT a training variable.
//
// Output (all under TMVA-output/):
//   models/<dataset>/weights/TMVAClassification_<method>.weights.xml
//                          (read by applyTMVA)
//   models/<dataset>/train_config.txt  (everything needed to reproduce the
//                          training)
//   training/tmva_training<tag>.root   (open with TMVAGui)
// The variable ranking is printed in the TMVA log ("Ranking result").
// Training again with the same --dataset overwrites that model.

namespace {

// "(a)&&(b)", skipping empty parts
std::string andCuts(const std::vector<std::string>& parts) {
    std::string out;
    for (const std::string& p : parts) {
        if (p.empty())
            continue;
        if (!out.empty())
            out += "&&";
        out += "(" + p + ")";
    }
    return out;
}

std::string num(double x) {
    std::ostringstream s;
    s << x;
    return s.str();
}

struct Source {
    TFile* file = nullptr;
    TTree* train = nullptr;   // even event numbers (in memory)
    TTree* test  = nullptr;   // odd event numbers  (in memory)
    bool hasIsSignal = false;
};

bool loadSource(const std::string& path, const std::vector<std::string>& vars,
                bool needXs, Source& src) {
    src.file = TFile::Open(path.c_str(), "READ");
    if (!src.file || src.file->IsZombie()) {
        std::cerr << "Could not open " << path << "\n";
        return false;
    }
    TTree* cand = nullptr;
    src.file->GetObject("Candidates", cand);
    if (!cand) {
        std::cerr << "No TTree 'Candidates' in " << path << "\n";
        return false;
    }

    std::vector<std::string> required = vars;
    required.push_back("event");
    required.push_back("mass");
    required.push_back("weight");
    if (needXs)
        required.push_back("xsWeight");
    for (const std::string& b : required)
        if (!cand->GetBranch(b.c_str())) {
            std::cerr << "Branch '" << b << "' not found in Candidates of "
                      << path << "\n";
            return false;
        }
    src.hasIsSignal = cand->GetBranch("isSignal") != nullptr;

    gROOT->cd();   // copies go to memory, not into the (read-only) input file
    src.train = cand->CopyTree("(event%2)==0");
    src.test  = cand->CopyTree("(event%2)==1");
    return true;
}

}  // namespace

int main(int argc, char* argv[]) {

    if (!checkOptions(argc, argv,
            {"method", "ntrees", "depth", "min-node", "ncuts", "separation",
             "boost-type", "ada-beta", "shrinkage", "bagging", "bag-fraction",
             "prune-method", "prune-strength", "bdt-options", "norm-mode",
             "train-weights", "presel", "cut", "sig-cut", "bkg-cut",
             "signal-file", "background-file", "vars", "drop", "drop-group",
             "dataset", "out", "tag"}))
        return 1;

    const std::string inputName =
        getInputFile(argc, argv, "data/candidates_inclusive.root");
    const std::string tag     = getTag(argc, argv, inputName);
    const std::string dataset = getOption(argc, argv, "dataset", "dataset");
    const std::string method  = getOption(argc, argv, "method", "BDT");
    const std::string outName =
        getOption(argc, argv, "out",
                  outDir("training") + "/tmva_training" + tag + ".root");

    if (!validDatasetName(dataset)) {
        std::cerr << "--dataset may contain only letters, digits, '_', '-' "
                     "and '.', and must not start with '.'\n";
        return 1;
    }

    const std::string signalFile = getOption(argc, argv, "signal-file", "");
    const std::string backgroundFile =
        getOption(argc, argv, "background-file", "");
    const bool twoFiles = !signalFile.empty() || !backgroundFile.empty();
    if (twoFiles && (signalFile.empty() || backgroundFile.empty())) {
        std::cerr << "--signal-file and --background-file must be given "
                     "together\n";
        return 1;
    }

    const std::string presel =
        getOption(argc, argv, "presel", getOption(argc, argv, "cut", ""));
    const std::string sigExtra = getOption(argc, argv, "sig-cut", "");
    const std::string bkgExtra = getOption(argc, argv, "bkg-cut", "");

    const std::string trainWeights =
        getOption(argc, argv, "train-weights", "flat");
    if (trainWeights != "flat" && trainWeights != "xs") {
        std::cerr << "--train-weights must be flat or xs\n";
        return 1;
    }
    const std::string weightExpr =
        (trainWeights == "xs") ? "weight*xsWeight" : "weight";
    const std::string normMode =
        getOption(argc, argv, "norm-mode", "EqualNumEvents");

    // ---- BDT parameters ----
    const int    nTrees   = static_cast<int>(getNumber(argc, argv, "ntrees", 400));
    const int    depth    = static_cast<int>(getNumber(argc, argv, "depth", 3));
    const double minNode  = getNumber(argc, argv, "min-node", 2.0);
    const int    nCuts    = static_cast<int>(getNumber(argc, argv, "ncuts", 20));
    const double adaBeta  = getNumber(argc, argv, "ada-beta", 0.5);
    const double shrink   = getNumber(argc, argv, "shrinkage", 0.1);
    const bool   bagging  = getNumber(argc, argv, "bagging", 1) != 0;
    const double bagFrac  = getNumber(argc, argv, "bag-fraction", 0.5);
    const std::string separation =
        getOption(argc, argv, "separation", "GiniIndex");
    const std::string pruneMethod =
        getOption(argc, argv, "prune-method", "NoPruning");
    const double pruneStrength = getNumber(argc, argv, "prune-strength", 0.0);
    const std::string rawOpts = getOption(argc, argv, "bdt-options", "");
    const std::string boostType = getOption(argc, argv, "boost-type",
        method == "BDTG" ? "Grad" : "AdaBoost");

    if (boostType != "AdaBoost" && boostType != "RealAdaBoost" &&
        boostType != "Grad" && boostType != "Bagging") {
        std::cerr << "--boost-type must be AdaBoost, RealAdaBoost, Grad or "
                     "Bagging\n";
        return 1;
    }

    std::string opts = "!H:!V:NTrees=" + std::to_string(nTrees) +
                       ":MaxDepth=" + std::to_string(depth) +
                       ":MinNodeSize=" + num(minNode) + "%" +
                       ":nCuts=" + std::to_string(nCuts) +
                       ":BoostType=" + boostType +
                       ":SeparationType=" + separation +
                       ":PruneMethod=" + pruneMethod;
    if (boostType == "AdaBoost")
        opts += ":AdaBoostBeta=" + num(adaBeta);
    if (boostType == "Grad")
        opts += ":Shrinkage=" + num(shrink);
    if (bagging)
        opts += ":UseBaggedBoost:BaggedSampleFraction=" + num(bagFrac);
    if (pruneMethod != "NoPruning")
        opts += ":PruneStrength=" + num(pruneStrength);
    if (!rawOpts.empty())
        opts += ":" + rawOpts;

    const std::vector<std::string> vars = getTMVAVars(argc, argv);
    if (vars.empty()) {
        std::cerr << "No input variables selected\n";
        return 1;
    }

    gSystem->mkdir(outDir("training").c_str(), true);
    gSystem->mkdir((modelDir(dataset) + "/weights").c_str(), true);

    // -----------------------------
    // Read candidates, split by event parity (trees live in memory)
    // -----------------------------
    const bool needXs = (trainWeights == "xs");
    const std::string sigPath = twoFiles ? signalFile : inputName;
    const std::string bkgPath = twoFiles ? backgroundFile : inputName;

    Source sigSrc, bkgSrc;
    if (!loadSource(sigPath, vars, needXs, sigSrc))
        return 1;
    if (twoFiles) {
        if (!loadSource(bkgPath, vars, needXs, bkgSrc))
            return 1;
    } else {
        bkgSrc = sigSrc;
    }
    if (!sigSrc.hasIsSignal) {
        std::cerr << "The signal file needs an isSignal branch\n";
        return 1;
    }

    const std::string sigSel =
        andCuts({"isSignal>0.5", presel, sigExtra});
    const std::string bkgSel =
        andCuts({bkgSrc.hasIsSignal ? "isSignal<0.5" : "", presel, bkgExtra});

    const Long64_t nSigTrain = sigSrc.train->GetEntries(sigSel.c_str());
    const Long64_t nSigTest  = sigSrc.test->GetEntries(sigSel.c_str());
    const Long64_t nBkgTrain = bkgSrc.train->GetEntries(bkgSel.c_str());
    const Long64_t nBkgTest  = bkgSrc.test->GetEntries(bkgSel.c_str());

    std::cout << "Training rows: signal " << nSigTrain
              << ", background " << nBkgTrain << "\n"
              << "Testing  rows: signal " << nSigTest
              << ", background " << nBkgTest << "\n";

    if (nSigTrain < 100 || nSigTest < 100 || nBkgTrain < 100 || nBkgTest < 100) {
        std::cerr << "Too few signal or background rows to train. Generate "
                     "more events (or loosen the selection).\n";
        return 1;
    }

    // -----------------------------
    // Record the configuration next to the weights
    // -----------------------------
    {
        std::ofstream cfg(modelDir(dataset) + "/train_config.txt");
        cfg << "command:";
        for (int i = 0; i < argc; ++i)
            cfg << " " << argv[i];
        cfg << "\nsignal source:     " << sigPath
            << "\nbackground source: " << bkgPath
            << "\nsignal selection:     " << sigSel
            << "\nbackground selection: " << bkgSel
            << "\nweights:           " << weightExpr
            << "\nnorm mode:         " << normMode
            << "\nmethod title:      " << method
            << "\nTMVA options:      " << opts
            << "\nrows (train/test): signal " << nSigTrain << "/" << nSigTest
            << ", background " << nBkgTrain << "/" << nBkgTest
            << "\nvariables (" << vars.size() << "):";
        for (const std::string& v : vars)
            cfg << " " << v;
        cfg << "\n";
    }

    // -----------------------------
    // TMVA setup
    // -----------------------------
    TMVA::Tools::Instance();

    TFile* outFile = TFile::Open(outName.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Could not create " << outName << "\n";
        return 1;
    }

    TMVA::Factory factory("TMVAClassification", outFile,
        "!V:!Silent:!Color:!DrawProgressBar:AnalysisType=Classification");

    // TMVA writes <dataset>/weights relative to the working directory, so
    // run it from TMVA-output/models/: the model then appears in place and
    // nothing has to be moved afterwards. The output and input files are
    // already open and the trees are in memory, so they are not affected.
    const std::string startDir = gSystem->WorkingDirectory();
    if (!gSystem->ChangeDirectory(outDir("models").c_str())) {
        std::cerr << "Could not enter " << outDir("models") << "\n";
        return 1;
    }

    TMVA::DataLoader* loader = new TMVA::DataLoader(dataset);

    for (const std::string& v : vars)
        loader->AddVariable(v.c_str(), 'F');
    loader->AddSpectator("mass", "mass");

    loader->AddSignalTree(sigSrc.train, 1.0, TMVA::Types::kTraining);
    loader->AddSignalTree(sigSrc.test,  1.0, TMVA::Types::kTesting);
    loader->AddBackgroundTree(bkgSrc.train, 1.0, TMVA::Types::kTraining);
    loader->AddBackgroundTree(bkgSrc.test,  1.0, TMVA::Types::kTesting);

    loader->SetSignalWeightExpression(weightExpr.c_str());
    loader->SetBackgroundWeightExpression(weightExpr.c_str());

    // NormMode: during training, signal and background get equal total
    // weight (the sample is ~1000:1 background). This does not change the
    // test-sample weights that applyTMVA uses for yields.
    loader->PrepareTrainingAndTestTree(TCut(sigSel.c_str()),
                                       TCut(bkgSel.c_str()),
                                       ("NormMode=" + normMode + ":!V").c_str());

    // BDTG is also TMVA's kBDT type; only the boost options differ.
    factory.BookMethod(loader, TMVA::Types::kBDT, method.c_str(), opts.c_str());

    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();

    gSystem->ChangeDirectory(startDir.c_str());

    outFile->Close();
    sigSrc.file->Close();
    if (twoFiles)
        bkgSrc.file->Close();

    std::cout << "\n===== trainTMVA summary =====\n"
              << "Signal:   " << sigPath << "  (" << sigSel << ")\n"
              << "Background: " << bkgPath << "  (" << bkgSel << ")\n"
              << "Method:   " << method << "\n"
              << "Options:  " << opts << "\n"
              << "Training weights: " << weightExpr
              << ", norm mode " << normMode << "\n"
              << "Variables (" << vars.size() << "):";
    for (const std::string& v : vars) std::cout << " " << v;
    std::cout << "\nWeights:  " << weightFilePath(dataset, method) << "\n"
              << "Config:   " << modelDir(dataset) << "/train_config.txt\n"
              << "TMVA file: " << outName
              << "   (root -l -e 'TMVA::TMVAGui(\"" << outName << "\")')\n"
              << "Next:     ./applyTMVA " << inputName
              << " --dataset=" << dataset << " --method=" << method << "\n"
              << "=============================\n";
    return 0;
}
