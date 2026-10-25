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

#include <iostream>
#include <string>

// trainTMVA: train a boosted decision tree (TMVA) to separate true
// D0 -> K pi candidates (isSignal = 1) from combinatorial pairs.
//
// Usage: ./trainTMVA [candidatesFile] [--option=value ...]
//   candidatesFile    default data/candidates_inclusive.root
//
// Options (defaults in brackets):
//   --method=BDT|BDTG   AdaBoost BDT or gradient-boost BDTG  [BDT]
//   --ntrees=N          number of trees                       [400]
//   --depth=N           max tree depth                        [3]
//   --min-node=X        min node size, percent of events      [2.0]
//   --train-weights=flat|xs
//                       flat: train with "weight" only (every generated
//                       event counts the same, so each pThat bin
//                       contributes by its number of events).
//                       xs: also multiply by xsWeight (physical mixture of
//                       pThat bins; a few heavy low-pThat rows can dominate
//                       and make training noisy).               [flat]
//   --cut=EXPR          extra TCut on both classes, e.g. "pairPt>1.0"
//   --vars=a,b,c        use exactly these branches   (see tmvaCommon.h)
//   --drop=a,b          remove branches from the default list
//   --dataset=NAME      TMVA dataset / weight directory       [dataset]
//   --out=PATH          TMVA output file   [data/tmva_training<tag>.root]
//
// Train / test split: candidates from EVEN event numbers train, ODD event
// numbers test. Splitting by event (not by candidate) keeps all candidates
// of one event on the same side, so the event-level features (nCharged,
// sphericityT, ...) cannot leak between the two samples.
//
// Event weights come from the "weight" branch (1 for signal, 1/bkg-keep for
// background). "mass" is only a spectator: it is NOT a training variable.
//
// Output: <dataset>/weights/TMVAClassification_<method>.weights.xml, read
// by applyTMVA, plus the usual TMVA training file (open with TMVAGui).

int main(int argc, char* argv[]) {

    const std::string inputName =
        getInputFile(argc, argv, "data/candidates_inclusive.root");
    const std::string tag     = getTag(argc, argv, inputName);
    const std::string dataset = getOption(argc, argv, "dataset", "dataset");
    const std::string method  = getOption(argc, argv, "method", "BDT");
    const std::string outName =
        getOption(argc, argv, "out", "data/tmva_training" + tag + ".root");
    const std::string extraCut = getOption(argc, argv, "cut", "");
    const std::string trainWeights =
        getOption(argc, argv, "train-weights", "flat");

    if (trainWeights != "flat" && trainWeights != "xs") {
        std::cerr << "--train-weights must be flat or xs\n";
        return 1;
    }
    const std::string weightExpr =
        (trainWeights == "xs") ? "weight*xsWeight" : "weight";

    const int    nTrees  = static_cast<int>(getNumber(argc, argv, "ntrees", 400));
    const int    depth   = static_cast<int>(getNumber(argc, argv, "depth", 3));
    const double minNode = getNumber(argc, argv, "min-node", 2.0);

    if (method != "BDT" && method != "BDTG") {
        std::cerr << "--method must be BDT or BDTG\n";
        return 1;
    }

    const std::vector<std::string> vars = getTMVAVars(argc, argv);
    if (vars.empty()) {
        std::cerr << "No input variables selected\n";
        return 1;
    }

    gSystem->mkdir("data", true);

    // -----------------------------
    // Read candidates, split by event parity (trees live in memory)
    // -----------------------------
    TFile input(inputName.c_str(), "READ");
    if (input.IsZombie()) {
        std::cerr << "Could not open " << inputName << "\n";
        return 1;
    }
    TTree* cand = nullptr;
    input.GetObject("Candidates", cand);
    if (!cand) {
        std::cerr << "No TTree 'Candidates' in " << inputName << "\n";
        return 1;
    }
    if (trainWeights == "xs" && !cand->GetBranch("xsWeight")) {
        std::cerr << "No xsWeight branch: rebuild the candidates with the "
                     "current buildtree\n";
        return 1;
    }
    for (const std::string& v : vars)
        if (!cand->GetBranch(v.c_str())) {
            std::cerr << "Branch '" << v << "' not found in Candidates\n";
            return 1;
        }

    gROOT->cd();   // copies go to memory, not into the (read-only) input file
    TTree* trainTree = cand->CopyTree("(event%2)==0");
    TTree* testTree  = cand->CopyTree("(event%2)==1");

    const Long64_t nSigTrain = trainTree->GetEntries("isSignal==1");
    const Long64_t nBkgTrain = trainTree->GetEntries("isSignal==0");
    const Long64_t nSigTest  = testTree->GetEntries("isSignal==1");
    const Long64_t nBkgTest  = testTree->GetEntries("isSignal==0");

    std::cout << "Training rows: signal " << nSigTrain
              << ", background " << nBkgTrain << "\n"
              << "Testing  rows: signal " << nSigTest
              << ", background " << nBkgTest << "\n";

    if (nSigTrain < 100 || nSigTest < 100 || nBkgTrain < 100 || nBkgTest < 100) {
        std::cerr << "Too few signal or background rows to train. Generate "
                     "more events (or lower the preselection).\n";
        return 1;
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

    TMVA::DataLoader* loader = new TMVA::DataLoader(dataset);

    for (const std::string& v : vars)
        loader->AddVariable(v.c_str(), 'F');
    loader->AddSpectator("mass", "mass");

    loader->AddSignalTree(trainTree, 1.0, TMVA::Types::kTraining);
    loader->AddSignalTree(testTree,  1.0, TMVA::Types::kTesting);
    loader->AddBackgroundTree(trainTree, 1.0, TMVA::Types::kTraining);
    loader->AddBackgroundTree(testTree,  1.0, TMVA::Types::kTesting);

    loader->SetSignalWeightExpression(weightExpr.c_str());
    loader->SetBackgroundWeightExpression(weightExpr.c_str());

    TCut sigCut("isSignal>0.5");
    TCut bkgCut("isSignal<0.5");
    if (!extraCut.empty()) {
        sigCut = sigCut && TCut(extraCut.c_str());
        bkgCut = bkgCut && TCut(extraCut.c_str());
    }

    // NormMode: during training, signal and background get equal total
    // weight (the sample is ~1000:1 background). This does not change the
    // test-sample weights that applyTMVA uses for yields.
    loader->PrepareTrainingAndTestTree(sigCut, bkgCut,
        "NormMode=EqualNumEvents:!V");

    std::string opts;
    if (method == "BDT") {
        opts = "!H:!V:NTrees=" + std::to_string(nTrees) +
               ":MaxDepth=" + std::to_string(depth) +
               ":MinNodeSize=" + std::to_string(minNode) + "%" +
               ":BoostType=AdaBoost:AdaBoostBeta=0.5"
               ":UseBaggedBoost:BaggedSampleFraction=0.5"
               ":SeparationType=GiniIndex:nCuts=20";
    } else {
        opts = "!H:!V:NTrees=" + std::to_string(nTrees) +
               ":MaxDepth=" + std::to_string(depth) +
               ":MinNodeSize=" + std::to_string(minNode) + "%" +
               ":BoostType=Grad:Shrinkage=0.1"
               ":UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20";
    }

    // BDTG is also TMVA's kBDT type; only the boost options differ.
    factory.BookMethod(loader, TMVA::Types::kBDT, method.c_str(), opts.c_str());

    factory.TrainAllMethods();
    factory.TestAllMethods();
    factory.EvaluateAllMethods();

    outFile->Close();
    input.Close();

    std::cout << "\n===== trainTMVA summary =====\n"
              << "Input:    " << inputName << "\n"
              << "Method:   " << method << "  (" << nTrees << " trees, depth "
              << depth << ")\n"
              << "Training weights: " << weightExpr << "\n"
              << "Variables (" << vars.size() << "):";
    for (const std::string& v : vars) std::cout << " " << v;
    std::cout << "\nWeights:  " << dataset << "/weights/TMVAClassification_"
              << method << ".weights.xml\n"
              << "TMVA file: " << outName
              << "   (root -l -e 'TMVA::TMVAGui(\"" << outName << "\")')\n"
              << "Next:     ./applyTMVA " << inputName
              << " --method=" << method << "\n"
              << "=============================\n";
    return 0;
}
