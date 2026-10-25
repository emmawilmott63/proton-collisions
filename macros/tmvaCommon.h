#ifndef TMVA_COMMON_H
#define TMVA_COMMON_H

// Shared by trainTMVA.cc and applyTMVA.cc so that training and application
// use EXACTLY the same input variables in the same order (TMVA::Reader
// requires this).
//
//   --vars=a,b,c    use exactly these branches
//   --drop=a,b      remove these from the default list
//
// "mass" is deliberately never in the default list (see buildtree.cc):
// it is registered as a spectator only.

#include "analysisUtils.h"

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

inline std::vector<std::string> splitCommas(const std::string& s) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ','))
        if (!item.empty())
            out.push_back(item);
    return out;
}

inline std::vector<std::string> getTMVAVars(int argc, char* argv[]) {
    std::vector<std::string> vars = {
        // daughters / pair
        "ptK", "ptPi", "etaK", "etaPi", "dEta", "dPhi", "dR", "openAngle",
        "ptBalance", "pairPt", "pairEta", "cosThetaStar",
        // jet association
        "sameJet", "jetPt", "jetMass", "jetNConst", "zJet", "dRJet",
        // isolation
        "isoRel02", "isoRel03", "isoRel", "isoRel06", "isoNeutralRel",
        // event level
        "nCharged", "nJets", "sphericityT", "thrustT"
    };

    const std::string custom = getOption(argc, argv, "vars", "");
    if (!custom.empty())
        vars = splitCommas(custom);

    for (const std::string& d : splitCommas(getOption(argc, argv, "drop", "")))
        vars.erase(std::remove(vars.begin(), vars.end(), d), vars.end());

    return vars;
}

#endif
