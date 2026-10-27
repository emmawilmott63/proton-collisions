#ifndef TMVA_COMMON_H
#define TMVA_COMMON_H

// Shared by trainTMVA.cc and applyTMVA.cc so that training and application
// use EXACTLY the same input variables in the same order (TMVA::Reader
// requires this).
//
//   --vars=a,b,c        use exactly these branches
//   --drop=a,b          remove these from the default list
//   --drop-group=g,h    remove whole feature groups (see featureGroup below):
//                       daughters, pair, jet, isolation, event
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

// Feature groups, used for group permutation importance (importanceTMVA) and
// for ablation (--drop-group). Features that carry overlapping information
// are in the same group, so removing or shuffling a group measures what the
// BDT would lose without that kind of information.
inline std::string featureGroup(const std::string& v) {
    static const char* daughters[] = {"ptK", "ptPi", "etaK", "etaPi", "ptHard",
                                      "ptSoft", "ptSum", "ptBalance"};
    static const char* pair[]      = {"dEta", "dPhi", "dR", "openAngle",
                                      "pairPt", "pairEta", "cosThetaStar"};
    static const char* jet[]       = {"sameJet", "jetPt", "jetEta", "jetMass",
                                      "jetNConst", "zJet", "dRJet"};
    static const char* iso[]       = {"isoRel02", "isoRel03", "isoRel",
                                      "isoRel06", "isoNeutralRel"};
    static const char* evt[]       = {"nCharged", "nJets", "sphericityT",
                                      "thrustT"};
    for (const char* n : daughters) if (v == n) return "daughters";
    for (const char* n : pair)      if (v == n) return "pair";
    for (const char* n : jet)       if (v == n) return "jet";
    for (const char* n : iso)       if (v == n) return "isolation";
    for (const char* n : evt)       if (v == n) return "event";
    return "other";
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

    for (const std::string& g :
         splitCommas(getOption(argc, argv, "drop-group", "")))
        vars.erase(std::remove_if(vars.begin(), vars.end(),
                       [&g](const std::string& v) {
                           return featureGroup(v) == g;
                       }),
                   vars.end());

    return vars;
}

#endif
