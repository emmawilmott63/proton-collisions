#include "generateCommon.h"
#include "analysisUtils.h"

#include <cstdlib>
#include <iostream>
#include <string>

// INCLUSIVE ("normal") sample: all hard QCD processes, so charm
// appears at its natural, rare rate.
// Purpose: realistic signal-to-background for evaluating the analysis.
//
// Usage: ./jets-inc [nEvents] [seed] [outFile] [--pthat-min=X] [--pthat-max=X]
//   defaults: 100000  54321  data/d0_inclusive.root  pTHat 5 - 10 GeV/c
// (different default seed from the charm sample, so they are independent)

int main(int argc, char* argv[]) {

    const int nEvents =
        std::atoi(getPositional(argc, argv, 0, "5000000").c_str());
    const int seed =
        std::atoi(getPositional(argc, argv, 1, "54321").c_str());
    const std::string outName =
        getPositional(argc, argv, 2, "data/d0_inclusive.root");

    const double pTHatMin = getNumber(argc, argv, "pthat-min", 5.0);
    const double pTHatMax = getNumber(argc, argv, "pthat-max", 100.0);

    if (pTHatMax <= pTHatMin) {
        std::cerr << "--pthat-max must be larger than --pthat-min\n";
        return 1;
    }

    return generateSample(SampleMode::Inclusive, nEvents, seed, outName,
                          pTHatMin, pTHatMax);
}
