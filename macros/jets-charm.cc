#include "generateCommon.h"
#include "analysisUtils.h"

#include <cstdlib>
#include <iostream>
#include <string>

// Charm-ENRICHED sample: every event contains a c cbar hard process.
// Purpose: confirm the D0 -> K pi reconstruction works and that the
// charm signal is visible.
//
// Usage: ./jets-charm [nEvents] [seed] [outFile] [--pthat-min=X] [--pthat-max=X]
//   defaults: 100000  12345  data/d0_charm.root  pTHat 5 - 10 GeV/c

int main(int argc, char* argv[]) {

    const int nEvents =
        std::atoi(getPositional(argc, argv, 0, "100000").c_str());
    const int seed =
        std::atoi(getPositional(argc, argv, 1, "12345").c_str());
    const std::string outName =
        getPositional(argc, argv, 2, "data/d0_charm.root");

    const double pTHatMin = getNumber(argc, argv, "pthat-min", 5.0);
    const double pTHatMax = getNumber(argc, argv, "pthat-max", 10.0);

    if (pTHatMax <= pTHatMin) {
        std::cerr << "--pthat-max must be larger than --pthat-min\n";
        return 1;
    }

    return generateSample(SampleMode::Charm, nEvents, seed, outName,
                          pTHatMin, pTHatMax);
}
