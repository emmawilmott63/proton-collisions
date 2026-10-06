#include "generateCommon.h"

#include <cstdlib>
#include <string>

// INCLUSIVE ("normal") sample: all hard QCD processes, so charm
// appears at its natural, rare rate.
// Purpose: realistic signal-to-background for evaluating the analysis.
//
// Usage: ./generateD0Inclusive [nEvents] [seed] [outFile]
//   defaults: 100000  54321  data/d0_inclusive.root
// (different default seed from the charm sample, so they are independent)

int main(int argc, char* argv[]) {

    const int         nEvents = (argc > 1) ? std::atoi(argv[1]) : 100000;
    const int         seed    = (argc > 2) ? std::atoi(argv[2]) : 54321;
    const std::string outName = (argc > 3) ? argv[3] : "data/d0_inclusive.root";

    return generateSample(SampleMode::Inclusive, nEvents, seed, outName);
}
