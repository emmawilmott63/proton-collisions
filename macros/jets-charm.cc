#include "generateCommon.h"

#include <cstdlib>
#include <string>

// Charm-ENRICHED sample: every event contains a c cbar hard process.
// Purpose: confirm the D0 -> K pi reconstruction works and that the
// charm signal is visible.
//
// Usage: ./jets-charm [nEvents] [seed] [outFile]
//   defaults: 100000  12345  data/d0_charm.root

int main(int argc, char* argv[]) {

    const int         nEvents = (argc > 1) ? std::atoi(argv[1]) : 100000;
    const int         seed    = (argc > 2) ? std::atoi(argv[2]) : 12345;
    const std::string outName = (argc > 3) ? argv[3] : "data/d0_charm.root";

    return generateSample(SampleMode::Charm, nEvents, seed, outName);
}
