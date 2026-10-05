#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TTree.h"

int main() {
    Pythia8::Pythia p;
    p.readString("HardQCD:all = on");
    p.readString("Beams:eCM = 200.");
    p.readString("PhaseSpace:pTHatMin = 5.");
    p.readString("PhaseSpace:pTHatMax = 10.");
    p.init();

    TFile out("data/d0.root", "RECREATE");
    TTree tree("D0", "D0 candidates");

    int event, d0id;
    double d0px, d0py, d0pz, d0e;
    double kpx, kpy, kpz, ke;
    double pipx, pipy, pipz, pie;

    tree.Branch("event", &event);
    tree.Branch("d0id", &d0id);
    tree.Branch("d0px", &d0px);
    tree.Branch("d0py", &d0py);
    tree.Branch("d0pz", &d0pz);
    tree.Branch("d0e", &d0e);
    tree.Branch("kpx", &kpx);
    tree.Branch("kpy", &kpy);
    tree.Branch("kpz", &kpz);
    tree.Branch("ke", &ke);
    tree.Branch("pipx", &pipx);
    tree.Branch("pipy", &pipy);
    tree.Branch("pipz", &pipz);
    tree.Branch("pie", &pie);

    for (event = 0; event < 10000; ++event) {
        if (!p.next()) continue;

        for (int i = 0; i < p.event.size(); ++i) {
            if (p.event[i].id() != 421) continue;

            int k = p.event[i].daughter1();
            int pi = p.event[i].daughter2();

            if (k <= 0 || pi <= 0) continue;
            if (p.event[k].id() != -321) continue;
            if (p.event[pi].id() != 211) continue;

            d0id = 421;

            d0px = p.event[i].px();
            d0py = p.event[i].py();
            d0pz = p.event[i].pz();
            d0e  = p.event[i].e();

            kpx = p.event[k].px();
            kpy = p.event[k].py();
            kpz = p.event[k].pz();
            ke  = p.event[k].e();

            pipx = p.event[pi].px();
            pipy = p.event[pi].py();
            pipz = p.event[pi].pz();
            pie  = p.event[pi].e();

            tree.Fill();
        }
    }

    tree.Write();
    out.Close();

    return 0;
}
