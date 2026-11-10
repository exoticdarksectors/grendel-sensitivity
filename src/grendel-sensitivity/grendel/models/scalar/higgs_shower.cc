// Parton shower of the BC5 Higgs LHE events with Pythia 8 (ISR + FSR, no MPI,
// no hadronisation, Higgs stable). Streams the final Higgs of every event as
// "H pT y E px py pz" to stdout, then the run summary to stderr.
//   higgs_shower <plain LHE file> <seed> [max events, 0 = all]
#include "Pythia8/Pythia.h"
#include <cstdio>
#include <string>
using namespace Pythia8;

int main(int argc, char* argv[]) {
  if (argc < 3) { fprintf(stderr, "usage: higgs_shower <lhe> <seed> [max_events]\n"); return 2; }
  std::string lhe = argv[1];
  long seed = atol(argv[2]);
  long maxEvents = (argc > 3) ? atol(argv[3]) : 0;

  Pythia pythia("", false);
  for (const char* s : {"Beams:frameType = 4", "PartonLevel:ISR = on", "PartonLevel:FSR = on",
                        "PartonLevel:MPI = off", "HadronLevel:all = off", "25:mayDecay = off",
                        "ProcessLevel:resonanceDecays = off", "Check:event = off",
                        "Next:numberCount = 0", "Next:numberShowEvent = 0", "Next:numberShowInfo = 0",
                        "Next:numberShowProcess = 0", "Next:numberShowLHA = 0",
                        "Init:showChangedSettings = off", "Init:showChangedParticleData = off",
                        "Print:quiet = on", "Random:setSeed = on"})
    pythia.readString(s);
  pythia.readString("Beams:LHEF = " + lhe);
  pythia.readString("Random:seed = " + std::to_string(seed));
  if (!pythia.init()) { fprintf(stderr, "pythia init failed on %s\n", lhe.c_str()); return 1; }

  long nDone = 0, nHiggs = 0;
  while (maxEvents <= 0 || nDone < maxEvents) {
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) break;
      continue;
    }
    ++nDone;
    int iH = -1;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].id() == 25 && pythia.event[i].isFinal()) iH = i;
    if (iH < 0) continue;
    const Particle& h = pythia.event[iH];
    printf("H %.8e %.8e %.8e %.8e %.8e %.8e\n", h.pT(), h.y(), h.e(), h.px(), h.py(), h.pz());
    ++nHiggs;
  }
  fprintf(stderr, "VERSION %.3f\n", pythia.parm("Pythia:versionNumber"));
  fprintf(stderr, "NSHOWERED %ld\n", nDone);
  fprintf(stderr, "NHIGGS %ld\n", nHiggs);
  fprintf(stderr, "SIGMAGEN_MB %.10e\n", pythia.info.sigmaGen());
  return 0;
}
