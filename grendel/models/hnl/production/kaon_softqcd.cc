// Pythia8 SoftQCD charged-kaon spectrum for the GRENDEL HNL kaon gate.
// pp @ 14 TeV, inelastic. Records final-state K+- (pT, y); prints the
// inelastic cross section and <n_K+->/event, then a (pT,y) list to stdout.
#include "Pythia8/Pythia.h"
#include <cstdio>
using namespace Pythia8;

int main(int argc, char* argv[]) {
  long nEvents = (argc > 1) ? atol(argv[1]) : 20000;
  long seed    = (argc > 2) ? atol(argv[2]) : 42;

  Pythia pythia;
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 14000.");
  pythia.readString("SoftQCD:inelastic = on");
  pythia.readString("Print:quiet = on");
  pythia.readString("Next:numberCount = 0");
  pythia.readString("Random:setSeed = on");
  char sbuf[64]; snprintf(sbuf, sizeof(sbuf), "Random:seed = %ld", seed);
  pythia.readString(sbuf);
  if (!pythia.init()) { fprintf(stderr, "pythia init failed\n"); return 1; }

  long nKaon = 0;
  // stream K+- (pT, y) to stdout after a header; cap |y|<8, pT<10 (soft window)
  for (long iEv = 0; iEv < nEvents; ++iEv) {
    if (!pythia.next()) continue;
    for (int i = 0; i < pythia.event.size(); ++i) {
      Particle& p = pythia.event[i];
      if (!p.isFinal()) continue;
      if (p.idAbs() != 321) continue;    // K+-
      double pT = p.pT();
      double y  = p.y();
      if (pT > 10.0 || fabs(y) > 8.0) continue;
      ++nKaon;
      printf("K %.5f %.5f\n", pT, y);
    }
  }
  double sigmaInel_mb = pythia.info.sigmaGen();       // mb
  double nKaonPerEvt  = double(nKaon) / double(nEvents);
  // header on stderr so stdout stays pure data
  fprintf(stderr, "SIGMA_INEL_MB %.6f\n", sigmaInel_mb);
  fprintf(stderr, "N_EVENTS %ld\n", nEvents);
  fprintf(stderr, "N_KAON %ld\n", nKaon);
  fprintf(stderr, "NKAON_PER_EVT %.6f\n", nKaonPerEvt);
  // sigma_kaon in pb = sigma_inel(mb) * 1e9 (pb/mb) * <n_K>
  fprintf(stderr, "SIGMA_KAON_PB %.6e\n", sigmaInel_mb * 1.0e9 * nKaonPerEvt);
  return 0;
}
