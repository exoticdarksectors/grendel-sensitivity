# Vendored FairShip HNL decay modules

These are the FairShip (github.com/ShipSoft/FairShip) Heavy-Neutral-Lepton
physics modules, used here to sample **flavor-dependent** HNL rest-frame decays
(the visible final states differ for Ue / Umu / Utau). They were imported via
the project's `hnl_alaship` package, which vendored them to avoid a full
FairShip install.

## Files

| file | role |
|------|------|
| `hnl.py` | `HNLbranchings` / `HNL`: partial + total widths and lifetime for arbitrary couplings `[U2e, U2mu, U2tau]` (analytic; no Pythia) |
| `readDecayTable.py` | parse the decay-selection config and register HNL decay channels |
| `pythia8_conf_utils.py` | Pythia8 configuration helpers (needs numpy + scipy) |
| `shipunit.py` | FairShip unit constants |
| `alpha_s.dat` | running strong-coupling table used by `hnl.py` |

The decay-channel selection config consumed by these modules lives next to the
driver at `hnl/analysis/fairship_decay_selection.conf`.

## How they are used

`hnl/analysis/fairship_decay.py` drives these modules through `ROOT.TPythia8`
to sample rest-frame decay templates; `hnl/analysis/generate_decay_templates.py`
caches them to `.npz`. This **template-generation stage requires PyROOT +
Pythia8** and runs in a dedicated Homebrew-ROOT venv (NOT the conda `hnl` env);
see `hnl/README.md`. The downstream acceptance (`decay_reco_acceptance.py`)
consumes the cached `.npz` in the conda env with no ROOT dependency.

The modules are unmodified copies. Lifetime (`computeNLifetime`) from `hnl.py`
is also the flavor-aware source for the HNL ctau used by the analysis.
