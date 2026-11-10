# SensCalc quartic-coupling production tables (dark scalar, PBC BC5)

Files copied verbatim (2026-09-25) from SensCalc v1.3.3, byte-identical in v1.4
(commit 4f7210b, 2026-10-11)
(M. Ovchynnikov, J.-L. Tastet, O. Mikulenko, K. Bondarenko, arXiv:2305.13383;
<https://github.com/maksymovchynnikov/SensCalc>, BSD-3-Clause, `LICENSE`),
directory `phenomenology/Scalar/Production probabilities/`:

- `BplustoXsSStotal.dat` — BR(B+ -> X_s S S), the sum over the strange final
  states (K, K*, K0*, K1, K2*, ...) of the quartic-coupling channel;
- `BstoSS.dat` — BR(B_s -> S S);
- `BplustoKplusSS.dat` — the exclusive K+ channel, kept only as the
  normalisation check against arXiv:1904.10447 Table 4 (1.4e-10) and against
  the eq. (G.8) integral implemented in `grendel.models.scalar.model`
  (`tests/unit/scalar/test_bc5.py`).

Columns: `m_S [GeV]`, branching ratio at the quartic coupling `alpha = 1 GeV`
of Boiarska et al. arXiv:1904.10447 eq. (2.14), `L = (alpha/2) S^2 h`, for
which BR(h -> SS) ~ 0.02 (their eq. 2.15). The rates scale as `alpha^2`;
`model.alpha_quartic` fixes `alpha` from BR(h -> SS) = 0.01, the PBC BC5
convention. Zero beyond the last tabulated mass (the channel is closed).
