"""Sensitivity of the pool CSV bytes to 1-ulp differences in libm results (cos, sin, cosh, sinh, exp)."""
import sys
import numpy as np
mode, out = sys.argv[1], sys.argv[2]
if mode == "perturb":
    rng = np.random.default_rng(7)          # separate stream: does not touch the production RNG
    def wrap(f):
        def g(x, *a, **k):
            r = f(x, *a, **k)
            if isinstance(r, np.ndarray) and r.dtype == np.float64 and r.size > 1:
                flip = rng.random(r.shape) < 0.1          # 10 % of the results move by one ulp
                r = np.where(flip, np.nextafter(r, np.where(rng.random(r.shape) < 0.5, np.inf, -np.inf)), r)
            return r
        return g
    import grendel.production.fonll.meson_sampler as ms, grendel.production.decay_engine.kinematics as kin
    for mod in (ms, kin):
        mod.np = type("np_proxy", (), {k: getattr(np, k) for k in dir(np) if not k.startswith("__")})
        for name in ("cos", "sin", "cosh", "sinh", "exp"):
            setattr(mod.np, name, staticmethod(wrap(getattr(np, name))))
from grendel.models.alp.production import main
main(["--mass", "0.22", "1.24", "3.42", "--n-pool", "1200000", "--high-pt-tilt-scale", "5", "--out-dir", out])
