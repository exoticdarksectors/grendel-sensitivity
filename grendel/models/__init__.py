"""The benchmark models.

``hnl``     BC6-BC8 heavy neutral leptons (mixing |U|^2 with e, mu or tau)
``scalar``  BC4 Higgs-portal dark scalar (mixing sin^2 theta)
``alp``     BC10 fermiophilic axion-like particle (1/f)

Each package has a ``spec`` (its ``ModelSpec`` for ``grendel.scan``), a
``scan`` command-line entry point, its physics (widths, branching ratios,
production) and the generators for its rest-frame decay templates.
"""
