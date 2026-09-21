"""The digitisation pipelines behind ``curves/data``.

Every module is a command (``python -m curves.tools.<name>``) that rewrites
its part of ``curves/data`` from a published source: a figure's PDF path
operators (``vector_*``), a rendered raster (``raster_bc7``), the HNLimits
package (``hnlimits_pull``, ``hnlimits_process``) or an experiment's own
data files (``original_sources``). The source PDFs are not tracked; see
``curves/sources/README.md`` for where each one comes from.

``hnlimits_process --check`` is the one that runs in CI: it regenerates the
transformed HNLimits regions and compares them byte for byte with the
tracked files.
"""
