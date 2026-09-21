"""Batch ray-cast of particle directions through the detector mesh.

Every model needs, per four-vector, whether it crosses the fiducial air
volume and at which distances it enters and leaves. That depends only on the
direction, so the result is cached per mass point as an ``.npz`` next to the
model's other outputs; the cache is invalidated when the source four-vectors
are newer than it. Geometry-code changes are not tracked by mtime: pass
``force=True`` after editing the mesh.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from ..constants import CMS_ORIGIN

# intersects_location allocates per-ray x per-triangle intermediates, so
# casting every candidate at once spikes to several GB on the highest-
# multiplicity points. Rays are independent: batching is exact.
RAY_CHUNK = 25000


def get_mesh():
    """The fiducial detector mesh (built on first import of the geometry)."""
    from .grendel_geometry import mesh_fiducial
    return mesh_fiducial


def directions_from_eta_phi(eta, phi):
    """(eta, phi) arrays -> (N, 3) unit direction vectors."""
    theta = 2.0 * np.arctan(np.exp(-eta))
    dx = np.sin(theta) * np.cos(phi)
    dy = np.sin(theta) * np.sin(phi)
    dz = np.cos(theta)
    return np.column_stack([dx, dy, dz])


def compute_geometry(eta, phi, mesh, origin=CMS_ORIGIN, batch_label=""):
    """Ray-cast (eta, phi) directions against the mesh.

    Returns ``(hits, entry_d, exit_d)``: a boolean mask of directions that
    cross the volume, and the distances from ``origin`` to the first two
    intersections (NaN where there is no hit).
    """
    n = len(eta)
    origin_arr = np.array(origin, dtype=np.float64)
    hits = np.zeros(n, dtype=bool)
    entry_d = np.full(n, np.nan)
    exit_d = np.full(n, np.nan)

    directions = directions_from_eta_phi(eta, phi)
    candidates = np.where(directions[:, 1] > 0.01)[0]
    n_cand = len(candidates)
    if n_cand == 0:
        return hits, entry_d, exit_d

    print(f"  Batch ray-casting {n_cand}/{n} candidates {batch_label}...",
          flush=True)

    cand_dirs = directions[candidates]
    loc_parts, rid_parts = [], []
    for cs in range(0, n_cand, RAY_CHUNK):
        ce = min(cs + RAY_CHUNK, n_cand)
        chunk_dirs = cand_dirs[cs:ce]
        chunk_origins = np.tile(origin_arr, (len(chunk_dirs), 1))
        loc_c, rid_c, _ = mesh.ray.intersects_location(
            ray_origins=chunk_origins, ray_directions=chunk_dirs)
        if len(loc_c):
            loc_parts.append(loc_c)
            rid_parts.append(rid_c + cs)  # ray ids are chunk-local

    if not loc_parts:
        print(f"  0/{n} events hit detector", flush=True)
        return hits, entry_d, exit_d
    locations = np.concatenate(loc_parts)
    ray_ids = np.concatenate(rid_parts)

    dists = np.linalg.norm(locations - origin_arr, axis=1)
    order = np.argsort(ray_ids)
    sorted_ray_ids = ray_ids[order]
    sorted_dists = dists[order]

    unique_rays, start_idx, counts = np.unique(
        sorted_ray_ids, return_index=True, return_counts=True)

    valid = counts >= 2
    valid_rays = unique_rays[valid]
    valid_starts = start_idx[valid]
    valid_counts = counts[valid]

    for i in range(len(valid_rays)):
        ray_local = valid_rays[i]
        orig_idx = candidates[ray_local]
        s = valid_starts[i]
        e = s + valid_counts[i]
        ray_dists = sorted_dists[s:e]
        ray_dists.sort()
        hits[orig_idx] = True
        entry_d[orig_idx] = ray_dists[0]
        exit_d[orig_idx] = ray_dists[1]

    n_hits = int(hits.sum())
    print(f"  {n_hits}/{n} events hit detector ({n_hits / n * 100:.2f}%)",
          flush=True)
    if n_hits > 0:
        path_lens = exit_d[hits] - entry_d[hits]
        print(f"  Mean path length: {path_lens.mean():.2f} m", flush=True)
    return hits, entry_d, exit_d


def load_or_compute_geometry(cache_path, eta, phi, mesh, *, origin=CMS_ORIGIN,
                             force=False, source_mtime=None, batch_label=""):
    """Load the cached ray-cast at ``cache_path`` or compute and cache it.

    The cache is served only if it exists, ``force`` is not set, and it is not
    older than ``source_mtime`` (the four-vector file it was computed from).
    The file is written atomically so an interrupted run never leaves a
    truncated cache behind.
    """
    cache_path = Path(cache_path)
    fresh = cache_path.exists() and not force
    if fresh and source_mtime is not None:
        fresh = cache_path.stat().st_mtime >= source_mtime
    if fresh:
        data = np.load(cache_path)
        return data["hits"].astype(bool), data["entry_d"], data["exit_d"]

    hits, entry_d, exit_d = compute_geometry(eta, phi, mesh, origin,
                                             batch_label=batch_label)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache_path.with_suffix(cache_path.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, hits=hits, entry_d=entry_d, exit_d=exit_d)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, cache_path)
    return hits, entry_d, exit_d
