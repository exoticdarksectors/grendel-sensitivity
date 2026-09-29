"""Batch ray-cast of particle directions through the detector mesh."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

import numpy as np

from ..constants import CMS_ORIGIN

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
    """Ray-cast (eta, phi) directions against the mesh."""
    n = len(eta)
    origin_arr = np.array(origin, dtype=np.float64)
    directions = directions_from_eta_phi(eta, phi)
    candidates = np.where(directions[:, 1] > 0.01)[0]
    n_cand = len(candidates)
    if n_cand == 0:
        return _segments(n, np.empty(0, dtype=int), np.empty(0))

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
            rid_parts.append(rid_c + cs)

    if not loc_parts:
        print(f"  0/{n} events hit detector", flush=True)
        return _segments(n, np.empty(0, dtype=int), np.empty(0))
    locations = np.concatenate(loc_parts)
    ray_ids = candidates[np.concatenate(rid_parts)]
    dists = np.linalg.norm(locations - origin_arr, axis=1)
    hits, entry_d, exit_d = _segments(n, ray_ids, dists)

    n_hits = int(hits.sum())
    print(f"  {n_hits}/{n} events hit detector ({n_hits / n * 100:.2f}%)",
          flush=True)
    if n_hits > 0:
        path_lens = path_length(entry_d[hits], exit_d[hits])
        n_multi = int((np.isfinite(entry_d[hits, 1:])).any(axis=1).sum()) if entry_d.shape[1] > 1 else 0
        print(f"  Mean path length: {path_lens.mean():.2f} m"
              + (f" ({n_multi} rays cross the volume more than once)" if n_multi else ""),
              flush=True)
    return hits, entry_d, exit_d


def _segments(n, ray_ids, dists):
    """Pair each ray's sorted intersections into in-volume segments."""
    entry_d = np.full((n, 1), np.nan)
    exit_d = np.full((n, 1), np.nan)
    hits = np.zeros(n, dtype=bool)
    if len(ray_ids) == 0:
        return hits, entry_d, exit_d
    order = np.lexsort((dists, ray_ids))
    ray_ids, dists = ray_ids[order], dists[order]
    rays, start, counts = np.unique(ray_ids, return_index=True, return_counts=True)
    pairs = counts // 2
    k_max = max(1, int(pairs.max()))
    entry_d = np.full((n, k_max), np.nan)
    exit_d = np.full((n, k_max), np.nan)
    rank = np.arange(len(ray_ids)) - np.repeat(start, counts)
    keep = rank < 2 * np.repeat(pairs, counts)
    rid, seg, d, is_exit = ray_ids[keep], rank[keep] // 2, dists[keep], rank[keep] % 2 == 1
    entry_d[rid[~is_exit], seg[~is_exit]] = d[~is_exit]
    exit_d[rid[is_exit], seg[is_exit]] = d[is_exit]
    hits[rays[pairs >= 1]] = True
    return hits, entry_d, exit_d


def path_length(entry_d, exit_d):
    """Total in-volume path length per ray from ``(n,)`` or ``(n, n_segments)`` entry/exit arrays
    (unused segment slots count zero)."""
    entry_d = np.asarray(entry_d, dtype=float)
    exit_d = np.asarray(exit_d, dtype=float)
    if entry_d.ndim == 1:
        return exit_d - entry_d
    return np.nansum(exit_d - entry_d, axis=1)


def sample_decay_distances(entry_d, exit_d, n_samples, rng):
    """``(n, n_samples)`` decay distances uniform over each ray's in-volume path."""
    entry_d = np.asarray(entry_d, dtype=float)
    exit_d = np.asarray(exit_d, dtype=float)
    if entry_d.ndim == 1:
        entry_d, exit_d = entry_d[:, None], exit_d[:, None]
    d = rng.uniform(entry_d[:, :1], exit_d[:, :1], size=(len(entry_d), n_samples))
    if entry_d.shape[1] == 1:
        return d
    lengths = np.nan_to_num(exit_d - entry_d)
    multi = (lengths[:, 1:] > 0).any(axis=1)
    if multi.any():
        e, L = entry_d[multi], lengths[multi]
        u = (d[multi] - e[:, :1]) / L[:, :1] * L.sum(axis=1, keepdims=True)
        cum = np.cumsum(L, axis=1)
        k = np.minimum((u[:, :, None] >= cum[:, None, :]).sum(axis=2), L.shape[1] - 1)
        before = cum - L
        d[multi] = (np.take_along_axis(e, k, axis=1)
                    + u - np.take_along_axis(before, k, axis=1))
    return d


def _fingerprint(eta, phi, origin, mesh):
    """SHA-256 of everything a cached ray-cast depends on."""
    h = hashlib.sha256(b"grendel-raycast-v2")
    for arr in (np.asarray(eta, dtype=np.float64), np.asarray(phi, dtype=np.float64),
                np.asarray(origin, dtype=np.float64), np.asarray(mesh.vertices, dtype=np.float64),
                np.asarray(mesh.faces, dtype=np.int64)):
        h.update(np.ascontiguousarray(arr).tobytes())
    return h.hexdigest()


def _write_cache(cache_path, hits, entry_d, exit_d, fingerprint):
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache_path.with_suffix(cache_path.suffix + f".tmp.{os.getpid()}")
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, hits=hits, entry_d=entry_d, exit_d=exit_d,
                            fingerprint=np.array(fingerprint))
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, cache_path)


def load_or_compute_geometry(cache_path, eta, phi, mesh, *, origin=CMS_ORIGIN,
                             force=False, source_mtime=None, batch_label=""):
    """Load the cached ray-cast at ``cache_path`` or compute and cache it."""
    cache_path = Path(cache_path)
    fingerprint = _fingerprint(eta, phi, origin, mesh)
    fresh = cache_path.exists() and not force
    if fresh and source_mtime is not None:
        fresh = cache_path.stat().st_mtime >= source_mtime
    if fresh:
        with np.load(cache_path) as data:
            if "fingerprint" in data.files and str(data["fingerprint"]) == fingerprint:
                return data["hits"].astype(bool), data["entry_d"], data["exit_d"]

    hits, entry_d, exit_d = compute_geometry(eta, phi, mesh, origin,
                                             batch_label=batch_label)
    _write_cache(cache_path, hits, entry_d, exit_d, fingerprint)
    return hits, entry_d, exit_d
