"""HNL panels: HNLimits single-flavour exclusions overlaid with GRENDEL.

The experimental bounds are the HNLimits curves transformed with the
metadata of the upstream workbook ``local_HNL_database.xlsx`` (mass units,
CL normalisation, Dirac/Majorana conversion); the raw ``.dat`` files alone
lose those transformations. The workbook is parsed with the standard
library, so ``openpyxl`` is not needed. The cosmological bound replaces
HNLimits' with a single-Majorana BBN calculation (``COSMO_PER_FLAVOUR_FILES``).
``processed_curves`` returns the regions drawn; ``curves.tools.hnlimits_process``
writes them out as the tracked ``.dat`` files under ``curves/data/hnl``.
"""
from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import numpy as np

from .common import (GRENDEL_DARK_RED, GRENDEL_LW, GRENDEL_RED, GRENDEL_THRESHOLD_LW, PROJECTION_LW,
                     PROJECTION_STYLES, SIGNAL_THRESHOLD, parse_csv_bool, parse_csv_float)
from .data import HNL_DATA, HNLIMITS_DATABASE, HNLIMITS_DATABASE_URL

X_RANGE = (1.0e-3, 1.0e2)
ZOOM_X_RANGE = (0.3, 8.0)
ZOOM_X_TICKS = (0.3, 0.5, 1.0, 2.0, 5.0, 8.0)
Y_RANGE = (1.0e-10, 1.0e-1)
OPEN_UPPER = Y_RANGE[1]
GRENDEL_RUN_LABEL = (r"HL-LHC, $\mathcal{L}=3\ \mathrm{ab}^{-1}$" "\n"
                     r"$pp$, $\sqrt{s}=14\ \mathrm{TeV}$" "\n"
                     r"$N_{\mathrm{sig}}\geq 3$ (bkg-free)")

SCENARIOS = {
    "100": {"sheet": "Ue4", "flavor": "Ue", "latex": r"$|U_{eN}|^2$",
            "description": "single electron-flavor dominance"},
    "010": {"sheet": "Umu4", "flavor": "Umu", "latex": r"$|U_{\mu N}|^2$",
            "description": "single muon-flavor dominance"},
    "001": {"sheet": "Utau4", "flavor": "Utau", "latex": r"$|U_{\tau N}|^2$",
            "description": "single tau-flavor dominance"},
}

# HNLimits' `cosmo` row is the BBN bound of Sabti et al. (arXiv:2006.07387),
# computed for a Dirac HNL = two degenerate Majorana HNLs; its Sec. 5.2 says
# "the bounds do not apply to models with single Majorana HNLs", which is the
# paper's benchmark. (The row also points every flavour at Ue4_*_Cosmo.dat,
# whose `_bottom` is a flat two-point line, and tags the curve Dirac and
# 2 sigma, so HNLimits' rules would scale it by 1/2 x sqrt(2.71/4) = 0.41 on
# top of a curve that is already per Majorana state.) Each panel instead
# draws the single-Majorana BBN bound of Chen and Zhang (arXiv:2410.07343,
# 95% CL) for its own flavour, vector-extracted by curves.tools.vector_bbn_hnl
# and drawn as published: no CL or Dirac/Majorana rule. Paths are relative
# to curves/data/hnl.
COSMO_PER_FLAVOUR_FILES = {flavor: f"cosmology/BBN_ChenZhang2025_{flavor}.dat" for flavor in ("Ue", "Umu", "Utau")}


def hnlimits_data_root() -> Path:
    import HNLimits
    return Path(HNLimits.__file__).resolve().parent / "include" / "data"


# ---------------------------------------------------------------- workbook --

def _xlsx_col_index(cell_ref: str) -> int:
    letters = "".join(ch for ch in cell_ref if ch.isalpha())
    idx = 0
    for ch in letters:
        idx = idx * 26 + ord(ch.upper()) - ord("A") + 1
    return idx - 1


def read_xlsx_sheet(path: Path, sheet_name: str) -> list[dict[str, object]]:
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    rel_ns = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
    rels_ns = "{http://schemas.openxmlformats.org/package/2006/relationships}"
    with ZipFile(path) as zf:
        strings: list[str] = []
        shared = ET.fromstring(zf.read("xl/sharedStrings.xml"))
        for item in shared.findall(f"{ns}si"):
            strings.append("".join(t.text or "" for t in item.iter(f"{ns}t")))
        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        relmap = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels.findall(f"{rels_ns}Relationship")}
        target = None
        for sheet in workbook.find(f"{ns}sheets"):
            if sheet.attrib["name"] == sheet_name:
                target = relmap[sheet.attrib[f"{rel_ns}id"]]
                break
        if target is None:
            raise ValueError(f"Sheet {sheet_name!r} not found in {path}")
        sheet_xml = ET.fromstring(zf.read(f"xl/{target}"))
    rows: list[list[object]] = []
    for row in sheet_xml.findall(f".//{ns}row"):
        values: dict[int, object] = {}
        for cell in row.findall(f"{ns}c"):
            idx = _xlsx_col_index(cell.attrib.get("r", ""))
            value = cell.find(f"{ns}v")
            if value is None:
                values[idx] = None
                continue
            if cell.attrib.get("t") == "s":
                values[idx] = strings[int(value.text)]
            elif cell.attrib.get("t") == "b":
                values[idx] = value.text == "1"
            else:
                text = value.text or ""
                try:
                    number = float(text)
                    values[idx] = int(number) if number.is_integer() else number
                except ValueError:
                    values[idx] = text
        if values:
            rows.append([values.get(i) for i in range(max(values) + 1)])
    header = [str(v).strip() if v is not None else "" for v in rows[0]]
    records: list[dict[str, object]] = []
    for row in rows[1:]:
        row = row + [None] * (len(header) - len(row))
        record = {key: row[i] for i, key in enumerate(header) if key}
        if record.get("id"):
            records.append(record)
    return records


def clean_cell(value: object) -> object | None:
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    if isinstance(value, str):
        value = value.strip()
        return value if value else None
    return value


def load_metadata(sheet_name: str, database: Path = HNLIMITS_DATABASE) -> list[dict[str, object | None]]:
    if not database.is_file():
        raise FileNotFoundError(f"HNLimits metadata workbook not found: {database}\n"
                                f"Download it from {HNLIMITS_DATABASE_URL} or set GRENDEL_HNLIMITS_DATABASE.")
    return [{key: clean_cell(v) for key, v in row.items()} for row in read_xlsx_sheet(database, sheet_name)]


# ------------------------------------------------------------------ curves --

def data_file(rel_path: str, root: Path | None = None) -> Path:
    """The HNLimits data file a workbook row names.

    A few rows spell a directory or file in a different letter case from the
    shipped tree; a case-insensitive file system resolves them silently and a
    case-sensitive one does not, so the lookup falls back to a case-insensitive
    match component by component. The exact path is returned when it exists
    or when no single match is found (the caller reports the missing file).
    """
    root = hnlimits_data_root() if root is None else Path(root)
    exact = root / rel_path
    if exact.is_file():
        return exact
    current = root
    for part in Path(rel_path).parts:
        if (current / part).exists():
            current = current / part
            continue
        candidates = [c for c in current.iterdir() if c.name.lower() == part.lower()] if current.is_dir() else []
        if len(candidates) != 1:
            return exact
        current = candidates[0]
    return current


def load_curve(rel_path: object | None, record: dict[str, object | None]) -> np.ndarray | None:
    if rel_path is None:
        return None
    path = data_file(str(rel_path))
    if not path.is_file():
        raise FileNotFoundError(f"{record['id']}: missing HNLimits data file {path}")
    data = np.genfromtxt(path, comments="#")
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError(f"{record['id']}: expected two columns in {path}")
    return data[:, :2]


def converted_curve(rel_path: object | None, record: dict[str, object | None],
                    nature: str = "majorana") -> np.ndarray | None:
    from HNLimits.hnl_tools import cl90_dict, dirac_to_majorana_dic, unit_dict
    data = load_curve(rel_path, record)
    if data is None:
        return None
    mass = np.asarray(data[:, 0], dtype=float)
    u2 = np.asarray(data[:, 1], dtype=float)
    cl = record.get("CL")
    if cl in cl90_dict:
        u2 = u2 * cl90_dict[cl]
    if nature == "majorana" and record.get("hnl_type") == "Dirac":
        u2 = u2 * dirac_to_majorana_dic[str(record.get("type"))]
    units = record.get("units")
    if units not in unit_dict:
        raise ValueError(f"{record['id']}: unsupported mass unit {units!r}")
    mass = mass * unit_dict[units] / 1000.0
    good = np.isfinite(mass) & np.isfinite(u2) & (mass > 0) & (u2 > 0)
    return np.column_stack([mass[good], u2[good]])


def sorted_curve(curve: np.ndarray) -> np.ndarray:
    return curve[np.argsort(curve[:, 0])]


def log_interp(x: np.ndarray, y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    order = np.argsort(x)
    return 10.0 ** np.interp(np.log10(x_new), np.log10(x[order]), np.log10(y[order]))


def load_region_file(path: Path) -> dict[str, np.ndarray]:
    """The ``# bottom`` and ``# top`` sections of a local region file (GeV, U^2)."""
    sections: dict[str, list[list[float]]] = {}
    current = None
    for line in Path(path).read_text().splitlines():
        text = line.strip()
        if text in ("# bottom", "# top"):
            current = sections.setdefault(text[2:], [])
        elif text and not text.startswith("#") and current is not None:
            current.append([float(v) for v in text.split()[:2]])
    return {name: np.asarray(rows, dtype=float) for name, rows in sections.items()}


def processed_dir(scenario: str, config: dict[str, str]) -> Path:
    return HNL_DATA / f"hnlimits_{scenario}_{config['flavor']}_majorana"


def processed_curves(scenario: str, config: dict[str, str],
                     records: list[dict[str, object | None]]) -> list[dict[str, object]]:
    """The regions one scenario draws: HNLimits' transformed, and the local cosmology."""
    processed = []
    flavor = str(config["flavor"])
    for record in records:
        if record.get("id") == "cosmo" and flavor in COSMO_PER_FLAVOUR_FILES:
            rel_path = COSMO_PER_FLAVOUR_FILES[flavor]
            region = load_region_file(HNL_DATA / rel_path)
            bottom, top = region.get("bottom"), region.get("top")
            file_top = file_bottom = f"curves/data/hnl/{rel_path}"
        else:
            file_top, file_bottom = record.get("file_top"), record.get("file_bottom")
            bottom = converted_curve(file_bottom, record)
            top = converted_curve(file_top, record)
        if bottom is None and top is None:
            continue
        category = "cosmology" if record["id"] == "cosmo" else "experimental"
        processed.append({"record": record, "bottom": bottom, "top": top, "category": category,
                          "files": (file_top, file_bottom),
                          "path": processed_dir(scenario, config) / f"{record['id']}_{config['flavor']}_majorana.dat"})
    return processed


def processed_file_text(scenario: str, config: dict[str, str], item: dict[str, object],
                        database: Path = HNLIMITS_DATABASE) -> str:
    """The ``.dat`` record of one transformed region, with its provenance."""
    record = item["record"]
    file_top, file_bottom = item["files"]
    replaced = (file_top, file_bottom) != (record.get("file_top"), record.get("file_bottom"))
    lines = [f"# {record['plot_label']}  (HNLimits id: {record['id']})",
             f"# Scenario: {scenario}, {config['description']}",
             "# Nature: Majorana",
             f"# Category: {item['category']}",
             "# Source metadata: local_HNL_database.xlsx",
             f"# Source metadata SHA256: {hashlib.sha256(database.read_bytes()).hexdigest()}",
             f"# Source URL: {HNLIMITS_DATABASE_URL}",
             "# Transformations: " + ("none, drawn as published" if replaced
                                      else "HNLimits unit, CL, and Dirac/Majorana rules")]
    if replaced:
        lines += ["# Curves REPLACED locally: the HNLimits row is the Sabti et al.",
                  "#   (arXiv:2006.07387) BBN bound for two degenerate Majorana HNLs, which",
                  "#   does not apply to a single Majorana HNL. Drawn instead: the single-",
                  "#   Majorana BBN bound of Chen and Zhang (arXiv:2410.07343), 95% CL, with",
                  "#   no CL or Dirac/Majorana rescaling. See COSMO_PER_FLAVOUR_FILES in",
                  "#   curves/plot/hnl_panel.py and curves.tools.vector_bbn_hnl.",
                  f"# File used: {file_top}"]
    lines.append("# mass_GeV  U2")
    out = "\n".join(lines) + "\n"
    if item["bottom"] is not None:
        out += "# bottom\n" + "".join(f"{m:.8e}  {u:.8e}\n" for m, u in item["bottom"])
    if item["top"] is not None:
        out += "\n# top\n" + "".join(f"{m:.8e}  {u:.8e}\n" for m, u in item["top"])
    return out


# ------------------------------------------------------------------ drawing --

def plot_limit(ax, item: dict[str, object], label: str | None = None) -> None:
    from HNLimits import plot_tools
    record, bottom, top, category = item["record"], item["bottom"], item["top"], item["category"]
    fill = "#9EC9E2" if category == "cosmology" else "0.62"
    alpha = 0.46 if category == "cosmology" else 0.38
    zorder = 1 if category == "cosmology" else 0
    same_file = record.get("file_top") is not None and record.get("file_top") == record.get("file_bottom")
    if same_file and bottom is not None:
        xs, ys = plot_tools.get_ordered_closed_region((bottom[:, 0], bottom[:, 1]), logy=True)
        ax.fill(xs, ys, facecolor=fill, edgecolor="none", alpha=alpha, zorder=zorder, label=label)
        return
    if bottom is None:
        return
    b = sorted_curve(bottom)
    if top is None:
        x = np.geomspace(max(X_RANGE[0], b[:, 0].min()), min(X_RANGE[1], b[:, 0].max()), 800)
        y_low = log_interp(b[:, 0], b[:, 1], x)
        y_high = np.full_like(x, Y_RANGE[1])
    else:
        t = sorted_curve(top)
        lo = max(X_RANGE[0], b[:, 0].min(), t[:, 0].min())
        hi = min(X_RANGE[1], b[:, 0].max(), t[:, 0].max())
        if not lo < hi:
            return
        x = np.geomspace(lo, hi, 800)
        y_low = log_interp(b[:, 0], b[:, 1], x)
        y_high = log_interp(t[:, 0], t[:, 1], x)
    lower = np.minimum(y_low, y_high)
    upper = np.maximum(y_low, y_high)
    ax.fill_between(x, lower, upper, facecolor=fill, edgecolor="none", alpha=alpha, zorder=zorder, label=label)


def append_high_mass_tip(rows: list[dict[str, object]], threshold: float | None = None) -> list[dict[str, object]]:
    """Close the island where the optimised yield falls below threshold.

    The sensitivity CSV keeps rows above the last excluded point with
    ``has_sensitivity = False`` and finite ``peak_N``/``peak_u2``. The physical
    high-mass endpoint is where the two roots merge at the maximum of the
    yield scan, i.e. ``u2_min = u2_max = peak_u2`` and ``peak_N`` = threshold.
    """
    threshold = SIGNAL_THRESHOLD if threshold is None else float(threshold)
    if not rows or "peak_N" not in rows[0] or "peak_u2" not in rows[0]:
        return rows
    rows = sorted((row.copy() for row in rows), key=lambda row: float(row["mass_GeV"]))
    sensitive = [bool(row["has_sensitivity"]) for row in rows]
    if not any(sensitive):
        return []
    last_pos = max(i for i, s in enumerate(sensitive) if s)
    closed = [row.copy() for row in rows if row["has_sensitivity"]]
    if last_pos + 1 >= len(rows):
        return closed
    left = rows[last_pos]
    below = [row for row in rows[last_pos + 1:]
             if np.isfinite(row["peak_N"]) and np.isfinite(row["peak_u2"])
             and row["peak_N"] > 0 and row["peak_u2"] > 0 and row["peak_N"] < threshold]
    if not below:
        return closed
    right = below[0]
    if not (np.isfinite(left["peak_N"]) and np.isfinite(left["peak_u2"])
            and left["peak_N"] > threshold and left["peak_u2"] > 0 and right["mass_GeV"] > left["mass_GeV"]):
        return closed
    log_left = np.log10(float(left["peak_N"]))
    log_right = np.log10(float(right["peak_N"]))
    log_threshold = np.log10(threshold)
    if log_left == log_right:
        return closed
    frac = float(np.clip((log_threshold - log_left) / (log_right - log_left), 0.0, 1.0))
    tip_mass = float(left["mass_GeV"] + frac * (right["mass_GeV"] - left["mass_GeV"]))
    tip_u2 = 10.0 ** (np.log10(float(left["peak_u2"]))
                      + frac * (np.log10(float(right["peak_u2"])) - np.log10(float(left["peak_u2"]))))
    tip = left.copy()
    tip.update(mass_GeV=tip_mass, u2_min=tip_u2, u2_max=tip_u2, u2_min_open=False, u2_max_open=False,
               peak_N=threshold, peak_u2=tip_u2, has_sensitivity=True)
    return sorted([*closed, tip], key=lambda row: float(row["mass_GeV"]))


def load_flavor_rows(path: Path, flavor: str, close_tip: bool = True,
                     threshold: float | None = None) -> list[dict[str, object]]:
    """The rows of one flavor from a sensitivity CSV, tip closed at ``threshold``."""
    rows: list[dict[str, object]] = []
    with Path(path).open(newline="") as fh:
        for raw in csv.DictReader(fh):
            if raw.get("flavor") != flavor:
                continue
            rows.append({"flavor": raw["flavor"], "mass_GeV": parse_csv_float(raw.get("mass_GeV")),
                         "u2_min": parse_csv_float(raw.get("u2_min")), "u2_max": parse_csv_float(raw.get("u2_max")),
                         "u2_min_open": parse_csv_bool(raw.get("u2_min_open")),
                         "u2_max_open": parse_csv_bool(raw.get("u2_max_open")),
                         "peak_N": parse_csv_float(raw.get("peak_N")), "peak_u2": parse_csv_float(raw.get("peak_u2")),
                         "has_sensitivity": parse_csv_bool(raw.get("has_sensitivity"))})
    rows = sorted(rows, key=lambda row: float(row["mass_GeV"]))
    if close_tip:
        return append_high_mass_tip(rows, threshold=threshold)
    return [row for row in rows if row["has_sensitivity"]]


def plot_grendel(ax, rows: list[dict[str, object]], label: str = "GRENDEL\nthis work") -> None:
    if not rows:
        return
    x = np.asarray([row["mass_GeV"] for row in rows], dtype=float)
    lo = np.asarray([row["u2_min"] for row in rows], dtype=float)
    hi = np.asarray([row["u2_max"] for row in rows], dtype=float)
    max_open = np.asarray([row["u2_max_open"] for row in rows], dtype=bool)
    hi_fill = np.where(np.isfinite(hi), hi, np.where(max_open, OPEN_UPPER, np.nan))
    hi_line = np.where(max_open, np.nan, hi)
    ax.fill_between(x, lo, hi_fill, color=GRENDEL_RED, alpha=0.18, lw=0, zorder=8, label=label)
    ax.plot(x, lo, color=GRENDEL_RED, lw=GRENDEL_LW, zorder=9)
    ax.plot(x, hi_line, color=GRENDEL_RED, lw=GRENDEL_LW, zorder=9)


def plot_grendel_threshold(ax, rows: list[dict[str, object]], label: str, color: str = GRENDEL_DARK_RED,
                           lw: float = GRENDEL_THRESHOLD_LW, ls: object = "-") -> None:
    """The band re-solved at a stricter threshold, an unfilled solid dark-red
    contour on top of the baseline band it is nested inside."""
    if not rows:
        return
    x = np.asarray([row["mass_GeV"] for row in rows], dtype=float)
    lo = np.asarray([row["u2_min"] for row in rows], dtype=float)
    hi = np.asarray([row["u2_max"] for row in rows], dtype=float)
    max_open = np.asarray([row["u2_max_open"] for row in rows], dtype=bool)
    hi_line = np.where(max_open, np.nan, hi)
    ax.plot(x, lo, color=color, lw=lw, ls=ls, zorder=10, label=label)
    ax.plot(x, hi_line, color=color, lw=lw, ls=ls, zorder=10)


def load_projection_segments(path: Path) -> list[np.ndarray]:
    from .common import load_polygon_segments
    return load_polygon_segments(path)


def _drop_top_closure(segment, rtol: float = 0.02):
    """Split a contour at its flat top edge, yielding the parts below it.

    Some HNLimits projection contours are closed along the top of their
    source figure's frame at a round value such as 1e-3. That horizontal run
    is a frame edge, not a sensitivity boundary, so it is not drawn.
    """
    y = segment[:, 1]
    ymax = float(np.nanmax(y))
    keep = y < ymax * (1.0 - rtol)
    if keep.all():
        return [segment]
    parts, start = [], None
    for i, k in enumerate(keep):
        if k and start is None:
            start = i
        elif not k and start is not None:
            parts.append(segment[start:i])
            start = None
    if start is not None:
        parts.append(segment[start:])
    return [q for q in parts if len(q) > 1]


def plot_hnlimits_projection(ax, flavor: str, name: str, label: str, color: str, linestyle: object,
                             linewidth: float = PROJECTION_LW, drop_top_closure: bool = False) -> None:
    curve_path = HNL_DATA / f"hnlimits_projections_{flavor}" / f"{name}_{flavor}.dat"
    if not curve_path.is_file():
        return
    segments = []
    for segment in load_projection_segments(curve_path):
        if len(segment) == 0:
            continue
        segments.extend(_drop_top_closure(segment) if drop_top_closure else [segment])
    for idx, segment in enumerate(segments):
        ax.plot(segment[:, 0], segment[:, 1], color=color, lw=linewidth, ls=linestyle, zorder=7,
                label=label if idx == 0 else None)


SHIP_6E20 = {"Ue": ("SHiP_bb_Ue.dat", "arXiv:2511.03883"), "Umu": ("SHiP_bb_Umu.dat", "arXiv:2511.03883"),
             "Utau": ("SHiP_ecn3_Utau.dat", "arXiv:2310.17726")}


def plot_ship_projection(ax, flavor: str) -> None:
    # 6e20 protons on target as on the BC4/BC5/BC10 panels (curves.tools.vector_ship_hnl): the ESPP 2026
    # Briefing Book for e and mu, the PBC ECN3 report for tau. The HNLimits curves (arXiv:1811.00930) are 2e20.
    style = PROJECTION_STYLES["SHiP"]
    name, source = SHIP_6E20[flavor]
    for idx, segment in enumerate(load_projection_segments(HNL_DATA / "vector" / name)):
        ax.plot(segment[:, 0], segment[:, 1], color=style["color"], lw=PROJECTION_LW, ls=style["ls"], zorder=7,
                label=f"SHiP ($6\\times10^{{20}}$ PoT)\n{source}" if idx == 0 else None)


def plot_faser2_projection(ax, flavor: str) -> None:
    # arXiv:2405.07330 Fig. 9 (FPF geometry, kaon production included; curves.tools.raster_faser2_hnl), which
    # supersedes the 2018 strawman curves of HNLimits (arXiv:1811.12522). The line leaves its source frame at 1e-3.
    style = PROJECTION_STYLES["FASER2"]
    for idx, segment in enumerate(load_projection_segments(HNL_DATA / "raster" / f"FASER2_2405_{flavor}.dat")):
        ax.plot(segment[:, 0], segment[:, 1], color=style["color"], lw=PROJECTION_LW, ls=style["ls"], zorder=7,
                label="FASER2 (3 ab$^{-1}$)\narXiv:2405.07330" if idx == 0 else None)


def plot_anubis_projection(ax, flavor: str) -> None:
    curve_path = HNL_DATA / "original_sources" / "ANUBIS" / f"ANUBIS_2026_{flavor}.dat"
    if not curve_path.is_file():
        return
    for idx, segment in enumerate(load_projection_segments(curve_path)):
        if len(segment) == 0:
            continue
        ax.plot(segment[:, 0], segment[:, 1], color=PROJECTION_STYLES["ANUBIS"]["color"], lw=PROJECTION_LW,
                ls=PROJECTION_STYLES["ANUBIS"]["ls"], zorder=7,
                label="ANUBIS (3 ab$^{-1}$)\narXiv:2606.26862" if idx == 0 else None)


def plot_anubis_meson_projection(ax, flavor: str) -> None:
    if flavor != "Ue":
        return
    curve_path = HNL_DATA / "original_sources" / "ANUBIS" / "ANUBIS_WangZhang_2025_ceiling_Ue.dat"
    if not curve_path.is_file():
        return
    for idx, segment in enumerate(load_projection_segments(curve_path)):
        if len(segment) == 0:
            continue
        ax.plot(segment[:, 0], segment[:, 1], color="#8c564b", lw=PROJECTION_LW, ls=(0, (2, 2)), zorder=7,
                label="ANUBIS D/B meson\narXiv:2512.13011" if idx == 0 else None)


def plot_codexb_projection(ax, flavor: str) -> None:
    # Vector-extracted from the CODEX-b EoI's own per-flavour panels
    # (curves.tools.vector_codexb_hnl); column 2 is already converted from
    # the EoI's Dirac convention to Majorana (lower edge / sqrt2, upper / 2),
    # NaN on the frame vertex.
    curve_path = HNL_DATA / "original_sources" / "CODEX-b" / f"CODEX-b_2019_{flavor}.dat"
    if not curve_path.is_file():
        return
    for idx, segment in enumerate(load_projection_segments(curve_path)):
        segment = segment[np.isfinite(segment[:, 1])]
        if len(segment) == 0:
            continue
        ax.plot(segment[:, 0], segment[:, 1], color=PROJECTION_STYLES["CODEX-b"]["color"], lw=PROJECTION_LW,
                ls=PROJECTION_STYLES["CODEX-b"]["ls"], zorder=7,
                label="CODEX-b (300 fb$^{-1}$)\narXiv:1911.00481" if idx == 0 else None)
