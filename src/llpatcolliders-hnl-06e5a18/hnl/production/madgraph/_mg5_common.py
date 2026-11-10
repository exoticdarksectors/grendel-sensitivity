"""Shared MG5 executable, LHAPDF, and macOS rpath helpers."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
PYTHON_EXE = Path(sys.executable)


def _resolve_mg5_exe():
    """Resolve mg5_aMC from env, vendored, shared, or sibling installs."""
    env_path = os.environ.get("HNL_MG5_EXE")
    if env_path:
        return Path(env_path)

    vendored = PROJECT_ROOT / "vendored" / "MG5_aMC_v3_6_6" / "bin" / "mg5_aMC"
    if vendored.exists():
        return vendored

    shared = PROJECT_ROOT.parent.parent / "vendored" / "MG5_aMC_v3_6_6" / "bin" / "mg5_aMC"
    if shared.exists():
        return shared

    rel = Path("llpatcolliders_FONLL") / "vendored" / "MG5_aMC_v3_6_6" / "bin" / "mg5_aMC"
    for base in (PROJECT_ROOT.parent, PROJECT_ROOT.parent.parent):
        candidate = base / rel
        if candidate.exists():
            return candidate
    return PROJECT_ROOT.parent.parent / rel


MG5_EXE = _resolve_mg5_exe()


def _resolve_lhapdf_config():
    """Resolve lhapdf-config from $HNL_LHAPDF_CONFIG or the active conda env.

    No sibling-directory fallback: if neither resolves, the returned path will
    not exist and the MG5 runners report a precise error before launching.
    """
    env_path = os.environ.get("HNL_LHAPDF_CONFIG")
    if env_path:
        return Path(env_path)
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        return Path(conda_prefix) / "bin" / "lhapdf-config"
    found = shutil.which("lhapdf-config")
    return Path(found) if found else Path("lhapdf-config")


LHAPDF_CONFIG = _resolve_lhapdf_config()
LHAPDF_LIBDIR = LHAPDF_CONFIG.parent.parent / "lib"

# PDF data set (NNPDF40_nlo_as_01180, ~340 MB). Defaults to the active conda
# env's LHAPDF data dir. Provision it once with
# `lhapdf install NNPDF40_nlo_as_01180`, or point $HNL_LHAPDF_DATA at an
# existing LHAPDF data directory (e.g. a prior FONLL install).
LHAPDF_DATA_DIR = Path(os.environ.get(
    "HNL_LHAPDF_DATA",
    LHAPDF_CONFIG.parent.parent / "share" / "LHAPDF",
))


def mg5_subprocess_env():
    """Inject LHAPDF libdir + PDF data path for MG5 subprocesses."""
    env = os.environ.copy()
    for key in ("DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"):
        existing = env.get(key, "")
        env[key] = f"{LHAPDF_LIBDIR}:{existing}" if existing else str(LHAPDF_LIBDIR)
    existing_data = env.get("LHAPDF_DATA_PATH", "")
    env["LHAPDF_DATA_PATH"] = (
        f"{LHAPDF_DATA_DIR}:{existing_data}" if existing_data else str(LHAPDF_DATA_DIR)
    )
    return env


def patch_me5_configuration(cards_dir):
    """Disable browser auto-open and point MadEvent at LHAPDF."""
    cfg = cards_dir / "me5_configuration.txt"
    if not cfg.exists():
        return
    text = cfg.read_text()
    new_text = text.replace("# web_browser = None", "web_browser = None")

    def set_option(config_text, key, value):
        lines = config_text.splitlines()
        updated = False
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith(f"{key} =") or stripped.startswith(f"# {key} ="):
                lines[i] = f"{key} = {value}"
                updated = True
        if not updated:
            lines.append(f"{key} = {value}")
        return "\n".join(lines) + "\n"

    new_text = set_option(new_text, "automatic_html_opening", "False")
    new_text = set_option(new_text, "notification_center", "False")

    lhapdf_line = f"lhapdf_py3 = {LHAPDF_CONFIG}"
    if "# lhapdf_py3 = lhapdf-config" in new_text:
        new_text = new_text.replace("# lhapdf_py3 = lhapdf-config", lhapdf_line)
    elif "lhapdf_py3 =" not in new_text:
        new_text += f"\n{lhapdf_line}\n"

    if new_text != text:
        cfg.write_text(new_text)


def force_compile_subprocesses(work_subdir, log_path, timeout=900):
    """Compile MG5 subprocesses before event generation via direct make.

    Bypasses the interactive bin/madevent shell (which in MG5 3.6.6 routes
    `compile` through the launch flow's switches/cards dialogs and chokes
    on piped stdin). We build the Source/ libs once, then `make madevent`
    in every SubProcesses/P*_*/ directly.

    A fresh MG5 output tree ships only matrix*_orig.f; the *_optim.f files
    are generated lazily by the helicity-recycling step inside the launch
    flow. We skip that optimization and pass MATRIX=matrix*_orig.o
    explicitly so the link picks up the original matrix elements (the
    speed cost is negligible for our event counts).
    """
    env = mg5_subprocess_env()
    with open(log_path, "w") as log:
        log.write("=== make -C Source ===\n")
        log.flush()
        result = subprocess.run(
            ["make"],
            cwd=work_subdir / "Source",
            env=env,
            stdout=log, stderr=subprocess.STDOUT, timeout=timeout,
        )
        if result.returncode != 0:
            return False

        any_built = False
        for sub in sorted((work_subdir / "SubProcesses").glob("P*_*")):
            if not sub.is_dir() or not (sub / "Makefile").exists():
                continue
            orig_objs = " ".join(sorted(
                p.with_suffix(".o").name for p in sub.glob("matrix*_orig.f")
            ))
            if not orig_objs:
                log.write(f"\n=== skip {sub.name}: no matrix*_orig.f ===\n")
                continue
            log.write(f"\n=== make -C {sub.name} madevent MATRIX={orig_objs} ===\n")
            log.flush()
            result = subprocess.run(
                ["make", "madevent", f"MATRIX={orig_objs}"],
                cwd=sub,
                env=env,
                stdout=log, stderr=subprocess.STDOUT, timeout=timeout,
            )
            if result.returncode != 0 or not (sub / "madevent").exists():
                return False
            any_built = True
    return any_built


def write_process_block(out_file, proc_lines):
    """Copy MG5 define/generate/add-process lines from a process card."""
    for line in proc_lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if (
            stripped.startswith("define ")
            or stripped.startswith("generate ")
            or stripped.startswith("add process ")
        ):
            out_file.write(line)


def has_five_flavor_proton(work_subdir):
    """Return True iff the cached MG5 process dir was built with b in p."""
    card = work_subdir / "Cards" / "proc_card_mg5.dat"
    if not card.exists():
        return False
    for line in card.read_text().splitlines():
        stripped = line.strip()
        if not stripped.startswith("define p"):
            continue
        if "=" not in stripped:
            continue
        tokens = stripped.split("=", 1)[1].split()
        return "b" in tokens
    return False


def patch_rpath_for_lhapdf(work_subdir):
    """macOS dyld workaround: symlink LHAPDF + OpenMP next to MG5 binaries.

    The MG5 madevent binary references @rpath/libLHAPDF.dylib and
    @rpath/libomp.dylib (when LHAPDF was installed via conda-forge, which
    ships libomp alongside libLHAPDF and the linker picks it up over GCC's
    libgomp). DYLD_LIBRARY_PATH is unreliable across MG5's grandchild
    Fortran subprocesses on macOS SIP, so we symlink each required dylib
    into the binary's own directory — dyld resolves @rpath through
    @loader_path and picks them up.
    """
    if sys.platform != "darwin":
        return
    sources = [LHAPDF_LIBDIR / name for name in ("libLHAPDF.dylib", "libomp.dylib")]
    sources = [s for s in sources if s.exists()]
    if not sources:
        return
    for sub in (work_subdir / "SubProcesses").glob("P*_*"):
        if not sub.is_dir():
            continue
        for src in sources:
            link = sub / src.name
            if link.exists() or link.is_symlink():
                continue
            link.symlink_to(src)
