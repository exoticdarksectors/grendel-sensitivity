"""Reusable MG5 process-runner helpers for HNL production."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from production.madgraph._mg5_common import (
    MG5_EXE,
    PYTHON_EXE,
    force_compile_subprocesses,
    has_five_flavor_proton,
    mg5_subprocess_env,
    patch_me5_configuration,
    patch_rpath_for_lhapdf,
    write_process_block,
)


def ensure_process_dir(
    *,
    label: str,
    model_import: str,
    proc_card: Path,
    work_dir: Path,
    process_dir: Path,
    generation_timeout: int = 600,
    compile_timeout: int = 900,
) -> Path | None:
    """Build or reuse a cached MG5 process directory."""
    if (
        (process_dir / "bin" / "generate_events").exists()
        and has_five_flavor_proton(process_dir)
    ):
        return process_dir

    if process_dir.exists():
        shutil.rmtree(process_dir, ignore_errors=True)

    if not proc_card.exists():
        raise FileNotFoundError(f"Process card not found: {proc_card}")

    work_dir.mkdir(parents=True, exist_ok=True)
    cmd_file = work_dir / f"mg5_gen_{label}.txt"
    log_file = work_dir / f"mg5_gen_{label}.log"
    compile_log = work_dir / f"mg5_compile_{label}.log"

    proc_lines = proc_card.read_text().splitlines(keepends=True)
    with open(cmd_file, "w") as f:
        f.write(f"{model_import}\n\n")
        f.write("set automatic_html_opening False\n")
        write_process_block(f, proc_lines)
        f.write(f"\noutput {process_dir} -nojpeg\n")
        f.write("quit\n")

    with open(log_file, "w") as log:
        result = subprocess.run(
            [str(PYTHON_EXE), str(MG5_EXE), str(cmd_file)],
            stdout=log,
            stderr=subprocess.STDOUT,
            cwd=work_dir,
            timeout=generation_timeout,
            env=mg5_subprocess_env(),
        )

    if result.returncode != 0 or not (process_dir / "bin" / "generate_events").exists():
        print(f"    FAILED: process generation (see {log_file})")
        return None

    if not force_compile_subprocesses(process_dir, compile_log, timeout=compile_timeout):
        print(f"    FAILED: SubProcess pre-compile (see {compile_log})")
        return None

    patch_rpath_for_lhapdf(process_dir)
    cmd_file.unlink(missing_ok=True)
    return process_dir


def write_run_card(work_subdir: Path, cards_dir: Path, n_events: int) -> None:
    """Write run_card.dat into an MG5 process directory."""
    dest_cards = work_subdir / "Cards"
    dest_cards.mkdir(exist_ok=True)

    run_content = (cards_dir / "run_card_template.dat").read_text()
    run_content = run_content.replace("N_EVENTS_PLACEHOLDER", str(n_events))
    (dest_cards / "run_card.dat").write_text(run_content)
    patch_me5_configuration(dest_cards)


def run_events(
    work_subdir: Path,
    run_name: str,
    *,
    nb_core: int = 1,
    timeout: int = 3600,
) -> Path | None:
    """Run MG5 generate_events and return the generated LHE path."""
    log_file = work_subdir / f"generate_events_{run_name}.log"
    cmd = [
        str(PYTHON_EXE),
        "bin/generate_events",
        run_name,
        "-f",
        "--laststep=parton",
    ]
    if nb_core > 1:
        cmd += ["--multicore", f"--nb_core={nb_core}"]
    else:
        cmd += ["--nb_core=1"]

    with open(log_file, "w") as log:
        result = subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            cwd=work_subdir,
            timeout=timeout,
            env=mg5_subprocess_env(),
        )

    if result.returncode != 0:
        print(f"    FAILED: event generation (see {log_file})")
        return None

    run_dir = work_subdir / "Events" / run_name
    for lhe in (run_dir / "unweighted_events.lhe.gz", run_dir / "unweighted_events.lhe"):
        if lhe.exists():
            return lhe
    return None
