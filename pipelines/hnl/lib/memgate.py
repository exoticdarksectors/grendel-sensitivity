"""Memory gate for the fine-grid scan jobs (lib/scan_point.py fine): at most HNL_FINE_SCAN_MAX of them run at once on
this machine (default 3), and one starts only when the memory it will need is free.

A fine job holds four 1M pools (B, B_c, Lambda_b, W/Z): the scan peaks at 3.07 GiB RSS (Umu 3 GeV, `time -v`;
~2.4-2.6 GiB for most of the ray casting), against 0.3-0.7 GiB for a coarse scan and 0.1 GiB for a W/Z job. A job
passes the gate when

    MemAvailable - sum over the running fine jobs of max(0, PEAK - their RSS now)  >=  PEAK + MARGIN

(GiB; PEAK = HNL_FINE_SCAN_PEAK_GB, default 3.2; MARGIN = HNL_FINE_SCAN_MARGIN_GB, default 1.0), so the growth
still to come of a job that has just started is counted. The margin does not cover memory-heavy jobs of other
pipelines that start later: with two fine scans and four BC4/BC5 scans (~1.4 GiB each) swap filled up, so while such
pipelines share the machine raise it (HNL_FINE_SCAN_MARGIN_GB=2 was used then).
Tokens are flocks on work/hnl/locks/fine_mem_token<i>.lock (the holder's PID is written into the file once it has
passed the gate; the lock vanishes with the process); the check runs under work/hnl/locks/fine_mem_start.lock so two
jobs never pass on the same free memory.

    acquire()      in the job, before its first heavy step; returns the open token (keep it until the job ends)
    would_pass()   for the driver, before it queues for a CPU slot (advisory: no token taken)
"""
from __future__ import annotations

import fcntl
import os
import time
from pathlib import Path

LOCKS = Path(os.environ["REPRO"]) / "work/hnl/locks"
GB = 1024 ** 3


def _cfg():
    return (max(1, int(os.environ.get("HNL_FINE_SCAN_MAX", "3"))),
            float(os.environ.get("HNL_FINE_SCAN_PEAK_GB", "3.2")) * GB,
            float(os.environ.get("HNL_FINE_SCAN_MARGIN_GB", "1.0")) * GB)


def mem_available() -> float:
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024.0
    raise RuntimeError("no MemAvailable in /proc/meminfo")


def tree_rss(pid: int) -> float:
    """RSS (bytes) of pid and all its descendants."""
    children: dict[int, list[int]] = {}
    for d in os.listdir("/proc"):
        if not d.isdigit():
            continue
        try:
            with open(f"/proc/{d}/stat") as f:
                ppid = int(f.read().rsplit(")", 1)[1].split()[1])
        except (OSError, ValueError, IndexError):
            continue
        children.setdefault(ppid, []).append(int(d))
    total, todo, page = 0.0, [pid], os.sysconf("SC_PAGE_SIZE")
    while todo:
        p = todo.pop()
        try:
            with open(f"/proc/{p}/statm") as f:
                total += int(f.read().split()[1]) * page
        except (OSError, ValueError, IndexError):
            continue
        todo += children.get(p, [])
    return total


def _running(skip: int | None = None) -> tuple[list[int], list[int]]:
    """(free token indices, PIDs of the jobs that hold a token and have passed the gate)."""
    n, _, _ = _cfg()
    free, pids = [], []
    LOCKS.mkdir(parents=True, exist_ok=True)
    for i in range(n):
        if i == skip:
            continue
        with open(LOCKS / f"fine_mem_token{i}.lock", "a+") as f:
            try:
                fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
                fcntl.flock(f, fcntl.LOCK_UN)
                free.append(i)
                continue
            except BlockingIOError:
                pass
            f.seek(0)
            txt = f.read().strip()
        if txt.isdigit() and os.path.exists(f"/proc/{txt}"):
            pids.append(int(txt))
    return free, pids


def _headroom(pids) -> float:
    _, peak, margin = _cfg()
    return mem_available() - sum(max(0.0, peak - tree_rss(p)) for p in pids) - (peak + margin)


def would_pass() -> bool:
    free, pids = _running()
    return bool(free) and _headroom(pids) >= 0


def acquire(log=print, poll: float = 30.0):
    n, peak, margin = _cfg()
    LOCKS.mkdir(parents=True, exist_ok=True)
    t0, said = time.time(), False
    while True:                                           # a token
        for i in range(n):
            tok = open(LOCKS / f"fine_mem_token{i}.lock", "a+")
            try:
                fcntl.flock(tok, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                tok.close()
                continue
            tok.seek(0)
            tok.truncate()
            tok.flush()
            break
        else:
            if not said:
                log(f"memgate: {n} fine jobs running, waiting for a token")
                said = True
            time.sleep(poll)
            continue
        break
    with open(LOCKS / "fine_mem_start.lock", "a+") as st:  # memory, one starter at a time
        fcntl.flock(st, fcntl.LOCK_EX)
        said = False
        while True:
            _, pids = _running(skip=i)
            h = _headroom(pids)
            if h >= 0:
                break
            if not said:
                log(f"memgate: waiting for memory (MemAvailable {mem_available() / GB:.1f} GB, {len(pids)} other "
                    f"fine job(s), need {(peak + margin) / GB:.1f} GB + their growth; short by {-h / GB:.1f} GB)")
                said = True
            time.sleep(poll)
        tok.seek(0)
        tok.truncate()
        tok.write(f"{os.getpid()}\n")
        tok.flush()
    log(f"memgate: token {i}, waited {time.time() - t0:.0f} s, MemAvailable {mem_available() / GB:.1f} GB")
    return tok
