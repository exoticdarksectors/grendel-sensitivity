#!/usr/bin/env python3
"""Run a command once enough CPU slots are free on this machine.

All pipelines share one pool of slots (default: number of cores minus one,
override with GRENDEL_SLOTS), so that several benchmarks can queue full-statistics
jobs at the same time without oversubscribing the machine. A slot is an flock on
$REPRO/work/slots/slot<i>.lock; locks vanish with the process, so a killed job
never leaves a slot taken.

    slot_run.py [-n K] [--label NAME] -- command args...

-n K takes K slots at once (all or nothing, so jobs never deadlock). The exit
status is the command's.

Requests are served strictly in arrival order (first come, first served): each
request takes a ticket (work/slots/queue/<ticket>.<pid>, flock-held while it
waits) and only the oldest live ticket may take slots. A job that needs K slots
therefore keeps the slots that free up until it has K, instead of losing them to
single-slot jobs that arrive later. A waiter that is killed leaves a ticket that is
no longer locked; the next waiter removes it.
"""
import argparse
import fcntl
import os
import random
import subprocess
import sys
import time
from pathlib import Path


def slot_dir() -> Path:
    repro = os.environ.get("REPRO")
    if not repro:
        repro = str(Path(__file__).resolve().parents[2])
    d = Path(repro) / "work" / "slots"
    (d / "queue").mkdir(parents=True, exist_ok=True)
    return d


def n_slots() -> int:
    env = os.environ.get("GRENDEL_SLOTS")
    if env:
        return max(1, int(env))
    return max(1, (os.cpu_count() or 2) - 1)


def try_acquire(k: int, total: int, d: Path):
    held = []
    order = list(range(total))
    random.shuffle(order)
    for i in order:
        f = open(d / f"slot{i}.lock", "a+")
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            f.close()
            continue
        held.append(f)
        if len(held) == k:
            return held
    for f in held:
        fcntl.flock(f, fcntl.LOCK_UN)
        f.close()
    return None


def take_ticket(d: Path):
    """Next ticket number (a counter under flock), as a locked file in the queue."""
    with open(d / "ticket.counter", "a+") as c:
        fcntl.flock(c, fcntl.LOCK_EX)
        c.seek(0)
        txt = c.read().strip()
        n = int(txt) + 1 if txt.isdigit() else 1
        c.seek(0)
        c.truncate()
        c.write(f"{n}\n")
        c.flush()
    name = f"{n:012d}.{os.getpid()}"
    tmp = d / "queue" / f".{name}.tmp"
    f = open(tmp, "a+")
    fcntl.flock(f, fcntl.LOCK_EX)
    os.rename(tmp, d / "queue" / name)  # visible only once locked
    return name, f


def at_head(d: Path, name: str) -> bool:
    """True when no live ticket is older than ours (dead tickets are removed)."""
    for other in sorted(p.name for p in (d / "queue").iterdir() if not p.name.startswith(".")):
        if other >= name:
            return True
        path = d / "queue" / other
        try:
            g = open(path, "a+")
        except FileNotFoundError:
            continue
        try:
            fcntl.flock(g, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            g.close()
            return False  # an older request is still waiting
        try:
            path.unlink()  # its process died while waiting
        except FileNotFoundError:
            pass
        g.close()
    return True


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-n", type=int, default=1, help="slots to take (default 1)")
    p.add_argument("--label", default="", help="written into the slot file while held")
    p.add_argument("cmd", nargs=argparse.REMAINDER)
    a = p.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd:
        p.error("no command given")
    total = n_slots()
    k = min(a.n, total)
    d = slot_dir()
    name, ticket = take_ticket(d)
    held = None
    try:
        while not at_head(d, name):
            time.sleep(0.5 + 0.5 * random.random())
        # queue.lock: the head-of-line lock of the previous version of this script, still taken so that
        # requests queued by it (already running processes) and by this version never take slots at once.
        head = open(d / "queue.lock", "a+")
        fcntl.flock(head, fcntl.LOCK_EX)
        try:
            while held is None:
                held = try_acquire(k, total, d)
                if held is None:
                    time.sleep(0.5 + 0.5 * random.random())
        finally:
            fcntl.flock(head, fcntl.LOCK_UN)
            head.close()
    finally:
        try:
            (d / "queue" / name).unlink()
        except FileNotFoundError:
            pass
        fcntl.flock(ticket, fcntl.LOCK_UN)
        ticket.close()
    try:
        for f in held:
            f.seek(0)
            f.truncate()
            f.write(f"{os.getpid()} {a.label or ' '.join(cmd)[:200]}\n")
            f.flush()
        return subprocess.call(cmd)
    finally:
        for f in held:
            f.seek(0)
            f.truncate()
            fcntl.flock(f, fcntl.LOCK_UN)
            f.close()


if __name__ == "__main__":
    sys.exit(main())
