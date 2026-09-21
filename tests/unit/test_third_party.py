"""third_party/fetch.py: the pins are well formed, and the fetch and verify
paths work against local git repositories and archives."""
import hashlib
import importlib.util
import json
import subprocess
import tarfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
FETCH = ROOT / "third_party" / "fetch.py"


def _load():
    spec = importlib.util.spec_from_file_location("grendel_fetch", FETCH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


fetch = _load()


def test_pins_are_well_formed():
    pins = fetch.load_pins()
    assert {"exHad", "HNLCalc", "FairShip", "SM_HeavyN_CKM_AllMasses_LO", "fonll", "MG5_aMC_v3_6_6",
            "pythia8315", "pythia8317"} <= set(pins)
    dests = set()
    for name, pin in pins.items():
        assert pin["kind"] in ("git", "archive"), name
        assert pin["url"].startswith("https://"), name
        assert pin["what"] and pin["license"], name
        if pin["kind"] == "git":
            assert len(pin["commit"]) == 40 and int(pin["commit"], 16) >= 0, name
        else:
            assert len(pin["sha256"]) == 64 and pin["archive_dir"], name
        dest = pin.get("dest", name)
        assert dest not in dests, f"duplicate destination {dest}"
        dests.add(dest)


def test_fonll_pin_matches_the_grid_generator_pins():
    fonll = fetch.load_pins()["fonll"]
    patches = json.loads((ROOT / "grendel" / "production" / "fonll_grids" / "patches" / "PINS.json").read_text())
    assert fonll["commit"] == patches["fonll"]["mirror_commit"]


def test_list_and_verify_on_an_empty_tree(tmp_path):
    assert subprocess.run(["python3", str(FETCH), "--list"], capture_output=True, text=True).returncode == 0
    result = subprocess.run(["python3", str(FETCH), "--verify", "--dest", str(tmp_path)], capture_output=True,
                            text=True)
    assert result.returncode == 1                       # everything missing, nothing mismatched
    assert result.stdout.count("missing") == len(fetch.load_pins())


def _git(*args, cwd):
    return subprocess.run(["git", *args], cwd=cwd, check=True, capture_output=True, text=True).stdout.strip()


@pytest.fixture
def local_repo(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    _git("init", "-q", cwd=src)
    _git("config", "user.email", "t@example.org", cwd=src)
    _git("config", "user.name", "t", cwd=src)
    _git("config", "uploadpack.allowAnySHA1InWant", "true", cwd=src)
    _git("config", "uploadpack.allowFilter", "true", cwd=src)
    (src / "python").mkdir()
    (src / "python" / "hnl.py").write_text("x = 1\n")
    (src / "big").mkdir()
    (src / "big" / "blob.bin").write_bytes(b"\0" * 1000)
    _git("add", "-A", cwd=src)
    _git("commit", "-q", "-m", "one", cwd=src)
    first = _git("rev-parse", "HEAD", cwd=src)
    (src / "python" / "hnl.py").write_text("x = 2\n")
    _git("commit", "-q", "-am", "two", cwd=src)
    return src, first


def test_fetch_git_checks_out_the_pinned_commit(local_repo, tmp_path):
    src, first = local_repo
    pin = {"kind": "git", "url": str(src), "commit": first}
    dest = fetch.fetch("thing", pin, tmp_path / "out", force=False)
    assert (dest / "python" / "hnl.py").read_text() == "x = 1\n"
    assert (dest / "UPSTREAM_COMMIT.txt").read_text().strip() == first
    assert fetch.verify("thing", pin, tmp_path / "out") == ("ok", f"{dest} @ {first[:12]}")
    # A second call leaves it alone; a changed pin is reported.
    fetch.fetch("thing", pin, tmp_path / "out", force=False)
    assert fetch.verify("thing", {**pin, "commit": "0" * 40}, tmp_path / "out")[0] == "mismatch"


def test_fetch_git_sparse_takes_only_the_listed_paths(local_repo, tmp_path):
    src, first = local_repo
    pin = {"kind": "git", "url": str(src), "commit": first, "sparse": ["python"]}
    dest = fetch.fetch("sparse", pin, tmp_path / "out", force=False)
    assert (dest / "python" / "hnl.py").is_file()
    assert not (dest / "big").exists()


def test_fetch_archive_checks_the_hash_and_unpacks(tmp_path):
    payload = tmp_path / "SM_Model"
    payload.mkdir()
    (payload / "__init__.py").write_text("# model\n")
    archive = tmp_path / "model.tgz"
    with tarfile.open(archive, "w:gz") as tar:
        tar.add(payload, arcname="SM_Model")
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    pin = {"kind": "archive", "url": archive.as_uri(), "sha256": digest, "archive_dir": "SM_Model"}
    dest = fetch.fetch("SM_Model", pin, tmp_path / "out", force=False)
    assert (dest / "__init__.py").read_text() == "# model\n"
    assert fetch.verify("SM_Model", pin, tmp_path / "out")[0] == "ok"
    with pytest.raises(RuntimeError, match="sha256"):
        fetch.fetch("bad", {**pin, "sha256": "0" * 64}, tmp_path / "out2", force=False)
