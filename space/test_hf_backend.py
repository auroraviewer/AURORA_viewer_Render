"""
Unit tests for the Hugging Face backend of SampleStore.

Network access is stubbed, so these pin the contract rather than the Hub: repo
layout -> sample ids, one sample spanning two files, the exact URL h5py is
handed, prefetch progress accounting, and that unknown samples fail first.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from samples import SampleStore  # noqa: E402

failures = []


def check(label, cond, extra=""):
    ok = bool(cond)
    print(("PASS  " if ok else "FAIL  ") + label + (f"  {extra}" if extra else ""))
    if not ok:
        failures.append(label)


REPO = "AURORAData/AURORA_prediction"
SUBFOLDER = "predictions_7um_28um/TCGA_LUAD"
HIGHRES = "TCGA-44-2668-01A-01-TS1-highres-pred.h5"
COARSE = "TCGA-44-2668-01A-01-TS1-pred.h5"


class FakeFS:
    def __init__(self, paths):
        self.paths = paths
        self.calls = []

    def find(self, root):
        self.root = root
        return [p for p in self.paths if p.startswith(root)]

    def open(self, url, mode="rb", block_size=None, **kw):
        self.calls.append({"url": url, "mode": mode, "block_size": block_size})
        return object()


class FakeH5:
    class File:
        def __init__(self, target, mode="r"):
            self.target = target
            self.mode = mode


def make_store(paths, **kw):
    store = SampleStore(backend="hf", repo_id=REPO, subfolder=SUBFOLDER,
                        block_size=4 << 20, **kw)
    fs = FakeFS(paths)
    store._fs = fs
    return store, fs


# --------------------------------------------------------------------------
# discovery and sample-id grouping
# --------------------------------------------------------------------------
paths = [
    f"datasets/{REPO}/{SUBFOLDER}/{HIGHRES}",
    f"datasets/{REPO}/{SUBFOLDER}/{COARSE}",
    f"datasets/{REPO}/{SUBFOLDER}/TCGA-05-4244-01A-01-BS1-highres-pred.h5",
    f"datasets/{REPO}/{SUBFOLDER}/README.md",
    f"datasets/{REPO}/other_folder/SHOULD_NOT_MATCH-highres-pred.h5",
]
store, fs = make_store(paths)
n = store.refresh()
check("finds prediction files under the subfolder", n == 2, n)
check("repo root is the subfolder", fs.root == f"datasets/{REPO}/{SUBFOLDER}", fs.root)
ids = sorted(s.sample_id for s in store.samples())
check("sample ids strip the schema suffix",
      ids == ["TCGA-05-4244-01A-01-BS1", "TCGA-44-2668-01A-01-TS1"], ids)

s = store.get("TCGA-44-2668-01A-01-TS1")
check("a sample merging two files exposes both",
      sorted(s.paths) == sorted([f"{SUBFOLDER}/{HIGHRES}", f"{SUBFOLDER}/{COARSE}"]),
      s.paths)
check("the highres file is the primary path", s.path == f"{SUBFOLDER}/{HIGHRES}", s.path)
check("kind reports highres", s.kind == "highres", s.kind)
check("backend is hf", s.backend == "hf", s.backend)

# comma-separated subfolders: both folders are searched and the paths stay
# relative to the repo root, so one URL format works for every file
multi = SampleStore(backend="hf", repo_id=REPO,
                    subfolder=f"{SUBFOLDER},other_folder", block_size=4 << 20)
multi._fs = FakeFS(paths)
nm = multi.refresh()
check("comma-separated subfolders search every folder", nm == 3, nm)
m = multi.get("TCGA-44-2668-01A-01-TS1")
check("paths stay repo-relative across folders",
      all(p.startswith(SUBFOLDER + "/") for p in m.paths), m.paths)

single = store.get("TCGA-05-4244-01A-01-BS1")
check("a single-file sample has one path", len(single.paths) == 1, single.paths)

try:
    store.get("NOPE")
    check("unknown sample raises KeyError", False)
except KeyError:
    check("unknown sample raises KeyError", True)

# --------------------------------------------------------------------------
# open() hands h5py one handle per file
# --------------------------------------------------------------------------
import types  # noqa: E402

real_h5py = sys.modules.get("h5py")
sys.modules["h5py"] = types.SimpleNamespace(File=FakeH5.File)
try:
    handles = store.open("TCGA-44-2668-01A-01-TS1")
    check("open returns one handle per file", len(handles) == 2, len(handles))
    check("handles are h5py.File objects", all(isinstance(h, FakeH5.File) for h in handles))
    check("one filesystem open per file", len(fs.calls) == 2, len(fs.calls))
    urls = [c["url"] for c in fs.calls]
    check("both URLs point inside the repo",
          all(u.startswith(f"datasets/{REPO}/{SUBFOLDER}/") for u in urls), urls)
    check("URLs cover both schema files",
          any(u.endswith(HIGHRES) for u in urls) and any(u.endswith(COARSE) for u in urls))
    check("opened in binary read mode", all(c["mode"] == "rb" for c in fs.calls))
    check("block size is forwarded",
          all(c["block_size"] == 4 << 20 for c in fs.calls))
finally:
    if real_h5py is not None:
        sys.modules["h5py"] = real_h5py
    else:
        del sys.modules["h5py"]

# --------------------------------------------------------------------------
# prefetch accounting
# --------------------------------------------------------------------------
store2, fs2 = make_store(paths, local_dir=None)
store2.refresh()
store2.cache_dir = os.path.join(HERE, "_prefetch_test")
import shutil  # noqa: E402

shutil.rmtree(store2.cache_dir, ignore_errors=True)
try:
    from samples import PrefetchJob

    job = PrefetchJob(sample_id="X", total=1000, downloaded=250)
    check("prefetch progress is a fraction", job.progress == 0.25, job.progress)
    check("prefetch reports a state", job.as_dict()["state"] == "pending")

    # a local-backend sample is already "done"
    local = SampleStore(backend="local", local_dir=os.path.dirname(HERE))
    local.refresh()
    real_id = local.samples()[0].sample_id
    st = local.prefetch_status(real_id)
    check("local samples report prefetch done", st["state"] == "done", st)

    # _download writes to the cache dir and tracks bytes
    class FakeResp:
        headers = {"content-length": "10"}

        def __init__(self):
            self.n = 0

        def raise_for_status(self):
            return None

        def iter_content(self, n):
            for _ in range(2):
                yield b"x" * 5

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    import requests  # noqa: E402

    real_get = requests.get
    requests.get = lambda *a, **k: FakeResp()
    try:
        sample = store2.get("TCGA-44-2668-01A-01-TS1")
        job2 = store2._jobs.setdefault(sample.sample_id,
                                       PrefetchJob(sample_id=sample.sample_id))
        path = store2._download(sample, sample.paths[0], job2)
        check("download writes the file into the cache dir",
              os.path.exists(path) and path.startswith(store2.cache_dir), path)
        check("download reports bytes received", job2.downloaded == 10, job2.downloaded)
        # a second call is a no-op that still accounts for the bytes
        before = job2.downloaded
        store2._download(sample, sample.paths[0], job2)
        check("re-downloading an existing file adds its size",
              job2.downloaded > before, (before, job2.downloaded))
    finally:
        requests.get = real_get
finally:
    shutil.rmtree(store2.cache_dir, ignore_errors=True)

# --------------------------------------------------------------------------
# local backend
# --------------------------------------------------------------------------
DATA = os.path.dirname(HERE)
local = SampleStore(backend="local", local_dir=DATA)
check("local backend finds the sample", local.refresh() >= 1, list(local._samples))
check("local backend needs no repo id", local.repo_id is None)

print()
if failures:
    print(f"{len(failures)} FAILURE(S): " + "; ".join(failures))
    sys.exit(1)
print("ALL HF BACKEND CHECKS PASSED")

