"""
Resolve sample IDs to prediction H5 files and open them.

Two backends:

``local``
    A directory of prediction ``.h5`` files, in either schema.  Used for
    development and for any deployment that mirrors the files on disk.

``hf``
    A Hugging Face Hub dataset repo, e.g. ``AURORAData/AURORA_prediction`` with
    subfolder ``predictions_7um_28um/TCGA_LUAD``.  A sample can either be
    *streamed* through ``HfFileSystem`` (HTTP range requests with an fsspec block
    cache) or *prefetched* to local disk.

Why prefetch matters
--------------------
``exp_pred_14`` is chunked ``(256 channels, 16384 spatial rows)``, so the
*channel* axis is split while every chunk spans the full spatial axis.  Reading
one gene's complete row therefore touches 24 chunks.  Measured on a 291 MB
sample over a normal connection:

    open the file                 2.7 s
    first 7 um gene (24 chunks)  50.8 s   <- unusable interactively
    second 7 um gene              0.2 s   <- most chunks already cached

So streaming is fine for a warm cache but poor for the first view of each new
gene.  ``prefetch`` downloads the whole file once (a few tens of seconds) and
then every view is local.
"""

from __future__ import annotations

import os
import posixpath
import threading
import time
from dataclasses import dataclass, field

H5_SUFFIXES = ("-highres-pred.h5", "-pred.h5", ".h5")


@dataclass
class Sample:
    """
    One biological sample.

    A sample can be spread over more than one file: the 7/28 um level lives in
    ``<id>-highres-pred.h5`` and the 112 um level in ``<id>-pred.h5``, which
    share a sample id.  ``paths`` lists every file contributing to it, ordered so
    the finest level's file comes first.
    """

    sample_id: str
    paths: list[str]
    backend: str                 # "local" | "hf"
    size: int = 0

    @property
    def path(self) -> str:
        """Primary file (the finest-resolution one)."""
        return self.paths[0]

    @property
    def kind(self) -> str:
        return "highres" if any("highres-pred" in p for p in self.paths) else "coarse"


@dataclass
class PrefetchJob:
    sample_id: str
    state: str = "pending"       # pending | running | done | error
    downloaded: int = 0
    total: int = 0
    path: str | None = None
    error: str | None = None
    started: float = field(default_factory=time.time)
    finished: float | None = None

    @property
    def progress(self) -> float:
        if self.total <= 0:
            return 0.0
        return min(1.0, self.downloaded / self.total)

    def as_dict(self) -> dict:
        return {
            "sample_id": self.sample_id,
            "state": self.state,
            "downloaded": self.downloaded,
            "total": self.total,
            "progress": round(self.progress, 4),
            "path": self.path,
            "error": self.error,
            "elapsed": round((self.finished or time.time()) - self.started, 1),
        }


def hf_token():
    return os.environ.get("HF_TOKEN") or os.environ.get("HUGGING_FACE_HUB_TOKEN")


class SampleStore:
    def __init__(
        self,
        backend: str | None = None,
        local_dir: str | None = None,
        repo_id: str | None = None,
        subfolder: str | None = None,
        revision: str | None = None,
        block_size: int | None = None,
    ):
        self.backend = (backend or os.environ.get("AURORA_BACKEND") or "local").lower()
        # Resolve relative paths now, while the working directory is the one the
        # user launched from.  `AURORA_LOCAL_DIR=../` is a common setting and it
        # means something different if the process later changes directory.
        local = local_dir or os.environ.get("AURORA_LOCAL_DIR")
        if local and not os.path.isabs(local):
            local = os.path.abspath(local)
        self.local_dir = local
        self.repo_id = repo_id or os.environ.get("AURORA_HF_REPO")
        # One or more comma-separated folders inside the repo.  The predictions
        # may be split (e.g. the 112 um files in one folder and the 7/28 um files
        # in another), so each is searched recursively and the results merged by
        # sample id.  An empty value searches the whole repo.
        raw = subfolder if subfolder is not None else os.environ.get("AURORA_HF_SUBFOLDER", "")
        self.subfolders = [s.strip("/") for s in (raw or "").split(",") if s.strip("/")]
        self.subfolder = ",".join(self.subfolders)      # kept for reporting
        self.revision = revision or os.environ.get("AURORA_HF_REVISION") or None
        self.block_size = int(block_size or os.environ.get("AURORA_HF_BLOCK") or 8 << 20)
        self.prefetch = os.environ.get("AURORA_PREFETCH") == "1"
        self.cache_dir = os.environ.get("AURORA_CACHE_DIR") or os.path.join(
            os.path.expanduser("~"), ".cache", "aurora_viewer")
        self._samples: dict[str, Sample] = {}
        self._lock = threading.Lock()
        self._fs = None
        self._jobs: dict[str, PrefetchJob] = {}

    # -- helpers -----------------------------------------------------------
    def _filesystem(self):
        if self._fs is None:
            from huggingface_hub import HfFileSystem

            self._fs = HfFileSystem(token=hf_token())
        return self._fs

    def _repo_root(self) -> str:
        if not self.repo_id:
            raise RuntimeError("AURORA_HF_REPO is not set")
        return f"datasets/{self.repo_id}"

    @staticmethod
    def _sample_id(filename: str) -> str:
        for suf in H5_SUFFIXES:
            if filename.endswith(suf):
                return filename[: -len(suf)]
        return posixpath.splitext(filename)[0]

    # -- listing -----------------------------------------------------------
    def refresh(self) -> int:
        return self._refresh_hf() if self.backend == "hf" else self._refresh_local()

    def _refresh_local(self) -> int:
        if not self.local_dir or not os.path.isdir(self.local_dir):
            raise RuntimeError(f"AURORA_LOCAL_DIR is not a directory: {self.local_dir!r}")
        grouped: dict[str, list[str]] = {}
        sizes: dict[str, int] = {}
        for root, _dirs, files in os.walk(self.local_dir):
            for fn in files:
                if not fn.endswith(".h5"):
                    continue
                sid = self._sample_id(fn)
                full = os.path.join(root, fn)
                grouped.setdefault(sid, []).append(full)
                sizes[sid] = sizes.get(sid, 0) + os.path.getsize(full)
        found = {
            sid: Sample(sid, self._order(paths), "local", sizes[sid])
            for sid, paths in grouped.items()
        }
        with self._lock:
            self._samples = found
        return len(found)

    @staticmethod
    def _order(paths: list[str]) -> list[str]:
        """Finest-resolution file first, so `path` is the richest one."""
        def rank(p: str) -> int:
            return 0 if "highres-pred" in p else 1
        return sorted(paths, key=rank)

    def _refresh_hf(self) -> int:
        fs = self._filesystem()
        root = self._repo_root()
        grouped: dict[str, list[str]] = {}
        prefixes = [f"{root}/{s}" for s in self.subfolders] or [root]
        for prefix in prefixes:
            for p in fs.find(prefix):
                if not p.endswith(".h5"):
                    continue
                fn = posixpath.basename(p)
                sid = self._sample_id(fn)
                # stored relative to the REPO root (not the subfolder), so
                # open()/_download() can build one URL regardless of which
                # folder the file came from
                rel = p[len(root):].lstrip("/") if p.startswith(root) else fn
                grouped.setdefault(sid, []).append(rel)
        found = {
            sid: Sample(sid, self._order(paths), "hf")
            for sid, paths in grouped.items()
        }
        with self._lock:
            self._samples = found
        return len(found)

    # -- access ------------------------------------------------------------
    def samples(self) -> list[Sample]:
        if not self._samples:
            self.refresh()
        with self._lock:
            return sorted(self._samples.values(), key=lambda s: s.sample_id)

    def get(self, sample_id: str) -> Sample:
        if not self._samples:
            self.refresh()
        with self._lock:
            s = self._samples.get(sample_id)
        if s is None:
            raise KeyError(sample_id)
        return s

    def open(self, sample_id: str):
        """
        Open every file contributing to a sample.

        Returns a list of open ``h5py.File`` objects (the 7/28 um file and the
        112 um file when both exist); ``AuroraPredictions`` merges their levels.
        """
        import h5py

        s = self.get(sample_id)
        if s.backend == "local":
            return [h5py.File(p, "r") for p in s.paths]

        handles = []
        for rel in s.paths:
            if self.prefetch:
                local = self._local_path(rel)
                if not os.path.exists(local):
                    self._download(s, rel)
                handles.append(h5py.File(local, "r"))
            else:
                fs = self._filesystem()
                url = f"{self._repo_root()}/{rel}"
                handles.append(h5py.File(fs.open(url, "rb", block_size=self.block_size), "r"))
        return handles

    # -- prefetch ----------------------------------------------------------
    def _local_path(self, rel: str) -> str:
        return os.path.join(self.cache_dir, os.path.basename(rel))

    def prefetch_status(self, sample_id: str) -> dict:
        s = self.get(sample_id)
        if s.backend == "local":
            return {"state": "done", "local": True, "path": s.path,
                    "total": s.size, "downloaded": s.size, "progress": 1.0}

        job = self._jobs.get(sample_id)
        local = self._local_path(s.path)
        if os.path.exists(local) and job is None:
            size = os.path.getsize(local)
            return {"state": "done", "local": True, "path": local,
                    "total": size, "downloaded": size, "progress": 1.0}
        if job is None:
            job = PrefetchJob(sample_id=sample_id, total=s.size)
            self._jobs[sample_id] = job
            threading.Thread(target=self._download_all, args=(s,), daemon=True).start()
        return job.as_dict()

    def _download_all(self, s: Sample):
        """Download every file of a sample, tracking aggregate progress."""
        job = self._jobs.setdefault(s.sample_id,
                                    PrefetchJob(sample_id=s.sample_id, total=s.size))
        job.state = "running"
        try:
            for rel in s.paths:
                self._download(s, rel, job)
            job.state = "done"
            job.path = self._local_path(s.path)
        except Exception as exc:
            job.state = "error"
            job.error = f"{type(exc).__name__}: {exc}"
        finally:
            job.finished = time.time()

    def _download(self, s: Sample, rel: str, job: PrefetchJob | None = None) -> str:
        import requests

        if job is None:
            job = self._jobs.setdefault(s.sample_id,
                                        PrefetchJob(sample_id=s.sample_id, total=s.size))
        os.makedirs(self.cache_dir, exist_ok=True)
        target = self._local_path(rel)
        if os.path.exists(target):
            job.downloaded += os.path.getsize(target)
            return target
        tmp = target + ".part"
        url = f"https://huggingface.co/datasets/{self.repo_id}/resolve/main/{rel}"
        headers = {}
        token = hf_token()
        if token:
            headers["Authorization"] = f"Bearer {token}"
        try:
            with requests.get(url, stream=True, timeout=180, headers=headers) as r:
                r.raise_for_status()
                if not job.total:
                    job.total = int(r.headers.get("content-length") or 0)
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(1 << 20):
                        fh.write(chunk)
                        job.downloaded += len(chunk)
            os.replace(tmp, target)
            return target
        except Exception:
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
            except Exception:
                pass
            raise
