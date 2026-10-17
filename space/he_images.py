"""
Locate the H&E image that belongs to a sample.

The prediction files do not ship with images, and the repo the old viewer's
README pointed at (``AURORAData/prediction_plots``) is gone.  So the image source
is configurable and resolved by naming convention:

``AURORA_HE_DIR``
    A directory scanned recursively for ``<sample_id>.png|jpg|jpeg|tif|tiff|webp``
    (also ``<sample_id>_he.*``, ``<sample_id>-he.*``).

``AURORA_HE_REPO`` / ``AURORA_HE_SUBFOLDER``
    A Hugging Face **dataset** repo to fetch images from, same naming rules.

``AURORA_HE_PLACEHOLDER``
    A single image used for every sample.  Intended for development and for
    demoing the side-by-side view before the real images are available; the UI
    labels it as a placeholder so it is never mistaken for the real slide.

Nothing here resamples the image into the prediction coordinate space 閳?the
frontend places it using the layer extent, so the two stay aligned by
construction.  ``viewer_core.align_report`` can verify the pairing.
"""

from __future__ import annotations

import io
import json
import os
import posixpath
import re
import threading
from collections import OrderedDict

from PIL import Image
from fastapi import Response

# Rescaled whole-slide H&E images are routinely hundreds of megapixels (the
# workspace sample is 58152 x 9960 = 579 MP), which trips Pillow's
# decompression-bomb guard.  These are first-party, trusted images, so the guard
# is raised deliberately rather than worked around.
Image.MAX_IMAGE_PIXELS = None

IMAGE_EXTS = (".png", ".jpg", ".jpeg", ".tif", ".tiff", ".webp", ".bmp")

# Whole-slide formats (.svs, .ndpi, .mrxs, ...) are deliberately NOT listed:
# Pillow cannot decode them and the predictions were made against the *rescaled*
# 0.5 um/px image anyway, so a rescaled TIFF/PNG is what is wanted here.
NAME_PATTERNS = ("{id}", "{id}_he", "{id}-he", "{id}_HE", "{id}_he_image")

# Suffixes appended to the sample id by the AURORA preprocessing.  `_rescaled`
# is the 0.5 um/px image the predictions were actually made against, so it is
# the one we want.
ID_SUFFIXES = (
    "_rescaled", "-rescaled", "_rescaled_image",
    "_he", "-he", "_HE", "_he_image",
)

# Directories that hold viewer chrome rather than slide images.
SKIP_DIRS = {"images", "assets", "static", "icons", "templates"}

# An H&E file's basename should look like a specimen id.  This keeps the index
# from picking up things like ``zoomout_rest.png`` from a slider skin.
SAMPLE_RE = re.compile(r"^(TCGA|TARGET|BWH|STANFORD|[A-Z]{2,6}[-_])", re.IGNORECASE)


class HEImageStore:
    def __init__(self, directory=None, repo_id=None, subfolder=None,
                 placeholder=None, revision=None):
        self.directory = directory or os.environ.get("AURORA_HE_DIR")
        self.repo_id = repo_id or os.environ.get("AURORA_HE_REPO")
        self.subfolder = (subfolder or os.environ.get("AURORA_HE_SUBFOLDER") or "").strip("/")
        self.placeholder = placeholder or os.environ.get("AURORA_HE_PLACEHOLDER")
        self.revision = revision or os.environ.get("AURORA_HE_REVISION") or None
        self._index: dict[str, str] | None = None
        self._lock = threading.Lock()
        # encoded, downscaled copies keyed by (path, max_width)
        self._photo_cache: OrderedDict = OrderedDict()
        self._photo_cache_max = 8
        # decoded source images keyed by path (decoding is the slow part)
        self._source_cache: OrderedDict = OrderedDict()
        self._source_cache_max = 3
        self._photo_lock = threading.Lock()

    # -- discovery ---------------------------------------------------------
    def _build_index(self) -> dict[str, str]:
        index: dict[str, str] = {}
        if not self.directory or not os.path.isdir(self.directory):
            return index
        for root, dirs, files in os.walk(self.directory):
            # prune viewer-chrome folders in place so os.walk does not descend
            dirs[:] = [d for d in dirs if d.lower() not in SKIP_DIRS]
            for fn in files:
                stem, ext = os.path.splitext(fn)
                if ext.lower() not in IMAGE_EXTS:
                    continue
                key = stem
                for pat in ID_SUFFIXES:
                    if key.endswith(pat):
                        key = key[: -len(pat)]
                        break
                if not SAMPLE_RE.match(key):
                    continue
                index.setdefault(key, os.path.join(root, fn))
        return index

    def index(self) -> dict[str, str]:
        with self._lock:
            if self._index is None:
                self._index = self._build_index()
            return self._index

    def refresh(self):
        with self._lock:
            self._index = None
        return len(self.index())

    # -- lookup ------------------------------------------------------------
    def path_for(self, sample_id: str) -> str | None:
        hit = self.index().get(sample_id)
        if hit:
            return hit
        if self.placeholder and os.path.exists(self.placeholder):
            return self.placeholder
        if self.repo_id:
            # fall through to the Hub on demand
            return self._download(sample_id)
        return None

    def is_placeholder(self, sample_id: str) -> bool:
        return self.index().get(sample_id) is None and bool(self.placeholder)

    def _download(self, sample_id: str) -> str | None:
        try:
            from huggingface_hub import hf_hub_download

            for pat in NAME_PATTERNS:
                for ext in IMAGE_EXTS:
                    rel = "/".join(x for x in (self.subfolder, pat.format(id=sample_id) + ext) if x)
                    try:
                        return hf_hub_download(
                            repo_id=self.repo_id, filename=rel, repo_type="dataset",
                            revision=self.revision,
                            token=os.environ.get("HF_TOKEN"),
                        )
                    except Exception:
                        continue
        except Exception:
            return None
        return None

    def describe(self, sample_id: str) -> dict:
        path = self.path_for(sample_id)
        if path is None:
            # a pyramid in the HF repo counts as available even though there is
            # no single flat file to open
            pyr = self.hf_pyramid(sample_id) if self.repo_id else None
            if pyr:
                return {
                    "available": True,
                    "placeholder": False,
                    "basename": f"{sample_id}_pyramid",
                    "url": f"/api/he/{sample_id}",
                    "image_size": pyr["native"],
                    "pyramid": True,
                }
            # same for a local pyramid: read its manifest (or infer from the
            # level name) so a pyramid alone is enough - no giant flat file needed
            ldir = self.pyramid_dir(sample_id)
            if ldir:
                native = self._pyramid_native(ldir)
                return {
                    "available": True,
                    "placeholder": False,
                    "basename": os.path.basename(ldir),
                    "url": f"/api/he/{sample_id}",
                    "image_size": native,
                    "pyramid": True,
                }
            return {"available": False,
                    "hint": "set AURORA_HE_DIR (or AURORA_HE_PLACEHOLDER) to enable "
                            "the side-by-side H&E view"}
        # The served image is downscaled, so report the NATIVE pixel size: the
        # client needs it to place the image in the same coordinate space as the
        # predictions, whose coordinates are in native image pixels.
        size = None
        try:
            with Image.open(path) as im:
                size = [int(im.width), int(im.height)]
        except Exception:
            size = None
        return {
            "available": True,
            "placeholder": self.is_placeholder(sample_id),
            "basename": os.path.basename(path),
            "url": f"/api/he/{sample_id}",
            "image_size": size,
        }

    def status(self) -> dict:
        return {
            "directory": self.directory,
            "repo": self.repo_id,
            "subfolder": self.subfolder,
            "placeholder": bool(self.placeholder),
            "indexed": len(self.index()),
        }

    # -- serving -----------------------------------------------------------
    def pyramid_dir(self, sample_id: str) -> str | None:
        """
        A pre-built pyramid for this sample, if one is present.

        Building a pyramid is a one-off job (`build_he_pyramid.py`); serving it is
        then a file copy.  Without one the server has to decode the whole slide on
        first request, which is ~13 s for a 579 MP TIFF.
        """
        for base in (self.directory, os.environ.get("AURORA_HE_PYRAMID_DIR")):
            if not base:
                continue
            for cand in (os.path.join(base, sample_id),
                         os.path.join(base, "he", sample_id),
                         os.path.join(base, f"{sample_id}_pyramid")):
                if os.path.isdir(cand) and self._pyramid_levels(cand):
                    return cand
        return None

    @staticmethod
    def _pyramid_levels(directory: str) -> list[tuple[int, str]]:
        """[(width, path)] for the level images in a pyramid directory, ascending."""
        out = []
        try:
            names = os.listdir(directory)
        except OSError:
            return out
        for fn in names:
            stem, ext = os.path.splitext(fn)
            if ext.lower() not in (".webp", ".jpg", ".jpeg", ".png"):
                continue
            if not stem.startswith("L"):
                continue
            p = os.path.join(directory, fn)
            try:
                with Image.open(p) as im:
                    out.append((im.width, p))
            except Exception:
                continue
        out.sort()
        return out

    def _pyramid_native(self, directory: str) -> list[int] | None:
        """
        Native image size in prediction pixels, from a pyramid directory.

        Prefers ``manifest.json`` (``{"size": [width, height]}``); otherwise it
        infers from the finest level's name: ``L<k>`` is the native image halved
        ``k`` times, so native = level * 2**k.  The batch converter writes the
        manifest, making the inference a fallback only.
        """
        mf = os.path.join(directory, "manifest.json")
        if os.path.exists(mf):
            try:
                with open(mf, "r", encoding="utf-8") as fh:
                    size = json.load(fh).get("size")
                if isinstance(size, (list, tuple)) and len(size) == 2:
                    return [int(size[0]), int(size[1])]
            except Exception:
                pass
        levels = self._pyramid_levels(directory)
        if not levels:
            return None
        w, path = levels[-1]                     # widest = finest = least halved
        m = re.match(r"L(\d+)", os.path.basename(path))
        k = int(m.group(1)) if m else 0
        try:
            with Image.open(path) as im:
                h = im.height
        except Exception:
            h = 0
        return [w * (2 ** k), h * (2 ** k)]

    def serve_pyramid(self, directory: str, max_width: int) -> Response | None:
        """
        Serve from a pre-built pyramid: the smallest level that still covers
        `max_width`, downscaled only if it is much larger than needed.
        """
        levels = self._pyramid_levels(directory)
        if not levels:
            return None
        pick = None
        for w, p in levels:            # ascending, so the first fit is smallest
            if w >= max_width:
                pick = (w, p)
                break
        if pick is None:
            pick = levels[-1]          # nothing big enough; use the largest
        w, path = pick
        if w > max_width * 2:
            # still far too big - fall back to resizing it
            return self.serve(path, max_width)
        try:
            with open(path, "rb") as fh:
                data = fh.read()
        except OSError:
            return None
        mime = "image/webp" if path.lower().endswith(".webp") else (
            "image/jpeg" if path.lower().endswith((".jpg", ".jpeg")) else "image/png")
        return Response(content=data, media_type=mime,
                        headers={"Cache-Control": "public, max-age=86400",
                                 "X-Cache": "pyramid", "X-Level-Width": str(w)})

    def _hf_fs(self):
        from huggingface_hub import HfFileSystem

        return HfFileSystem(token=os.environ.get("HF_TOKEN")
                            or os.environ.get("HUGGING_FACE_HUB_TOKEN"))

    def _hf_root(self) -> str:
        return f"datasets/{self.repo_id}" + (f"/{self.subfolder}" if self.subfolder else "")

    # -- HF pyramid ---------------------------------------------------------
    def hf_pyramid(self, sample_id: str) -> dict | None:
        """
        A pre-built pyramid for this sample, stored in the HF repo.

        Returns ``{"levels": [(width, height, rel_path), ...], "native": [w, h]}``
        or None.  Levels are discovered in ``<subfolder>/<id>_pyramid`` (or
        ``<id>``/``he/<id>``), and the native image size is read from a
        ``manifest.json`` when present, else inferred from the level name:
        ``L<k>`` is the native image halved ``k`` times, so native = level * 2**k.
        """
        if not self.repo_id:
            return None
        try:
            fs = self._hf_fs()
        except Exception:
            return None
        root = f"datasets/{self.repo_id}"
        for cand in (f"{self.subfolder}/{sample_id}_pyramid",
                     f"{self.subfolder}/he/{sample_id}_pyramid",
                     f"{self.subfolder}/{sample_id}"):
            url = "/".join(x for x in (root, cand) if x)
            try:
                entries = fs.ls(url)
                # fsspec/HfFileSystem ls returns FileInfo dicts (or strings)
                names = [posixpath.basename(e.get("name") if isinstance(e, dict) else e)
                         for e in entries]
            except Exception:
                continue
            if not any(n.endswith((".webp", ".jpg", ".jpeg", ".png")) and n.startswith("L")
                       for n in names):
                continue
            levels = []
            native = None
            # manifest.json carries the exact native size when the builder wrote it
            if "manifest.json" in names:
                try:
                    with fs.open(f"{url}/manifest.json", "rb") as fh:
                        native = json.loads(fh.read())["size"]
                except Exception:
                    native = None
            for fn in sorted(names):
                stem, ext = posixpath.splitext(fn)
                if ext.lower() not in (".webp", ".jpg", ".jpeg", ".png") or not stem.startswith("L"):
                    continue
                rel = f"{cand}/{fn}" if cand else fn
                try:
                    with fs.open(f"{root}/{rel}", "rb") as fh:
                        with Image.open(fh) as im:
                            w, h = im.width, im.height
                    levels.append((w, h, rel))
                except Exception:
                    continue
            if not levels:
                continue
            levels.sort()
            if native is None:
                # infer from the finest (largest) level's name
                m = re.match(r"L(\d+)", posixpath.basename(levels[-1][2]))
                scale = 2 ** int(m.group(1)) if m else 1
                w, h = levels[-1][0], levels[-1][1]
                native = [w * scale, h * scale]
            return {"levels": levels, "native": native}
        return None

    def _hf_fetch(self, rel: str) -> bytes:
        fs = self._hf_fs()
        with fs.open(f"datasets/{self.repo_id}/{rel}", "rb") as fh:
            return fh.read()

    def serve_hf_pyramid(self, sample_id: str, max_width: int) -> Response | None:
        info = self.hf_pyramid(sample_id)
        if not info:
            return None
        levels = info["levels"]
        pick = None
        for w, h, rel in levels:            # ascending; first that covers is smallest
            if w >= max_width:
                pick = (w, h, rel)
                break
        if pick is None:
            pick = levels[-1]
        w, h, rel = pick
        try:
            data = self._hf_fetch(rel)
        except Exception:
            return None
        if w > max_width * 2:
            # still far bigger than asked; downscale and re-encode once
            try:
                img = Image.open(io.BytesIO(data))
                img = img.resize((max_width, max(1, int(img.height * max_width / img.width))),
                                 Image.LANCZOS)
                buf = io.BytesIO()
                img.convert("RGB").save(buf, format="JPEG", quality=88, optimize=True)
                data = buf.getvalue()
                mime = "image/jpeg"
            except Exception:
                return None
        else:
            mime = "image/webp" if rel.lower().endswith(".webp") else (
                "image/jpeg" if rel.lower().endswith((".jpg", ".jpeg")) else "image/png")
        return Response(content=data, media_type=mime,
                        headers={"Cache-Control": "public, max-age=86400",
                                 "X-Cache": "pyramid", "X-Level-Width": str(w)})

    def serve(self, path: str, max_width: int = 4096) -> Response:
        """
        Serve an H&E image, downscaled and re-encoded for the web.

        H&E is photographic, so JPEG is used for colour images: a 4096 px PNG of
        a 40 MB source weighs ~25 MB and is unusable in a browser, while JPEG
        lands around 1 MB with no visible loss at viewing scale.  Greyscale
        sources (masks, single-channel TIFFs) stay PNG.

        The encoded result is cached: decoding a 40 MB source PNG takes ~20 s.
        """
        key = (path, int(max_width))
        with self._photo_lock:
            hit = self._photo_cache.get(key)
            if hit is not None:
                self._photo_cache.move_to_end(key)
                data, mime = hit
                return Response(content=data, media_type=mime,
                                headers={"Cache-Control": "public, max-age=86400",
                                         "X-Cache": "hit"})

        try:
            img = self._open_cached(path)
        except Exception as exc:
            return Response(status_code=500, content=f"cannot read image: {exc}")

        if img.width > max_width:
            ratio = max_width / img.width
            img = img.resize((max_width, max(1, int(img.height * ratio))),
                             Image.LANCZOS)

        buf = io.BytesIO()
        if img.mode in ("L", "1"):
            img.save(buf, format="PNG", compress_level=3)
            mime = "image/png"
        else:
            img = img.convert("RGB")
            img.save(buf, format="JPEG", quality=88, optimize=True, progressive=True)
            mime = "image/jpeg"

        data = buf.getvalue()
        with self._photo_lock:
            self._photo_cache[key] = (data, mime)
            while len(self._photo_cache) > self._photo_cache_max:
                self._photo_cache.popitem(last=False)
        return Response(content=data, media_type=mime,
                        headers={"Cache-Control": "public, max-age=86400",
                                 "X-Cache": "miss"})

    def _open_cached(self, path: str) -> Image.Image:
        """
        Decode a source image, memoised.

        Decoding a 40 MB PNG or a rescaled TIFF takes ~20 s, which would
        otherwise be paid on every sample switch.
        """
        with self._photo_lock:
            hit = self._source_cache.get(path)
            if hit is not None:
                self._source_cache.move_to_end(path)
                return hit
        img = Image.open(path)
        img.load()
        if img.mode not in ("RGB", "RGBA", "L", "1"):
            img = img.convert("RGB")
        with self._photo_lock:
            self._source_cache[path] = img
            while len(self._source_cache) > self._source_cache_max:
                self._source_cache.popitem(last=False)
        return img


