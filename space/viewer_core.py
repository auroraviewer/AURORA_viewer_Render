"""
Framework-independent core of the viewer service.

Decides *what* gets served; ``app.py`` is a thin HTTP wrapper.  Everything here
is testable without FastAPI.

The renderer produces **two** rasters per view:

* ``values``   - uint8 expression (or cell-type proportion) code, 0 where no
                 patch exists.  The browser maps this through the colormap, so
                 an expression of 0 is drawn as the colormap's low end, *not* as
                 background.
* ``presence`` - uint8, 255 where a patch exists, 0 otherwise.  This is what
                 makes "no patch" transparent in the browser.

Values are drawn exactly as stored.  Nothing is masked, hatched or annotated:
whatever the file contains is what appears.

Guarding against a real bug this fixes: the earlier version multiplied values by
the presence mask, which made "expression 0" and "no data" indistinguishable.
"""

from __future__ import annotations

import io
import os
import threading
import time
from collections import OrderedDict

import numpy as np
from PIL import Image

# whole-slide H&E images exceed Pillow's default bomb guard; see he_images.py
Image.MAX_IMAGE_PIXELS = None

from aurora_h5 import CELL_MICRONS, AuroraPredictions

# Target width for a served raster.  The native 7 um grid is 1040 cells wide and
# a 1600 px viewport shows the whole slide at ~1.5 screen px per cell, so 1x is
# already generous; 2x keeps it crisp under zoom.
TARGET_WIDTH = 2080
MAX_ZOOM = 8
PNG_COMPRESS_LEVEL = 1


def decimate_to(a: np.ndarray, target_w: int) -> np.ndarray:
    """Integer max-pooling to at most `target_w`, preserving hot spots."""
    h, w = a.shape[:2]
    if w <= target_w:
        return a
    k = int(np.ceil(w / target_w))
    h2, w2 = h // k * k, w // k * k
    if h2 == 0 or w2 == 0:
        return a
    if a.ndim == 2:
        return a[:h2, :w2].reshape(h2 // k, k, w2 // k, k).max(axis=(1, 3))
    return a[:h2, :w2].reshape(h2 // k, k, w2 // k, k, a.shape[2]).max(axis=(1, 3))


def _png_l(arr: np.ndarray) -> bytes:
    buf = io.BytesIO()
    Image.fromarray(np.ascontiguousarray(arr, dtype=np.uint8), mode="L").save(
        buf, format="PNG", compress_level=PNG_COMPRESS_LEVEL
    )
    return buf.getvalue()


def _png_rgb(values: np.ndarray, presence: np.ndarray) -> bytes:
    """
    Two channels of information in an RGB PNG:

    ``R`` value     the expression (or proportion) code, 0-255
    ``B`` presence  255 where a patch exists, 0 where there is no tissue

    The presence channel is what makes "no patch" genuinely transparent.  Without
    it the client cannot tell background from a legitimate expression of 0, and
    the raster is drawn as a rectangle of near-black over the whole image.

    Values are shown exactly as stored.  No cell is masked, hatched or otherwise
    annotated: whatever the file contains is what gets drawn.
    """
    rgb = np.zeros(values.shape + (3,), dtype=np.uint8)
    rgb[..., 0] = values
    rgb[..., 2] = presence
    buf = io.BytesIO()
    Image.fromarray(rgb, mode="RGB").save(buf, format="PNG",
                                          compress_level=PNG_COMPRESS_LEVEL)
    return buf.getvalue()


def stretch_present_to_codes(unit: np.ndarray, present: np.ndarray
                             ) -> tuple[np.ndarray, float, float]:
    """
    Map raw probabilities onto colour codes over their **observed range**.

    Cell-type proportions are softmax outputs, so the nine values at every spot
    sum to 1 and each type sits somewhere near 1/9.  On the workspace sample all
    nine lie between about 0.093 and 0.243 - only the bottom quarter of a 0-1
    scale - so mapping them straight onto the colormap wastes three quarters of
    the ramp and every map renders as near-uniform blue.

    Stretching each cell type over the range it actually occupies uses the whole
    ramp, and the value behind the colour is still the raw probability.  The
    bounds are returned so the legend can quote the real numbers.

    Cells with no patch hold 0 and are excluded from the bounds, then forced back
    to code 0, so the background stays "no data" instead of dragging the lower
    bound down and squashing the map.
    """
    if unit.size == 0 or not present.any():
        return np.zeros(unit.shape, dtype=np.uint8), 0.0, 0.0
    vals = unit[present]
    lo, hi = float(vals.min()), float(vals.max())
    if hi <= lo:
        # a flat map: sit in the middle rather than dividing by zero
        out = np.full(unit.shape, 128, dtype=np.uint8)
    else:
        out = np.rint((unit.astype(np.float32) - lo) * (255.0 / (hi - lo))
                      ).clip(0, 255).astype(np.uint8)
    out = np.where(present, out, 0).astype(np.uint8)
    return out, lo, hi


def render_view(layer, kind: str, index: int, target_w: int = TARGET_WIDTH):
    """
    Rasterise one gene or cell type at one resolution.

    Returns ``(png, info)``.  ``kind`` is ``"gene"`` or ``"celltype"``.
    """
    t0 = time.time()
    flat, W, H, _xs, _ys = layer.grid()

    raw_lo = raw_hi = None
    if kind == "gene":
        codes = layer.exp_codes([index])[:, 0]
        raw_unit = "expression (dequantised)"
        raw_max = None
    elif kind == "celltype":
        # keep the float probability.  Round-tripping through the 8-bit code
        # loses the low end (code 24/255 = 0.094, which would read as 0.0).
        props = layer.proportions([index])[:, 0]
        codes = np.rint(props * 255.0).clip(0, 255).astype(np.uint8)
        unit_flat = props.astype(np.float32)
        raw_unit = "fraction of cells"
        raw_max = float(props.max()) if props.size else 0.0
    else:
        raise ValueError(f"unknown kind {kind!r}")

    values = np.zeros(W * H, dtype=np.uint8)
    values[flat] = codes
    values = values.reshape(H, W)

    unit = None
    if kind == "celltype":
        unit = np.zeros(W * H, dtype=np.float32)
        unit[flat] = unit_flat
        unit = unit.reshape(H, W)

    # presence: 255 where a patch exists, 0 where there is no tissue.  This is
    # what makes the background of the prediction transparent rather than being
    # confused with a legitimate expression of 0.
    presence = np.zeros(W * H, dtype=np.uint8)
    presence[flat] = 255
    presence = presence.reshape(H, W)

    if target_w and W > target_w:
        # Integer max-pooling keeps the aspect ratio exact, so the decimated
        # raster still maps linearly onto the layer's extent.  Cell types use max
        # pooling too, which for a proportion means "the most abundant of the
        # pooled cells" - the right summary for a peak-preserving display.
        values_d = decimate_to(values, target_w)
        presence_d = decimate_to(presence, target_w)
    else:
        values_d, presence_d = values, presence

    # Cell types are stretched over the range they actually occupy.  Done on the
    # decimated raster so the reported bounds match exactly what was drawn.
    # Genes already fill the colormap because exp_max is their quantisation max.
    if kind == "celltype":
        unit_d = (decimate_to(unit, target_w) if (target_w and W > target_w)
                  else unit)
        values_d, raw_lo, raw_hi = stretch_present_to_codes(unit_d, presence_d > 0)
        scale = 255.0 / (raw_hi - raw_lo) if raw_hi > raw_lo else 1.0
    else:
        scale = 1.0

    t_raster = time.time() - t0

    t0 = time.time()
    png = _png_rgb(values_d, presence_d)
    mode = "RGB"
    t_png = time.time() - t0

    # The prediction does NOT fill the H&E frame: its own extent is narrower, and
    # it may start at a negative coordinate.  The client needs both so it can
    # place the raster rather than stretch it over the whole image.  Take the
    # extent from the layer, which measures it on the raster's coordinate axis -
    # computing it as grid x cell would close up any hole in the patch lattice.
    origin = list(layer.origin())
    extent = list(layer.extent())

    return png, {
        "kind": kind,
        "index": index,
        "cell_px": layer.cell_px,
        "microns": layer.microns,
        "native_grid": [W, H],
        "raster": [values_d.shape[1], values_d.shape[0]],
        "origin": origin,            # world-pixel position of raster cell (0, 0)
        "extent": extent,            # world-pixel size of the predicted area
        "world_bounds": [origin[0], origin[1],
                         origin[0] + extent[0], origin[1] + extent[1]],
        "mode": mode,
        # so the legend can state what the colour actually means
        "raw_unit": raw_unit,
        "raw_max": raw_max,
        "raw_min": raw_lo,
        "raw_lo": raw_lo,
        "raw_hi": raw_hi,
        "normalized": bool(kind == "celltype"),
        "scale": round(scale, 6),
        "bytes": len(png),
        "raster_seconds": round(t_raster, 4),
        "png_seconds": round(t_png, 4),
    }


class ViewerCache:
    """One open prediction file plus cached rasters and encoded PNGs."""

    def __init__(self, store=None, max_pngs: int | None = None, max_meta: int = 24):
        self.store = store
        self._lock = threading.RLock()
        self._sample_id: str | None = None
        self._h5: AuroraPredictions | None = None
        self._pngs: OrderedDict = OrderedDict()
        # small on purpose: PNGs are ~200-500 KB each and a 512 MB host has no
        # room for 64 of them; exact repeats still hit the cache
        self._max_pngs = int(max_pngs if max_pngs is not None
                             else os.environ.get("AURORA_PNG_CACHE", "16"))
        self._meta: OrderedDict = OrderedDict()
        self._max_meta = max_meta
        self.stats = {"requests": 0, "png_hits": 0, "renders": 0,
                      "render_seconds": 0.0}

    # -- sample ------------------------------------------------------------
    def _close(self):
        if self._h5 is not None:
            try:
                self._h5.close()
            except Exception:
                pass
        self._h5 = None
        self._sample_id = None
        self._pngs.clear()

    def close(self):
        with self._lock:
            self._close()

    @property
    def open_sample(self):
        return self._sample_id

    def handle(self, sample_id: str) -> AuroraPredictions:
        with self._lock:
            if self._sample_id != sample_id:
                if self.store is not None:
                    sample = self.store.get(sample_id)   # KeyError if unknown
                    h5 = self.store.open(sample_id)      # merged multi-file view
                else:
                    h5 = sample_id                       # allow a direct path
                self._close()
                self._h5 = AuroraPredictions(h5)
                self._sample_id = sample_id
            return self._h5

    # -- metadata ----------------------------------------------------------
    def meta(self, sample_id: str) -> dict:
        with self._lock:
            hit = self._meta.get(sample_id)
            if hit is not None:
                self._meta.move_to_end(sample_id)
                return hit
        h = self.handle(sample_id)
        with self._lock:
            info = h.describe()
            info["sample_id"] = sample_id
            self._meta[sample_id] = info
            while len(self._meta) > self._max_meta:
                self._meta.popitem(last=False)
            return info

    # -- renders -----------------------------------------------------------
    def render(self, sample_id: str, cell_px: int, kind: str, name: str,
               target_w: int = TARGET_WIDTH):
        """(png, info) for a gene or cell type, cached."""
        h = self.handle(sample_id)
        try:
            layer = h.layer(cell_px)
        except KeyError:
            # a bad resolution is a client error, not a missing sample
            raise ValueError(
                f"no {cell_px} px layer for {sample_id}; available: "
                f"{sorted(h.layers)}"
            )
        if kind == "gene":
            idx = h.gene_index()
        elif kind == "celltype":
            idx = h.cell_type_index()
        else:
            raise ValueError(f"unknown kind {kind!r}")
        if name not in idx:
            raise KeyError(name)
        index = idx[name]

        key = (sample_id, cell_px, kind, index, target_w)
        with self._lock:
            self.stats["requests"] += 1
            hit = self._pngs.get(key)
            if hit is not None:
                self._pngs.move_to_end(key)
                self.stats["png_hits"] += 1
                return hit[0], dict(hit[1], cached=True)

        png, info = render_view(layer, kind, index, target_w)
        info["name"] = name
        with self._lock:
            self.stats["renders"] += 1
            self.stats["render_seconds"] += info["raster_seconds"] + info["png_seconds"]
            self._pngs[key] = (png, info)
            while len(self._pngs) > self._max_pngs:
                self._pngs.popitem(last=False)
        return png, info


# --------------------------------------------------------------------------
# H&E images
# --------------------------------------------------------------------------
def tissue_mask_from_image(lum: np.ndarray, white: int = 220) -> np.ndarray:
    """
    Rough tissue mask from an H&E luminance raster: tissue is anything that is
    not near-white.  Used to align an H&E image with the prediction extent.
    """
    return lum < white


def align_report(img: Image.Image, layer, white: int = 220) -> dict:
    """
    Compare the prediction footprint with the tissue in an H&E image.

    Both are evaluated on the layer's **image-pixel extent** (width x height in
    rescaled-image pixels), not on the compressed patch raster: the H&E is
    resized to that extent and the footprint is rasterised onto the same grid.
    Doing it on the patch raster distorts the tissue mask whenever the lattice
    has gaps, which understates the agreement.

    A correct pairing scores a high IoU; a mismatch does not.

    The image is downscaled to the model grid (a separate pre-generated pyramid is
    used for display), then **cropped** to the prediction extent.  Cropping
    matters: the prediction covers only part of the slide, and stretching the
    whole image onto the prediction extent would compress the tissue and
    understate the agreement.
    """
    he_w, he_h = img.size
    _flat, W, H, xs, ys = layer.grid()
    max_dim = 1400                      # downsample for speed; ratios are what matter
    scale = min(1.0, max_dim / max(W, H))
    gw, gh = max(1, int(round(W * scale))), max(1, int(round(H * scale)))

    lum = np.asarray(img.convert("L").resize((gw, gh), Image.BILINEAR),
                     dtype=np.float32)
    tissue = tissue_mask_from_image(lum, white=white)

    # prediction extent in downscaled image pixels, clipped to the image
    x0, y0 = layer.origin()
    cell = layer.cell_px
    px0 = max(0, int(round(x0 * gw / he_w)))
    py0 = max(0, int(round(y0 * gh / he_h)))
    px1 = min(gw, int(round((x0 + W * cell) * gw / he_w)))
    py1 = min(gh, int(round((y0 + H * cell) * gh / he_h)))
    if px1 <= px0 or py1 <= py0:
        return {"iou": 0.0, "dice": 0.0, "he_tissue_fraction": 0.0,
                "pred_tissue_fraction": 0.0, "grid": [0, 0],
                "error": "prediction extent does not overlap the image"}

    cw, ch = px1 - px0, py1 - py0
    tissue_c = tissue[py0:py1, px0:px1]

    # Rasterise the footprint directly in image coordinates.  NOTE: `grid()`'s
    # xs/ys hold the real image coordinates of the compressed raster's columns and
    # rows - they are NOT 0..W, so they must be mapped through the image scale.
    # One downscaled pixel per cell is enough for an IoU; no loop needed.
    fx = np.rint(xs * gw / he_w - px0).astype(np.int64)
    fy = np.rint(ys * gh / he_h - py0).astype(np.int64)
    inx = (fx >= 0) & (fx < cw)
    iny = (fy >= 0) & (fy < ch)
    # the presence mask is indexed [row, col] = [y, x]
    present = layer.mask()[np.ix_(iny, inx)]
    foot = np.zeros((ch, cw), dtype=bool)
    foot[np.ix_(fy[iny], fx[inx])] = present

    inter = int((tissue_c & foot).sum())
    union = int((tissue_c | foot).sum())
    return {
        "iou": round(inter / union, 4) if union else 0.0,
        "dice": round(2 * inter / (tissue_c.sum() + foot.sum()), 4)
        if (tissue_c.sum() + foot.sum()) else 0.0,
        "he_tissue_fraction": round(float(tissue_c.mean()), 4),
        "pred_tissue_fraction": round(float(foot.mean()), 4),
        "grid": [cw, ch],
    }






