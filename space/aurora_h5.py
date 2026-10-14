"""
AURORA prediction readers.

Two on-disk schemas exist and are handled by one interface:

``v2`` (``predictions_7um_28um``, ``*-highres-pred.h5``)
    * ``exp_pred_<cell>`` : (genes, patches) uint8, zstd + ``spatial_gene_delta_v1``
    * ``prop_pred_<cell>``: (celltypes, patches) uint8, log-scale, same encoding
    * ``patch_loc_<cell>``: (N, 2) int32 **(x, y)** top-left corners
    * ``row_order_<cell>``: stored row -> original row (not needed to render)
    * cell size 14 px = 7 um, 56 px = 28 um

``v1`` (``predictions_112um``, ``*-pred.h5``)
    * ``exp_pred_224``  : (patches, genes) uint8, **no delta encoding**
    * ``prop_pred_224`` : (patches, celltypes) uint8, log-scale
    * ``patch_loc_224`` : (N, 2) int16 **(x, y)**, stored with a negative offset
    * cell size 224 px = 112 um

The datasets in v2 are compressed with a third-party zstd HDF5 filter (id 32015)
that h5py has no plugin for, so reads go through raw chunk decoding.
"""

from __future__ import annotations

import os

import numpy as np
import h5py

try:
    from imagecodecs import zstd_decode as _zstd_decode

    _BACKEND = "imagecodecs"
except Exception:  # pragma: no cover
    from numcodecs.zstd import Zstd

    _zc = Zstd(level=9)
    _zstd_decode = lambda buf: _zc.decode(np.frombuffer(buf, dtype=np.uint8))
    _BACKEND = "numcodecs"


# cell size (px, at 0.5 um/px) -> microns
CELL_MICRONS = {14: 7, 56: 28, 224: 112}

DELTA_BLOCK_DEFAULT = 16384


# --------------------------------------------------------------------------
# filter chain reversal (v2)
# --------------------------------------------------------------------------
def unshuffle(buf: bytes, elem_size: int) -> bytes:
    """Inverse of the HDF5 shuffle filter."""
    if elem_size <= 1:
        return buf
    a = np.frombuffer(buf, dtype=np.uint8)
    n = a.size // elem_size
    if n * elem_size != a.size:
        raise ValueError("buffer is not a multiple of the element size")
    return a.reshape(elem_size, n).T.copy().tobytes()


def _chain_of(dset) -> list[tuple[str, int]]:
    dcpl = dset.id.get_create_plist()
    chain = []
    for i in range(dcpl.get_nfilters()):
        fid, _flags, cd_values, name = dcpl.get_filter(i)
        nm = name if isinstance(name, str) else name.decode()
        if fid == 2:
            chain.append(("shuffle", int(cd_values[0]) if len(cd_values) else 1))
        elif fid == 32015:
            chain.append(("zstd", int(cd_values[0]) if len(cd_values) else 0))
        else:
            chain.append((f"filter{fid}", 0))
    return chain


def _decode_chunk(raw: bytes, chain) -> bytes:
    out = raw
    for name, param in reversed(chain):
        if name == "zstd":
            out = _zstd_decode(out)
            if hasattr(out, "tobytes"):
                out = out.tobytes()
        elif name == "shuffle":
            out = unshuffle(out, param)
        else:
            raise NotImplementedError(f"unsupported HDF5 filter: {name}")
    return out


CACHE_MAX_ENTRIES = 64
# The shared decode cache only holds small 1-D arrays (gene maxima etc.); the
# big per-gene chunk slices are deliberately NOT cached (see read_array) because
# a full channel block is 4 MB and a few views used to accumulate ~390 MB, which
# is fatal on a 512 MB host.  The PNG cache in viewer_core covers exact repeats.
CACHE_MAX_BYTES = int(os.environ.get("AURORA_DECODE_CACHE_MB", "32")) << 20


def _cache_bytes(cache):
    return sum(int(v.nbytes) for v in cache.values())


def _cache_put(cache, key, arr):
    cache[key] = arr
    while len(cache) > CACHE_MAX_ENTRIES or _cache_bytes(cache) > CACHE_MAX_BYTES:
        if len(cache) <= 1:
            break
        cache.pop(next(iter(cache)))


def _chunk(dset, chain, off0, off1, cache, store):
    key = (dset.name, off0, off1)
    hit = cache.get(key)
    if hit is not None:
        return hit
    shape, chunks, dtype = dset.shape, dset.chunks, dset.dtype
    off = (off0, off1) if len(shape) > 1 else (off0,)
    _, raw = dset.id.read_direct_chunk(off)
    buf = _decode_chunk(raw, chain)
    nr = min(off0 + chunks[0], shape[0]) - off0
    if len(shape) == 1:
        # flat layout: the logical values are the first `nr` bytes, padding after
        arr = np.frombuffer(buf, dtype=dtype, count=nr)
    else:
        nc = min(off1 + chunks[1], shape[1]) - off1
        total = len(buf)
        # THE CRITICAL PART: a partially written chunk still occupies its full
        # physical extent on disk, and the buffer is laid out with the PHYSICAL
        # strides.  Reshaping by the logical extent (nr x nc) therefore shifts
        # every channel after the first into the padding - each gene after the
        # first in the final block read garbage.  Reshape by the physical extent
        # and slice the logical part out.
        phys_elems = chunks[0] * chunks[1]
        phys_bytes = phys_elems * dtype.itemsize
        if total == phys_bytes:
            full = np.frombuffer(buf, dtype=dtype, count=phys_elems)
            arr = np.ascontiguousarray(full.reshape(chunks[0], chunks[1])[:nr, :nc])
        elif (total // dtype.itemsize) % nr == 0:
            # the buffer holds only logical rows; second axis is not padded
            rows = (total // dtype.itemsize) // nr
            arr = np.frombuffer(buf, dtype=dtype, count=rows * nr).reshape(nr, rows)
            arr = np.ascontiguousarray(arr[:, :nc])
        else:
            raise ValueError(
                f"unexpected chunk buffer {total} bytes for chunk {off} "
                f"({chunks=}, dtype={dtype}, logical {nr}x{nc})")
    if store:
        _cache_put(cache, key, arr)
    return arr


def _as_indexer(sel, n: int) -> np.ndarray:
    if isinstance(sel, slice):
        return np.arange(*sel.indices(n), dtype=np.int64)
    a = np.asarray(sel)
    if a.ndim == 0:
        a = a.reshape(1)
    a = a.astype(np.int64, copy=False)
    if a.size and (a.min() < 0 or a.max() >= n):
        raise IndexError(f"index out of range for {n}")
    return a


def read_array(dset, rows=None, cols=None, cache=None) -> np.ndarray:
    """Read ``dset[rows][:, cols]``, decoding only intersecting chunks."""
    shape, dtype = dset.shape, dset.dtype
    view = cache if cache is not None else {}
    store = cache is not None

    if dset.chunks is None:
        r = slice(None) if rows is None else np.asarray(rows)
        if len(shape) == 1:
            return dset[r]
        c = slice(None) if cols is None else np.asarray(cols)
        return dset[r, :][:, c]

    chain = _chain_of(dset)

    if len(shape) == 1:
        idx = (np.arange(shape[0], dtype=np.int64) if rows is None
               else np.atleast_1d(np.asarray(rows, dtype=np.int64)))
        out = np.empty(idx.size, dtype=dtype)
        if idx.size == 0:
            return out
        c0 = dset.chunks[0]
        order = np.argsort(idx, kind="stable")
        s = idx[order]
        blocks = s // c0
        for grp in np.split(np.arange(s.size), np.flatnonzero(np.diff(blocks)) + 1):
            if grp.size == 0:
                continue
            rb = int(blocks[grp[0]])
            arr = _chunk(dset, chain, rb * c0, None, view, store)
            out[grp] = arr[s[grp] - rb * c0]
        res = np.empty_like(out)
        res[order] = out
        return res

    r_idx = (np.arange(shape[0], dtype=np.int64) if rows is None
             else np.atleast_1d(np.asarray(rows, dtype=np.int64)))
    c_idx = (np.arange(shape[1], dtype=np.int64) if cols is None
             else np.atleast_1d(np.asarray(cols, dtype=np.int64)))
    out = np.empty((r_idx.size, c_idx.size), dtype=dtype)
    if r_idx.size == 0 or c_idx.size == 0:
        return out

    c0, c1 = dset.chunks[0], dset.chunks[1]
    a0, a1 = c0 < shape[0], c1 < shape[1]
    row_blocks = np.unique(r_idx // c0) if a0 else np.array([0])
    col_blocks = np.unique(c_idx // c1) if a1 else np.array([0])
    # 2-D reads are never written to the shared chunk cache: a single-gene view
    # still decodes whole channel blocks (256 x 16384 = 4 MB each), and caching
    # those pushed memory past 390 MB on a 512 MB host for almost no reuse (the
    # request-specific slice is what's needed, and exact repeats are already
    # covered by the PNG cache).  Higher-level caches (patch_loc, grid) memoise
    # the expensive full reads themselves.
    for rb in row_blocks:
        rmask = (r_idx // c0) == rb if a0 else np.ones(r_idx.size, bool)
        r_local = r_idx[rmask] - rb * c0 if a0 else r_idx[rmask]
        for cb in col_blocks:
            cmask = (c_idx // c1) == cb if a1 else np.ones(c_idx.size, bool)
            c_local = c_idx[cmask] - cb * c1 if a1 else c_idx[cmask]
            arr = _chunk(dset, chain,
                         int(rb * c0) if a0 else 0,
                         int(cb * c1) if a1 else 0,
                         {}, False)
            out[np.ix_(rmask, cmask)] = arr[np.ix_(r_local, c_local)]
    return out


def delta_decode_u8(enc: np.ndarray, block_cols: int) -> np.ndarray:
    """Inverse of ``spatial_gene_delta_v1`` (uint8 modulo-256 block cumsum)."""
    if enc.dtype != np.uint8:
        enc = enc.astype(np.uint8)
    out = np.empty_like(enc)
    n = enc.shape[1]
    for start in range(0, n, block_cols):
        stop = min(start + block_cols, n)
        out[:, start:stop] = np.cumsum(enc[:, start:stop], axis=1, dtype=np.uint8)
    return out


def _s(x):
    return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)


# --------------------------------------------------------------------------
# per-resolution view
# --------------------------------------------------------------------------
class Layer:
    """
    One resolution level of a sample.

    Exposes a uniform interface over the two on-disk schemas: gene and cell-type
    codes as ``(patches, channels)`` uint8, plus the patch grid.
    """

    def __init__(self, owner: "AuroraPredictions", cell_px: int):
        self.owner = owner
        self.f = owner.f
        self.cell_px = cell_px
        self.microns = CELL_MICRONS.get(cell_px, cell_px // 2)
        self.suffix = str(cell_px)
        self.is_v1 = owner.schema == 1

        self._exp = f"exp_pred_{self.suffix}"
        self._prop = f"prop_pred_{self.suffix}"
        self._loc = f"patch_loc_{self.suffix}"
        if self._exp not in self.f:
            raise KeyError(f"no dataset {self._exp}")

        d = self.f[self._exp]
        if self.is_v1:
            # (patches, genes), absolute codes
            self.n_patches, self.n_genes = int(d.shape[0]), int(d.shape[1])
            self._block = None
        else:
            # (genes, patches)
            self.n_genes, self.n_patches = int(d.shape[0]), int(d.shape[1])
            self._block = int(d.attrs.get("delta_block_rows", DELTA_BLOCK_DEFAULT))

        self._loc_cache = None
        self._grid_cache = None

    # -- geometry ---------------------------------------------------------
    @property
    def name(self) -> str:
        return f"{self.microns}um"

    @property
    def cell_um(self) -> float:
        return self.cell_px / 2.0        # 0.5 um per pixel

    def patch_loc(self) -> np.ndarray:
        """
        (N, 2) int32 patch top-left corners in the rescaled-image pixel space.

        Stored as **(y, x)**: the raw column order is vertical-then-horizontal,
        which was settled empirically by matching the prediction footprint
        against the tissue in a real rescaled H&E (IoU 0.83 for (y, x) versus
        0.19 for (x, y) - see DATA_REFERENCE.md).

        Coordinates are used **as stored**; no shift is applied.  A patch may
        start at a negative coordinate, in which case only the part inside the
        image is visible, and :meth:`grid` clips it accordingly.
        """
        if self._loc_cache is None:
            self._loc_cache = read_array(self.f[self._loc],
                                         cache=self.owner._cache).astype(np.int64)
        return self._loc_cache

    def grid(self):
        """
        (flat_index, W, H, x_axis, y_axis) placing every patch on a raster.

        Raster cell (row, col) = (y, x) with x = column 1 and y = column 0 of the
        stored coordinates.  Raster origins are the **minimum** stored
        coordinates, so a patch starting at a negative coordinate simply lands at
        the raster edge and only its visible part is represented - no shifting,
        which would misalign the prediction with the H&E.

        **Grid cells are spaced by the lattice pitch, not by one another.**  A
        patch lattice can have a hole in it: this sample has a ~14,000 px stretch
        of slide in the middle with no patches at all, and at 28 um the x axis is
        also 2 cells short at the left.  Laying the patches out contiguously would
        close that hole, slide everything to its right across by the width of the
        gap, and butt two separate pieces of tissue together - which reads as a
        seam with a different colour range either side.  Placing each patch at
        ``round((x - x0) / pitch)`` keeps the hole a hole, so the raster matches
        the image.  Reference behaviour: the plotting code draws every patch at
        its own stored coordinate and never builds a raster at all.

        ``x_axis`` / ``y_axis`` therefore hold the true image coordinate of each
        raster column and row.

        Memoised: it reads the full ``patch_loc`` array and sorts it, which costs
        far more than the per-view work that follows.
        """
        if self._grid_cache is not None:
            return self._grid_cache
        loc = self.patch_loc()
        y = loc[:, 0]
        x = loc[:, 1]

        def axis(values):
            """(coordinate per grid index, index per value) on the lattice pitch."""
            u = np.unique(values)
            if u.size == 1:
                return u, np.zeros(values.shape, dtype=np.int64)
            steps = np.diff(u)
            pitch = int(steps[steps > 0].min()) if (steps > 0).any() else 1
            idx = np.rint((values - u.min()) / pitch).astype(np.int64)
            n = int(idx.max()) + 1
            # the coordinate that each grid index stands for
            coord = u.min() + np.arange(n, dtype=np.int64) * pitch
            return coord, idx

        xs, xpos = axis(x)
        ys, ypos = axis(y)
        W = xs.size
        flat = ypos * W + xpos
        self._grid_cache = (flat, W, ys.size, xs, ys)
        return self._grid_cache

    def origin(self) -> tuple[int, int]:
        """
        (x0, y0) of the top-left corner of raster cell (0, 0), in image pixels.

        Normally (0, 0), but a grid whose smallest stored coordinate is negative
        starts left of / above the image origin, so the rendered raster has to be
        placed at that offset to line up with the H&E.
        """
        loc = self.patch_loc()
        return int(loc[:, 1].min()), int(loc[:, 0].min())

    def extent(self) -> tuple[int, int]:
        """
        (width, height) of the predicted area in rescaled-image pixels.

        Measured on the raster's own coordinate axis, so a grid starting at a
        negative coordinate is reported at its full size and positioned with
        :meth:`origin`, and any hole in the patch lattice is included in the span
        rather than being closed up.
        """
        _flat, W, H, xs, ys = self.grid()
        w = int(xs[-1] - xs[0]) + self.cell_px
        h = int(ys[-1] - ys[0]) + self.cell_px
        return w, h

    def bounds(self) -> tuple[int, int, int, int]:
        """
        (x0, y0, x1, y1) of the prediction in image pixels, clipped to the part
        that lies inside the image.  Patches starting at a negative coordinate
        contribute only their visible portion.
        """
        x0, y0 = self.origin()
        w, h = self.extent()
        return x0, y0, x0 + w, y0 + h

    # -- values -----------------------------------------------------------
    def _idx(self, sel, n: int, names: list[str] | None = None) -> np.ndarray:
        """
        Normalise a channel selector to indices.

        Accepts None (every channel), strings (matched against `names`), or
        anything array-like of integers.
        """
        if sel is None:
            return np.arange(n, dtype=np.int64)
        if isinstance(sel, str):
            sel = [sel]
        seq = list(sel) if not isinstance(sel, (int, np.integer)) else [sel]
        if not seq:
            return np.arange(n, dtype=np.int64)
        if names is not None and any(isinstance(x, str) for x in seq):
            lookup = {nm: i for i, nm in enumerate(names)}
            out = []
            for x in seq:
                if isinstance(x, str):
                    if x not in lookup:
                        raise KeyError(x)
                    out.append(lookup[x])
                else:
                    out.append(int(x))
            return np.asarray(out, dtype=np.int64)
        return np.atleast_1d(np.asarray(seq, dtype=np.int64))

    def exp_codes(self, genes=None, rows=None) -> np.ndarray:
        """(len(rows), len(genes)) absolute uint8 codes.  Genes may be names."""
        d = self.f[self._exp]
        g = self._idx(genes, self.n_genes, self.owner.genes)
        if self.is_v1:
            if rows is None:
                return np.ascontiguousarray(d[:, g])
            return np.ascontiguousarray(d[np.atleast_1d(np.asarray(rows, np.int64)), :][:, g])
        block = self._block
        if rows is None:
            enc = read_array(d, g, None, cache=self.owner._cache)
            return delta_decode_u8(enc, block).T
        want = np.atleast_1d(np.asarray(rows, dtype=np.int64))
        if want.size == 0:
            return np.empty((0, g.size), dtype=np.uint8)
        lo = int(want.min()) // block * block
        hi = min(int(want.max()) + 1, self.n_patches)
        enc = read_array(d, g, np.arange(lo, hi), cache=self.owner._cache)
        dec = delta_decode_u8(enc, block)
        return dec[:, want - lo].T

    def prop_codes(self, celltypes=None, rows=None) -> np.ndarray:
        """(len(rows), len(celltypes)) absolute uint8 codes (log-scale values)."""
        if self._prop not in self.f:
            raise KeyError(f"no dataset {self._prop}")
        d = self.f[self._prop]
        c = self._idx(celltypes, len(self.owner.cell_types), self.owner.cell_types)
        if self.is_v1:
            if rows is None:
                return np.ascontiguousarray(d[:, c])
            return np.ascontiguousarray(d[np.atleast_1d(np.asarray(rows, np.int64)), :][:, c])
        block = int(d.attrs.get("delta_block_rows", DELTA_BLOCK_DEFAULT))
        if rows is None:
            enc = read_array(d, c, None, cache=self.owner._cache)
            return delta_decode_u8(enc, block).T
        want = np.atleast_1d(np.asarray(rows, dtype=np.int64))
        if want.size == 0:
            return np.empty((0, c.size), dtype=np.uint8)
        lo = int(want.min()) // block * block
        hi = min(int(want.max()) + 1, self.n_patches)
        enc = read_array(d, c, np.arange(lo, hi), cache=self.owner._cache)
        dec = delta_decode_u8(enc, block)
        return dec[:, want - lo].T

    def exp_values(self, genes=None, rows=None) -> np.ndarray:
        """Expression in the model's own units, float32."""
        codes = self.exp_codes(genes, rows)
        mx = self.exp_max(self._idx(genes, self.n_genes, self.owner.genes))
        return codes.astype(np.float32) / 255.0 * mx[None, :]

    def proportions(self, celltypes=None, rows=None) -> np.ndarray:
        """
        Cell-type proportions as probabilities that sum to 1 per patch.

        ``prop_pred_*`` stores **log** values, so we exponentiate and renormalise
        (see DATA_REFERENCE.md).

        The softmax is always taken over **every** cell type and the result is
        subset afterwards.  Normalising over a subset would divide a single
        channel by itself, turning every value into 1.0 鈥?which silently produced
        blank cell-type maps.
        """
        # all channels, normalised together
        all_codes = self.prop_codes(None, rows)
        mx = self.prop_max(None)
        logits = all_codes.astype(np.float32) / 255.0 * mx[None, :]
        logits -= logits.max(axis=1, keepdims=True)
        p = np.exp(logits)
        p /= p.sum(axis=1, keepdims=True)
        sel = self._idx(celltypes, len(self.owner.cell_types), self.owner.cell_types)
        return np.ascontiguousarray(p[:, sel])

    def exp_max(self, genes=None) -> np.ndarray:
        """Per-gene quantisation maxima; accepts indices or gene names."""
        mx = self.f[f"exp_max_{self.suffix}"][()]
        if genes is None:
            return mx
        return mx[self._idx(genes, self.n_genes, self.owner.genes)]

    def prop_max(self, celltypes=None) -> np.ndarray:
        """Per-cell-type maxima; accepts indices or cell-type names."""
        mx = self.f[f"prop_max_{self.suffix}"][()]
        if celltypes is None:
            return mx
        return mx[self._idx(celltypes, len(self.owner.cell_types),
                            self.owner.cell_types)]

    # -- rasters ----------------------------------------------------------
    def mask(self) -> np.ndarray:
        """(H, W) bool: True where a patch exists, on the dense raster."""
        flat, W, H, _xs, _ys = self.grid()
        out = np.zeros(W * H, dtype=bool)
        out[flat] = True
        return out.reshape(H, W)

    def lattice_is_uniform(self) -> bool:
        """
        True when every occupied grid column and row holds patches.

        The grid itself always steps by the lattice pitch, so its axis is uniform
        by construction; what can vary is occupancy.  This sample has a ~14 kpx
        stretch of slide in the middle with no patches, which shows up as a run of
        blank columns rather than as a larger axis step.
        """
        return not self.lattice_gaps()["blank_columns"] and \
            not self.lattice_gaps()["blank_rows"]

    def lattice_gaps(self) -> dict:
        """
        Holes in the patch lattice, in image pixels.

        ``x_max_step`` / ``y_max_step`` are the largest jump between consecutive
        **occupied** coordinates, so a hole shows up here rather than as a larger
        grid step.  ``blank_columns`` / ``blank_rows`` count grid lines with no
        patches at all.
        """
        loc = self.patch_loc()
        _flat, W, H, _xs, _ys = self.grid()
        ux = np.unique(loc[:, 1])
        uy = np.unique(loc[:, 0])
        return {
            "x_max_step": int(np.diff(ux).max()) if ux.size > 1 else 0,
            "y_max_step": int(np.diff(uy).max()) if uy.size > 1 else 0,
            "blank_columns": int(W - ux.size),
            "blank_rows": int(H - uy.size),
            "cell_px": self.cell_px,
        }


# --------------------------------------------------------------------------
# file wrapper
# --------------------------------------------------------------------------
class AuroraPredictions:
    """
    Reader for both AURORA prediction schemas.

    ``schema`` is 2 for ``*-highres-pred.h5`` (7/28 um delta files) and 1 for
    ``*-pred.h5`` (112 um plain files); it is detected from the datasets present.
    """

    def __init__(self, path):
        if isinstance(path, h5py.File):
            sources = [path]
        elif isinstance(path, (list, tuple)):
            sources = list(path)
        else:
            sources = [path]

        self.readers: list[AuroraPredictions] = []
        for src in sources:
            if isinstance(src, AuroraPredictions):
                self.readers.append(src)
            else:
                self.readers.append(_SingleFile(src))

        if not self.readers:
            raise ValueError("AuroraPredictions needs at least one source")

        # merged view: every resolution level across every file
        merged: dict[int, Layer] = {}
        for r in self.readers:
            for cell_px, layer in r.layers.items():
                merged.setdefault(cell_px, layer)
        self.layers = merged
        if not self.layers:
            raise ValueError(f"{path}: no prediction datasets found")

        self.schema = self.readers[0].schema
        self.path = " + ".join(str(r.path) for r in self.readers)
        self._cache: dict = {}
        # a shared chunk cache so a merged sample does not double-cache
        for r in self.readers:
            r._cache = self._cache
        self._genes = None
        self._cell_types = None

    # -- lifecycle --------------------------------------------------------
    def close(self):
        for r in self.readers:
            r.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __repr__(self):
        return (f"AuroraPredictions({self.path!r}, "
                f"{len(self.genes)} genes, levels={sorted(self.layers)})")

    # -- metadata ---------------------------------------------------------
    @property
    def genes(self) -> list[str]:
        if self._genes is None:
            self._genes = self.readers[0].genes
        return self._genes

    @property
    def cell_types(self) -> list[str]:
        if self._cell_types is None:
            self._cell_types = self.readers[0].cell_types
        return self._cell_types

    def gene_index(self) -> dict[str, int]:
        return {g: i for i, g in enumerate(self.genes)}

    def cell_type_index(self) -> dict[str, int]:
        return {c: i for i, c in enumerate(self.cell_types)}

    def layer(self, cell_px: int) -> Layer:
        if cell_px not in self.layers:
            raise KeyError(f"no layer for {cell_px} px; have {sorted(self.layers)}")
        return self.layers[cell_px]

    def default_layer(self) -> Layer:
        return self.layers[min(self.layers)]

    def describe(self) -> dict:
        out = {
            "schema": f"v{self.schema}",
            "sources": [os.path.basename(str(r.path)) for r in self.readers],
            "genes": self.genes,
            "cell_types": self.cell_types,
            "layers": [],
        }
        for cell_px, layer in sorted(self.layers.items()):
            extent = list(layer.extent())
            x0, y0 = layer.origin()
            _flat, gw, gh, _xs, _ys = layer.grid()
            out["layers"].append({
                "cell_px": cell_px,
                "microns": layer.microns,
                "name": layer.name,
                "patches": layer.n_patches,
                # the raster's own size, which is larger than patches/cell when the
                # lattice has a hole in it
                "grid": [gw, gh],
                "extent": extent,
                # where the grid's top-left corner sits in image pixels; non-zero
                # when patches start at negative coordinates
                "origin": [x0, y0],
                "uniform_lattice": layer.lattice_is_uniform(),
                "source": os.path.basename(str(layer.owner.path)),
            })
        return out


class _SingleFile:
    """One prediction file (either schema)."""

    def __init__(self, path):
        if isinstance(path, h5py.File):
            self.f = path
            self.path = getattr(path, "filename", "<open h5py.File>")
            self._owns = False
        else:
            self.path = str(path)
            self.f = h5py.File(self.path, "r")
            self._owns = True
        self._cache: dict = {}
        self._genes = None
        self._cell_types = None
        self.schema = self._detect_schema()
        self.layers = self._discover_layers()
        if not self.layers:
            raise ValueError(f"{self.path}: no prediction datasets found")

    def _detect_schema(self) -> int:
        if "exp_pred_224" in self.f and "exp_pred_14" not in self.f:
            return 1
        return 2

    def _discover_layers(self) -> dict[int, Layer]:
        out = {}
        for cell_px in sorted(CELL_MICRONS):
            if f"exp_pred_{cell_px}" in self.f:
                try:
                    out[cell_px] = Layer(self, cell_px)
                except Exception:
                    continue
        return out

    def close(self):
        if self._owns:
            self.f.close()
        else:
            try:
                self.f.id._close()
            except Exception:
                pass

    @property
    def genes(self) -> list[str]:
        if self._genes is None:
            self._genes = [_s(x) for x in self.f["genes"][()]]
        return self._genes

    @property
    def cell_types(self) -> list[str]:
        if self._cell_types is None:
            self._cell_types = [_s(x) for x in self.f["cell_types"][()]]
        return self._cell_types


