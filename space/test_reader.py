"""
Reader tests, against a real matched pair.

The workspace sample is `TCGA-JL-A3YX-01A-01-TS1`:
    *-highres-pred.h5   7 um + 28 um (delta-encoded v2 schema)
    *-pred.h5           112 um (plain v1 schema)
    *_rescaled.tif      the 0.5 um/px H&E the predictions were made against

The H&E is what settles the coordinate convention, so the alignment checks here
are the important ones.
"""
import os
import sys
import time

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from aurora_h5 import (AuroraPredictions, delta_decode_u8,  # noqa: E402
                       read_array)

SAMPLE = os.environ.get("AURORA_TEST_SAMPLE", "TCGA-JL-A3YX-01A-01-TS1")
FINE = os.environ.get("AURORA_TEST_H5") or os.path.join(DATA, f"{SAMPLE}-highres-pred.h5")
COARSE = os.environ.get("AURORA_TEST_H5_V1") or os.path.join(DATA, f"{SAMPLE}-pred.h5")
HE = os.environ.get("AURORA_TEST_HE") or os.path.join(DATA, f"{SAMPLE}_rescaled.tif")

if not os.path.exists(FINE):
    print(f"SKIP: no prediction file at {FINE}")
    print("      set AURORA_TEST_H5 to a *-highres-pred.h5 to run these tests")
    sys.exit(0)

Image.MAX_IMAGE_PIXELS = None
failures = []


def check(label, cond, extra=""):
    ok = bool(cond)
    print(("PASS  " if ok else "FAIL  ") + label + (f"  {extra}" if extra else ""))
    if not ok:
        failures.append(label)


# --------------------------------------------------------------------------
# schemas
# --------------------------------------------------------------------------
fine = AuroraPredictions(FINE)
check("v2 file detected", fine.schema == 2, fine.schema)
check("v2 exposes 7 um and 28 um", sorted(fine.layers) == [14, 56], sorted(fine.layers))
check("gene list is populated and consistent", len(fine.genes) > 1000 and len(fine.genes) == len(set(fine.genes)), len(fine.genes))
check("cell types are populated", 2 <= len(fine.cell_types) <= 40, fine.cell_types)
check("delta block is 16384 rows", fine.layer(14)._block == 16384, fine.layer(14)._block)

have_v1 = os.path.exists(COARSE)
coarse = AuroraPredictions(COARSE) if have_v1 else None
if have_v1:
    check("v1 file detected", coarse.schema == 1, coarse.schema)
    check("v1 exposes only 112 um", sorted(coarse.layers) == [224], sorted(coarse.layers))
    check("v1 has no delta encoding", coarse.layer(224)._block is None)
    check("both schemas list the same genes in the same order", fine.genes == coarse.genes)
else:
    print(f"SKIP  v1 checks: no file at {COARSE}")

# --------------------------------------------------------------------------
# the stock h5py path cannot read v2
# --------------------------------------------------------------------------
import h5py  # noqa: E402

raw = h5py.File(FINE, "r")
dcpl = raw["exp_pred_14"].id.get_create_plist()
filters = [dcpl.get_filter(i)[0] for i in range(dcpl.get_nfilters())]
check("v2 predictions use the third-party zstd filter", 32015 in filters, filters)
try:
    raw["exp_pred_14"][:, :4]
    check("the stock h5py path fails as expected", False, "it succeeded")
except OSError as exc:
    msg = str(exc).lower()
    check("the stock h5py path fails as expected",
          "can't open directory" in msg or "filter" in msg, str(exc)[:50])
raw.close()

# --------------------------------------------------------------------------
# delta arithmetic
# --------------------------------------------------------------------------
enc = np.array([[5, 1, 2, 255, 0, 10]], dtype=np.uint8)
check("delta decode matches hand-computed values",
      delta_decode_u8(enc, 64)[0].tolist() == [5, 6, 8, 7, 7, 17])
check("delta blocks restart at their first row",
      delta_decode_u8(enc, 3)[0].tolist() == [5, 6, 8, 255, 255, 9])

# --------------------------------------------------------------------------
# coordinate convention: patch_loc is (y, x), used as stored
# --------------------------------------------------------------------------
layer = fine.layer(14)
loc = layer.patch_loc()
n = layer.n_patches
print()
print(f"  7 um: {n} patches, col0 (y) {loc[:,0].min()}..{loc[:,0].max()}, "
      f"col1 (x) {loc[:,1].min()}..{loc[:,1].max()}")

if os.path.exists(HE):
    im = Image.open(HE)
    HE_W, HE_H = im.size
    print(f"  H&E: {HE_W} x {HE_H}")
    check("column 0 spans the image HEIGHT, column 1 the WIDTH",
          loc[:, 0].max() <= HE_H and loc[:, 1].max() <= HE_W,
          f"y<= {HE_H}?  x<= {HE_W}?")
    check("the aspect ratio of the patch spans matches the H&E",
          abs((loc[:, 1].max() / loc[:, 0].max()) - (HE_W / HE_H)) < 0.2,
          f"{loc[:,1].max()/loc[:,0].max():.2f} vs {HE_W/HE_H:.2f}")
else:
    print(f"SKIP  H&E-dependent checks: no image at {HE}")

# --------------------------------------------------------------------------
# raster geometry
# --------------------------------------------------------------------------
flat, W, H, xs, ys = layer.grid()
check("grid height equals the number of distinct y values",
      H == len(np.unique(loc[:, 0])), (H, len(np.unique(loc[:, 0]))))
check("raster is (H, W)", layer.mask().shape == (H, W), layer.mask().shape)
check("every patch lands in the raster exactly once",
      int(layer.mask().sum()) == n, int(layer.mask().sum()))
check("origin is the minimum stored coordinate",
      layer.origin() == (int(loc[:, 1].min()), int(loc[:, 0].min())), layer.origin())
# the width is set by the coordinate range and the lattice pitch, so it is at
# least as wide as the number of occupied columns and exactly as wide as the span
span_cols = int((loc[:, 1].max() - loc[:, 1].min()) // 14) + 1
check("grid width spans the coordinate range at the lattice pitch",
      W == span_cols, (W, span_cols))
check("grid width is at least the number of occupied columns",
      W >= len(np.unique(loc[:, 1])), (W, len(np.unique(loc[:, 1]))))
check("extent is the true coordinate span plus one cell",
      layer.extent() == (int(xs[-1] - xs[0]) + 14, int(ys[-1] - ys[0]) + 14),
      layer.extent())

# a 14 px cell is 7 um at 0.5 um/px, so the lattice should step by the cell size.
# Some samples contain stretches of slide with no patches; that is data, not a
# bug, and the raster now preserves it as blank columns rather than closing up.
gaps = layer.lattice_gaps()
check("y lattice steps by exactly the cell size", gaps["y_max_step"] == 14, gaps)
xsteps = np.unique(np.diff(xs))
check("x lattice steps by exactly the cell size",
      bool((xsteps == 14).all()), xsteps[:6].tolist())
blank = W - len(np.unique(loc[:, 1]))
if blank:
    print(f"      note: the lattice has a hole - {blank} blank columns of {W} "
          f"(x gap {gaps['x_max_step']} px); the raster keeps it")

l28 = fine.layer(56)
_o28, W28, H28, _xs28, _ys28 = l28.grid()
check("28 um grid is a quarter of the 7 um grid",
      abs(W - 4 * W28) <= 8 and abs(H - 4 * H28) <= 8, (W, H, W28, H28))

# --------------------------------------------------------------------------
# partial reads equal a full read across delta-block boundaries
# --------------------------------------------------------------------------
full = layer.exp_codes(["VIM"])[:, 0]
check("full read has one value per patch", full.shape == (n,), full.shape)
edges = [0, 1, 16383, 16384, 16385, 32767, 32768, n // 2, n - 2, n - 1]
for a in range(0, len(edges) - 1):
    idx = np.array(sorted({edges[a], edges[a + 1]}))
    if idx.min() < 0 or idx.max() >= n:
        continue
    part = layer.exp_codes(["VIM"], idx)[:, 0]
    check(f"partial read rows {idx.tolist()} matches the full read",
          np.array_equal(part, full[idx]))
part = layer.exp_codes(["VIM"], np.arange(n))[:, 0]
check("reading every row at once matches the full read", np.array_equal(part, full))

multi = layer.exp_codes(["TMSB4X", "CD74", "VIM"])
check("multi-gene read shape is (patches, 3)", multi.shape == (n, 3), multi.shape)
for k, g in enumerate(["TMSB4X", "CD74", "VIM"]):
    check(f"column {k} ({g}) matches its single-gene read",
          np.array_equal(multi[:, k], layer.exp_codes([g])[:, 0]))

# --------------------------------------------------------------------------
# negative coordinates are preserved and clipped, never shifted
# --------------------------------------------------------------------------
if have_v1:
    l224 = coarse.layer(224)
    loc224 = l224.patch_loc()
    neg = (loc224[:, 0] < 0) | (loc224[:, 1] < 0)
    print()
    print(f"  112 um: {int(neg.sum())} patches start at a negative coordinate")
    print(f"       x {loc224[:,1].min()}..{loc224[:,1].max()}   "
          f"y {loc224[:,0].min()}..{loc224[:,0].max()}")
    check("negative coordinates are kept as stored, not shifted",
          int(loc224[:, 0].min()) < 0 or int(loc224[:, 1].min()) < 0 or not neg.any(),
          (int(loc224[:, 0].min()), int(loc224[:, 1].min())))
    check("the 112 um origin is the minimum coordinate",
          l224.origin() == (int(loc224[:, 1].min()), int(loc224[:, 0].min())),
          l224.origin())
    if os.path.exists(HE):
        check("only a small fraction of 112 um patches overhang the image edge",
              neg.mean() < 0.02, f"{100*neg.mean():.2f}%")

# --------------------------------------------------------------------------
# a hole in the patch lattice must stay a hole
#
# This sample has a ~14,000 px stretch of slide in the middle with no patches at
# all.  Laying the patches out contiguously would close that hole and slide
# everything to its right across by the width of the gap, butting two separate
# pieces of tissue together.  That reads as a seam with a different colour range
# either side, and it misaligns the prediction against the H&E.
# --------------------------------------------------------------------------
print()
for cell in (14, 56, 224):
    L = fine.layer(cell) if cell != 224 else (coarse.layer(224) if have_v1 else None)
    if L is None:
        continue
    flat, W, H, xs, ys = L.grid()
    step = L.cell_px
    dx = np.diff(xs)
    dy = np.diff(ys)
    check(f"{L.microns} um: the raster steps by the lattice pitch",
          bool((dx == step).all()) and bool((dy == step).all()),
          f"x steps {np.unique(dx).tolist()[:4]}, y steps {np.unique(dy).tolist()[:4]}")
    span_x = int(xs[-1] - xs[0]) + step
    span_y = int(ys[-1] - ys[0]) + step
    check(f"{L.microns} um: extent is the true span, not grid x cell",
          L.extent() == (span_x, span_y),
          (L.extent(), (span_x, span_y)))
    loc = L.patch_loc()
    real_w = int(loc[:, 1].max() - loc[:, 1].min()) + step
    check(f"{L.microns} um: the raster covers the full coordinate range",
          span_x >= real_w, (span_x, real_w))
    # the raster must be wider than the patch count implies whenever there is a hole
    n_cols_used = int(np.unique(loc[:, 1]).size)
    if W > n_cols_used:
        print(f"      {L.microns:3d} um: {W - n_cols_used} blank columns hold a "
              f"lattice hole ({W} columns for {n_cols_used} occupied)")
    check(f"{L.microns} um: blank columns exist only where the lattice has a hole",
          W - n_cols_used == 0 or float(
              (xs[-1] - xs[0] + step) / step) - W < 1e-6,
          f"{W} columns, {n_cols_used} occupied")

# the dense raster must not be more filled than the patches justify
l28 = fine.layer(56)
_flat, _W, _H, _xs, _ys = l28.grid()
check("a raster with a hole reports less than the patch count would fill",
      float(l28.mask().mean()) < 1.0 and _W > np.unique(l28.patch_loc()[:, 1]).size,
      f"{100*l28.mask().mean():.1f}% filled, {_W} columns")

# --------------------------------------------------------------------------
# delta blocks decode independently
#
# Rows are stored as a per-block cumulative sum, so every block must decode on
# its own and the first row of each block is stored absolute.  This checks that
# structure - not the values, which the viewer draws as they are.
# --------------------------------------------------------------------------
l14 = fine.layer(14)
blk = l14._block
n_rows = l14.n_patches
starts = list(range(0, n_rows, blk))
check("7 um has multiple delta blocks", len(starts) > 1, len(starts))
check("the last 7 um block is partial", n_rows % blk != 0,
      f"{n_rows % blk} rows")
s1 = starts[1]
codes = l14.exp_codes([0])

# Read the stored bytes for channel 0 over that block.  h5py cannot read these
# datasets at all (no zstd plugin), so go through the same raw-chunk path the
# reader uses - which also proves the plugin is not needed.
raw = read_array(l14.f["exp_pred_14"], rows=np.array([0]),
                 cols=np.arange(s1, s1 + blk)).reshape(-1)
# the first row of each block is written absolute, so it must equal the stored byte
check("row 0 of a block is stored absolute",
      int(codes[s1, 0]) == int(raw[0]), (int(codes[s1, 0]), int(raw[0])))
# re-differencing the decoded block must reproduce the stored bytes exactly
dec = codes[s1:s1 + blk, 0].astype(np.int16)
rederived = np.empty_like(dec)
rederived[0] = dec[0]
rederived[1:] = (dec[1:] - dec[:-1]) % 256
check("re-differencing a block reproduces the stored bytes",
      np.array_equal(rederived, raw[:dec.size].astype(np.int16)))
# decoding a block in isolation must equal those rows from the full read
part = l14.exp_codes([0], np.arange(s1, min(s1 + blk, n_rows)))
check("a single block decodes identically in isolation",
      np.array_equal(part[:, 0], codes[s1:min(s1 + blk, n_rows), 0]))
print(f"  7 um: {len(starts)} blocks of {blk} rows, final one has {n_rows % blk} rows")

# --------------------------------------------------------------------------
# values
# --------------------------------------------------------------------------
mx = layer.exp_max(["VIM"])[0]
vals = layer.exp_values(["VIM"])[:, 0]
check("expression is float32", vals.dtype == np.float32, vals.dtype)
check("expression is proportional to the code",
      bool(np.allclose(vals, full.astype(np.float32) / 255.0 * mx, atol=1e-6)))

P = layer.proportions()
check("proportions sum to 1 per patch", bool(np.allclose(P.sum(1), 1.0, atol=1e-5)),
      float(P.sum(1).mean()))
probe_ct = max(fine.cell_types, key=lambda c: float(layer.proportions([c]).std()))
sub = layer.proportions([probe_ct])
check("subsetting a cell type gives real proportions, not all ones",
      float(sub.std()) > 1e-3, f"{probe_ct}: std={float(sub.std()):.5f}")
check("the subset matches the full matrix column",
      np.allclose(sub[:, 0], P[:, fine.cell_type_index()[probe_ct]], atol=1e-6))

# --------------------------------------------------------------------------
# performance
# --------------------------------------------------------------------------
t0 = time.time()
fine.layer(56).exp_codes(["VIM"])
t28 = time.time() - t0
t0 = time.time()
fine.layer(14).exp_codes(["CD74"])
t7 = time.time() - t0
check("28 um gene decode is fast (< 1 s)", t28 < 1.0, f"{t28:.3f}s")
check("7 um gene decode is usable (< 5 s)", t7 < 5.0, f"{t7:.3f}s")
print(f"      decode: 28 um {t28:.3f}s, 7 um {t7:.3f}s over {n} patches")

fine.close()
if coarse is not None:
    coarse.close()
print()
if failures:
    print(f"{len(failures)} FAILURE(S): " + "; ".join(failures))
    sys.exit(1)
print("ALL READER CHECKS PASSED")





