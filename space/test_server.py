"""
Serving-layer tests: read a real prediction file through ViewerCache and check
the metadata, the raster placement, the rendered PNGs, cell types and error
paths.  Needs a real *-highres-pred.h5; set AURORA_TEST_H5, or keep the sample
next to this folder.

Point AURORA_TEST_SAMPLE at another sample id to run against a different file.
"""
import io
import os
import sys
import time

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from samples import SampleStore           # noqa: E402
import viewer_core as C                   # noqa: E402

SAMPLE = os.environ.get("AURORA_TEST_SAMPLE", "TCGA-JL-A3YX-01A-01-TS1")
Image.MAX_IMAGE_PIXELS = None

failures = []


def check(label, cond, extra=""):
    ok = bool(cond)
    print(("PASS  " if ok else "FAIL  ") + label + (f"  {extra}" if extra else ""))
    if not ok:
        failures.append(label)


print("=" * 78)
print(f"sample {SAMPLE}")
print("=" * 78)

t0 = time.time()
store = SampleStore(backend="local", local_dir=DATA)
store.refresh()
cache = C.ViewerCache(store=store)
t_open = time.time() - t0
print(f"      store built in {t_open:.2f}s: {len(store.samples())} sample(s)")

# --------------------------------------------------------------------------
# metadata
# --------------------------------------------------------------------------
t0 = time.time()
meta = cache.meta(SAMPLE)
t_meta = time.time() - t0
genes = meta["genes"]
cell_types = meta["cell_types"]
by_cell = {L["cell_px"]: L for L in meta["layers"]}

check("metadata lists genes", len(genes) > 100, len(genes))
check("metadata lists cell types", len(cell_types) >= 2, len(cell_types))
check("metadata lists three resolutions", len(by_cell) == 3, sorted(by_cell))
check("every level reports an extent and an origin",
      all("extent" in L and "origin" in L for L in by_cell.values()))
check("the coarsest level keeps a negative origin is allowed",
      all(isinstance(v, int) for L in by_cell.values() for v in L["origin"]))
print(f"      metadata built in {t_meta:.2f}s: {len(genes)} genes, "
      f"{len(cell_types)} cell types")
for c in sorted(by_cell):
    L = by_cell[c]
    print(f"      {L['microns']:3d} um: grid {L['grid'][0]}x{L['grid'][1]} "
          f"origin {L['origin']} extent {L['extent']} patches {L['patches']}")

# --------------------------------------------------------------------------
# placement: the raster must be POSITIONED, not stretched over the image
#
# This was a real bug: the client drew the raster across the whole H&E extent.
# The prediction is narrower than the image and a layer may start at a negative
# coordinate, so the render has to report origin and extent, and decimation must
# keep the aspect ratio.
# --------------------------------------------------------------------------
for cell in (14, 56, 224):
    _png, info = cache.render(SAMPLE, cell, "gene", genes[0])
    L = by_cell[cell]
    check(f"{info['microns']} um: render reports an origin",
          "origin" in info and len(info["origin"]) == 2, info.get("origin"))
    check(f"{info['microns']} um: render origin matches the layer",
          tuple(info["origin"]) == tuple(L["origin"]),
          (info["origin"], L["origin"]))
    check(f"{info['microns']} um: extent equals grid x cell size",
          info["extent"] == [info["native_grid"][0] * info["cell_px"],
                             info["native_grid"][1] * info["cell_px"]],
          (info["extent"], info["native_grid"], info["cell_px"]))
    nw, nh = info["raster"]
    gw, gh = info["native_grid"]
    check(f"{info['microns']} um: decimation preserves the aspect ratio",
          abs((nw / nh) - (gw / gh)) < 0.02,
          f"raster {nw}x{nh} ({nw/nh:.3f}) vs grid {gw}x{gh} ({gw/gh:.3f})")

he_size = None
for cand in (f"{SAMPLE}_rescaled.tif", f"{SAMPLE}.png", f"{SAMPLE}_he.png"):
    p = os.path.join(DATA, cand)
    if os.path.exists(p):
        with Image.open(p) as im:
            he_size = [im.width, im.height]
        break

if he_size:
    print(f"      H&E native size {he_size[0]}x{he_size[1]}")
    for cell in (14, 56, 224):
        L = by_cell[cell]
        x0, y0 = L["origin"]
        w, h = L["extent"]
        slack = 2 * L["cell_px"]
        check(f"{L['microns']} um: prediction sits inside the H&E frame",
              x0 >= -slack and y0 >= -slack
              and x0 + w <= he_size[0] + slack and y0 + h <= he_size[1] + slack,
              f"x {x0}..{x0+w} of {he_size[0]}, y {y0}..{y0+h} of {he_size[1]}")
    check("the 7 um prediction is narrower than the H&E, so stretching was wrong",
          by_cell[14]["extent"][0] < he_size[0],
          f"prediction {by_cell[14]['extent'][0]} vs image {he_size[0]}")
else:
    print("      (no H&E image found; skipping frame checks)")

# --------------------------------------------------------------------------
# rendering a gene at every resolution
# --------------------------------------------------------------------------
for cell in (14, 56, 224):
    t0 = time.time()
    png, info = cache.render(SAMPLE, cell, "gene", "VIM")
    t_cold = time.time() - t0
    t0 = time.time()
    png2, info2 = cache.render(SAMPLE, cell, "gene", "VIM")
    t_warm = time.time() - t0
    img = Image.open(io.BytesIO(png))
    check(f"{info['microns']} um: PNG decodes", img.size == tuple(info["raster"]),
          (img.size, info["raster"]))
    check(f"{info['microns']} um: second render is served from cache",
          info2.get("cached") is True and png == png2)
    check(f"{info['microns']} um: payload is web-sized", info["bytes"] < 3_000_000,
          f"{info['bytes']/1024:.0f} KB")
    check(f"{info['microns']} um: mode is RGB", info["mode"] == "RGB", info["mode"])
    print(f"      {info['microns']:3d} um: cold {t_cold:.2f}s warm {t_warm:.3f}s "
          f"{info['raster'][0]}x{info['raster'][1]} {info['mode']} "
          f"{info['bytes']/1024:.1f} KB")

# --------------------------------------------------------------------------
# value / presence channels
#
# presence is what makes "no patch" transparent.  Without it the client cannot
# tell background from a legitimate expression of 0, and the whole raster is
# drawn as a rectangle of near-black over the image - a bug that made the
# prediction panel look empty.
# --------------------------------------------------------------------------
for cell in (14, 56, 224):
    png, info = cache.render(SAMPLE, cell, "gene", "VIM")
    img = Image.open(io.BytesIO(png))
    arr = np.asarray(img)
    check(f"{info['microns']} um: PNG is RGB", img.mode == "RGB", img.mode)
    values, presence = arr[..., 0], arr[..., 2]
    frac = float((presence > 0).mean())
    check(f"{info['microns']} um: presence marks part of the raster",
          0.01 < frac < 0.99, f"{100*frac:.1f}%")
    check(f"{info['microns']} um: background carries no value",
          bool((values[presence == 0] == 0).all()))
    check(f"{info['microns']} um: values are drawn not all-zero",
          int(values[presence > 0].max()) >= 200,
          int(values[presence > 0].max()))

# --------------------------------------------------------------------------
# cell types
# --------------------------------------------------------------------------
layer14 = cache.handle(SAMPLE).layer(14)
P = layer14.proportions()
check("proportions sum to 1 per patch", bool(np.allclose(P.sum(1), 1.0, atol=1e-5)),
      float(P.sum(1).mean()))
stds = [float(P[:, k].std()) for k in range(P.shape[1])]
check("several cell types vary", sum(1 for s in stds if s > 1e-4) >= 2,
      [round(s, 5) for s in stds])
probe = cell_types[int(np.argmax(stds))]
sub = layer14.proportions([probe])
check(f"a single cell-type map varies ({probe})", float(sub.std()) > 1e-3,
      f"std={float(sub.std()):.5f}")
for cell in (56, 224):
    png, info = cache.render(SAMPLE, cell, "celltype", probe)
    check(f"{info['microns']} um cell-type render", info["kind"] == "celltype"
          and info["bytes"] > 0, f"{info['bytes']} bytes")

# --------------------------------------------------------------------------
# cell-type colour spans the colormap
#
# Proportions are softmax probabilities; on this sample every cell type lies
# between about 0.095 and 0.24.  Scaling by 255 alone caps the code at ~60 and
# paints every cell type in the dark end of the colormap, so all maps look blank
# and identical.  The viewer stretches each type over the range it occupies and
# reports those bounds so the legend can state the real numbers.
# --------------------------------------------------------------------------
raw_max = float(cache.handle(SAMPLE).layer(14).proportions().max())
check("proportion codes would NOT fill the colormap unnormalised",
      raw_max * 255 < 200, f"raw max {raw_max:.4f} -> code {raw_max*255:.0f}")

for cell in (14, 56, 224):
    for name in cell_types:
        png, info = cache.render(SAMPLE, cell, "celltype", name)
        arr = np.asarray(Image.open(io.BytesIO(png)))
        vals = arr[..., 0][arr[..., 2] > 0]
        check(f"{info['microns']} um '{name}': colour fills the ramp",
              int(vals.min()) == 0 and int(vals.max()) == 255,
              f"codes {int(vals.min())}..{int(vals.max())}")
        # max-pooling preserves the maximum exactly, so the top bound must match
        col = cache.handle(SAMPLE).layer(cell).proportions([name])[:, 0]
        check(f"{info['microns']} um '{name}': raw upper bound is exact",
              abs(info["raw_hi"] - float(col.max())) < 1e-6,
              f"reported {info['raw_hi']:.6f}, actual {float(col.max()):.6f}")
        check(f"{info['microns']} um '{name}': raw bounds are a fraction",
              0.0 < info["raw_lo"] < info["raw_hi"] <= 1.0,
              f"{info['raw_lo']:.4f}..{info['raw_hi']:.4f}")

# genes already fill the colormap, so they must NOT be rescaled
for name in genes[:3]:
    png, info = cache.render(SAMPLE, 14, "gene", name)
    check(f"gene '{name}' is not rescaled",
          not info["normalized"] and info["scale"] == 1.0,
          f"normalized={info['normalized']} scale={info['scale']}")

# --------------------------------------------------------------------------
# error paths
# --------------------------------------------------------------------------
try:
    cache.render(SAMPLE, 14, "gene", "NOT_A_REAL_GENE")
    check("unknown gene raises KeyError", False)
except KeyError:
    check("unknown gene raises KeyError", True)

try:
    cache.render(SAMPLE, 99, "gene", "VIM")
    check("unknown resolution raises ValueError", False)
except ValueError as exc:
    check("unknown resolution raises ValueError", True, str(exc)[:44])

try:
    cache.meta("NOT_A_SAMPLE")
    check("unknown sample raises KeyError", False)
except KeyError:
    check("unknown sample raises KeyError", True)

try:
    cache.render(SAMPLE, 14, "celltype", "NotACellType")
    check("unknown cell type raises KeyError", False)
except KeyError:
    check("unknown cell type raises KeyError", True)

# --------------------------------------------------------------------------
# raster semantics: expression 0 is distinguishable from "no patch"
# --------------------------------------------------------------------------
layer = cache.handle(SAMPLE).layer(14)
flat, W, H, _xs, _ys = layer.grid()
mask = layer.mask()
check("the raster is (H, W) and matches the grid", mask.shape == (H, W),
      (mask.shape, (H, W)))
# every patch must land on a cell the mask marks present
check("every patch lands on a cell the mask marks present",
      bool(mask.reshape(-1)[flat].all()))
# expression 0 and "no patch" stay separable: the mask is the discriminator
check("the two concepts are separable via the mask",
      bool(mask.any()) and bool((~mask).any()),
      f"{100*mask.mean():.1f}% of cells hold a patch")

print()
print(f"stats: {cache.stats}")
if failures:
    print(f"{len(failures)} FAILURE(S): " + "; ".join(failures))
    sys.exit(1)
print("ALL SERVING CHECKS PASSED")


