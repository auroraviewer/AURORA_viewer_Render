# H&E images: format and coordinate convention

Two things the user asked about, both settled by measurement against a real
matched pair (`TCGA-JL-A3YX-01A-01-TS1`).

The source image is `*_rescaled.tif`: **58152 x 9960** (579 MP), RGB, LZW, 736 MB
on disk, ~9 s to decode. That decode time is the reason a container format must
never be served per request.

---

## 1. Coordinate convention: `patch_loc` is (y, x)

`patch_loc_14` / `_56` / `_224` store each patch's top-left corner as
**column 0 = y (vertical), column 1 = x (horizontal)** — i.e. the raw order is
(y, x), not (x, y).

Settled empirically, twice over:

* **Against the stored order.** Column 1 is monotonically non-decreasing through
  the file (0 violations), column 0 is not. The writer sorts spatially by the
  primary axis, so column 1 is the outer/slow axis.
* **Against the H&E.** Score the prediction footprint against the tissue in the
  real rescaled image under every plausible mapping:

  | mapping | IoU |
  | --- | --- |
  | **(y, x)**, no flip | **0.79** |
  | (y, x) with both axes flipped | 0.54 |
  | (x, y), any flip | ≤ 0.19 |

  With the correct convention the patch spans recover the H&E dimensions exactly:
  x 26..58140 in a 58152-wide image, y 10..9950 in a 9960-tall image, span ratio
  5.84 against an H&E aspect of 5.84.

**Coordinates are used exactly as stored — no shift is applied.** An earlier
build added a constant offset to the 112 µm positions to bring them into range;
that misaligned the prediction against the H&E and has been removed.

## 2. Negative coordinates are clipped, not shifted

A patch may start at a negative coordinate, meaning part of it lies outside the
rescaled image. In this sample 17 of 4451 112 µm patches do (0.38%): x from −44,
y from −60.

Handling: the raster's origin is the **minimum** stored coordinate, so a patch
starting at −44 simply lands at the raster edge and only its visible portion is
represented. The renderer is told that origin via `Layer.origin()`, so the raster
is placed at the right offset in the image. Nothing is shifted, which is what
keeps the overlay aligned.

`Layer.bounds()` returns the clipped `(x0, y0, x1, y1)` for callers that need it.

## 3. Storage on the Hub: a WebP pyramid, ~8 MB per sample

**Recommendation: ship a pre-built image pyramid of WebP files, not a TIFF.**

Measured on the real image (`space/build_he_pyramid.py`):

| level | size | WebP q85 |
| --- | --- | --- |
| 0 | 58152 x 9960 | *exceeds WebP's 16383 px limit* |
| 1 | 29076 x 4980 | *exceeds the limit* |
| **2** | **14538 x 2490** | **5.9 MB** |
| **3** | **7269 x 1245** | **1.7 MB** |
| **4** | **3634 x 622** | **0.4 MB** |
| | **total** | **8.1 MB** |

Why this is the right trade:

* **Nothing pre-generated per gene.** One pyramid per sample serves every gene and
  every cell type, because the colour is applied in the browser.
* **~8 MB per sample** versus 736 MB for the source TIFF — about 90x smaller, and
  ~4 GB for 500 samples.
* **Detail is sufficient.** At the finest WebP level one 7 um cell is 4.6 image
  pixels wide, so the heatmap stays smooth well past any realistic zoom. Beyond
  that the browser interpolates rather than revealing more data than exists.
* **No codec limitation fights you.** WebP cannot encode above 16383 px per side.
  Trying to include levels 0 and 1 produced 762 MB of PNG tiles and a 770 MB
  pyramid — the opposite of fast loading. Capping at level 2 avoids the problem
  entirely. (JPEG can encode the full size: level 0 is 69 MB, level 1 22 MB, for
  a ~100 MB pyramid. Only worth it if pixel-level zoom is genuinely needed.)

### Answering "should I rescale by another 0.1?"

Not necessary, and it would cost detail. The pyramid already solves loading
speed: the browser fetches only the visible tiles at the level it needs. A 0.1x
rescale (5815 x 996) is simply the overview level that the pyramid already
contains — as its own file it would be ~1.3 MB, useful as a first-paint
thumbnail, but it does not replace the pyramid.

### Suggested layout on the Hub

The viewer already prefers a pyramid when it finds one. Put it next to the
predictions (or point `AURORA_HE_PYRAMID_DIR` at it), named either
`<sample_id>/`, `he/<sample_id>/` or `<sample_id>_pyramid/`:

```
<sample_id>/L2.webp      finest level (14538 x 2490)  -- 5.9 MB
<sample_id>/L3.webp                                    -- 1.7 MB
<sample_id>/L4.webp                                    -- 0.4 MB
```

Levels are picked by width, so the names' numbers only need to be ordered. A
request for a 2500 px wide image is then served straight from `L3.webp` in **under
half a second**, against ~13 s to decode the source TIFF — measured.

`space/build_he_pyramid.py` builds this from a rescaled TIFF:

```bash
python build_he_pyramid.py --tif <sample>_rescaled.tif --out <sample>_pyramid --levels 5
# then delete the L0*/L1* PNG tiles it had to emit for the oversized levels
```

### If you want to keep it simple

A single **14538 x 2490 WebP** (5.9 MB) is already enough for a good viewer:
whole-slide overview plus useful zoom, one file per sample, no tiling logic. The
extra 2.2 MB for levels 3 and 4 buys faster zoomed-out rendering and quicker
first paint.
