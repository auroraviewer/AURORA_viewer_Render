---
sdk: docker
app_port: 7860
title: AURORA H&E + Gene Expression Viewer
---

# AURORA H&E + Gene Expression Viewer

Interactive viewer for AURORA high-resolution gene-expression predictions.
Pick a **sample**, a **gene** and a **resolution**; the
heatmap is decoded on demand from the original prediction files and shown next to
the H&E image with **synchronised pan and zoom**. No images are pre-generated.

---

## Status of the five reported problems

| # | Report | Status |
| --- | --- | --- |
| 1 | 28 µm broken, bottom part a mess | **Root-caused and fixed - it was a decoder bug, not the data.** A partially written final HDF5 chunk keeps the *physical* row stride, and the decoder reshaped it by the *logical* extent, so every gene after the first in the final block read padding as values. `test_partial_chunk.py` verifies the decode against the native h5py read. The prediction files are correct as written. |
| 2 | H&E side by side + synced zoom | **Implemented and verified** against a real matched pair. `AURORA_HE_DIR` finds `<sample_id>_rescaled.tif` automatically; alignment scores IoU 0.79 (7 µm) / 0.78 (28 µm) / 0.98 (112 µm). See [HE_FORMAT.md](../HE_FORMAT.md) for the storage recommendation. |
| 3 | 112 µm predictions | **Implemented.** Different schema, handled by the same reader; all three resolutions appear as buttons. See [DATA_REFERENCE.md](../DATA_REFERENCE.md). |
| 4 | Cell-type visualisation | **Removed from the UI.** Decoding the cell-type proportions needs a softmax over the whole patch array and overflowed the 512 MB limit of a free hosting instance. The backend endpoints remain, but the picker is gone so nothing decodes that data. |
| 5 | No list of genes / cell types | **Implemented** — searchable gene picker read from the file. Gene counts differ between samples (3000/3002), so nothing is hard-coded. |

A rendering bug was also found and fixed: prediction code `0` (real zero
expression) was drawn as background, which made low-expression regions look
empty. Values and "no patch" are now tracked separately through an explicit
presence channel.

The coordinate convention was corrected too: **`patch_loc` is (y, x)**, and
coordinates are used exactly as stored — an earlier build applied a constant
offset to the 112 µm positions, which misaligned the overlay. Patches starting at
a negative coordinate are clipped, not shifted.

## Run it

From `space/`, this is enough when the prediction files sit next to it — which is
the repository layout — because the sample directory is then auto-detected:

```powershell
cd space
python -m uvicorn app:app --port 7860
```

then open <http://127.0.0.1:7860>.

To point it somewhere else:

```powershell
cd space
$env:AURORA_BACKEND   = 'local'
$env:AURORA_LOCAL_DIR = 'C:\path\to\folder\with\prediction\h5'
$env:AURORA_HE_DIR    = 'C:\path\to\folder\with\H&E\images'   # optional
python -m uvicorn app:app --port 7860
```

Only unset variables are filled in, so anything you set explicitly still wins.
If a page says no sample is loaded, check `GET /api/health` — it reports the
directory actually in use:

```json
{ "backend": "local", "samples_loaded": 1,
  "sample_dir": "...\\AURORA_viewer-gh-pages", "sample_dir_exists": true,
  "listing_error": null }
```

Setting the H&E path is covered in detail under
[Setting up the H&E images](#setting-up-the-he-images) below.

Against the Hub (see the hosting note in the repository README):

```powershell
$env:AURORA_BACKEND    = 'hf'
$env:AURORA_HF_REPO    = 'AURORAData/AURORA_prediction'
$env:AURORA_HF_SUBFOLDER = 'predictions_7um_28um/TCGA_LUAD'
$env:AURORA_PREFETCH   = '1'      # strongly recommended, see below
python -m uvicorn app:app --port 7860
```

### Why `AURORA_PREFETCH=1` matters

`exp_pred_14` is chunked `(256 channels, 16384 spatial rows)` — the *channel* axis
is split and every chunk spans the full spatial axis — so reading **one** gene's
row touches 24 chunks. Measured over a normal connection:

| | |
| --- | --- |
| open a sample | 2.7 s |
| first 7 µm gene (24 chunks) | **50.8 s** |
| second 7 µm gene | 0.2 s |

Locally the same operations take 0.1 s and 0.01 s. So streaming is only usable
with a warm cache; `AURORA_PREFETCH=1` downloads the sample once (a few tens of
seconds) and makes every view local.

### Setting up the H&E images

The prediction files contain **no images**, and the repo the old viewer's README
pointed at (`AURORAData/prediction_plots`) is gone. You must supply them, and this
is the only thing the viewer cannot work out for itself.

You need the **rescaled 0.5 µm/px image** — the same one the predictions were made
against, normally named `<sample_id>_rescaled.tif`. Using the original whole-slide
`.svs` will not work (Pillow cannot read it, and it is at the wrong scale anyway).

#### Option A — a folder on disk (simplest)

```powershell
cd space
$env:AURORA_HE_DIR = 'D:\aurora\he'      # a folder you control
python -m uvicorn app:app --port 7860
```

`AURORA_HE_DIR` is scanned **recursively**, so subfolders are fine. Point it at a
parent directory if you like:

```
D:\aurora\he\
  TCGA-JL-A3YX-01A-01-TS1_rescaled.tif
  TCGA-05-4244-01A-01-BS1_rescaled.tif
  another_batch\TCGA-44-2668-01A-01-TS1_rescaled.tif
```

#### Filenames it recognises

The filename must be the sample id, optionally with a known suffix. All of these
resolve to sample `TCGA-JL-A3YX-01A-01-TS1` (`_rescaled` is the one to use):

| filename | |
| --- | --- |
| `<id>_rescaled.tif` | ✅ preferred — what the pipeline produces |
| `<id>-rescaled.tif`, `<id>_rescaled_image.tif` | ✅ |
| `<id>.png`, `<id>.jpg`, `<id>.tif`, `<id>.webp` | ✅ |
| `<id>_he.png`, `<id>-he.jpg`, `<id>_HE.tif`, `<id>_he_image.png` | ✅ |
| `<id>.svs`, `.ndpi`, `.mrxs` | ❌ whole-slide formats, not supported |
| `logo.png`, `zoomout_rest.png` | ❌ ignored (not specimen ids) |

Files whose name does not start like a specimen id are skipped, so a folder
holding viewer icons alongside slides is fine. `images/`, `assets/` and `static/`
subtrees are skipped entirely.

> If two files would resolve to the same sample id, whichever `os.walk` reaches
> first wins. Keep one image per sample.

#### Option B — straight from a Hugging Face dataset repo

Good when the images already live on the Hub and you do not want them on local
disk. Sampling is by range request, so the viewer downloads only what it serves.

```powershell
$env:AURORA_HE_REPO      = 'YourOrg/aurora-he'
$env:AURORA_HE_SUBFOLDER = 'rescaled'          # optional
# $env:HF_TOKEN          = 'hf_...'            # only for private/gated repos
```

It looks for the same names as above, in this order: `<id>`, `<id>_he`,
`<id>-he`, `<id>_HE`, `<id>_he_image`, each with every supported extension.

#### Option C — a pre-built pyramid (fastest, recommended for real use)

Decoding a 579 MP rescaled TIFF takes ~13 s, and that is paid on the first request
for each image. A pyramid turns it into a file read. See
[HE_FORMAT.md](../HE_FORMAT.md) for the measurements and the reasoning.

```powershell
# one-off, per sample
python build_he_pyramid.py --tif D:\aurora\he\<id>_rescaled.tif `
                           --out D:\aurora\he\<id>_pyramid --levels 5
# delete the L0*/L1* PNG tiles it had to emit for the oversized levels
```

The viewer finds a pyramid automatically if it sits next to the image as `<id>/`,
`he/<id>/` or `<id>_pyramid/`, or point `AURORA_HE_PYRAMID_DIR` at wherever they
live. Measured effect on the sample above: a 2500 px request went from **13.5 s to
441 ms**.

#### Option D — a placeholder, for demoing without real images

```powershell
$env:AURORA_HE_PLACEHOLDER = 'C:\path\to\any_slide.png'
```

Used for every sample. The UI labels it **"placeholder image"** so it can never be
mistaken for the real slide, and the side-by-side view can be exercised before the
real images are to hand.

#### All the variables

| Variable | Meaning |
| --- | --- |
| `AURORA_HE_DIR` | folder scanned recursively for `<sample_id>[suffix][ext]` |
| `AURORA_HE_REPO` | Hugging Face **dataset** repo to fetch images from |
| `AURORA_HE_SUBFOLDER` | folder inside that repo |
| `AURORA_HE_REVISION` | optional branch/tag/commit within that repo |
| `AURORA_HE_PYRAMID_DIR` | where pre-built pyramids live, if not beside the images |
| `AURORA_HE_PLACEHOLDER` | one image used for every sample; dev only |
| `HF_TOKEN` / `HUGGING_FACE_HUB_TOKEN` | only for private/gated repos |

#### Checking it worked

Start the server and open a sample. If the left panel says *"No H&E image for this
sample"*, the lookup failed — these tell you why:

```powershell
# what got indexed, and how many
curl http://127.0.0.1:7860/api/health

# is a specific sample matched, and is it a placeholder?
curl http://127.0.0.1:7860/api/samples/TCGA-JL-A3YX-01A-01-TS1/meta
#   -> "he": {"available": true, "placeholder": false,
#             "basename": "TCGA-JL-A3YX-01A-01-TS1_rescaled.tif", ...}

# does the image actually line up with the predictions?
curl "http://127.0.0.1:7860/api/he/TCGA-JL-A3YX-01A-01-TS1/align?cell=14"
```

The `align` score is the real check: it compares the tissue in the image against
the prediction footprint, so it catches a mismatched or wrongly scaled image that
would otherwise look plausible. Measured on the sample here: **IoU 0.79 (7 µm),
0.78 (28 µm), 0.98 (112 µm)**. If you get something below ~0.5, the pairing is
wrong — most often the image is a different sample, or a full-resolution slide
rather than the 0.5 µm/px rescale.

#### How images are served

Downscaled on the way out and cached: from a pyramid level it is a file read;
otherwise the source is decoded once (~13 s for a 579 MP TIFF) and re-encoded as
JPEG (a 4096 px PNG of a 40 MB source is 25 MB; JPEG is ~4.6 MB). The decoded
source and recent encodes are kept in memory, so repeat views are milliseconds.

## Layout

```
space/                     the viewer service
  app.py                   FastAPI routes
  viewer_core.py           decode + raster + PNG, no web framework
  aurora_h5.py             readers for both on-disk schemas
  samples.py               sample discovery, HF streaming, prefetch
  he_images.py             H&E discovery and serving
  build_he_pyramid.py      one-off: build a web pyramid from a rescaled TIFF
  download_samples.py      mirror samples from the Hub to a local folder
  static/                  front end: dual canvas, one shared camera
  test_reader.py           both schemas, delta arithmetic, defect detector
  test_server.py           rendering, flagging, cell types, error paths
  test_hf_backend.py       HF backend contract (Hub stubbed)
  test_he_setup.py         the H&E lookup rules documented above
  test_frontend.js         front-end behaviours, by inspection
  test_boot.js             runs app.js against a mock DOM and checks start-up
  test_partial_chunk.py    partial HDF5 chunk decode vs native h5py

DATA_REFERENCE.md          what each dataset contains, per schema
HE_FORMAT.md               how to store the H&E images, with measurements
gene-names.txt             legacy gene list (stale, see DATA_REFERENCE)
index.html, assets/...     the previous static site, kept for reference
```

## API

| Route | Purpose |
| --- | --- |
| `GET /` | the viewer |
| `GET /api/samples` | list sample ids |
| `GET /api/samples/{id}/meta` | genes, cell types, layers, extents, H&E availability |
| `POST /api/samples/{id}/prefetch` | start/report a background download |
| `GET /api/gene.png?sample=&gene=&cell=&width=` | gene raster |
| `GET /api/celltype.png?sample=&celltype=&cell=&width=` | cell-type raster |
| `GET /api/view/stats` | value range for a view |
| `GET /api/he/{id}` | the H&E image |
| `GET /api/he/{id}/align` | H&E vs prediction footprint IoU |
| `POST /api/refresh` · `GET /api/cache/clear` | rescan / drop caches |

`cell` is the cell size in rescaled-image pixels: `14` (7 µm), `56` (28 µm) or
`224` (112 µm).

Rasters are sent as RGB PNGs and the browser applies the colormap, so colour never
round-trips through the server. The channels are:

| channel | meaning |
| --- | --- |
| `R` | the value code, 0-255 |
| `B` | presence: 255 where a patch exists, 0 where there is none |

Presence is what makes "no patch" transparent: without it the client cannot tell
background from a legitimate expression of 0, and the raster is drawn as a
rectangle of near-black over the whole image. Values are drawn exactly as stored —
nothing is masked, hatched or annotated.

## Tests

```powershell
cd space
python test_reader.py       # readers, delta arithmetic, block boundaries
python test_server.py       # rendering, flagging, cell types, errors
python test_hf_backend.py   # HF backend contract
python test_he_setup.py     # the H&E lookup rules documented above
node test_frontend.js       # front-end behaviours, by inspection
node test_boot.js           # actually RUNS app.js against a mock DOM
python test_partial_chunk.py # partial HDF5 chunk decode vs native h5py
```

The Python suites expect a matched pair on disk. They default to
`TCGA-JL-A3YX-01A-01-TS1` in the repository root; override with any of:

| Variable | Meaning |
| --- | --- |
| `AURORA_TEST_SAMPLE` | the sample id to test |
| `AURORA_TEST_H5` | path to the `*-highres-pred.h5` (7/28 um) |
| `AURORA_TEST_H5_V1` | path to the `*-pred.h5` (112 um); its checks are skipped if absent |
| `AURORA_TEST_HE` | path to the rescaled H&E; its checks are skipped if absent |

A suite skips cleanly rather than failing when a file is missing, so
`test_reader.py` runs even without the 112 um file or the H&E.
