# AURORA H&E + Gene Expression Viewer

An interactive viewer for AURORA high-resolution gene-expression predictions.
Pick a **sample**, a **gene** and a **resolution** and the heatmap is decoded on
demand from the original prediction files. **No images are pre-generated**, so a
new gene or a new sample costs nothing to add.

The viewer runs as a small web service (`space/`) that can be deployed to a free
Hugging Face Space while the predictions stay in the dataset repo.

---

## Repository layout

```
space/                   the viewer service (deploy this)
  app.py                 FastAPI routes
  viewer_core.py         decode + raster + PNG encoding, cache, no web framework
  aurora_h5.py           readers for both on-disk schemas
  samples.py             sample discovery (local dir or Hugging Face repo)
  he_images.py           H&E discovery (folder / Hub / pyramid / placeholder)
  build_he_pyramid.py    one-off: build a web pyramid from a rescaled TIFF
  static/                front end: two canvases, one shared camera
  test_reader.py         both schemas, delta arithmetic, block boundaries
  test_server.py         rendering, flagging, cell types, error paths
  test_hf_backend.py     Hugging Face backend (Hub stubbed)
  test_he_setup.py       the H&E lookup rules documented in space/README.md
  test_frontend.js       headless checks of the browser code

DATA_REFERENCE.md        what each dataset contains, per schema
HE_FORMAT.md             how to store the H&E images (measured recommendation)

index.html, assets/,     the previous static site (kept for reference); it
images/, README.md       required manually uploading PNGs per gene
gene-names.txt           legacy gene list (stale, see DATA_REFERENCE.md)
```

## Running it locally

From `space/`, with the prediction files next to it (the repository layout), no
configuration is needed — the sample directory is auto-detected:

```powershell
cd space
python -m uvicorn app:app --port 7860
```

then open <http://127.0.0.1:7860>. To use a different folder:

```powershell
cd space
$env:AURORA_BACKEND   = 'local'
$env:AURORA_LOCAL_DIR = 'C:\path\to\folder\with\*-highres-pred.h5'
$env:AURORA_HE_DIR    = 'C:\path\to\folder\with\H&E\images'   # optional
python -m uvicorn app:app --port 7860
```

Only unset variables are filled in, so an explicit setting always wins. If the
page reports no sample loaded, `GET /api/health` reports the directory actually in
use, which is normally all you need to see.

**The H&E images are the one thing you must supply yourself** — the prediction
files contain none, and the repo the old viewer pointed at is gone. Setting
`AURORA_HE_DIR` is covered step by step, with the accepted filenames and a
troubleshooting section, in
[space/README.md](space/README.md#setting-up-the-he-images).

## Deploying to a Hugging Face Space

> **A free account cannot do this.** Docker and Gradio Spaces now require a paid
> plan; only *Static* Spaces are free. See the hosting discussion below.

1. Create a Space with the **Docker** SDK and copy the contents of `space/` into it.
2. Set these variables in **Settings -> Variables and secrets**:

   | Name | Value |
   | --- | --- |
   | `AURORA_BACKEND` | `hf` |
   | `AURORA_HF_REPO` | `AURORAData/AURORA_prediction` |
   | `AURORA_HF_SUBFOLDER` | `TCGA_LUAD_website,predictions_7um_28um/TCGA_LUAD` |
   | `AURORA_HE_REPO` | `AURORAData/AURORA_prediction` |
   | `AURORA_HE_SUBFOLDER` | `HE_image_website` |
   | `AURORA_DEFAULT_SAMPLE` | `TCGA-44-5645-01B-04-BS4` |
   | `AURORA_DEFAULT_GENE` | `EPCAM` |
   | `AURORA_PREFETCH` | `1` — strongly recommended, see below |

   `AURORA_HF_SUBFOLDER` accepts **comma-separated folders**: the LUAD
   predictions are split, with the 112 um files in `TCGA_LUAD_website` and the
   7/28 um files in `predictions_7um_28um/TCGA_LUAD`, and the viewer merges them
   by sample id.

   Add `HF_TOKEN` as a **secret** only if a repo is private or gated.
   `AURORA_HF_REVISION` pins a branch/tag/commit; `AURORA_HF_BLOCK` tunes the
   range-request block size (default 8 MiB); `AURORA_CACHE_DIR` sets where
   prefetched samples are stored.
3. `AURORA_PREFETCH=1` caches each whole sample locally the first time it is
   opened. Slower for the very first view, much faster for every gene afterwards:
   a cold 7 um gene over HTTP costs ~40 MB of chunk reads.

`space/README.md` carries the Space card metadata (`sdk: docker`,
`app_port: 7860`).

## Controls

* **Sample** - every prediction file found in the configured location.
* **Gene** - type to search; the list is read from the file, so it always matches
  the predictions that exist (counts differ between samples: 3000 vs 3002).
* **Resolution** - `7 um` (14 px cells), `28 um` (56 px) or `112 um` (224 px).
* **Colormap** - turbo (the historical choice), viridis, or gray.

Pan by dragging, zoom with the wheel or a double-click to fit, `n`/`p` (or the
arrow keys) to step through genes. The URL carries `#sample=...&gene=...&res=...`
so a view can be linked or reloaded.

The H&E panel carries a **scale bar** that updates with the zoom (the slide is
0.5 µm per pixel, so the bar reads µm at fit zoom and switches to mm when zoomed
out). Inside the prediction raster, cells with no patch are painted **white**, so
the missing points look like the H&E's white glass; the dark matte remains only
around the raster itself.

---

## How the data is read

The prediction datasets are compressed with a third-party **zstd HDF5 filter
(id 32015)**, sometimes in a `shuffle -> zstd` chain. `h5py` ships no plugin for
it, so a normal read fails with a misleading

```
OSError: Can't synchronously read data (can't open directory)
```

`aurora_h5.py` ignores HDF5's filter pipeline for those datasets: it reads raw
chunks with `Dataset.id.read_direct_chunk()` and reverses the `shuffle -> zstd`
chain in numpy. It also implements the `spatial_gene_delta_v1` encoding, which is
a uint8 modulo-256 first difference along the spatially sorted rows, restarted
every 16384 rows. Decoding is therefore a uint8 cumulative sum per block, and a
partial read must be widened back to its enclosing block - the tests in
`test_server.py` check exactly that.

Layout facts verified against real files:

| Fact | Value |
| --- | --- |
| Schema | `spatial_gene_delta_v1`, uint8 modulo-256 deltas, 16384-row blocks |
| Genes | 3000-3002, depending on the sample |
| Cell types | 9-10, with sample-specific names |
| Resolutions | 7 um (14 px cells), 28 um (56 px), 112 um (224 px) |
| `patch_loc_*` order | **(y, x)** - see below |
| Coordinate space | the 0.5 um/px rescaled image (the sample here is 58152 x 9960) |
| Colormap | turbo, per-gene normalisation, as in the old PNGs |

### Three things worth knowing

**`patch_loc_*` is `(y, x)`, not `(x, y)`.** Column 0 is vertical, column 1
horizontal. Settled empirically: column 1 is monotonic through the stored file
(the writer sorts by the primary axis), and scoring the prediction footprint
against tissue in a real rescaled H&E gives IoU 0.79 for (y, x) versus 0.19 or
less for (x, y). Coordinates are used **exactly as stored** - a patch may begin at
a negative coordinate, in which case only its visible part is drawn. Some comments
in the producer code (`plot.py`, `predict_highres_TCGA_optimized.py`) claim the
opposite; they are stale. Details in [HE_FORMAT.md](HE_FORMAT.md).

**The 28 um level is not an average of the 7 um level.** They are independent
predictions, so the coarse level cannot be dropped and derived on the fly.

**The number of genes and the cell-type names differ between samples.** Read them
from the file; nothing in the viewer hard-codes them.

### Why the current file format is fine to keep

The corpus does not need regenerating for this viewer to work. Measured on a
real sample: the payload compresses to about **0.51 bytes per gene-patch value**
(4.06 bits), which is close to the lossless floor for this signal - the field is
smooth (neighbouring cells differ by ~3 of 255 codes) but not low-entropy, so
zstd only reaches ~2.1x. Re-arranging the data does not shrink it.

There is one real inefficiency, and it only matters over a network:
`exp_pred_14` is chunked `(256 channels, 16384 rows)`, so the *channel* axis is
split while every chunk spans the full spatial axis. Reading **one gene's**
complete row therefore touches 24 chunks - about **36-50 MB of compressed data
for 0.38 MB of useful values**. Locally that is invisible (a 7 um gene decodes in
~0.5 s). Over HTTP it is the dominant cost, which is why `AURORA_PREFETCH=1`
exists.

If bandwidth ever becomes the bottleneck, the fix is a derived copy of the
corpus re-chunked one gene per chunk (~0.18 MB per gene, same total size), which
would turn each gene into a single small fetch. That is an additive build step,
not a change to the predictions.

### Values and colour

Expression targets are per-sample, per-gene **min-max normalised** during
training, so the natural colour coordinate is `uint8 / 255` - which is what the
viewer uses, matching the historical turbo plots. The stored `exp_max_*` arrays
are the quantisation maxima, needed only to recover absolute values (range
`[0, ~0.71]`); the viewer does not display them.

`gene-names.txt` (3026 entries) is stale: it also lists the 10 cell types and
marker genes such as `MKI67` and `PDCD1` that are not in the predictions. The
gene picker loads its list from `space/static/LUAD_gene.csv`, and any gene that a
loaded sample's H5 does not contain is shown dimmed and cannot be selected — so
the viewer never offers a gene it cannot render.

---

## Tests

```powershell
cd space
python test_server.py        # reader + serving layer, against a real H5
python test_reader.py     # delta encoding, partial reads
python -m pytest -q          # if you prefer
node test_frontend.js        # browser code, headless
python test_hf_backend.py    # Hugging Face backend, Hub stubbed
```

`test_server.py` and `test_reader.py` need a real `*-highres-pred.h5`; point
`AURORA_LOCAL_DIR` at the folder containing one, or drop it next to the tests.

