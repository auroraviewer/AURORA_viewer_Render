# AURORA Viewer — Setup & Run Instructions

Interactive H&E + gene-expression viewer that reads AURORA prediction files
(`*-highres-pred.h5` and `*-pred.h5`) directly, with no pre-rendering. It
decodes the HDF5 chunks on demand, serves them as PNGs, and the browser applies
the colormap. No GPU is required.

---

## 1. What you need

- **Python 3.10+** (3.11 or 3.12 recommended).
- Your **prediction H5 files** — either on disk, or in a Hugging Face dataset repo.
- Your **H&E images** (optional but recommended) — already rescaled to the
  0.5 µm/px size the predictions were made against (e.g. `<id>_rescaled.tif`).

## 2. Files in this package

| path | purpose |
| --- | --- |
| `app.py` | FastAPI HTTP service (routes, entry point) |
| `samples.py` | finds and opens prediction files (local dir or HF repo) |
| `aurora_h5.py` | reads/decodes the prediction HDF5 (incl. the zstd chunks) |
| `viewer_core.py` | rasterises a gene/cell-type view to a PNG |
| `he_images.py` | finds and serves H&E images / pyramids |
| `build_he_pyramid.py` | builds one H&E pyramid from a rescaled image |
| `batch_he_pyramid.py` | **batch**-converts a folder of H&E images into pyramids |
| `download_samples.py` | mirrors prediction H5 files from a HF dataset to disk |
| `requirements.txt` | Python dependencies |
| `Dockerfile` | container image (for HF Spaces / Render) |
| `static/` | the browser UI (`index.html`, `app.js`, `colormaps.js`, `LUAD_gene.csv`, `logo.png`) |
| `test_*.py` / `test_*.js` | self-checks (optional) |
| `render.yaml` | Render blueprint (only if you deploy on Render) |

You do **not** need the large `.h5` / `.tif` data files — those are yours.

---

## 3. Install

```bash
cd space
python -m pip install -r requirements.txt
```

## 4. Run locally (files on disk)

Put your prediction files in one folder, e.g. `D:\aurora`:

```
D:\aurora\
    TCGA-...-highres-pred.h5     # 7 um + 28 um
    TCGA-...-pred.h5             # 112 um
    TCGA-..._rescaled.tif        # (optional) H&E
```

```bash
# Windows PowerShell
$env:AURORA_BACKEND = "local"
$env:AURORA_LOCAL_DIR = "D:\aurora"
$env:AURORA_HE_DIR    = "D:\aurora"     # optional, for H&E
python -m uvicorn app:app --port 7860
```

Open **http://127.0.0.1:7860/** in a browser.

> If `AURORA_LOCAL_DIR` is not set, the app auto-detects the parent folder
> (the repo root) — useful if you just drop the H5 files next to `space/`.

## 5. H&E images (make the side-by-side view fast)

The viewer can serve a plain rescaled image, but a **pyramid** (halved WebP
levels) loads much faster. Convert one image:

```bash
python build_he_pyramid.py --tif D:\aurora\TCGA-..._rescaled.tif --out D:\aurora\TCGA-..._pyramid
```

Or convert a whole folder at once:

```bash
python batch_he_pyramid.py --src D:\he_images --out D:\aurora
# add --resize 58152x9960 if the source is not already at 0.5 um/px
```

This writes `<id>_pyramid/` (with `L*.webp` + `manifest.json`) next to the data.
Point `AURORA_HE_DIR` (local) or `AURORA_HE_SUBFOLDER` (HF) at that directory.

Whole-slide formats (`.svs`, `.ndpi`, `.mrxs`) are not decoded — rescale those
first.

## 6. Stream from a Hugging Face dataset (no local files)

If the predictions live in a HF **dataset** repo — this is the configuration used
for the LUAD website:

```bash
$env:AURORA_BACKEND = "hf"
$env:AURORA_HF_REPO = "AURORAData/AURORA_prediction"
# the LUAD predictions are split across two folders (comma-separated):
#   112 um      -> TCGA_LUAD_website
#   7 um / 28um -> predictions_7um_28um/TCGA_LUAD
$env:AURORA_HF_SUBFOLDER = "TCGA_LUAD_website,predictions_7um_28um/TCGA_LUAD"
$env:AURORA_HE_REPO = "AURORAData/AURORA_prediction"
$env:AURORA_HE_SUBFOLDER = "HE_image_website"
$env:AURORA_DEFAULT_SAMPLE = "TCGA-44-5645-01B-04-BS4"
$env:AURORA_PREFETCH = "1"          # download each sample once, then serve locally
$env:HF_TOKEN = "hf_..."            # only for private/gated repos
python -m uvicorn app:app --port 7860
```

`AURORA_HF_SUBFOLDER` accepts **comma-separated folders**; files are merged by
sample id, so a sample whose 112 um and 7/28 um files live in different folders
still appears once, with all three resolutions.

Set `AURORA_PREFETCH=1` so the first view downloads the sample to disk and every
later view is fast; a LUAD sample is a few hundred MB to ~1.5 GB.

## 7. Deploy

- **Render (free/cheap):** push this repo to GitHub and connect it with a
  **Blueprint** — `render.yaml` does the rest. Set the variables above in the
  dashboard. Free instances sleep when idle (cold start on the next visit).
- **Hugging Face Space (Docker SDK):** copy the contents of `space/` into a
  Docker Space, set the same variables, and set `sdk: docker`, `app_port: 7860`
  in the Space `README.md`.

## 8. Verify

```bash
python test_reader.py       # H5 decode vs the reference reader
python test_server.py       # serving / rendering
python test_he_setup.py     # H&E discovery claims
node test_frontend.js       # UI logic (needs node)
node test_boot.js           # UI boots without throwing (needs node)
```

## 9. Troubleshooting

- **"No sample is loaded"** — check `http://127.0.0.1:7860/api/health`: it shows
  which directory/repo was searched and how many samples were found.
- **Changes don't appear** — restart the server (Python caches imported modules)
  and hard-refresh the browser (`Ctrl+Shift+R`).
- **Port busy** — another process holds 7860; kill it or use `--port 8000`.
