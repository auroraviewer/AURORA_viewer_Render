"""
AURORA gene-expression viewer - HTTP service.

Serves gene and cell-type maps for any (sample, resolution) straight from the
prediction H5 files, plus the H&E image for a sample so the two can be shown side
by side.  Nothing is pre-rendered: a view is decoded on request and the browser
applies the colormap.

Run locally:

    uvicorn app:app --port 7860

That is enough when the prediction files sit next to this folder (the repository
layout): the sample directory is then auto-detected.  Otherwise point it at them
explicitly, or use the Hub:

    set AURORA_BACKEND=local
    set AURORA_LOCAL_DIR=C:\\path\\to\\prediction\\h5

    set AURORA_BACKEND=hf
    set AURORA_HF_REPO=AURORAData/AURORA_prediction
    set AURORA_HF_SUBFOLDER=predictions_7um_28um/TCGA_LUAD
    set AURORA_PREFETCH=1

Serving logic lives in ``viewer_core.py`` and is testable without FastAPI; the
on-disk formats are documented in ``../DATA_REFERENCE.md``.
"""

from __future__ import annotations

import hashlib
import logging
import os

import numpy as np
from fastapi import FastAPI, HTTPException, Query, Request, Response
from fastapi.responses import FileResponse, JSONResponse

from he_images import HEImageStore
from samples import SampleStore
from viewer_core import TARGET_WIDTH, ViewerCache

HERE = os.path.dirname(os.path.abspath(__file__))
STATIC_DIR = os.path.join(HERE, "static")
REPO_ROOT = os.path.dirname(HERE)

# The sample the viewer opens on load.  Overridable so a deployment can point at
# a different specimen without touching the frontend.
DEFAULT_SAMPLE = os.environ.get("AURORA_DEFAULT_SAMPLE", "TCGA-44-5645-01B-04-BS4")

# The gene the viewer preselects.  Only used when the loaded sample actually has
# it; otherwise the first available gene is used.
DEFAULT_GENE = os.environ.get("AURORA_DEFAULT_GENE", "EPCAM")

log = logging.getLogger("aurora")


def _autodetect_dirs() -> None:
    """
    Fill in the usual locations when nothing has been configured.

    The prediction files and the rescaled H&E normally sit next to ``space/``, so
    ``uvicorn app:app --port 7860`` should work with no environment at all.  Only
    set a variable when it is unset: an explicit setting always wins, and a bad
    one should fail loudly rather than be silently replaced.
    """
    if os.environ.get("AURORA_BACKEND", "local").lower() != "local":
        return
    for name, default in (("AURORA_LOCAL_DIR", REPO_ROOT),
                          ("AURORA_HE_DIR", REPO_ROOT)):
        if not os.environ.get(name) and os.path.isdir(default):
            os.environ[name] = default


_autodetect_dirs()

app = FastAPI(title="AURORA gene expression viewer", version="3.0")
store = SampleStore()
cache = ViewerCache(store=store)
he_store = HEImageStore()


@app.on_event("startup")
def _report_configuration():
    """
    Say out loud where the samples were looked for and how many were found.

    Without this, a misconfigured directory produces a server that starts fine and
    an empty viewer, which is a confusing way to fail.
    """
    target = store.local_dir if store.backend == "local" else store.repo_id
    try:
        n = len(store.samples())
        log.warning("AURORA viewer: backend=%s  %s=%s  samples=%d  H&E=%d",
                    store.backend,
                    "dir" if store.backend == "local" else "repo", target, n,
                    he_store.status()["indexed"])
        if n == 0:
            log.warning("AURORA viewer: NO SAMPLES FOUND in %s - set "
                        "AURORA_LOCAL_DIR to the folder holding the .h5 files, or "
                        "check GET /api/health", target)
    except Exception as exc:
        log.warning("AURORA viewer: cannot list samples from %s: %s", target, exc)


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def _meta(sample_id: str) -> dict:
    try:
        return cache.meta(sample_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"unknown sample: {sample_id}")
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"{type(exc).__name__}: {exc}")


def _render(kind: str, sample: str, name: str, cell_px: int, width: int):
    try:
        return cache.render(sample, cell_px, kind, name, target_w=width)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=f"unknown {kind}: {exc}")
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"{type(exc).__name__}: {exc}")


def _image_response(request: Request, png: bytes, info: dict, tag: str):
    # ETag is the CONTENT hash, not the URL.  An earlier version keyed the ETag
    # on sample|cell|gene|width, so after the lattice code was fixed the browser
    # kept serving the stale PNG it had cached (max-age was 86400) for any gene
    # viewed before the fix - which looked exactly like "EPCAM is broken but
    # TMSB4X is fine".  Hashing the bytes makes any change to the render produce
    # a new ETag, and a short max-age bounds how long a stale copy can live.
    etag = '"' + hashlib.sha256(png).hexdigest() + '"'
    headers = {
        "Cache-Control": "public, max-age=300",
        "ETag": etag,
        "X-Render-Ms": str(int((info.get("raster_seconds", 0)
                                + info.get("png_seconds", 0)) * 1000)),
        "X-Cache": "hit" if info.get("cached") else "miss",
    }
    if request.headers.get("if-none-match") == etag:
        return Response(status_code=304, headers=headers)
    return Response(content=png, media_type="image/png", headers=headers)


# --------------------------------------------------------------------------
# API
# --------------------------------------------------------------------------
@app.get("/api/health")
def health():
    """
    Status, including **where the server looked for samples**.

    That last part is the point: "no sample is loaded" is nearly always a path
    problem, and the directory actually in use is the first thing to check.
    """
    try:
        found = len(store.samples())
        listing_error = None
    except Exception as exc:
        found = 0
        listing_error = f"{type(exc).__name__}: {exc}"
    return {
        "ok": True,
        "backend": store.backend,
        "samples_loaded": found,
        "sample_dir": store.local_dir,
        "sample_dir_exists": bool(store.local_dir and os.path.isdir(store.local_dir)),
        "repo": store.repo_id,
        "subfolder": store.subfolder,
        "listing_error": listing_error,
        "open_sample": cache.open_sample,
        "he_images": he_store.status(),
        "prefetch": store.prefetch,
        "stats": cache.stats,
    }


@app.get("/api/samples")
def list_samples(q: str | None = None):
    try:
        items = store.samples()
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"cannot list samples: {exc}")
    out = [{"sample_id": s.sample_id, "size": s.size} for s in items]
    if q:
        ql = q.lower()
        out = [s for s in out if ql in s["sample_id"].lower()]
    ids = [s["sample_id"] for s in out]
    # the sample the UI should open first; falls back to the first listed one
    default = DEFAULT_SAMPLE if DEFAULT_SAMPLE in ids else (ids[0] if ids else None)
    return {"count": len(out), "samples": out, "default": default,
            "default_gene": DEFAULT_GENE}


@app.post("/api/refresh")
def refresh():
    try:
        n = store.refresh()
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc))
    return {"ok": True, "count": n}


@app.post("/api/samples/{sample_id}/prefetch")
def prefetch(sample_id: str):
    """
    Start (or report) a background download of a whole sample.

    Reading one 7 um gene over HTTP needs ~40 MB of chunks because the file is
    chunked by channel, so a cold gene takes tens of seconds.  Downloading the
    file once makes every later view local.
    """
    try:
        store.get(sample_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"unknown sample: {sample_id}")
    return store.prefetch_status(sample_id)


@app.get("/api/samples/{sample_id}/meta")
def sample_meta(sample_id: str):
    m = _meta(sample_id)
    m["he"] = he_store.describe(sample_id)
    return JSONResponse(m, headers={"Cache-Control": "private, max-age=300"})


@app.get("/api/samples/{sample_id}/genes.txt")
def sample_genes(sample_id: str):
    m = _meta(sample_id)
    return Response("\n".join(m["genes"]), media_type="text/plain",
                    headers={"Cache-Control": "private, max-age=600"})


@app.get("/api/gene.png")
def gene_png(
    request: Request,
    sample: str = Query(...),
    gene: str = Query(...),
    cell: int = Query(14, description="cell size in rescaled-image px: 14, 56 or 224"),
    width: int = Query(TARGET_WIDTH, ge=64, le=4096),
):
    png, info = _render("gene", sample, gene, cell, width)
    return _image_response(request, png, info, f"{sample}|{cell}|gene|{gene}|{width}")


@app.get("/api/celltype.png")
def celltype_png(
    request: Request,
    sample: str = Query(...),
    celltype: str = Query(...),
    cell: int = Query(14),
    width: int = Query(TARGET_WIDTH, ge=64, le=4096),
):
    png, info = _render("celltype", sample, celltype, cell, width)
    return _image_response(request, png, info,
                           f"{sample}|{cell}|celltype|{celltype}|{width}")


@app.get("/api/view/stats")
def api_view_stats(sample: str, name: str, kind: str = "gene", cell: int = 14,
                   width: int = TARGET_WIDTH):
    """
    Value range for a view, plus where its raster sits in image-pixel space.

    The client needs `origin` / `extent` to place the raster correctly: the
    prediction covers only part of the slide, and a layer may start at a negative
    coordinate, so stretching it over the whole image would misalign it.
    """
    _meta(sample)
    h = cache.handle(sample)
    try:
        layer = h.layer(cell)
    except KeyError as exc:
        raise HTTPException(status_code=400, detail=str(exc))

    if kind == "gene":
        idx = h.gene_index()
        values = layer.exp_values([idx[name]])[:, 0] if name in idx else None
    elif kind == "celltype":
        idx = h.cell_type_index()
        values = layer.proportions([idx[name]])[:, 0] if name in idx else None
    else:
        raise HTTPException(status_code=400, detail=f"unknown kind {kind}")
    if values is None:
        raise HTTPException(status_code=404, detail=f"unknown {kind}: {name}")

    v = values

    # Render (from cache) so this reports the *same* bounds the PNG was drawn
    # with.  The client labels the legend from these, and recomputing them here
    # would let the label drift away from the colours.
    png, info = _render(kind, sample, name, cell, width)
    return {
        "sample": sample, "kind": kind, "name": name,
        "cell_px": cell, "microns": layer.microns,
        "raster": info["raster"],
        "native_grid": info["native_grid"],
        "origin": info["origin"],
        "extent": info["extent"],
        "world_bounds": info["world_bounds"],
        "mode": info["mode"],
        # what the colour means: for cell types the ramp spans raw_lo..raw_hi of
        # the raw probability; for genes it spans 0..exp_max
        "raw_unit": info["raw_unit"],
        "raw_max": info["raw_max"],
        "raw_lo": info["raw_lo"],
        "raw_hi": info["raw_hi"],
        "normalized": info["normalized"],
        "scale": info["scale"],
        "bytes": info["bytes"],
        "n": int(values.size),
        "min": float(v.min()) if v.size else None,
        "max": float(v.max()) if v.size else None,
        "mean": float(v.mean()) if v.size else None,
    }


@app.get("/api/he/{sample_id}")
def he_image(sample_id: str, max_width: int = Query(4096, ge=256, le=20000)):
    # prefer a pre-built pyramid: serving it is a file read, whereas the raw
    # slide costs ~13 s to decode on a cold cache.  Look on local disk first,
    # then in the HF repo.
    pyr = he_store.pyramid_dir(sample_id)
    if pyr:
        resp = he_store.serve_pyramid(pyr, max_width)
        if resp is not None:
            return resp
    resp = he_store.serve_hf_pyramid(sample_id, max_width)
    if resp is not None:
        return resp
    path = he_store.path_for(sample_id)
    if path is None:
        raise HTTPException(status_code=404, detail=f"no H&E image for {sample_id}")
    return he_store.serve(path, max_width)


@app.get("/api/he/{sample_id}/align")
def he_align(sample_id: str, cell: int = 14, white: int = Query(220, ge=100, le=255)):
    """
    Check that an H&E image and the prediction extent describe the same region.

    Scales the image to the prediction grid and compares its tissue mask with the
    prediction footprint.  A correct pairing scores a high IoU; a mismatch shows
    up immediately, which matters because the H&E images are supplied separately
    from the prediction files.
    """
    import numpy as np
    from PIL import Image

    from viewer_core import align_report

    path = he_store.path_for(sample_id)
    if path is None:
        raise HTTPException(status_code=404, detail=f"no H&E image for {sample_id}")
    _meta(sample_id)
    h = cache.handle(sample_id)
    try:
        layer = h.layer(cell)
    except KeyError as exc:
        raise HTTPException(status_code=400, detail=str(exc))
    try:
        img = Image.open(path)
        img.load()
    except Exception as exc:
        raise HTTPException(status_code=500, detail=f"cannot read image: {exc}")

    report = align_report(img, layer, white=white)
    report.update({
        "sample": sample_id,
        "cell_px": cell,
        "microns": layer.microns,
        "he_size": list(img.size),
        "extent": list(layer.extent()),
        "placeholder": he_store.is_placeholder(sample_id),
    })
    return report


@app.get("/api/cache/clear")
def clear_cache():
    cache.close()
    return {"ok": True}


# --------------------------------------------------------------------------
# static front end
# --------------------------------------------------------------------------
@app.get("/")
def index():
    idx = os.path.join(STATIC_DIR, "index.html")
    if not os.path.exists(idx):
        return JSONResponse({"detail": "static/index.html is missing"}, status_code=500)
    return FileResponse(idx)


if os.path.isdir(STATIC_DIR):
    from fastapi.staticfiles import StaticFiles

    app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")



