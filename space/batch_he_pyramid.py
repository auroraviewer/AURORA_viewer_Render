"""
Batch-convert H&E images into the pyramid format the viewer serves.

For every image it finds it writes::

    <sample_id>_pyramid/
        L2.webp          # halved WebP levels; L<k> is the native image / 2**k
        L3.webp
        L4.webp
        manifest.json    # {"size": [native_w, native_h], ...}

The viewer serves a pyramid directly (no re-decode on every request) and reads
``manifest.json`` for the native size, so the prediction and H&E stay aligned.
Levels wider/taller than 16383 px are skipped because WebP cannot encode them;
their *names* are still counted from the native image, so L2 is always a quarter
of the native size and the viewer can still reconstruct the native dimensions.

The input image must be at the **rescaled 0.5 um/px** size the predictions were
made against (the viewer does not rescale; it places the H&E by its native size).
Use ``--resize WxH`` to convert a source at a different resolution.

Only ``<sample>_rescaled.tif`` images are converted.  A matching mask file
(``<sample>_mask_rescaled.tif`` or ``<sample>_HistoSweep_mask_rescaled.tif``) is
**applied** rather than converted: everything outside the tissue is painted
white, matching the reference figures.  Disable with ``--no-mask``.

Whole-slide formats (.svs/.ndpi/.mrxs) are not decoded by Pillow and are skipped
with a warning - rescale those first (e.g. with the AURORA preprocessing).

Examples
--------
    # images are already rescaled
    python batch_he_pyramid.py --src ./he_images --out ./he_out

    # also resize every image to 58152x9960 first
    python batch_he_pyramid.py --src ./he_images --out ./he_out --resize 58152x9960

    # see what would happen
    python batch_he_pyramid.py --src ./he_images --out ./he_out --dry-run
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

from PIL import Image

Image.MAX_IMAGE_PIXELS = None

IMAGE_EXTS = (".tif", ".tiff", ".png", ".jpg", ".jpeg", ".webp", ".bmp")
UNSUPPORTED_EXTS = (".svs", ".ndpi", ".mrxs", ".czi", ".vms", ".vmu", ".scn", ".dcm")

# same suffixes the viewer strips when matching an image to a sample id
ID_SUFFIXES = (
    "_rescaled", "-rescaled", "_rescaled_image",
    "_he", "-he", "_HE", "_he_image",
)

# same filter the viewer uses to recognise a slide id; keeps logo.png and
# prediction PNGs out of the batch
SAMPLE_RE = re.compile(r"^(TCGA|TARGET|BWH|STANFORD|[A-Z]{2,6}[-_])", re.IGNORECASE)

# WebP refuses canvases whose longest side exceeds this
WEBP_MAX_SIDE = 16383


def sample_id_from(path: str) -> str:
    stem = os.path.splitext(os.path.basename(path))[0]
    for suf in ID_SUFFIXES:
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
            break
    return stem


def find_mask(he_path: str, sample_id: str) -> str | None:
    """A ``<sample>*mask*`` image next to the H&E, preferring HistoSweep."""
    d = os.path.dirname(he_path)
    cands = []
    for fn in sorted(os.listdir(d)):
        stem, ext = os.path.splitext(fn)
        if ext.lower() not in IMAGE_EXTS:
            continue
        if not stem.lower().startswith(sample_id.lower()):
            continue
        if "mask" not in stem.lower():
            continue
        cands.append(os.path.join(d, fn))
    if not cands:
        return None

    def rank(p):
        s = os.path.basename(p).lower()
        if "histosweep" in s:
            return 0
        return 1
    return sorted(cands, key=rank)[0]


def apply_mask(im: Image.Image, mask_path: str) -> Image.Image:
    """Paint everything outside the tissue mask white (reference-figure style)."""
    from PIL import ImageOps
    mask = Image.open(mask_path).convert("L")
    if mask.size != im.size:
        mask = mask.resize(im.size, Image.NEAREST)
    inv = ImageOps.invert(mask)                    # 255 = background, 0 = tissue
    white = Image.new("RGB", im.size, (255, 255, 255))
    return Image.composite(white, im, inv)


def discover(src: str, any_name: bool = False) -> list[str]:
    out = []
    for root, dirs, files in os.walk(src):
        # don't descend into already-built pyramids (their L*.webp would be
        # mistaken for whole-slide images named "L2" etc.)
        dirs[:] = sorted(d for d in dirs
                         if d.lower() not in {"images", "assets", "static", "icons",
                                              "templates"}
                         and not d.endswith("_pyramid"))
        for fn in sorted(files):
            p = os.path.join(root, fn)
            ext = os.path.splitext(fn)[1].lower()
            if "mask" in os.path.splitext(fn)[0].lower():
                continue                 # masks are applied, never converted
            if ext in UNSUPPORTED_EXTS:
                print(f"[skip] {fn}: whole-slide format, rescale it first", file=sys.stderr)
                continue
            if ext in IMAGE_EXTS:
                if any_name or SAMPLE_RE.match(sample_id_from(p)):
                    out.append(p)
                else:
                    print(f"[skip] {fn}: name does not look like a slide id", file=sys.stderr)
    return out


def build_pyramid(im: Image.Image, out_dir: str, n_levels: int, quality: int,
                  fmt: str) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    native_w, native_h = im.size
    manifest = {"size": [native_w, native_h], "levels": []}

    cur = im
    for level in range(n_levels):
        w, h = cur.size
        saved = None
        if fmt == "webp" and max(w, h) <= WEBP_MAX_SIDE:
            import io
            buf = io.BytesIO()
            cur.save(buf, format="WEBP", quality=quality, method=6)
            fn = f"L{level}.webp"
            with open(os.path.join(out_dir, fn), "wb") as fh:
                fh.write(buf.getvalue())
            saved = fn
        elif fmt == "png" and max(w, h) <= WEBP_MAX_SIDE:
            fn = f"L{level}.png"
            cur.save(os.path.join(out_dir, fn), format="PNG", compress_level=3)
            saved = fn
        manifest["levels"].append({"level": level, "size": [w, h],
                                   "file": saved, "skipped": saved is None})
        if min(w, h) <= 512:
            break
        cur = cur.resize((max(1, w // 2), max(1, h // 2)), Image.LANCZOS)

    with open(os.path.join(out_dir, "manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    return manifest


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", required=True, help="directory of H&E images")
    ap.add_argument("--out", required=True, help="where to write <id>_pyramid/ dirs")
    ap.add_argument("--resize", metavar="WxH",
                    help="resize every image to WxH before building (e.g. 58152x9960)")
    ap.add_argument("--levels", type=int, default=5)
    ap.add_argument("--quality", type=int, default=85)
    ap.add_argument("--format", default="webp", choices=("webp", "png"))
    ap.add_argument("--force", action="store_true", help="rebuild existing pyramids")
    ap.add_argument("--no-mask", action="store_true",
                    help="do not apply a matching tissue mask")
    ap.add_argument("--any", action="store_true",
                    help="convert images whose names do not look like slide ids")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    if a.resize:
        try:
            tw, th = (int(v) for v in a.resize.lower().split("x", 1))
        except ValueError:
            ap.error("--resize must look like 58152x9960")

    paths = discover(a.src, a.any)
    if not paths:
        print(f"no images found under {a.src}")
        return

    print(f"found {len(paths)} images")
    done = skipped = failed = 0
    for p in paths:
        sid = sample_id_from(p)
        out_dir = os.path.join(a.out, f"{sid}_pyramid")
        if os.path.exists(os.path.join(out_dir, "manifest.json")) and not a.force:
            print(f"[skip] {sid}: pyramid already built")
            skipped += 1
            continue
        if a.dry_run:
            print(f"[dry-run] {sid}: {p} -> {out_dir}")
            continue

        t0 = time.time()
        try:
            im = Image.open(p).convert("RGB")
            im.load()
            if a.resize:
                im = im.resize((tw, th), Image.LANCZOS)
            mask_note = ""
            if not a.no_mask:
                mask_path = find_mask(p, sid)
                if mask_path:
                    im = apply_mask(im, mask_path)
                    mask_note = f"  masked with {os.path.basename(mask_path)}"
            manifest = build_pyramid(im, out_dir, a.levels, a.quality, a.format)
            kept = [lv for lv in manifest["levels"] if lv["file"]]
            print(f"[ ok ] {sid}: {manifest['size'][0]}x{manifest['size'][1]} -> "
                  f"{len(kept)} levels ({', '.join(lv['file'] for lv in kept)})"
                  f"{mask_note}  in {time.time()-t0:.1f}s")
            done += 1
        except Exception as exc:
            print(f"[fail] {p}: {type(exc).__name__}: {exc}", file=sys.stderr)
            failed += 1

    print(f"\ndone: {done} built, {skipped} skipped, {failed} failed")


if __name__ == "__main__":
    main()
