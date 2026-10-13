"""
Build and measure a real H&E pyramid, to decide how to store it on the Hub.

Produces the actual files (WebP per level) and reports exact totals, so the
answer is a measurement rather than an estimate.

Usage:  python build_he_pyramid.py --tif <rescaled.tif> --out <dir> [--levels 5]
"""
import argparse
import io
import os
import sys
import time

from PIL import Image

Image.MAX_IMAGE_PIXELS = None

# WebP refuses very large canvases; past this the levels are split into tiles.
WEBP_MAX_SIDE = 16383


def build(tif_path, out_dir, n_levels=5, quality=85, fmt="webp", tile=0):
    os.makedirs(out_dir, exist_ok=True)
    t0 = time.time()
    im = Image.open(tif_path).convert("RGB")
    im.load()
    print(f"decoded {im.size} in {time.time()-t0:.1f}s")

    results = []
    cur = im
    for level in range(n_levels):
        w, h = cur.size
        t0 = time.time()
        if fmt == "webp" and max(w, h) <= WEBP_MAX_SIDE:
            buf = io.BytesIO()
            cur.save(buf, format="WEBP", quality=quality, method=6)
            data = buf.getvalue()
            name = f"L{level}.webp"
            with open(os.path.join(out_dir, name), "wb") as fh:
                fh.write(data)
            results.append((level, (w, h), len(data), name))
        else:
            # too large for the codec: split into tiles at this level
            nt = 0
            total = 0
            step = tile or 4096
            for ty in range(0, h, step):
                for tx in range(0, w, step):
                    box = (tx, ty, min(tx + step, w), min(ty + step, h))
                    b = io.BytesIO()
                    cur.crop(box).save(b, format="PNG", compress_level=3)
                    fn = f"L{level}_{tx}_{ty}.png"
                    with open(os.path.join(out_dir, fn), "wb") as fh:
                        fh.write(b.getvalue())
                    total += b.tell()
                    nt += 1
            results.append((level, (w, h), total, f"{nt} PNG tiles"))
        print(f"  L{level}: {w}x{h}  {results[-1][2]/1e6:7.2f} MB  "
              f"({time.time()-t0:.1f}s)  {results[-1][3]}")
        if min(w, h) <= 512:
            break
        cur = cur.resize((max(1, w // 2), max(1, h // 2)), Image.LANCZOS)

    total = sum(r[2] for r in results)
    print()
    print(f"TOTAL: {total/1e6:.2f} MB across {len(results)} levels")
    return results, total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tif", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--levels", type=int, default=5)
    ap.add_argument("--quality", type=int, default=85)
    ap.add_argument("--format", default="webp", choices=("webp", "png"))
    ap.add_argument("--tile", type=int, default=0,
                    help="tile size for levels too large for the codec")
    a = ap.parse_args()
    build(a.tif, a.out, a.levels, a.quality, a.format, a.tile)


if __name__ == "__main__":
    main()
