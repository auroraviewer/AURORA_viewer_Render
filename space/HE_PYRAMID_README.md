# H&E pyramid converter

Converts H&E slide images into the small, fast "pyramid" format used by the
AURORA viewer: a `<sample_id>_pyramid/` folder of halved WebP levels plus a
`manifest.json` with the native size.

## Requirements

- Python 3.8+
- Pillow

```bash
pip install Pillow
```

## Run

```bash
python batch_he_pyramid.py --src <folder-of-images> --out <output-folder>
```

It converts **only** the `<sample>_rescaled.tif` images. If a matching tissue mask
(`<sample>_mask_rescaled.tif` or `<sample>_HistoSweep_mask_rescaled.tif`) sits next
to an image, it is **applied automatically** — everything outside the tissue is
painted white, like the reference figures. Skip that with `--no-mask`.

Options:

| flag | meaning |
| --- | --- |
| `--resize WxH` | resize each image first, e.g. `--resize 58152x9960` |
| `--no-mask` | do not apply a matching tissue mask |
| `--levels N` | number of levels (default 5) |
| `--quality N` | WebP quality 0–100 (default 85) |
| `--format png` | write PNG levels instead of WebP |
| `--force` | rebuild pyramids that already exist |
| `--dry-run` | list what would happen without writing |

## Notes

- Only `<sample>_rescaled.tif` images are converted; `*mask*` files are consumed
  as masks, not converted as separate images.
- The image must be at the **0.5 µm/px rescaled size** the predictions use, or
  pass `--resize WxH` to set that size.
- Whole-slide formats (`.svs`, `.ndpi`, `.mrxs`) are skipped — rescale them first.
- Output looks like:

```
<sample_id>_pyramid/
    L2.webp        # native / 4  (levels >16383 px are skipped: WebP's limit)
    L3.webp        # native / 8
    L4.webp        # native / 16
    manifest.json  # {"size": [native_w, native_h], ...}
```
