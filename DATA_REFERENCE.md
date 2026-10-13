# Data reference: what each AURORA product contains

Verified against real files, not just the producer code.

## `predictions_7um_28um` - fine + coarse (schema v2)

One `*-highres-pred.h5` per sample (115 MB - 1.8 GB; 885 files for TCGA_LUAD).

| dataset | shape | notes |
| --- | --- | --- |
| `exp_pred_14` | (genes, patches) uint8 | 7 um |
| `exp_pred_56` | (genes, patches) uint8 | 28 um |
| `prop_pred_14` / `_56` | (celltypes, patches) uint8 | cell-type proportions, **log** scale |
| `exp_max_14` / `_56` | (genes,) float32 | per-gene quantisation maximum |
| `prop_max_14` / `_56` | (celltypes,) float32 | per-cell-type maximum |
| `patch_loc_14` / `_56` | (N, 2) | **(y, x)** top-left corners |
| `genes`, `cell_types` | object | names |

* Chunked `(256 channels, 16384 spatial rows)`, zstd (HDF5 filter id 32015), with
  `shuffle -> zstd` for the auxiliary 1-D/2-D arrays. **h5py cannot read these
  without a plugin** - see `space/aurora_h5.py` for the raw-chunk workaround.
* Encoding `spatial_gene_delta_v1`: uint8 modulo-256 first difference along the
  spatially sorted rows, restarted every `delta_block_rows` (16384) rows.
* Grid stride 14 px = 7 um and 56 px = 28 um, in the 0.5 um/px rescaled image.
* Gene counts and cell-type names **differ between samples**: files seen so far
  give 3002 or 3000 genes, and 10 or 9 cell types with different labels
  (`B_cell` vs `B-cells`, `Malignant` vs `Cancer Epithelial`). Read them from the
  file; never hard-code.
* **No H&E images are in this repo.**

## `predictions_112um` - coarse, different schema (v1)

One `*.h5` per sample inside `TCGA_LUAD_predictions_V2.tar.gz` (9.4 GB, gzip, so
whole-file download only). Per-sample files are small (4-20 MB) and are laid out
as plain contiguous arrays - **no delta encoding, no compression, no chunking**:

| dataset | shape | dtype |
| --- | --- | --- |
| `exp_pred_224` | (patches, genes) | uint8 |
| `prop_pred_224` | (patches, celltypes) | uint8 |
| `patch_loc_224` | (patches, 2) | int16 |
| `exp_max_224` | (genes,) | float32 |
| `prop_max_224` | (celltypes,) | float32 |
| `resolutions` | (1,) | int64 = [224] |
| `genes` / `cell_types` | object | names |

Key differences, each of which the reader handles:

1. **Axis order is transposed.** `exp_pred_224` is (patches, genes); the fine
   files store (genes, patches).
2. **No delta encoding** - the values are absolute uint8 codes already.
3. **No `row_order`** dataset: rows are not spatially sorted into a raster; the
   order is `patch_loc_224` order.
4. **224 px patches = 112 um** at 0.5 um/px, i.e. 16 x 16 fine (7 um) cells.

## Coordinate convention

`patch_loc_*` stores **(y, x)**, not (x, y): column 0 is vertical, column 1 is
horizontal. Settled empirically two ways - column 1 is monotonic through the
stored file (the writer sorts by the primary axis), and scoring the prediction
footprint against tissue in a real rescaled H&E gives IoU 0.79 for (y, x) versus
0.19 or less for (x, y).

Coordinates are used **exactly as stored**. A patch may start at a negative
coordinate, meaning part of it lies outside the image; the raster origin is then
the minimum stored coordinate and only the visible portion is represented. See
`HE_FORMAT.md` for the details and for how the H&E itself should be stored.

## What is *not* on the Hub

* **H&E images.** The old viewer's README pointed at
  `AURORAData/prediction_plots`, which now returns 404. `AURORA_prediction`
  contains zero image files; `AURORAData` has no other repos. The rescaled H&E has
  to be regenerated or supplied, as the user did for `TCGA-JL-A3YX-01A-01-TS1`.
* **A per-sample 112 um file.** Only the 9.4 GB gzip tar exists, and gzip has no
  random access, so a single member cannot be fetched by range request. Extract
  it, or upload per-sample files.
