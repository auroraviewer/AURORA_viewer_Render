"""
Regression test for the partial-chunk stride bug.

exp_pred_* is chunked (256 channels, 16384 rows).  A partially written final
chunk still occupies its full physical extent on disk, and the decompressed
buffer keeps the PHYSICAL row stride.  Reshaping by the logical extent therefore
shifted every channel after the first into the padding - the final block of
every level decoded as garbage after gene 0, which looked exactly like a
corrupt prediction block.

This test compares the reader's output for the final partial chunk against the
native h5py read (needs hdf5plugin, which registers the zstd filter); when the
plugin is absent the test is skipped rather than failed, because the comparison
is the only independent ground truth available.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from aurora_h5 import AuroraPredictions, read_array, delta_decode_u8  # noqa: E402

failures = []


def check(label, cond, extra=""):
    ok = bool(cond)
    print(("PASS  " if ok else "FAIL  ") + label + (f"  {extra}" if extra else ""))
    if not ok:
        failures.append(label)


S = "TCGA-JL-A3YX-01A-01-TS1"
FINE = os.path.join(DATA, f"{S}-highres-pred.h5")

try:
    import hdf5plugin  # noqa: F401
    import h5py
    HAVE_PLUGIN = True
except ImportError:
    HAVE_PLUGIN = False

print("=" * 78)
print("partial final chunk: reader vs native h5py")
print("=" * 78)

f = AuroraPredictions([FINE, os.path.join(DATA, f"{S}-pred.h5")])
sel = np.asarray(list(range(0, len(f.genes), 29)))

for cell in (14, 56):
    L = f.layer(cell)
    blk = L._block
    row0 = (L.n_patches // blk) * blk
    nblk = L.n_patches - row0

    mine = L.exp_codes(sel.tolist(), np.arange(row0, L.n_patches))  # (nblk, k)

    if HAVE_PLUGIN:
        with h5py.File(FINE, "r") as h:
            d = h[f"exp_pred_{cell}"]
            stored = np.ascontiguousarray(d[:, row0:L.n_patches])     # (3000, nblk)
            gt = np.cumsum(stored[sel], axis=1, dtype=np.uint8).T     # (nblk, k)
        same = np.array_equal(mine, gt)
        check(f"{cell}px final block matches the native h5py decode", same)
        if not same:
            diff = np.nonzero(mine.astype(np.int16) - gt.astype(np.int16))
            print(f"      {diff[0].size} mismatching cells; first at "
                  f"row {diff[0][0]}, gene-sel {diff[1][0]}")
    else:
        print(f"SKIP  {cell}px comparison (hdf5plugin not installed)")

    # the second gene must not decode to zeros just because the chunk is partial:
    # channel k >= 1 sat in the padding under the old reshape
    col1 = mine[:, 1]
    nz = int((col1 != 0).sum())
    check(f"{cell}px second sampled gene has real values in the final block",
          nz > 0.5 * nblk, f"{nz} nonzero of {nblk}")

    # and re-differencing must reproduce the stored bytes, whatever they are
    red = np.empty_like(mine)
    red[0] = mine[0]
    red[1:] = (mine[1:].astype(np.int16) - mine[:-1].astype(np.int16)) % 256
    back = red.astype(np.uint8)
    enc = read_array(f.layer(cell).f[f"exp_pred_{cell}"], sel,
                     np.arange(row0, L.n_patches))
    check(f"{cell}px final block re-differences to the stored bytes",
          np.array_equal(back.T, enc))

f.close()

print()
if failures:
    print(f"{len(failures)} FAILURE(S): " + "; ".join(failures))
    sys.exit(1)
print("PARTIAL CHUNK CHECKS PASSED")
