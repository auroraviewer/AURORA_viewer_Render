"""
Live end-to-end check against a real Hugging Face dataset.

Skipped unless AURORA_LIVE_HF=1, because it needs the network and a specific
repo.  Run it manually before deploying a Space:

    $env:AURORA_LIVE_HF = '1'
    python test_hf_live.py

Covers the exact layout the Space expects:
  <repo>/predictions_7um_28um/<id>-highres-pred.h5
  <repo>/predictions_112um/<id>-pred.h5
  <repo>/HE_image/<id>_pyramid/{L2,L3,L4}.webp
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

REPO = os.environ.get("AURORA_LIVE_HF_REPO", "AURORAData/AURORA_Viewer_test")
S = os.environ.get("AURORA_LIVE_HF_SAMPLE", "TCGA-JL-A3YX-01A-01-TS1")

if os.environ.get("AURORA_LIVE_HF") != "1":
    print("SKIP  set AURORA_LIVE_HF=1 to run the live Hugging Face check")
    sys.exit(0)

failures = []


def check(label, cond, extra=""):
    ok = bool(cond)
    print(("PASS  " if ok else "FAIL  ") + label + (f"  {extra}" if extra else ""))
    if not ok:
        failures.append(label)


from samples import SampleStore
from he_images import HEImageStore

print("=" * 78)
print("PREDICTIONS (hf backend)")
print("=" * 78)
store = SampleStore(backend="hf", repo_id=REPO, subfolder="")
store.refresh()
ids = [s.sample_id for s in store.samples()]
print(f"  samples: {ids}")
check("finds the sample", S in ids)
s = store.get(S)
check("merges 7/28 and 112 um files", len(s.paths) == 2, s.paths)
handles = store.open(S)
check("opens the sample", len(handles) >= 1)
for h in handles:
    h.close()

print()
print("=" * 78)
print("H&E PYRAMID (hf repo)")
print("=" * 78)
he = HEImageStore(repo_id=REPO, subfolder="HE_image")
pyr = he.hf_pyramid(S)
check("discovers the pyramid", pyr is not None)
if pyr:
    print(f"  levels: {[(w, h) for w, h, _ in pyr['levels']]}")
    print(f"  native: {pyr['native']}")
    check("infers the native size", pyr["native"] == [58152, 9960], pyr["native"])
    check("describe reports available", he.describe(S)["available"])
    for mw in (1200, 3000, 8000):
        resp = he.serve_hf_pyramid(S, mw)
        check(f"serves max_width={mw}", resp is not None and resp.status_code == 200
              and len(resp.body) > 0,
              resp.headers.get("X-Level-Width") if resp else None)

print()
if failures:
    print(f"{len(failures)} FAILURE(S): " + "; ".join(failures))
    sys.exit(1)
print("LIVE HF CHECK PASSED")
