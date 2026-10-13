"""
Verify the README's H&E instructions are literally true.

Each check mirrors a documented instruction, so a passing run means a reader can
follow the README and succeed.
"""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from he_images import HEImageStore  # noqa: E402

SID = "TCGA-JL-A3YX-01A-01-TS1"
failures = []


def check(label, cond, extra=""):
    ok = bool(cond)
    print(("PASS  " if ok else "FAIL  ") + label + (f"  {extra}" if extra else ""))
    if not ok:
        failures.append(label)


print("=" * 84)
print("README claim: AURORA_HE_DIR is scanned recursively (subfolders are fine)")
print("=" * 84)
tmp = os.path.join(DATA, "_readme_probe")
shutil.rmtree(tmp, ignore_errors=True)
deep = os.path.join(tmp, "another_batch")
os.makedirs(deep)
img = os.path.join(deep, f"{SID}_rescaled.tif")
shutil.copyfile(
    os.path.join(DATA, f"{SID}_rescaled.tif") if os.path.exists(
        os.path.join(DATA, f"{SID}_rescaled.tif")) else
    os.path.join(DATA, "TCGA-05-4244-01A-01-BS1.png"),
    img)
st = HEImageStore(directory=tmp)
check("a file in a subfolder is found", st.path_for(SID) == img, st.path_for(SID))
check("describe() reports it as available and not a placeholder",
      st.describe(SID)["available"] and not st.describe(SID)["placeholder"])
shutil.rmtree(tmp, ignore_errors=True)

print()
print("=" * 84)
print("README claim: these filenames resolve to the sample id")
print("=" * 84)
GOOD = [f"{SID}_rescaled.tif", f"{SID}-rescaled.tif", f"{SID}_rescaled_image.tif",
        f"{SID}.png", f"{SID}.jpg", f"{SID}.tif", f"{SID}.webp",
        f"{SID}_he.png", f"{SID}-he.jpg", f"{SID}_HE.tif", f"{SID}_he_image.png"]
BAD = [f"{SID}.svs", f"{SID}.ndpi", "logo.png", "zoomout_rest.png"]
tmp = os.path.join(DATA, "_readme_probe")
for fn in GOOD + BAD:
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(tmp)
    open(os.path.join(tmp, fn), "wb").close()
    got = HEImageStore(directory=tmp).path_for(SID)
    if fn in GOOD:
        check(f"{fn} resolves", got is not None)
    else:
        check(f"{fn} is ignored", got is None)
shutil.rmtree(tmp, ignore_errors=True)

print()
print("=" * 84)
print("README claim: the workspace sample works with AURORA_HE_DIR=repo root")
print("=" * 84)
st = HEImageStore(directory=DATA)
check("the rescaled TIFF in the repo root is found",
      st.path_for(SID) is not None, st.describe(SID).get("basename"))
check("viewer chrome under images/ is not indexed as a sample",
      st.path_for("zoomout_rest") is None)
idx = st.index()
check("the H&E index holds few entries (only specimen ids)",
      len(idx) <= 10, f"{len(idx)}: {sorted(idx)}")

print()
print("=" * 84)
print("README claim: pyramid layouts the viewer looks for")
print("=" * 84)
for layout in ("{sid}", "he/{sid}", "{sid}_pyramid"):
    tmp = os.path.join(DATA, "_readme_probe")
    shutil.rmtree(tmp, ignore_errors=True)
    d = os.path.join(tmp, layout.format(sid=SID))
    os.makedirs(d)
    # a 1x1 webp stands in for a level image
    from PIL import Image
    Image.new("RGB", (100, 20)).save(os.path.join(d, "L2.webp"), format="WEBP")
    st = HEImageStore(directory=tmp)
    found = st.pyramid_dir(SID)
    # normalise: the layout string uses "/" but the platform separator may differ
    check(f"layout {layout} is discovered",
          found is not None and os.path.normpath(found) == os.path.normpath(d), found)
    shutil.rmtree(tmp, ignore_errors=True)

print()
print("=" * 84)
print("README claim: levels are picked by width")
print("=" * 84)
from PIL import Image  # noqa: E402
tmp = os.path.join(DATA, "_readme_probe")
shutil.rmtree(tmp, ignore_errors=True)
d = os.path.join(tmp, SID)
os.makedirs(d)
for name, w in (("L2.webp", 14538), ("L3.webp", 7269), ("L4.webp", 3634)):
    Image.new("RGB", (w, max(1, w // 6))).save(os.path.join(d, name), format="WEBP")
st = HEImageStore(directory=tmp)
# ask for 2500 -> L4 (3634) is the smallest that covers it
resp = st.serve_pyramid(d, 2500)
check("a 2500 px request is served from the pyramid", resp is not None)
if resp is not None:
    check("it picks L4 (3634 >= 2500)", resp.headers.get("X-Level-Width") == "3634",
          resp.headers.get("X-Level-Width"))
resp = st.serve_pyramid(d, 8000)
check("an 8000 px request picks L2 (14538)", resp is not None and
      resp.headers.get("X-Level-Width") == "14538",
      None if resp is None else resp.headers.get("X-Level-Width"))
shutil.rmtree(tmp, ignore_errors=True)

print()
print("=" * 84)
print("README claim: no configuration is needed when the files sit beside space/")
print("=" * 84)
for name in ("AURORA_LOCAL_DIR", "AURORA_HE_DIR", "AURORA_BACKEND"):
    os.environ.pop(name, None)
import importlib  # noqa: E402
import app as _app  # noqa: E402
importlib.reload(_app)
check("AURORA_LOCAL_DIR is auto-detected to the repository root",
      os.environ.get("AURORA_LOCAL_DIR") == DATA,
      os.environ.get("AURORA_LOCAL_DIR"))
check("AURORA_HE_DIR is auto-detected to the repository root",
      os.environ.get("AURORA_HE_DIR") == DATA,
      os.environ.get("AURORA_HE_DIR"))
check("the auto-detected directory yields the sample",
      len(_app.store.samples()) >= 1, len(_app.store.samples()))
_h = _app.health()
check("health reports the directory in use and no listing error",
      _h["sample_dir"] == DATA and _h["listing_error"] is None
      and _h["samples_loaded"] >= 1,
      (_h["sample_dir"], _h["listing_error"], _h["samples_loaded"]))

# an explicit setting must win over auto-detection
os.environ["AURORA_LOCAL_DIR"] = os.path.join(DATA, "space")
importlib.reload(_app)
check("an explicit AURORA_LOCAL_DIR is not overridden",
      os.environ["AURORA_LOCAL_DIR"] == os.path.join(DATA, "space"),
      os.environ["AURORA_LOCAL_DIR"])
os.environ.pop("AURORA_LOCAL_DIR", None)

print()
if failures:
    print(f"{len(failures)} FAILURE(S): " + "; ".join(failures))
    sys.exit(1)
print("ALL README H&E CLAIMS VERIFIED")
