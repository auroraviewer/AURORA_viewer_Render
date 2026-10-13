"""
Download prediction H5 files from a Hugging Face dataset for local use.

Useful when you would rather mirror the corpus on disk than stream it, e.g.:

    python download_samples.py --repo AURORAData/AURORA_prediction \\
        --subfolder predictions_7um_28um/TCGA_LUAD --out D:\\aurora_h5

    python download_samples.py --repo ... --subfolder ... --out ... --limit 5
    python download_samples.py --repo ... --subfolder ... --out ... --sample TCGA-44-2668-01A-01-TS1

Then run the viewer against the folder:

    $env:AURORA_BACKEND='local'; $env:AURORA_LOCAL_DIR='D:\\aurora_h5'

Note each sample is roughly 600 MB, so budget disk accordingly.  A file that is
already present locally is skipped unless --overwrite is given.
"""

from __future__ import annotations

import argparse
import os
import sys


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--repo", required=True, help="dataset repo id, e.g. ORG/NAME")
    p.add_argument("--subfolder", default="", help="folder inside the repo")
    p.add_argument("--out", required=True, help="local directory to write into")
    p.add_argument("--revision", default=None, help="branch, tag or commit")
    p.add_argument("--limit", type=int, default=None, help="stop after N files")
    p.add_argument("--sample", action="append", default=None,
                   help="only this sample id (repeatable)")
    p.add_argument("--overwrite", action="store_true", help="re-download existing files")
    p.add_argument("--token", default=None, help="HF token (else HF_TOKEN env)")
    p.add_argument("--list-only", action="store_true", help="show what would be fetched")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)

    from huggingface_hub import HfApi, hf_hub_download

    token = args.token or os.environ.get("HF_TOKEN") or os.environ.get("HUGGING_FACE_HUB_TOKEN")
    api = HfApi(token=token)

    sub = args.subfolder.strip("/")
    print(f"listing {args.repo}" + (f"/{sub}" if sub else "") + " ...")
    entries = api.list_repo_tree(
        args.repo, path_in_repo=sub or None, recursive=True,
        repo_type="dataset", revision=args.revision,
    )
    paths = [e.path for e in entries
             if getattr(e, "path", "").endswith(".h5") and "-highres-pred.h5" in e.path]

    if args.sample:
        wanted = set(args.sample)
        paths = [p for p in paths
                 if os.path.basename(p)[: -len("-highres-pred.h5")] in wanted]
    paths.sort()
    if args.limit:
        paths = paths[: args.limit]

    total = len(paths)
    print(f"{total} prediction file(s) selected")
    if not total:
        print("nothing to do")
        return 1
    if args.list_only:
        for p in paths:
            print("  ", p)
        return 0

    os.makedirs(args.out, exist_ok=True)
    failures = 0
    for i, path in enumerate(paths, 1):
        name = os.path.basename(path)
        target = os.path.join(args.out, name)
        if os.path.exists(target) and not args.overwrite:
            print(f"[{i}/{total}] skip (exists) {name}")
            continue
        print(f"[{i}/{total}] downloading {name} ...", flush=True)
        try:
            local = hf_hub_download(
                repo_id=args.repo,
                filename=path,
                repo_type="dataset",
                revision=args.revision,
                token=token,
            )
            # hf_hub_download returns a cache path; mirror it under --out so the
            # folder can be used directly as AURORA_LOCAL_DIR.
            import shutil

            if os.path.abspath(local) != os.path.abspath(target):
                shutil.copyfile(local, target)
            size = os.path.getsize(target)
            print(f"          -> {target}  ({size/1e6:.1f} MB)")
        except Exception as exc:
            failures += 1
            print(f"          FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)

    print()
    print(f"done: {total - failures}/{total} ok"
          + (f", {failures} failed" if failures else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
