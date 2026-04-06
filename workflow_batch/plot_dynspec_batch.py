#!/usr/bin/env python3
"""
Batch plotting helper for workflow_batch outputs.

All plotting inputs are configured in this file (no command-line arguments required).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = PROJECT_ROOT / "workflow_batch" / "outputs"
PLOT_DIR = PROJECT_ROOT / "workflow_batch" / "plots"


# ---------------------------------------------------------------------------
# User configuration (edit this block)
# ---------------------------------------------------------------------------

# Exact FITS filenames to plot (relative to INPUT_DIR).
# Example:
# INPUT_FILES = [
#     "ts2vphima__spectraUVES__FeII5469__... .fits",
#     "ts2vphima__spectraUVES__FeII6318__... .fits",
# ]
# If this is not empty, only these files are plotted.
INPUT_FILES: list[str] = []

# Optional glob-based selection (used only when INPUT_FILES is empty).
# Example:
# INPUT_GLOBS = ["*FeII5469*.fits", "*FeII6318*.fits"]
INPUT_GLOBS: list[str] = ['*5018*.fits']

# Only image-like dynspec outputs are plotted by default.
ALLOWED_PREFIXES = ("ts2ima__", "ts2vima__", "ts2phima__", "ts2vphima__", "ts2cov__")

PLOT_EXTENSION = ".png"
OVERWRITE = True

# plot_dynspec.py options
PAGE = "letter"
ORIENTATION = "landscape"
DPI = 180
CMAP = "viridis"
USE_HEADER_CUTS = True
CLIP_BELOW_ZERO = False
INFO = False
NAN_COLOR = "black"
NAN_ALPHA = 1.0
# Color-scale pairs. These lists must have the same length.
# Each (VMIN_LIST[i], VMAX_LIST[i]) pair creates one output plot per FITS file.
# Use None for auto-limits.
VMIN_LIST: list[float | None] = [0.78007]#[0.9695]
VMAX_LIST: list[float | None] = [None]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def collect_input_fits() -> list[Path]:
    if INPUT_FILES:
        return [INPUT_DIR / name for name in INPUT_FILES]

    if INPUT_GLOBS:
        found: dict[Path, None] = {}
        for pat in INPUT_GLOBS:
            for p in sorted(INPUT_DIR.glob(pat)):
                found[p] = None
        return list(found.keys())

    fits_files = sorted(INPUT_DIR.glob("*.fits"))
    return [p for p in fits_files if p.name.startswith(ALLOWED_PREFIXES)]


def build_plot_command(in_fits: Path, out_plot: Path, vmin: float | None, vmax: float | None) -> list[str]:
    cmd = [
        sys.executable,
        str(PROJECT_ROOT / "plot_dynspec.py"),
        str(in_fits),
        str(out_plot),
        "--page",
        PAGE,
        "--orientation",
        ORIENTATION,
        "--dpi",
        str(DPI),
        "--cmap",
        CMAP,
        "--nan-color",
        NAN_COLOR,
        "--nan-alpha",
        str(NAN_ALPHA),
    ]
    if USE_HEADER_CUTS:
        cmd.append("--use-header-cuts")
    if CLIP_BELOW_ZERO:
        cmd.append("--clip-below-zero")
    if INFO:
        cmd.append("--info")
    if vmin is not None:
        cmd += ["--vmin", str(vmin)]
    if vmax is not None:
        cmd += ["--vmax", str(vmax)]
    return cmd


def run_one(cmd: list[str]) -> int:
    proc = subprocess.run(cmd, cwd=PROJECT_ROOT, text=True, capture_output=True)
    if proc.returncode != 0:
        print("FAILED:", " ".join(cmd))
        if proc.stdout.strip():
            print(proc.stdout.strip().splitlines()[-20:])
        if proc.stderr.strip():
            print(proc.stderr.strip().splitlines()[-20:])
    return int(proc.returncode)


def scale_token(value: float | None) -> str:
    if value is None:
        return "auto"
    s = f"{value:.6g}"
    return s.replace("-", "m").replace(".", "p")


def main() -> int:
    if not INPUT_DIR.exists():
        raise FileNotFoundError(f"Input directory does not exist: {INPUT_DIR}")

    if len(VMIN_LIST) != len(VMAX_LIST):
        raise ValueError("VMIN_LIST and VMAX_LIST must have the same length.")
    scale_pairs = list(zip(VMIN_LIST, VMAX_LIST))

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    in_files = collect_input_fits()
    if not in_files:
        print("No FITS files found to plot.")
        return 0

    n_ok = 0
    n_skip = 0
    n_fail = 0

    print(f"Input directory: {INPUT_DIR}")
    print(f"Plot directory : {PLOT_DIR}")
    print(f"Total files    : {len(in_files)}")
    print(f"Scale pairs    : {len(scale_pairs)}")
    for i, (vmin, vmax) in enumerate(scale_pairs, start=1):
        print(f"  pair {i}: vmin/vmax = {vmin} / {vmax}")

    for in_fits in in_files:
        if not in_fits.exists():
            print(f"MISSING: {in_fits}")
            n_fail += 1
            continue

        for vmin, vmax in scale_pairs:
            suffix = f"__vmin{scale_token(vmin)}__vmax{scale_token(vmax)}"
            out_plot = PLOT_DIR / f"{in_fits.stem}{suffix}{PLOT_EXTENSION}"
            if out_plot.exists() and not OVERWRITE:
                print(f"SKIP {out_plot.name}")
                n_skip += 1
                continue

            print(f"PLOT {in_fits.name} -> {out_plot.name}")
            cmd = build_plot_command(in_fits, out_plot, vmin, vmax)
            rc = run_one(cmd)
            if rc == 0:
                n_ok += 1
            else:
                n_fail += 1

    print("")
    print(f"Done. ok={n_ok} skip={n_skip} fail={n_fail}")
    return 0 if n_fail == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
