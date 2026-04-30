#!/usr/bin/env python3

from __future__ import annotations

import os
import re
import subprocess
import sys
import tempfile
from fnmatch import fnmatch
from glob import glob
from pathlib import Path

import numpy as np
from astropy.io import fits as astrofits
from scipy.ndimage import gaussian_filter1d

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from spectral_line_catalog import wavelength_for


# ============================================================
# --- INPUTS ---
# ============================================================

# Input spectra: glob patterns matched inside spectra/ directory.
# Ignored if SPECTRA_LIST_FILE is set.
FITS_GLOBS = [
    "*CHIRON*.fits",
    # "*UVES*.fits",
    # "*BeSS*VIS_6*.fits",
]

# Optional: path to a hand-crafted .list file (one absolute path per line).
# When set, FITS_GLOBS is ignored and smoothing is skipped.
# Set to None to use FITS_GLOBS.
SPECTRA_LIST_FILE: str | None = "spectra.list"
# SPECTRA_LIST_FILE = "/home/rklement/Documents/dynspec/my_spectra.list"

# Tag used in output filenames (replaces the old list-file stem)
OUTPUT_TAG = "CHIRON_UVES"

# Per-instrument Gaussian pre-smoothing: {glob_pattern: fwhm_angstrom}
# Applied per-spectrum (before dynspec stacking) to files whose basename matches.
# {} to disable.
SMOOTH_FWHM = {
    "*CHIRON*.fits": 0.3,
}

# Output directory for plots
# OUTPUT_DIR = ROOT / "outputs"   # uncomment to also save FITS files
PLOT_DIR = ROOT / "plots"

# ---- WHAT TO RUN ----

# Spectral lines to process (names from spectral_line_catalog.py)
# (name, vmin, vmax[, title]) — vmin/vmax/title are None for auto/default
LINE_NAMES = [
    ("Hbeta", 0.17508, 1.45129, r"H$\beta$"),
    # HeI triplet lines
    # ("HeI 3889", None, None),
    # ("HeI 4471", None, None),
    # ("HeI 4713", 0.985, 1.01),
    # ("HeI 5876", 0.93827, 1.0331),
    # ("HeI 7065", 0.975, 1.01),
    # ("HeI 10830", None, None),
    # HeI singlet lines
    # ("HeI 3965", None, None),
    # ("HeI 4922", None, None),
    # ("HeI 5016", None, None),
    # ("HeI 6678", None, None),
    # ("HeI 7281", None, None),
    # HeII
    # ("HeII 4686", 0.994, 1.005),
    # Fe
    # ("FeII 4629", None, None),
    # ("FeII 4924", 0.75, 1.02),
    ("FeII 5018", 0.82, 1.06),
    ("FeIII 5127", 0.96872, 1.02314),
    ("FeII 5317", 0.89967, 1.06723),
    # ("NaI 5890", None, None),
    # ("NaI 5896", None, None),
]

# Modes to run (subset of: ts2ima, ts2vima, ts2phima, ts2vphima, ts2tab, ts2vtab, ts2cov)
MODES = [
    # "ts2vphima",
    "ts2vima",
    # "ts2ima",
    # "ts2phima",
    # "ts2tab",
    # "ts2vtab",
    # "ts2cov",
]

# ---- WAVELENGTH / VELOCITY WINDOW ----

WAVE_HALF_WIDTH_A = 15.0    # ±Angstrom window around line center
LAMSTEP_A         =  0.1    # wavelength step (for wavelength-axis modes)

VLO_KMS   = -400.0          # velocity lower limit (km/s)
VHI_KMS   =  400.0          # velocity upper limit (km/s)
VSTEP_KMS =    5.0          # velocity step (km/s)

# ---- TIME / PHASE ----

TSTEP_D = 4               # time step in days (for time-axis modes)

PHSTEP          = 0.0625    # phase step
OVERLAP_PERCENT = 0.0       # phase bin overlap in percent

ORBITAL_SOLUTIONS = [
    {"tag": "orb_period", "period": 87.742, "t0_mjd": 59988.22},
]

# ---- PROCESSING ----

RENORM        = False
MEDIAN_FILTER = False

# ---- BEHAVIOR ----

OVERWRITE     = True    # overwrite existing FITS/plots
STOP_ON_ERROR = False   # stop if a FITS generation job fails

# ---- PLOT SETTINGS ----

FIGURE_WIDTH_IN  = 2.0
FIGURE_HEIGHT_IN = 2.5
FONT_SIZE = 7.0
DPI       = 300
# CMAP      = "RdBu_r"
CMAP      = "nipy_spectral"  


PHASE_LO        = -0.5      # lower phase limit for phase plots
PHASE_HI        =  1.5      # upper phase limit for phase plots
USE_HEADER_CUTS = False     # use VCUTLO/VCUTHI from FITS header
CLIP_BELOW_ZERO = False     # hard-clip values below vmin
NAN_COLOR       = "black"
NAN_ALPHA       = 1.0
OUTPUT_FORMATS: list[str] = ["png", "pdf"]   # output formats: any combo of "png", "pdf"
PLOT_TITLE: str | None  = None   # plot title; None = no title
HIDE_XLABEL      = False   # hide x-axis label
HIDE_YLABEL      = True   # hide y-axis label
HIDE_CBAR_LABEL  = True   # hide colorbar label
INFO            = True     # print image stats before each plot


# ============================================================
# --- CODE ---
# ============================================================

PHASE_MODES = {"ts2phima", "ts2vphima"}
WAVE_MODES  = {"ts2ima", "ts2phima", "ts2tab"}
VELO_MODES  = {"ts2vima", "ts2vphima", "ts2vtab", "ts2cov"}


def slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "", text)


def tok(value: float, ndp: int) -> str:
    s = f"{value:.{ndp}f}".rstrip("0")
    if s.endswith("."):
        s += "0"
    return s.replace("-", "m").replace(".", "p")


def scale_token(value: float | None) -> str:
    if value is None:
        return "auto"
    return f"{value:.6g}".replace("-", "m").replace(".", "p")


def make_fits_name(mode, line_key, lamlo, lamhi, orb, output_tag):
    parts = [mode, slug(output_tag), slug(line_key)]
    if mode in WAVE_MODES:
        parts += [f"w{tok(lamlo, 2)}to{tok(lamhi, 2)}", f"dw{tok(LAMSTEP_A, 3)}"]
    if mode in VELO_MODES:
        parts += [f"v{tok(VLO_KMS, 1)}to{tok(VHI_KMS, 1)}", f"dv{tok(VSTEP_KMS, 3)}"]
    if mode in {"ts2ima", "ts2vima"}:
        parts.append(f"dt{tok(TSTEP_D, 3)}")
    if mode in PHASE_MODES and orb is not None:
        parts += [
            f"P{tok(float(orb['period']), 5)}",
            f"T0{tok(float(orb['t0_mjd']), 5)}",
            f"dph{tok(PHSTEP, 4)}",
        ]
    if mode == "ts2vphima":
        parts.append(f"ov{tok(OVERLAP_PERCENT, 3)}")
    parts.append("renorm" if RENORM else "norenorm")
    parts.append("medfilt" if MEDIAN_FILTER else "nomedfilt")
    return "_".join(parts) + ".fits"


def build_fits_cmd(mode, out_path, lamc, lamlo, lamhi, orb, list_path):
    cmd = [
        sys.executable, str(ROOT / f"{mode}.py"),
        "--out", str(out_path),
        "--list", str(list_path),
        "--extrapolate-x",
    ]
    if mode in WAVE_MODES:
        cmd += ["--lamlo", f"{lamlo}", "--lamhi", f"{lamhi}", "--lamstep", f"{LAMSTEP_A}"]
    if mode in VELO_MODES:
        cmd += ["--lamc", f"{lamc}", "--vlo", f"{VLO_KMS}", "--vhi", f"{VHI_KMS}", "--vstep", f"{VSTEP_KMS}"]
    if mode in {"ts2ima", "ts2vima"}:
        cmd += ["--tstep", f"{TSTEP_D}"]
    if mode in PHASE_MODES and orb is not None:
        cmd += ["--period", f"{orb['period']}", "--t0-mjd", f"{orb['t0_mjd']}", "--phstep", f"{PHSTEP}"]
    if mode == "ts2vphima":
        cmd += ["--overlap", f"{OVERLAP_PERCENT}"]
    cmd.append("--renorm" if RENORM else "--no-renorm")
    cmd.append("--filter" if MEDIAN_FILTER else "--no-filter")
    return cmd


def build_plot_cmd(in_fits, out_plot, vmin, vmax, title=None):
    cmd = [
        sys.executable, str(ROOT / "plot_dynspec.py"),
        str(in_fits), str(out_plot),
        "--figsize", str(FIGURE_WIDTH_IN), str(FIGURE_HEIGHT_IN),
        "--fontsize", str(FONT_SIZE),
        "--dpi", str(DPI),
        "--cmap", CMAP,
        "--nan-color", NAN_COLOR,
        "--nan-alpha", str(NAN_ALPHA),
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

    cmd += ["--phase-lo", str(PHASE_LO), "--phase-hi", str(PHASE_HI)]
    if title:
        cmd += ["--title", title]
    if HIDE_XLABEL:
        cmd.append("--hide-xlabel")
    if HIDE_YLABEL:
        cmd.append("--hide-ylabel")
    if HIDE_CBAR_LABEL:
        cmd.append("--hide-cbar-label")
    return cmd


def run_cmd(cmd):
    proc = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True)
    if proc.returncode != 0:
        print("FAILED:", " ".join(cmd))
        for line in (proc.stdout + proc.stderr).strip().splitlines()[-20:]:
            print(" ", line)
    return proc.returncode


def main():

    smooth_tmp_files = []
    owns_list = False

    if SPECTRA_LIST_FILE is not None:
        list_path = Path(SPECTRA_LIST_FILE)
        if not list_path.exists():
            raise FileNotFoundError(f"SPECTRA_LIST_FILE not found: {list_path}")
        n = sum(1 for ln in list_path.read_text().splitlines() if ln.strip())
        print(f"Using list file: {list_path} ({n} spectra)")
    else:
        _seen = set()
        spectra_files = []
        for pat in FITS_GLOBS:
            for f in sorted(glob(str(ROOT / "spectra" / pat))):
                if f not in _seen:
                    _seen.add(f)
                    spectra_files.append(f)
        if not spectra_files:
            raise FileNotFoundError(f"No FITS files found in spectra/ with globs {FITS_GLOBS}")
        print(f"Found {len(spectra_files)} spectra")

        final_files = []
        for f in spectra_files:
            fname = os.path.basename(f)
            fwhm = next((v for pat, v in SMOOTH_FWHM.items() if fnmatch(fname, pat)), 0)
            if fwhm > 0:
                with astrofits.open(f) as hdul:
                    cdelt1 = hdul[0].header["CDELT1"]
                    sigma_px = (fwhm / 2.3548) / cdelt1
                    smoothed = gaussian_filter1d(hdul[0].data.astype(np.float64), sigma=sigma_px)
                    hdul[0].data = smoothed.astype(hdul[0].data.dtype)
                    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".fits")
                    os.close(tmp_fd)
                    hdul.writeto(tmp_path, overwrite=True)
                smooth_tmp_files.append(tmp_path)
                final_files.append(tmp_path)
                print(f"  smoothed {fname} (FWHM={fwhm} Å, sigma={sigma_px:.2f} px)")
            else:
                final_files.append(f)

        tmp_list_fd, tmp_list_path = tempfile.mkstemp(suffix=".list")
        os.close(tmp_list_fd)
        list_path = Path(tmp_list_path)
        list_path.write_text("\n".join(final_files) + "\n")
        owns_list = True

    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    n_ok = n_skip = n_fail = 0
    p_ok = p_skip = p_fail = 0

    for entry in LINE_NAMES:
        line_key, vmin, vmax = entry[0], entry[1], entry[2]
        line_title = entry[3] if len(entry) > 3 else None
        lamc = float(wavelength_for(line_key))
        lamlo = lamc - WAVE_HALF_WIDTH_A
        lamhi = lamc + WAVE_HALF_WIDTH_A
        for mode in MODES:
            orb_list = ORBITAL_SOLUTIONS if mode in PHASE_MODES else [None]
            for orb in orb_list:
                plot_stem = Path(make_fits_name(mode, line_key, lamlo, lamhi, orb, OUTPUT_TAG)).stem
                cmap_slug = re.sub(r"[^A-Za-z0-9]", "", CMAP)
                base_name = f"{plot_stem}_vmin{scale_token(vmin)}_vmax{scale_token(vmax)}_{cmap_slug}"
                all_outputs = [
                    PLOT_DIR / f"{base_name}.{fmt}"
                    for fmt in OUTPUT_FORMATS
                ]
                if all(p.exists() for p in all_outputs) and not OVERWRITE:
                    print(f"SKIP {plot_stem}")
                    n_skip += 1
                    continue

                tmp_fd, tmp_path = tempfile.mkstemp(suffix=".fits")
                os.close(tmp_fd)
                tmp_fits = Path(tmp_path)
                print(f"RUN  [{mode}] {plot_stem}")
                rc = run_cmd(build_fits_cmd(mode, tmp_fits, lamc, lamlo, lamhi, orb, list_path))
                if rc != 0:
                    n_fail += 1
                    if STOP_ON_ERROR:
                        return 2
                    continue
                n_ok += 1

                for fmt in OUTPUT_FORMATS:
                    out_plot = PLOT_DIR / f"{base_name}.{fmt}"
                    if out_plot.exists() and not OVERWRITE:
                        print(f"SKIP plot {out_plot.name}")
                        p_skip += 1
                        continue
                    print(f"PLOT -> {out_plot.name}")
                    rc = run_cmd(build_plot_cmd(tmp_fits, out_plot, vmin, vmax, title=PLOT_TITLE or line_title or line_key))
                    if rc == 0:
                        p_ok += 1
                    else:
                        p_fail += 1

                tmp_fits.unlink(missing_ok=True)

    if owns_list:
        list_path.unlink(missing_ok=True)
    for f in smooth_tmp_files:
        Path(f).unlink(missing_ok=True)
    print(f"\nDone.  FITS ok={n_ok} skip={n_skip} fail={n_fail}  |  plots ok={p_ok} skip={p_skip} fail={p_fail}")
    return 0 if (n_fail + p_fail) == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
