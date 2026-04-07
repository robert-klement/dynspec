#!/usr/bin/env python3
"""
Central batch runner for the dynspec suite.

All inputs are configured in this file (no command-line arguments required).
This script launches the existing working tools in the parent directory and
writes outputs with deterministic names so files do not overwrite each other.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from spectral_line_catalog import wavelength_for  # noqa: E402


# ---------------------------------------------------------------------------
# User configuration (edit this block)
# ---------------------------------------------------------------------------

LIST_FILE = PROJECT_ROOT / "spectra_CHIRON.list"
OUTPUT_DIR = PROJECT_ROOT / "workflow_batch" / "outputs"

# Any subset of:
#   ts2ima, ts2vima, ts2phima, ts2vphima, ts2tab, ts2vtab, ts2cov
MODES = [
    # "ts2vima",
    "ts2vphima",
    # "ts2ima",
    # "ts2phima",
    # "ts2tab",
    # "ts2vtab",
    # "ts2cov",
]

# Use names or aliases from spectral_line_catalog.py
LINE_NAMES = [
    "HeI 6678",
    # "FeII 5018",
    # "FeII 5469",
]

# Wavelength-region controls
WAVE_HALF_WIDTH_A = 15.0
LAMSTEP_A = 0.1  # coarse by default for fast workflow

# Velocity-region controls
VLO_KMS = -600.0
VHI_KMS = 600.0
VSTEP_KMS = 5  # coarse by default for fast workflow

# Time/phase controls
TSTEP_D = 5
PHSTEP =  0.0625 # 1/16 in phase, coarse by default for fast workflow
OVERLAP_PERCENT = 0.0

# Multiple orbital solutions for phase modes
ORBITAL_SOLUTIONS = [
    {"tag": "orb_period", "period": 87.49, "t0_mjd": 59988.29},
    # {"tag": "VR_cycle", "period": 1800, "t0_mjd": 55000.00},
]

# Optional processing flags
RENORM = False
MEDIAN_FILTER = False

# Behavior
OVERWRITE = True
STOP_ON_ERROR = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

PHASE_MODES = {"ts2phima", "ts2vphima"}
WAVE_MODES = {"ts2ima", "ts2phima", "ts2tab"}
VELO_MODES = {"ts2vima", "ts2vphima", "ts2vtab", "ts2cov"}


def slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "", text)


def tok(value: float, ndp: int) -> str:
    s = f"{value:.{ndp}f}".rstrip("0").rstrip(".")
    return s.replace("-", "m").replace(".", "p")


def make_output_name(
    mode: str,
    line_key: str,
    lamc: float,
    lamlo: float,
    lamhi: float,
    orb: dict[str, float] | None,
) -> str:
    parts = [mode, slug(LIST_FILE.stem), slug(line_key), f"lam{tok(lamc, 3)}"]

    if mode in WAVE_MODES:
        parts.append(f"w{tok(lamlo, 2)}to{tok(lamhi, 2)}")
        parts.append(f"dw{tok(LAMSTEP_A, 3)}")
    if mode in VELO_MODES:
        parts.append(f"v{tok(VLO_KMS, 1)}to{tok(VHI_KMS, 1)}")
        parts.append(f"dv{tok(VSTEP_KMS, 3)}")
    if mode in {"ts2ima", "ts2vima"}:
        parts.append(f"dt{tok(TSTEP_D, 3)}")
    if mode in PHASE_MODES and orb is not None:
        parts.append(slug(str(orb["tag"])))
        parts.append(f"P{tok(float(orb['period']), 5)}")
        parts.append(f"T0{tok(float(orb['t0_mjd']), 5)}")
        parts.append(f"dph{tok(PHSTEP, 4)}")
    if mode == "ts2vphima":
        parts.append(f"ov{tok(OVERLAP_PERCENT, 3)}")

    parts.append("renorm" if RENORM else "norenorm")
    parts.append("medfilt" if MEDIAN_FILTER else "nomedfilt")
    return "__".join(parts) + ".fits"


def build_command(
    mode: str,
    out_path: Path,
    lamc: float,
    lamlo: float,
    lamhi: float,
    orb: dict[str, float] | None,
) -> list[str]:
    cmd = [
        sys.executable,
        str(PROJECT_ROOT / f"{mode}.py"),
        "--out",
        str(out_path),
        "--list",
        str(LIST_FILE),
        "--extrapolate-x",
    ]

    if mode in WAVE_MODES:
        cmd += ["--lamlo", f"{lamlo}", "--lamhi", f"{lamhi}", "--lamstep", f"{LAMSTEP_A}"]
    if mode in VELO_MODES:
        cmd += [
            "--lamc",
            f"{lamc}",
            "--vlo",
            f"{VLO_KMS}",
            "--vhi",
            f"{VHI_KMS}",
            "--vstep",
            f"{VSTEP_KMS}",
        ]
    if mode in {"ts2ima", "ts2vima"}:
        cmd += ["--tstep", f"{TSTEP_D}"]
    if mode in PHASE_MODES and orb is not None:
        cmd += [
            "--period",
            f"{orb['period']}",
            "--t0-mjd",
            f"{orb['t0_mjd']}",
            "--phstep",
            f"{PHSTEP}",
        ]
    if mode == "ts2vphima":
        cmd += ["--overlap", f"{OVERLAP_PERCENT}"]

    cmd.append("--renorm" if RENORM else "--no-renorm")
    cmd.append("--filter" if MEDIAN_FILTER else "--no-filter")
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


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not LIST_FILE.exists():
        raise FileNotFoundError(f"List file not found: {LIST_FILE}")

    jobs: list[tuple[str, Path, list[str]]] = []
    for line_key in LINE_NAMES:
        lamc = float(wavelength_for(line_key))
        lamlo = lamc - WAVE_HALF_WIDTH_A
        lamhi = lamc + WAVE_HALF_WIDTH_A

        for mode in MODES:
            if mode in PHASE_MODES:
                for orb in ORBITAL_SOLUTIONS:
                    out_name = make_output_name(mode, line_key, lamc, lamlo, lamhi, orb)
                    out_path = OUTPUT_DIR / out_name
                    cmd = build_command(mode, out_path, lamc, lamlo, lamhi, orb)
                    jobs.append((mode, out_path, cmd))
            else:
                out_name = make_output_name(mode, line_key, lamc, lamlo, lamhi, None)
                out_path = OUTPUT_DIR / out_name
                cmd = build_command(mode, out_path, lamc, lamlo, lamhi, None)
                jobs.append((mode, out_path, cmd))

    n_ok = 0
    n_skip = 0
    n_fail = 0

    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Total jobs: {len(jobs)}")

    for mode, out_path, cmd in jobs:
        if out_path.exists() and not OVERWRITE:
            print(f"SKIP [{mode}] {out_path.name}")
            n_skip += 1
            continue

        print(f"RUN  [{mode}] {out_path.name}")
        rc = run_one(cmd)
        if rc == 0:
            n_ok += 1
        else:
            n_fail += 1
            if STOP_ON_ERROR:
                break

    print("")
    print(f"Done. ok={n_ok} skip={n_skip} fail={n_fail}")
    return 0 if n_fail == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
