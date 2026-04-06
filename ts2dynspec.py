#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ts2dynspec.py — Build dynamical spectra (time/phase x velocity/wavelength) from FITS 1D spectra.

Modes:
  - TS2VIMA   : --x velocity   --y time
  - TS2IMA    : --x wavelength --y time
  - TS2VPHIMA : --x velocity   --y phase
  - TS2PHIMA  : --x wavelength --y phase

Key fixes in this version
-------------------------
1) FITS header comments use ASCII only ("Angstrom", not "Å"), because FITS header comments
   must be standard printable ASCII.
2) Y-resampling is gap-aware and does NOT extrapolate outside the support of each segment.
   This avoids huge spline blow-ups when using coarse --tstep or --phstep.
3) Very narrow Y-segments (narrower than one output bin) are collapsed to the nearest output
   row instead of disappearing or exploding.
4) Phase modes use C-style bin averaging in phase; ts2vphima supports overlap in percent.

Dependencies: numpy, scipy, astropy
"""

from __future__ import annotations

import argparse
import logging
import math
import warnings
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional, Sequence, Tuple, Literal

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from astropy.constants import c
from astropy import units as u
from scipy.interpolate import CubicSpline, interp1d

C_KMS = float(c.to_value("km/s"))
log = logging.getLogger("ts2dynspec")

XMode = Literal["velocity", "wavelength"]
YMode = Literal["time", "phase"]


@dataclass(frozen=True)
class Params:
    # Mode selection
    xmode: XMode = "velocity"
    ymode: YMode = "time"

    # Velocity mode
    lamc: Optional[float] = None   # Angstrom
    vlo: Optional[float] = -1000    # km/s
    vhi: Optional[float] = 1000     # km/s
    vstep: float = 5.0             # km/s

    # Wavelength mode
    lamlo: Optional[float] = None   # Angstrom
    lamhi: Optional[float] = None   # Angstrom
    lamstep: Optional[float] = None # Angstrom

    # Time mode
    tstep: float = 1.0              # days

    # Phase mode
    period: Optional[float] = None  # days
    t0_mjd: Optional[float] = None  # MJD
    phstep: float = 0.01
    overlap: float = 0.0

    # Processing
    renorm: bool = False
    median_filter: bool = False

    # MIDAS-like defaults
    edge_win: int = 25
    half_med: int = 1
    gap_factor: float = 6.0

    # Numerics
    extrapolate_x: bool = True     # X-direction rebinning can extrapolate slightly if desired
    fill_value: float = np.nan
    average_duplicates: bool = True


def configure_logging(verbosity: int) -> None:
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG
    logging.basicConfig(level=level, format="%(levelname)s: %(message)s")


# -------------------------------------------------------------------
# FITS reading
# -------------------------------------------------------------------

def read_time_mjd(header: fits.Header) -> float:
    """Extract observation time in MJD from common FITS keywords."""
    for key in ("MJD-OBS", "MJD_OBS", "MJDOBS", "MJD"):
        if key in header:
            return float(header[key])

    for key in ("JD-OBS", "JD_OBS", "JD"):
        if key in header:
            return float(Time(float(header[key]), format="jd", scale="utc").mjd)

    for key in ("DATE-OBS", "DATEOBS"):
        if key in header:
            s = str(header[key]).strip()
            for fmt in ("isot", "fits", "iso"):
                try:
                    return float(Time(s, format=fmt, scale="utc").mjd)
                except Exception:
                    pass
            return float(Time(s, scale="utc").mjd)

    raise KeyError("No recognizable time keyword found (tried MJD-OBS/JD-OBS/DATE-OBS).")


def wave_to_angstrom(wave: np.ndarray, unit_raw, source: Path) -> np.ndarray:
    """
    Convert wavelength axis to Angstrom with a guard against common bad CUNIT metadata.

    Some spectra are tagged as nm while CRVAL/CDELT are already in Angstrom-scale values
    (~5000-10000). Converting those would shift them by x10 and miss optical lines.
    """
    if unit_raw is None:
        return np.asarray(wave, dtype=np.float64)

    unit_str = str(unit_raw).strip()
    if not unit_str:
        return np.asarray(wave, dtype=np.float64)

    try:
        unit = u.Unit(unit_str)
    except Exception:
        log.warning("%s: unrecognized wavelength unit %r; assuming Angstrom.", source.name, unit_str)
        return np.asarray(wave, dtype=np.float64)

    # Heuristic for mis-labeled XSHOOTER-like products:
    # CUNIT1='nm' but CRVAL/CDELT are already in Angstrom scale.
    if unit == u.nm:
        finite = np.isfinite(wave)
        if np.any(finite):
            med = float(np.nanmedian(np.asarray(wave, dtype=np.float64)[finite]))
            if med > 3000.0:
                log.info(
                    "%s: CUNIT1=%r but wavelength values look Angstrom-like; using raw values.",
                    source.name,
                    unit_str,
                )
                return np.asarray(wave, dtype=np.float64)

    return np.asarray((wave * unit).to_value(u.AA), dtype=np.float64)


def read_wave_flux(path: Path) -> Tuple[np.ndarray, np.ndarray, fits.Header]:
    """
    Read wavelength [Angstrom], flux, and the source header from a FITS file.

    Supported:
      1) 1D image HDU with linear WCS (CRVAL1/CDELT1, optional CRPIX1/CUNIT1)
      2) Table HDU with common wavelength/flux column names
    """
    with fits.open(path, memmap=False) as hdul:
        # Prefer a 1D image spectrum
        for hdu in hdul:
            if isinstance(hdu, (fits.PrimaryHDU, fits.ImageHDU)) and hdu.data is not None:
                data = np.asarray(hdu.data)
                if data.ndim == 1 and data.size > 0:
                    flux = np.asarray(data, dtype=np.float64).copy()
                    hdr = hdu.header

                    if "CRVAL1" not in hdr or "CDELT1" not in hdr:
                        raise KeyError(f"{path}: missing CRVAL1/CDELT1 for wavelength axis")

                    crval = float(hdr["CRVAL1"])
                    cdelt = float(hdr["CDELT1"])
                    crpix = float(hdr.get("CRPIX1", 1.0))

                    pix_fits = np.arange(flux.size, dtype=np.float64) + 1.0
                    wave = crval + (pix_fits - crpix) * cdelt

                    cunit = hdr.get("CUNIT1", None)
                    if cunit:
                        wave = wave_to_angstrom(wave, cunit, path)

                    return np.asarray(wave, dtype=np.float64), flux, hdr

        # Otherwise try a table
        for hdu in hdul:
            if isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                tab = Table(hdu.data)
                hdr = hdu.header
                colmap = {name.upper(): name for name in tab.colnames}

                wcol = None
                for cand in ("WAVE", "WAVELENGTH", "LAMBDA", "WL", "LAMBDA_A", "LAMBDA_ANG", "LAMBDAANG"):
                    if cand in colmap:
                        wcol = colmap[cand]
                        break

                fcol = None
                for cand in ("FLUX", "SPEC", "F", "Y", "NORMFLUX"):
                    if cand in colmap:
                        fcol = colmap[cand]
                        break

                if wcol and fcol:
                    wave = np.asarray(tab[wcol], dtype=np.float64)
                    flux = np.asarray(tab[fcol], dtype=np.float64)
                    wu = getattr(tab[wcol], "unit", None)
                    if wu:
                        wave = wave_to_angstrom(wave, wu, path)
                    return np.asarray(wave, dtype=np.float64), np.asarray(flux, dtype=np.float64), hdr

    raise ValueError(f"{path}: could not find 1D image spectrum or wavelength/flux columns in a table")


# -------------------------------------------------------------------
# Axis helpers
# -------------------------------------------------------------------

def velocity_from_lambda(wave_aa: np.ndarray, lamc_aa: float) -> np.ndarray:
    """Non-relativistic Doppler velocity relative to lamc."""
    return C_KMS * (wave_aa / lamc_aa - 1.0)


def infer_lamstep_from_first_spectrum(files: Sequence[Path]) -> float:
    """
    Infer default wavelength step from the first usable input spectrum.
    Mirrors C behavior where wavelength modes use the first frame step unless overridden.
    """
    for fp in files:
        if not fp.exists():
            continue
        try:
            wave, _flux, _hdr = read_wave_flux(fp)
        except Exception:
            continue

        w = np.asarray(wave, dtype=np.float64)
        m = np.isfinite(w)
        if np.count_nonzero(m) < 2:
            continue

        dw = np.diff(w[m])
        dw = dw[np.isfinite(dw)]
        dw = np.abs(dw[dw != 0.0])
        if dw.size == 0:
            continue

        step = float(np.nanmedian(dw))
        if np.isfinite(step) and step > 0:
            return step

    raise RuntimeError(
        "Could not infer default wavelength step from input spectra; please provide --lamstep."
    )


def make_x_grid(params: Params) -> Tuple[np.ndarray, str, str]:
    """Build the common X grid and return (x_grid, CTYPE1, CUNIT1)."""
    if params.xmode == "velocity":
        if params.lamc is None or params.vlo is None or params.vhi is None:
            raise ValueError("velocity mode requires --lamc, --vlo, and --vhi")
        v_lo, v_hi = (params.vlo, params.vhi) if params.vhi >= params.vlo else (params.vhi, params.vlo)
        n = int(math.floor((v_hi - v_lo) / params.vstep)) + 1
        x = v_lo + np.arange(n, dtype=np.float64) * params.vstep
        return x, "VELO", "km/s"

    # wavelength mode
    if params.lamlo is None or params.lamhi is None:
        # Convenience: derive wavelength limits from lamc + velocity limits if available
        if params.lamc is None or params.vlo is None or params.vhi is None:
            raise ValueError("wavelength mode requires --lamlo/--lamhi OR --lamc/--vlo/--vhi")
        w1 = params.lamc * (1.0 + params.vlo / C_KMS)
        w2 = params.lamc * (1.0 + params.vhi / C_KMS)
        lamlo, lamhi = (w1, w2) if w1 < w2 else (w2, w1)
    else:
        lamlo, lamhi = (params.lamlo, params.lamhi) if params.lamhi >= params.lamlo else (params.lamhi, params.lamlo)

    if params.lamstep is None:
        if params.lamc is not None:
            lamstep = params.lamc * (params.vstep / C_KMS)
        else:
            lamstep = 0.05
        log.warning("No --lamstep given; using %.6f Angstrom.", lamstep)
    else:
        lamstep = float(params.lamstep)

    n = int(math.floor((lamhi - lamlo) / lamstep)) + 1
    x = lamlo + np.arange(n, dtype=np.float64) * lamstep
    return x, "WAVE", "Angstrom"


def make_y_grid(params: Params, mjd: np.ndarray) -> Tuple[np.ndarray, str, str]:
    """Build the output Y grid and return (y_grid, CTYPE2, CUNIT2)."""
    if params.ymode == "time":
        t0 = float(np.min(mjd))
        t1 = float(np.max(mjd))
        n = int(math.floor((t1 - t0) / params.tstep)) + 1
        y = t0 + np.arange(n, dtype=np.float64) * params.tstep
        return y, "MJD", "d"

    if params.period is None or params.t0_mjd is None:
        raise ValueError("phase mode requires --period and --t0-mjd")

    n = int(math.floor(1.0 / params.phstep))
    y = np.arange(n, dtype=np.float64) * params.phstep
    return y, "PHASE", ""


def compute_y_samples(params: Params, mjd: np.ndarray) -> np.ndarray:
    """Compute the Y coordinate for each observation: MJD or orbital phase in [0,1)."""
    if params.ymode == "time":
        return mjd.copy()
    return np.remainder((mjd - float(params.t0_mjd)) / float(params.period), 1.0)


# -------------------------------------------------------------------
# Interpolation / filtering helpers
# -------------------------------------------------------------------

def cubic_spline_1d(x: np.ndarray, y: np.ndarray, x_new: np.ndarray, extrapolate: bool) -> np.ndarray:
    """
    X-direction cubic spline interpolation using bc_type='clamped' (zero first derivative at ends).

    We sort x, drop duplicate x points, and ignore non-finite data.
    """
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < 2:
        return np.full_like(x_new, np.nan, dtype=np.float64)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    uniq = np.concatenate(([True], np.diff(x) != 0))
    x = x[uniq]
    y = y[uniq]
    if x.size < 2:
        return np.full_like(x_new, np.nan, dtype=np.float64)

    cs = CubicSpline(x, y, bc_type="clamped", extrapolate=extrapolate)
    return cs(x_new)


def renormalize_edges(img: np.ndarray, edge_win: int) -> np.ndarray:
    """
    Renormalize each row by a linear continuum defined by edge medians.
    """
    out = img.copy()
    n_rows, n_x = out.shape
    w = int(edge_win)
    if n_x < 2 * w + 2:
        return out

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        left = np.nanmedian(out[:, :w], axis=1)
        right = np.nanmedian(out[:, -w:], axis=1)

    a = (right - left) / (n_x - w)
    b = left - a * (w / 2.0)

    x = np.arange(n_x, dtype=np.float64)[None, :]
    cont = a[:, None] * x + b[:, None]
    cont[~np.isfinite(cont) | (cont == 0)] = np.nan
    return out / cont


def median_filter_y(img: np.ndarray, half_med: int) -> np.ndarray:
    """
    Median filter along Y with window size (2*half_med+1), NaN-safe.
    """
    out = img.copy()
    h = int(half_med)
    if h <= 0 or img.shape[0] < 2 * h + 1:
        return out

    windows = np.lib.stride_tricks.sliding_window_view(img, window_shape=(2 * h + 1), axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        med = np.nanmedian(windows, axis=-1)
    out[h:-h, :] = med
    return out


def segments_from_gaps(y: np.ndarray, y_step: float, gap_factor: float):
    """Return contiguous segments [i0, i1] that do not cross large gaps."""
    if y.size == 0:
        return []
    gaps = np.diff(y)
    breaks = np.where(gaps > gap_factor * y_step)[0]
    segs = []
    start = 0
    for b in breaks:
        segs.append((start, b))
        start = b + 1
    segs.append((start, y.size - 1))
    return segs


def average_duplicate_y(
    y: np.ndarray,
    rows: np.ndarray,
    files: Sequence[Path],
    average: bool,
) -> Tuple[np.ndarray, np.ndarray, list[Path]]:
    """
    Ensure Y is strictly increasing for spline interpolation.
    Average exact duplicates if requested.
    """
    if y.size <= 1:
        return y, rows, list(files)

    dup = np.where(np.diff(y) == 0)[0]
    if dup.size == 0:
        return y, rows, list(files)

    if not average:
        i = int(dup[0])
        raise ValueError(f"Duplicate y-samples at indices {i} and {i+1} (y={y[i]}).")

    uniq, inv = np.unique(y, return_inverse=True)
    out_rows = np.empty((uniq.size, rows.shape[1]), dtype=np.float64)
    out_files: list[Path] = []

    for k in range(uniq.size):
        sel = np.where(inv == k)[0]
        out_rows[k] = np.nanmean(rows[sel, :], axis=0)
        out_files.append(files[sel[0]])

    return uniq, out_rows, out_files


def phase_bin_average(
    phase: np.ndarray,
    rows: np.ndarray,
    *,
    phstep: float,
    fill_value: float,
    overlap_percent: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    C-style phase binning:
    - ts2phima: average into phase bins on [0,1).
    - ts2vphima: optionally extend by overlap and duplicate edge samples.
    """
    if phstep <= 0:
        raise ValueError("--phstep must be > 0 for phase modes")

    phase_bin = int(round(1.0 / phstep))
    if phase_bin < 1:
        raise ValueError("invalid phase binning: phstep produced <1 bins")

    pstep = 1.0 / float(phase_bin)
    over = max(0.0, float(overlap_percent)) / 100.0

    pstart = 0.0 - over
    if over > 0.0:
        nbin = int(phase_bin * (1.0 + 2.0 * over))
    else:
        nbin = phase_bin
    nbin = max(1, nbin)

    y_grid = pstart + np.arange(nbin, dtype=np.float64) * pstep

    sums = np.zeros((nbin, rows.shape[1]), dtype=np.float64)
    counts = np.zeros((nbin, rows.shape[1]), dtype=np.int32)

    for t in range(phase.size):
        p = float(phase[t])
        row = rows[t]
        finite = np.isfinite(row)
        if not np.any(finite):
            continue

        b = int((p - pstart) / pstep)
        if 0 <= b < nbin:
            sums[b, finite] += row[finite]
            counts[b, finite] += 1

        if over > 0.0:
            b2: Optional[int] = None
            if p >= (1.0 - over):
                b2 = b - phase_bin
            elif p <= over:
                b2 = b + phase_bin

            if b2 is not None and 0 <= b2 < nbin:
                sums[b2, finite] += row[finite]
                counts[b2, finite] += 1

    img = np.full((nbin, rows.shape[1]), float(fill_value), dtype=np.float64)
    good = counts > 0
    img[good] = sums[good] / counts[good]
    return img, y_grid


def resample_rows_to_uniform_y(
    y_samples: np.ndarray,
    rows: np.ndarray,
    y_grid: np.ndarray,
    *,
    y_step: float,
    gap_factor: float,
    fill_value: float,
) -> np.ndarray:
    """
    Resample rows defined at y_samples onto a uniform y_grid.

    Important behavior:
    - No Y extrapolation is performed. This avoids huge spline blow-ups when a segment
      is much narrower than the chosen output step (e.g. very coarse --tstep).
    - If a segment is narrower than one output bin and contains no y_grid points inside it,
      the segment is collapsed to the nearest grid row using the mean row.
    - Large gaps remain empty.
    """
    img = np.full((y_grid.size, rows.shape[1]), fill_value, dtype=np.float64)
    segs = segments_from_gaps(y_samples, y_step, gap_factor)

    for i0, i1 in segs:
        ys = y_samples[i0:i1 + 1]
        rs = rows[i0:i1 + 1, :]

        # Grid rows that actually fall inside the support of this segment
        inside = np.where((y_grid >= ys[0]) & (y_grid <= ys[-1]))[0]

        # Single sample: put it on the nearest output row
        if ys.size == 1:
            j = int(np.argmin(np.abs(y_grid - ys[0])))
            img[j, :] = rs[0, :]
            continue

        # Narrow segment (no output row falls inside the support): collapse to nearest row
        if inside.size == 0:
            j = int(np.argmin(np.abs(y_grid - np.mean(ys))))
            img[j, :] = np.nanmean(rs, axis=0)
            continue

        yg = y_grid[inside]

        # C-style interpolation choice in time direction:
        #   2 points -> linear
        #   >2 points -> cubic spline
        # We evaluate per X-column to stay robust if some columns contain NaNs.
        vals = np.full((yg.size, rs.shape[1]), np.nan, dtype=np.float64)
        use_linear = ys.size == 2
        for x in range(rs.shape[1]):
            col = rs[:, x]
            good = np.isfinite(col)
            ng = int(np.count_nonzero(good))
            if ng == 0:
                continue
            if ng == 1:
                j = int(np.argmin(np.abs(yg - ys[good][0])))
                vals[j, x] = col[good][0]
                continue

            ysg = ys[good]
            colg = col[good]

            if use_linear or ng == 2:
                f = interp1d(
                    ysg,
                    colg,
                    kind="linear",
                    bounds_error=False,
                    fill_value=np.nan,
                )
                vals[:, x] = f(yg)
            else:
                cs = CubicSpline(ysg, colg, bc_type="clamped", extrapolate=False)
                vals[:, x] = cs(yg)

        good = np.isfinite(vals)
        for k, j in enumerate(inside):
            if np.any(good[k]):
                img[j, good[k]] = vals[k, good[k]]

    return img


# -------------------------------------------------------------------
# Main builder
# -------------------------------------------------------------------

def build_dynspec(files: Sequence[Path], params: Params) -> Tuple[np.ndarray, Table, fits.Header]:
    """
    Build the final dynamical spectrum image and return:
      img, obs_table, output_header
    """
    if params.xmode == "wavelength" and params.lamstep is None:
        inferred_lamstep = infer_lamstep_from_first_spectrum(files)
        params = replace(params, lamstep=inferred_lamstep)
        log.info(
            "No --lamstep provided; using first-spectrum wavelength step %.8g Angstrom.",
            inferred_lamstep,
        )

    x_grid, ctype1, cunit1 = make_x_grid(params)

    # Wavelength pre-slicing limits for efficiency/stability
    if params.xmode == "velocity":
        lamc = float(params.lamc)
        w1 = lamc * (1.0 + float(np.min(x_grid)) / C_KMS)
        w2 = lamc * (1.0 + float(np.max(x_grid)) / C_KMS)
        wlo, whi = (w1, w2) if w1 < w2 else (w2, w1)
    else:
        wlo, whi = float(np.min(x_grid)), float(np.max(x_grid))

    mjd_list: list[float] = []
    rows: list[np.ndarray] = []
    used: list[Path] = []

    # Read and rebin each input spectrum
    for fp in files:
        if not fp.exists():
            log.warning("Missing file: %s (skipping)", fp)
            continue

        wave, flux, hdr = read_wave_flux(fp)
        mjd = read_time_mjd(hdr)

        m = np.isfinite(wave) & np.isfinite(flux) & (wave >= wlo) & (wave <= whi)
        n_in_window = int(np.count_nonzero(m))
        if n_in_window < 4:
            log.warning("%s: too few points in window [%.3f, %.3f] Angstrom; skipping.", fp.name, wlo, whi)
            continue

        wave_sel = wave[m]
        flux_sel = flux[m]

        if params.xmode == "velocity":
            x_sel = velocity_from_lambda(wave_sel, float(params.lamc))
        else:
            x_sel = wave_sel

        row = cubic_spline_1d(x_sel, flux_sel, x_grid, extrapolate=params.extrapolate_x)

        mjd_list.append(float(mjd))
        rows.append(row)
        used.append(fp)

    if not rows:
        raise RuntimeError("No usable spectra. Check list paths and wavelength window.")

    mjd = np.array(mjd_list, dtype=np.float64)
    rows_arr = np.vstack(rows).astype(np.float64)

    # Compute Y samples
    y_samples = compute_y_samples(params, mjd)
    if params.ymode == "phase":
        ctype2, cunit2 = "PHASE", ""
        y_grid = np.empty(0, dtype=np.float64)
    else:
        y_grid, ctype2, cunit2 = make_y_grid(params, mjd)
    y_step = float(params.tstep)

    # Sort by Y
    order = np.argsort(y_samples)
    y_samples = y_samples[order]
    mjd = mjd[order]
    rows_arr = rows_arr[order, :]
    used = [used[i] for i in order]

    # ------------------------------------------------------------------
    # Build the OBS provenance table NOW, before phase-wrap duplication
    # and before averaging duplicate Y values.
    # ------------------------------------------------------------------
    obs = Table({
        "mjd": mjd.copy(),
        "y": y_samples.copy(),
        "filename": [str(p) for p in used],
    })
    obs.rename_column("y", "phase_sorted" if params.ymode == "phase" else "mjd_sorted")

    if params.ymode == "phase":
        # C order for phase modes: sort by phase, optional renorm/filter, then phase-bin averaging.
        if params.renorm:
            rows_arr = renormalize_edges(rows_arr, params.edge_win)
        if params.median_filter:
            rows_arr = median_filter_y(rows_arr, params.half_med)

        over = float(params.overlap) if params.xmode == "velocity" else 0.0
        img, y_grid = phase_bin_average(
            y_samples,
            rows_arr,
            phstep=float(params.phstep),
            fill_value=float(params.fill_value),
            overlap_percent=over,
        )
    else:
        # Time modes keep interpolation/resampling behavior unchanged.
        y_samples, rows_arr, used = average_duplicate_y(y_samples, rows_arr, used, params.average_duplicates)

        img = resample_rows_to_uniform_y(
            y_samples,
            rows_arr,
            y_grid,
            y_step=y_step,
            gap_factor=float(params.gap_factor),
            fill_value=float(params.fill_value),
        )

        # Convert non-NaN fill values to NaN for robust stats/filtering, then restore.
        restore_fill = not np.isnan(params.fill_value)
        if restore_fill:
            img = img.copy()
            img[np.isclose(img, params.fill_value)] = np.nan

        if params.renorm:
            img = renormalize_edges(img, params.edge_win)

        if params.median_filter:
            img = median_filter_y(img, params.half_med)

        if restore_fill:
            img = np.where(np.isfinite(img), img, float(params.fill_value))

    # Output FITS header (ASCII-only comments)
    hdr = fits.Header()

    hdr["CTYPE1"] = (ctype1, "X axis type")
    hdr["CUNIT1"] = (cunit1, "X axis unit")
    hdr["CRPIX1"] = (1.0, "Reference pixel 1-based")
    hdr["CRVAL1"] = (float(x_grid[0]), "X at reference pixel")
    hdr["CDELT1"] = (float(x_grid[1] - x_grid[0]) if x_grid.size > 1 else 1.0, "X step")

    hdr["CTYPE2"] = (ctype2, "Y axis type")
    hdr["CUNIT2"] = (cunit2, "Y axis unit")
    hdr["CRPIX2"] = (1.0, "Reference pixel 1-based")
    hdr["CRVAL2"] = (float(y_grid[0]), "Y at reference pixel")
    hdr["CDELT2"] = (float(y_grid[1] - y_grid[0]) if y_grid.size > 1 else 1.0, "Y step")

    hdr["XMODE"] = (params.xmode, "velocity or wavelength")
    hdr["YMODE"] = (params.ymode, "time or phase")
    hdr["RENORM"] = (bool(params.renorm), "Edge renormalization applied")
    hdr["MEDFILT"] = (bool(params.median_filter), "Median filter along Y applied")
    hdr["GAPFAC"] = (float(params.gap_factor), "No interp across gaps > GAPFAC*step")
    hdr["AVGDUP"] = (bool(params.average_duplicates), "Average identical y samples")
    hdr["N_SPEC"] = (int(len(obs)), "Spectra used after processing")

    if params.xmode == "velocity":
        hdr["LAMC"] = (float(params.lamc), "Central wavelength Angstrom")
        hdr["VLO"] = (float(params.vlo), "Lower velocity limit km/s")
        hdr["VHI"] = (float(params.vhi), "Upper velocity limit km/s")
        hdr["VSTEP"] = (float(params.vstep), "Velocity step km/s")
    else:
        hdr["LAMLO"] = (float(np.min(x_grid)), "Lower wavelength bound Angstrom")
        hdr["LAMHI"] = (float(np.max(x_grid)), "Upper wavelength bound Angstrom")
        hdr["LAMSTP"] = (float(x_grid[1] - x_grid[0]) if x_grid.size > 1 else 0.0, "Wavelength step Angstrom")

    if params.ymode == "phase":
        hdr["PERIOD"] = (float(params.period), "Orbital period days")
        hdr["T0MJD"] = (float(params.t0_mjd), "Reference epoch MJD")
        hdr["PHSTEP"] = (float(y_grid[1] - y_grid[0]) if y_grid.size > 1 else float(params.phstep), "Phase step")
        hdr["PHLO"] = (float(y_grid[0]), "Phase lower bound")
        hdr["PHHI"] = (float(y_grid[-1]), "Phase upper bound")
        if params.xmode == "velocity":
            hdr["OVERLAP"] = (float(params.overlap), "Phase overlap percent")
    else:
        hdr["TSTEP"] = (float(params.tstep), "Time step days")

    # Suggested display cuts: robust percentiles, lower cut forced >= 0
    finite = np.isfinite(img)
    if np.any(finite):
        vals = img[finite]
        p1, p99 = np.percentile(vals, [1, 99])
        hdr["VCUTLO"] = (max(0.0, float(p1)), "Suggested lower display cut")
        hdr["VCUTHI"] = (float(p99), "Suggested upper display cut")

    return img, obs, hdr


def write_fits(outpath: Path, img: np.ndarray, hdr: fits.Header, obs: Table) -> None:
    """Write the final image + OBS table."""
    hdus = [
        fits.PrimaryHDU(data=img.astype(np.float32), header=hdr),
        fits.BinTableHDU(obs, name="OBS"),
    ]
    fits.HDUList(hdus).writeto(outpath, overwrite=True)


# -------------------------------------------------------------------
# Validation / self-test / CLI
# -------------------------------------------------------------------

def read_list_file(list_path: Path) -> list[Path]:
    """Read list file; relative paths are resolved relative to the list file directory."""
    base = list_path.parent
    files: list[Path] = []
    for line in list_path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        p = Path(s).expanduser()
        if not p.is_absolute():
            p = (base / p).resolve()
        files.append(p)
    return files


def validate_files(files: Sequence[Path], limit: Optional[int] = None) -> int:
    """Validate that each file can be read and has a usable time keyword."""
    n_ok = 0
    n_fail = 0
    subset = files[: (limit or len(files))]
    for fp in subset:
        try:
            if not fp.exists():
                raise FileNotFoundError(str(fp))
            wave, flux, hdr = read_wave_flux(fp)
            _ = read_time_mjd(hdr)
            if wave.size != flux.size or wave.size < 2:
                raise ValueError(f"wave/flux mismatch or too short: {wave.size} vs {flux.size}")
            n_ok += 1
        except Exception as e:
            n_fail += 1
            print(f"FAIL: {fp} -> {type(e).__name__}: {e}")
    print(f"Validated {n_ok+n_fail} file(s): {n_ok} OK, {n_fail} FAIL")
    return 0 if n_fail == 0 else 2


def self_test() -> int:
    """
    Basic smoke test:
      - builds synthetic datasets from a template FITS spectrum in the current directory
      - exercises all four modes
    """
    configure_logging(2)
    work = Path("_ts2dynspec_selftest")
    work.mkdir(exist_ok=True)

    template = Path("spectra/ANCol_BeSSEch_VIS_60555_2930.fits")
    if not template.exists():
        cand = sorted(Path("spectra").glob("*.fits"))
        if not cand:
            log.error("Self-test needs any template *.fits spectrum in the spectra/ directory.")
            return 2
        template = cand[0]
        log.warning("Default template not found; using %s", template)

    mjds = [60555.10, 60555.60, 60556.10, 60556.60, 60560.10, 60560.60]
    files = []
    for i, mjd in enumerate(mjds):
        out = work / f"spec_{i:02d}.fits"
        with fits.open(template, memmap=False) as hdul:
            hdul[0].header["MJD-OBS"] = float(mjd)
            hdul.writeto(out, overwrite=True)
        files.append(out)

    tests = [
        ("time_velo", Params(xmode="velocity", ymode="time", lamc=6562.8, vlo=-400, vhi=400, vstep=5, tstep=0.5)),
        ("time_wave", Params(xmode="wavelength", ymode="time", lamlo=5000, lamhi=5050, lamstep=0.05, tstep=0.5)),
        ("phase_velo", Params(xmode="velocity", ymode="phase", lamc=6562.8, vlo=-400, vhi=400, vstep=5,
                              period=5.0, t0_mjd=60555.0, phstep=0.02)),
        ("phase_wave", Params(xmode="wavelength", ymode="phase", lamlo=5000, lamhi=5050, lamstep=0.05,
                              period=5.0, t0_mjd=60555.0, phstep=0.02)),
    ]

    for name, params in tests:
        img, obs, hdr = build_dynspec(files, params)
        assert img.ndim == 2
        assert np.any(np.isfinite(img)), f"{name}: no finite pixels"
        write_fits(work / f"{name}.fits", img, hdr, obs)
        log.info("Self-test wrote %s", work / f"{name}.fits")

    log.info("Self-test OK.")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Build dynamical spectra (time/phase x velocity/wavelength) from FITS 1D spectra."
    )
    ap.add_argument("--out", type=Path, help="Output FITS file.")
    ap.add_argument("--list", type=Path, help="Text file with one FITS spectrum path per line.")

    ap.add_argument("--x", dest="xmode", choices=["velocity", "wavelength"], default="velocity",
                    help="X axis mode. Default=velocity.")
    ap.add_argument("--y", dest="ymode", choices=["time", "phase"], default="time",
                    help="Y axis mode. Default=time.")

    ap.add_argument("--validate-only", action="store_true",
                    help="Only validate reading FITS and extracting wavelength/time.")
    ap.add_argument("--validate-limit", type=int, default=None,
                    help="If set, only validate the first N spectra.")
    ap.add_argument("--self-test", action="store_true",
                    help="Run the built-in smoke test and exit.")

    # Velocity mode args
    ap.add_argument("--lamc", type=float, default=None, help="Central wavelength Angstrom (velocity mode).")
    ap.add_argument("--vlo", type=float, default=None, help="Lower velocity limit km/s (velocity mode).")
    ap.add_argument("--vhi", type=float, default=None, help="Upper velocity limit km/s (velocity mode).")
    ap.add_argument("--vstep", type=float, default=5.0, help="Velocity step km/s. Default=5.")

    # Wavelength mode args
    ap.add_argument("--lamlo", type=float, default=None, help="Lower wavelength bound Angstrom (wavelength mode).")
    ap.add_argument("--lamhi", type=float, default=None, help="Upper wavelength bound Angstrom (wavelength mode).")
    ap.add_argument(
        "--lamstep",
        type=float,
        default=None,
        help="Wavelength step Angstrom. Default: inferred from first usable spectrum.",
    )

    # Time / phase args
    ap.add_argument("--tstep", type=float, default=0.5, help="Time step days. Default=0.5.")
    ap.add_argument("--period", type=float, default=None, help="Orbital period days (phase mode).")
    ap.add_argument("--t0-mjd", dest="t0_mjd", type=float, default=None, help="Reference epoch MJD (phase mode).")
    ap.add_argument("--phstep", type=float, default=0.01, help="Phase step. Default=0.01.")
    ap.add_argument("--overlap", type=float, default=0.0,
                    help="Phase overlap in percent for velocity+phase mode (ts2vphima). Default=0.")

    # Processing flags
    ap.add_argument("--renorm", action=argparse.BooleanOptionalAction, default=False,
                    help="Apply edge renormalization (changes absolute flux scale). Default=False.")
    ap.add_argument("--filter", dest="median_filter", action=argparse.BooleanOptionalAction, default=False,
                    help="Apply 3-point median filter along Y (smooths along time/phase). Default=False.")

    # Numerics
    ap.add_argument("--gap-factor", type=float, default=6.0,
                    help="Do not interpolate across gaps > gap-factor*step. Default=6.")
    ap.add_argument("--fill", default="nan",
                    help="Fill value for missing pixels: 'nan' or a number. Default=nan.")
    ap.add_argument("--extrapolate-x", action=argparse.BooleanOptionalAction, default=True,
                    help="Allow slight X-direction spline extrapolation. Default=True.")
    ap.add_argument("--average-duplicates", action=argparse.BooleanOptionalAction, default=True,
                    help="Average identical y-samples instead of error. Default=True.")

    ap.add_argument("-v", "--verbose", action="count", default=0,
                    help="Increase verbosity (-v, -vv).")

    args = ap.parse_args(argv)
    configure_logging(args.verbose)

    if args.self_test:
        return self_test()

    if args.list is None:
        ap.error("Missing required argument: --list")

    files = read_list_file(args.list)

    if args.validate_only:
        return validate_files(files, limit=args.validate_limit)

    if args.out is None:
        ap.error("Missing required argument: --out")

    fill_value = np.nan if str(args.fill).lower() == "nan" else float(args.fill)

    params = Params(
        xmode=args.xmode,
        ymode=args.ymode,
        lamc=args.lamc,
        vlo=args.vlo,
        vhi=args.vhi,
        vstep=args.vstep,
        lamlo=args.lamlo,
        lamhi=args.lamhi,
        lamstep=args.lamstep,
        tstep=args.tstep,
        period=args.period,
        t0_mjd=args.t0_mjd,
        phstep=args.phstep,
        overlap=args.overlap,
        renorm=bool(args.renorm),
        median_filter=bool(args.median_filter),
        gap_factor=float(args.gap_factor),
        fill_value=fill_value,
        extrapolate_x=bool(args.extrapolate_x),
        average_duplicates=bool(args.average_duplicates),
    )

    img, obs, hdr = build_dynspec(files, params)
    write_fits(args.out, img, hdr, obs)
    log.info("Wrote %s", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
