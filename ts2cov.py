#!/usr/bin/env python3
"""ts2cov.py — covariance matrix product (TS2COV analogue)."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import numpy as np
from astropy.io import fits

from ts2dynspec import (
    C_KMS,
    Params,
    configure_logging,
    cubic_spline_1d,
    make_x_grid,
    median_filter_y,
    read_list_file,
    read_time_mjd,
    read_wave_flux,
    renormalize_edges,
    velocity_from_lambda,
)


def build_velocity_rows(files: Sequence[Path], params: Params) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_grid, _ctype1, _cunit1 = make_x_grid(params)

    lamc = float(params.lamc)
    w1 = lamc * (1.0 + float(np.min(x_grid)) / C_KMS)
    w2 = lamc * (1.0 + float(np.max(x_grid)) / C_KMS)
    wlo, whi = (w1, w2) if w1 < w2 else (w2, w1)

    mjd_list: list[float] = []
    rows: list[np.ndarray] = []

    for fp in files:
        if not fp.exists():
            continue

        wave, flux, hdr = read_wave_flux(fp)
        mjd = read_time_mjd(hdr)

        m = np.isfinite(wave) & np.isfinite(flux) & (wave >= wlo) & (wave <= whi)
        if int(np.count_nonzero(m)) < 4:
            continue

        x_sel = velocity_from_lambda(wave[m], lamc)
        row = cubic_spline_1d(x_sel, flux[m], x_grid, extrapolate=params.extrapolate_x)
        mjd_list.append(float(mjd))
        rows.append(row)

    if not rows:
        raise RuntimeError("No usable spectra. Check list paths and velocity window.")

    mjd = np.asarray(mjd_list, dtype=np.float64)
    mat = np.vstack(rows).astype(np.float64)

    # Match C expectation: input must already be strictly ascending in time.
    if np.any(np.diff(mjd) <= 0.0):
        raise ValueError("Input list must be strictly ascending in time for ts2cov.")

    if params.renorm:
        mat = renormalize_edges(mat, params.edge_win)
    if params.median_filter:
        mat = median_filter_y(mat, params.half_med)

    return x_grid, mjd, mat


def covariance_like_c(mat: np.ndarray) -> np.ndarray:
    """
    Build covariance-like matrix following the TS2COV post-transform:
    - covariance from row samples
    - off-diagonals normalized by sqrt(var_i + var_j)
    - diagonals set to sqrt(var/2)
    """
    ny, nx = mat.shape
    if ny < 2:
        raise ValueError("Need at least 2 spectra to compute covariance.")

    finite = np.isfinite(mat)
    mean = np.nanmean(mat, axis=0)
    centered = np.where(finite, mat - mean[None, :], 0.0)

    counts = (finite.astype(np.int32).T @ finite.astype(np.int32)).astype(np.int32)
    sumprod = centered.T @ centered

    cov = np.full((nx, nx), np.nan, dtype=np.float64)
    valid = counts > 1
    cov[valid] = sumprod[valid] / (counts[valid] - 1.0)

    out = cov.copy()
    for i in range(nx):
        for j in range(i + 1, nx):
            denom = np.sqrt(cov[i, i] + cov[j, j]) if np.isfinite(cov[i, i]) and np.isfinite(cov[j, j]) else np.nan
            if np.isfinite(cov[i, j]) and np.isfinite(denom) and denom > 0.0:
                val = cov[i, j] / denom
            else:
                val = np.nan
            out[i, j] = val
            out[j, i] = val

    for i in range(nx):
        if np.isfinite(cov[i, i]) and cov[i, i] >= 0.0:
            out[i, i] = np.sqrt(cov[i, i] / 2.0)
        else:
            out[i, i] = np.nan

    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Compute covariance-like velocity matrix from a time series.")
    ap.add_argument("--out", type=Path, required=True, help="Output FITS image file.")
    ap.add_argument("--list", type=Path, required=True, help="Text file with one FITS spectrum path per line.")
    ap.add_argument("--lamc", type=float, required=True, help="Central wavelength Angstrom.")
    ap.add_argument("--vlo", type=float, required=True, help="Lower velocity limit km/s.")
    ap.add_argument("--vhi", type=float, required=True, help="Upper velocity limit km/s.")
    ap.add_argument("--vstep", type=float, default=5.0, help="Velocity step km/s. Default=5.")
    ap.add_argument("--renorm", action=argparse.BooleanOptionalAction, default=False,
                    help="Apply edge renormalization. Default=False.")
    ap.add_argument("--filter", dest="median_filter", action=argparse.BooleanOptionalAction, default=False,
                    help="Apply 3-point median filter along time. Default=False.")
    ap.add_argument("--extrapolate-x", action=argparse.BooleanOptionalAction, default=True,
                    help="Allow slight X-direction spline extrapolation. Default=True.")
    ap.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity.")
    args = ap.parse_args(argv)

    configure_logging(args.verbose)
    files = read_list_file(args.list)

    params = Params(
        xmode="velocity",
        ymode="time",
        lamc=float(args.lamc),
        vlo=float(args.vlo),
        vhi=float(args.vhi),
        vstep=float(args.vstep),
        renorm=bool(args.renorm),
        median_filter=bool(args.median_filter),
        extrapolate_x=bool(args.extrapolate_x),
    )

    x_grid, _mjd, mat = build_velocity_rows(files, params)
    cov = covariance_like_c(mat).astype(np.float32)

    h = fits.Header()
    h["CTYPE1"] = ("VELO", "x axis type")
    h["CUNIT1"] = ("km/s", "x axis unit")
    h["CRPIX1"] = (1.0, "reference pixel")
    h["CRVAL1"] = (float(x_grid[0]), "velocity at reference pixel")
    h["CDELT1"] = (float(x_grid[1] - x_grid[0]) if x_grid.size > 1 else 1.0, "velocity step")
    h["CTYPE2"] = ("VELO", "y axis type")
    h["CUNIT2"] = ("km/s", "y axis unit")
    h["CRPIX2"] = (1.0, "reference pixel")
    h["CRVAL2"] = (float(x_grid[0]), "velocity at reference pixel")
    h["CDELT2"] = (float(x_grid[1] - x_grid[0]) if x_grid.size > 1 else 1.0, "velocity step")
    h["LAMC"] = (float(params.lamc), "Central wavelength Angstrom")
    h["VLO"] = (float(np.min(x_grid)), "Lower velocity limit")
    h["VHI"] = (float(np.max(x_grid)), "Upper velocity limit")
    h["N_SPEC"] = (int(mat.shape[0]), "Number of spectra used")
    h["RENORM"] = (bool(params.renorm), "Edge renorm applied")
    h["MEDFILT"] = (bool(params.median_filter), "Median filter applied")

    fits.PrimaryHDU(data=cov, header=h).writeto(args.out, overwrite=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
