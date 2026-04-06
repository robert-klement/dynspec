#!/usr/bin/env python3
"""ts2vtab.py — velocity table product (TS2VTAB analogue)."""
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


def build_table_rows(files: Sequence[Path], params: Params) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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
        raise ValueError("Input list must be strictly ascending in time for ts2vtab.")

    if params.renorm:
        mat = renormalize_edges(mat, params.edge_win)
    if params.median_filter:
        mat = median_filter_y(mat, params.half_med)

    return x_grid, mjd, mat


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Convert a time series of spectra to velocity table output.")
    ap.add_argument("--out", type=Path, required=True, help="Output FITS table file.")
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

    x_grid, mjd, mat = build_table_rows(files, params)
    ny, nx = mat.shape

    cols = [fits.Column(name="jd24", format="D", array=mjd.astype(np.float64))]
    for i in range(nx):
        cols.append(fits.Column(name=f"l{i + 1:04d}", format="E", array=mat[:, i].astype(np.float32)))

    tab = fits.BinTableHDU.from_columns(cols, name="TS2VTAB")
    tab.header["VNPIX"] = (int(nx), "Number of velocity bins")
    tab.header["VSTART"] = (float(x_grid[0]), "Velocity start km/s")
    tab.header["VSTEP"] = (float(x_grid[1] - x_grid[0]) if nx > 1 else 0.0, "Velocity step km/s")
    tab.header["LAMC"] = (float(params.lamc), "Central wavelength Angstrom")
    tab.header["VLO"] = (float(np.min(x_grid)), "Lower velocity limit")
    tab.header["VHI"] = (float(np.max(x_grid)), "Upper velocity limit")
    tab.header["N_SPEC"] = (int(ny), "Number of spectra used")
    tab.header["RENORM"] = (bool(params.renorm), "Edge renorm applied")
    tab.header["MEDFILT"] = (bool(params.median_filter), "Median filter applied")

    fits.HDUList([fits.PrimaryHDU(), tab]).writeto(args.out, overwrite=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
