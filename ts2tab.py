#!/usr/bin/env python3
"""ts2tab.py — wavelength table product (TS2TAB analogue)."""
from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
from typing import Sequence

import numpy as np
from astropy.io import fits

from ts2dynspec import (
    Params,
    configure_logging,
    cubic_spline_1d,
    infer_lamstep_from_first_spectrum,
    make_x_grid,
    median_filter_y,
    read_list_file,
    read_time_mjd,
    read_wave_flux,
    renormalize_edges,
)


def build_table_rows(files: Sequence[Path], params: Params) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_grid, _ctype1, _cunit1 = make_x_grid(params)
    wlo = float(np.min(x_grid))
    whi = float(np.max(x_grid))

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

        row = cubic_spline_1d(wave[m], flux[m], x_grid, extrapolate=params.extrapolate_x)
        mjd_list.append(float(mjd))
        rows.append(row)

    if not rows:
        raise RuntimeError("No usable spectra. Check list paths and wavelength window.")

    mjd = np.asarray(mjd_list, dtype=np.float64)
    mat = np.vstack(rows).astype(np.float64)

    # Match C expectation: input must already be strictly ascending in time.
    if np.any(np.diff(mjd) <= 0.0):
        raise ValueError("Input list must be strictly ascending in time for ts2tab.")

    if params.renorm:
        mat = renormalize_edges(mat, params.edge_win)
    if params.median_filter:
        mat = median_filter_y(mat, params.half_med)

    return x_grid, mjd, mat


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Convert a time series of spectra to wavelength table output.")
    ap.add_argument("--out", type=Path, required=True, help="Output FITS table file.")
    ap.add_argument("--list", type=Path, required=True, help="Text file with one FITS spectrum path per line.")
    ap.add_argument("--lamlo", type=float, required=True, help="Lower wavelength bound Angstrom.")
    ap.add_argument("--lamhi", type=float, required=True, help="Upper wavelength bound Angstrom.")
    ap.add_argument(
        "--lamstep",
        type=float,
        default=None,
        help="Wavelength step Angstrom. Default: inferred from first usable spectrum.",
    )
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
        xmode="wavelength",
        ymode="time",
        lamlo=float(args.lamlo),
        lamhi=float(args.lamhi),
        lamstep=args.lamstep,
        renorm=bool(args.renorm),
        median_filter=bool(args.median_filter),
        extrapolate_x=bool(args.extrapolate_x),
    )

    if params.lamstep is None:
        params = replace(params, lamstep=infer_lamstep_from_first_spectrum(files))

    x_grid, mjd, mat = build_table_rows(files, params)
    ny, nx = mat.shape

    cols = [fits.Column(name="jd24", format="D", array=mjd.astype(np.float64))]
    for i in range(nx):
        cols.append(fits.Column(name=f"l{i + 1:04d}", format="E", array=mat[:, i].astype(np.float32)))

    tab = fits.BinTableHDU.from_columns(cols, name="TS2TAB")
    tab.header["LNPIX"] = (int(nx), "Number of wavelength bins")
    tab.header["LSTART"] = (float(x_grid[0]), "Wavelength start Angstrom")
    tab.header["LSTEP"] = (float(x_grid[1] - x_grid[0]) if nx > 1 else 0.0, "Wavelength step Angstrom")
    tab.header["XLO"] = (float(np.min(x_grid)), "Lower wavelength limit")
    tab.header["XHI"] = (float(np.max(x_grid)), "Upper wavelength limit")
    tab.header["N_SPEC"] = (int(ny), "Number of spectra used")
    tab.header["RENORM"] = (bool(params.renorm), "Edge renorm applied")
    tab.header["MEDFILT"] = (bool(params.median_filter), "Median filter applied")

    fits.HDUList([fits.PrimaryHDU(), tab]).writeto(args.out, overwrite=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
