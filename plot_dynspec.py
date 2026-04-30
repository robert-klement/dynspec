#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_dynspec.py — Plot a dynamical spectrum FITS (time/phase x velocity/wavelength)
with a fixed page-sized figure.

Features
--------
- Reads linear axis metadata from FITS header.
- Uses a fixed physical figure size (Letter/A4).
- Uses imshow(aspect="auto"), so the plot height stays fixed.
- Makes NaN pixels visible (black by default) so gaps are obvious.
- Forces the color bar lower limit to start at 0 by default.
- Optional hard clipping of values below vmin using Normalize(..., clip=True).
- For phase plots, displays a default range of -0.5..1.5 cycles.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter1d


def axis_from_header(hdr: fits.Header, axis: int, n: int) -> np.ndarray:
    """
    Reconstruct a linear axis from FITS-style CRVAL/CDELT/CRPIX.
    """
    crval = float(hdr.get(f"CRVAL{axis}"))
    cdelt = float(hdr.get(f"CDELT{axis}"))
    crpix = float(hdr.get(f"CRPIX{axis}", 1.0))
    pix_fits = np.arange(n, dtype=float) + 1.0
    return crval + (pix_fits - crpix) * cdelt



def label_from_ctype(ctype: str, cunit: str, axis: int) -> str:
    """
    Make a reasonable axis label from CTYPEn/CUNITn.
    """
    ctype_u = (ctype or "").strip().upper()
    unit = (cunit or "").strip()

    if axis == 1:
        if "VELO" in ctype_u:
            return "RV (km/s)" if unit else "Velocity"
        if "WAVE" in ctype_u or "LAMB" in ctype_u:
            return f"Wavelength ({unit})" if unit else "Wavelength"
        return f"X ({ctype})" if ctype else "X"

    if "MJD" in ctype_u:
        return "MJD"
    if "PHASE" in ctype_u:
        return "Phase"
    return f"Y ({ctype})" if ctype else "Y"


def tile_phase_for_plot(img: np.ndarray, y: np.ndarray, lo: float = -0.5, hi: float = 1.5) -> Tuple[np.ndarray, np.ndarray, bool]:
    """
    Tile phase rows by -1/0/+1 so the default view can cover [-0.5, 1.5].
    """
    if y.size == 0:
        return img, y, False

    chunks_img = []
    chunks_y = []

    for shift in (-1.0, 0.0, 1.0):
        ys = y + shift
        ys_lo = float(np.min(ys))
        ys_hi = float(np.max(ys))
        if ys_hi < lo or ys_lo > hi:
            continue
        chunks_img.append(img)
        chunks_y.append(ys)

    if not chunks_img:
        return img, y, False

    return np.vstack(chunks_img), np.concatenate(chunks_y), True


def gauss_filter_x(img: np.ndarray, fwhm_pix: float) -> np.ndarray:
    """NaN-aware Gaussian filter along the x-axis (wavelength direction, axis=1)."""
    sigma = fwhm_pix / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    nan_mask = ~np.isfinite(img)
    filled = np.where(nan_mask, 0.0, img)
    weight = np.where(nan_mask, 0.0, 1.0)
    smooth = gaussian_filter1d(filled, sigma, axis=1)
    wsmooth = gaussian_filter1d(weight, sigma, axis=1)
    result = np.where(wsmooth > 0, smooth / wsmooth, np.nan)
    result[nan_mask] = np.nan
    return result


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Plot a dynamical spectrum FITS with fixed page size."
    )
    ap.add_argument("input_fits", type=Path, help="Input FITS (2D image).")
    ap.add_argument("output", type=Path, help="Output plot (.pdf, .png, etc.).")

    # Figure size / layout controls
    ap.add_argument("--figsize", type=float, nargs=2, metavar=("WIDTH", "HEIGHT"),
                    default=[5.0, 4.0],
                    help="Figure width and height in inches. Default=5 4.")
    ap.add_argument("--dpi", type=int, default=200,
                    help="DPI for raster output (PNG). Default=200.")
    ap.add_argument("--fontsize", type=float, default=7.0,
                    help="Base font size in points (scales title, labels, ticks). Default=7.")

    # Color/display controls
    ap.add_argument("--cmap", default="viridis",
                    help="Matplotlib colormap. Default=viridis.")
    ap.add_argument("--vmin", type=float, default=None,
                    help="Color scale minimum. Default: max(0, VCUTLO or 0).")
    ap.add_argument("--vmax", type=float, default=None,
                    help="Color scale maximum. Default: image maximum (finite pixels).")
    ap.add_argument("--use-header-cuts", action="store_true",
                    help="Use header VCUTLO/VCUTHI defaults instead of full data range.")
    ap.add_argument("--clip-below-zero", action="store_true",
                    help="Hard-clip values below vmin to the minimum color.")
    ap.add_argument("--nan-color", default="black",
                    help="Color used for NaN pixels. Default=black.")
    ap.add_argument("--nan-alpha", type=float, default=1.0,
                    help="Alpha for NaN pixels (0=transparent). Default=1.")

    # Gaussian smoothing in wavelength direction
    ap.add_argument("--gauss-fwhm", type=float, default=None,
                    help="FWHM of Gaussian smoothing kernel along wavelength axis, in Angstroms.")

    # Phase display range
    ap.add_argument("--phase-lo", type=float, default=-0.5,
                    help="Lower phase limit for phase plots. Default=-0.5.")
    ap.add_argument("--phase-hi", type=float, default=1.5,
                    help="Upper phase limit for phase plots. Default=1.5.")

    # Label visibility
    ap.add_argument("--hide-xlabel", action="store_true", help="Hide x-axis label.")
    ap.add_argument("--hide-ylabel", action="store_true", help="Hide y-axis label.")
    ap.add_argument("--hide-cbar-label", action="store_true", help="Hide colorbar label.")

    # Diagnostics / cosmetics
    ap.add_argument("--title", default=None,
                    help="Optional plot title. Default: input FITS filename.")
    ap.add_argument("--info", action="store_true",
                    help="Print a quick summary before plotting.")

    args = ap.parse_args(argv)

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    plt.rcParams["font.size"] = args.fontsize

    with fits.open(args.input_fits, memmap=False) as hdul:
        img = np.asarray(hdul[0].data, dtype=float)
        hdr = hdul[0].header

    if img.ndim != 2:
        raise SystemExit(f"Expected a 2D image in the primary HDU; got shape {img.shape}")

    ny, nx = img.shape
    x = axis_from_header(hdr, axis=1, n=nx)
    y = axis_from_header(hdr, axis=2, n=ny)
    is_phase = "PHASE" in str(hdr.get("CTYPE2", "")).strip().upper()
    phase_default_view = False
    if is_phase:
        img, y, phase_default_view = tile_phase_for_plot(img, y, lo=args.phase_lo, hi=args.phase_hi)
        ny, nx = img.shape

    xlab = label_from_ctype(hdr.get("CTYPE1", ""), hdr.get("CUNIT1", ""), axis=1)
    ylab = label_from_ctype(hdr.get("CTYPE2", ""), hdr.get("CUNIT2", ""), axis=2)

    if args.gauss_fwhm is not None:
        cdelt1 = abs(float(hdr.get("CDELT1", 1.0)))
        xmode = str(hdr.get("XMODE", hdr.get("CTYPE1", ""))).strip().upper()
        if "VELO" in xmode or xmode == "VELOCITY":
            lamc = float(hdr.get("LAMC", 0.0))
            if lamc <= 0.0:
                raise SystemExit("--gauss-fwhm requires LAMC in the FITS header for velocity-axis files")
            C_KMS = 299792.458
            fwhm_pix = (args.gauss_fwhm / lamc) * C_KMS / cdelt1
        else:
            fwhm_pix = args.gauss_fwhm / cdelt1
        img = gauss_filter_x(img, fwhm_pix)

    if args.info:
        if args.gauss_fwhm is not None:
            print(f"Gaussian smoothing: FWHM={args.gauss_fwhm} AA ({fwhm_pix:.2f} pix)")
        finite = np.isfinite(img)
        finite_rows = int(np.sum(np.any(finite, axis=1)))
        finite_frac = float(np.mean(finite))
        print(f"Image shape: {img.shape} (Y rows x X cols)")
        print(f"Finite rows: {finite_rows}/{ny}  (finite fraction: {finite_frac:.4f})")
        print(f"X range: {x[0]:.6g} .. {x[-1]:.6g}   (Dx={hdr.get('CDELT1')})")
        print(f"Y range: {y[0]:.6g} .. {y[-1]:.6g}   (Dy={hdr.get('CDELT2')})")
        if "N_SPEC" in hdr:
            print(f"N_SPEC (spectra used): {hdr['N_SPEC']}")
        if "VCUTLO" in hdr and "VCUTHI" in hdr:
            print(f"Header display cuts: VCUTLO={hdr['VCUTLO']}, VCUTHI={hdr['VCUTHI']}")

    # Choose color scale
    finite = np.isfinite(img)
    finite_vals = img[finite] if np.any(finite) else np.array([0.0], dtype=float)

    if args.vmin is not None:
        vmin = max(0.0, float(args.vmin))
    elif args.use_header_cuts and hdr.get("VCUTLO") is not None:
        vmin = max(0.0, float(hdr.get("VCUTLO")))
    else:
        vmin = max(0.0, float(np.nanmin(finite_vals)))

    if args.vmax is not None:
        vmax = float(args.vmax)
    elif args.use_header_cuts and hdr.get("VCUTHI") is not None:
        vmax = float(hdr.get("VCUTHI"))
    else:
        vmax = float(np.nanmax(finite_vals))

    if vmax <= vmin:
        vmax = vmin + 1e-6

    if args.info:
        print(f"Display scale: vmin={vmin:.6g}, vmax={vmax:.6g}")

    w_in, h_in = args.figsize
    fig, ax = plt.subplots(figsize=(w_in, h_in), constrained_layout=True)

    # Make NaNs visible
    cmap = plt.get_cmap(args.cmap).copy()
    cmap.set_bad(color=args.nan_color, alpha=float(args.nan_alpha))

    extent = [float(x[0]), float(x[-1]), float(y[0]), float(y[-1])]

    if args.clip_below_zero:
        # Hard clipping below vmin (typically 0)
        norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        im = ax.imshow(
            img,
            origin="lower",
            aspect="auto",
            extent=extent,
            interpolation="nearest",
            cmap=cmap,
            norm=norm,
        )
    else:
        im = ax.imshow(
            img,
            origin="lower",
            aspect="auto",
            extent=extent,
            interpolation="nearest",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )

    if not args.hide_xlabel:
        ax.set_xlabel(xlab)
    if not args.hide_ylabel:
        ax.set_ylabel(ylab)
    if args.title:
        ax.set_title(args.title, fontsize=7)
    if phase_default_view:
        ax.set_ylim(args.phase_lo, args.phase_hi)

    # Inward ticks on all four sides (no secondary axes, which confuse constrained_layout)
    ax.minorticks_on()
    ax.tick_params(which="both", top=True, right=True, direction="in", color="white")

    cbar = fig.colorbar(im, ax=ax, location="top")
    if not args.hide_cbar_label:
        cbar.set_label("Flux")
    cbar_ticks = np.linspace(vmin, vmax, 3)
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([f"{t:.3f}" for t in cbar_ticks])
    cbar.ax.tick_params(which="both", direction="in")

    # Do not use bbox_inches="tight" — it conflicts with constrained_layout
    fig.savefig(args.output, dpi=args.dpi)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
