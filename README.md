# Dynspec Python Suite

This suite builds and plots dynamical spectra from a time series of 1D FITS spectra.

## Getting Started

Edit the inputs at the top of `run_batch.py` and run it. That's the main script — it generates FITS files and PNG plots in one go.

Outputs:
- FITS files → `outputs/`
- PNG plots → `plots/`

Set `SAVE_FITS = False` to delete FITS files after plotting and keep only the PNGs.

## Scripts

- `run_batch.py` — main batch runner: generates FITS and plots them
- `ts2dynspec.py` — core image-building library
- `ts2ima.py` — time × wavelength image
- `ts2vima.py` — time × velocity image
- `ts2phima.py` — phase × wavelength image
- `ts2vphima.py` — phase × velocity image
- `ts2tab.py` — wavelength table output
- `ts2vtab.py` — velocity table output
- `ts2cov.py` — velocity covariance-like matrix
- `plot_dynspec.py` — plot a single FITS image

Naming: `ts2` = "time series to". Suffixes: `ima` = image, `vima` = velocity image, `phima` = phase image, `vphima` = velocity-phase image, `tab` = table, `vtab` = velocity table, `cov` = covariance matrix.

## Input List Format

`LIST_FILE` must point to a text file with one FITS path per line (paths relative to the list file location, or absolute).

Example:
```text
spectra/ANCol_UVES_580U_60219_2032.fits
spectra/ANCol_UVES_580U_60225_2518.fits
spectra/ANCol_UVES_580U_60229_3508.fits
```

## Quick Validation

```bash
python ts2dynspec.py --validate-only --list spectra_UVES.list
```

## Behavior Notes

- Phase modes use C-style phase-bin averaging.
- `ts2vphima` supports `OVERLAP_PERCENT` as percent of one phase bin.
- If `--lamstep` is omitted in wavelength image modes, it is inferred from the first usable spectrum.
- Time interpolation in image modes:
  - 1 point: no interpolation
  - 2 points: linear interpolation
  - more than 2 points: cubic spline

## Understanding `ts2cov` output

`ts2cov` writes a 2D FITS image with shape `N_v × N_v`, where each pixel `(i, j)` compares variability at velocity `v_i` with velocity `v_j`.

Processing: spectra rebinned to velocity grid → optional renorm/filter → covariance matrix → legacy TS2COV transform:

- Off-diagonal: `M_ij = C_ij / sqrt(C_ii + C_jj)`
- Diagonal: `M_ii = sqrt(C_ii / 2)`

The matrix is symmetric. Positive off-diagonal = correlated variability; negative = anticorrelated. Input list must be strictly ascending in time.

Quick plotting snippet for `ts2cov`:

```python
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

with fits.open("ex_ts2cov.fits", memmap=False) as hdul:
    img = np.asarray(hdul[0].data, dtype=float)
    h = hdul[0].header

nx = img.shape[1]
x = h["CRVAL1"] + (np.arange(nx) + 1 - h.get("CRPIX1", 1.0)) * h["CDELT1"]
extent = [float(x[0]), float(x[-1]), float(x[0]), float(x[-1])]

absmax = np.nanpercentile(np.abs(img[np.isfinite(img)]), 99)
norm = TwoSlopeNorm(vmin=-absmax, vcenter=0.0, vmax=absmax)

plt.figure(figsize=(8, 6))
im = plt.imshow(img, origin="lower", aspect="auto", extent=extent, cmap="RdBu_r", norm=norm)
plt.xlabel("Velocity (km/s)")
plt.ylabel("Velocity (km/s)")
plt.colorbar(im, label="Covariance-like value")
plt.tight_layout()
plt.savefig("ex_ts2cov.png", dpi=200)
plt.close()
```
