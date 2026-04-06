# Dynspec Python Suite (Current Working Version)

This directory contains the working Python tools for building and plotting dynamical spectra:

- `ts2ima.py` (time x wavelength image)
- `ts2vima.py` (time x velocity image)
- `ts2phima.py` (phase x wavelength image)
- `ts2vphima.py` (phase x velocity image)
- `ts2tab.py` (wavelength table output)
- `ts2vtab.py` (velocity table output)
- `ts2cov.py` (velocity covariance-like matrix)
- `plot_dynspec.py` (plot FITS image outputs)

Naming note: `ts2` means `time series to`. Common suffixes are `ima` = image, `vima` = velocity image, `phima` = phase image, `vphima` = velocity-phase image, `tab` = table, `vtab` = velocity table, and `cov` = covariance-like matrix.

Core image-building logic is in `ts2dynspec.py`.

## Input List Format

`--list` must point to a text file containing one FITS path per line.

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

## Example Input Commands (Image Modes)

Use these as paste-ready examples.

### 1) `ts2ima` (time x wavelength)

```bash
python ts2ima.py \
  --out ex_ts2ima.fits \
  --list spectra_UVES.list \
  --lamlo 6200 --lamhi 6400 --lamstep 1.0 \
  --tstep 0.5
```

### 2) `ts2vima` (time x velocity)

```bash
python ts2vima.py \
  --out ex_ts2vima.fits \
  --list spectra_UVES.list \
  --lamc 6319 --vlo -400 --vhi 400 --vstep 2.0 \
  --tstep 0.5
```

### 3) `ts2phima` (phase x wavelength)

```bash
python ts2phima.py \
  --out ex_ts2phima.fits \
  --list spectra_UVES.list \
  --lamlo 6200 --lamhi 6400 --lamstep 1.0 \
  --period 87.49 --t0-mjd 59988.29 \
  --phstep 0.05
```

### 4) `ts2vphima` (phase x velocity, overlap supported)

```bash
python ts2vphima.py \
  --out ex_ts2vphima.fits \
  --list spectra_UVES.list \
  --lamc 6319 --vlo -400 --vhi 400 --vstep 2.0 \
  --period 87.49 --t0-mjd 59988.29 \
  --phstep 0.05 \
  --overlap 10
```

## Example Input Commands (Table / Covariance Modes)

### 5) `ts2tab` (wavelength table)

```bash
python ts2tab.py \
  --out ex_ts2tab.fits \
  --list spectra_UVES.list \
  --lamlo 6200 --lamhi 6400 --lamstep 1.0
```

### 6) `ts2vtab` (velocity table)

```bash
python ts2vtab.py \
  --out ex_ts2vtab.fits \
  --list spectra_UVES.list \
  --lamc 6319 --vlo -400 --vhi 400 --vstep 2.0
```

### 7) `ts2cov` (velocity covariance-like matrix)

```bash
python ts2cov.py \
  --out ex_ts2cov.fits \
  --list spectra_UVES.list \
  --lamc 6319 --vlo -400 --vhi 400 --vstep 2.0
```

#### Understanding `ts2cov` output

`ts2cov` writes a single 2D FITS image (Primary HDU) with shape `N_v x N_v`, where `N_v` is the number of velocity bins from `--vlo` to `--vhi` sampled by `--vstep`.

- `CTYPE1 = VELO`, `CUNIT1 = km/s`
- `CTYPE2 = VELO`, `CUNIT2 = km/s`

So each pixel `(i, j)` compares variability at velocity `v_i` with variability at velocity `v_j`.

Processing steps:

1. Each input spectrum is rebinned to the common velocity grid.
2. Optional preprocessing is applied in time order:
   - edge renormalization (`--renorm`)
   - 3-point median filter along time (`--filter`)
3. A covariance matrix across time is computed.
4. A legacy TS2COV transform is applied to the covariance matrix entries.

If `C_ij` is covariance and `C_ii`, `C_jj` are variances:

- Off-diagonal (`i != j`):
  - `M_ij = C_ij / sqrt(C_ii + C_jj)`
- Diagonal (`i == j`):
  - `M_ii = sqrt(C_ii / 2)`

Interpretation:

- Positive off-diagonal values: the two velocity bins tend to vary together.
- Negative off-diagonal values: anticorrelated variability.
- Near-zero off-diagonal values: weak shared variability.
- Diagonal values: variability amplitude proxy at each velocity bin (after TS2COV transform).

Notes:

- This is a covariance-like diagnostic, not a Pearson correlation matrix.
- The matrix is symmetric by construction.
- Input list order must be strictly ascending in time for `ts2cov`.

Quick plotting example for `ts2cov` (diverging colormap centered at 0):

```python
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

in_fits = "ex_ts2cov.fits"
out_png = "ex_ts2cov.png"

with fits.open(in_fits, memmap=False) as hdul:
    img = np.asarray(hdul[0].data, dtype=float)
    h = hdul[0].header

nx = img.shape[1]
ny = img.shape[0]
x = h["CRVAL1"] + (np.arange(nx) + 1 - h.get("CRPIX1", 1.0)) * h["CDELT1"]
y = h["CRVAL2"] + (np.arange(ny) + 1 - h.get("CRPIX2", 1.0)) * h["CDELT2"]
extent = [float(x[0]), float(x[-1]), float(y[0]), float(y[-1])]

finite = np.isfinite(img)
absmax = np.nanpercentile(np.abs(img[finite]), 99) if np.any(finite) else 1.0
norm = TwoSlopeNorm(vmin=-absmax, vcenter=0.0, vmax=absmax)

plt.figure(figsize=(8, 6))
im = plt.imshow(img, origin="lower", aspect="auto", extent=extent, cmap="RdBu_r", norm=norm)
plt.xlabel("Velocity (km/s)")
plt.ylabel("Velocity (km/s)")
plt.title("TS2COV matrix")
plt.colorbar(im, label="Covariance-like value")
plt.tight_layout()
plt.savefig(out_png, dpi=200)
plt.close()
```

## Plot Commands

```bash
python plot_dynspec.py ex_ts2ima.fits ex_ts2ima.png
python plot_dynspec.py ex_ts2vima.fits ex_ts2vima.png
python plot_dynspec.py ex_ts2phima.fits ex_ts2phima.png
python plot_dynspec.py ex_ts2vphima.fits ex_ts2vphima.png
```

For phase images, plotting defaults to a displayed phase range of `-0.5 .. 1.5`.

## Behavior Notes (Current Version)

- Phase modes use C-style phase-bin averaging.
- `ts2vphima` supports `--overlap` as percent of one phase bin.
- If `--lamstep` is omitted in wavelength image modes, default is inferred from the first usable spectrum.
- Time interpolation in image modes uses:
  - 1 point: no interpolation
  - 2 points: linear interpolation
  - more than 2 points: cubic spline

## Practical Tip

For fast smoke tests, use coarse sampling (`--lamstep 1.0`, `--vstep 2.0` or `5.0`) and narrow windows.
