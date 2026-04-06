# Workflow Batch Subdirectory

This subdirectory is a separate workflow layer and does not modify the main working scripts.

Files:

- `run_dynspec_batch.py`
- `plot_dynspec_batch.py`

## 1) Batch generation

Edit configuration inside `run_dynspec_batch.py`:

- `LIST_FILE`
- `LINE_NAMES` (read from `spectral_line_catalog.py`)
- `MODES`
- `ORBITAL_SOLUTIONS` (for phase modes)
- sampling/window parameters (`LAMSTEP_A`, `VSTEP_KMS`, etc.)

Run:

```bash
python workflow_batch/run_dynspec_batch.py
```

Outputs go to:

`workflow_batch/outputs/`

Output FITS names include mode + list + line + region + period/T0 (for phase) so files do not overwrite each other.

## 2) Batch plotting

Edit configuration inside `plot_dynspec_batch.py`:

- `INPUT_DIR`, `PLOT_DIR`
- `INPUT_FILES` (optional explicit list)
- plot style options (`DPI`, `CMAP`, etc.)

Run:

```bash
python workflow_batch/plot_dynspec_batch.py
```

Plots go to:

`workflow_batch/plots/`
