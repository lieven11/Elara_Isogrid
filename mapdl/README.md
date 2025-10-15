MAPDL Automation

Overview
- This folder contains a small scaffold to run parametric MAPDL (APDL) jobs for isogrid shells, separate from the heatmap workflow.
- Parameters align with existing heatmap code: `R, L, a, b, t` (and `n_theta` as an alternative to `a`).
- Outputs are written under `mapdl/runs/<case_id>/` to keep results separate.

Contents
- `templates/isogrid_0deg_tri_template.inp`: APDL template with placeholders for geometry and material parameters.
- `cases/`: Organised case libraries:
  - `samples/`: quick single cases (baseline per material).
  - `sweeps/`: 1D sweeps for thickness, length, and cell pitch covering requested ranges (`t = 0.5…4 mm`, `L = 0.25…2.0 m`, `a = 0.05…0.60 m`).
  - `grids/`: coarse 3×3×3 design grid so you can launch multi-parameter batches.
  - `generated/`: optional output location when you build custom libraries with the generator.
- `specs/`: Sweep templates (e.g. `tri0_dense_template.json`) describing ranges/steps.
- `geometry_preview/`: Static Matplotlib previews for lattice geometry sanity checks (`geometry_preview/tri0_preview.py`).
- `case_generator.py`: Expand a spec into a case-library JSON (thousands of combos if needed).
- `summarize_runs.py`: Collate per-case summary CSVs into one `summary.csv` after MAPDL finishes.
- `runner.py`: Expand cases into APDL input files and optionally run MAPDL if available.
- `utils.py`: Helpers to pull materials and derive basic geometry values.

Usage
1) Pick a case-library JSON (e.g. `mapdl/cases/sweeps/thickness.json`) **or** generate your own from a spec:
   - Thickness sweep (inputs only, macOS/Linux): `python3 mapdl/runner.py --cases mapdl/cases/sweeps/thickness.json --dry-run`
   - Thickness sweep (run MAPDL, Windows with Ansys): `python mapdl/runner.py --cases mapdl/cases/sweeps/thickness.json`
   - Length sweep (inputs only, macOS/Linux): `python3 mapdl/runner.py --cases mapdl/cases/sweeps/length.json --dry-run`
   - Length sweep (run MAPDL, Windows with Ansys): `python mapdl/runner.py --cases mapdl/cases/sweeps/length.json`
   - Cell-pitch sweep (inputs only, macOS/Linux): `python3 mapdl/runner.py --cases mapdl/cases/sweeps/cellsize.json --dry-run`
   - Cell-pitch sweep (run MAPDL, Windows with Ansys): `python mapdl/runner.py --cases mapdl/cases/sweeps/cellsize.json`
   - Dense multi-parameter spec → library: `python3 mapdl/case_generator.py --spec mapdl/specs/tri0_dense_template.json --out mapdl/cases/generated/tri0_dense.json`
2) Generate APDL input files (no solve) for the chosen library (example shown for the dense library):
   `python3 mapdl/runner.py --cases mapdl/cases/generated/tri0_dense.json --dry-run`
   (use `python` instead of `python3` on Windows if that’s how Python is installed)
3) Run MAPDL when available (drop `--dry-run`, using the same library as above):
   `python mapdl/runner.py --cases mapdl/cases/generated/tri0_dense.json`
4) After runs finish, gather all per-case metrics into a single table:
   `python3 mapdl/summarize_runs.py`

Command examples
- Dense 3-parameter library (no solve): `python mapdl/runner.py --cases mapdl/cases/generated/tri0_dense.json --dry-run`
- Dense 3-parameter library (solve if MAPDL available): `python mapdl/runner.py --cases mapdl/cases/generated/tri0_dense.json`
- Rib thickness sweep only: `python mapdl/runner.py --cases mapdl/cases/sweeps/thickness.json [--dry-run]`
- Segment length sweep only: `python mapdl/runner.py --cases mapdl/cases/sweeps/length.json [--dry-run]`
- Cell pitch sweep only: `python mapdl/runner.py --cases mapdl/cases/sweeps/cellsize.json [--dry-run]`
- Coarse 3-3-3 design grid: `python mapdl/runner.py --cases mapdl/cases/grids/tri0_design_grid.json [--dry-run]`
- Geometry preview (interactive window): `py -3.12 mapdl/geometry_preview/tri0_preview.py`
- Geometry preview to PNG (headless): `py -3.12 mapdl/geometry_preview/tri0_preview.py --save preview.png`

Post-release workflow
1. Pull the latest code and review the changelog or git log for MAPDL-related updates.
2. Regenerate any derived case libraries if specs changed (for example: `python3 mapdl/case_generator.py --spec ... --out ...`).
3. Rebuild APDL input decks without solving to confirm nothing fails: `python3 mapdl/runner.py --cases <library>.json --dry-run`.
4. Launch the full solves on the target machine (drop `--dry-run` and prefer the Windows launcher if MAPDL is local).
5. Refresh the aggregated outputs once solves finish: `python3 mapdl/summarize_runs.py`.
6. Inspect `mapdl/runs/<case_id>/` for spot checks and open `mapdl/runs/summary.csv` (or your merged file) to compare against the previous release.

Updating case libraries
- Re-run the generator with the same `--out` path to overwrite a library that should pick up new logic: `python3 mapdl/case_generator.py --spec mapdl/specs/tri0_dense_template.json --out mapdl/cases/generated/tri0_dense.json`.
- Add `--dry-run` first if you want to confirm the case count before writing.
- To refresh every spec in one go on PowerShell:
  ```powershell
  Get-ChildItem mapdl/specs -Filter *.json | ForEach-Object {
    $spec = $_.FullName
    $out = "mapdl/cases/generated/$($_.BaseName).json"
    py -3.12 mapdl/case_generator.py --spec $spec --out $out
  }
  ```
- Adjust the launcher (`python3`, `python`, or `py -3.12`) to match your environment.

Python 3.12 (Windows launcher)
- If you manage Python via the Windows launcher, replace the `python`/`python3` calls above with `py -3.12` (for example: `py -3.12 mapdl/runner.py --cases mapdl/cases/sweeps/thickness.json --dry-run`).

Sometimes it doesnt work: run
1. $Env:ANSYSLMD_LICENSE_FILE = '1055@WIN-ELARA-GPU'
2. echo $Env:ANSYSLMD_LICENSE_FILE
1055@WIN-ELARA-GPU
3. py -3.12 -c "from ansys.mapdl.core import launch_mapdl; mapdl = launch_mapdl(); print(mapdl.version); mapdl.exit()"
>> 
24.2
4. py -3.12 mapdl/runner.py --cases mapdl/cases/sweeps/thickness.json 
(für 4.: oder anderer command für andere Fälle)

Notes
- The APDL template is a starting point for a 0°-family triangular isogrid using BEAM188. It parameterizes lattice pitch `a`, rib width `b`, thickness `t`, length `L`, and radius `R`.
- If `a` is not provided, it is computed from `R` and `n_theta` (cells around the circumference).
- Material options mirror the dashboard/heatmap scripts: AISI 304, AISI 316L, Al 6061-T6, Al 2024-T3, Al 7075-T6 (see `isogrid/materials.py`).
- `case_generator.py` respects range/step definitions in the spec. Tightening step sizes (e.g. 0.01 m on `L`, 0.002 m on `a`, 0.00025 m on `t`) quickly scales to tens of thousands of permutations—check the console summary before writing.
- Each MAPDL run writes `<case_id>_summary.csv` inside its run folder with the key outputs (tip deflection, linear buckling factor, etc.), so `summarize_runs.py` can stitch everything into one overview CSV for post-processing.
- Heatmap outputs remain under the existing `out/` folder for now; MAPDL outputs are under `mapdl/runs/`. If you want, we can later consolidate heatmap outputs under `out_heatmap/` and update references in your scripts.
- On macOS you typically only run the generator/runner in `--dry-run` mode (MAPDL isn’t available); execute the same command without `--dry-run` on a Windows/Linux machine with Ansys installed.
