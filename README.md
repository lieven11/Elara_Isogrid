# Elara Isogrid Workspace

This repository now separates the interactive Python tooling from the MAPDL batch
workflows so it is clearer where each set of assets lives.

- `heatmap/` – Plotly dashboards, sweep scripts, generated HTML outputs (`heatmap/out`),
  and archived dashboard variants (`heatmap/archive`). The package is importable as
  `heatmap.scripts.*`.
- `mapdl/` – Input templates, case generators, run utilities, and generated MAPDL
  results (`mapdl/runs`).
- `isogrid/` – Shared analytical helpers (geometry, materials, mechanics) used by both
  the heatmap tooling and the MAPDL automation.

Each script that writes files ensures the target directory exists, so you can still
run them from any working directory. The default output destinations were updated to
reflect the new layout.

