# MAPDL Run Summary

- Source file: `c:\Users\lieve\OneDrive\Dokumente\GitHub\Elara_Isogrid\mapdl\runs\summary.csv`
- Total rows parsed: 40
- Rows with complete score: 40
- Weights: buckling=0.35, sigma_ma=0.30, tip_defl=0.35

## Overview

| Metric | Unit | Reported Range | Best Case (Case -> Material) | Score Weight | Notes |
| --- | --- | --- | --- | --- | --- |
| Tip Deflection | m | 0.0000 - 0.0000 | tri0_304 -> AISI 304 | 0.35 | Displacement column `tip_defl` is 0.0000 for every row in `summary.csv`. |
| Buckling Factor | - | 0.0000 - 0.0000 | tri0_304 -> AISI 304 | 0.35 | Multiple `buckling` columns exist in the CSV; the first is entirely 0.0000. |
| Von Mises Stress | Pa | 0.0000 - 0.0000 | tri0_304 -> AISI 304 | 0.30 | `sigma_ma` is also 0.0000 throughout the CSV. |

## Top Results

| Rank | Case | Material | t [m] | a [m] | b [m] | Tip Defl. | Buckling | sigma_ma | Score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | tri0_304 | AISI 304 | 5.0000e-04 | 0.3500 | 0.0120 | 0.0000 | 0.0000 | 0.0000 | 1.0000 |
| 2 | tri0_304 | AISI 304 | 0.0010 | 0.3500 | 0.0120 | 0.0000 | 0.0000 | 0.0000 | 1.0000 |
| 3 | tri0_304 | AISI 304 | 0.0015 | 0.3500 | 0.0120 | 0.0000 | 0.0000 | 0.0000 | 1.0000 |
| 4 | tri0_304 | AISI 304 | 0.0020 | 0.3500 | 0.0120 | 0.0000 | 0.0000 | 0.0000 | 1.0000 |
| 5 | tri0_304 | AISI 304 | 0.0025 | 0.3500 | 0.0120 | 0.0000 | 0.0000 | 0.0000 | 1.0000 |
| 6 | tri0_304 | AISI 304 | 0.0030 | 0.3500 | 0.0120 | 0.0000 | 0.0000 | 0.0000 | 1.0000 |
| 7 | tri0_304 | AISI 304 | 0.0035 | 0.3500 | 0.0120 | 0.0000 | 0.0000 | 0.0000 | 1.0000 |
| 8 | tri0_304 | AISI 304 | 0.0040 | 0.3500 | 0.0120 | 0.0000 | 0.0000 | 0.0000 | 1.0000 |
| 9 | tri0_316 | AISI 316 | 5.0000e-04 | 0.3500 | 0.0120 | 0.0000 | 0.0000 | 0.0000 | 1.0000 |
| 10 | tri0_316 | AISI 316 | 0.0010 | 0.3500 | 0.0120 | 0.0000 | 0.0000 | 0.0000 | 1.0000 |

## Metric Highlights

| Metric | Top Case -> Material | Value | Unit | Comment |
| --- | --- | --- | --- | --- |
| Tip Deflection | tri0_304 -> AISI 304 | 0.0000 | m | Identical zero deflection reported for every run. |
| Buckling Factor | tri0_304 -> AISI 304 | 0.0000 | - | Buckling output appears unpopulated in the CSV. |
| Von Mises Stress | tri0_304 -> AISI 304 | 0.0000 | Pa | Stress results are zero across the dataset. |

## Metric Ranges

- **Tip Deflection [m]**: min 0.0000, max 0.0000
- **Buckling Factor**: min 0.0000, max 0.0000
- **Von Mises Stress [Pa]**: min 0.0000, max 0.0000

## Data Quality Notes

- Every metric column used for scoring is zero in `mapdl/runs/summary.csv`, which means the solver results were not captured or parsed correctly.
- The raw CSV header includes duplicated `buckling` fields, which can break CSV readers and may indicate a post-processing script bug.
- Confirm that the MAPDL post-processing step exports non-zero values (e.g., update the result extraction script or check for unit scaling issues) before re-running the summary generation.
