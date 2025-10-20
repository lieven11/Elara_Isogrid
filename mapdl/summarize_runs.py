#!/usr/bin/env python3
"""Collect per-case summary CSVs into a single table."""

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


SUMMARY_FILE_SUFFIX = "_summary.csv"

CANONICAL_HEADER: Sequence[str] = (
    "case_id",
    "material",
    "mat_code",
    "R_m",
    "L_m",
    "a_m",
    "b_m",
    "t_m",
    "n_theta",
    "tip_defl_m",
    "buckling_factor",
    "sigma_max_pa",
    "axial_load_face_N",
    "face_pressure_pa",
    "load_pair_total_N",
    "buckling_face_N",
    "buckling_pair_N",
    "total_mass_kg",
    "buckling_per_mass",
    "tip_per_load_m_per_N",
    "sigma_per_load_pa_per_N",
)
LEGACY_HEADER_ORDER: Sequence[str] = CANONICAL_HEADER[:-2]
LEGACY_ADDED_COLUMNS: Sequence[str] = CANONICAL_HEADER[-2:]
LEGACY_MARKERS = {"load_pai", "axial_lo"}

COLUMN_DESCRIPTIONS = {
    "case_id": "Unique case identifier combining lattice family, material, and geometry inputs.",
    "material": "Human-readable material name pulled from the MAPDL material library.",
    "mat_code": "Short material code (e.g. 304, 316L) used by MAPDL input files.",
    "R_m": "Shell mid-surface radius in metres.",
    "L_m": "Shell axial length in metres.",
    "a_m": "Isogrid cell pitch along the circumferential direction in metres.",
    "b_m": "Rib (strut) width in metres.",
    "t_m": "Rib (strut) thickness in metres.",
    "n_theta": "Number of lattice bays around the circumference (blank if pitch a_m was supplied).",
    "tip_defl_m": "Tip deflection under the reference load case in metres.",
    "buckling_factor": "Linear buckling eigenvalue relative to the reference axial load (1.0 means equal to reference load).",
    "sigma_max_pa": "Maximum von Mises stress reported by MAPDL in pascals.",
    "axial_load_face_N": "Axial force applied to a single face sheet in newtons.",
    "face_pressure_pa": "Distributed pressure applied to the face sheet in pascals.",
    "load_pair_total_N": "Total axial load applied to the opposing faces combined in newtons.",
    "buckling_face_N": "Estimated linear buckling capacity for a single face in newtons.",
    "buckling_pair_N": "Estimated linear buckling capacity for the face pair in newtons.",
    "total_mass_kg": "Total mass of the analysed shell segment in kilograms.",
    "buckling_per_mass": "Buckling capacity per unit mass (buckling_pair_N / total_mass_kg) in N/kg.",
    "tip_per_load_m_per_N": "Tip compliance per unit load in metres per newton.",
    "sigma_per_load_pa_per_N": "Stress amplification per unit load in pascals per newton.",
}

PREVIEW_COLUMNS: Sequence[str] = (
    "case_id",
    "material",
    "R_m",
    "L_m",
    "t_m",
    "a_m",
    "b_m",
    "tip_defl_m",
    "buckling_factor",
    "sigma_max_pa",
    "load_pair_total_N",
    "buckling_pair_N",
    "total_mass_kg",
    "buckling_per_mass",
)
PREVIEW_ROW_LIMIT = 12


def _read_summary_file(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    """Return (fieldnames, rows) parsed from a per-case summary CSV.

    Some MAPDL summary exports include a trailing comma in the header, which
    pushes the final column name onto its own line. We normalise that here so
    downstream consumers see a single, clean header row.
    """
    with path.open(newline="") as f:
        raw_rows = list(csv.reader(f))

    if not raw_rows:
        return [], []

    header = [cell.strip() for cell in raw_rows[0] if cell.strip()]
    data_rows = raw_rows[1:]

    # Fold in header fragments that end up on their own line (single non-empty cell)
    while data_rows:
        first = data_rows[0]
        non_empty = [cell.strip() for cell in first if cell.strip()]
        if len(non_empty) == 1 and all(not cell.strip() for cell in first[1:]):
            header.append(non_empty[0])
            data_rows.pop(0)
        else:
            break

    legacy_format = _is_legacy_header(header)
    rows: List[Dict[str, str]] = []
    for raw in data_rows:
        if not any(cell.strip() for cell in raw):
            continue
        values = [cell.strip() for cell in raw]
        if legacy_format:
            row = _build_legacy_row(values)
            rows.append(row)
            continue
        if len(values) < len(header):
            values.extend([""] * (len(header) - len(values)))
        elif len(values) > len(header):
            values = values[: len(header)]
        row = {header[i]: values[i] for i in range(len(header))}
        rows.append(row)

    if legacy_format:
        return list(CANONICAL_HEADER), rows
    return header, rows


def _is_legacy_header(header: Sequence[str]) -> bool:
    if not header:
        return False
    if any(marker in header for marker in LEGACY_MARKERS):
        return True
    if "tip_defl_m" not in header and "tip_defl" in header:
        return True
    if header.count("buckling") > 1:
        return True
    return False


def _build_legacy_row(values: Sequence[str]) -> Dict[str, str]:
    normalised: Dict[str, str] = {}
    trimmed = list(values)[: len(LEGACY_HEADER_ORDER)]
    if len(trimmed) < len(LEGACY_HEADER_ORDER):
        trimmed.extend([""] * (len(LEGACY_HEADER_ORDER) - len(trimmed)))
    for name, value in zip(LEGACY_HEADER_ORDER, trimmed):
        normalised[name] = value
    for name in LEGACY_ADDED_COLUMNS:
        normalised[name] = ""
    return normalised


def _merge_fieldnames(existing: List[str], incoming: Iterable[str]) -> None:
    for name in incoming:
        if name not in existing:
            existing.append(name)


def gather_rows(runs_dir: Path) -> Tuple[List[Dict[str, str]], List[str]]:
    rows: List[Dict[str, str]] = []
    fieldnames: List[str] = []
    for subdir in sorted(runs_dir.iterdir()):
        if not subdir.is_dir():
            continue
        case_id = subdir.name
        summary_path = subdir / f"{case_id}{SUMMARY_FILE_SUFFIX}"
        if not summary_path.exists():
            continue
        header, parsed_rows = _read_summary_file(summary_path)
        if not parsed_rows:
            continue
        if not fieldnames:
            fieldnames = list(header)
        else:
            _merge_fieldnames(fieldnames, header)
        for row in parsed_rows:
            _merge_fieldnames(fieldnames, row.keys())
            rows.append(row)
    return rows, fieldnames


def _format_markdown_table(rows: Sequence[Dict[str, str]], columns: Sequence[str], limit: int) -> Sequence[str]:
    preview = list(rows[:limit])
    usable_columns = [name for name in columns if any(row.get(name) for row in rows)]
    if not preview or not usable_columns:
        return ()
    widths = {name: len(name) for name in usable_columns}
    for row in preview:
        for name in usable_columns:
            value = (row.get(name) or "").strip() or "N/A"
            widths[name] = max(widths[name], len(value))
    header = "| " + " | ".join(name.ljust(widths[name]) for name in usable_columns) + " |"
    divider = "| " + " | ".join("-" * widths[name] for name in usable_columns) + " |"
    lines = [header, divider]
    for row in preview:
        line = "| " + " | ".join(
            ((row.get(name) or "").strip() or "N/A").ljust(widths[name]) for name in usable_columns
        ) + " |"
        lines.append(line)
    return lines


def _render_comment_block(rows: Sequence[Dict[str, str]], fieldnames: Sequence[str]) -> str:
    lines = []
    total_rows = len(rows)
    lines.append("# Consolidated MAPDL run summary")
    lines.append("# Lines starting with '#' are comments for human readers; CSV parsers should ignore them.")
    lines.append(f"# Rows reported below: {total_rows}")
    lines.append("# Values filled with '*' come directly from MAPDL and typically indicate missing data.")
    lines.append("#")
    if fieldnames:
        name_width = max(len(name) for name in fieldnames)
        lines.append("# Column glossary:")
        for name in fieldnames:
            description = COLUMN_DESCRIPTIONS.get(name, "Additional MAPDL output column.")
            lines.append(f"# - {name.ljust(name_width)} : {description}")
        lines.append("#")
    preview_lines = _format_markdown_table(rows, PREVIEW_COLUMNS, PREVIEW_ROW_LIMIT)
    if preview_lines:
        lines.append("# Preview table (first rows):")
        for line in preview_lines:
            lines.append(f"# {line}")
        if total_rows > PREVIEW_ROW_LIMIT:
            lines.append(f"# ... truncated after {PREVIEW_ROW_LIMIT} rows in the preview (full data below in CSV).")
        lines.append("#")
    else:
        lines.append("# No rows available to preview.")
        lines.append("#")
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--runs-dir",
        type=Path,
        default=Path(__file__).parent / "runs",
        help="Directory containing per-case run folders",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).parent / "runs" / "summary.csv",
        help="Output CSV path",
    )
    args = ap.parse_args()

    if not args.runs_dir.exists():
        raise SystemExit(f"Runs directory '{args.runs_dir}' does not exist")

    rows, fieldnames = gather_rows(args.runs_dir)
    if not rows:
        raise SystemExit("No summary CSVs found. Run MAPDL cases first")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        comment_block = _render_comment_block(rows, fieldnames)
        f.write(comment_block)
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})

    print(f"Wrote consolidated summary with {len(rows)} rows to {args.out}")


if __name__ == "__main__":  # pragma: no cover
    main()
