#!/usr/bin/env python3
"""Collect per-case summary CSVs into a single table."""

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


SUMMARY_FILE_SUFFIX = "_summary.csv"
DEFAULT_REPORT_PATH = Path(__file__).parent / "runs" / "summary.md"

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
    "result_status",
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
    "result_status": "Computation status for this case: computed, zero_output, or error (see run logs).",
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
    "result_status",
)

PREVIEW_ROW_LIMIT = 12

MAPDL_ROOT = Path(__file__).resolve().parent
if str(MAPDL_ROOT) not in sys.path:
    sys.path.insert(0, str(MAPDL_ROOT))

import report_best_runs as best_report


ERROR_MARKERS = (
    "*** ERROR ***",
    "*** FATAL ERROR ***",
    "***FATAL ERROR***",
    "*** SEVERE ERROR ***",
    "SOLUTION ABORTED",
)
STATUS_NUMERIC_FIELDS: Sequence[str] = (
    "tip_defl_m",
    "buckling_factor",
    "sigma_max_pa",
    "buckling_pair_N",
    "buckling_per_mass",
)


def _safe_float(raw: Optional[str]) -> Optional[float]:
    if raw is None:
        return None
    text = raw.strip()
    if not text or set(text) == {"*"}:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _has_meaningful_value(raw: Optional[str], tol: float = 1e-12) -> bool:
    value = _safe_float(raw)
    if value is None:
        return False
    return abs(value) > tol


def _infer_result_status(row: Dict[str, str], run_dir: Path) -> str:
    error_text: List[str] = []
    for name in ("file0.err", "file1.err", "file0.out", "file1.out"):
        path = run_dir / name
        if not path.exists():
            continue
        try:
            error_text.append(path.read_text(encoding="utf-8", errors="ignore"))
        except Exception:
            continue
    blob = "\n".join(error_text)
    if any(marker in blob for marker in ERROR_MARKERS):
        return "error"
    if any(_has_meaningful_value(row.get(field)) for field in STATUS_NUMERIC_FIELDS):
        return "computed"
    return "zero_output"


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
            row["result_status"] = row.get("result_status", "")
            row["result_status"] = _infer_result_status(row, subdir)
            _merge_fieldnames(fieldnames, ("result_status",))
            _merge_fieldnames(fieldnames, row.keys())
            rows.append(row)
    
    return rows, fieldnames


def _read_consolidated_summary(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    if not path.exists():
        return [], []
    with path.open(newline="") as fp:
        data_lines = [line for line in fp if not line.lstrip().startswith("#")]
    if not data_lines:
        return [], []
    reader = csv.DictReader(data_lines)
    fieldnames = list(reader.fieldnames or [])
    rows: List[Dict[str, str]] = []
    for raw in reader:
        row: Dict[str, str] = {}
        for key, value in raw.items():
            if key is None:
                continue
            row[key] = (value or "").strip()
        rows.append(row)
    return fieldnames, rows


def _normalise_rows_for_comparison(
    rows: Sequence[Dict[str, str]],
    columns: Sequence[str],
) -> List[Tuple[str, ...]]:
    normalised: List[Tuple[str, ...]] = []
    for row in rows:
        normalised.append(tuple((row.get(column, "") or "").strip() for column in columns))
    normalised.sort()
    return normalised


def _compute_digest(rows: Sequence[Dict[str, str]], fieldnames: Sequence[str]) -> str:
    canonical = sorted(set(fieldnames) | {"case_id"})
    serialisable = [
        [(row.get(column, "") or "").strip() for column in canonical]
        for row in rows
    ]
    serialisable.sort()
    payload = json.dumps(serialisable, ensure_ascii=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _summaries_match(
    new_rows: Sequence[Dict[str, str]],
    new_fields: Sequence[str],
    existing_rows: Sequence[Dict[str, str]],
    existing_fields: Sequence[str],
) -> bool:
    if not existing_rows and new_rows:
        return False
    columns = sorted(set(new_fields) | set(existing_fields) | {"case_id"})
    return _normalise_rows_for_comparison(new_rows, columns) == _normalise_rows_for_comparison(existing_rows, columns)



def _render_full_markdown_table(
    rows: Sequence[Dict[str, str]],
    fieldnames: Sequence[str],
) -> str:
    if not rows or not fieldnames:
        return "_No data available to render._"

    usable_columns = [name for name in fieldnames if any((row.get(name) or "").strip() for row in rows)]
    if not usable_columns:
        return "_Column set produced no displayable data._"

    def _alignment_marker(width: int, align: str) -> str:
        width = max(width, 3)
        if align == "right":
            return "-" * (width - 1) + ":"
        if align == "center":
            if width == 3:
                return ":-:"
            return ":" + "-" * (width - 2) + ":"
        return ":" + "-" * (width - 1)

    def _format_numeric(value: float) -> str:
        if value == 0:
            return "0"
        abs_value = abs(value)
        if abs_value >= 1e3 or abs_value < 1e-3:
            return f"{value:.3e}"
        if abs_value >= 1:
            text = f"{value:.4f}"
        else:
            text = f"{value:.4g}"
        return text.rstrip("0").rstrip(".")

    column_is_numeric: Dict[str, bool] = {name: True for name in usable_columns}
    for row in rows:
        for name in usable_columns:
            raw = (row.get(name, "") or "").strip()
            if not raw or set(raw) == {"*"}:
                continue
            try:
                float(raw)
            except ValueError:
                column_is_numeric[name] = False

    formatted_rows: List[List[str]] = []
    widths = {name: len(name) for name in usable_columns}
    for row in rows:
        formatted_row = []
        for name in usable_columns:
            raw = (row.get(name, "") or "").strip()
            if not raw:
                value = "-"
            elif set(raw) == {"*"}:
                value = raw
            elif column_is_numeric[name]:
                try:
                    value = _format_numeric(float(raw))
                except ValueError:
                    value = raw
            else:
                value = raw
            widths[name] = max(widths[name], len(value))
            formatted_row.append(value)
        formatted_rows.append(formatted_row)

    header_parts: List[str] = []
    divider_parts: List[str] = []
    for name in usable_columns:
        align = "right" if column_is_numeric[name] else "left"
        header_parts.append(name.rjust(widths[name]) if align == "right" else name.ljust(widths[name]))
        divider_parts.append(_alignment_marker(widths[name], align))
    header = "| " + " | ".join(header_parts) + " |"
    divider = "| " + " | ".join(divider_parts) + " |"
    lines = [header, divider]
    for formatted_row in formatted_rows:
        cells = []
        for idx, name in enumerate(usable_columns):
            align = "right" if column_is_numeric[name] else "left"
            value = formatted_row[idx]
            cells.append(value.rjust(widths[name]) if align == "right" else value.ljust(widths[name]))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def _build_combined_markdown(
    summary_path: Path,
    runs_dir: Path,
    digest: str,
    rows: Sequence[Dict[str, str]],
    fieldnames: Sequence[str],
    ranked: Sequence[best_report.Row],
    weights: Dict[str, float],
    stats: Dict[str, best_report.MetricStats],
    top_n: int,
) -> str:
    total_rows = len(rows)
    ranked_count = len(ranked)
    selected = list(enumerate(ranked, start=1))[:top_n]

    lines: List[str] = []
    lines.append("# MAPDL Run Summary\n")
    lines.append(f"- Source summary: `{summary_path}`")
    lines.append(f"- Runs directory: `{runs_dir}`")
    lines.append(f"- Cases aggregated: {total_rows}")
    lines.append(f"- Ranked rows: {ranked_count}")
    lines.append(f"- Digest: `{digest}`")
    weight_text = ", ".join(f"{name}={weights[name]:.2f}" for name in sorted(weights))
    lines.append(f"- Weights: {weight_text}\n")

    lines.append("## Top Results\n")
    if selected:
        lines.append(best_report.render_table(selected))
    else:
        lines.append("No rows contained enough data to compute a score.")

    lines.append("\n## Metric Highlights\n")
    lines.append(best_report.render_metric_highlights(ranked, stats))
    lines.append("\n## Metric Ranges\n")
    lines.append(best_report.render_stats_section(stats))
    lines.append("\n## Full Dataset\n")
    lines.append(_render_full_markdown_table(rows, fieldnames))
    lines.append("\n")

    return "\n".join(lines).rstrip() + "\n"

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





def _write_full_markdown_table(
    path: Path,
    rows: Sequence[Dict[str, str]],
    fieldnames: Sequence[str],
) -> None:
    content = "# Consolidated MAPDL run summary (all rows)\n\n" + _render_full_markdown_table(rows, fieldnames) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _render_comment_block(
    rows: Sequence[Dict[str, str]],
    fieldnames: Sequence[str],
    digest: str,
    report_path: Optional[Path],
    top_n: int,
) -> str:
    lines = []
    total_rows = len(rows)
    lines.append("# Consolidated MAPDL run summary")
    lines.append("# Lines starting with '#' are comments for human readers; CSV parsers should ignore them.")
    lines.append(f"# Rows reported below: {total_rows}")
    lines.append("# Values filled with '*' come directly from MAPDL and typically indicate missing data.")
    lines.append(f"# Digest: {digest}")
    lines.append("# Use `python summarize_runs.py --check-only` to verify this summary before relying on cached outputs.")
    lines.append("#")
    if fieldnames:
        name_width = max(len(name) for name in fieldnames)
        lines.append("# Column glossary:")
        for name in fieldnames:
            description = COLUMN_DESCRIPTIONS.get(name, "Additional MAPDL output column.")
            lines.append(f"# - {name.ljust(name_width)} : {description}")
        lines.append("#")
    if report_path is not None:
        lines.append(f"# Combined Markdown summary (top {top_n} + full table) written to: {report_path}")
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
    ap.add_argument(
        "--report",
        "--markdown",
        dest="report",
        type=Path,
        default=DEFAULT_REPORT_PATH,
        help="Path for the combined Markdown summary (use '-' to skip).",
    )
    ap.add_argument(
        "--top",
        type=int,
        default=10,
        help="Number of high-scoring rows to include at the top of the Markdown summary.",
    )
    ap.add_argument(
        "--weight",
        action="append",
        default=[],
        metavar="METRIC=VALUE",
        help="Override default metric weights, e.g. --weight tip_per_load=0.5.",
    )
    ap.add_argument(
        "--force",
        action="store_true",
        help="Rebuild outputs even if the existing summary matches the current runs.",
    )
    ap.add_argument(
        "--check-only",
        action="store_true",
        help="Only report whether the consolidated summary is up to date without writing files.",
    )
    args = ap.parse_args()

    if not args.runs_dir.exists():
        raise SystemExit(f"Runs directory '{args.runs_dir}' does not exist")

    rows, fieldnames = gather_rows(args.runs_dir)
    if not rows:
        raise SystemExit("No summary CSVs found. Run MAPDL cases first")

    digest = _compute_digest(rows, fieldnames)

    if args.out.exists():
        existing_fieldnames, existing_rows = _read_consolidated_summary(args.out)
        summary_matches = _summaries_match(rows, fieldnames, existing_rows, existing_fieldnames)
    else:
        summary_matches = False

    if args.check_only:
        if summary_matches:
            print(f"Summary at {args.out} matches the current run data ({len(rows)} rows).")
            return
        raise SystemExit(f"Summary at {args.out} is missing or stale. Re-run without --check-only to rebuild.")

    report_path: Optional[Path]
    if str(args.report) == "-":
        report_path = None
    else:
        report_path = args.report

    weights = best_report.parse_weight_overrides(args.weight)

    needs_write = args.force or not summary_matches
    if not needs_write and report_path is not None and not report_path.exists():
        needs_write = True

    if not needs_write:
        print(f"Summary already up to date at {args.out}")
        if report_path is not None and report_path.exists():
            print(f"Combined Markdown summary already up to date at {report_path}")
        return

    args.out.parent.mkdir(parents=True, exist_ok=True)
    comment_block = _render_comment_block(rows, fieldnames, digest, report_path, args.top)
    with args.out.open("w", newline="") as f:
        f.write(comment_block)
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})

    print(f"Wrote consolidated summary with {len(rows)} rows to {args.out}")

    if report_path is not None:
        ranking_rows = best_report.rows_from_dicts(rows)
        stats = best_report.compute_metric_stats(ranking_rows)
        best_report.assign_scores(ranking_rows, weights, stats)
        ranked = best_report.sorted_rows(ranking_rows)
        combined_markdown = _build_combined_markdown(
            summary_path=args.out,
            runs_dir=args.runs_dir,
            digest=digest,
            rows=rows,
            fieldnames=fieldnames,
            ranked=ranked,
            weights=weights,
            stats=stats,
            top_n=args.top,
        )
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(combined_markdown, encoding="utf-8")
        print(f"Wrote Markdown summary with {min(args.top, len(ranked))} ranked rows to {report_path}")


if __name__ == "__main__":  # pragma: no cover
    main()
