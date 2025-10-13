#!/usr/bin/env python3
"""Collect per-case summary CSVs into a single table."""

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


SUMMARY_FILE_SUFFIX = "_summary.csv"


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

    rows: List[Dict[str, str]] = []
    for raw in data_rows:
        if not any(cell.strip() for cell in raw):
            continue
        values = [cell.strip() for cell in raw]
        if len(values) < len(header):
            values.extend([""] * (len(header) - len(values)))
        elif len(values) > len(header):
            values = values[: len(header)]
        row = {header[i]: values[i] for i in range(len(header))}
        rows.append(row)

    return header, rows


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
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})

    print(f"Wrote consolidated summary with {len(rows)} rows to {args.out}")


if __name__ == "__main__":  # pragma: no cover
    main()

