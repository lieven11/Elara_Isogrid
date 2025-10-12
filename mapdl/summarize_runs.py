#!/usr/bin/env python3
"""Collect per-case summary CSVs into a single table."""

import argparse
import csv
from pathlib import Path
from typing import Dict, List


SUMMARY_FILE_SUFFIX = "_summary.csv"


def gather_rows(runs_dir: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for subdir in sorted(runs_dir.iterdir()):
        if not subdir.is_dir():
            continue
        case_id = subdir.name
        summary_path = subdir / f"{case_id}{SUMMARY_FILE_SUFFIX}"
        if not summary_path.exists():
            continue
        with summary_path.open() as f:
            reader = csv.DictReader(f)
            for row in reader:
                rows.append(row)
    return rows


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

    rows = gather_rows(args.runs_dir)
    if not rows:
        raise SystemExit("No summary CSVs found. Run MAPDL cases first")

    fieldnames = list(rows[0].keys())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote consolidated summary with {len(rows)} rows to {args.out}")


if __name__ == "__main__":  # pragma: no cover
    main()

