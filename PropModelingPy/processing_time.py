"""Measure the processing time of BET, BEMT and Vortex."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
from pathlib import Path
from time import perf_counter

import numpy as np

from apc_data import (
    discover_apc_cases,
    load_geometry,
    load_performance,
    select_representative_rpm,
)
from bemt import run_bemt
from bet import run_bet
from polar_lookup import PolarDatabase
from vortex import run_vortex


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
POLAR = DATA / "ClarkY" / "polares_clarky_expan.mat"
MODELS = (("BET", run_bet), ("BEMT", run_bemt), ("Vortex", run_vortex))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    """Write a table without requiring pandas."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def measure_case(case: object) -> list[dict[str, object]]:
    """Load one propeller and time only the three model calls."""
    geometry = load_geometry(case.geometry_path)
    blocks = load_performance(case.performance_path)
    rpm, _, _ = select_representative_rpm(blocks, geometry.diameter_in)
    reference = blocks[rpm]
    polar = PolarDatabase(POLAR, mode="mature")
    rows = []
    for name, model in MODELS:
        start = perf_counter()
        result = model(geometry, reference.j, rpm, polar)
        elapsed = perf_counter() - start
        rows.append(
            {
                "propeller": case.name,
                "rpm": rpm,
                "model": name,
                "seconds": elapsed,
                "j_points": reference.j.size,
                "radial_stations": geometry.x.size,
                "finite_ct_points": int(np.count_nonzero(np.isfinite(result.ct))),
            }
        )
    return rows


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Calculate mean, median, spread and the cost relative to BET."""
    summary = []
    for name, _ in MODELS:
        values = np.array([row["seconds"] for row in rows if row["model"] == name])
        summary.append(
            {
                "model": name,
                "propellers": values.size,
                "mean_seconds": float(np.mean(values)),
                "median_seconds": float(np.median(values)),
                "std_seconds": float(np.std(values, ddof=1)) if values.size > 1 else 0.0,
                "min_seconds": float(np.min(values)),
                "max_seconds": float(np.max(values)),
            }
        )
    bet_median = summary[0]["median_seconds"]
    for row in summary:
        row["ratio_to_bet"] = float(row["median_seconds"] / bet_median)
    return summary


def run(limit: int | None, jobs: int, output_dir: Path, include_marine: bool) -> None:
    """Measure all confirmed Clark-Y propellers or a limited sample."""
    cases = [
        case
        for case in discover_apc_cases(DATA, include_marine=include_marine)
        if case.status == "included"
    ]
    if limit is not None:
        cases = cases[:limit]

    rows: list[dict[str, object]] = []
    if jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            iterator = executor.map(measure_case, cases)
            for index, case_rows in enumerate(iterator, 1):
                rows.extend(case_rows)
                print(f"[{index}/{len(cases)}] benchmark complete", flush=True)
    else:
        for index, case in enumerate(cases, 1):
            rows.extend(measure_case(case))
            print(f"[{index}/{len(cases)}] benchmark complete", flush=True)

    summary = summarize(rows)
    write_csv(output_dir / "processing_time_by_propeller.csv", rows)
    write_csv(output_dir / "processing_time_summary.csv", summary)
    print("\nModel      median [s]   ratio to BET")
    for row in summary:
        print(f"{row['model']:7s} {row['median_seconds']:12.4f} {row['ratio_to_bet']:14.2f}x")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--include-marine", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "Results" / "ProcessingTime")
    args = parser.parse_args()
    run(args.limit, max(1, args.jobs), args.output_dir, args.include_marine)


if __name__ == "__main__":
    main()
