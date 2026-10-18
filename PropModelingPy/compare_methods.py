"""Compare BET, BEMT and Vortex with one APC propeller."""

from __future__ import annotations

import argparse
from pathlib import Path

from apc_data import load_geometry, load_performance
from bemt import run_bemt
from bet import run_bet
from metrics import curve_metrics
from plots import plot_individual
from polar_lookup import PolarDatabase
from vortex import run_vortex


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
POLAR = DATA / "ClarkY" / "polares_clarky_expan.mat"


def performance_file(propeller: str) -> Path:
    """Find the APC performance file in the aircraft or marine folder."""
    base = DATA / "APC" / "performance"
    matches = sorted(base.rglob(f"PER3_{propeller}.dat"))
    if not matches:
        raise FileNotFoundError(f"APC performance file not found for {propeller}")
    return matches[0]


def compare(propeller: str, requested_rpm: int, output: Path) -> None:
    """Run the same geometry, operating point and polar in all three models."""
    geometry_path = DATA / "APC" / "geometry" / f"{propeller}-PERF.PE0"
    geometry = load_geometry(geometry_path)
    blocks = load_performance(performance_file(propeller))
    rpm = min(blocks, key=lambda value: abs(value - requested_rpm))
    reference = blocks[rpm]
    polar = PolarDatabase(POLAR, mode="mature")

    # Each method receives the same J values, geometry, RPM and Clark-Y polar.
    results = {
        "BET": run_bet(geometry, reference.j, rpm, polar),
        "BEMT": run_bemt(geometry, reference.j, rpm, polar),
        "Vortex": run_vortex(geometry, reference.j, rpm, polar),
    }

    print(f"APC {propeller} at {rpm} rpm")
    print("Model     sMAPE CT   sMAPE CP   sMAPE eta    r2 CT    r2 CP   r2 eta")
    for name, result in results.items():
        ct = curve_metrics(reference.ct, result.ct)
        cp = curve_metrics(reference.cp, result.cp)
        eta = curve_metrics(reference.eta, result.eta)
        print(
            f"{name:7s} {ct.smape:10.3f} {cp.smape:10.3f} {eta.smape:11.3f}"
            f" {ct.pearson_r2:8.4f} {cp.pearson_r2:8.4f} {eta.pearson_r2:8.4f}"
        )
    plot_individual(propeller, rpm, reference, results, output)
    print(f"Figure saved to {output.resolve()}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--propeller", default="22x10")
    parser.add_argument("--rpm", type=int, default=5000)
    parser.add_argument("--output", type=Path, default=ROOT / "Results" / "comparison.png")
    args = parser.parse_args()
    compare(args.propeller, args.rpm, args.output)


if __name__ == "__main__":
    main()
