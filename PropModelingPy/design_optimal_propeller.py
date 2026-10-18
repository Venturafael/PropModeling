"""Reproduce the TCC Minimum Energy Loss design point from the command line."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, replace
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from polar_lookup import PolarDatabase
from apc_data import load_geometry
from mel import MELInputs, design_sweep, solve_power


ROOT = Path(__file__).resolve().parent
DATA_ROOT = ROOT / "Data"
POLAR_PATH = DATA_ROOT / "ClarkY" / "polares_clarky_expan.mat"
RESULTS_ROOT = ROOT / "Results"


def require_packaged_data() -> None:
    """Stop early when the local APC or polar data are missing."""
    if not (DATA_ROOT / "APC" / "geometry").is_dir():
        raise FileNotFoundError(DATA_ROOT / "APC" / "geometry")
    if not POLAR_PATH.is_file():
        raise FileNotFoundError(POLAR_PATH)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        if fields:
            writer.writeheader()
            writer.writerows(rows)


def _jsonable(value: object) -> object:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def radial_rows(state: dict[str, object]) -> list[dict[str, object]]:
    """Flatten radial arrays into an auditable station table."""

    fields = (
        "xi",
        "c",
        "beta",
        "alpha",
        "alpha_target",
        "re",
        "mach",
        "f",
        "g",
        "a",
        "ap",
        "dct_dxi",
        "dcp_dxi",
    )
    return [
        {field: float(np.asarray(state[field])[index]) for field in fields}
        for index in range(np.asarray(state["xi"]).size)
    ]


def plot_design(state: dict[str, object], output: Path) -> None:
    """Plot geometry, aerodynamic state and radial loading."""

    xi = np.asarray(state["xi"])
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4), constrained_layout=True)
    axes[0].plot(xi, 1000.0 * np.asarray(state["c"]), lw=2.2, color="#1f77b4")
    twin = axes[0].twinx()
    twin.plot(xi, state["beta"], lw=2.0, color="#d95f02")
    axes[0].set_ylabel("Chord [mm]", color="#1f77b4")
    twin.set_ylabel(r"Blade angle, $\beta$ [deg]", color="#d95f02")
    axes[0].set_title("Generated geometry")

    axes[1].plot(xi, np.asarray(state["re"]) / 1e3, lw=2.2, color="#1f77b4")
    twin_re = axes[1].twinx()
    twin_re.plot(xi, state["mach"], lw=2.0, color="#d95f02")
    axes[1].set_ylabel("Reynolds number [thousands]", color="#1f77b4")
    twin_re.set_ylabel("Local Mach", color="#d95f02")
    axes[1].set_title("Section operating state")

    axes[2].plot(xi, state["dct_dxi"], lw=2.2, label=r"$dC_T/d(r/R)$")
    axes[2].plot(xi, state["dcp_dxi"], lw=2.2, label=r"$dC_P/d(r/R)$")
    axes[2].set_ylabel("Coefficient density")
    axes[2].set_title("Radial loading")
    axes[2].legend(frameon=False)
    for axis in axes:
        axis.set_xlabel("Normalized radius, r/R")
        axis.grid(True, alpha=0.22)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Adkins–Liebeck Minimum Energy Loss design")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=190, bbox_inches="tight")
    plt.close(fig)


def run(
    *,
    cl_target: float,
    alpha_bias_deg: float,
    stations: int,
    sweep: bool,
    output_dir: Path,
) -> dict[str, object]:
    """Run the reproducibility case and write all numerical intermediates."""

    require_packaged_data()
    if stations < 4:
        raise ValueError("stations must be at least four")
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    reference_geometry = load_geometry(
        DATA_ROOT / "APC" / "geometry" / "22x10-PERF.PE0"
    )
    xi_root = max(0.10, float(np.min(reference_geometry.x)))
    xi = np.unique(
        np.clip(
            np.r_[np.linspace(xi_root, 0.95, stations - 1), 0.99],
            xi_root,
            0.99,
        )
    )
    inputs = replace(MELInputs(), alpha_bias_deg=alpha_bias_deg)
    # Legacy-linear interpolation is intentional: it validates the selected
    # MATLAB implementation. Analysis models use the mature common lookup.
    polar = PolarDatabase(POLAR_PATH, mode="legacy_linear")
    if sweep:
        sweep_records, state = design_sweep(inputs, polar, xi)
        _write_csv(output_dir / "mel_cl_sweep.csv", sweep_records)
    else:
        state = solve_power(cl_target, xi, inputs, polar)

    _write_csv(output_dir / "mel_radial_geometry_and_state.csv", radial_rows(state))
    _write_csv(
        output_dir / "mel_power_closure_history.csv",
        list(state["iteration_history"]),
    )
    plot_design(state, output_dir / "mel_design.png")
    summary: dict[str, object] = {
        "inputs": asdict(inputs),
        "selected_cl": state["cl_target"],
        "zeta": state["zeta"],
        "advance_ratio": inputs.advance_ratio,
        "ct": state["ct"],
        "cp": state["cp"],
        "eta": state["eta"],
        "thrust_n": state["thrust_n"],
        "power_w": state["power_w"],
        "torque_nm": state["torque_nm"],
        "activity_factor": state["activity_factor"],
        "polar_file": str(POLAR_PATH),
        "polar_mode": "legacy_linear",
        "alpha_bias_caveat": (
            "The -3.5 deg TCC default changes the alpha used for Cd while chord "
            "and epsilon retain target Cl. Preserved for reproduction; validate "
            "aerodynamic consistency before claiming a physical optimum."
        ),
    }
    (output_dir / "mel_summary.json").write_text(
        json.dumps(summary, indent=2, default=_jsonable), encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cl-target", type=float, default=0.76)
    parser.add_argument(
        "--alpha-bias-deg",
        type=float,
        default=-3.5,
        help="Preserved MATLAB/TCC default; use 0 only as an explicit sensitivity case.",
    )
    parser.add_argument("--stations", type=int, default=15)
    parser.add_argument("--sweep", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=RESULTS_ROOT / "mel")
    args = parser.parse_args()
    summary = run(
        cl_target=args.cl_target,
        alpha_bias_deg=args.alpha_bias_deg,
        stations=args.stations,
        sweep=args.sweep,
        output_dir=args.output_dir,
    )
    print(json.dumps(summary, indent=2, default=_jsonable))


if __name__ == "__main__":
    main()
