"""Run the auditable global APC Clark--Y comparison."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import json
from pathlib import Path
import traceback

import numpy as np
from scipy.stats import rankdata, spearmanr

from polar_lookup import PolarDatabase
from apc_data import discover_apc_cases, load_geometry, load_performance, select_representative_rpm
from metrics import average_metric_rank, curve_metrics
from bemt import run_bemt
from bet import run_bet
from vortex import run_vortex
from plots import plot_aggregate, plot_individual


ROOT = Path(__file__).resolve().parent
DATA_ROOT = ROOT / "Data"
POLAR_PATH = DATA_ROOT / "ClarkY" / "polares_clarky_expan.mat"
RESULTS_ROOT = ROOT / "Results"
REPOSITORY_ROOT = ROOT


def require_packaged_data() -> None:
    """Stop early when the local APC or polar data are missing."""
    if not (DATA_ROOT / "APC" / "geometry").is_dir():
        raise FileNotFoundError(DATA_ROOT / "APC" / "geometry")
    if not POLAR_PATH.is_file():
        raise FileNotFoundError(POLAR_PATH)


_WORKER_POLAR: PolarDatabase | None = None


def display_path(path: Path) -> str:
    """Prefer repository-relative paths without breaking external datasets."""

    try:
        return str(path.relative_to(REPOSITORY_ROOT))
    except ValueError:
        return str(path)


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        if fieldnames:
            writer.writeheader()
            writer.writerows(rows)


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def inventory_cases(cases: list, max_tip_mach: float) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in cases:
        row: dict[str, object] = {
            "propeller": case.name,
            "geometry_file": display_path(case.geometry_path),
            "performance_file": display_path(case.performance_path) if case.performance_path else "",
            "airfoil_metadata": " | ".join(case.airfoils),
            "status": case.status,
            "reason": case.reason,
        }
        try:
            geometry = load_geometry(case.geometry_path)
            row.update({
                "diameter_in": geometry.diameter_in,
                "nominal_pitch_in": geometry.nominal_pitch_in,
                "blades": geometry.blades,
                "radial_stations": geometry.x.size,
            })
            if case.performance_path:
                blocks = load_performance(case.performance_path)
                rpms = sorted(blocks)
                chosen, tip_mach, selection = select_representative_rpm(
                    blocks, geometry.diameter_in, max_tip_mach=max_tip_mach
                )
                all_j = np.concatenate([block.j for block in blocks.values()])
                row.update({
                    "rpm_available_count": len(rpms), "rpm_min": min(rpms), "rpm_max": max(rpms),
                    "rpm_available": ";".join(map(str, rpms)), "selected_rpm": chosen,
                    "selected_static_tip_mach": tip_mach, "rpm_selection_rule": selection,
                    "j_min_all_rpms": float(np.min(all_j)), "j_max_all_rpms": float(np.max(all_j)),
                })
        except Exception as exc:
            row["inventory_error"] = f"{type(exc).__name__}: {exc}"
        rows.append(row)
    return rows


def metrics_for_case(propeller: str, rpm: int, reference, results: dict) -> list[dict[str, object]]:
    rows = []
    for model, result in results.items():
        met_ct = curve_metrics(reference.ct, result.ct)
        met_cp = curve_metrics(reference.cp, result.cp)
        met_eta = curve_metrics(reference.eta, result.eta)
        rows.append({
            "propeller": propeller, "rpm": rpm, "model": model,
            "n_ct": met_ct.n, "n_cp": met_cp.n, "n_eta": met_eta.n,
            "smape_ct": met_ct.smape, "smape_cp": met_cp.smape, "smape_eta": met_eta.smape,
            "r2_ct": met_ct.pearson_r2, "r2_cp": met_cp.pearson_r2, "r2_eta": met_eta.pearson_r2,
        })
    return rows


def _initialize_worker() -> None:
    global _WORKER_POLAR
    _WORKER_POLAR = PolarDatabase(POLAR_PATH, mode="mature")


def _process_case(payload: tuple[object, float, bool, Path]) -> tuple[list, list, list, dict]:
    case, max_tip_mach, skip_plots, output_dir = payload
    global _WORKER_POLAR
    polar_mature = _WORKER_POLAR or PolarDatabase(POLAR_PATH, mode="mature")
    geometry = load_geometry(case.geometry_path)
    blocks = load_performance(case.performance_path)
    rpm, static_tip_mach, selection_rule = select_representative_rpm(
        blocks, geometry.diameter_in, max_tip_mach=max_tip_mach
    )
    reference = blocks[rpm]
    results = {}
    issue_rows: list[dict[str, object]] = []
    for model, runner in (("BET", run_bet), ("BEMT", run_bemt), ("Vortex", run_vortex)):
        try:
            result = runner(geometry, reference.j, rpm, polar_mature)
            results[model] = result
            for issue in result.issues:
                issue_rows.append({"propeller": case.name, "rpm": rpm, **issue})
            stats = result.polar_stats
            clamp_total = sum(value for key, value in stats.items() if key != "queries")
            if clamp_total:
                issue_rows.append({
                    "propeller": case.name, "rpm": rpm, "model": model, "j": "ALL",
                    "category": "polar_domain", "reason": json.dumps(stats, sort_keys=True),
                    "treatment": "coordinate clipped to nearest database boundary; no extrapolation",
                })
        except Exception as exc:
            issue_rows.append({
                "propeller": case.name, "rpm": rpm, "model": model, "j": "ALL",
                "category": "fatal_numerical", "reason": f"{type(exc).__name__}: {exc}",
                "treatment": "model omitted for this propeller", "traceback": traceback.format_exc(),
            })

    metric_rows = metrics_for_case(case.name, rpm, reference, results)
    model_for_regime = next(iter(results.values()), None)
    regime_row = {
        "propeller": case.name, "diameter_in": geometry.diameter_in,
        "nominal_pitch_in": geometry.nominal_pitch_in, "blades": geometry.blades,
        "rpm": rpm, "rpm_selection_rule": selection_rule, "static_tip_mach": static_tip_mach,
        "tip_mach": float(np.nanmedian(model_for_regime.mach_tip)) if model_for_regime else np.nan,
        "representative_reynolds": float(np.nanmedian(model_for_regime.reynolds_representative)) if model_for_regime else np.nan,
        "apc_reynolds_75pct_median": float(np.nanmedian(reference.reynolds_apc)),
        "apc_mach_median": float(np.nanmedian(reference.mach_apc)),
        "j_min": float(np.min(reference.j)), "j_max": float(np.max(reference.j)), "j_mean": float(np.mean(reference.j)),
        "ct_reference_mean": float(np.mean(reference.ct)), "cp_reference_mean": float(np.mean(reference.cp)),
    }
    point_rows: list[dict[str, object]] = []
    for point_index, j in enumerate(reference.j):
        point_rows.append({
            "propeller": case.name, "rpm": rpm, "model": "APC reference", "j": j,
            "ct": reference.ct[point_index], "cp": reference.cp[point_index], "eta": reference.eta[point_index],
        })
        for model, result in results.items():
            point_rows.append({
                "propeller": case.name, "rpm": rpm, "model": model, "j": j,
                "ct": result.ct[point_index], "cp": result.cp[point_index], "eta": result.eta[point_index],
            })
    if not skip_plots:
        plot_individual(
            case.name,
            rpm,
            reference,
            results,
            output_dir / "individual" / f"{case.name}.png",
        )
    return metric_rows, point_rows, issue_rows, regime_row


def aggregate(metric_rows: list[dict[str, object]], regime_rows: list[dict[str, object]]) -> tuple[list, list, list, list, list]:
    models = ["BET", "BEMT", "Vortex"]
    outputs = ["ct", "cp", "eta"]
    summary: list[dict[str, object]] = []
    for model in models:
        subset = [row for row in metric_rows if row["model"] == model]
        for output in outputs:
            for metric in ("smape", "r2"):
                values = np.array([safe_float(row[f"{metric}_{output}"]) for row in subset])
                values = values[np.isfinite(values)]
                summary.append({
                    "model": model, "output": output, "metric": metric, "n": values.size,
                    "mean": float(np.mean(values)) if values.size else np.nan,
                    "median": float(np.median(values)) if values.size else np.nan,
                    "std": float(np.std(values, ddof=1)) if values.size > 1 else np.nan,
                    "min": float(np.min(values)) if values.size else np.nan,
                    "max": float(np.max(values)) if values.size else np.nan,
                })

    wins: list[dict[str, object]] = []
    props = sorted({str(row["propeller"]) for row in metric_rows})
    for output in outputs:
        for criterion, key, sign in (("lowest_smape", f"smape_{output}", 1), ("highest_r2", f"r2_{output}", -1)):
            count = {model: 0.0 for model in models}
            eligible = 0
            for prop in props:
                subset = {str(row["model"]): row for row in metric_rows if row["propeller"] == prop}
                values = np.array([safe_float(subset.get(model, {}).get(key)) for model in models])
                if not np.all(np.isfinite(values)):
                    continue
                eligible += 1
                best = np.min(sign * values)
                tied = np.isclose(sign * values, best, rtol=1e-10, atol=1e-12)
                for model, is_best in zip(models, tied):
                    count[model] += float(is_best) / np.count_nonzero(tied)
            for model in models:
                wins.append({"criterion": criterion, "output": output, "model": model, "wins": count[model], "eligible_propellers": eligible})

    rank_rows: list[dict[str, object]] = []
    for prop in props:
        subset = {str(row["model"]): row for row in metric_rows if row["propeller"] == prop}
        ranks = {model: [] for model in models}
        if len(subset) == len(models):
            for output in outputs:
                smape = np.array([safe_float(subset[m][f"smape_{output}"]) for m in models])
                r2 = np.array([safe_float(subset[m][f"r2_{output}"]) for m in models])
                if np.all(np.isfinite(smape)):
                    for model, rank in zip(models, rankdata(smape, method="average")): ranks[model].append(float(rank))
                if np.all(np.isfinite(r2)):
                    for model, rank in zip(models, rankdata(-r2, method="average")): ranks[model].append(float(rank))
        for model in models:
            rank_rows.append({"propeller": prop, "model": model, "metrics_ranked": len(ranks[model]), "mean_rank": float(np.mean(ranks[model])) if ranks[model] else np.nan})

    global_ranks = average_metric_rank(metric_rows)
    global_rank_rows = [{"model": model, "global_mean_rank": global_ranks.get(model, np.nan)} for model in models]

    regime_by_prop = {str(row["propeller"]): row for row in regime_rows}
    trend_rows: list[dict[str, object]] = []
    predictors = ["diameter_in", "rpm", "representative_reynolds", "tip_mach", "j_mean", "ct_reference_mean"]
    for model in models:
        subset = [row for row in metric_rows if row["model"] == model]
        for output in outputs:
            y_all = np.array([safe_float(row[f"smape_{output}"]) for row in subset])
            for predictor in predictors:
                x_all = np.array([safe_float(regime_by_prop.get(str(row["propeller"]), {}).get(predictor)) for row in subset])
                valid = np.isfinite(x_all) & np.isfinite(y_all)
                if np.count_nonzero(valid) >= 3 and np.unique(x_all[valid]).size > 1 and np.unique(y_all[valid]).size > 1:
                    rho, p_value = spearmanr(x_all[valid], y_all[valid])
                else:
                    rho = p_value = np.nan
                trend_rows.append({"model": model, "output": output, "predictor": predictor, "n": int(np.count_nonzero(valid)), "spearman_rho": rho, "p_value": p_value})
    return summary, wins, rank_rows, global_rank_rows, trend_rows


def run(
    limit: int | None,
    only_case: str | None,
    skip_plots: bool,
    max_tip_mach: float,
    jobs: int,
    *,
    output_dir: Path = RESULTS_ROOT,
    include_marine: bool = False,
) -> None:
    """Run the auditable comparison and write tables plus optional figures.

    The default reproduces the aircraft-data scope (224 Clark--Y cases).
    Marine PER3 files are retained in the repository but are included only
    with an explicit ``include_marine=True``.
    """
    require_packaged_data()
    output_dir = Path(output_dir).resolve()
    cases = discover_apc_cases(DATA_ROOT, include_marine=include_marine)
    inventory = inventory_cases(cases, max_tip_mach)
    write_csv(output_dir / "tables" / "propeller_inventory.csv", inventory)
    included = [case for case in cases if case.status == "included"]
    if only_case:
        included = [case for case in included if case.name.lower() == only_case.lower()]
    if limit is not None:
        included = included[:limit]

    metric_rows: list[dict[str, object]] = []
    point_rows: list[dict[str, object]] = []
    issue_rows: list[dict[str, object]] = []
    regime_rows: list[dict[str, object]] = []

    payloads = [(case, max_tip_mach, skip_plots, output_dir) for case in included]
    if jobs > 1:
        executor = ProcessPoolExecutor(max_workers=jobs, initializer=_initialize_worker)
        iterator = executor.map(_process_case, payloads, chunksize=1)
    else:
        executor = None
        _initialize_worker()
        iterator = map(_process_case, payloads)
    for index, (case, output) in enumerate(zip(included, iterator), 1):
        print(f"[{index}/{len(included)}] {case.name}", flush=True)
        new_metrics, new_points, new_issues, regime = output
        metric_rows.extend(new_metrics)
        point_rows.extend(new_points)
        issue_rows.extend(new_issues)
        regime_rows.append(regime)
        write_csv(output_dir / "tables" / "metrics_by_propeller.csv", metric_rows)
        write_csv(output_dir / "tables" / "issues.csv", issue_rows)
    if executor is not None:
        executor.shutdown()

    write_csv(output_dir / "tables" / "curves_all.csv", point_rows)
    write_csv(output_dir / "tables" / "operating_regimes.csv", regime_rows)
    summary, wins, rank_rows, global_rank_rows, trends = aggregate(metric_rows, regime_rows)
    write_csv(output_dir / "aggregate" / "method_summary.csv", summary)
    write_csv(output_dir / "aggregate" / "wins.csv", wins)
    write_csv(output_dir / "aggregate" / "rank_by_propeller.csv", rank_rows)
    write_csv(output_dir / "aggregate" / "global_rank.csv", global_rank_rows)
    write_csv(output_dir / "aggregate" / "regime_trends_spearman.csv", trends)
    if not skip_plots:
        plot_aggregate(metric_rows, wins, rank_rows, output_dir / "aggregate")
    metadata = {
        "included_propellers_requested": len(included), "completed_propellers": len(regime_rows),
        "polar_file": display_path(POLAR_PATH), "polar_domain": PolarDatabase(POLAR_PATH).domain,
        "max_static_tip_mach_for_rpm_selection": max_tip_mach,
        "parallel_jobs": jobs,
        "include_marine": include_marine,
    }
    (output_dir / "aggregate").mkdir(parents=True, exist_ok=True)
    (output_dir / "aggregate" / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int)
    parser.add_argument("--case")
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--max-tip-mach", type=float, default=0.86)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--include-marine", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=RESULTS_ROOT)
    args = parser.parse_args()
    run(
        args.limit,
        args.case,
        args.skip_plots,
        args.max_tip_mach,
        max(1, args.jobs),
        output_dir=args.output_dir,
        include_marine=args.include_marine,
    )


if __name__ == "__main__":
    main()
