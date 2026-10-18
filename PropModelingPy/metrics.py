"""Metrics used in the SETI article and transparent aggregation helpers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import rankdata


@dataclass(frozen=True)
class CurveMetrics:
    n: int
    smape: float
    pearson_r2: float


def curve_metrics(reference: np.ndarray, predicted: np.ndarray) -> CurveMetrics:
    reference = np.asarray(reference, dtype=float).ravel()
    predicted = np.asarray(predicted, dtype=float).ravel()
    valid = np.isfinite(reference) & np.isfinite(predicted)
    y = reference[valid]
    yhat = predicted[valid]
    if y.size == 0:
        return CurveMetrics(0, float("nan"), float("nan"))
    denominator = np.abs(yhat) + np.abs(y)
    term = np.zeros_like(y)
    nonzero = denominator > np.finfo(float).eps
    term[nonzero] = 2.0 * np.abs(yhat[nonzero] - y[nonzero]) / denominator[nonzero]
    smape = 100.0 * float(np.mean(term))
    if y.size < 2 or np.isclose(np.var(y), 0.0) or np.isclose(np.var(yhat), 0.0):
        pearson_r2 = float("nan")
    else:
        pearson_r2 = float(np.corrcoef(y, yhat)[0, 1] ** 2)
    return CurveMetrics(int(y.size), smape, pearson_r2)


def average_metric_rank(rows: list[dict[str, object]]) -> dict[str, float]:
    """Average per-propeller rank across 3 sMAPE and 3 r2 metrics."""
    models = sorted({str(row["model"]) for row in rows})
    propellers = sorted({str(row["propeller"]) for row in rows})
    ranks: dict[str, list[float]] = {model: [] for model in models}
    for propeller in propellers:
        subset = {str(row["model"]): row for row in rows if row["propeller"] == propeller}
        if len(subset) != len(models):
            continue
        for metric in ("ct", "cp", "eta"):
            smape = np.array([float(subset[m][f"smape_{metric}"]) for m in models])
            r2 = np.array([float(subset[m][f"r2_{metric}"]) for m in models])
            if np.all(np.isfinite(smape)):
                for model, rank in zip(models, rankdata(smape, method="average")):
                    ranks[model].append(float(rank))
            if np.all(np.isfinite(r2)):
                for model, rank in zip(models, rankdata(-r2, method="average")):
                    ranks[model].append(float(rank))
    return {model: float(np.mean(values)) if values else float("nan") for model, values in ranks.items()}

