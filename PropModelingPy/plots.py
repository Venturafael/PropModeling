"""Plotting with fixed colors, names, and English labels."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from apc_data import PerformanceBlock
from model_common import ModelResult


COLORS = {"APC reference": "#111111", "BET": "#2a9d68", "BEMT": "#2457a6", "Vortex": "#c43c35"}
LABELS = {"BET": "BET", "BEMT": "BEMT", "Vortex": "Vortex Theory"}
STYLES = {"APC reference": ":", "BET": "-", "BEMT": "-.", "Vortex": "--"}


def _finish_axis(axis: plt.Axes, ylabel: str) -> None:
    axis.set_xlabel("Advance ratio, J")
    axis.set_ylabel(ylabel)
    axis.grid(True, alpha=0.22)
    axis.spines[["top", "right"]].set_visible(False)


def plot_individual(
    propeller: str,
    rpm: int,
    reference: PerformanceBlock,
    results: dict[str, ModelResult],
    output_path: str | Path,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3), constrained_layout=True)
    series = (("ct", r"$C_T$"), ("cp", r"$C_P$"), ("eta", r"$\eta$"))
    for axis, (field, ylabel) in zip(axes, series):
        axis.plot(
            reference.j, getattr(reference, field), color=COLORS["APC reference"],
            linestyle=STYLES["APC reference"], marker="o", markersize=3.2,
            linewidth=2.0, label="APC reference",
        )
        for model in ("BET", "BEMT", "Vortex"):
            result = results.get(model)
            if result is None:
                continue
            axis.plot(
                result.j, getattr(result, field), color=COLORS[model],
                linestyle=STYLES[model], linewidth=2.0, label=LABELS[model],
            )
        _finish_axis(axis, ylabel)
        if field == "eta":
            axis.set_ylim(bottom=0.0)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.04),
        ncol=4, frameon=False,
    )
    fig.suptitle(f"APC {propeller} — {rpm:,} rpm".replace(",", " "), y=1.08, fontsize=14)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_aggregate(
    metric_rows: list[dict[str, object]],
    win_rows: list[dict[str, object]],
    rank_rows: list[dict[str, object]],
    output_directory: str | Path,
) -> None:
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    models = ["BET", "BEMT", "Vortex"]
    outputs = ["ct", "cp", "eta"]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
    for axis, output in zip(axes, outputs):
        data = [
            [float(row[f"smape_{output}"]) for row in metric_rows if row["model"] == model and np.isfinite(float(row[f"smape_{output}"]))]
            for model in models
        ]
        boxes = axis.boxplot(data, labels=[LABELS[m] for m in models], patch_artist=True, showfliers=False)
        for patch, model in zip(boxes["boxes"], models):
            patch.set_facecolor(COLORS[model]); patch.set_alpha(0.72)
        axis.set_title(output.upper() if output != "eta" else r"$\eta$")
        axis.set_ylabel("sMAPE [%]")
        axis.grid(axis="y", alpha=0.22)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Error distribution across confirmed APC Clark–Y propellers")
    fig.savefig(output_directory / "smape_boxplots.png", dpi=190, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    x = np.arange(len(models))
    for axis, output in zip(axes, outputs):
        values = []
        for model in models:
            match = next((row for row in win_rows if row["criterion"] == "lowest_smape" and row["output"] == output and row["model"] == model), None)
            values.append(int(match["wins"]) if match else 0)
        bars = axis.bar(x, values, color=[COLORS[m] for m in models])
        axis.bar_label(bars, fontsize=9)
        axis.set_xticks(x, [LABELS[m] for m in models])
        axis.set_ylabel("Propeller wins")
        axis.set_title(output.upper() if output != "eta" else r"$\eta$")
        axis.grid(axis="y", alpha=0.22)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Per-propeller wins by lowest sMAPE")
    fig.savefig(output_directory / "smape_wins.png", dpi=190, bbox_inches="tight")
    plt.close(fig)

    if rank_rows:
        props = sorted({str(row["propeller"]) for row in rank_rows})
        matrix = np.full((len(models), len(props)), np.nan)
        for row in rank_rows:
            matrix[models.index(str(row["model"])), props.index(str(row["propeller"]))] = float(row["mean_rank"])
        fig_width = max(14, len(props) * 0.07)
        fig, axis = plt.subplots(figsize=(fig_width, 3.2), constrained_layout=True)
        image = axis.imshow(matrix, aspect="auto", cmap="RdYlGn_r", vmin=1, vmax=3)
        axis.set_yticks(np.arange(len(models)), [LABELS[m] for m in models])
        step = max(1, len(props) // 25)
        ticks = np.arange(0, len(props), step)
        axis.set_xticks(ticks, [props[i] for i in ticks], rotation=70, ha="right", fontsize=7)
        axis.set_title("Mean rank per propeller (six metrics; 1 is best)")
        fig.colorbar(image, ax=axis, label="Mean rank")
        fig.savefig(output_directory / "rank_by_propeller.png", dpi=190, bbox_inches="tight")
        plt.close(fig)
