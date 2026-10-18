"""Shared model types and numerical utilities."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from polar_lookup import PolarDatabase


@dataclass(frozen=True)
class Atmosphere:
    density_kg_m3: float = 1.225
    dynamic_viscosity_pa_s: float = 1.81e-5
    speed_of_sound_m_s: float = 343.0


@dataclass
class ModelResult:
    model: str
    j: np.ndarray
    ct: np.ndarray
    cp: np.ndarray
    eta: np.ndarray
    reynolds_representative: np.ndarray
    mach_tip: np.ndarray
    issues: list[dict[str, Any]] = field(default_factory=list)
    polar_stats: dict[str, int] = field(default_factory=dict)


def local_reynolds_mach(
    radius_in: float,
    chord_in: float,
    j: float,
    rpm: float,
    prop_radius_in: float,
    atmosphere: Atmosphere,
) -> tuple[float, float, float]:
    inch = 0.0254
    radius_m = radius_in * inch
    chord_m = chord_in * inch
    diameter_m = 2.0 * prop_radius_in * inch
    rotations_s = rpm / 60.0
    omega = 2.0 * np.pi * rotations_s
    freestream = j * rotations_s * diameter_m
    tangential = omega * radius_m
    relative = float(np.hypot(tangential, freestream))
    reynolds = atmosphere.density_kg_m3 * relative * chord_m / atmosphere.dynamic_viscosity_pa_s
    mach = relative / atmosphere.speed_of_sound_m_s
    return float(reynolds), float(mach), relative


def finalize_result(
    model: str,
    j: np.ndarray,
    ct: np.ndarray,
    cp: np.ndarray,
    reynolds_rep: np.ndarray,
    mach_tip: np.ndarray,
    issues: list[dict[str, Any]],
    polar: PolarDatabase,
) -> ModelResult:
    eta = np.full_like(ct, np.nan, dtype=float)
    propulsive = np.isfinite(ct) & np.isfinite(cp) & (ct > 0.0) & (cp > 0.0)
    eta[propulsive] = j[propulsive] * ct[propulsive] / cp[propulsive]
    for idx in np.where(~propulsive)[0]:
        issues.append({
            "model": model,
            "j": float(j[idx]),
            "category": "physical_regime",
            "reason": "efficiency undefined because CT<=0, CP<=0, or coefficient non-finite",
            "treatment": "eta retained as NaN; CT and CP retained",
        })
    return ModelResult(
        model=model,
        j=j.copy(), ct=ct, cp=cp, eta=eta,
        reynolds_representative=reynolds_rep,
        mach_tip=mach_tip,
        issues=issues,
        polar_stats=polar.stats.as_dict(),
    )


def representative_station_index(x: np.ndarray, target: float = 0.75) -> int:
    return int(np.argmin(np.abs(x - target)))


def polar_slope(
    polar: PolarDatabase,
    alpha_deg: float,
    reynolds: float,
    mach: float,
    delta_deg: float = 0.5,
) -> float:
    high = polar.lookup(alpha_deg + delta_deg, reynolds, mach, record=False).cl
    low = polar.lookup(alpha_deg - delta_deg, reynolds, mach, record=False).cl
    return float((high - low) / (2.0 * delta_deg))


def prandtl_factor(blades: int, x: float, phi_rad: float) -> float:
    denominator = x * np.sin(phi_rad)
    if not np.isfinite(denominator) or denominator <= 0.0:
        return float("nan")
    exponent = -(blades / 2.0) * (1.0 - x) / denominator
    argument = float(np.clip(np.exp(exponent), 0.0, 1.0))
    return float(max(2.0 / np.pi * np.arccos(argument), 0.05))

