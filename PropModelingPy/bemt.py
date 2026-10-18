"""Lowry Blade Element Momentum Theory (BEMT).

BEMT requires each annulus to satisfy two descriptions simultaneously: the
forces predicted from the local airfoil polar and the momentum change of the
streamtube.  The unknown section angle of attack is therefore solved from a
nonlinear residual.  Prandtl's factor weakens loading near the finite blade
tip.  Compared with BET, this closure usually improves force magnitude at the
cost of a root solve at every radial station and advance ratio.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from polar_lookup import PolarDatabase
from apc_data import Geometry
from model_common import (
    Atmosphere,
    ModelResult,
    finalize_result,
    local_reynolds_mach,
    prandtl_factor,
    representative_station_index,
)


def _legacy_aero_coefficients(alpha_deg: float) -> tuple[float, float]:
    c1 = 0.09458 * (alpha_deg + 4.8)
    cl = (1.5 + c1 - np.sqrt((1.5 - c1) ** 2 + 0.04)) / 2.0
    cd = 0.015 + 0.03 * ((4.0 * c1 / 3.0) - 1.0) ** 4
    return float(cl), float(cd)


def _residual(
    alpha: float,
    beta: float,
    j: float,
    x: float,
    sigma: float,
    blades: int,
    reynolds: float,
    mach: float,
    polar: PolarDatabase,
    legacy_solver_polar: bool,
) -> float:
    """Difference between annular momentum and blade-element advance ratio."""
    phi = beta - alpha
    phi_rad = np.radians(phi)
    if legacy_solver_polar:
        cl, cd = _legacy_aero_coefficients(alpha)
    else:
        query = polar.lookup(alpha, reynolds, mach, record=False)
        cl, cd = query.cl, query.cd
    f = prandtl_factor(blades, x, phi_rad)
    if not np.all(np.isfinite([cl, cd, f])):
        return float("nan")
    lambda_t = cl * np.cos(phi_rad) - cd * np.sin(phi_rad)
    lambda_p = cl * np.sin(phi_rad) + cd * np.cos(phi_rad)
    numerator = np.pi * x * (4.0 * f * np.sin(phi_rad) ** 2 - lambda_t * sigma)
    denominator = 4.0 * f * np.sin(phi_rad) * np.cos(phi_rad) + lambda_p * sigma
    return float(numerator / denominator - j) if abs(denominator) > 1e-14 else float("nan")


def _solve_alpha(
    beta: float,
    j: float,
    x: float,
    sigma: float,
    blades: int,
    reynolds: float,
    mach: float,
    polar: PolarDatabase,
    legacy_solver_polar: bool,
) -> float:
    """Bracket and solve the physically coupled sectional angle of attack."""
    grid = np.arange(-10.0, 25.0001, 0.5)
    phi = beta - grid
    phi_rad = np.radians(phi)
    if legacy_solver_polar:
        c1 = 0.09458 * (grid + 4.8)
        cl = (1.5 + c1 - np.sqrt((1.5 - c1) ** 2 + 0.04)) / 2.0
        cd = 0.015 + 0.03 * ((4.0 * c1 / 3.0) - 1.0) ** 4
    else:
        cl, cd = polar.lookup_many(grid, reynolds, mach)
    denominator_f = x * np.sin(phi_rad)
    f = np.full_like(grid, np.nan)
    valid_f = denominator_f > 0.0
    argument = np.clip(np.exp(-(blades / 2.0) * (1.0 - x) / denominator_f[valid_f]), 0.0, 1.0)
    f[valid_f] = np.maximum(2.0 / np.pi * np.arccos(argument), 0.05)
    lambda_t = cl * np.cos(phi_rad) - cd * np.sin(phi_rad)
    lambda_p = cl * np.sin(phi_rad) + cd * np.cos(phi_rad)
    numerator = np.pi * x * (4.0 * f * np.sin(phi_rad) ** 2 - lambda_t * sigma)
    denominator = 4.0 * f * np.sin(phi_rad) * np.cos(phi_rad) + lambda_p * sigma
    values = numerator / denominator - j
    for idx in range(grid.size - 1):
        if not np.all(np.isfinite(values[idx : idx + 2])):
            continue
        if values[idx] == 0.0:
            return float(grid[idx])
        if values[idx] * values[idx + 1] < 0.0:
            try:
                return float(brentq(
                    _residual, grid[idx], grid[idx + 1],
                    args=(beta, j, x, sigma, blades, reynolds, mach, polar, legacy_solver_polar),
                    xtol=1e-10, rtol=1e-12, maxiter=100,
                ))
            except (ValueError, RuntimeError):
                continue
    return float("nan")


def run_bemt(
    geometry: Geometry,
    j_values: np.ndarray,
    rpm: float,
    polar: PolarDatabase,
    atmosphere: Atmosphere = Atmosphere(),
    *,
    legacy_solver_polar: bool = False,
) -> ModelResult:
    """Evaluate BEMT with the shared polar also inside the alpha solver.

    ``legacy_solver_polar=True`` exists only for MATLAB migration validation;
    it reproduces the hidden ``aeroCoeffs.m`` dependency and is never used by
    the global comparison.
    """
    j_values = np.asarray(j_values, dtype=float)
    x_values = geometry.x
    diameter_in = geometry.diameter_in
    ct = np.full(j_values.size, np.nan)
    cp = np.full(j_values.size, np.nan)
    reynolds_rep = np.full(j_values.size, np.nan)
    mach_tip = np.full(j_values.size, np.nan)
    issues: list[dict[str, object]] = []
    rep_i = representative_station_index(x_values)
    polar.reset_stats()

    for jj, advance_ratio in enumerate(j_values):
        if advance_ratio < 0.05:
            issues.append({
                "model": "BEMT", "j": float(advance_ratio), "category": "model_domain",
                "reason": "legacy BEMT validity cutoff J<0.05",
                "treatment": "CT, CP, eta retained as NaN",
            })
            continue
        dct_dx = np.zeros(x_values.size)
        dcp_dx = np.zeros(x_values.size)
        valid = np.ones(x_values.size, dtype=bool)
        for i, (x, radius, chord, pitch) in enumerate(zip(
            x_values, geometry.station_in, geometry.chord_in, geometry.pitch_le_te_in
        )):
            beta = np.degrees(np.arctan(pitch / (2.0 * np.pi * radius)))
            sigma = geometry.blades * chord / (2.0 * np.pi * x * geometry.radius_in)
            re_i, ma_i, _ = local_reynolds_mach(
                radius, chord, advance_ratio, rpm, geometry.radius_in, atmosphere
            )
            if i == rep_i:
                reynolds_rep[jj] = re_i
            if i == x_values.size - 1:
                mach_tip[jj] = ma_i
            # This is the central BEMT closure: alpha changes Cl/Cd, those
            # forces change induction, and induction changes the inflow angle.
            alpha = _solve_alpha(
                beta, advance_ratio, x, sigma, geometry.blades, re_i, ma_i,
                polar, legacy_solver_polar,
            )
            if not np.isfinite(alpha):
                valid[i] = False
                continue
            phi = beta - alpha
            phi_rad = np.radians(phi)
            query = polar.lookup(alpha, re_i, ma_i)
            f = prandtl_factor(geometry.blades, x, phi_rad)
            sin_phi = np.sin(phi_rad)
            if not np.all(np.isfinite([query.cl, query.cd, f])) or abs(sin_phi) < 1e-9:
                valid[i] = False
                continue
            lambda_t = query.cl * np.cos(phi_rad) - query.cd * sin_phi
            lambda_p = query.cl * sin_phi + query.cd * np.cos(phi_rad)
            if abs(lambda_t * sigma) < 1e-14:
                valid[i] = False
                continue
            a_ax = 1.0 / (4.0 * f * sin_phi**2 / (lambda_t * sigma) - 1.0)
            factor = (1.0 + a_ax) ** 2 / sin_phi**2
            dct_dx[i] = 0.25 * advance_ratio**2 * factor * chord * geometry.blades * lambda_t / diameter_in
            dcp_dx[i] = 0.25 * np.pi * x * advance_ratio**2 * factor * chord * geometry.blades * lambda_p / diameter_in
            if not np.all(np.isfinite([dct_dx[i], dcp_dx[i]])):
                valid[i] = False
        failed = int(np.count_nonzero(~valid))
        if failed:
            issues.append({
                "model": "BEMT", "j": float(advance_ratio), "category": "numerical",
                "reason": f"alpha solver or local equation failed at {failed}/{valid.size} radial stations",
                "treatment": "legacy treatment retained: failed station contribution set to zero",
            })
        ct[jj] = np.trapz(dct_dx, x_values)
        cp[jj] = np.trapz(dcp_dx, x_values)
    return finalize_result("BEMT", j_values, ct, cp, reynolds_rep, mach_tip, issues, polar)
