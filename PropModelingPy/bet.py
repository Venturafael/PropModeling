"""Lowry simple Blade Element Theory (BET).

Each radial station is treated as an independent two-dimensional airfoil.
The inflow angle follows directly from advance ratio and radius; BET does not
feed the blade loading back into an induced velocity.  Lift and drag are
resolved into axial/tangential components and integrated to obtain ``CT`` and
``CP``.  This makes BET inexpensive but leaves wake contraction, induction and
tip-loss coupling outside the model.
"""

from __future__ import annotations

import numpy as np

from polar_lookup import PolarDatabase
from apc_data import Geometry
from model_common import (
    Atmosphere,
    ModelResult,
    finalize_result,
    local_reynolds_mach,
    representative_station_index,
)


def run_bet(
    geometry: Geometry,
    j_values: np.ndarray,
    rpm: float,
    polar: PolarDatabase,
    atmosphere: Atmosphere = Atmosphere(),
) -> ModelResult:
    """Evaluate the selected Lowry MATLAB equations without induced flow."""
    j_values = np.asarray(j_values, dtype=float)
    x = geometry.x
    diameter_in = geometry.diameter_in
    n_j = j_values.size
    ct = np.full(n_j, np.nan)
    cp = np.full(n_j, np.nan)
    reynolds_rep = np.full(n_j, np.nan)
    mach_tip = np.full(n_j, np.nan)
    issues: list[dict[str, object]] = []
    rep_i = representative_station_index(x)
    polar.reset_stats()

    for jj, advance_ratio in enumerate(j_values):
        if advance_ratio <= 0.0:
            issues.append({
                "model": "BET", "j": float(advance_ratio), "category": "model_domain",
                "reason": "Lowry J-form equations are singular at J=0",
                "treatment": "CT, CP, eta retained as NaN",
            })
            continue
        dct_dx = np.full(x.size, np.nan)
        dcp_dx = np.full(x.size, np.nan)
        for i, (xi, radius, chord, pitch) in enumerate(zip(
            x, geometry.station_in, geometry.chord_in, geometry.pitch_le_te_in
        )):
            # Geometric pitch angle minus purely kinematic inflow gives the
            # sectional angle of attack used by the common polar database.
            beta = np.degrees(np.arctan(pitch / (2.0 * np.pi * radius)))
            phi = np.degrees(np.arctan(advance_ratio / (np.pi * xi)))
            alpha = beta - phi
            re_i, ma_i, _ = local_reynolds_mach(
                radius, chord, advance_ratio, rpm, geometry.radius_in, atmosphere
            )
            query = polar.lookup(alpha, re_i, ma_i)
            if i == rep_i:
                reynolds_rep[jj] = re_i
            if i == x.size - 1:
                mach_tip[jj] = ma_i
            sin_phi = np.sin(np.radians(phi))
            if abs(sin_phi) < 1e-12 or not np.all(np.isfinite([query.cl, query.cd])):
                continue
            # Resolve sectional lift/drag into thrust and torque directions.
            lambda_t = query.cl * np.cos(np.radians(phi)) - query.cd * sin_phi
            lambda_p = query.cl * sin_phi + query.cd * np.cos(np.radians(phi))
            dct_dx[i] = (
                geometry.blades * chord / (4.0 * diameter_in)
                * advance_ratio**2 * lambda_t / sin_phi**2
            )
            dcp_dx[i] = (
                np.pi * geometry.blades * chord * xi / (4.0 * diameter_in)
                * advance_ratio**2 * lambda_p / sin_phi**2
            )
        valid = np.isfinite(dct_dx) & np.isfinite(dcp_dx)
        if np.count_nonzero(valid) < 2:
            issues.append({
                "model": "BET", "j": float(advance_ratio), "category": "numerical",
                "reason": "fewer than two valid radial stations",
                "treatment": "CT, CP, eta retained as NaN",
            })
            continue
        if not np.all(valid):
            issues.append({
                "model": "BET", "j": float(advance_ratio), "category": "numerical",
                "reason": f"{np.count_nonzero(~valid)} invalid radial stations",
                "treatment": "invalid stations omitted from radial trapezoidal integration",
            })
        ct[jj] = np.trapz(dct_dx[valid], x[valid])
        cp[jj] = np.trapz(dcp_dx[valid], x[valid])
    return finalize_result("BET", j_values, ct, cp, reynolds_rep, mach_tip, issues, polar)

