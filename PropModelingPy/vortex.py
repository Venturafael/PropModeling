"""McCormick lifting-line / helical-wake Vortex Theory implementation.

The rotating bound circulation sheds a helical trailing wake.  McCormick's
relations iterate tangential induced velocity ``wt`` so that the velocity
induced back at each blade station is consistent with sectional circulation.
Because the wake geometry changes with advance ratio, the method retains more
of the curve-shape physics than a purely local force sum, but its nested
fixed-point iterations and repeated polar queries make it the most expensive
of the three packaged analysis models.
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
    polar_slope,
    representative_station_index,
)


def run_vortex(
    geometry: Geometry,
    j_values: np.ndarray,
    rpm: float,
    polar: PolarDatabase,
    atmosphere: Atmosphere = Atmosphere(),
    *,
    max_iterations: int = 500,
    tolerance: float = 1e-5,
    use_pe0_thickness: bool = True,
) -> ModelResult:
    """Port the selected McCormick MATLAB core with an auditable iteration."""
    j_values = np.asarray(j_values, dtype=float)
    x_values = geometry.x
    beta_deg = np.degrees(np.arctan2(
        geometry.pitch_le_te_in, 2.0 * np.pi * geometry.station_in
    ))
    # The MATLAB loader names this field THICKNESS_RATIO, while the Vortex
    # script checks T_OVER_C/TC/t_over_c and therefore silently falls back to
    # 0.103.  The global analysis fixes that wiring and uses the PE0 geometry;
    # validation can request the exact legacy fallback.
    thickness_values = (
        geometry.thickness_ratio
        if use_pe0_thickness
        else np.full_like(geometry.thickness_ratio, 0.103)
    )
    tip_index = int(np.argmax(x_values))
    tc_tip_index = max(tip_index - 1, 0)
    phi_tip_corrected = np.radians(
        beta_deg[tip_index] + 41.0 * thickness_values[tc_tip_index]
    )
    sin_phi_tip = max(float(np.sin(phi_tip_corrected)), 1e-8)

    ct = np.full(j_values.size, np.nan)
    cp = np.full(j_values.size, np.nan)
    reynolds_rep = np.full(j_values.size, np.nan)
    mach_tip = np.full(j_values.size, np.nan)
    issues: list[dict[str, object]] = []
    rep_i = representative_station_index(x_values)
    polar.reset_stats()

    for jj, advance_ratio in enumerate(j_values):
        lam = advance_ratio / np.pi
        dct_dx = np.zeros(x_values.size)
        dcp_dx = np.zeros(x_values.size)
        failed = 0

        for k, (x, radius, chord, thickness, beta) in enumerate(zip(
            x_values,
            geometry.station_in,
            geometry.chord_in,
            thickness_values,
            beta_deg,
        )):
            sigma_mean = geometry.blades * chord / (np.pi * geometry.radius_in)
            wt = 0.010
            converged = False
            cl = cd = phi_relative = float("nan")
            re_raw, ma_raw, _ = local_reynolds_mach(
                radius, chord, advance_ratio, rpm, geometry.radius_in, atmosphere
            )
            if k == rep_i:
                reynolds_rep[jj] = re_raw
            if k == tip_index:
                mach_tip[jj] = ma_raw

            # Fixed point: wake-induced tangential velocity changes inflow;
            # inflow changes circulation, which in turn changes the wake.
            for _ in range(max_iterations):
                radicand = lam**2 + 4.0 * wt * (x - wt)
                if radicand <= 0.0:
                    break
                wa = 0.5 * (-lam + np.sqrt(radicand))
                phi_relative = np.arctan2(lam + wa, x - wt)
                # Keep McCormick/MATLAB's one-argument atan convention here.
                # atan2 would add pi when x-2*wt is negative and materially
                # change the empirical Delta-theta correction.
                denominator_theta = x - 2.0 * wt
                if abs(denominator_theta) < 1e-14:
                    break
                delta_theta = np.degrees(
                    np.arctan((lam + wa) / denominator_theta)
                    - np.arctan((lam + wa) / x)
                )
                delta_alpha = (
                    4.0 * lam * sigma_mean / (15.0 * (lam**2 + x**2))
                    * thickness * 180.0 / np.pi
                )
                alpha_geometric = beta - np.degrees(phi_relative) - delta_alpha
                query = polar.lookup(alpha_geometric, re_raw, ma_raw)
                slope = polar_slope(polar, alpha_geometric, re_raw, ma_raw)
                if not np.isfinite(slope) or slope <= 0.0:
                    slope = 0.1 * (1.0 + thickness)
                cl = query.cl - slope * delta_theta / 4.0
                cd = query.cd
                relative_non_dimensional = np.hypot(lam + wa, x - wt)
                f_tip = (geometry.blades / 2.0) * (1.0 - x) / sin_phi_tip
                prandtl = max(2.0 / np.pi * np.arccos(np.exp(-max(f_tip, 0.0))), 0.05)
                wt_new = sigma_mean * cl * relative_non_dimensional / (8.0 * x * prandtl)
                wt_new = float(np.clip(wt_new, -0.9 * x, 0.9 * x))
                wt_new = 0.55 * wt_new + 0.45 * wt
                if abs(wt_new - wt) < tolerance:
                    wt = wt_new
                    converged = True
                    break
                wt = wt_new

            if not converged or not np.all(np.isfinite([cl, cd, phi_relative])):
                failed += 1
                continue
            speed_squared_term = advance_ratio**2 + (np.pi * x) ** 2
            dct_dx[k] = (
                np.pi / 8.0 * speed_squared_term * sigma_mean
                * (cl * np.cos(phi_relative) - cd * np.sin(phi_relative))
            )
            dcp_dx[k] = (
                np.pi / 8.0 * (np.pi * x) * speed_squared_term * sigma_mean
                * (cl * np.sin(phi_relative) + cd * np.cos(phi_relative))
            )

        if failed:
            issues.append({
                "model": "Vortex", "j": float(advance_ratio), "category": "convergence",
                "reason": f"fixed-point solver failed at {failed}/{x_values.size} radial stations",
                "treatment": "legacy treatment retained: failed station contribution set to zero",
            })
        ct[jj] = np.trapz(dct_dx, x_values)
        cp[jj] = np.trapz(dcp_dx, x_values)
    return finalize_result("Vortex", j_values, ct, cp, reynolds_rep, mach_tip, issues, polar)
