"""Adkins--Liebeck Minimum Energy Loss (MEL) design solver.

This module is a direct, modular migration of ``adkins_uai_TUNADO.m`` and
the ``mel_*`` MATLAB helpers selected during the source audit.  The default
inputs reproduce the TCC design point; they are not universal defaults for
all propellers.

Important reproducibility caveat
--------------------------------
The MATLAB/TCC implementation applies ``alpha_bias_deg = -3.5`` after
finding the angle that provides the requested design lift coefficient.  The
biased angle is then used to query drag, while chord and the loss ratio retain
the requested ``Cl``.  This empirical choice is preserved verbatim so that
the port remains traceable.  It creates an aerodynamic-state inconsistency
that should be investigated before interpreting the result as a new validated
physical optimum.  Set ``alpha_bias_deg=0`` only as an explicit sensitivity
case, not as an undocumented correction.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from polar_lookup import PolarDatabase
from apc_data import Geometry


@dataclass(frozen=True)
class MELInputs:
    """Operating/design inputs used by the TCC MEL implementation."""

    blades: int = 2
    diameter_m: float = 0.559
    density_kg_m3: float = 1.225
    kinematic_viscosity_m2_s: float = 1.47e-5
    velocity_m_s: float = 11.9
    rpm: float = 4300.0
    shaft_power_w: float = 500.0
    alpha_bias_deg: float = -3.5
    cd_factor: float = 1.0
    activity_factor_min: float = 70.0
    speed_of_sound_m_s: float = 343.0
    prandtl_min: float = 0.05
    axial_induction_min: float = -0.20
    axial_induction_max: float = 0.60
    tangential_induction_min: float = -0.50
    tangential_induction_max: float = 0.50
    local_tolerance: float = 1e-6

    @property
    def rotations_s(self) -> float:
        return self.rpm / 60.0

    @property
    def radius_m(self) -> float:
        return self.diameter_m / 2.0

    @property
    def advance_ratio(self) -> float:
        return self.velocity_m_s / (self.rotations_s * self.diameter_m)

    @property
    def lambda_value(self) -> float:
        return self.advance_ratio / np.pi

    @property
    def power_coefficient_target(self) -> float:
        return (
            2.0
            * self.shaft_power_w
            / (
                self.density_kg_m3
                * self.velocity_m_s**3
                * np.pi
                * self.radius_m**2
            )
        )


def inverse_alpha(
    cl_target: float,
    reynolds: float,
    mach: float,
    polar: PolarDatabase,
) -> float:
    """Return the global least-squares inverse of ``Cl(alpha)``.

    The coarse search followed by bounded minimization reproduces the branch
    choice of MATLAB ``fminbnd`` while avoiding thousands of MAT-file queries.
    """

    grid = np.linspace(polar.alpha_min, polar.alpha_max, 241)
    cl_grid, _ = polar.lookup_many(grid, reynolds, mach)
    objective = (cl_grid - cl_target) ** 2
    index = int(np.nanargmin(objective))
    lower = float(grid[max(index - 1, 0)])
    upper = float(grid[min(index + 1, grid.size - 1)])
    if lower == upper:
        return float(grid[index])
    result = minimize_scalar(
        lambda alpha: (
            polar.lookup(float(alpha), reynolds, mach, record=False).cl
            - cl_target
        )
        ** 2,
        bounds=(lower, upper),
        method="bounded",
        options={"xatol": 1e-6, "maxiter": 40},
    )
    return float(result.x)


def mel_angles(
    zeta: float,
    radius_fraction: np.ndarray,
    inputs: MELInputs,
) -> dict[str, np.ndarray]:
    """Evaluate helical-flow angles and Prandtl/MEL loading functions."""

    xi = np.asarray(radius_fraction, dtype=float)
    lam = inputs.lambda_value
    x = xi / lam
    phi_tip = np.arctan(lam * (1.0 + 0.5 * zeta))
    f_tip = (inputs.blades / 2.0) * (1.0 - xi) / max(np.sin(phi_tip), 1e-12)
    f = np.maximum(
        (2.0 / np.pi) * np.arccos(np.exp(-np.maximum(f_tip, 0.0))),
        inputs.prandtl_min,
    )
    phi = np.arctan(np.tan(phi_tip) / np.maximum(xi, 1e-8))
    sin_phi = np.maximum(np.sin(phi), 1e-8)
    cos_phi = np.cos(phi)
    tan_phi = np.tan(phi)
    g = f * x * cos_phi * sin_phi
    return {
        "x": x,
        "phi": phi,
        "sin_phi": sin_phi,
        "cos_phi": cos_phi,
        "tan_phi": tan_phi,
        "f": f,
        "g": g,
    }


def mel_state(
    zeta: float,
    cl_target: float,
    radius_fraction: np.ndarray,
    inputs: MELInputs,
    polar: PolarDatabase,
) -> dict[str, object]:
    """Evaluate one fully coupled MEL state for a prescribed ``zeta``."""

    xi = np.asarray(radius_fraction, dtype=float)
    if xi.ndim != 1 or xi.size < 2 or not np.all(np.diff(xi) > 0.0):
        raise ValueError("radius_fraction must be a strictly increasing 1-D array")
    if np.any((xi <= 0.0) | (xi > 1.0)):
        raise ValueError("radius_fraction entries must satisfy 0 < r/R <= 1")
    if cl_target <= 0.0:
        raise ValueError("cl_target must be positive")

    angles = mel_angles(zeta, xi, inputs)
    wc = (
        4.0
        * np.pi
        * inputs.lambda_value
        * angles["g"]
        * inputs.velocity_m_s
        * inputs.radius_m
        * zeta
        / (cl_target * inputs.blades)
    )
    names = (
        "a",
        "ap",
        "w",
        "c",
        "re",
        "mach",
        "alpha",
        "alpha_target",
        "cd",
    )
    arrays = {name: np.zeros(xi.size) for name in names}

    for k in range(xi.size):
        a_value = 0.05
        ap = relative = chord = reynolds = mach = alpha_target = alpha = cd = np.nan
        for _ in range(100):
            relative = inputs.velocity_m_s * (1.0 + a_value) / angles["sin_phi"][k]
            chord = wc[k] / max(relative, 1e-12)
            reynolds = relative * chord / max(inputs.kinematic_viscosity_m2_s, 1e-12)
            mach = relative / inputs.speed_of_sound_m_s
            alpha_target = inverse_alpha(cl_target, reynolds, mach, polar)
            # Preserved TCC/MATLAB empirical convention; see module caveat.
            alpha = inputs.alpha_bias_deg + alpha_target
            cd = inputs.cd_factor * polar.lookup(
                alpha, reynolds, mach, record=False
            ).cd
            epsilon = max(cd / max(cl_target, 1e-12), 1e-12)
            a_new = (
                0.5
                * zeta
                * angles["cos_phi"][k] ** 2
                * (1.0 - epsilon * angles["tan_phi"][k])
            )
            ap = (
                0.5
                * zeta
                * angles["cos_phi"][k]
                * angles["sin_phi"][k]
                * (1.0 + epsilon / max(angles["tan_phi"][k], 1e-12))
                / max(angles["x"][k], 1e-12)
            )
            a_new = np.clip(
                a_new, inputs.axial_induction_min, inputs.axial_induction_max
            )
            ap = np.clip(
                ap,
                inputs.tangential_induction_min,
                inputs.tangential_induction_max,
            )
            updated = 0.7 * a_value + 0.3 * a_new
            if abs(updated - a_value) < inputs.local_tolerance:
                a_value = float(updated)
                break
            a_value = float(updated)
        arrays["a"][k] = a_value
        arrays["ap"][k] = ap
        arrays["w"][k] = relative
        arrays["c"][k] = chord
        arrays["re"][k] = reynolds
        arrays["mach"][k] = mach
        arrays["alpha"][k] = alpha
        arrays["alpha_target"][k] = alpha_target
        arrays["cd"][k] = cd

    epsilon = np.maximum(arrays["cd"] / max(cl_target, 1e-12), 1e-12)
    i1p = 4.0 * xi * angles["g"] * (1.0 - epsilon * angles["tan_phi"])
    i2p = (
        inputs.lambda_value
        * (i1p / np.maximum(2.0 * xi, 1e-12))
        * (1.0 + epsilon / np.maximum(angles["tan_phi"], 1e-12))
        * angles["sin_phi"]
        * angles["cos_phi"]
    )
    j1p = 4.0 * xi * angles["g"] * (
        1.0 + epsilon / np.maximum(angles["tan_phi"], 1e-12)
    )
    j2p = (
        0.5
        * j1p
        * (1.0 - epsilon * angles["tan_phi"])
        * angles["cos_phi"] ** 2
    )
    dxi = np.gradient(xi)
    i1, i2, j1, j2 = [
        float(np.sum(integrand * dxi)) for integrand in (i1p, i2p, j1p, j2p)
    ]
    tc = i1 * zeta - i2 * zeta**2
    pc = j1 * zeta + j2 * zeta**2
    disk_area = np.pi * inputs.radius_m**2
    thrust = 0.5 * inputs.density_kg_m3 * inputs.velocity_m_s**2 * disk_area * tc
    aerodynamic_power = (
        0.5
        * inputs.density_kg_m3
        * inputs.velocity_m_s**3
        * disk_area
        * pc
    )
    torque = aerodynamic_power / (2.0 * np.pi * inputs.rotations_s)
    ct = thrust / (
        inputs.density_kg_m3 * inputs.rotations_s**2 * inputs.diameter_m**4
    )
    cp = aerodynamic_power / (
        inputs.density_kg_m3 * inputs.rotations_s**3 * inputs.diameter_m**5
    )
    efficiency = inputs.advance_ratio * ct / max(cp, 1e-12)
    beta = np.degrees(angles["phi"]) + arrays["alpha"]
    activity_factor = float(
        (100000.0 / 16.0)
        * inputs.blades
        * np.sum((arrays["c"] / inputs.diameter_m) * xi**3 * dxi)
    )
    kct = (
        0.5
        * inputs.velocity_m_s**2
        * disk_area
        / (inputs.rotations_s**2 * inputs.diameter_m**4)
    )
    kcp = (
        0.5
        * inputs.velocity_m_s**3
        * disk_area
        / (inputs.rotations_s**3 * inputs.diameter_m**5)
    )
    return {
        "zeta": float(zeta),
        "cl_target": float(cl_target),
        "xi": xi,
        **angles,
        **arrays,
        "beta": beta,
        "epsilon": epsilon,
        "i1p": i1p,
        "i2p": i2p,
        "j1p": j1p,
        "j2p": j2p,
        "dct_dxi": kct * (i1p * zeta - i2p * zeta**2),
        "dcp_dxi": kcp * (j1p * zeta + j2p * zeta**2),
        "thrust_n": float(thrust),
        "power_w": float(aerodynamic_power),
        "torque_nm": float(torque),
        "ct": float(ct),
        "cp": float(cp),
        "eta": float(efficiency),
        "activity_factor": activity_factor,
    }


def solve_power(
    cl_target: float,
    radius_fraction: np.ndarray,
    inputs: MELInputs,
    polar: PolarDatabase,
    zeta_initial: float = 0.1,
) -> dict[str, object]:
    """Close the MEL design on the requested shaft power."""

    xi = np.asarray(radius_fraction, dtype=float)
    zeta = zeta_initial
    history: list[dict[str, float]] = []
    for _ in range(60):
        state = mel_state(zeta, cl_target, xi, inputs, polar)
        gradient = np.gradient(xi)
        j1 = float(np.sum(np.asarray(state["j1p"]) * gradient))
        j2 = max(
            float(np.sum(np.asarray(state["j2p"]) * gradient)),
            np.finfo(float).eps,
        )
        zeta_new = (
            -j1
            + np.sqrt(j1**2 + 4.0 * j2 * inputs.power_coefficient_target)
        ) / (2.0 * j2)
        history.append(
            {
                "iteration": len(history) + 1,
                "zeta_in": float(zeta),
                "power_w_at_zeta_in": float(state["power_w"]),
                "zeta_update": float(zeta_new),
            }
        )
        if (
            abs(zeta_new - zeta) < 1e-3
            or abs(zeta_new - zeta) / max(abs(zeta_new), 1.0) < 5e-4
        ):
            zeta = float(zeta_new)
            break
        zeta = float(zeta_new)
    final = mel_state(zeta, cl_target, xi, inputs, polar)
    final["iteration_history"] = history
    return final


def design_sweep(
    inputs: MELInputs,
    polar: PolarDatabase,
    radius_fraction: np.ndarray,
    *,
    cl_values: np.ndarray | None = None,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    """Sweep design lift coefficient and select maximum feasible thrust."""

    values = (
        np.arange(0.50, 1.4001, 0.02)
        if cl_values is None
        else np.asarray(cl_values, dtype=float)
    )
    records: list[dict[str, object]] = []
    states: list[dict[str, object]] = []
    for cl_target in values:
        state = solve_power(float(cl_target), radius_fraction, inputs, polar)
        states.append(state)
        records.append(
            {
                "alpha_bias_deg": inputs.alpha_bias_deg,
                "cl_target": float(cl_target),
                "zeta": state["zeta"],
                "thrust_n": state["thrust_n"],
                "power_w": state["power_w"],
                "ct": state["ct"],
                "cp": state["cp"],
                "eta": state["eta"],
                "activity_factor": state["activity_factor"],
                "mean_chord_mm": 1000.0 * float(np.mean(state["c"])),
                "max_chord_mm": 1000.0 * float(np.max(state["c"])),
                "mean_alpha_deg": float(np.mean(state["alpha"])),
                "min_alpha_deg": float(np.min(state["alpha"])),
                "max_alpha_deg": float(np.max(state["alpha"])),
            }
        )
    feasible = [
        index
        for index, record in enumerate(records)
        if float(record["activity_factor"]) >= inputs.activity_factor_min
    ]
    if feasible:
        selected = max(feasible, key=lambda index: float(records[index]["thrust_n"]))
    else:
        selected = int(np.argmax([float(record["activity_factor"]) for record in records]))
    return records, states[selected]


def candidate_geometry(state: dict[str, object], inputs: MELInputs) -> Geometry:
    """Convert a MEL state into the common geometry representation."""

    xi = np.asarray(state["xi"], dtype=float)
    chord_m = np.asarray(state["c"], dtype=float)
    beta_deg = np.asarray(state["beta"], dtype=float)
    radius_in = inputs.radius_m / 0.0254
    station_in = xi * radius_in
    chord_in = chord_m / 0.0254
    pitch_in = 2.0 * np.pi * station_in * np.tan(np.radians(beta_deg))
    return Geometry(
        name="MEL_Cl076",
        station_in=station_in,
        chord_in=chord_in,
        pitch_quoted_in=pitch_in,
        pitch_le_te_in=pitch_in,
        thickness_ratio=np.full_like(xi, 0.117),
        radius_in=radius_in,
        blades=inputs.blades,
        airfoils=("CLARK-Y", "CLARK-Y"),
        path=Path("generated_from_MEL"),
    )
