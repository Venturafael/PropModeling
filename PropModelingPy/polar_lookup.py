"""Single Clark--Y polar lookup used by BET, BEMT, and Vortex Theory.

The mature mode follows ``busca_polares_corretas.m`` where defensible:

* PCHIP interpolation in angle of attack;
* bilinear interpolation in log10(Re) and Mach;
* no extrapolation: each coordinate is clipped and every clip is counted.

Unlike the MATLAB routine, this implementation is stateless and order invariant.
The MATLAB persistent monotonic bridge has no effect for ascending J with the
shared kinematic Reynolds/Mach calculation, but it makes a lookup depend on
which point happened to run before it.  Removing that state is necessary for a
reproducible global batch.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.io import loadmat


@dataclass
class PolarStats:
    queries: int = 0
    alpha_clamped: int = 0
    reynolds_clamped_low: int = 0
    reynolds_clamped_high: int = 0
    mach_clamped_low: int = 0
    mach_clamped_high: int = 0

    def as_dict(self) -> dict[str, int]:
        return dict(vars(self))


@dataclass(frozen=True)
class PolarQuery:
    cl: float
    cd: float
    cm: float
    alpha_used: float
    reynolds_used: float
    mach_used: float
    alpha_clamped: bool
    reynolds_clamped: bool
    mach_clamped: bool


class PolarDatabase:
    """Load and query the expanded Clark--Y XFOIL database."""

    def __init__(
        self,
        mat_path: str | Path,
        mode: Literal["mature", "legacy_linear"] = "mature",
    ) -> None:
        self.path = Path(mat_path).resolve()
        self.mode = mode
        raw = loadmat(self.path, struct_as_record=False, squeeze_me=True)
        try:
            grid = raw["ans"].perfil
        except (KeyError, AttributeError) as exc:
            raise ValueError(f"Unsupported polar structure in {self.path}") from exc
        if grid.ndim != 2:
            raise ValueError(f"Expected a two-dimensional polar grid, got {grid.shape}")

        re_axis = np.array([float(grid[i, 0].Re) for i in range(grid.shape[0])])
        mach_axis = np.array([float(grid[0, j].Mach) for j in range(grid.shape[1])])
        i_order = np.argsort(re_axis)
        j_order = np.argsort(mach_axis)
        self.reynolds_axis = re_axis[i_order]
        self.mach_axis = mach_axis[j_order]
        self.grid = grid[np.ix_(i_order, j_order)]
        self._curves: dict[tuple[int, int], tuple[np.ndarray, ...]] = {}
        self._pchips: dict[tuple[int, int], PchipInterpolator] = {}
        self.stats = PolarStats()

        alpha_min = []
        alpha_max = []
        alpha_lengths = []
        for cell in self.grid.flat:
            alpha = np.asarray(cell.Alpha, dtype=float).ravel()
            alpha_min.append(float(np.min(alpha)))
            alpha_max.append(float(np.max(alpha)))
            alpha_lengths.append(alpha.size)
        self.alpha_min = float(min(alpha_min))
        self.alpha_max = float(max(alpha_max))
        self.safe_alpha_min = float(max(alpha_min))
        self.safe_alpha_max = float(min(alpha_max))
        self.alpha_lengths = tuple(sorted(set(alpha_lengths)))

    def reset_stats(self) -> None:
        self.stats = PolarStats()

    @property
    def domain(self) -> dict[str, object]:
        return {
            "reynolds": [float(self.reynolds_axis[0]), float(self.reynolds_axis[-1])],
            "mach": [float(self.mach_axis[0]), float(self.mach_axis[-1])],
            "alpha_union_deg": [self.alpha_min, self.alpha_max],
            "alpha_common_deg": [self.safe_alpha_min, self.safe_alpha_max],
            "alpha_point_counts": list(self.alpha_lengths),
            "shape": list(self.grid.shape),
            "mode": self.mode,
        }

    @staticmethod
    def _bracket(value: float, axis: np.ndarray) -> tuple[int, int, float]:
        hi = int(np.searchsorted(axis, value, side="left"))
        if hi <= 0:
            return 0, 0, 0.0
        if hi >= axis.size:
            last = axis.size - 1
            return last, last, 0.0
        lo = hi - 1
        fraction = (value - axis[lo]) / (axis[hi] - axis[lo])
        return lo, hi, float(np.clip(fraction, 0.0, 1.0))

    def _cell_arrays(self, i: int, j: int) -> tuple[np.ndarray, ...]:
        key = (i, j)
        if key not in self._curves:
            cell = self.grid[i, j]
            alpha = np.asarray(cell.Alpha, dtype=float).ravel()
            cl = np.asarray(cell.Cl, dtype=float).ravel()
            cd = np.asarray(cell.Cd, dtype=float).ravel()
            cm = np.asarray(getattr(cell, "Cm", np.zeros_like(alpha)), dtype=float).ravel()
            order = np.argsort(alpha)
            alpha, cl, cd, cm = alpha[order], cl[order], cd[order], cm[order]
            unique, unique_idx = np.unique(alpha, return_index=True)
            self._curves[key] = (unique, cl[unique_idx], cd[unique_idx], cm[unique_idx])
        return self._curves[key]

    def _cell_value(self, i: int, j: int, alpha: float) -> tuple[float, float, float, float]:
        arrays = self._cell_arrays(i, j)
        a = float(np.clip(alpha, arrays[0][0], arrays[0][-1]))
        if self.mode == "legacy_linear":
            return (
                float(np.interp(a, arrays[0], arrays[1])),
                float(np.interp(a, arrays[0], arrays[2])),
                float(np.interp(a, arrays[0], arrays[3])),
                a,
            )
        key = (i, j)
        if key not in self._pchips:
            self._pchips[key] = PchipInterpolator(
                arrays[0], np.column_stack(arrays[1:]), axis=0, extrapolate=False
            )
        values = np.asarray(self._pchips[key](a), dtype=float)
        return float(values[0]), float(values[1]), float(values[2]), a

    def lookup_many(
        self,
        alpha_deg: np.ndarray,
        reynolds: float,
        mach: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Fast Cl/Cd lookup for a vector of alpha at one Re/Mach condition.

        This is used only to locate BEMT root brackets; the accepted solution
        is queried again through :meth:`lookup` and is therefore audited.
        """
        alpha_deg = np.asarray(alpha_deg, dtype=float)
        re_clipped = float(np.clip(reynolds, self.reynolds_axis[0], self.reynolds_axis[-1]))
        ma_clipped = float(np.clip(mach, self.mach_axis[0], self.mach_axis[-1]))
        re_coord = np.log10(re_clipped) if self.mode == "mature" else re_clipped
        re_axis = np.log10(self.reynolds_axis) if self.mode == "mature" else self.reynolds_axis
        i1, i2, tr = self._bracket(re_coord, re_axis)
        j1, j2, tm = self._bracket(ma_clipped, self.mach_axis)
        cl = np.zeros_like(alpha_deg)
        cd = np.zeros_like(alpha_deg)
        total_weight = 0.0
        for i, j, weight in (
            (i1, j1, (1 - tr) * (1 - tm)),
            (i1, j2, (1 - tr) * tm),
            (i2, j1, tr * (1 - tm)),
            (i2, j2, tr * tm),
        ):
            if weight <= 0.0:
                continue
            arrays = self._cell_arrays(i, j)
            alpha = np.clip(alpha_deg, arrays[0][0], arrays[0][-1])
            if self.mode == "legacy_linear":
                cl_value = np.interp(alpha, arrays[0], arrays[1])
                cd_value = np.interp(alpha, arrays[0], arrays[2])
            else:
                key = (i, j)
                if key not in self._pchips:
                    self._pchips[key] = PchipInterpolator(
                        arrays[0], np.column_stack(arrays[1:]), axis=0, extrapolate=False
                    )
                values = np.asarray(self._pchips[key](alpha), dtype=float)
                cl_value, cd_value = values[:, 0], values[:, 1]
            cl += weight * cl_value
            cd += weight * cd_value
            total_weight += weight
        if total_weight == 0.0:  # both coordinates exactly on a lower boundary
            arrays = self._cell_arrays(i1, j1)
            alpha = np.clip(alpha_deg, arrays[0][0], arrays[0][-1])
            if self.mode == "legacy_linear":
                return np.interp(alpha, arrays[0], arrays[1]), np.interp(alpha, arrays[0], arrays[2])
            key = (i1, j1)
            if key not in self._pchips:
                self._pchips[key] = PchipInterpolator(
                    arrays[0], np.column_stack(arrays[1:]), axis=0, extrapolate=False
                )
            values = np.asarray(self._pchips[key](alpha), dtype=float)
            return values[:, 0], values[:, 1]
        return cl / total_weight, cd / total_weight

    def lookup(
        self,
        alpha_deg: float,
        reynolds: float,
        mach: float,
        *,
        record: bool = True,
    ) -> PolarQuery:
        """Return Cl/Cd/Cm without extrapolating any database coordinate."""
        if not np.all(np.isfinite([alpha_deg, reynolds, mach])):
            return PolarQuery(*(float("nan"),) * 6, False, False, False)

        re_min, re_max = self.reynolds_axis[[0, -1]]
        ma_min, ma_max = self.mach_axis[[0, -1]]
        re_clipped = float(np.clip(reynolds, re_min, re_max))
        ma_clipped = float(np.clip(mach, ma_min, ma_max))

        re_coord = np.log10(re_clipped) if self.mode == "mature" else re_clipped
        re_axis = np.log10(self.reynolds_axis) if self.mode == "mature" else self.reynolds_axis
        i1, i2, tr = self._bracket(re_coord, re_axis)
        j1, j2, tm = self._bracket(ma_clipped, self.mach_axis)

        corners = []
        weights = []
        for i, j, weight in (
            (i1, j1, (1 - tr) * (1 - tm)),
            (i1, j2, (1 - tr) * tm),
            (i2, j1, tr * (1 - tm)),
            (i2, j2, tr * tm),
        ):
            if weight > 0.0 or (i1 == i2 and j1 == j2 and not corners):
                corners.append(self._cell_value(i, j, alpha_deg))
                weights.append(weight if weight > 0.0 else 1.0)
        w = np.asarray(weights, dtype=float)
        w /= np.sum(w)
        values = np.asarray(corners, dtype=float)
        cl, cd, cm, alpha_used = np.sum(values * w[:, None], axis=0)

        alpha_clamped = not np.allclose(values[:, 3], alpha_deg, atol=1e-12, rtol=0)
        re_clamped = not np.isclose(re_clipped, reynolds)
        ma_clamped = not np.isclose(ma_clipped, mach)
        if record:
            self.stats.queries += 1
            self.stats.alpha_clamped += int(alpha_clamped)
            self.stats.reynolds_clamped_low += int(reynolds < re_min)
            self.stats.reynolds_clamped_high += int(reynolds > re_max)
            self.stats.mach_clamped_low += int(mach < ma_min)
            self.stats.mach_clamped_high += int(mach > ma_max)
        return PolarQuery(
            float(cl), float(cd), float(cm), float(alpha_used), re_clipped, ma_clipped,
            alpha_clamped, re_clamped, ma_clamped,
        )
