"""Readers and auditable discovery for APC PE0/PER3 data files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np


@dataclass(frozen=True)
class Geometry:
    name: str
    station_in: np.ndarray
    chord_in: np.ndarray
    pitch_quoted_in: np.ndarray
    pitch_le_te_in: np.ndarray
    thickness_ratio: np.ndarray
    radius_in: float
    blades: int
    airfoils: tuple[str, ...]
    path: Path

    @property
    def diameter_in(self) -> float:
        return 2.0 * self.radius_in

    @property
    def nominal_pitch_in(self) -> float:
        finite = self.pitch_quoted_in[np.isfinite(self.pitch_quoted_in)]
        return float(np.median(finite)) if finite.size else float("nan")

    @property
    def x(self) -> np.ndarray:
        return self.station_in / self.radius_in


@dataclass(frozen=True)
class PerformanceBlock:
    rpm: int
    j: np.ndarray
    eta: np.ndarray
    ct: np.ndarray
    cp: np.ndarray
    mach_apc: np.ndarray
    reynolds_apc: np.ndarray


@dataclass(frozen=True)
class PropellerCase:
    name: str
    geometry_path: Path
    performance_path: Path | None
    airfoils: tuple[str, ...]
    status: str
    reason: str


def _float_rows(lines: list[str], expected: int) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in lines:
        try:
            values = [float(token) for token in line.replace(",", ".").split()]
        except ValueError:
            continue
        if len(values) == expected:
            rows.append(values)
    return rows


def load_geometry(path: str | Path) -> Geometry:
    path = Path(path).resolve()
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    header = next((i for i, line in enumerate(lines) if "STATION" in line and "CHORD" in line), None)
    radius_line = next((i for i, line in enumerate(lines) if "RADIUS:" in line), None)
    if header is None or radius_line is None:
        raise ValueError(f"Geometry table not found in {path}")
    rows = _float_rows(lines[header + 1 : radius_line], 13)
    if not rows:
        raise ValueError(f"No 13-column geometry rows in {path}")
    data = np.asarray(rows, dtype=float)
    radius_match = re.search(r"RADIUS:\s*([0-9Ee+.-]+)", lines[radius_line])
    blades_match = re.search(r"BLADES:\s*(\d+)", text)
    if not radius_match or not blades_match:
        raise ValueError(f"RADIUS/BLADES metadata missing in {path}")
    airfoils = tuple(
        value.strip().upper()
        for value in re.findall(r"AIRFOIL\d+:\s*[^,]*,\s*([^\r\n(]+)", text, re.I)
    )
    name = path.name.removesuffix("-PERF.PE0")
    return Geometry(
        name=name,
        station_in=data[:, 0],
        chord_in=data[:, 1],
        pitch_quoted_in=data[:, 2],
        pitch_le_te_in=data[:, 3],
        thickness_ratio=data[:, 6],
        radius_in=float(radius_match.group(1)),
        blades=int(blades_match.group(1)),
        airfoils=airfoils,
        path=path,
    )


def load_performance(path: str | Path) -> dict[int, PerformanceBlock]:
    path = Path(path).resolve()
    lines = path.read_text(errors="replace").splitlines()
    starts: list[tuple[int, int]] = []
    for index, line in enumerate(lines):
        match = re.search(r"PROP RPM\s*=\s*(\d+)", line)
        if match:
            starts.append((index, int(match.group(1))))
    blocks: dict[int, PerformanceBlock] = {}
    for pos, (start, rpm) in enumerate(starts):
        end = starts[pos + 1][0] if pos + 1 < len(starts) else len(lines)
        numeric: list[list[float]] = []
        for line in lines[start + 1 : end]:
            try:
                values = [float(token) for token in line.split()]
            except ValueError:
                continue
            if len(values) >= 15:
                numeric.append(values)
        if not numeric:
            continue
        data = np.asarray(numeric, dtype=float)
        blocks[rpm] = PerformanceBlock(
            rpm=rpm,
            j=data[:, 1], eta=data[:, 2], ct=data[:, 3], cp=data[:, 4],
            mach_apc=data[:, 12], reynolds_apc=data[:, 13],
        )
    if not blocks:
        raise ValueError(f"No performance blocks found in {path}")
    return blocks


def select_representative_rpm(
    blocks: dict[int, PerformanceBlock],
    diameter_in: float,
    speed_of_sound: float = 343.0,
    max_tip_mach: float = 0.86,
) -> tuple[int, float, str]:
    """Select the highest available RPM inside the adopted tip-Mach ceiling."""
    diameter_m = diameter_in * 0.0254
    candidates = []
    for rpm in sorted(blocks):
        mach = np.pi * diameter_m * rpm / (60.0 * speed_of_sound)
        if mach <= max_tip_mach:
            candidates.append((rpm, mach))
    if candidates:
        rpm, mach = candidates[-1]
        return rpm, float(mach), "highest_rpm_with_static_tip_mach_le_0.86"
    rpm = min(blocks)
    mach = np.pi * diameter_m * rpm / (60.0 * speed_of_sound)
    return rpm, float(mach), "minimum_rpm_used_but_static_tip_mach_above_0.86"


def _apc_sources(root: Path) -> list[Path]:
    """Return existing APC roots in deterministic priority order.

    The first entry is the clean public-repository layout.  Historical layouts
    remain supported so the same reader can audit the unmodified TCC tree.
    """

    candidates = [root / "APC", root / "APC (2)", root]
    sources: list[Path] = []
    for candidate in candidates:
        geometry_exists = any(
            (candidate / directory).is_dir()
            for directory in ("geometry", "Geometria")
        )
        performance_exists = any(
            (candidate / directory).is_dir()
            for directory in ("performance", "Performance")
        )
        if geometry_exists and performance_exists and candidate not in sources:
            sources.append(candidate)
    return sources


def discover_apc_cases(
    root: str | Path,
    *,
    include_marine: bool = False,
) -> list[PropellerCase]:
    """Inventory and pair APC PE0 geometry with PER3 performance files.

    Parameters
    ----------
    root:
        Either the repository ``data`` directory, the ``data/APC`` directory,
        or the historical TCC root containing ``APC`` / ``APC (2)`` mirrors.
    include_marine:
        Also match files below ``PERFILES2-MARINE``.  It defaults to ``False``
        to reproduce the 224-propeller aircraft-data study.  The marine files
        are still retained in the repository and can be audited explicitly.

    Notes
    -----
    All downloaded APC files may remain in the dataset.  Inclusion in a model
    batch is stricter: the PE0 file must explicitly identify Clark--Y at every
    declared airfoil station and must have a matching PER3 file.
    """
    root = Path(root).resolve()
    sources = _apc_sources(root)
    if not sources:
        raise FileNotFoundError(
            f"No APC geometry/performance directory pair found below {root}"
        )
    geometries: dict[str, Path] = {}
    performances: dict[str, Path] = {}
    for source in sources:
        geometry_directory = next(
            directory
            for directory in (source / "geometry", source / "Geometria")
            if directory.is_dir()
        )
        performance_directory = next(
            directory
            for directory in (source / "performance", source / "Performance")
            if directory.is_dir()
        )
        for path in sorted(geometry_directory.iterdir()):
            if not path.is_file() or path.suffix.upper() != ".PE0":
                continue
            name = re.sub(r"-PERF\.PE0$", "", path.name, flags=re.I)
            geometries.setdefault(name, path.resolve())
        search_root = performance_directory
        standard_directory = performance_directory / "PERFILES2"
        if standard_directory.is_dir() and not include_marine:
            search_root = standard_directory
        for path in sorted(search_root.rglob("*")):
            if not path.is_file():
                continue
            if not include_marine and "PERFILES2-MARINE" in path.parts:
                continue
            match = re.fullmatch(r"PER3_(.+)\.dat", path.name, flags=re.I)
            if match:
                performances.setdefault(match.group(1), path.resolve())

    cases: list[PropellerCase] = []
    for name, geometry_path in sorted(geometries.items()):
        try:
            geometry = load_geometry(geometry_path)
            airfoils = geometry.airfoils
        except Exception as exc:  # inventory must retain unreadable cases
            cases.append(PropellerCase(name, geometry_path, performances.get(name), (), "excluded", f"geometry_parse_error: {exc}"))
            continue
        performance_path = performances.get(name)
        if not airfoils:
            status, reason = "uncertain", "airfoil metadata absent"
        elif not all(airfoil == "CLARK-Y" for airfoil in airfoils):
            status, reason = "excluded", "PE0 airfoil metadata is not exclusively CLARK-Y"
        elif performance_path is None:
            status, reason = "excluded", "Clark-Y confirmed but matching PER3 performance file absent"
        else:
            status, reason = "included", "both PE0 airfoil stations explicitly identify CLARK-Y"
        cases.append(PropellerCase(name, geometry_path, performance_path, airfoils, status, reason))
    return cases

