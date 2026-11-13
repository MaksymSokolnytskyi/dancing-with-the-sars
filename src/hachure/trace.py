"""Hachure placement: strokes between successive contours, along the line of steepest descent.

For each contour band [L - I, L]:
  1. seed candidates on the band's mid-contour (priority), then on quarter contours (gap fill);
  2. trace every candidate up to L and down to L - I with vectorized RK2 on the unit
     downhill field;
  3. accept strokes greedily against an occupancy grid (Jobard & Lefer 1997 idea):
     a seed must be >= `spacing` from accepted strokes, a stroke is cut where it comes
     closer than `min_gap` (convergence in hollows).
Rows of strokes per band give the engraved "terrace" rhythm of Lehmann/Dufour maps.
"""
from __future__ import annotations

from dataclasses import dataclass

import contourpy
import numpy as np
from scipy.ndimage import map_coordinates

from .field import Field


@dataclass
class Stroke:
    pts: np.ndarray   # (k, 2) array of (col, row), from top (uphill) to bottom
    dark: float       # 0..1


def _sample(arr: np.ndarray, pts: np.ndarray) -> np.ndarray:
    return map_coordinates(arr, [pts[:, 1], pts[:, 0]], order=1, mode="nearest")


def _resample_line(line: np.ndarray, spacing: float, offset: float) -> np.ndarray:
    seg = np.hypot(*np.diff(line, axis=0).T)
    cum = np.concatenate([[0], np.cumsum(seg)])
    if cum[-1] < spacing:
        return np.empty((0, 2))
    t = np.arange(offset % spacing, cum[-1], spacing)
    return np.column_stack([np.interp(t, cum, line[:, 0]), np.interp(t, cum, line[:, 1])])


def _trace(f: Field, seeds: np.ndarray, sign: float, z_lo: float, z_hi: float,
           step: float, max_steps: int, min_slope: float) -> np.ndarray:
    """Vectorized RK2. Returns (max_steps+1, n, 2); finished paths repeat their last point."""
    h, w = f.z.shape
    pos = seeds.copy()
    active = np.ones(len(seeds), bool)
    out = [pos.copy()]
    for _ in range(max_steps):
        if not active.any():
            break
        p = pos[active]
        v1 = np.column_stack([_sample(f.ux, p), _sample(f.uy, p)])
        mid = p + 0.5 * step * sign * v1
        v2 = np.column_stack([_sample(f.ux, mid), _sample(f.uy, mid)])
        n2 = np.hypot(v2[:, 0], v2[:, 1])[:, None]
        new = p + step * sign * v2 / np.where(n2 > 1e-9, n2, 1)
        zz = _sample(f.z, new)
        sl = _sample(f.slope, new)
        ok = ((zz >= z_lo) & (zz <= z_hi) & (sl >= min_slope) & (n2[:, 0] > 1e-6)
              & (new[:, 0] >= 0) & (new[:, 0] <= w - 1) & (new[:, 1] >= 0) & (new[:, 1] <= h - 1))
        idx = np.flatnonzero(active)
        pos[idx[ok]] = new[ok]
        active[idx[~ok]] = False
        out.append(pos.copy())
    return np.stack(out)


def _dedupe(path: np.ndarray) -> np.ndarray:
    keep = np.concatenate([[True], np.any(np.diff(path, axis=0) != 0, axis=1)])
    return path[keep]


class _Occupancy:
    def __init__(self, shape, cell: float):
        self.cell = cell
        self.grid = np.full((int(shape[0] / cell) + 3, int(shape[1] / cell) + 3), -1, np.int32)

    def _ij(self, p):
        return int(p[1] / self.cell) + 1, int(p[0] / self.cell) + 1

    def free(self, p, radius_cells: int, owner: int) -> bool:
        i, j = self._ij(p)
        r = radius_cells
        win = self.grid[max(i - r, 0):i + r + 1, max(j - r, 0):j + r + 1]
        return not np.any((win >= 0) & (win != owner))

    def mark(self, pts, owner: int):
        for p in pts:
            i, j = self._ij(p)
            self.grid[i, j] = owner


def hachures(f: Field, *, interval_m: float = 50.0, spacing_px: float = 1.6,
             min_gap_px: float | None = None, step_px: float = 0.4,
             max_len_px: float = 12.0, min_len_px: float = 0.8,
             min_slope_deg: float = 2.0, min_dark: float = 0.04) -> list[Stroke]:
    """Generate hachure strokes. Units are grid pixels except interval (field units)."""
    min_gap = min_gap_px if min_gap_px is not None else 0.4 * spacing_px
    cell = min_gap / 2
    gap_r = 2

    def seed_radius(d: float) -> int:
        # Lehmann: dark slopes -> dense strokes, light slopes -> sparse
        return max(1, int(round(spacing_px * (1.5 - 0.8 * d) / cell)))
    max_steps = int(max_len_px / step_px)

    gen = contourpy.contour_generator(z=f.z, line_type=contourpy.LineType.Separate)
    z_min, z_max = float(np.nanmin(f.z)), float(np.nanmax(f.z))
    levels = np.arange(np.ceil(z_min / interval_m) * interval_m, z_max + interval_m, interval_m)

    strokes: list[Stroke] = []
    owner = 0
    for top in levels:
        lo, hi = top - interval_m, top
        occ = _Occupancy(f.z.shape, cell)   # per band: neighbouring bands may touch
        # mid contour first, then quarter contours to fill gaps on spurs / in hollows
        for frac in (0.5, 0.25, 0.75):
            lvl = hi - frac * interval_m
            lines = gen.lines(lvl)
            if not lines:
                continue
            seeds = np.concatenate([_resample_line(ln, 0.5 * spacing_px,
                                                   offset=frac * spacing_px)
                                    for ln in lines if len(ln) > 1] or [np.empty((0, 2))])
            if len(seeds) == 0:
                continue
            dark = _sample(f.dark, seeds)
            slope = _sample(f.slope, seeds)
            keep = (dark >= min_dark) & (slope >= min_slope_deg)
            seeds, dark = seeds[keep], dark[keep]
            if len(seeds) == 0:
                continue
            up = _trace(f, seeds, -1.0, lo, hi, step_px, max_steps // 2, min_slope_deg)
            down = _trace(f, seeds, +1.0, lo, hi, step_px, max_steps // 2, min_slope_deg)
            for k in range(len(seeds)):
                if not occ.free(seeds[k], seed_radius(dark[k]), -2):
                    continue
                # grow from seed outward on both sides, cut at first conflict
                halves = []
                for path in (up[:, k], down[:, k]):
                    path = _dedupe(path)
                    cut = len(path)
                    for m in range(1, len(path)):
                        if not occ.free(path[m], gap_r, owner):
                            cut = m
                            break
                    halves.append(path[:cut])
                pts = np.concatenate([halves[0][::-1], halves[1][1:]])
                length = np.hypot(*np.diff(pts, axis=0).T).sum() if len(pts) > 1 else 0
                if length < min_len_px:
                    continue
                occ.mark(pts, owner)
                strokes.append(Stroke(pts=pts, dark=float(dark[k])))
                owner += 1
    return strokes
