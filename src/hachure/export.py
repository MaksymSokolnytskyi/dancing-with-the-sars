"""Exporters: SVG (plotter/print), PNG preview, GeoJSON (web map).

Multiple inks = multiple Layers (e.g. grey terrain under red subsidence),
the way engraved maps printed rock, ice and water from separate plates.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from rasterio.warp import transform as warp_transform

from .trace import Stroke

PAPER = "#f3eedf"
INK = "#1b1a17"


@dataclass
class Layer:
    strokes: list[Stroke]
    ink: str = INK
    w_min: float = 0.12
    w_max: float = 1.6
    name: str = "hachures"


def _width(dark: float, w_min: float, w_max: float) -> float:
    return w_min + (w_max - w_min) * dark


def _as_layers(obj, w_min, w_max, ink) -> list[Layer]:
    if obj and isinstance(obj[0], Layer):
        return obj
    return [Layer(obj, ink=ink, w_min=w_min, w_max=w_max)]


def to_svg(layers, shape, path: str | Path, *, scale: float = 2.0,
           w_min: float = 0.12, w_max: float = 1.6, classes: int = 10,
           ink: str = INK, paper: str = PAPER, simplify: int = 2) -> Path:
    """`layers`: list[Layer] or plain list[Stroke]. One <path> per width class per layer."""
    h, w = shape
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w * scale:.0f} {h * scale:.0f}" '
             f'width="{w * scale:.0f}" height="{h * scale:.0f}">',
             f'<rect width="100%" height="100%" fill="{paper}"/>']
    for layer in _as_layers(layers, w_min, w_max, ink):
        buckets: dict[int, list[str]] = {}
        for s in layer.strokes:
            c = min(int(s.dark * classes), classes - 1)
            p = s.pts[::simplify] if len(s.pts) > 2 * simplify else s.pts
            if len(p) < 2:
                p = s.pts[[0, -1]]
            d = "M" + "L".join(f"{x * scale:.1f},{y * scale:.1f}" for x, y in p)
            buckets.setdefault(c, []).append(d)
        parts.append(f'<g id="{layer.name}" fill="none" stroke="{layer.ink}" stroke-linecap="round">')
        for c in sorted(buckets):
            sw = _width((c + 0.5) / classes, layer.w_min, layer.w_max)
            parts.append(f'<path stroke-width="{sw:.2f}" d="{"".join(buckets[c])}"/>')
        parts.append("</g>")
    parts.append("</svg>")
    path = Path(path)
    path.write_text("\n".join(parts), encoding="utf-8")
    return path


def to_png(layers, shape, path: str | Path, *, scale: float = 2.0,
           w_min: float = 0.12, w_max: float = 1.6, ink: str = INK, paper: str = PAPER,
           dpi: int = 100) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    h, w = shape
    fig = plt.figure(figsize=(w * scale / dpi, h * scale / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor(paper)
    fig.patch.set_facecolor(paper)
    for layer in _as_layers(layers, w_min, w_max, ink):
        if not layer.strokes:
            continue
        # widths are in output px; matplotlib wants points (1pt = dpi/72 px)
        widths = [_width(s.dark, layer.w_min, layer.w_max) * 72 / dpi for s in layer.strokes]
        ax.add_collection(LineCollection([s.pts for s in layer.strokes], linewidths=widths,
                                         colors=layer.ink, capstyle="round"))
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")
    fig.savefig(path, dpi=dpi, facecolor=paper)
    plt.close(fig)
    return Path(path)


def to_geojson(strokes: list[Stroke], transform, crs, path: str | Path,
               simplify: int = 2, extra: dict | None = None) -> Path:
    """Pixel (col,row) -> projected via affine -> WGS84 lon/lat."""
    feats = []
    for s in strokes:
        p = s.pts[::simplify] if len(s.pts) > 2 * simplify else s.pts
        xs = transform.c + (p[:, 0] + 0.5) * transform.a
        ys = transform.f + (p[:, 1] + 0.5) * transform.e
        lon, lat = warp_transform(crs, "EPSG:4326", xs.tolist(), ys.tolist())
        feats.append({
            "type": "Feature",
            "properties": {"dark": round(s.dark, 3), **(extra or {})},
            "geometry": {"type": "LineString",
                         "coordinates": [[round(a, 6), round(b, 6)] for a, b in zip(lon, lat)]},
        })
    path = Path(path)
    path.write_text(json.dumps({"type": "FeatureCollection", "features": feats}), encoding="utf-8")
    return path
