"""Draw the Edelweiss emblem as an engraving (SVG, ink lines only) -> web/edelweiss.svg.

    uv run python scripts/draw_edelweiss.py

The same hand as the plates: a form is shaped by strokes, not by fills. Each woolly bract gets an outline, a midrib and
lengthwise hatching that crowds towards its shadow side (light from the upper left); the flower heads in the middle are
small rings with a stipple; the stem and leaves are hatched the same way. Deterministic (fixed seed), so the file is
stable in git.
"""
from __future__ import annotations

import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INK = "#3a2a1c"
PAPER = "#ece2c9"
LIGHT = (-0.7, -0.7)                       # light comes from the upper left (screen coordinates, y down)
rnd = random.Random(7)


def f(x):
    return f"{x:.1f}"


def path(pts, close=False):
    return "M" + " L".join(f"{f(x)},{f(y)}" for x, y in pts) + (" Z" if close else "")


def bract(cx, cy, ang, length, width, bend):
    """Lanceolate woolly bract from (cx, cy) along angle `ang` (radians). Returns svg elements."""
    ux, uy = math.cos(ang), math.sin(ang)
    nx, ny = -uy, ux

    def half(t):                          # half-width along the axis: narrow base, broad belly, pointed tip
        return width * (t ** 0.45) * ((1 - t) ** 0.6) * 1.55

    def at(t, u):                         # point at axis fraction t, across offset u in [-1, 1] of the half-width
        b = bend * math.sin(math.pi * t) * length      # a gentle curve of the whole bract
        x = cx + ux * length * t + nx * (u * half(t) + b)
        y = cy + uy * length * t + ny * (u * half(t) + b)
        return x, y

    ts = [i / 40 for i in range(41)]
    outline = [at(t, 1) for t in ts] + [at(t, -1) for t in reversed(ts)]
    # woolly edge: a slight random wobble on the outline
    outline = [(x + rnd.uniform(-0.6, 0.6), y + rnd.uniform(-0.6, 0.6)) for x, y in outline]
    els = [f'<path d="{path(outline, True)}" fill="{PAPER}" stroke="{INK}" stroke-width="1.1" stroke-linejoin="round"/>']
    # shadow side: the half whose normal faces away from the light
    shade = 1 if (nx * LIGHT[0] + ny * LIGHT[1]) < 0 else -1
    for k in range(1, 12):
        u = k / 12
        for side in (1, -1):
            dark = side == shade
            if not dark and k % 3:        # the lit half gets every third line only
                continue
            t0, t1 = 0.07 + 0.05 * u, 0.96 - 0.5 * u ** 1.6
            pts = [at(t0 + (t1 - t0) * i / 24, side * u * (0.97 + rnd.uniform(-0.02, 0.02))) for i in range(25)]
            w = (0.55 if dark else 0.35) * (1.15 - 0.6 * u)
            els.append(f'<path d="{path(pts)}" fill="none" stroke="{INK}" stroke-width="{w:.2f}" stroke-linecap="round"/>')
    mid = [at(0.04 + 0.9 * i / 20, 0) for i in range(21)]
    els.append(f'<path d="{path(mid)}" fill="none" stroke="{INK}" stroke-width="0.8" stroke-linecap="round"/>')
    # wool: short ticks off the outline near the tip
    for _ in range(14):
        t = rnd.uniform(0.35, 0.9)
        s = rnd.choice((1, -1))
        x, y = at(t, s)
        ex, ey = x + nx * s * rnd.uniform(1.5, 3.2), y + ny * s * rnd.uniform(1.5, 3.2)
        els.append(f'<path d="M{f(x)},{f(y)} L{f(ex)},{f(ey)}" stroke="{INK}" stroke-width="0.45"/>')
    return els


def leaf(x0, y0, ang, length, width):
    ux, uy = math.cos(ang), math.sin(ang)
    nx, ny = -uy, ux

    def at(t, u):
        h = width * math.sin(math.pi * t ** 0.8) * (1 - 0.3 * t)
        return x0 + ux * length * t + nx * u * h, y0 + uy * length * t + ny * u * h

    ts = [i / 30 for i in range(31)]
    els = [f'<path d="{path([at(t, 1) for t in ts] + [at(t, -1) for t in reversed(ts)], True)}" fill="{PAPER}" '
           f'stroke="{INK}" stroke-width="1"/>']
    els.append(f'<path d="{path([at(t, 0) for t in ts])}" fill="none" stroke="{INK}" stroke-width="0.7"/>')
    for i in range(3, 28, 2):            # oblique hatching from the vein to the lower edge
        t = i / 30
        a, b = at(t, 0), at(min(1, t + 0.06), -0.85)
        els.append(f'<path d="M{f(a[0])},{f(a[1])} L{f(b[0])},{f(b[1])}" stroke="{INK}" stroke-width="0.4"/>')
    return els


def main():
    els = []
    # stem: two lines, hatched on the right (shadow) side
    stem = [(0 + 6 * math.sin(i / 20 * 1.3), 18 + i * 6.2) for i in range(21)]
    els.append(f'<path d="{path(stem)}" fill="none" stroke="{INK}" stroke-width="2.2" stroke-linecap="round"/>')
    for x, y in stem[2::2]:
        els.append(f'<path d="M{f(x + 1.5)},{f(y)} L{f(x + 4.5)},{f(y + 3)}" stroke="{INK}" stroke-width="0.45"/>')
    els += leaf(stem[9][0], stem[9][1], math.radians(-58), 50, 6)
    els += leaf(stem[12][0], stem[12][1], math.radians(-128), 46, 5.5)
    els += leaf(stem[16][0], stem[16][1], math.radians(-62), 36, 5)
    # big bracts, then small ones in the gaps
    n = 9
    for i in range(n):
        a = -math.pi / 2 + 2 * math.pi * i / n + rnd.uniform(-0.08, 0.08)
        els += bract(0, 0, a, rnd.uniform(62, 92), rnd.uniform(20, 25), rnd.uniform(-0.06, 0.06))
    for i in range(n):
        a = -math.pi / 2 + 2 * math.pi * (i + 0.5) / n + rnd.uniform(-0.1, 0.1)
        els += bract(0, 0, a, rnd.uniform(40, 52), rnd.uniform(13, 16), rnd.uniform(-0.05, 0.05))
    # the flower heads: a disc of small rings, each with a dot, and stipple between them
    heads = [(0, 0)] + [(12 * math.cos(2 * math.pi * k / 5 + 0.4), 12 * math.sin(2 * math.pi * k / 5 + 0.4)) for k in range(5)]
    for hx, hy in heads:                  # each head: a paper dome with a ring, crowded with florets
        els.append(f'<circle cx="{f(hx)}" cy="{f(hy)}" r="7.2" fill="{PAPER}" stroke="{INK}" stroke-width="0.8"/>')
        for _ in range(13):
            r = 5.6 * math.sqrt(rnd.random())
            a = rnd.uniform(0, 2 * math.pi)
            x, y = hx + r * math.cos(a), hy + r * math.sin(a)
            els.append(f'<circle cx="{f(x)}" cy="{f(y)}" r="1.5" fill="{PAPER}" stroke="{INK}" stroke-width="0.6"/>')
            els.append(f'<circle cx="{f(x + 0.3)}" cy="{f(y + 0.3)}" r="0.45" fill="{INK}"/>')
    for _ in range(110):                  # stipple, denser on the shadow (lower right) side
        r = 20 * math.sqrt(rnd.random())
        a = rnd.uniform(0, 2 * math.pi)
        x, y = r * math.cos(a), r * math.sin(a)
        if rnd.random() < 0.5 + 0.4 * (x + y) / 37:
            els.append(f'<circle cx="{f(x)}" cy="{f(y)}" r="0.42" fill="{INK}"/>')
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="-115 -110 230 255" role="img" aria-label="Edelweiss">'
           f'<title>Edelweiss</title><g>{"".join(els)}</g></svg>\n')
    out = ROOT / "web/edelweiss.svg"
    out.write_text(svg, encoding="utf-8")
    print(f"{out.relative_to(ROOT)}: {len(els)} elements, {len(svg) / 1e3:.0f} kB")


if __name__ == "__main__":
    main()
