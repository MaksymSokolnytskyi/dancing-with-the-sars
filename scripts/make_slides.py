"""Build the 7-slide deck (web/demo/edelweiss-slides.pdf) from the demo cards and plate screenshots.

    uv run python scripts/make_slides.py <montage_dir>

montage_dir: the same folder scripts/make_demo.py reads (title.jpg, end.jpg, m0*..m4*). Numbers on the slides are the
plates' own numbers, with their errors and the failed tests (see notes/submission-2026.md).
"""
from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
W, H = 1280, 720
PAPER, INK, MUTED, RED = (239, 230, 210), (40, 28, 18), (110, 92, 70), (150, 52, 36)


def font(size, italic=False):
    return ImageFont.truetype("C:/Windows/Fonts/georgiai.ttf" if italic else "C:/Windows/Fonts/georgia.ttf", size)


def shot(p):
    return Image.open(p).convert("RGB").resize((W, H), Image.LANCZOS)


def plate(p, title, lines):
    """A plate screenshot with a paper band at the foot: the claim and its numbers."""
    im = shot(p)
    d = ImageDraw.Draw(im)
    band = 46 + 34 * len(lines)
    d.rectangle([0, H - band, W, H], fill=PAPER)
    d.line([0, H - band, W, H - band], fill=INK, width=1)
    d.text((40, H - band + 12), title, font=font(30), fill=INK)
    for i, t in enumerate(lines):
        d.text((40, H - band + 52 + 34 * i), t, font=font(23, i == len(lines) - 1), fill=INK if i < len(lines) - 1 else MUTED)
    return im


def main():
    m = Path(sys.argv[1])
    pick = lambda k, i: sorted(m.glob(f"m{k}*.jpg"))[i]  # noqa: E731
    duf = Image.new("RGB", (W, H), PAPER)
    duf.paste(Image.open(ROOT / "notes" / "img" / "dufour_jungfrau_z13.png").convert("RGB").resize((560, 560)), (90, 80))
    d = ImageDraw.Draw(duf)
    for text, y, size, it, col in [("The idea", 150, 30, True, MUTED), ("Radar data drawn like", 200, 44, False, INK),
                                   ("Dufour's engraved maps", 256, 44, False, INK), ("hachures, ink washes, parcels and signs,", 360, 26, False, INK),
                                   ("not rainbow fringes: readable in seconds", 396, 26, False, INK),
                                   ("NISAR L-band radar, every 12 days, 2025–2026", 470, 24, True, MUTED),
                                   ("Dufour map: © swisstopo", 650, 18, False, MUTED)]:
        d.text((710, y), text, font=font(size, it), fill=col)
    slides = [
        shot(m / "title.jpg"),
        duf,
        plate(pick(0, -2), "Reykjanes, Iceland: the dome over the magma rises",
              ["up to 23.5 ± 1.0 cm a year; the towns on its flank tilt",
               "only one of four GPS stations agrees, and the plate says so"]),
        plate(pick(1, -2), "Mexico City sinks into its old lakebed",
              ["up to about 30 cm a year, along the shore of the vanished lake; GPS −27.0 cm/yr, radar −27.8",
               "the town map first, then the moving ground"]),
        plate(pick(2, -1), "Ostrava–Karviná: coal under two countries",
              ["22 sinking patches, 12 pass a rule fixed before drawing; ~1,140 buildings on the Polish ones",
               "early data; the check against Europe's ground-motion map came out only partial"]),
        plate(pick(4, 1), "Fields through the year: Central Valley and Wisconsin",
              ["the farm calendar in radar coherence; corn 3–8 dB brighter than soybean; the harvest is seen",
               "failed tests (rain, lake ice, frost heave, swamp forest, dry afternoons) stay on the plates"]),
        shot(m / "end.jpg"),
    ]
    out = ROOT / "web" / "demo" / "edelweiss-slides.pdf"
    slides[0].save(out, save_all=True, append_images=slides[1:], resolution=96)
    print(f"{out.relative_to(ROOT)}: {len(slides)} slides, {out.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
