"""Build the 30-second demo video (web/demo/edelweiss-demo.mp4) from plate screenshots, a Dufour sample and a narration.

    uv run python scripts/make_demo.py <montage_dir> <voice_dir>

montage_dir: title.jpg and end.jpg (web/demo/cards.html), then m0*..m4* = Reykjanes, Mexico City, Ostrava-Karvina,
Central Valley, Wisconsin, screenshots of the plates without the frame, each opening on the town map (OSM) and then
the data playing through the year, with a zoom at the end;
voice_dir: c1..c7.wav, the narration lines in order (placed at the times in VOICE_AT).
The narration is read by the Windows speech synthesizer (System.Speech); the script was drafted with AI (Claude),
which the end card says. The Dufour sample is (c) swisstopo (free use with credit).
"""
from __future__ import annotations

import subprocess
import sys
import wave
from pathlib import Path

import imageio_ffmpeg
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
W, H, FPS, DUR = 1280, 720, 24, 30.0
PAPER, INK, MUTED, RED = (239, 230, 210), (40, 28, 18), (110, 92, 70), (150, 52, 36)
FONT = "C:/Windows/Fonts/georgia.ttf"
FONT_I = "C:/Windows/Fonts/georgiai.ttf"
# when each narration line starts (s): title and Dufour, Reykjanes, then a few words per plate, and the close
VOICE_AT = [0.2, 4.9, 11.3, 17.0, 21.0, 24.6, 27.1]
# seconds per part: title, Dufour, then the five plates, then the end card
TITLE_S, DUFOUR_S, PLATE_S, END_S = 2.2, 2.6, (5.6, 6.0, 4.2, 3.6, 3.2), 2.6


def font(size, italic=False):
    return ImageFont.truetype(FONT_I if italic else FONT, size)


def card(lines):
    im = Image.new("RGB", (W, H), PAPER)
    d = ImageDraw.Draw(im)
    for text, y, size, italic, col in lines:
        f = font(size, italic)
        d.text(((W - d.textlength(text, font=f)) / 2, y), text, font=f, fill=col)
    return im


def spaced(t):
    return " ".join(t)


def main():
    duf = card([])
    dimg = Image.open(ROOT / "notes" / "img" / "dufour_jungfrau_z13.png").convert("RGB").resize((560, 560))
    duf.paste(dimg, (90, 80))
    d = ImageDraw.Draw(duf)
    for text, y, size, italic, col in [("Inspired by", 210, 34, True, MUTED), ("Dufour's engraved", 260, 46, False, INK),
                                        ("maps of the Alps", 318, 46, False, INK), ("1845â€“1865", 380, 30, True, MUTED),
                                        ("radar data drawn as ink hachures,", 470, 26, False, INK),
                                        ("not rainbow fringes", 506, 26, False, INK), ("Dufour map: Â© swisstopo", 650, 18, False, MUTED)]:
        d.text((720, y), text, font=font(size, italic), fill=col)
    mdir = Path(sys.argv[1])
    shot = lambda p: Image.open(p).convert("RGB").resize((W, H), Image.LANCZOS)  # noqa: E731
    title, end = shot(mdir / "title.jpg"), shot(mdir / "end.jpg")
    seq = [(title, TITLE_S), (duf, DUFOUR_S)]
    for k, secs in enumerate(PLATE_S):                        # each plate's frames in turn: the town, the year playing, a zoom
        fs = sorted(mdir.glob(f"m{k}*.jpg"))
        seq += [(shot(f), secs / len(fs)) for f in fs]
    seq += [(end, END_S)]
    out = ROOT / "web" / "demo"
    out.mkdir(exist_ok=True)
    tmp = out / "_frames"
    tmp.mkdir(exist_ok=True)
    n, t_edges, acc = int(DUR * FPS), [], 0.0
    for _, s in seq:
        t_edges.append((acc, acc + s)); acc += s
    for k in range(n):
        t = k / FPS
        i = next(j for j, (a, b) in enumerate(t_edges) if t < b or j == len(t_edges) - 1)
        a, b = t_edges[i]
        im = seq[i][0]
        x = min(0.35, (b - a) * 0.6)                          # short frames of a playing year blend most of the way
        if i + 1 < len(seq) and t > b - x:
            im = Image.blend(im, seq[i + 1][0], (t - (b - x)) / x)
        im.save(tmp / f"f{k:04d}.jpg", quality=90)
    vdir = Path(sys.argv[2])
    clips = sorted(vdir.glob("c*.wav"), key=lambda q: int(q.stem[1:]))
    for c, at in zip(clips, VOICE_AT):
        with wave.open(str(c)) as w:
            print(f"{c.name} at {at:.1f} s, ends {at + w.getnframes() / w.getframerate():.1f} s")
    ff = imageio_ffmpeg.get_ffmpeg_exe()
    mp4 = out / "edelweiss-demo.mp4"
    cmd = [ff, "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", str(tmp / "f%04d.jpg")]
    for c in clips:
        cmd += ["-i", str(c)]
    mix = "".join(f"[{k + 1}:a]adelay={int(at * 1000)}|{int(at * 1000)}[a{k}];" for k, at in enumerate(VOICE_AT[:len(clips)]))
    mix += "".join(f"[a{k}]" for k in range(len(clips))) + f"amix=inputs={len(clips)}:normalize=0,apad[aout]"
    cmd += ["-filter_complex", mix, "-map", "0:v", "-map", "[aout]", "-t", str(DUR), "-c:v", "libx264", "-pix_fmt", "yuv420p",
            "-crf", "20", "-c:a", "aac", "-b:a", "128k", "-metadata", "comment=Narration: synthetic voice; script drafted with AI (Claude)",
            "-movflags", "+faststart", str(mp4)]
    subprocess.run(cmd, check=True)
    for p in tmp.iterdir():
        p.unlink()
    tmp.rmdir()
    print(f"{mp4.relative_to(ROOT)}: {mp4.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
