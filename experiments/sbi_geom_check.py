"""D diagnostic (blind to USGS): does the real pair fit better under another LOS geometry?

    uv run --group ml python scripts/sbi_geom_check.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import torch

sys.path[:0] = [str(Path(__file__).parent), str(Path(__file__).parents[1] / "scripts")]
from sbi_fault import HEADING, Sim, load_event, lsq  # noqa: E402


def main():
    torch.manual_seed(0)
    for event in ("nagqu", "kumamoto"):
        ev = load_event(event)
        sim = Sim(ev)
        th = torch.acos(sim.lv[2])
        for name, az, sgn in [("left-looking (current)", HEADING + 90, 1), ("right-looking", HEADING - 90, 1),
                              ("left, LOS sign flipped", HEADING + 90, -1), ("vertical only", None, 1)]:
            if az is None:
                sim.lv = torch.stack([torch.zeros_like(th), torch.zeros_like(th), torch.ones_like(th)])
            else:
                a = math.radians(az)
                sim.lv = sgn * torch.stack([torch.sin(th) * math.sin(a), torch.sin(th) * math.cos(a), torch.cos(th)])
            best, rms, ystd = lsq(sim, starts=96, iters=1200)
            print(f"  {event:8s} {name:24s} rms {rms * 100:5.1f} cm (data {ystd * 100:.1f})  strike {best['strike'] % 360:.0f}"
                  f" dip {best['dip']:.0f} rake {(best['rake'] + 180) % 360 - 180:.0f} slip {best['slip']:.1f}", flush=True)


if __name__ == "__main__":
    main()
