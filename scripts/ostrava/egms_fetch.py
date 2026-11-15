"""Download the EGMS L3 Ortho vertical tiles over the Ostrava–Karviná box (for test O3b, notes/experiments.md).

    EGMS_TOKEN=... uv run python scripts/ostrava/egms_fetch.py

EGMS (European Ground Motion Service, Copernicus Land Monitoring Service) needs a free account; the download links in
the EGMS Explorer carry a personal `id=` key. Put it in the environment as EGMS_TOKEN; it is never printed or stored.
Writes data/ostrava/egms/<tile>.zip and the unpacked GeoTIFF (ignored by git). Skips tiles already there.

Product: L3 Ortho, vertical component (U), release 2019–2023, 100 m grid, 100 km tiles in EPSG:3035.
The box 18.13–18.83 E, 49.68–50.07 N falls on tiles E49N29 and E49N30.
"""
from __future__ import annotations

import os
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "ostrava" / "egms"
URL = "https://egms.land.copernicus.eu/insar-api/archive/download/{name}.zip?id={token}"
TILES = ["E49N29", "E49N30"]
RELEASE = "2019_2023_1"


def main() -> int:
    token = os.environ.get("EGMS_TOKEN", "").strip()
    if not token:
        print("EGMS_TOKEN is not set (the id= key from an EGMS Explorer download link)")
        return 1
    OUT.mkdir(parents=True, exist_ok=True)
    for tile in TILES:
        name = f"EGMS_L3_{tile}_100km_U_{RELEASE}"
        dst = OUT / f"{name}.zip"
        if not dst.exists():
            req = urllib.request.Request(URL.format(name=name, token=token), headers={"User-Agent": "space-apps-research"})
            tmp = dst.with_suffix(".part")
            with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as f:
                while chunk := r.read(1 << 20):
                    f.write(chunk)
            tmp.rename(dst)
        with zipfile.ZipFile(dst) as z:
            tifs = [m for m in z.namelist() if m.lower().endswith((".tif", ".tiff"))]
            for m in tifs:
                if not (OUT / Path(m).name).exists():
                    (OUT / Path(m).name).write_bytes(z.read(m))
        print(f"{name}: {dst.stat().st_size / 1e6:.1f} MB, {', '.join(Path(m).name for m in tifs) or 'no GeoTIFF inside'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
