# Reykjanes research pipeline

The analysis behind the Reykjanes plate (`web/reykjanes-scroll.html`) and its math window (`docs/uk/reykjanes-math.html`):
NISAR uplift over Svartsengi, the source (Mogi, ellipsoid, layered-crust FE), MCMC vs NPE, the 10 Nov 2023 dike from
Sentinel-1 + GNSS, FIRMS / HLS / coherence checks, and the Kīlauea polarization study. Not part of the `sars` CLI yet.

Run every script from the repo root, e.g. `uv run python scripts/reykjanes/build_breathing_data.py`. Scripts in this folder
import each other by module name, so run them by path, not as a package.

## Paths and secrets

- Inputs go to `data/` and results to `output/` at the repo root; both are ignored by git.
- Secrets are read from the environment only (`_secrets.py`): `EARTHDATA_TOKEN` or the file `~/.edl_token` (NASA Earthdata),
  `FIRMS_MAP_KEY` (NASA FIRMS area API). Nothing is read from files inside the repo.
- Extra packages used only here: `uv sync --group ml --group research` installs `torch` and `zuko` (NPE and real-time dike
  scripts) and `shapely` and `osmium` (OSM extract). Name every group you need in one `uv sync`: groups left out are removed.

## Order

| step | scripts | writes |
|---|---|---|
| 1. NISAR windows | `download_site.py reykjanes` | `data/reykjanes/gunw/` |
| 2. Context | GNSS `.tenv3` (NGL) into `data/gnss/`; Copernicus DEM tile into `data/dem/`; `firms_hotspots.py`; `osm_extract.py` (needs `data/osm/iceland-latest.osm.pbf` from Geofabrik) | `data/gnss/`, `data/dem/`, `data/firms/`, `data/osm/reykjanes_layers.json` |
| 3. Uplift and source | `reykjanes_analysis.py`, `reykjanes_mogi.py`, `reykjanes_ellipsoid.py`, `reykjanes_gnss.py`, `build_breathing_data.py`, `build_breathing_app.py` | `output/reykjanes/*.npz`, `breathing_*.{npz,json}`, `map_layers.json` |
| 4. Inference | `run_pipeline.py` (observation cache, FE table and MCMC, NPE flow vs MCMC, real-time dike, FE NPE) | `output/reykjanes/ml/` |
| 5. Real dike (Sentinel-1) | `run_s1_dike.py` / `hyp3_wait.py` (HyP3 jobs), `s1_grindavik_dike.py`, `s1_dike_npe_v2.py`, `s1_dike_npe_v3.py`, `s1_dike_npe_v4.py`, `s1_dike_distributed.py`, `coherence_fissures.py`, `hls_*.py`, `results_v4_section.py` | `output/reykjanes/ml/s1_*.json`, `coherence_fissures.json`, `firms_hotspots.json` |
| 6. Pages | `build_results_page.py` (report), `../build_reykjanes_scroll.py` (the plate) | `output/reykjanes/results.html`, `web/reykjanes-scroll.html` |

Kīlauea: `download_site.py kilauea`, then `kilauea_episodes.py`, `kilauea_polarization.py`, `kilauea_tephra_figure.py`.

`progress.py` is a small helper that writes task status to `output/progress/`; the status server itself is not included.
Paths were made repo-relative when the scripts moved here; the full chain has not yet been re-run on a clean clone.
Tests and their verdicts are logged in `notes/experiments.md`, block R (a post-hoc record, plus pre-registered R5).
