# Data on disk (pre-hackathon sandbox), 2026-11-10

| What | Where | For |
|---|---|---|
| NISAR GUNW windows, Mexico City, descending track 113 frame 79: 27-pair chain, 2025-10-25 → 2026-10-20 | `data/raw/nisar_cache/bbox_.../` | subsidence map, animation, noise statistics |
| Same windows, ascending track 19 frame 11: 29 pairs (download started 2026-11-10) | same cache, output `data/processed/mexico-city-asc/` | asc+desc split into vertical and east-west, second independent noise sample (morning vs evening pass) |
| GNSS daily positions (Nevada Geodetic Lab, IGS20 frame), 22 stations in the Mexico City box | `data/raw/gnss/*.tenv3` | ground truth |
| Earthquake windows (co-seismic pairs + 2 pre-event pairs each) | `data/raw/nisar_events/<name>/` | Okada inversion and simulation-based inference tests |

## GNSS stations
- Only MMX1 (airport) has data covering the whole NISAR year.
- ICMX, UNIP, UNVA, TNGF, MXTX and MXAJ have long histories that end 2022–2025. Use them for long-term rates only, since the subsidence is close to linear.

## Earthquakes
- Nagqu, Tibet: M7.0, 2026-10-22, 18 km deep. Co-seismic pair 10-19 → 10-31 on ascending track 112 frame 18. About 1 m of line-of-sight motion and clean fringes.
- Kumamoto, Japan: M6.8, 2026-07-28, 8 km deep. Download started 2026-11-10.
- Other candidates from the USGS scan, shallow M ≥ 5.8 events spanned by NISAR pairs:
  - Panama, M7.7, 2026-10-09
  - Venezuela, M7.5, 2026-06-25
  - Indonesia (Flores), M7.8, 2026-08-15
  - Hubbard Glacier, M7.0, 2025-12-06

## Scripts
| Script | What it does |
|---|---|
| `scripts/fetch_event.py` | Earthquake windows plus quick-look images |
| `scripts/atmo_stats.py` | Noise statistics; output in `data/processed/mexico-city/atmo_stats.{json,png}` |
| `scripts/gnss_check.py` | NISAR vs GNSS comparison |

## Environment
- PyTorch 2.14.1 + CUDA 12.6 is installed in the project env (`uv sync --group ml`).
- It works on the Quadro RTX 5000.
