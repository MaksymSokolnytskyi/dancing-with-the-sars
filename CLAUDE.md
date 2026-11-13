# CLAUDE.md: context for Claude Code sessions in this repo

NASA Space Apps Challenge 2026, challenge **Dancing with the SARs**. Team from Ukraine (Dnipro location). Several team members
work in this repo through their own Claude Code sessions. Read this file, then [README.md](README.md).

## What we build

An interactive atlas of surface change seen by the **NISAR** radar (L-band, 12-day repeat). Each **plate** is one place and
one kind of change (a city sinking, a debris flow, fields through the year, an earthquake). The plates are drawn as 19th-century
engraved hachures (Dufour, Lehmann) instead of rainbow fringes, and each one tells its story as a few short pages.
Each plate also says how sure it is.

- Challenge page: https://www.spaceappschallenge.org/2026/challenges/dancing-with-the-sars/
- Deadline: 15 Nov 2026, 23:55 local. Everything shown to the jury is in English.

## Map of the repo

| Path | What | Entry point |
|---|---|---|
| `src/hachure/cli.py` | the `sars` command: doctor, sites, search, fetch, lock, build, story, pack, serve, new | `uv run sars --help` |
| `src/hachure/nisar.py` | NISAR GUNW: ASF search, signed URL, HDF5 window read, `RawWindow` cache, `process` to vertical motion | `read_raw`, `process` |
| `src/hachure/gcov.py` | NISAR GCOV (backscatter HH/HV): search, window read, cache | `search`, `fetch`, `Window` |
| `src/hachure/engine.py` | site config to web layers: velocity stack, trust, hachures, context, manifest | `load_site`, `cached_chain` |
| `src/hachure/field.py`, `trace.py` | the hachure algorithm (fields, then stroke tracing) | `field.prepare`, `trace.hachures` |
| `src/hachure/trust.py` | per-pixel error budget behind "how sure" | `audit` |
| `sites/<plate>.json` | one plate: bbox, direction, layers, context, pages (`.story.json`), pinned pairs (`.lock.json`) | |
| `scripts/` | pipeline steps the CLI calls (`build_site`, `story`, `pack_scene`, fetchers, `export_quake`, `sbi_*`) | |
| `web/` | the map app (MapLibre + deck.gl, no build step); `web/data/<plate>/` = built plates, in git | `web/README.md` |
| `experiments/` | one-off pre-registered tests behind `notes/experiments.md` | |
| `notes/` | plans, research, experiment log; index in `notes/README.md` | |
| `docs/` | data guide; `docs/uk/` = Ukrainian onboarding | `docs/data-guide.md` |
| `tests/smoke.py` | offline smoke test (no network, no token) | `uv run python tests/smoke.py` |

## Commands

```bash
uv sync                                   # add --group ml for PyTorch (SBI / quake experiments only)
uv run sars doctor                        # environment, token, network, path length
uv run sars fetch <plate> --dry-run       # what would be downloaded
uv run sars fetch <plate>                 # DEM, NISAR windows (pinned by the lock), OSM, WorldCover, context, GNSS
uv run sars build <plate>                 # -> web/data/<plate>/
uv run sars serve                         # http://localhost:8765/?scene=<plate>&theme=scroll
uv run python tests/smoke.py              # before every push
```

Datasets, licences, cache formats and per-plate needs are in [docs/data-guide.md](docs/data-guide.md).

## Rules for every session

1. **Secrets.** The Earthdata token lives in `~/.edl_token` or `EARTHDATA_TOKEN`, never in the repo, a chat, an issue or a
   doc. Never print it; `sars doctor` reports only its length. Never create accounts or type passwords for a teammate.
2. **Data stays out of git.** `data/` is ignored. Fetch with `sars fetch`, never by hand-copying big files. Read NISAR only by
   window (MBs, not the 2 GB files). Before any new download larger than ~1 GB, say the size and ask the human.
3. **Pre-register.** Before a test that supports a claim on a plate, write the prediction and the pass/fail thresholds in
   `notes/experiments.md`, then run it. Never move a threshold after seeing the result. Label anything found after the fact
   as **post-hoc**.
4. **Honest plates.** Weak signal is drawn faded or dashed, numbers carry their error, and failed tests are reported.
   Do not write a claim on a plate that a test did not support.
5. **English** in code, comments, commit messages, docs and everything the jury sees. The web UI is bilingual (`en`, `uk`)
   through `{ "en": ..., "uk": ... }` labels; keep both when you add text. Ukrainian docs for the team live in `docs/uk/`.
6. **No personal data** (emails, phone numbers) in committed files or in requests to outside services. Use the User-Agent
   `space-apps-research`.
7. **Small, reviewable changes.** Branch, run the smoke test, open a PR (see [CONTRIBUTING.md](CONTRIBUTING.md)).
   Do not push to `main` directly unless the team agrees.
8. **Credits.** A new dataset goes into the credits in `README.md` and into `docs/data-guide.md` with its licence.

## Team

| Name | GitHub | Telegram |
|---|---|---|
| Andrii Bilych (team owner) | `Bodomerka` | @vox_astutiae |
| Oleh Paziuka | `ProrabVasili` | @prorabvasili |
| Maksym Sokolnytskyi (repo owner) | `MaksymSokolnytskyi` | @sokolnytskyi |
| Yehor Matiiash | `edoTenseii` | @kiber_atlett |
| Illia Bondar | `kirillovw` | @kirillovw |
| Yelyzaveta Razumovska | `geekata` | |

Use the GitHub login (with `@`) only to request a review or assign a PR; GitHub notifies that person. The team's
brainstorm doc and the Claude-to-Claude thread protocol are described in `notes/claude-threads-protocol.md`.
