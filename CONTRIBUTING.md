# Contributing

Short rules so that six people and their Claude Code sessions can work in one repo for 48 hours without stepping on each other.

## Setup (once)

See [README.md](README.md#quick-start-the-same-for-everyone): install uv, clone into a short folder, run `uv sync`, save
your own Earthdata token, then run `uv run sars doctor`.

## Daily loop

1. `git pull` on `main`, then branch: `git switch -c <your-name>/<topic>` (e.g. `oleh/gcov-plate`).
2. Work. Fetch data only with `uv run sars fetch <plate>` (see [docs/data-guide.md](docs/data-guide.md)).
3. Before pushing, run `uv run python tests/smoke.py` (offline, ~15 s). If you changed a plate, also run
   `uv run sars build <plate>` and look at it with `uv run sars serve`.
4. Push the branch and open a pull request. Ask one teammate to look at it. Merge when the smoke test passes and someone
   has looked.

## What goes where

| You are adding | Put it in |
|---|---|
| a new place | `uv run sars new <id> --bbox W S E N ...` creates `sites/<id>.json`; pages go in `sites/<id>.story.json` |
| engine code | `src/hachure/`; the CLI step in `src/hachure/cli.py` |
| a pipeline step the CLI calls | `scripts/` |
| a one-off test | `experiments/`, with its pre-registration in `notes/experiments.md` |
| a look (colours, fonts) | `web/themes/<name>.json` (see `web/README.md`) |
| notes, plans, research | `notes/` (English); Ukrainian team docs go in `docs/uk/` |

## Commits and pull requests

- English, imperative mood, one topic per commit: `Add GCOV water mask to the Houston plate`.
- Commit the lock file (`sites/<id>.lock.json`) when a plate's pairs change. Commit `web/data/<id>/` when a plate is
  rebuilt for the demo.
- Never commit `data/`, tokens, `.env` files, videos, or files with personal emails. `.gitignore` covers the usual cases;
  check `git status` anyway.
- Large files: keep each committed file under ~10 MB (the biggest built layer today is 5.4 MB).

## Science rules

- **Pre-register**: write the prediction and thresholds in `notes/experiments.md` before running a test that supports a
  claim on a plate. Do not move a threshold after seeing the data. Label anything found later as **post-hoc**.
- **Report failures**: a failed test goes in the log and, if it touches a plate, on the plate.
- **Credit data**: every dataset is listed with its licence in `README.md` and `docs/data-guide.md`.

## Working with Claude Code

- Every session reads `CLAUDE.md` first; keep it current when you change the layout or the rules.
- Give your session a narrow task and a plate or folder to work in. Two sessions should not edit the same file at the same
  time.
- Your session must never see or print your Earthdata token. If it asks for the token, refuse, and point it at
  `~/.edl_token`.
