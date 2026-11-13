# Experiments

One-off scripts behind the entries in `notes/experiments.md` (each test was pre-registered with its thresholds before it
ran; post-hoc findings are labelled there). They read the same caches as the pipeline (`data/raw/nisar_cache/`), so run
`uv run sars fetch <plate>` first. Notes written before 13.11.2026 call them `scripts/<name>.py`; they now live here.

    uv run python experiments/<name>.py
