# NASA Space Apps 2026 — "Dancing with the SARs"

Context for any Claude session working in this folder.

## Challenge
Build an interactive, open-source app that uses **NISAR** radar data to track and visualize
one or more types of surface change (floods, deforestation, glaciers, subsidence, quakes...).
Page: https://www.spaceappschallenge.org/2026/challenges/dancing-with-the-sars/
Official resources tab opens **Fri 2026-11-13**. Until then: prototype algorithms on
substitute data so we can switch to NISAR immediately.

## Team (Ukraine)
- Andrii Bilych (@vox_astutiae) — team owner
- Oleh Paziuka (@prorabvasili)
- Maksym Sokolnytskyi (@sokolnytskyi)
- Yehor Matiiash (@kiber_atlett)
- Illia Bondar (@kirillovw)
- + Maksym's wife (joining later)

## Shared spaces
- Brainstorm doc: https://docs.google.com/document/d/1Aq228TEg3yEg1NPnly5zXZHRSBcp_cJrYWedwIv15ac/edit
  - Readable without a connector at `.../mobilebasic` (link-shared).
  - Tab "Claude threads" → "main room": Claude-to-Claude threads, protocol in `notes/claude-threads-protocol.md`.
- GitHub: a fresh shared repo will be created when the hackathon starts. This folder is the pre-hackathon sandbox.

## Core idea (working hypothesis)
Render SAR-derived change as **hachures / streamlines** (Dufour/Lehmann-style engraving),
not rainbow fringes. See `notes/hachures-x-nisar.md`.

## Layout
- `notes/` — thinking, decisions, protocols
- `data/raw/`, `data/processed/` — downloads (gitignored later)
- `src/` — Python pipeline (DEM/SAR → vector field → streamlines → GeoJSON/SVG)
- `notebooks/` — experiments
- `web/` — interactive front-end (MapLibre / deck.gl)

## Tooling on this machine
Python via `py` (3.12) and `uv`; Node 26; `gh` logged in as MaksymSokolnytskyi.
