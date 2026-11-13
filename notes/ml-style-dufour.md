# ML "engraver's hand" for the grey relief (pet project, not for the hackathon)

Idea: an ML model draws the static grey relief in the style of the Dufour map. The red data
hachures stay algorithmic, which keeps them exact and lets them animate without flicker.
"ML for beauty, physics for truth."

## Data recon (2026-11-10)
- Layer `ch.swisstopo.hiks-dufour` (Dufour map, first edition, 1:100 000, 1845-1865, 25 sheets).
  Free WMTS tiles, no key:
  `https://wmts.geo.admin.ch/1.0.0/ch.swisstopo.hiks-dufour/default/current/3857/{z}/{x}/{y}.png`
  - The tiles are already georeferenced in Web Mercator.
  - z13 (~13 m/px at 46°N) looks best.
  - From z15 the tiles are upscaled; the native scan resolution is around z14.
- Samples: `img/dufour_jungfrau_z13.png` and `img/dufour_lauterbrunnen_z15.png`.
- License: swisstopo data are free for any use with the credit "© swisstopo".
  - Georeferenced files of the old maps (as opposed to tiles) are paid, on request via geodata@swisstopo.ch.
  - Check the exact wording before publishing.
- Modern DEM for pairs: swissALTI3D (free).
  - Caveat: the 19th-century relief is off by tens to hundreds of metres, so pixel pairs are inaccurate.

## What the tiles show (and what must be cleaned)
- Lots of text: place names and spot heights.
- Forest dots, roads, rivers.
- High mountains drawn differently: rock drawing and white glaciers.
- Visible seams between sheets, with different print darkness.
- For training, keep only tiles that are mostly hachures, with no text.

## Plan, simplest first
1. Stable Diffusion 1.5 + a LoRA on 100-300 clean tiles + an off-the-shelf ControlNet (lineart or depth).
   - Input: our own hachure render or a hillshade.
   - LoRA training takes 1-3 h on the Quadro RTX 5000 (16 GB). Generation takes a few seconds per 512 px tile.
2. CycleGAN "our renders → Dufour scans" (no pairs needed). 1-2 days of GPU per attempt.
3. Our own ControlNet on DEM→Dufour pairs. Needs a cloud A100 and weeks of work. Not worth it for now.
