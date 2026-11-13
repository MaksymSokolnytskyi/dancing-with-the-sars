# Hachures × NISAR: breaking down the idea

## TL;DR
The idea is strong, but one thing in Gemini's advice needs fixing:
**InSAR (interferometry) gives a scalar, not a vector**: displacement along the line of sight (LOS).
Arrows "from LOS" are artificial vectors. Real vector fields from SAR do exist, and there are three types:

| Vector source | What it is | NISAR product | Where it works best |
|---|---|---|---|
| **Pixel offsets** (amplitude correlation) | true 2D displacement (range + azimuth), meters | **GOFF** (L2), ROFF (L1) | glaciers, large ground displacements, coseismic displacements |
| **Gradient of a scalar field**: the classic Lehmann/Dufour hachure | −∇f of any raster | GUNW (LOS displacement), GCOV (backscatter), DEM | subsidence bowls, "slopes" of change |
| **Asc + desc decomposition** | E–W + Up from two geometries | 2× GUNW | urban subsidence, volcanoes |

So the hachure is not decoration but an **operator**: take any field (relief, displacement, backscatter
change) → compute the gradient → trace lines of steepest descent. The same engine
for all layers.

## What your examples show
- **Everest (Washburn / Swissair map)**: rock is gray hachures, **ice is blue lines,
  and flow lines are visible on the Kangshung glacier**. This is exactly our visual grammar:
  - gray stroke = static relief (DEM),
  - **blue stroke = ice/water motion (NISAR offsets / water mask)**,
  - a third color (red/ochre) = change (Δ backscatter, subsidence).
- **Dufour Blatt XVII**: northwest illumination creates volume without any color →
  color stays free as a data channel. **This is where the information value lies**:
  monochrome hachures encode direction + magnitude, and color encodes a second variable.

## Answer to "is this just about beauty?"
No, if every visual parameter of the stroke is tied to physics:
- **direction** = direction of motion / of the gradient,
- **length/thickness** = speed / steepness (Lehmann's law),
- **density** = InSAR **coherence**. Low coherence → we do not draw → a white patch
  honestly says "the radar is not confident here". This is the radar version of Dufour's generalization:
  he removed noise deliberately, we do it by a quality metric.
- **animation** = time series (particles moving along flow lines; date on a slider).

## Study area candidates
1. **Glacier** (Greenland: Sermeq Kujalleq/Jakobshavn; Alaska: Malaspina):
   GOFF gives true 2D vectors → fits the flow-line idea best. Control:
   NASA ITS_LIVE (ready-made velocity fields).
2. **Mexico City: subsidence** ~30+ cm/year: visible even over short intervals; the subsidence
   bowl = an ideal "upside-down mountain" for hachures.
3. **Former Kakhovka Reservoir**: a personal Ukrainian story: the bed has been growing into forest since 2023.
   L-band is sensitive to woody vegetation → backscatter change (HV) as
   a "relief" that grows. Strong storytelling, but NISAR coverage needs checking.
4. Flood / high water: a combination of a water mask from GCOV + hydrological hachures from the DEM
   (where the water drains). Your experience with NLCD masks applies directly here
   (ESA WorldCover / NLCD as context).

## What to do BEFORE 13.11 (while the resources are closed)
- [x] Hachure engine on a DEM (Copernicus GLO-30): gradient → seeding (Jobard–Lefer
      evenly-spaced streamlines) → RK2 → Lehmann thickness + Lambert shading → GeoJSON/SVG.
      (2026-11-10, `src/hachure/`, examples: Dents du Midi, Everest)
- [ ] The same engine on a vector field (ITS_LIVE as a stand-in for GOFF).
- [ ] Check what is already in ASF for NISAR (`asf_search`, dataset NISAR): formats
      HDF5, structure of GCOV/GUNW/GOFF.
- [ ] Front-end prototype: MapLibre + deck.gl, a line layer + particle animation + before/after swipe.

## Small corrections to the Gemini text
- Do not swap NISAR for Sentinel-1/ALOS in the final version: the challenge is specifically about NISAR.
  S1/ALOS only for comparing frequencies (C vs L) or as a fallback.
- The algorithm is Jobard & **Lefer** (1997), not "Leferve".
- The claims about Synthetic Vision in avionics "as hachures" are doubtful; do not use them in the pitch.
