# Step texts for the "Honest map" (draft 11.11, for editing)

How to read: each step = what is visible on the map + on-screen text (1–2 sentences, no jargon) + how we check it against reality.
Marks: ✅ data already exists · 🔧 needs to be done on the map · 🧪 needs to be checked with a test (pre-registered) · ❓ fact needs to be checked against a source.

The through-line: **"The ground does not stand still"**. Slowly, where people pump water (Mexico City, California). Suddenly, an earthquake and a collapse (Panama, Nepal). Each time the map says where it is sure and where it is not.

---

## 0. Intro (over a world map)

**0.1** Visible: a globe or world map, four place dots.
> The NISAR satellite scans the whole Earth with radar every 12 days: through clouds, day and night. By comparing two passes, it notices that the ground has moved by a few centimeters.

**0.2** Visible: the relief strokes of Mexico City, large, like an engraving.
> We draw this motion with strokes, the way 19th-century cartographers drew mountains: the denser the strokes, the larger the motion. And where we are not confident, the strokes fade or disappear.

---

## 1. Mexico City: a sinking city

**1.1** Visible: the whole city, relief only, names of a few districts. 🔧 (names)
> Mexico City stands on the bed of a former lake. Under it is soft clay, and for decades the city has been drinking water from under its own feet.

**1.2** Visible: ▶, the strokes grow from October 2025 to October 2026. ✅
> In one year, the ground at the deepest spot sank by **28 cm** ❓. Each new ring of strokes is another 2 cm.

**1.3** Visible: zoom to station MMX1, a station marker. 🔧 (marker)
> Can we trust this? There is a GPS station here, MMX1. It says: −27.0 cm per year. The radar says: −26.9. ✅

**1.4** Visible: pale strokes plus cells with series breaks. ✅
> Pale strokes are where the signal is weak. We do not hide uncertainty, we draw it.

❓ **Clarify the number before the demo.** The map says "deepest 28 cm" (smoothed picture), but individual strokes have rates down to −46 cm/year (near −99.05, 19.27, southeast of the city). One cell shows −72 cm/year with an error of 7; this looks like an artifact. Until we sort out which number to quote, we do not put it on screen. Also check against publications how much Mexico City sinks per year.

---

## 2. Central Valley, California: farms and water

Two options; we choose after the test.

### Option A (stronger, but needs a test 🧪): "the ground breathes with irrigation"

**2.1** Visible: the valley between mountains, city names (Fresno, Hanford, Corcoran). 🔧
> The Central Valley is one of America's main food-growing regions ❓ (check the share of the harvest). When there is little rain, the farms pump groundwater.

**2.2** Visible: the subsidence strokes over the year.
> The ground under the fields is sinking. GPS stations in the valley show up to **−17 cm per year** ✅ (from the station data; we put up a radar number only after the test 🧪).

**2.3** Visible: a month slider; the curve of one station, for example HANF or CAHA.
> The ground breathes with irrigation: in summer they pump, and it sinks by 5–6 cm; in winter it almost stands still. ✅ (GPS: May–September −5…−6 cm, November–March ≈ 0)

**2.4** (if X-season passes 🧪)
> Part of this subsidence will not come back: compressed clay no longer holds water.

⚠ **Trap.** Our radar series in the valley has a known winter bias: about +10 cm of fake uplift due to the ionosphere correction (P-port). It overlaps exactly with the seasonal signal. The test must show that the radar (with the "plane" correction) sees the same summer subsidence as GPS. Otherwise we do not show option A.

### Option B (fallback, data exists ✅): "when the fields are plowed"

> The radar notices when a field changes: plowed, irrigated, harvested. In summer the fields "flicker", in winter they stay calm.

Weaker, because without a crop map (USDA CDL) we do not know what exactly happened to each field.

---

## 3. Panama: an earthquake

**3.1** Visible: the coast, the peninsula, village names, the epicenter. 🔧 (coast and names)
> 9 October 2026. A magnitude 7.7 earthquake, 10 km from the village of Pitaloza Arriba ✅ (USGS) on the Azuero Peninsula ❓. More than 2000 people reported to USGS that they felt it ✅.

**3.2** Visible: strokes of the ground jump (ring = N cm), a single "before and after" image instead of a year. 🔧
> The radar compared images before and after. The ground moved within seconds: on the morning image up to ~40 cm ✅, and near the epicenter even more ❓ (check against the data).

**3.3** Visible: the fault line under the strokes.
> From this picture, the engine reconstructs the fault underground in a few seconds: about 135 km long, at a depth of ~11 km.

**3.4** Visible: a card with the "NOT CONFIDENT" flag.
> But the engine itself says: "I am not confident". The fault direction matched the seismologists (285° vs 288°), but the dip angle and the strength did not. On test earthquakes, this flag separates correct answers from wrong ones. ✅

**3.5** (optional) Visible: the evening image.
> Two and a half hours later another shock hit, M6.6, 76 km to the west ✅. The evening image sees both together, which is why the picture is so tangled.

❓ On the evening image there are values over a meter near the epicenter. Check whether this is real displacement or a phase unwrapping error before showing it.

---

## 4. Langtang, Nepal: a mountain collapse

**4.1** Visible: relief in strokes, label "Langtang Lirung, 7227 m". 🔧 (label)
> The Langtang valley in Nepal, north of Kathmandu ❓ (check the distance).

**4.2** Visible: a mark at the place it broke off from (~5000–5500 m).
> 26 August 2026, 08:37. About 67 million cubic meters of ice and rock collapsed from the wall of Langtang Lirung ✅ (ICIMOD). The avalanche turned into a debris flow and went down the valley.

**4.3** Visible: the debris-flow scar as engraving hatching instead of an orange patch. 🔧 (hatching style)
> The radar saw where the debris flow went: a 4.8 km² strip, almost 6 km down the valley, from 4000 to 2650 m ✅. After the event it simply "did not recognize" this ground.

**4.4** Visible: the same measure on pairs before the event, empty.
> A week and a month before, the same measure showed nothing ✅. This is not noise, it is an event.

**4.5** Visible: the collapse source, in gray.
> Honestly: our radar could not see that the mountain was getting ready to fall. Ice at 5000 m "flickers" every summer and does not allow measuring motion ✅ (L-retro). For that, more precise NISAR products are needed.

**4.6** (desirable) Visible: Sentinel-2 satellite photo before and after, side by side. 🔧 (needs an OK for the download)
> And here is how it looks on an ordinary photo.

⚠ For us, not for the screen: the same rule gives a false patch in September (autumn "flicker"). Before making this into a general layer, we need a seasonal safeguard.

---

## 5. Finale

> Slowly or suddenly, the ground moves. NISAR sees it every 12 days across the whole planet. We turn this into a map that anyone can read, and we honestly mark where we are not sure. A new place is one settings file.

---

## What to decide

1. Order: Mexico City → California → Panama → Nepal (slow → sudden)? Or start with Panama, because it is the most dramatic?
2. California: do we run the test for option A, or take B right away?
3. Do we keep step 4.5 "what we could not do"? I like it: it is honesty, not weakness.
4. Downloads: coast, rivers and names (Natural Earth / OSM); Sentinel-2 photos for Nepal. Both are public, no login, but your OK is needed.

## What to build after the edits (briefly)

- Story mode: steps in `sites/*.json` (`story: [{view, layers, time, text}]`), "Next" and "Explore on your own" buttons.
- A names and coastline layer in engraving style; station markers.
- The ground jump for Panama as strokes (like Mexico City, but one pair).
- The debris-flow scar as hatching.
- Intro over a world map.
