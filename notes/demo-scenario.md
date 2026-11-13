# Demo scenario: "Honest map" (draft, 11.11, for team discussion)

The idea in one sentence: **the NISAR radar sees how the ground moves, and our map shows it with strokes, like old engravings, and
honestly says where it is not confident.**

Places and colors may change: the scenario rests on four "skills" of the engine, not on specific cities.
Each skill = one layer type (see `web/README.md`). To replace a place you need a new `sites/*.json` and a rebuild.

## 3 minutes

| time | what is on screen | what we say (gist) |
|---|---|---|
| 0:00–0:20 | Mexico City, strokes have not grown yet | The ground under millions of people sinks by centimeters per year. The NISAR satellite measures this every 12 days. Ordinary radar maps look like a rainbow that only specialists can read, and they never say where they are wrong. |
| 0:20–1:10 | We press ▶: the strokes grow over a year (Oct 2025 → Oct 2026) | Each ring of strokes = 2 cm of subsidence, the deepest is 28 cm. Pale strokes = weak signal (1–3σ); we do not draw noise at all. We turn on "chain breaks": only ~9% of points have an unbroken series, and the map shows it. A "radar only" badge; next to it, a check against GNSS station MMX1: −26.9 radar vs −27.0 cm/year. |
| 1:10–1:50 | Panama M7.7: the earthquake card | From two passes, the engine produces the fault and the magnitude in a few seconds (12 runs of ~0.5 s each). The key thing on the card is the **"NOT CONFIDENT"** flag. Strike 285° (USGS 288°) matched; dip and magnitude did not, and the card warned about that. On synthetic tests, the flag separates correct cards from wrong ones (AUC 0.93). |
| 1:50–2:30 | Central Valley: "coherence on a date", ▶ | This is the "dance": fields in different months "flicker" differently for the radar (plowing, irrigation, harvest). We switch to "4 rhythm types" and click a field to see its curve against the mean curve of its type. |
| 2:30–2:45 | Langtang, Nepal: debris-flow scar | The same measure, but for sudden changes: a 4.8 km² patch where coherence dropped in a single pass. As a control, the same measure before the event gives 0.1 km². Honest label: the rule was found on this event and has not yet been tested blind. |
| 2:45–3:00 | Style switch (engraving → night), places menu | The engine is open. A new place = one settings file, a new style = one theme file. From 13.11 we plug in the official NISAR data. |

Backup contrast, if there is time: Houston, descending pass (errors are the same for all points on a date → "why
you cannot trust a single pair"). There is no ready scene for it yet.

## What we need to have before the demo

- [ ] Four scenes on official NISAR data from the challenge resources (for now, data we found ourselves).
- [ ] Check the Space Apps rules on work done before the hackathon (Participant FAQ). How we present the engine depends on this.
- [ ] A 30 s backup video: screen recording with ▶ on Mexico City and on the fields (in case the internet is slow).
- [ ] One page "how to read the strokes" (a legend for the judges).
- [ ] Hosting: the static `web/` folder → GitHub Pages in the team's new repository.

## Friday, 13.11: who does what (a proposal, the team decides)

| role | what exactly | who |
|---|---|---|
| NISAR data | open the resources tab, find GUNW for our places, fill the cache, rebuild the scenes | Maksym (already knows the engine) |
| New places | pick 1–2 places with a strong story (flood? glacier? Ukraine?), write `sites/*.json` | ? |
| Look | theme, fonts, mobile version, test the legend on "people from the street" | ? |
| Story and submission | script, submission text, video, slides | ? |
| Verification | GNSS, USGS, comparison with published maps; someone who looks for where we are wrong | ? |
| Repository | new shared GitHub, move `src/`, `scripts/`, `sites/`, `web/`, license, README | ? |

Open questions for the team: do we keep the debris flow (post-hoc for now)? do we do a blind earthquake live if an M6+ happens?
what order for the acts: Mexico City first or the fields first (the "dance" hook)?
