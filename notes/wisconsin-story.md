# Wisconsin story: a year of a northern landscape at L-band (`edoTenseii`)

Pages for the Wisconsin plate, written for the atlas page format (`web/README.md`, "Atlas pages"). Every claim here is a
test that passed in `notes/experiments.md` (sections W, WS, TW). Failed tests are told as failures. The layer and data
names below are what `experiments/ws_tests.py` writes into `data/processed/wisconsin/ws/` (not in git; regenerate with
the commands at the end).

Window: northern Wisconsin, W −91.00, S 45.02, E −90.56, N 45.28 (forest 67 %, fields 20 %, lakes, swamps).
Radar: NISAR GCOV (band B 80 m HH/HV, band A 20 m full quad-pol), GUNW, four tracks (dawn A76, A148; dusk D41, D113),
Oct 2025 – Nov 2026 (no data Feb – mid-June 2026).

## Pages

1. **intro: "A year in the north woods, heard at 24 cm"**
   - en: NISAR's L-band wave is 24 cm long. It passes through leaves and is scattered by branches, stalks and the ground.
     Here it watches a year of forest, fields, swamps and lakes.
   - uk: Хвиля NISAR L-діапазону має довжину 24 см. Вона проходить крізь листя і розсіюється на гілках, стеблах і ґрунті.
     Тут вона спостерігає рік лісу, полів, боліт і озер.
   - map: classes (`story_cells.json`, column `class`), time slider over all dates.
   - facts: 52 dates · 4 tracks · 24 cm.

2. **what: "Leaves are invisible" (W2, W3)**
   - en: When the maples drop their leaves, the forest's cross-polarised echo hardly moves (0.8 – 2.8 dB over the year). The wave
     is scattered by branches comparable to its length; a leaf is too thin to matter.
   - uk: Коли клени скидають листя, крос-поляризоване відлуння лісу майже не змінюється (0.8 – 2.8 дБ за рік). Хвилю
     розсіюють гілки, порівнянні з її довжиною, а лист для неї надто тонкий.
   - map: `story_cells.json` mode `hv_anom_db` (time), forest cells stay pale through October.

3. **what: "Corn stands up" (WS1b; WS1a failed)**
   - en: In late October 2025, corn fields were 3 – 8 dB brighter in HV than soybean fields: standing stalks and tall stubble make
     a volume, soybean stubble does not. The gap shrank week by week as the corn was harvested. The co-pol ratio test (HH/VV)
     passed on only 3 of 5 dates.
   - uk: Наприкінці жовтня 2025 року поля кукурудзи були на 3 – 8 дБ яскравіші в HV, ніж поля сої: стебла й висока стерня
     утворюють об'єм, а соєва стерня ні. Різниця танула щотижня в міру збирання кукурудзи. Тест співвідношення HH/VV пройшов
     лише на 3 датах із 5.
   - facts: +8.4 dB (22 Oct) → +3.0 dB (20 Nov).

4. **what: "Harvest, seen twice" (WS4)**
   - en: Sentinel-2 sees a field turn yellow; the radar sees its stalks cut. For 93 % of the 1 km crop blocks, the two harvest dates
     agree within 24 days, the radar about two weeks after the optical drop.
   - uk: Sentinel-2 бачить, як поле жовтіє, а радар бачить, як стебла зрізають. Для 93 % кілометрових блоків полів дві дати
     жнив збігаються в межах 24 днів; радар фіксує жнива приблизно на два тижні пізніше за оптику.
   - map: `story_harvest.json` (radar date per block); play Aug – Oct 2026.

5. **what: "Rain wakes the bare fields" (WS5)**
   - en: After harvest, bare soil is a mirror whose brightness follows its water. Rain in the three days before a pass raises the
     fields' HH echo (ρ = 0.64): wet soil has a higher permittivity.
   - uk: Після жнив голий ґрунт працює як дзеркало, яскравість якого залежить від води. Дощ за три дні до прольоту
     піднімає відлуння полів у HH (ρ = 0.64), бо вологий ґрунт має вищу діелектричну проникність.

6. **what: "The forest freezes" (W4) and "snow, not frost, blurs the winter" (WS6 failed)**
   - en: When trees and soil freeze, water turns to ice (ε 80 → 3) and HH drops by about 2 dB. We expected frozen ground to keep the
     interferometric fingerprint; it did not. Every frozen pair was a snowy pair, and the forest lost coherence under fresh snow.
   - uk: Коли дерева й ґрунт замерзають, вода стає льодом (ε 80 → 3), і HH падає приблизно на 2 дБ. Ми чекали, що мерзла земля
     збереже інтерферометричний «відбиток», але цього не сталося: кожна мерзла пара була сніжною, і ліс під свіжим снігом
     втрачав когерентність.
   - map: `story_cells.json` mode `hh_anom_db`, Dec – Jan.

7. **why: "The radar measures the sky twice" (WS8, WS9a; WS9b failed)**
   - en: On its way down and back, the wave crosses the ionosphere. Its polarisation turns (Faraday rotation), and its phase is
     delayed differently at the band's two ends. The two measurements agree (ρ = 0.59 on 45 pairs, and their ratio is
     within ~10 % of theory). The rotation also follows the global TEC maps (ρ = 0.81): up to 5° on summer evenings, less at
     dawn. In the 2026 processing the sign agrees too (32 of 33 scenes). The 2025 releases did not.
   - uk: Дорогою вниз і назад хвиля двічі проходить іоносферу. Її поляризація повертається (ефект Фарадея), а фаза
     запізнюється по-різному на двох краях діапазону. Обидва виміри узгоджуються (ρ = 0.59 на 45 парах, а їхнє
     співвідношення в межах ~10 % від теорії). Поворот іде і за глобальними картами TEC (ρ = 0.81): до 5° влітку ввечері,
     менше на світанку. У обробці 2026 року збігається й знак (32 з 33 знімків), у релізах 2025 року — ні.
   - map: a bars strip per date: measured Ω and predicted Ω (`ws_results.json` → `WS9.points`).

8. **sure: "How sure, and what failed"**
   - en: Failed and stated: early lake ice is dark, not bright (WS3); no frost heave above the noise (WS7); swamp forest is not
     brighter in HH (WS2a). The product's co-pol phase is not calibrated (forest at −120° or +140° instead of 0°); 2025
     releases have a ~1 dB HV/VH imbalance (W5a) and a Faraday offset; the noise layer is often broken (W1a). Summer dusk
     HV of the forest repeats within 0.1 dB (TW).
   - uk: Провали чесно: ранній лід на озерах темний, а не яскравий (WS3); морозного здимання понад шум не видно (WS7); болотний ліс
     не яскравіший у HH (WS2a). Фаза між HH і VV у продукті не відкалібрована (ліс на −120° або +140° замість 0°); у релізах
     2025 року є дисбаланс HV/VH ~1 дБ (W5a) і зсув фарадеївського кута; шар шуму часто зламаний (W1a). Вечірній HV лісу
     влітку повторюється з точністю 0.1 дБ (TW).

9. **next: monthly polarimetry**
   - en: Entropy, alpha and Pauli colour of the fields month by month (`ws10_SUB2.png`): fields switch from surface to volume
     scattering and back.
   - uk: Ентропія, кут α і колір Паулі полів місяць за місяцем (`ws10_SUB2.png`): поля перемикаються з поверхневого
     розсіювання на об'ємне і назад.

## Data for the plate (regenerate)

```bash
uv run python experiments/w_fetch.py        # band B windows, class grid (section W)
uv run python experiments/ws_refs.py        # CDL 2025, NWI, SUB2
uv run python experiments/ws_fetch.py weather
uv run python experiments/ws_fetch.py gim
uv run python experiments/ws_fetch.py gcov
uv run python experiments/ws_fetch.py gunw
uv run python experiments/ws_fetch.py s2
uv run python experiments/ws_tests.py       # -> data/processed/wisconsin/ws/
uv run python experiments/ws_plot.py
```

| file | what | for |
|---|---|---|
| `story_cells.json` | 1 km blocks: lon, lat, class, HH and HV per date relative to the block's own summer on its own track; dates, tracks, freeze state | a `cells` layer with `time: true` modes |
| `story_harvest.json` | per crop block: radar harvest date per track, optical date | a `cells` or `polygons` layer by date |
| `ws_results.json` | every per-date row (band B, band A, GUNW) and every test | bars strips (Ω, coherence, class medians) |
| `ws_tests.png`, `ws10_SUB.png`, `ws10_SUB2.png` | figures | the "sure" and "next" pages |
