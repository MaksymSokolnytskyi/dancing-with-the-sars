"""HTML section "NPE validation on the real dike" for results.html (reads the JSON outputs of the validation steps)."""
import json
from pathlib import Path

ML = Path(__file__).resolve().parents[2] / "output" / "reykjanes" / "ml"
OUT = ML.parent


def j(name):
    p = ML / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def validation_html(img_file, table):
    v2, ap, wide = j("s1_npe_v2.json"), j("s1_npe_v2_apply.json"), j("s1_npe_v2_wide.json")
    v3, ap3 = j("s1_npe_v3.json"), j("s1_npe_v2_apply_v3.json")
    dg, dn, dm = j("s1_dike_distributed.json"), j("s1_dike_distributed_nognss.json"), j("s1_dike_distributed_mogi.json")
    nq, nt, gn = j("s1_noise_quiet.json"), j("s1_noise_total.json"), json.loads((ML.parents[2] / "data" / "gnss" / "s1_gnss_2023.json").read_text())
    seng = json.loads((OUT / "seng_imo_compare.json").read_text()) if (OUT / "seng_imo_compare.json").exists() else None
    if not (v2 and ap and dg and dn):
        return ""
    f = lambda v: f"{v[1]:.0f} [{v[0]:.0f}–{v[2]:.0f}]"
    f5 = lambda v: f"{v[2]:.0f} [{v[1]:.0f}–{v[3]:.0f}]"
    mg, mi = v2["mcmc"]["gnss"], v2["mcmc"]["insar_only"]
    sb, pp = v2["sbc"], v2["ppc"]
    n_ok = sum(p > 0.01 for p in sb["p_uniform"].values())
    w68 = lambda v: v[3] - v[1]
    # ---- comparison table
    best = lambda d: min(d["depth_profile"], key=lambda p: p["chi2"])
    rows = ([["NPE v3 (рівень шуму невідомий, InSAR + GNSS)", f5(v3["real"]["volume"]), f(v3["real"]["pct"]["src depth km"]), f(v3["real"]["pct"]["src dV A016 1e6 m3"]),
              f"{v3['timing']['infer_real']:.0f} с"]] if v3 else []) + \
           [["NPE v2 (фіксований шум, InSAR + GNSS)", f5(v2["real"]["volume"]), f(v2["real"]["pct"]["src depth km"]), f(v2["real"]["pct"]["src dV A016 1e6 m3"]),
             f"{v2['timing'].get('infer_real', 0):.0f} с"],
            ["MCMC, та сама модель, InSAR + GNSS", f5(mg["volume"]), f(mg["pct"]["src depth km"]), f(mg["pct"]["src dV A016 1e6 m3"]), f"{mg['seconds'] / 60:.0f} хв"],
            ["MCMC, та сама модель, лише InSAR", f5(mi["volume"]), f(mi["pct"]["src depth km"]), f(mi["pct"]["src dV A016 1e6 m3"]), f"{mi['seconds'] / 60:.0f} хв"],
            ["розподілене розкриття, InSAR + GNSS, МСЕ", f"{dg['volume_1e6m3']:.0f}", f"{dg['geometry']['src_depth']:.1f}", f"{dg['src_dV_1e6m3']['A016']:.0f}", f"{dg['seconds_one_solve']:.1f} с"],
            ["розподілене розкриття, InSAR + GNSS, Моґі", f"{dm['volume_1e6m3']:.0f}" if dm else "—", f"{dm['geometry']['src_depth']:.1f}" if dm else "—",
             f"{dm['src_dV_1e6m3']['A016']:.0f}" if dm else "—", f"{dm['seconds_one_solve']:.1f} с" if dm else "—"],
            ["розподілене розкриття, лише InSAR, МСЕ", f"{dn['volume_1e6m3']:.0f}", f"{dn['geometry']['src_depth']:.1f} (слабко)", f"{dn['src_dV_1e6m3']['A016']:.0f}", f"{dn['seconds_one_solve']:.1f} с"],
            ["Sigmundsson et al. 2024 (~20 станцій GNSS, 4 супутники, зсуви пікселів)", "130–139", "4.0–4.2 (Моґі)", "−76…−82", "дні"]]
    comp = table(["метод", "об’єм дайки, млн м³", "глибина Свартсенгі, км", "dV Свартсенгі (до 12.11), млн м³", "час"], rows)
    # ---- events
    lab = {"quiet_2023": "тиша, 7/5 → 19/17.10.2023", "inflation_2024": "лише інфляція, 21/19.07 → 2.08/31.07.2024",
           "event_20240208": "дайка й виверження 08.02.2024", "event_20240302": "дайка без виверження 02.03.2024",
           "event_20240316": "виверження 16.03.2024", "event_20240529": "виверження 29.05.2024"}
    lit = {"quiet_2023": "немає дайки", "inflation_2024": "немає дайки; інфляція", "event_20240208": "—", "event_20240302": "~1.3 (з резервуару)",
           "event_20240316": "~17–23 (накопичення перед подією)", "event_20240529": "13–19 (дайка + 1-й тиждень)"}
    A3 = ap3 or ap
    ev_rows = [[lab[k], f"{ap[k]['p_nodike']:.2f} / {A3[k]['p_nodike']:.2f}", f"{f(ap[k]['volume'])} / {f(A3[k]['volume'])}",
                f"{A3[k]['pct']['src dV A016 1e6 m3'][1]:+.0f}", lit[k], ", ".join(f"{d:.1f}" for d in A3[k].get("vent_distance_km", [])) or "—"]
               for k in lab if k in ap and k in A3]
    ev = table(["пара", "P(без дайки): v2 / v3", "об’єм дайки, млн м³: v2 / v3", "dV Св. (v3), млн м³", "література", "до жерел HLS (v3), км"], ev_rows)
    # ---- hypotheses
    q_rms = ", ".join(f"{t} {v['rms_cm']:.1f}" for t, v in nq["tracks"].items()) if nq else "—"
    m_sill = ", ".join(f"{t} {v['sill2_cm2']:.0f}" for t, v in nt["tracks"].items()) if nt else "—"
    H = [["H1", "GNSS прибирає компроміс «об’єм дайки ↔ джерело Свартсенгі»",
          f"68 % інтервал V (MCMC): {w68(mi['volume']):.0f} → {w68(mg['volume']):.0f} млн м³ (−{100 * (1 - w68(mg['volume']) / w68(mi['volume'])):.0f} %); "
          f"глибина джерела визначена (Δχ² різкий); але corr(V, глибина) {mi['corr_V_depth']:+.2f} → {mg['corr_V_depth']:+.2f}; NPE–MCMC: "
          f"{v2['real']['volume'][2]:.0f} проти {mg['volume'][2]:.0f}", "частково"],
         ["H2", "розбіжність NPE ↔ MCMC — від моделі, а не від NPE",
          f"зі змінним розкриттям і маскою: медіани в межах 1σ; але χ² зразків NPE {pp['chi2_median']:.0f} проти {-2 * mg['loglike_median']:.0f} у MCMC, "
          f"залишок уздовж напрямків дайки {pp['dike_projection_median']:.0f} (шум ≤ {pp['dike_projection_noise_95']:.0f})", "ні"],
         ["H3", "коваріацію шуму завищено", f"тиха пара: rms {q_rms} см; похибка моделі (sill) {m_sill} см²; із сумою χ² {dg['chi2']:.0f} / {dg['n_data']}", "так"],
         ["H4", "NPE не вигадує дайок", f"v3: тиша P(без дайки) {A3['quiet_2023']['p_nodike']:.2f}; інфляція {A3['inflation_2024']['p_nodike']:.2f}, dV {A3['inflation_2024']['pct']['src dV A016 1e6 m3'][1]:+.0f} "
          f"(v2: {ap['quiet_2023']['p_nodike']:.2f} і {ap['inflation_2024']['p_nodike']:.2f})", "так (v3)" if A3['inflation_2024']['p_nodike'] >= 0.8 and A3['quiet_2023']['p_nodike'] >= 0.8 else "здебільшого"],
         ["H5", "NPE відкалібрований (SBC)", (f"v3 (випадковий рівень шуму): {sum(p > 0.01 for p in v3['sbc']['p_uniform'].values())} із {len(v3['sbc']['p_uniform'])} параметрів p > 0.01, об’єм p = {v3['sbc']['p_uniform_volume']:.2f}; "
          f"малі дайки 1–30 млн м³ за тихого шуму виявлено в {v3['sbc']['small_dikes']['detected_frac_quiet'] * 100:.0f} % випадків. " if v3 else "") +
          f"v2: {n_ok} із {len(sb['p_uniform'])}", "здебільшого"],
         ["H6", "сигнал не від апріорного", f"широке апріорне (мала мережа, {wide['n_sims'] // 1000 if wide else '—'} тис.): дайку бачить (P(без) {wide['real']['p_nodike'] if wide else '—'}), але простягання {wide['real']['pct']['strike deg'][1]:.0f}° [{wide['real']['pct']['strike deg'][0]:.0f}–{wide['real']['pct']['strike deg'][2]:.0f}]" if wide else "—", "не підтверджено (замало навчання)"],
         ["H7", "NPE переноситься на події 2024 р.", f"v3: 02.03 — {f(A3['event_20240302']['volume'])} (опубл. ~1.3, в інтервалі); 29.05 — {f(A3['event_20240529']['volume'])} проти 13–19 "
          f"(v2: {ap['event_20240529']['volume'][1]:.0f}); 08.02 і 16.03 — P(без дайки) ~{A3['event_20240208']['p_nodike']:.2f}, невпевнено", "частково (v3)"],
         ["H8", "шарувата кора узгоджується зі статтею", f"Моґі: {dm['geometry']['src_depth']:.1f} км, МСЕ: {dg['geometry']['src_depth']:.1f} км (на ~1 км глибше, як для NISAR); стаття — 4.0–4.2 км (Моґі)" if dm else "—", "так (зсув ~1 км)"]]
    hyp = table(["", "гіпотеза", "результат", "вердикт"], H)
    s = gn.get("SENG", {}).get("windows", {}).get("A016", {}).get("enu_cm", [0, 0, 0])
    return f"""
<section id="validation">
  <h2>5. Перевірка NPE на справжній дайці 10.11.2023</h2>
  <p>Перша інверсія дайки за Sentinel-1 дала розбіжність: NPE (12 км, 114 млн м³) і MCMC (8.8 км, 61 млн м³) з тією самою простою моделлю. Тому модель і перевірку посилено:</p>
  <ul>
    <li><b>дайка зі змінним розкриттям</b> (8 сегментів, 5 вузлів) і розподілене розкриття (патчі 1 × 1 км, лапласіан);</li>
    <li><b>маска грабену</b> ±1 км навколо сліду;</li>
    <li><b>GNSS:</b> 4 станції NGL з даними восени 2023 р. — SENG, NYLA, VOGS, REYK. SENG зсунувся на {s[0]:.0f} см на схід (тобто на захід) і {s[2]:.0f} см по вертикалі, як на рис. 3A статті. На VOGS знайдено й виключено апаратний стрибок 0.97 м;</li>
    <li><b>шум:</b> атмосфера з тихої пари плюс похибка моделі;</li>
    <li><b>клас «без дайки»:</b> 25 % симуляцій.</li>
  </ul>
  <p>Еталон — Sigmundsson et al. 2024 (Science 383:1228).</p>
  <h3>Гіпотези і вердикти</h3>
  {hyp}
  <h3>Об’єм дайки і джерело Свартсенгі різними методами</h3>
  {comp}
  <p class="note">Медіана [16–84 %]; для NPE і MCMC об’єм — медіана [16–84 %] з 95-відсоткового набору. dV Свартсенгі — вікно треку A016 (31.10 → 12.11). Час — один ноутбук.</p>
  <p><b>Що це означає.</b></p>
  <ul>
    <li><b>Лише InSAR:</b> об’єм дайки стійкий (~100–120 млн м³ за будь-якої глибини Свартсенгі), але джерело Свартсенгі майже не визначене.</li>
    <li><b>З GNSS:</b> глибина джерела визначається різко, dV Свартсенгі (−83…−94 млн м³) наближається до опублікованого, а об’єм дайки зростає до 150–170 млн м³.</li>
    <li><b>Збіг зі статтею:</b> за глибини джерела 4.0–4.2 км (як у статті) наш Моґі дає ~138 млн м³ — збіг з опублікованими 130–139. Глибше джерело в нас — наслідок того, що маємо лише 4 станції GNSS, і одна з них (SENG) стоїть за 2 км від дайки, де проста модель найслабша.</li>
  </ul>
  {img_file(ML / 's1_tradeoff.png', 'Компроміс об’єм дайки — джерело Свартсенгі з GNSS і без')}
  {img_file(ML / 's1_dike_distributed.png', 'Розподілене розкриття дайки, InSAR + GNSS, маска грабену')}
  <p class="note">Розподілене розкриття: максимум {dg['max_opening_m']:.1f} м на глибині ~{dg['depth_of_max_opening_km']:.1f} км. У статті — до ~8 м на 2–4 км. Різниця від того, що в нас немає даних ближнього поля: грабен замасковано, а зсувів пікселів (pixel offsets) немає. Згладжування (вибране перехресною перевіркою) розмазує максимум.</p>
  <h3>Чи не вигадує NPE дайку: нульові пари і події 2024 р. (без перенавчання)</h3>
  {ev}
  {img_file(ML / ('s1_null_tests_v3.png' if ap3 else 's1_null_tests.png'), 'Нульовий тест NPE v3')}
  {img_file(ML / ('s1_events_2024_v3.png' if ap3 else 's1_events_2024.png'), 'NPE v3 на подіях 2024 року')}
  {f'''<p><b>Виправлення (NPE v3).</b> Залишки 2023 р. не ростуть із сигналом: rms ~4 см і там, де деформація 0.4 см, і там, де 17 см. Отже, «похибка моделі» — це переважно шум саме цієї листопадової пари (зимова атмосфера плюс неописана деформація). Шум відрізняється між парами в ~10 разів: у тихій парі 0.3–0.4 см. Тому v3 навчено з невідомим рівнем шуму: кожна симуляція має свій рівень (0.1–40 см²) і кореляційну довжину (3–16 км), і мережа сама оцінює шум за даними.</p>
  <ul>
    <li><b>10.11.2023:</b> v3 дає {f5(v3['real']['volume'])} млн м³ (MCMC {f5(mg['volume'])}) і dV Свартсенгі {f(v3['real']['pct']['src dV A016 1e6 m3'])} млн м³ — як у статті (−76…−82).</li>
    <li><b>Нульові пари:</b> P(без дайки) {A3['quiet_2023']['p_nodike']:.2f} і {A3['inflation_2024']['p_nodike']:.2f}.</li>
    <li><b>Події 2024 р.:</b> 29.05 — {f(A3['event_20240529']['volume'])} млн м³ (v2: {f(ap['event_20240529']['volume'])}; опубл. 13–19).</li>
    <li><b>Що лишилося:</b> окремі зразки NPE пояснюють дані гірше за MCMC — χ² у медіані {v3['ppc']['chi2_at_median']:.0f} проти {-2 * mg['loglike_median']:.0f}. Тому NPE — швидка перша оцінка, а точну відповідь дає MCMC зі стартом від NPE.</li>
  </ul>''' if v3 and ap3 else ''}
  <p><b>Висновок про v2.</b> На тихій парі NPE відповідає «дайки немає» (P = {ap['quiet_2023']['p_nodike']:.2f}), а для інфляції дає додатне dV джерела. Але події 2024 р. NPE майже не бачить: P(без дайки) 0.5–0.7, а 29.05 він недооцінює. Причина — у моделі шуму. Похибку моделі (~31 см²) оцінено на величезній дайці 2023 р., і для слабких сигналів 2024 р. (χ² {ap['event_20240302']['chi2']:.0f}–{ap['event_20240529']['chi2']:.0f} на {ap['event_20240529']['n_data']} даних) шум завищено в рази. Похибка моделі має масштабуватися з сигналом — наступний крок. Відстані до жерел HLS малі ({', '.join(f'{d:.1f}' for d in ap['event_20240316'].get('vent_distance_km', []))} км), але слід лишається в коридорі 2023 р., тож це не доказ.</p>
  <h3>Калібрування (SBC), передбачення (PPC), апріорне</h3>
  {img_file(ML / ('s1_sbc_v3.png' if v3 else 's1_sbc.png'), 'SBC: ранги правди серед зразків NPE v3')}
  {img_file(ML / 's1_ppc.png', 'PPC: залишки для 10.11.2023')}
  {img_file(ML / 's1_prior_sensitivity.png', 'Чутливість до апріорного розподілу')}
  <p class="note">SBC: 500 симуляцій, з них 25 % без дайки; ранги правди серед 1000 зразків NPE.</p>
  <p class="note">PPC: 200 зразків NPE → пряма модель → вибілені залишки. Виходить, що NPE центрований правильно (медіана об’єму в межах 1σ від MCMC), але окремі зразки пояснюють дані гірше, ніж MCMC: розподіл NPE ширший за справжній.</p>
  <p class="note">Широке апріорне: мала мережа (30 тис. симуляцій, 12 епох) дайку виявляє, але не локалізує. Для чесної перевірки H6 потрібен двораундовий NPE.</p>
  <h3>Незалежні перевірки: оптика і GPS Свартсенгі</h3>
  {img_file(ML / 'hls_dike_optical.png', 'Оптична перевірка дайки: HLS SWIR, нова лава, жерла')}
  {img_file(ML / 'hls_timelapse.png', 'Часова серія SWIR кратерного ряду Сундгнукур 2024–2025')}
  {img_file(OUT / 'seng_imo_compare.png', 'SENG: GPS і NISAR, як на графіку Ісландського метеобюро')}
  <p class="note">SENG, 2025–2026: NISAR {seng['stats']['up']['nisar_end_mm']:.0f} мм вгору (до {seng['stats']['up']['nisar_last_date']}), GPS {seng['stats']['up']['gps_end_mm']:.0f} мм, IMO ~{seng['imo_end']['up']} мм (до 09.11.2026). Розбіжність NISAR − GPS: {seng['stats']['up']['rms_diff_mm']:.0f} мм rms.</p>
</section>
""" if seng else ""
