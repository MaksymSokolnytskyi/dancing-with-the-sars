"""Results of the Reykjanes / Kilauea study on one readable page, without the map.

Same numbers as earth_breathes.html (read from the same json/npz outputs), laid out as a report:
uplift and GPS check, source inference (MCMC vs NPE, with timings), layered-crust FE, real-time
dike, Kilauea tephra / lava and polarizations, data sources and limitations.

Writes output/reykjanes/results.html (self-contained: figures embedded).
Usage: python scripts/reykjanes/build_results_page.py
"""
import base64
import html
import io
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from build_breathing_app import ML, OUT, expensive_payload, flow_timings, methods_payload
from results_validation_section import validation_html
from results_v4_section import v4_html

KIL = OUT.parent / "kilauea"
MNAME = {"mcmc_screens": "MCMC, справжні екрани атмосфери (еталон)", "npe_flow": "NPE-потік (вибілений, 2 раунди)",
         "npe_raw": "NPE, справжні екрани", "npe_gp": "NPE, екрани за GP", "mcmc": "MCMC, гауссів шум",
         "npe_gauss2": "NPE, гауссів шум"}
ORDER = ["mcmc_screens", "npe_flow", "npe_raw", "npe_gp", "mcmc", "npe_gauss2"]
COL = {"mcmc_screens": "#1f5fa8", "npe_flow": "#d0632b", "npe_raw": "#3c9a5f", "npe_gp": "#8a63b8", "mcmc": "#7a8794", "npe_gauss2": "#b8a03a"}


def esc(s):
    return html.escape(str(s))


def img_file(p, alt):
    if not p.exists():
        return ""
    return f'<figure><img src="data:image/png;base64,{base64.b64encode(p.read_bytes()).decode()}" alt="{esc(alt)}"></figure>'


def img_fig(fig, alt):
    b = io.BytesIO()
    fig.savefig(b, format="png", dpi=110, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return f'<figure><img src="data:image/png;base64,{base64.b64encode(b.getvalue()).decode()}" alt="{esc(alt)}"></figure>'


def f(v):
    return f"{v:.0f}" if abs(v) >= 10 else f"{v:.2f}" if abs(v) < 1 else f"{v:.1f}"


def ci(s):
    return f'{f(s[2])} <span class="ci">[{f(s[1])}…{f(s[3])}]</span>' if s else "—"


def dur(s):
    if s is None:
        return "—"
    return f"{s:.1f} с" if s < 1 else f"{s:.0f} с" if s < 90 else f"{s / 60:.0f} хв" if s < 5400 else f"{s / 3600:.1f} год"


def table(head, rows, cls=""):
    h = "".join(f"<th>{c}</th>" for c in head)
    b = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<div class="tw"><table class="{cls}"><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>'


def methods_section(M, T):
    out = []
    for kind, title in (("mogi", "Точкове джерело (Моґі)"), ("ellipsoid", "Видовжений резервуар (еліпсоїд)")):
        for ds in ("NISAR + GNSS", "NISAR only"):
            rows = []
            for m in ORDER:
                S = (M.get(kind, {}).get(m) or {}).get(ds)
                if not S:
                    continue
                sm = S["summary"]
                c = S.get("calibration_real") or S.get("calibration_flow")
                cal = (f'{min(c["68"])}–{max(c["68"])} %' + ("*" if "calibration_flow" in S and "calibration_real" not in S else "")) if c else \
                    (f'R̂ {S["mcmc"]["rhat"]:.2f}' if S.get("mcmc") else "—")
                tm = (T.get(kind, {}).get(ds) or {}).get(m) or {}
                train = "не потрібне" if tm and not tm.get("train_s") else (f'{dur(tm["train_s"])} <span class="ci">({tm["n_sims"] // 1000} тис. сим.)</span>' if tm else "—")
                infer = f'{dur(tm["infer_s"])}' if tm else "—"
                r = [f'<span class="sw" style="background:{COL[m]}"></span>{MNAME[m]}', ci(sm["depth"]), ci(sm["volume"])]
                if kind == "ellipsoid":
                    r += [ci([2 * x for x in sm["half_length"]]), ci(sm["strike"]), f'{round(S["p_long"] * 100)} %' if S.get("p_long") is not None else "—"]
                r += [cal, train, infer]
                rows.append(r)
            if not rows:
                continue
            head = ["метод", "глибина, км", "об’єм, млн м³/рік"] + (["довжина 2a, км", "вісь, °", "P(&gt;3 км)"] if kind == "ellipsoid" else []) + \
                   ["68 % інтервал містить правду", "навчання (раз)", "новий набір даних"]
            out.append(f'<h3>{title} · {"NISAR + GPS" if ds == "NISAR + GNSS" else "лише NISAR"}</h3>' + table(head, rows))
    return "".join(out)


def methods_figure(M):
    """Posterior histograms of the key parameters, every method on one axis (NISAR + GNSS)."""
    panels = [("mogi", "depth", "Моґі: глибина, км", (2.5, 5.5)), ("mogi", "volume", "Моґі: об’єм, млн м³/рік", (8, 26)),
              ("ellipsoid", "depth", "еліпсоїд: глибина, км", (1.5, 6)), ("ellipsoid", "half_length", "еліпсоїд: довжина 2a, км", (0, 12)),
              ("ellipsoid", "strike", "еліпсоїд: вісь, ° від півночі", (0, 180)), ("ellipsoid", "volume", "еліпсоїд: об’єм, млн м³/рік", (5, 45))]
    idx = {"mogi": ("east", "north", "depth", "volume"), "ellipsoid": ("east", "north", "depth", "half_length", "strike", "volume", "sill_share")}
    fig, ax = plt.subplots(2, 3, figsize=(13, 6.4))
    for a, (kind, p, ttl, rng) in zip(ax.ravel(), panels):
        i = idx[kind].index(p)
        for m in ORDER:
            S = (M.get(kind, {}).get(m) or {}).get("NISAR + GNSS")
            if not S:
                continue
            v = np.array(S["samples"])[:, i] * (2 if p == "half_length" else 1)
            lw = 2.4 if m in ("mcmc_screens", "npe_flow") else 1.1
            a.hist(v, bins=60, range=rng, density=True, histtype="step", lw=lw, color=COL[m], label=MNAME[m])
        a.set_title(ttl, fontsize=10)
        a.set_yticks([])
        a.spines[["top", "right", "left"]].set_visible(False)
    ax[0, 0].legend(fontsize=7, loc="upper left", frameon=False)
    fig.suptitle("Апостеріорні розподіли параметрів джерела, NISAR + GPS: еталон MCMC (синій) і NPE-потік (помаранчевий) товсті", fontsize=11)
    fig.tight_layout()
    return img_fig(fig, "Апостеріорні розподіли глибини, об’єму, довжини й осі для всіх методів")


def main():
    vec = json.loads((OUT / "breathing_vectors.json").read_text(encoding="utf-8"))
    hl = vec["headline"]
    M = methods_payload()
    T = flow_timings(json.loads((ML / "timings.json").read_text(encoding="utf-8")))
    X = expensive_payload()
    fm = X["fem_mcmc"]["layered + topography"]["summary"]
    fn = X.get("fem_npe")
    mref = X["mogi_ref"]
    D, V = "depth km b.s.l.", "dV 1e6 m3/yr"
    ell = M["ellipsoid"]["mcmc_screens"]["NISAR + GNSS"]["summary"]

    # ---- GPS check
    g_rows = []
    for st, g in vec["gnss"].items():
        v = g.get("val")
        if v:
            g_rows.append([st, f'{v["nisar"]:.1f} ± {v["nisar_sigma"]:.1f}', f'{v["gnss"]:.1f} ± {v["gnss_sigma"]:.1f}', f'{v["diff"]:+.1f}', f'{v["z"]:+.1f}'])

    # ---- FE
    fe_rows = [["однорідний пружний півпростір (Моґі), MCMC", ci(mref[D]), ci(mref[V]), "~40 с на набір даних"]]
    for lab, name in (("layered + topography", "шари + рельєф ASTER GDEM, MCMC на емуляторі МСЕ"), ("layered, flat", "шари без рельєфу, MCMC на емуляторі")):
        s = X["fem_mcmc"][lab]
        fe_rows.append([name, ci(s["summary"][D]), ci(s["summary"][V]), f'{dur(s["seconds"])} + 91 розрахунок МСЕ для таблиці'])
    if fn:
        tm = fn["timing"]
        fe_rows.append(["шари + рельєф, NPE на прямих розрахунках МСЕ", ci(fn["summary"][D]), ci(fn["summary"][V]),
                        f'{fn["n_solves"]} розрахунків {dur(tm["sim1"] + tm["sim2"])}, навчання {dur(tm["train1"] + tm["train2"])}; нові дані {tm["infer"]:.1f} с'])
        cmp_rows = [[nm, f'{c["shift_sd"]:+.2f}σ', f'×{c["width_ratio"]:.2f}'] for nm, c in fn["compare"].items()]

    # ---- flow vs MCMC detail
    def flow_rows(kind):
        j = json.loads((ML / f"npe_flow_{kind}_gnss.json").read_text(encoding="utf-8"))
        return [[nm, ci(r["mcmc"]), ci(r["flow_round2"]), f'{r["flow_round2_shift_sd"]:+.2f}σ', f'×{r["flow_round2_width_ratio"]:.2f}']
                for nm, r in j["params"].items()], j

    fr_m, jm = flow_rows("mogi")
    fr_e, je = flow_rows("ellipsoid")

    # ---- dike
    d = X["dike"]
    dt, de = d["timing"], d["median_abs_error"]
    c68 = [round(v * 100) for v in d["calibration"]["0.68"]]
    c95 = [round(v * 100) for v in d["calibration"]["0.95"]]
    dike_names = ("схід, км", "північ, км", "простягання, °", "довжина, км", "верх, км", "висота, км", "розкриття, м")
    dike_rows = [[n, f"{e:.2f}", f"{a} %", f"{b} %"] for n, e, a, b in zip(dike_names, de.values(), c68, c95)]

    # ---- the REAL 10 Nov 2023 dike (Sentinel-1)
    s1 = X.get("s1")
    s1_html = ""
    if s1:
        q = lambda v: f"{v[1]:.1f} [{v[0]:.1f}–{v[2]:.1f}]"
        Ns, Ms = s1["npe"], s1["mcmc"]
        s1_rows = [["довжина, км", q(Ns["length km"]), q(Ms["length km"]), "~15"],
                   ["об’єм дайки, млн м³", q(s1["npe_volume"]), q(s1["mcmc_volume"]), "130–139"],
                   ["розкриття, м", q(s1["opening_m"]), q([10 ** v for v in Ms["log10 opening m"]]), "до ~8 (макс.)"],
                   ["простягання, °", q(Ns["strike deg"]), q(Ms["strike deg"]), "—"],
                   ["Свартсенгі dV, трек A016, млн м³", q(Ns["Mogi dV A016 1e6 m3"]), q(Ms["Mogi dV A016 1e6 m3"]), "−76…−82 (до 12.11)"],
                   ["глибина джерела Свартсенгі, км н.р.м.", q(Ns["Mogi depth km"]), q(Ms["Mogi depth km"]), "—"],
                   ["лог-правдоподібність у медіані", f"{s1['loglike']['npe']:.0f}", f"{s1['loglike']['mcmc']:.0f}", "—"],
                   ["час відповіді", f"{s1['timing']['infer_real']:.1f} с", dur(s1["timing"]["mcmc"]), "—"]]
        c68s = [round(v * 100) for v in s1["calibration"]["0.68"]]
        mo = s1.get("mogi_r1")
        s1_html = f"""
  <h3>Справжня дайка 10.11.2023 під Ґріндавіком (Sentinel-1, ASF HyP3)</h3>
  <p>NISAR ще не бачив дайки (дані з листопада 2025 р.), тому для перевірки на справжній події взято дві інтерферограми Sentinel-1, замовлені в ASF HyP3 (INSAR_GAMMA, 40 м): висхідний трек 16 (31.10 → 12.11.2023) і низхідний трек 155 (10.11 о 07:59 UTC, за години до дайки → 22.11). Дані зведено до клітинок 1 км ({s1['bins']['A016']} і {s1['bins']['D155']}), коваріацію шуму оцінено за залишками. Модель: дайка за Окадою (однорідний півпростір, однакове розкриття) + джерело під Свартсенгі в <b>шаруватій корі з рельєфом</b> (функції Гріна МСЕ, ASTER GDEM v3), окремий об’єм для кожного треку. NPE: два раунди по 100 тис. симуляцій, 256 головних компонент.</p>
  {table(["", "NPE", "MCMC", "опубліковано*"], s1_rows)}
  <p class="note">* Sigmundsson et al. 2024, Science 383:1228 (doi 10.1126/science.adn2838); EGU24-14829 — також у моделі однорідного пружного середовища. Медіана [16–84 %].</p>
  <p><b>Чесний підсумок.</b> NPE відповідає за {s1['timing']['infer_real']:.1f} с, на симуляціях відкалібрований ({min(c68s)}–{max(c68s)} % для 68-відсоткових інтервалів), і його довжина й об’єм близькі до опублікованих. Але MCMC знаходить значно кращу підгонку тієї самої моделі — і зі старту з апріорного розподілу, і зі старту з відповіді NPE: коротшу дайку на північ від Сундгнукура з удвічі меншим об’ємом. Отже одна прямокутна дайка з однаковим розкриттям погано описує справжню (розкриття змінне, до ~8 м; частина під морем і в зоні грабену, де немає когерентності), і на такій невідповідній моделі NPE і MCMC розходяться; збіг NPE з літературою тут не є доказом точності. Шарувата кора замість Моґі покращила узгодженість NPE з даними (лог-правдоподібність {f"{mo['loglike']['npe']:.0f}" if mo else '—'} → {s1['loglike']['npe']:.0f}). Наступний крок — дайка зі змінним розкриттям: для неї MCMC дорожчає в рази, а NPE — ні.</p>
  {img_file(ML / 's1_grindavik_dike.png', 'Справжня дайка 10.11.2023: дані Sentinel-1, модель NPE, залишки, можливі дайки')}"""

    # ---- polarization
    pz =json.loads((KIL / "polarization.json").read_text(encoding="utf-8"))
    longi = [i for i, s in enumerate(pz["span_days"]) if s >= 20]
    pol_rows = []
    for c in pz["channels"]:
        S = pz["summary"][c]
        sw = pz["per_pair"][c]["mean_change_dB"]["SW"]
        pol_rows.append([c, f'{np.median(pz["per_pair"][c]["noise_sd_dB"]):.2f}', " / ".join(f"{sw[i]:+.1f}" for i in longi),
                         f'{S["SW"]["detect_median"]:.0f} %', f'{S["caldera"]["cnr_abs_median"]:.1f} · {S["caldera"]["sign_consistent"]}/10',
                         f'{S["NE"]["detect_median"]:.0f} %'])
    tj = json.loads((KIL / "tephra_evidence.json").read_text(encoding="utf-8"))

    body = f"""
<header>
  <p class="kicker">NASA Space Apps 2026 · Dancing with the SARs · NISAR L-band</p>
  <h1>Земля дихає: результати</h1>
  <p class="lede">Під Свартсенгі (Рейк’янес, Ісландія) накопичується магма, і земля над нею піднімається. Тут — усі числові результати без карти: швидкість підняття, перевірка GPS, де і яке джерело, як швидко його можна оцінювати, і що NISAR бачить у радарній яскравості під час вивержень Кілауеа.</p>
  <div class="kpis">
    <div class="kpi"><b>+{hl['rate']:.1f} ± {hl['sigma']:.1f}</b><span>см/рік — найшвидше підняття (щонайменше +{hl['at_least']:.1f})</span></div>
    <div class="kpi"><b>{hl['area_km2']:.0f} км²</b><span>надійно піднімається</span></div>
    <div class="kpi"><b>{fm[D][2]:.1f} км</b><span>глибина джерела з шаруватою корою (Моґі: {mref[D][2]:.1f})</span></div>
    <div class="kpi"><b>~{2 * ell['half_length'][2]:.0f} км</b><span>довжина резервуару, вісь {ell['strike'][2]:.0f}° від півночі</span></div>
    <div class="kpi"><b>{dt['infer_one']:.2f} с</b><span>інверсія дайки після нової пари (MCMC ~{dur(dt['mcmc_one'])})</span></div>
  </div>
  <nav><a href="#uplift">Підняття і GPS</a><a href="#source">Джерело: MCMC проти NPE</a><a href="#fem">Шари кори</a><a href="#dike">Дайка за секунди</a><a href="#validation">Перевірка NPE</a><a href="#v4">Сліпий тест і незалежні перевірки</a><a href="#kilauea">Кілауеа і поляризації</a><a href="#data">Дані</a><a href="#limits">Обмеження</a></nav>
</header>

<section id="uplift">
  <h2>1. Підняття і незалежна перевірка GPS</h2>
  <p>Два треки NISAR (159A і 082D), листопад 2025 – жовтень 2026, розкладені на вертикаль і схід. Найшвидше точка піднімається на <b>{hl['rate']:.1f} ± {hl['sigma']:.1f} см/рік</b>, накопичено ~{hl['cum']:.0f} см. Порівняння з GPS-станціями, яких інверсія не бачила при оцінці швидкості:</p>
  {table(["станція", "NISAR, см/рік", "GPS, см/рік", "різниця", "z"], g_rows)}
  <p class="note">z = різниця / сумарна похибка; |z| &lt; 2 — збіг у межах похибок. NISAR — медіана 3 × 3 клітинок біля станції; жодна станція не використана ні у відліку, ні в поправках.</p>
  <p><b>Збіг є лише в центрі підняття:</b> SENG (Свартсенгі) — {g_rows[0][1] if g_rows else '—'} проти {g_rows[0][2] if g_rows else '—'} см/рік. На трьох віддалених станціях, де рух повільний, NISAR і GPS розходяться на 4–7 см/рік (|z| 2.7–4.0), і знак розбіжності різний. Причину ми не встановили: це можуть бути атмосферні залишки на великих відстанях, вибір точки відліку або локальні рухи біля станцій. Тому надійний висновок — швидкість у центрі підняття, а не дрібні швидкості на краях.</p>
</section>

<section id="source">
  <h2>2. Джерело: MCMC проти NPE</h2>
  <p>Та сама фізика (точкове джерело Моґі або видовжений резервуар), різні способи виводу. <b>MCMC</b> щоразу заново проганяє модель сотні тисяч разів; <b>NPE</b> (нейронна оцінка апостеріорного розподілу) навчається один раз, а потім відповідає за секунди. Еталон — MCMC зі справжньою коваріацією атмосферних екранів NISAR.</p>
  {methods_section(M, T)}
  <p class="note">«68 % інтервал містить правду»: як часто на не бачених при навчанні симуляціях 68-відсоткові інтервали мережі містять справжнє значення (ідеал 68 %); * — NPE-потік перевірено на симуляціях з параметричною моделлю справжніх екранів; для MCMC — діагностика збіжності R̂. Час — один ноутбук, 8 потоків.</p>
  {methods_figure(M)}
  <h3>NPE-потік проти MCMC параметр за параметром (NISAR + GPS)</h3>
  <div class="two">
    <div><h4>Моґі</h4>{table(["параметр", "MCMC", "NPE-потік", "зсув", "ширина"], fr_m)}</div>
    <div><h4>Видовжений резервуар</h4>{table(["параметр", "MCMC", "NPE-потік", "зсув", "ширина"], fr_e)}</div>
  </div>
  <p>Для Моґі потік конкурентний з MCMC: зсуви медіан ≤ 0.5σ, ширина ×1.01–1.14. Для резервуару ширина та сама, що в MCMC, довжина і об’єм збігаються, але глибина на ~0.4 км менша (−1.7σ), а вісь на ~6° ближча до півночі. Що зробило NPE конкурентним: дані вибілено тією самою коваріацією шуму, що в MCMC, і стиснено до головних компонент (шум стає N(0, I)); свіжий шум щоепохи; другий раунд симуляцій лише в околі відповіді першого. Для резервуару 48 компонент губили інформацію про глибину (48-ма була ще втричі сильніша за шум); 192 компоненти зменшили зсув глибини з 3.7σ до 1.7σ.</p>
</section>

<section id="fem">
  <h2>3. Дорога фізика: шарувата кора і рельєф</h2>
  <p>Осесиметрична скінченно-елементна модель (P2) пружної кори з профілем Vp за Růžek et al. 2023 до 2.5 км і типовою ісландською корою глибше (припущення), густина за Гарднером; рельєф NASA/METI ASTER GDEM v3 у наближенні змінної глибини.</p>
  {table(["модель і метод", "глибина нижче рівня моря, км", "об’єм, млн м³/рік", "час"], fe_rows)}
  <p><b>Висновок:</b> шари кори роблять джерело на ~{fm[D][2] - mref[D][2]:.1f} км глибшим і об’єм на ~{round((1 - fm[V][2] / mref[V][2]) * 100)} % меншим, ніж однорідний Моґі. Рельєф майже нічого не змінює ({fm[D][2] - X['fem_mcmc']['layered, flat']['summary'][D][2]:+.2f} км).</p>
  <h3>NPE на прямих розрахунках МСЕ проти MCMC на емуляторі</h3>
  {table(["параметр", "зсув медіани", "ширина інтервалу"], cmp_rows) if fn else ""}
  <p>MCMC прямо на МСЕ потребувало б ~256 000 розрахунків (~15 днів), тому для MCMC знадобився емулятор (таблиця функцій Гріна). NPE навчається на прямих розрахунках: кожен розрахунок задає форму поля для своєї глибини, а положення й об’єм її лише зсувають і масштабують, тож один розрахунок дає {fn['per_solve'] if fn else '—'} точних навчальних прикладів. Калібрування на не бачених розрахунках МСЕ: {', '.join(f'{round(v * 100)}' for v in fn['calibration']['0.68']) if fn and fn.get('calibration') else '—'} % (ідеал 68). NPE навчено ще з рельєфом Copernicus GLO-30 (різниця з ASTER у точках спостережень ~8 ± 11 м).</p>
</section>

<section id="dike">
  <h2>4. Дайка за секунди після нової пари NISAR</h2>
  <p>Сценарій «тисяч однакових задач»: кожні 12 днів приходить нова пара, і під час кризи потрібна відповідь одразу. Одна мережа, навчена раз ({dur(dt['simulate'] + dt['train'])} на {d['n_sims'] // 1000} тис. симуляцій дайки Окади), відповідає за <b>{dt['infer_one']:.2f} с</b> на подію; MCMC — ~{dur(dt['mcmc_one'])} на подію, ~{d['mcmc_for_all_events_hours']:.0f} год на {d['n_events']} подій. Перевірено на {d['n_events']} синтетичних дайках зі справжньою атмосферою NISAR (повернуті екрани).</p>
  {table(["параметр", "медіанна похибка", "у 68 % інтервалі", "у 95 % інтервалі"], dike_rows)}
  <p class="note">Події синтетичні: справжньої дайки в наших парах NISAR немає. На одній з двох контрольних подій MCMC застряг у хибному режимі, а NPE дав близьку відповідь.</p>
  {img_file(ML / 'dike_realtime.png', 'Інверсія дайки за секунди: приклад, 1000 подій, час')}
  {s1_html}
</section>

{validation_html(img_file, table)}

<section id="v4">
  <h2>6. Складніша кінематика, сліпий тест і незалежні перевірки</h2>
  {v4_html(img_file, table)}
</section>

<section id="kilauea">
  <h2>7. Кілауеа 2026: тефра, нова лава і поляризації</h2>
  <p>Епізоди фонтанування 41, 42, 43, 50 і 51 (24.01–16.07.2026), два треки (072D, 151A), знімки GCOV до і після кожного (100 м). Деформація: лопать 1–2 км на південний захід від Halemaʻumaʻu, −11.4 см в епізоді 42 і +8.4 см в епізоді 50; дно кальдери декорелює в усіх парах.</p>
  {img_file(KIL / 'episodes.png', 'Деформація і когерентність Кілауеа за епізодами')}
  <h3>Тефра і нова лава в радарній яскравості</h3>
  <p>Найтемніший напрям від жерла — на південний захід ({min(tj['darkest_sector_deg']):.0f}–{max(tj['darkest_sector_deg']):.0f}°; у 8 з 10 пар 187–262°), тобто за пасатом: тефра згладжує поверхню. Сильне затемнення — лише в довгих парах (−2.2 дБ, 24 дні; −1.2 дБ, 36 днів). Дно кальдери світлішає в 7 з 10 пар (+0.3…+0.9 дБ) — свіжа шорстка лава.</p>
  {img_file(KIL / 'tephra_evidence.png', 'Зміна радарної яскравості навколо Кілауеа після епізодів')}
  <h3>Які поляризації підкреслюють зміни</h3>
  <p>NISAR над Кілауеа знімає дві поляризації (DH): HH — розсіяння поверхнею, HV — об’ємне розсіяння. Комплексного члена HH·HV* у продуктах немає, тож повний поляриметричний розклад неможливий; порівнюємо HH, HV, відношення HV/HH і суму. Інтерферометрія (GUNW) у NISAR — лише HH.</p>
  {table(["канал", "шум фону, дБ", "тефра (SW): зміна в довгих парах, дБ", "тефра: виявлено клітинок", "кальдера: |CNR| · той самий знак", "контроль NE: виявлено"], pol_rows)}
  <ul>
    <li><b>Тефра — HH найпотужніший, HV/HH найчистіший.</b> HH падає на −2.5 і −1.3 дБ, HV у 2–7 разів слабше, тож відношення HV/HH <i>зростає</i> (+1.1 і +0.9 дБ): тефра пригнічує поверхневе розсіяння. Відношення має найменший шум (0.35 дБ проти ~0.5) і майже нуль на контрольному секторі NE, бо вологість і калібрування зсувають обидва канали разом.</li>
    <li><b>Нова лава в кальдері — HV.</b> Світлішає у 8 з 10 пар: шорстка лава додає об’ємного розсіяння.</li>
    <li><b>Час:</b> у 12-денних парах зміни слабкі (|CNR| &lt; 1); тефра помітна лише в парах 24–36 днів, що охоплюють кілька фонтанувань. Поріг виявлення — 5 % хибних тривог на фоні (&gt; 7 км від жерла).</li>
  </ul>
  {img_file(KIL / 'polarization.png', 'Зміни HH, HV і HV/HH навколо Кілауеа')}
</section>

<section id="data">
  <h2>8. Дані</h2>
  <ul>
    <li><b>NISAR L2 GUNW</b> (інтерферограми, HH) і <b>L2 GCOV</b> (яскравість, HH і HV) — NASA/ISRO через ASF DAAC / Earthdata.</li>
    <li><b>Sentinel-1</b> (ESA Copernicus) — інтерферограми листопада 2023 р., замовлені й оброблені в <b>ASF HyP3</b> (NASA ASF DAAC, INSAR_GAMMA): лише для перевірки NPE на справжній дайці.</li>
    <li><b>Sentinel-1, додаткові пари</b> через ASF HyP3: тиха пара жовтня 2023 р. (шум), інфляція липня 2024 р. і чотири події 2024 р. (08.02, 02.03, 16.03, 29.05) — для нульових тестів і крос-валідації.</li>
    <li><b>GNSS 2023</b> — Nevada Geodetic Laboratory (SENG, NYLA, VOGS, REYK; добові ряди IGS20, вилучено тренд за 2023 р.).</li>
    <li><b>NASA FIRMS</b> — теплові точки VIIRS 375 м і MODIS 1 км (архів стандартної обробки), 2023–2025 рр.</li>
    <li><b>NASA HLS</b> (Harmonized Landsat Sentinel-2, LP DAAC) — SWIR-знімки лави й жерел 2023–2025 рр.</li>
    <li><b>ASTER GDEM v3</b> (NASA/METI, LP DAAC) — рельєф для моделі МСЕ; NASADEM/SRTM до 63.9° пн. ш. не сягають.</li>
    <li><b>GPS</b> — Nevada Geodetic Laboratory (станції SENG, AUSV, GEVK, VOGC, NYLA); незалежна перевірка, не NASA.</li>
    <li><b>Швидкісний розріз кори</b> — Růžek et al. 2023 (EGUsphere-2023-2465); глибші шари — припущення.</li>
  </ul>
</section>

<section id="limits">
  <h2>9. Обмеження</h2>
  <ul>
    <li>На віддалених GPS-станціях (AUSV, GEVK, VOGC) NISAR розходиться з GPS на 4–7 см/рік; причину не встановлено.</li>
    <li>Атмосферний шум оцінено з самих даних NISAR (екрани з ~20 дат); коваріація параметрична (експоненційна + білий шум).</li>
    <li>NPE відповідає правильно лише для тих самих треків, геометрії і моделі шуму; для нового вікна чи моделі її треба перенавчати.</li>
    <li>Видовжений резервуар: NPE-потік досі дає глибину на ~0.4 км меншу за MCMC.</li>
    <li>Профіль кори глибше 2.5 км — припущення; шари роблять джерело глибшим, тож абсолютна глибина залежить від цього припущення.</li>
    <li>Дайки в тесті NISAR синтетичні. Справжню дайку 2023 р. перевірено на Sentinel-1 простою моделлю (одна прямокутна дайка з однаковим розкриттям, однорідний півпростір для самої дайки): NPE і MCMC тут розходяться, бо модель не описує змінного розкриття; коваріацію шуму оцінено із залишків, тож вона завищена.</li>
    <li>GNSS для дайки 2023 р. — лише 4 станції NGL (у статті ~20 станцій IMO); немає зсувів пікселів (pixel offsets), тож ближнє поле дайки (до ~8 м розкриття) ми не бачимо; маска грабену залежить від сліду, знайденого NPE.</li>
    <li>Похибку моделі оцінено на дайці 2023 р. і не масштабовано із сигналом, тому NPE v2 майже не бачить менших дайок 2024 р.; широке апріорне перевірено лише малою мережею.</li>
    <li>Поляризаційний аналіз — 10 пар над одним вулканом.</li>
  </ul>
</section>
"""
    css = """
:root{--bg:#f7f6f3;--card:#ffffff;--ink:#1d2329;--muted:#5d6873;--line:#dfe3e8;--accent:#1f5fa8;--accent2:#d0632b;--th:#eef1f5}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#14181c;--card:#1c2228;--ink:#e6e9ec;--muted:#9aa5b0;--line:#2e3740;--accent:#6ea8ef;--accent2:#f08a55;--th:#232b33}}
:root[data-theme="dark"]{--bg:#14181c;--card:#1c2228;--ink:#e6e9ec;--muted:#9aa5b0;--line:#2e3740;--accent:#6ea8ef;--accent2:#f08a55;--th:#232b33}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:16px/1.6 "Segoe UI",system-ui,-apple-system,sans-serif}
main{max-width:1080px;margin:0 auto;padding:24px 16px 64px}
header{padding:8px 0 16px}
.kicker{color:var(--muted);font-size:13px;letter-spacing:.04em;text-transform:uppercase;margin:0}
h1{font-size:34px;line-height:1.15;margin:6px 0 10px}
h2{font-size:23px;margin:0 0 10px}
h3{font-size:17px;margin:22px 0 8px}
h4{font-size:15px;margin:6px 0}
.lede{font-size:17px;max-width:780px}
section{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:20px 22px;margin:18px 0}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:10px;margin:16px 0}
.kpi{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px 14px}
.kpi b{display:block;font-size:24px;color:var(--accent);font-variant-numeric:tabular-nums}
.kpi span{font-size:13px;color:var(--muted)}
nav{display:flex;flex-wrap:wrap;gap:8px}
nav a{color:var(--accent);text-decoration:none;border:1px solid var(--line);border-radius:999px;padding:3px 12px;font-size:14px;background:var(--card)}
nav a:hover{border-color:var(--accent)}
.tw{overflow-x:auto;margin:6px 0}
table{border-collapse:collapse;width:100%;font-size:14px;font-variant-numeric:tabular-nums}
th{background:var(--th);text-align:left;font-weight:600;padding:7px 9px;border-bottom:1px solid var(--line);white-space:nowrap}
td{padding:6px 9px;border-bottom:1px solid var(--line);vertical-align:top}
tbody tr:hover{background:var(--th)}
.ci{color:var(--muted);font-size:12px;white-space:nowrap}
.sw{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:6px;vertical-align:middle}
.note{color:var(--muted);font-size:13px}
figure{margin:14px 0 4px}
figure img{width:100%;height:auto;border-radius:8px;background:#fff;border:1px solid var(--line)}
.two{display:grid;grid-template-columns:1fr 1fr;gap:16px}
@media (max-width:760px){.two{grid-template-columns:1fr}h1{font-size:27px}section{padding:16px 14px}}
ul{padding-left:20px}
li{margin:4px 0}
"""
    page = f"""<!doctype html>
<html lang="uk"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Земля дихає: результати</title><style>{css}</style></head>
<body><main>{body}</main></body></html>"""
    dst = OUT / "results.html"
    dst.write_text(page, encoding="utf-8")
    print("wrote", dst, f"{len(page) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
