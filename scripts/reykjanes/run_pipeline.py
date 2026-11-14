"""Run the remaining analysis steps one after another, each as a top-level task on the
status page (scripts/task_server.py); every script's own tasks hang under its step.

Usage: python scripts/reykjanes/run_pipeline.py [step ...]   (default: all, in order)
"""
import subprocess
import sys
import time
from pathlib import Path

from progress import Progress

HERE = Path(__file__).resolve().parent
PY = sys.executable
LOGS = HERE.parents[1] / "output" / "reykjanes" / "ml"
STEPS = [
    ("cache", "Кеш спостережень Рейк'янесу (з сирих GUNW)", ["-c", "import ml_npe_source as m; m.build_cache()"]),
    ("fem_table", "МСЕ: таблиця функцій Гріна (шарувата кора)", ["ml_fem_inference.py", "table", "--workers", "3"]),
    ("fem_mcmc", "МСЕ: MCMC на емуляторі (з рельєфом і без)", ["ml_fem_inference.py", "mcmc"]),
    ("flow_mogi", "NPE-потік проти MCMC: Моґі", ["ml_npe_flow.py", "mogi", "--n1", "60000", "--n2", "100000", "--ep1", "20", "--ep2", "30"]),
    ("flow_ellipsoid", "NPE-потік проти MCMC: видовжене тіло", ["ml_npe_flow.py", "ellipsoid", "--n1", "40000", "--n2", "80000", "--ep1", "20", "--ep2", "30", "--k", "192"]),
    ("dike", "Дайка в реальному часі: навчання і 1000 подій", ["ml_dike_realtime.py", "--sims", "100000", "--epochs", "30", "--mcmc-steps", "1500"]),
    ("fem_npe", "МСЕ: NPE на прямих симуляціях", ["ml_fem_inference.py", "npe", "--sims1", "300", "--sims2", "600", "--per", "60", "--workers", "3"]),
]


def main():
    want = sys.argv[1:] or [s[0] for s in STEPS]
    steps = [s for s in STEPS if s[0] in want]
    title = {"fem_table": "Конвеєр Б: МСЕ (паралельно)"}.get(want[0], "Конвеєр А: NPE проти MCMC, дайка") if sys.argv[1:] else "Конвеєр"
    with Progress(title, total=len(steps), unit="крок") as top:
        for key, name, args in steps:
            top.set(note=f"зараз: {name}")
            step = Progress(name, unit="крок", parent=top, note=f"лог: {key}_log.txt")
            t0 = time.time()
            with open(LOGS / f"{key}_log.txt", "w", encoding="utf-8") as log:
                r = subprocess.run([PY, *args], cwd=HERE, env={**step.child_env(), "PYTHONIOENCODING": "utf-8"},
                                   stdout=log, stderr=subprocess.STDOUT)
            tail = (LOGS / f"{key}_log.txt").read_text(encoding="utf-8", errors="replace").strip().splitlines()[-1:]
            step.close("done" if r.returncode == 0 else "failed",
                       note=f"{(time.time() - t0) / 60:.1f} хв · " + (tail[0][:160] if tail else ""))
            top.update(1)


if __name__ == "__main__":
    main()
