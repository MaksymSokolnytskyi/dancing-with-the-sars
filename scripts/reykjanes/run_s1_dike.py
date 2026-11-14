"""Driver for the Sentinel-1 Grindavik-dike demo with live progress on the task page (scripts/task_server.py).

Stages under one parent task: wait for the HyP3 jobs -> download -> bin -> NPE (simulations, training,
calibration, real-data posterior, MCMC). The NPE step runs as a child process whose own progress bars
(training epochs, MCMC steps) hang under the parent.
Usage: python scripts/reykjanes/run_s1_dike.py
"""
from _secrets import read_env
import json
import subprocess
import sys
import time
from pathlib import Path

from progress import Progress

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
S1 = ROOT / "data" / "sentinel1"


def job_status(tok, job_id):
    r = subprocess.run(["curl", "-s", "-H", f"Authorization: Bearer {tok}", f"https://hyp3-api.asf.alaska.edu/jobs/{job_id}"],
                       capture_output=True, text=True)
    return json.loads(r.stdout)["status_code"]


def main():
    jobs = json.loads((S1 / "jobs_submitted.json").read_text())["jobs"]
    env = read_env()
    tok = env["EARTHDATA_API_KEY"].strip().strip('"')
    with Progress("Справжня дайка 10.11.2023 (Sentinel-1 → NPE)", total=4, unit="етап") as top:
        with Progress("ASF HyP3: інтерферограми Sentinel-1", total=len(jobs), unit="робота", parent=top) as p:
            while True:
                st = {j["job_id"][:8]: job_status(tok, j["job_id"]) for j in jobs}
                done = sum(v == "SUCCEEDED" for v in st.values())
                p.set(n=done, note=", ".join(f"{k}: {v}" for k, v in st.items()))
                if any(v == "FAILED" for v in st.values()):
                    raise RuntimeError(f"HyP3 job failed: {st}")
                if done == len(jobs):
                    break
                time.sleep(60)
        del tok
        top.update(1, note="інтерферограми готові")
        for stage, note in (("download", "завантаження продуктів HyP3"), ("prep", "бінування 1 км, вектори огляду")):
            with Progress(note, total=1, unit="крок", parent=top) as p:
                subprocess.run([sys.executable, str(HERE / "s1_grindavik_dike.py"), stage], check=True, env=p.child_env(), cwd=HERE)
            top.update(1, note=note + ": готово")
        with Progress("NPE: симуляції, навчання, справжні дані, MCMC", total=1, unit="крок", parent=top) as p:
            subprocess.run([sys.executable, str(HERE / "s1_grindavik_dike.py"), "npe"], check=True, env=p.child_env(), cwd=HERE)
        top.update(1, note="готово")


if __name__ == "__main__":
    main()
