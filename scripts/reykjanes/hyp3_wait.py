"""Wait for HyP3 jobs (live progress on the task page), then download and unzip the finished products.

Usage: python scripts/reykjanes/hyp3_wait.py data/sentinel1/jobs_submitted3.json
The Earthdata token (EARTHDATA_TOKEN or ~/.edl_token) is read and sent only to hyp3-api.asf.alaska.edu; products are public URLs.
"""
from _secrets import read_env
import json
import subprocess
import sys
import time
import zipfile
from pathlib import Path

from progress import Progress

ROOT = Path(__file__).resolve().parents[2]
S1 = ROOT / "data" / "sentinel1"


def main():
    jobs = json.loads(Path(sys.argv[1]).read_text())["jobs"]
    env = read_env()
    tok = env["EARTHDATA_API_KEY"].strip().strip('"')
    done = set()
    with Progress(f"ASF HyP3: {len(jobs)} інтерферограм (тиша, інфляція, події 2024)", total=len(jobs), unit="продукт") as p:
        while len(done) < len(jobs):
            states = {}
            for j in jobs:
                if j["job_id"] in done:
                    states[j["job_id"]] = "DOWNLOADED"
                    continue
                r = subprocess.run(["curl", "-s", "-H", f"Authorization: Bearer {tok}", f"https://hyp3-api.asf.alaska.edu/jobs/{j['job_id']}"],
                                   capture_output=True, text=True)
                d = json.loads(r.stdout)
                states[j["job_id"]] = d["status_code"]
                if d["status_code"] == "SUCCEEDED":
                    f = d["files"][0]
                    z = S1 / f["filename"]
                    if not z.exists():
                        subprocess.run(["curl", "-s", "-L", "-o", str(z), f["url"]], check=True)
                    if z.suffix == ".zip":                    # InSAR products are zips; autoRIFT gives a netCDF
                        with zipfile.ZipFile(z) as zf:
                            zf.extractall(S1)
                        z.unlink()                            # keep only the extracted files (disk space)
                    done.add(j["job_id"])
                    p.update(1)
                elif d["status_code"] == "FAILED":
                    print("FAILED", j["name"], j["job_id"], flush=True)
                    done.add(j["job_id"])
                    p.update(1)
            summary = {s: list(states.values()).count(s) for s in set(states.values())}
            p.set(note=", ".join(f"{k}: {v}" for k, v in sorted(summary.items())))
            if len(done) < len(jobs):
                time.sleep(60)
    del tok
    print("all done")


if __name__ == "__main__":
    main()
