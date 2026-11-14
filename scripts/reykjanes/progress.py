"""tqdm-like progress reporting that a web page can watch (scripts/task_server.py).

Every task is a small JSON file in output/progress/: name, done / total, unit, status,
start / update time, process id, a note, and its parent task. Sub-tasks can live in other
processes: a script started with the environment variable PROGRESS_PARENT=<task id> hangs
its tasks under that parent, so a pipeline step and everything it does form one tree.

    from progress import Progress, track
    with Progress("MCMC ellipsoid", total=3000, unit="step") as p:
        for i in range(3000):
            ...
            p.update(1, note=f"acceptance {acc:.2f}")
    for x in track(items, "FE simulations"):
        ...
"""
import json
import os
import time
import uuid
from pathlib import Path

DIR = Path(__file__).resolve().parents[2] / "output" / "progress"


class Progress:
    def __init__(self, name, total=None, unit="it", parent=None, note=""):
        DIR.mkdir(parents=True, exist_ok=True)
        if isinstance(parent, Progress):
            parent = parent.id
        self.id = f"{int(time.time() * 1000)}_{uuid.uuid4().hex[:6]}"
        self.file = DIR / f"{self.id}.json"
        self.d = dict(id=self.id, name=name, total=total, n=0, unit=unit, status="running", note=note,
                      parent=parent if parent is not None else os.environ.get("PROGRESS_PARENT"),
                      pid=os.getpid(), started=time.time(), updated=time.time())
        self._last = 0.0
        self._write(force=True)

    def _write(self, force=False):
        now = time.time()
        if not force and now - self._last < 0.4:
            return
        self._last = now
        self.d["updated"] = now
        tmp = self.file.with_suffix(".tmp")
        try:
            tmp.write_text(json.dumps(self.d, ensure_ascii=False), encoding="utf-8")
            os.replace(tmp, self.file)
        except OSError:
            pass

    def update(self, k=1, note=None):
        self.d["n"] += k
        if note is not None:
            self.d["note"] = note
        # throttled; but never skip the step that completes the task
        self._write(force=bool(self.d["total"]) and self.d["n"] >= self.d["total"])

    def set(self, n=None, total=None, note=None):
        if n is not None:
            self.d["n"] = n
        if total is not None:
            self.d["total"] = total
        if note is not None:
            self.d["note"] = note
        self._write(force=True)          # rare, explicit changes are never throttled away

    def child_env(self):
        """Environment for a subprocess whose tasks should hang under this one."""
        env = dict(os.environ)
        env["PROGRESS_PARENT"] = self.id
        return env

    def close(self, status="done", note=None):
        if self.d["total"] and status == "done":
            self.d["n"] = self.d["total"]
        self.d["status"] = status
        if note is not None:
            self.d["note"] = note
        self._write(force=True)

    def __enter__(self):
        return self

    def __exit__(self, et, ev, tb):
        self.close("failed" if et else "done", note=f"{et.__name__}: {ev}"[:200] if et else None)
        return False


def track(it, name, total=None, unit="it", parent=None):
    if total is None:
        try:
            total = len(it)
        except TypeError:
            total = None
    with Progress(name, total, unit, parent) as p:
        for x in it:
            yield x
            p.update(1)
