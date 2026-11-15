"""Publish the committed web/ folder to the public GitHub Pages repository of the atlas.

    uv run python scripts/publish_site.py [--repo OWNER/NAME] [--ref HEAD] [--dry-run]

Only what is committed under web/ at --ref is published (git archive), so data/, notes/, CLAUDE.md and anything
uncommitted never leave this repository. The Pages repository gets the site at its root, a .nojekyll file and a short
README with the credits; its history is one commit per publish. Needs `git` with push access to --repo.
The site is then at https://<owner>.github.io/<name>/ (Pages: main branch, root folder).
"""
from __future__ import annotations

import argparse
import io
import re
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = "MaksymSokolnytskyi/edelweiss-atlas"
SKIP = {"README.md"}                              # web/README.md is for developers; the site repo gets its own


def git(*args, cwd=ROOT, out=False):
    r = subprocess.run(["git", *args], cwd=cwd, check=True, capture_output=True)
    return r.stdout if out else None


def site_readme():
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    m = re.search(r"^## Data and credits\n(.*?)(?=^## )", text, flags=re.S | re.M)
    credits = m.group(1).strip() if m else ""
    return ("# Edelweiss: an atlas of the moving ground\n\n"
            "NASA Space Apps Challenge 2026, challenge *Dancing with the SARs*. A team from Dnipro, Ukraine.\n\n"
            "Plates of surface change seen by the NISAR radar, drawn as 19th-century engraved hachures; each plate says how "
            "sure it is. Open the site and start from the title page.\n\n"
            "This repository only holds the built web pages. It is published from the team's working repository by "
            "`scripts/publish_site.py`; do not edit it by hand.\n\n"
            "## Data and credits\n\n" + credits + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo", default=REPO)
    ap.add_argument("--ref", default="HEAD")
    ap.add_argument("--dry-run", action="store_true", help="build the folder, show what would be pushed, push nothing")
    a = ap.parse_args()
    sha = git("rev-parse", "--short", a.ref, out=True).decode().strip()
    tar = tarfile.open(fileobj=io.BytesIO(git("archive", "--format=tar", a.ref, "web", out=True)))
    with tempfile.TemporaryDirectory() as tmp:
        site = Path(tmp) / "site"
        try:
            git("clone", "--depth", "1", f"https://github.com/{a.repo}.git", str(site), cwd=tmp)
        except subprocess.CalledProcessError:              # a new, empty repository
            site.mkdir()
            git("init", "-b", "main", cwd=site)
            git("remote", "add", "origin", f"https://github.com/{a.repo}.git", cwd=site)
        for p in site.iterdir():                           # the site is replaced as a whole
            if p.name != ".git":
                shutil.rmtree(p) if p.is_dir() else p.unlink()
        n = 0
        for m in tar.getmembers():
            if not m.isfile():
                continue
            rel = Path(m.name).relative_to("web")
            if str(rel) in SKIP:
                continue
            dst = site / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes(tar.extractfile(m).read())
            n += 1
        (site / ".nojekyll").write_text("", encoding="utf-8")
        (site / "README.md").write_text(site_readme(), encoding="utf-8", newline="\n")
        git("add", "-A", cwd=site)
        status = git("status", "--short", cwd=site, out=True).decode()
        size = sum(p.stat().st_size for p in site.rglob("*") if p.is_file() and ".git" not in p.parts)
        print(f"{n} files from web/ at {sha}, {size / 1e6:.1f} MB; {len(status.splitlines())} changed")
        if a.dry_run or not status.strip():
            print("nothing pushed" + (" (dry run)" if a.dry_run else " (no changes)"))
            return
        git("-c", "user.name=Edelweiss atlas", "-c", "user.email=noreply@users.noreply.github.com",
            "commit", "-q", "-m", f"Publish the atlas from {sha}", cwd=site)
        git("push", "-q", "origin", "HEAD:main", cwd=site)
        owner, name = a.repo.split("/")
        print(f"pushed; the site: https://{owner.lower()}.github.io/{name}/")


if __name__ == "__main__":
    main()
