"""Secrets for the analysis scripts, read from the environment only (never from files in the repo).

EARTHDATA_TOKEN (or the file ~/.edl_token) -> returned as EARTHDATA_API_KEY, the key name the scripts use;
FIRMS_MAP_KEY -> NASA FIRMS area API key. Missing keys raise a clear error; values are never printed.
"""
import os
from pathlib import Path


def read_env():
    env = dict(os.environ)
    tok = env.get("EARTHDATA_TOKEN")
    f = Path.home() / ".edl_token"
    if not tok and f.exists():
        tok = f.read_text().strip()
    if tok:
        env["EARTHDATA_API_KEY"] = tok
    return _Env(env)


class _Env(dict):
    def __missing__(self, k):
        hint = "set EARTHDATA_TOKEN or write ~/.edl_token" if k == "EARTHDATA_API_KEY" else f"set the environment variable {k}"
        raise SystemExit(f"{k} not found: {hint}")
