"""D: fault parameters from one NISAR co-seismic pair, by simulation-based inference (pre-registered).

Simulate on the GPU (Okada + real pre-event noise) -> train a CNN + mixture-density posterior on the
fly -> apply to the real pair. Classic multi-start least squares as the baseline. Blind: does not
read the USGS solution (scripts/sbi_compare.py does that afterwards).

    uv run --group ml python scripts/sbi_fault.py --event nagqu
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from rasterio.warp import transform as warp_transform

from hachure import nisar, okada

ROOT = Path(__file__).resolve().parents[1]
EVENTS = {
    "nagqu": dict(dir="nagqu-2026", co="co_*_112_A_*.npz", pre="pre_*_112_A_*.npz", epi=(93.7837, 32.3376)),
    "kumamoto": dict(dir="kumamoto-2026", co="co_*_024_A_*.npz", pre="pre_*_024_A_*.npz", epi=(130.6601, 32.6177)),
    "panama": dict(dir="panama-2026", co="co_*_076_A_*.npz", pre="pre_*_076_A_*.npz", epi=(-80.757, 7.607)),
    # D-3w diagnostic (unblinded): 2x wider window, 2 km blocks, magnitude-scaled length prior (env SBI_* below)
    "panama-wide": dict(dir="panama-2026-wide", co="co_*_076_A_*.npz", pre="pre_*_076_A_*.npz", epi=(-80.757, 7.607)),
    "kumamoto-desc": dict(dir="kumamoto-2026", co="co_*_003_D_*.npz", pre="pre_*_003_D_*.npz", epi=(130.6601, 32.6177),
                          heading=-168.0),
    "panama-wide-desc": dict(dir="panama-2026-wide", co="co_*_040_D_*.npz", pre="pre_*_040_D_*.npz", epi=(-80.757, 7.607),
                             heading=-169.0),
}
HEADING = -12.0          # NISAR ascending ground-track heading near 32 N (deg); left-looking
BLOCK = int(os.environ.get("SBI_BLOCK", 12))      # 12 x 80 m ~ 1 km (D-3w diagnostic: 24)
L_MAX = float(os.environ.get("SBI_LMAX", 80e3))    # length prior upper bound, m (D-3w: 200 km)
XY_MAX = float(os.environ.get("SBI_XYMAX", 25e3))  # centroid prior half-range, m (D-3w: 60 km)
K = 8                    # mixture components
dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
OUT = ROOT / "data" / "processed" / "events"


# ---------------------------------------------------------------- data
def los_map(path):
    r = nisar.RawWindow.load(path)
    ph = r.phase - (r.iono if r.iono is not None else 0)
    ok = np.isfinite(ph) & (r.comp > 0)
    return np.where(ok, nisar.LOS_SIGN * r.wavelength / (4 * np.pi) * ph, np.nan), r


def blocks(a):
    H, W = (a.shape[0] // BLOCK) * BLOCK, (a.shape[1] // BLOCK) * BLOCK
    v = a[:H, :W].reshape(H // BLOCK, BLOCK, W // BLOCK, BLOCK)
    n = np.isfinite(v).sum((1, 3))
    with np.errstate(all="ignore"):
        m = np.nanmean(v, (1, 3))
    m[n < 0.3 * BLOCK * BLOCK] = np.nan
    return m


def load_event(name):
    cfg = EVENTS[name]
    d = ROOT / "data" / "raw" / "nisar_events" / cfg["dir"]
    co = sorted(d.glob(cfg["co"]))[0]
    los, r = los_map(co)
    obs = blocks(los)
    inc = blocks(r.incidence.astype(float))
    H, W = obs.shape
    bx = r.x[: W * BLOCK].reshape(-1, BLOCK).mean(1)
    by = r.y[: H * BLOCK].reshape(-1, BLOCK).mean(1)
    ex, ey = warp_transform("EPSG:4326", f"EPSG:{r.epsg}", [cfg["epi"][0]], [cfg["epi"][1]])
    E, N = np.meshgrid(bx - ex[0], by - ey[0])
    noise = [blocks(los_map(p)[0]) for p in sorted(d.glob(cfg["pre"]))]
    print(f"{name}: {co.name[:62]}..., grid {obs.shape}, valid {np.isfinite(obs).mean():.2f}, "
          f"{len(noise)} noise maps, incidence {np.nanmedian(inc):.1f}", flush=True)
    return dict(obs=obs, inc=inc, E=E, N=N, noise=noise, co=co.name)


# ---------------------------------------------------------------- simulator (all on the GPU)
class Sim:
    def __init__(self, ev):
        self.mask = np.isfinite(ev["obs"]).ravel()
        self.shape = ev["obs"].shape
        self.m = torch.tensor(self.mask, device=dev)
        f = lambda a: torch.tensor(a.ravel()[self.mask], dtype=torch.float32, device=dev)   # noqa: E731
        self.e, self.n = f(ev["E"]), f(ev["N"])
        th = torch.deg2rad(f(np.where(np.isfinite(ev["inc"]), ev["inc"], np.nanmedian(ev["inc"]))))
        az = math.radians(HEADING + 90)                    # ground -> satellite azimuth (left-looking)
        self.lv = torch.stack([torch.sin(th) * math.sin(az), torch.sin(th) * math.cos(az), torch.cos(th)])
        self.A = torch.stack([torch.ones_like(self.e), self.e / 3e4, self.n / 3e4], 1)
        self.Apinv = torch.linalg.pinv(self.A)
        nz = []
        for nm in ev["noise"]:
            v = np.nan_to_num(nm - np.nanmedian(nm))
            for flip_ud in (False, True):
                for flip_lr in (False, True):
                    a = v[::-1] if flip_ud else v
                    a = a[:, ::-1] if flip_lr else a
                    nz.append(np.ascontiguousarray(a).ravel()[self.mask])
        self.noise = torch.tensor(np.array(nz), dtype=torch.float32, device=dev)
        self.coords = torch.tensor(np.stack([ev["E"], ev["N"]]) / 3e4, dtype=torch.float32, device=dev)
        self.obs = f(ev["obs"])

    def prior(self, B):
        u = lambda lo, hi: lo + (hi - lo) * torch.rand(B, device=dev)          # noqa: E731
        lu = lambda lo, hi: torch.exp(u(math.log(lo), math.log(hi)))            # noqa: E731
        L, W, dip = lu(5e3, L_MAX), lu(3e3, 30e3), u(10, 89)
        depth = torch.maximum(u(2e3, 20e3), W / 2 * torch.sin(torch.deg2rad(dip)) + 200)
        return dict(e0=u(-XY_MAX, XY_MAX), n0=u(-XY_MAX, XY_MAX), depth=depth, strike=u(0, 360), dip=dip,
                    rake=u(-180, 180), L=L, W=W, slip=lu(0.1, 8.0))

    def los(self, p):
        e = self.e[None, :] - p["e0"][:, None]
        n = self.n[None, :] - p["n0"][:, None]
        uE, uN, uZ = okada.displacement(e, n, p["depth"], p["strike"], p["dip"], p["L"], p["W"], p["rake"], p["slip"])
        return uE * self.lv[0] + uN * self.lv[1] + uZ * self.lv[2]

    def deplane(self, y):
        return y - (y @ self.Apinv.T) @ self.A.T

    def batch(self, B):
        p = self.prior(B)
        k = torch.randint(0, len(self.noise), (B,), device=dev)
        scale = (0.7 + 0.8 * torch.rand(B, 1, device=dev)) * torch.where(torch.rand(B, 1, device=dev) < 0.5, -1.0, 1.0)
        y = self.deplane(self.los(p) + scale * self.noise[k] + self.turbulence(B))
        return self.image(y), targets(p), p

    def turbulence(self, B):
        """Power-law (Kolmogorov-like, P ~ k^-8/3) random screens, std 0.3-3 cm, on the block grid."""
        H, W = self.shape
        ky = torch.fft.fftfreq(H, device=dev)[:, None]
        kx = torch.fft.rfftfreq(W, device=dev)[None, :]
        k = torch.sqrt(kx ** 2 + ky ** 2).clamp(min=1.0 / max(H, W))
        amp = k ** (-8 / 6)
        z = torch.randn(B, H, W // 2 + 1, device=dev) + 1j * torch.randn(B, H, W // 2 + 1, device=dev)
        f = torch.fft.irfft2(z * amp, s=(H, W)).reshape(B, -1)[:, self.m]
        f = f / f.std(1, keepdim=True)
        std = torch.exp(math.log(0.003) + (math.log(0.03) - math.log(0.003)) * torch.rand(B, 1, device=dev))
        return f * std

    def image(self, y):
        B = y.shape[0]
        img = torch.zeros(B, self.mask.size, device=dev)
        img[:, self.m] = torch.asinh(y / 0.05)
        mk = self.m.float().view(1, 1, *self.shape).expand(B, 1, *self.shape)
        return torch.cat([img.view(B, 1, *self.shape), mk, self.coords[None].expand(B, 2, *self.shape)], 1)


# ---------------------------------------------------------------- posterior network
NAMES = ["e0", "n0", "depth", "sin_st", "cos_st", "dip", "u1", "u2", "logL", "logW"]   # v2: signed slip components


def targets(p):
    st, rk = torch.deg2rad(p["strike"]), torch.deg2rad(p["rake"])
    return torch.stack([p["e0"] / 25e3, p["n0"] / 25e3, (p["depth"] / 1e3 - 10) / 6, torch.sin(st), torch.cos(st),
                        (p["dip"] - 50) / 23, torch.asinh(p["slip"] * torch.cos(rk) / 0.5),
                        torch.asinh(p["slip"] * torch.sin(rk) / 0.5), torch.log(p["L"] / 20e3) / 0.8,
                        torch.log(p["W"] / 10e3) / 0.7], 1)


def untarget(t):
    """-> physical units: km, deg, m."""
    return dict(e0=t[:, 0] * 25, n0=t[:, 1] * 25, depth=t[:, 2] * 6 + 10,
                strike=torch.rad2deg(torch.atan2(t[:, 3], t[:, 4])) % 360, dip=(t[:, 5] * 23 + 50).clamp(1, 90),
                rake=torch.rad2deg(torch.atan2(torch.sinh(t[:, 7]), torch.sinh(t[:, 6]))),
                L=torch.exp(t[:, 8] * 0.8) * 20, W=torch.exp(t[:, 9] * 0.7) * 10,
                slip=torch.hypot(torch.sinh(t[:, 6]), torch.sinh(t[:, 7])) * 0.5)


class Net(nn.Module):
    def __init__(self, D=len(NAMES)):
        super().__init__()
        c = lambda i, o: nn.Sequential(nn.Conv2d(i, o, 3, 2, 1), nn.GELU(), nn.Conv2d(o, o, 3, 1, 1), nn.GELU())  # noqa
        self.enc = nn.Sequential(c(4, 32), c(32, 64), c(64, 128), c(128, 128), nn.AdaptiveAvgPool2d(4), nn.Flatten())
        self.mlp = nn.Sequential(nn.Linear(128 * 16, 512), nn.GELU(), nn.Linear(512, 512), nn.GELU())
        self.D = D
        self.head = nn.Linear(512, K * (1 + 2 * D))

    def forward(self, x):
        h = self.head(self.mlp(self.enc(x)))
        logit, mu, ls = h[:, :K], h[:, K:K + K * self.D], h[:, K + K * self.D:]
        return logit, mu.view(-1, K, self.D), ls.view(-1, K, self.D).clamp(-7, 3)


def nll(net, x, t):
    logit, mu, ls = net(x)
    lp = -0.5 * (((t[:, None] - mu) / ls.exp()) ** 2).sum(-1) - ls.sum(-1)
    return -(torch.logsumexp(torch.log_softmax(logit, -1) + lp, -1)).mean()


@torch.no_grad()
def sample(net, x, n):
    logit, mu, ls = net(x)
    k = torch.distributions.Categorical(logits=logit[0]).sample((n,))
    return mu[0, k] + ls[0, k].exp() * torch.randn(n, mu.shape[-1], device=dev)


# ---------------------------------------------------------------- baseline: multi-start least squares
def lsq(sim, starts=128, iters=1500):
    p0 = sim.prior(starts)
    p0["ztop"] = (p0.pop("depth") - p0["W"] / 2 * torch.sin(torch.deg2rad(p0["dip"]))).clamp(min=100)
    lo = dict(e0=-40e3, n0=-40e3, ztop=50, strike=-720, dip=5, rake=-720, L=2e3, W=1e3, slip=0.02)
    hi = dict(e0=40e3, n0=40e3, ztop=30e3, strike=720, dip=89.5, rake=720, L=150e3, W=60e3, slip=30)
    z = {k: torch.logit(((v - lo[k]) / (hi[k] - lo[k])).clamp(1e-4, 1 - 1e-4)).clone().requires_grad_(True)
         for k, v in p0.items()}
    opt = torch.optim.Adam(z.values(), lr=0.05)
    yo = sim.deplane(sim.obs[None])[0]
    for it in range(iters):
        p = {k: lo[k] + (hi[k] - lo[k]) * torch.sigmoid(v) for k, v in z.items()}
        p["depth"] = p.pop("ztop") + p["W"] / 2 * torch.sin(torch.deg2rad(p["dip"]))
        r = sim.deplane(sim.los(p)) - yo
        loss = torch.nan_to_num((r ** 2).mean(1), nan=0.0)
        opt.zero_grad()
        loss.sum().backward()
        for v in z.values():
            v.grad.nan_to_num_(0.0)
        opt.step()
    with torch.no_grad():
        p = {k: lo[k] + (hi[k] - lo[k]) * torch.sigmoid(v) for k, v in z.items()}
        p["depth"] = p.pop("ztop") + p["W"] / 2 * torch.sin(torch.deg2rad(p["dip"]))
        loss = torch.nan_to_num(((sim.deplane(sim.los(p)) - yo) ** 2).mean(1), nan=1e9)
        b = int(torch.argmin(loss))
    return {k: float(v[b]) for k, v in p.items()}, float(loss[b].sqrt()), float(yo.std())


def physical(p):
    q = dict(p)
    for k in ("e0", "n0", "depth", "L", "W"):
        q[k] = p[k] / 1e3
    q["strike"] %= 360
    q["rake"] = (p["rake"] + 180) % 360 - 180
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event", required=True, choices=list(EVENTS))
    ap.add_argument("--steps", type=int, default=8000)
    ap.add_argument("--batch", type=int, default=256)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    torch.manual_seed(0)
    ev = load_event(args.event)
    sim = Sim(ev)

    t0 = time.time()
    x, _, _ = sim.batch(1024)
    torch.cuda.synchronize()
    print(f"simulator: 1024 noisy scenes in {time.time() - t0:.2f} s on {dev}", flush=True)

    net = Net().to(dev)
    opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.OneCycleLR(opt, max_lr=1e-3, total_steps=args.steps)
    t0 = time.time()
    for step in range(args.steps):
        x, t, _ = sim.batch(args.batch)
        loss = nll(net, x, t)
        opt.zero_grad()
        loss.backward()
        nn.utils.clip_grad_norm_(net.parameters(), 5.0)
        opt.step()
        sched.step()
        if step % 500 == 0 or step == args.steps - 1:
            print(f"  step {step:5d}  nll {loss.item():+.3f}  ({time.time() - t0:.0f} s)", flush=True)
    train_s = time.time() - t0
    net.eval()

    # calibration on held-out synthetics
    cover = {k: [] for k in ("e0", "n0", "depth", "dip", "L", "slip", "strike")}
    with torch.no_grad():
        for _ in range(4):
            x, t, p = sim.batch(250)
            for i in range(250):
                s = untarget(sample(net, x[i:i + 1], 400))
                truth = untarget(t[i:i + 1])
                for k in cover:
                    if k == "strike":
                        dlt = (s[k] - truth[k] + 180) % 360 - 180
                        lo, hi = torch.quantile(dlt, 0.16), torch.quantile(dlt, 0.84)
                        cover[k].append(bool(lo <= 0 <= hi))
                    else:
                        lo, hi = torch.quantile(s[k], 0.16), torch.quantile(s[k], 0.84)
                        cover[k].append(bool(lo <= truth[k][0] <= hi))
    cal = {k: float(np.mean(v)) for k, v in cover.items()}
    with torch.no_grad():
        x, t, p = sim.batch(300)
        sense = []
        for i in range(300):
            s = untarget(sample(net, x[i:i + 1], 400))
            sense.append(float((((s["rake"] - p["rake"][i] + 180) % 360 - 180).abs() < 90).float().mean()))
    cal["rake_sense_mass_median"] = float(np.median(sense))
    print(f"rake: posterior mass on the correct slip sense, median {np.median(sense):.2f} (v1 network: 0.50)")
    print("calibration, 68% interval coverage on 1000 held-out synthetics:",
          ", ".join(f"{k} {v * 100:.0f}%" for k, v in cal.items()))

    # the real pair
    yo = sim.deplane(sim.obs[None])
    xo = sim.image(yo)
    torch.cuda.synchronize()
    t0 = time.time()
    s = untarget(sample(net, xo, 20000))
    torch.cuda.synchronize()
    inf_s = time.time() - t0
    post = {}
    for k, v in s.items():
        v = v.cpu().numpy()
        if k in ("strike", "rake"):
            c = math.degrees(math.atan2(np.mean(np.sin(np.radians(v))), np.mean(np.cos(np.radians(v))))) % 360
            d = (v - c + 180) % 360 - 180
            post[k] = [c if k == "strike" else (c + 180) % 360 - 180, float(c + np.percentile(d, 16)),
                       float(c + np.percentile(d, 84))]
        else:
            post[k] = [float(np.median(v)), float(np.percentile(v, 16)), float(np.percentile(v, 84))]
    mw = np.array([okada.moment_magnitude(L * 1e3, W * 1e3, sl) for L, W, sl in
                   zip(s["L"].cpu().numpy(), s["W"].cpu().numpy(), s["slip"].cpu().numpy())])
    post["Mw"] = [float(np.median(mw)), float(np.percentile(mw, 16)), float(np.percentile(mw, 84))]
    print(f"\nSBI posterior for the real pair (inference {inf_s * 1000:.0f} ms for 20000 samples):")
    for k, (m, lo, hi) in post.items():
        print(f"  {k:7s} {m:8.2f}   [{lo:8.2f}, {hi:8.2f}]")

    t0 = time.time()
    best, rms, ystd = lsq(sim)
    torch.cuda.synchronize()
    lsq_s = time.time() - t0
    bp = physical(best)
    bp["Mw"] = okada.moment_magnitude(best["L"], best["W"], best["slip"])
    print(f"\nleast squares (128 starts, {lsq_s:.0f} s): rms misfit {rms * 100:.1f} cm (data std {ystd * 100:.1f} cm)")
    print("  " + ", ".join(f"{k} {v:.2f}" for k, v in bp.items()))

    # misfit of the SBI posterior median model
    with torch.no_grad():
        pm = {k: torch.tensor([post[k][0] * (1e3 if k in ("e0", "n0", "depth", "L", "W") else 1)],
                              device=dev, dtype=torch.float32) for k in ("e0", "n0", "depth", "strike", "dip", "rake",
                                                                         "L", "W", "slip")}
        model = sim.deplane(sim.los(pm))[0]
        rms_sbi = float(((model - yo[0]) ** 2).mean().sqrt())
    print(f"SBI median model rms misfit {rms_sbi * 100:.1f} cm")
    res = dict(event=args.event, pair=ev["co"], train_s=train_s, inference_ms=inf_s * 1000, lsq_s=lsq_s,
               calibration=cal, sbi=post, sbi_rms_cm=rms_sbi * 100, lsq=bp, lsq_rms_cm=rms * 100,
               data_std_cm=ystd * 100)
    json.dump(res, open(OUT / f"sbi_{args.event}.json", "w"), indent=1)
    np.savez_compressed(OUT / f"sbi_{args.event}_maps.npz", obs=yo[0].cpu().numpy(), sbi=model.cpu().numpy(),
                        mask=sim.mask, shape=np.array(sim.shape), E=ev["E"], N=ev["N"])
    torch.save(net.state_dict(), OUT / f"sbi_{args.event}_net.pt")


if __name__ == "__main__":
    main()
