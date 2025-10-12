
import argparse, numpy as np, pandas as pd
import mpmath as mp
import matplotlib.pyplot as plt

mp.mp.dps = 50

def hankel1(nu, z):
    return mp.besselj(nu, z) + 1j*mp.bessely(nu, z)

def find_hankel_zeros_lower(nu, n_roots=2, nmax=80):
    roots = []
    n = 1
    trials = 0
    while len(roots) < n_roots and trials < nmax:
        x0 = (n + nu/2 - 0.25) * mp.pi
        z0 = x0 - 0.25j - 0.10j*n
        z1 = x0 + 0.25 - 0.55j - 0.12j*n
        try:
            z = mp.findroot(lambda z: hankel1(nu, z), (z0, z1), tol=1e-20, maxsteps=100)
            if mp.im(z) < -1e-8 and all(abs(z - r) > 1e-2 for r in roots):
                roots.append(z)
        except:
            try:
                z0 = x0*1.01 - 0.2j - 0.08j*n
                z1 = x0*0.99 + 0.2 - 0.45j - 0.10j*n
                z = mp.findroot(lambda z: hankel1(nu, z), (z0, z1), tol=1e-20, maxsteps=100)
                if mp.im(z) < -1e-8 and all(abs(z - r) > 1e-2 for r in roots):
                    roots.append(z)
            except:
                pass
        n += 1
        trials += 1
    return roots

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--R", type=float, default=1.0)
    ap.add_argument("--m", type=int, nargs="+", default=[0,1])
    ap.add_argument("--alpha", type=float, nargs="+", default=[0.0,0.25])
    ap.add_argument("--xrange", type=float, nargs=2, default=[0.5, 12.0])
    ap.add_argument("--yrange", type=float, nargs=2, default=[-4.0, -0.05])
    ap.add_argument("--csv", type=str, required=True)
    ap.add_argument("--fig", type=str, default=None)
    args = ap.parse_args()

    rows = []
    for m in args.m:
        for a in args.alpha:
            nu = abs(m + a)
            zs = find_hankel_zeros_lower(nu, n_roots=2)
            for z in zs:
                xr, yi = float(mp.re(z)), float(mp.im(z))
                Gamma = -4.0 * xr * yi  # ħ^2/2μ=1
                rows.append({"m": m, "alpha": a, "nu": float(nu),
                             "Re(kR)": xr, "Im(kR)": yi, "Gamma": Gamma})
    df = pd.DataFrame(rows)
    Path(args.csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.csv, index=False)
    print(f"Wrote {args.csv}")

    if args.fig:
        fig = plt.figure(figsize=(6.8, 4.0))
        ax = fig.add_subplot(111)
        for (m,a), g in df.groupby(["m","alpha"]):
            ax.scatter(g["Re(kR)"], g["Im(kR)"], label=f"m={m}, alpha={a}")
        ax.axhline(0.0, lw=1, alpha=0.4)
        ax.set_xlim(args.xrange); ax.set_ylim(args.yrange)
        ax.set_xlabel(r"$\Re(kR)$"); ax.set_ylabel(r"$\Im(kR)$")
        ax.grid(True, alpha=0.3); ax.legend(loc="lower left", fontsize=8)
        Path(args.fig).parent.mkdir(parents=True, exist_ok=True)
        fig.tight_layout()
        fig.savefig(args.fig, dpi=180)
        plt.close(fig)
        print(f"Wrote {args.fig}")
