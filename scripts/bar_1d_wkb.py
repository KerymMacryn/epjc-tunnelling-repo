import argparse, numpy as np, pandas as pd
from pathlib import Path

def rectangular_transmission(V0, L, E):
    if E >= V0:
        return 1.0
    kappa = np.sqrt(V0 - E)
    # Transmisión barrera rectangular (ħ^2/2m = 1)
    T = 1.0 / (1.0 + (V0**2 * np.sinh(kappa*L)**2) / (4*E*(V0-E)))
    return float(T)

def agmon_action_rectangular(V0, L, E):
    return np.sqrt(max(V0 - E, 0.0)) * L

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--potential", default="rectangular", choices=["rectangular"])
    ap.add_argument("--V0", type=float, required=True)
    ap.add_argument("--L", type=float, required=True)
    ap.add_argument("--E", type=float, required=True)
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    if args.potential != "rectangular":
        raise SystemExit("Only --potential rectangular is implemented.")

    S = agmon_action_rectangular(args.V0, args.L, args.E)
    T_exact = rectangular_transmission(args.V0, args.L, args.E)
    T_wkb = np.exp(-2*S)

    df = pd.DataFrame([{
        "E": args.E,
        "S_Agmon": S,
        "T_exact": T_exact,
        "T_WKB": T_wkb,
        "rel_error": abs(T_wkb - T_exact)/max(T_exact, 1e-16)
    }])

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)  # <-- crea artifacts/data si no existe
    df.to_csv(out_path, index=False)
    print(f"Wrote {out_path.as_posix()}")
