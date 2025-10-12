
import argparse, numpy as np, pandas as pd

def rectangular_transmission(V0, L, E):
    if E >= V0:
        k = np.sqrt(E - V0)
        return 1.0
    kappa = np.sqrt(V0 - E)
    # Textbook formula for rectangular barrier (units: ħ^2/2m = 1)
    T = 1.0 / (1.0 + (V0**2 * np.sinh(kappa*L)**2) / (4*E*(V0-E)))
    return float(T)

def agmon_action_rectangular(V0, L, E):
    kappa = np.sqrt(max(V0 - E, 0.0))
    return kappa * L

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--potential", default="rectangular", choices=["rectangular"])
    ap.add_argument("--V0", type=float, required=True)
    ap.add_argument("--L", type=float, required=True)
    ap.add_argument("--E", type=float, required=True)
    ap.add_argument("--out", type=str, required=True)
    args = ap.parse_args()

    S = agmon_action_rectangular(args.V0, args.L, args.E)
    T_exact = rectangular_transmission(args.V0, args.L, args.E)
    T_wkb = np.exp(-2*S)

    df = pd.DataFrame([{
        "E": args.E, "S_Agmon": S, "T_exact": T_exact, "T_WKB": T_wkb,
        "rel_error": abs(T_wkb - T_exact)/max(T_exact, 1e-16)
    }])
    df.to_csv(args.out, index=False)
    print(f"Wrote {args.out}")
