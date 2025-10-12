
import argparse
from pathlib import Path
import numpy as np, pandas as pd

def grad_sup_gaussian(a): return abs(a)/0.2*np.exp(-0.5)   # σ=0.2
def grad_sup_cosine(a):   return abs(a)*np.pi               # φ(x)=a cos(πx)
def grad_sup_tophat(a):   return 1e-3*abs(a)                # regularized jump

def make_table(sector, order, Kstar, Cd, disc_err=0.01, domain_err=0.01, n=3):
    rows=[]
    for i in range(n):
        name = order[i % len(order)]
        a = float(0.5 + 0.5*i)  # 0.5, 1.0, 1.5
        gsup = {"gaussian":grad_sup_gaussian,"cosine":grad_sup_cosine,"tophat":grad_sup_tophat}[name](a)
        lb = Kstar - Cd*(gsup**2)
        rows.append({
            "sector": sector,
            "profile_id": i+1,
            "profile_name": name,
            "a_param": a,
            "Kstar_int": Kstar,
            "C_d_int": Cd,
            "grad_phi_sup": gsup,
            "lambda1_lower_bound": lb,
            "disc_error": disc_err,
            "domain_error": domain_err,
            "certificate_margin": lb - (disc_err+domain_err)
        })
    return pd.DataFrame(rows)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--su2_csv", required=True)
    ap.add_argument("--su3_csv", required=True)
    ap.add_argument("--Kstar_su2", type=float, default=0.50)
    ap.add_argument("--Kstar_su3", type=float, default=0.60)
    ap.add_argument("--Cd_su2",    type=float, default=0.75)
    ap.add_argument("--Cd_su3",    type=float, default=0.80)
    ap.add_argument("--disc_err",  type=float, default=0.01)
    ap.add_argument("--domain_err",type=float, default=0.01)
    args = ap.parse_args()

    su2 = make_table("SU(2)", ["gaussian","cosine","tophat"], args.Kstar_su2, args.Cd_su2, args.disc_err, args.domain_err)
    su3 = make_table("SU(3)", ["cosine","gaussian","tophat"], args.Kstar_su3, args.Cd_su3, args.disc_err, args.domain_err)

    Path(args.su2_csv).parent.mkdir(parents=True, exist_ok=True)
    su2.to_csv(args.su2_csv, index=False)
    Path(args.su3_csv).parent.mkdir(parents=True, exist_ok=True)
    su3.to_csv(args.su3_csv, index=False)
    print("Wrote", args.su2_csv, "and", args.su3_csv)
