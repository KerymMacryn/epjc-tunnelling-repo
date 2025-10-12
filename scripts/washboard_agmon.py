
import argparse, numpy as np, pandas as pd
import matplotlib.pyplot as plt

def washboard_potential(phi, EJ, eta):
    # U(phi) = -EJ (cos phi + eta * phi)  up to a constant
    return -EJ*(np.cos(phi) + eta*phi)

def find_min_and_barrier(EJ, eta, phi_window=(-6*np.pi, 6*np.pi), N=20000):
    phi = np.linspace(phi_window[0], phi_window[1], N)
    U = washboard_potential(phi, EJ, eta)
    # crude locate a local min and next barrier to the right
    i_min = np.argmin(U + (phi% (2*np.pi))*0)  # keep periodic structure
    # better: look around multiples of 2π
    i_min = np.argmin(U)
    # next barrier: first local max to the right
    for j in range(i_min+1, N-1):
        if U[j-1] < U[j] and U[j] > U[j+1]:
            i_bar = j
            break
    else:
        i_bar = min(i_min+2000, N-1)
    return phi[i_min], phi[i_bar], U[i_min], U[i_bar]

def omega_p(EJ, EC, phi_min, eta):
    # U''(phi) = EJ cos(phi); plasma freq (dimensionless units): sqrt(8*EC*U'')
    return np.sqrt(max(8.0*EC*EJ*np.cos(phi_min), 0.0))

def agmon_action(EJ, EC, E_level, phi1, phi2, n=4000):
    # S = ∫ sqrt((U(phi)-E_level)/(4 EC)) dphi over forbidden region
    import math
    phis = np.linspace(phi1, phi2, n)
    U = -EJ*(np.cos(phis) + eta*phis)  # eta is closed over main
    W = np.maximum(U - E_level, 0.0)
    return float(np.trapz(np.sqrt(W/(4.0*EC)), phis))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--EJ", type=float, required=True)
    ap.add_argument("--EC", type=float, required=True)
    ap.add_argument("--etas", type=float, nargs="+", required=True)
    ap.add_argument("--csv", type=str, required=True)
    ap.add_argument("--fig", type=str, default=None)
    args = ap.parse_args()

    rows = []
    for eta in args.etas:
        phi_min, phi_bar, U_min, U_bar = find_min_and_barrier(args.EJ, eta)
        om = omega_p(args.EJ, args.EC, phi_min, eta)
        # choose a metastable level not too close to the barrier
        E_level = U_min + 0.25*(U_bar - U_min)
        S = agmon_action(args.EJ, args.EC, E_level, min(phi_min, phi_bar), max(phi_min, phi_bar))
        Gamma_est = (om/(2*np.pi))*np.exp(-2*S)
        rows.append({
            "eta": eta, "phi_min": phi_min, "phi_barrier": phi_bar,
            "U_min": U_min, "U_barrier": U_bar, "dU": U_bar-U_min,
            "omega_p": om, "E_level": E_level, "S": S, "Gamma_est": Gamma_est
        })

    df = pd.DataFrame(rows)
    out = args.csv
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"Wrote {out}")

    if args.fig:
        # plot potential & annotate a representative case (first eta)
        eta0 = args.etas[0]
        phi = np.linspace(-3*np.pi, 3*np.pi, 2000)
        U = washboard_potential(phi, args.EJ, eta0)
        fig = plt.figure(figsize=(7.2, 4.2))
        plt.plot(phi, U, label=f"washboard (eta={eta0:.2f})")
        plt.xlabel(r"$\varphi$"); plt.ylabel(r"$U(\varphi)$")
        plt.title("Josephson washboard potential")
        plt.grid(True, alpha=0.3)
        Path(args.fig).parent.mkdir(parents=True, exist_ok=True)
        plt.tight_layout()
        fig.savefig(args.fig, dpi=180)
        plt.close(fig)
        print(f"Wrote {args.fig}")
