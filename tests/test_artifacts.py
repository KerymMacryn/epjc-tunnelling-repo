
import pandas as pd, numpy as np
from pathlib import Path

def test_ab_csv_schema():
    p = Path("artifacts/data/ab_shape_resonance_poles.csv")
    assert p.exists(), "AB CSV missing"
    df = pd.read_csv(p)
    for col in ["m","alpha","Re(kR)","Im(kR)","Gamma"]:
        assert col in df.columns, f"Missing column: {col}"
    assert (df["Im(kR)"]<0).all()

def test_washboard_csv_schema():
    p = Path("artifacts/data/washboard_escape_rates.csv")
    assert p.exists(), "washboard CSV missing"
    df = pd.read_csv(p)
    for col in ["eta","phi_min","phi_barrier","U_min","U_barrier","dU","omega_p","E_level","S","Gamma_est"]:
        assert col in df.columns, f"Missing column: {col}"
