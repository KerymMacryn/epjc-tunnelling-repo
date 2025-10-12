
import hashlib, sys
from pathlib import Path

FILES = [
    Path("artifacts/data/washboard_escape_rates.csv"),
    Path("artifacts/data/ab_shape_resonance_poles.csv"),
]

out = Path("artifacts/CHECKSUMS.txt")
lines = []
for p in FILES:
    if p.exists():
        h = hashlib.sha256()
        with p.open("rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                h.update(chunk)
        lines.append(f"{h.hexdigest()}  {p.as_posix()}")
    else:
        lines.append(f"<pending>  {p.as_posix()}")

out.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"Wrote {out}")
