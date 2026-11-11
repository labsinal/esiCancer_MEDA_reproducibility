"""Folder layout of the repository; every script takes its paths from here."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SIM = ROOT / "data" / "simulations"            # esiCancer per-run outputs (one folder per series)
GERMLINE = SIM / "germline"                    # inherited first-hit runs
TURNOVER = SIM / "phase2_turnover"             # matched-input turnover runs (when available)
EPI = ROOT / "data" / "epidemiology"           # published epidemiological inputs
OUT = ROOT / "output"                          # tables and text results
FIG = ROOT / "figures"                         # figures
for _d in (OUT, FIG):
    _d.mkdir(exist_ok=True)
