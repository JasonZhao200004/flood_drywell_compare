#!/usr/bin/env python3
"""Extract the closest available profile to 2880 h from the CSV."""

import csv
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
source = ROOT / "x3m_all_profiles_g_cm3.csv"
target = ROOT / "oxygen_initial_2880.dat"

with source.open(newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    columns = reader.fieldnames or []
    candidates = []
    for name in columns[1:]:
        match = re.search(r"O2_([0-9.]+)h_g_cm3", name)
        if match:
            candidates.append((abs(float(match.group(1)) - 2880.0), name))
    if not candidates:
        raise RuntimeError("No time-profile columns found")
    _, selected = min(candidates)
    rows = [(float(row[columns[0]]), float(row[selected])) for row in reader]

with target.open("w", encoding="ascii") as f:
    for elevation, concentration in rows:
        f.write(f"{elevation:.10e} {concentration:.10e}\n")

print(f"Selected {selected}")
print(f"Wrote {target} with {len(rows)} rows")
