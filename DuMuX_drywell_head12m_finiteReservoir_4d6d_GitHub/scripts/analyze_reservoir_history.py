#!/usr/bin/env python3
import csv
from pathlib import Path

here = Path(__file__).resolve().parents[1]
f = here / "results" / "diagnostic_240h" / "reservoir_production_240h_reservoir_history.csv"

with f.open(newline="") as fh:
    rows = list(csv.DictReader(fh))

for r in rows:
    for k in r:
        r[k] = float(r[k])

post = [r for r in rows if r["time_h"] > 96.0]

for target in [96, 96.1, 97, 98, 100, 102, 108, 120, 144, 168, 192, 216, 240]:
    r = min(rows, key=lambda x: abs(x["time_h"] - target))
    print(
        f'{target:6.1f} h : '
        f't={r["time_h"]:10.6f}  '
        f'H={r["reservoir_head_cm"]:10.3f} cm  '
        f'V={r["reservoir_storage_m3"]:10.5f} m3  '
        f'Q={r["reservoir_net_flux_m3_s"]:.6e} m3/s'
    )

positive = [r for r in post if r["reservoir_net_flux_m3_s"] > 0.0]
empty = [r for r in post if r["reservoir_head_cm"] <= 1e-9]

print()
print("positive soil->well flux rows =", len(positive))
if positive:
    print("first positive soil->well flux row =", positive[0])
print("reservoir empty by 240 h =", bool(empty))
print("final state =", rows[-1])

max_res = 0.0
sum_res = 0.0
n = 0
prev = None
for r in rows:
    if prev is not None and r["time_h"] > 96.0:
        dt = (r["time_h"] - prev["time_h"])*3600.0
        dV = r["reservoir_storage_m3"] - prev["reservoir_storage_m3"]
        pred = r["reservoir_net_flux_m3_s"]*dt
        res = dV - pred
        max_res = max(max_res, abs(res))
        sum_res += abs(res)
        n += 1
    prev = r

print()
print("mass-balance residual:")
print("  max abs  =", max_res, "m3")
print("  mean abs =", sum_res/n if n else 0.0, "m3")
