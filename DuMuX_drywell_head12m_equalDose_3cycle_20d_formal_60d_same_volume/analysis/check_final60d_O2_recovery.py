#!/usr/bin/env python3

from pathlib import Path
import sys
import math
import bisect
import xml.etree.ElementTree as ET

import numpy as np
import meshio


PVD = Path(sys.argv[1]).resolve()

# ------------------------------------------------------------
# FINAL 60-DAY FORMAL SCHEDULE
# ------------------------------------------------------------

CYCLES = [
    {
        "name": "C1",
        "start": 0.0,
        "wet_end": 176.490541000,
        "end": 480.0,
    },
    {
        "name": "C2",
        "start": 480.0,
        "wet_end": 725.850601953,
        "end": 960.0,
    },
    {
        "name": "C3",
        "start": 960.0,
        "wet_end": 1210.359144226,
        "end": 1440.0,
    },
]

PLUME_FRAC = 0.10
MAIN_Z_MAX = 27.0


# ============================================================
# READ PVD
# ============================================================

root = ET.parse(PVD).getroot()

records = []

for ds in root.findall(".//DataSet"):

    t = float(ds.attrib["timestep"])
    f = (PVD.parent / ds.attrib["file"]).resolve()

    records.append((t, f))

records.sort()

if records[-1][0] > 2000:
    records = [(t/3600.0, f) for t, f in records]

times = np.array([t for t, _ in records], dtype=float)


print("="*80)
print("FINAL 60-DAY RUN COMPLETENESS")
print("="*80)

missing = [f for _, f in records if not f.exists()]

print("PVD             :", PVD)
print("outputs         :", len(records))
print("first time [h]  :", times[0])
print("last time [h]   :", times[-1])
print("missing VTUs    :", len(missing))

if abs(times[-1] - 1440.0) < 0.01 and not missing:
    print("STATUS          : COMPLETE")
else:
    print("STATUS          : CHECK REQUIRED")


# ============================================================
# READ TRIANGLE GEOMETRY + O2
# ============================================================

cache = {}


def read_o2(path):

    path = Path(path)

    if path in cache:
        return cache[path]

    mesh = meshio.read(path)

    vals = []

    for bi, block in enumerate(mesh.cells):

        if block.type != "triangle":
            continue

        vals.append(
            np.asarray(
                mesh.cell_data["x^O2_liq"][bi],
                dtype=float
            ).reshape(-1)
        )

    vals = np.concatenate(vals)

    cache[path] = vals

    return vals


# Geometry from first file
mesh0 = meshio.read(records[0][1])

pts = np.asarray(mesh0.points, dtype=float)[:, :2]

tris = []

for block in mesh0.cells:
    if block.type == "triangle":
        tris.append(np.asarray(block.data, dtype=int))

tris = np.vstack(tris)

centers = pts[tris].mean(axis=1)

r = centers[:, 0]
z = centers[:, 1]

a = pts[tris[:, 0]]
b = pts[tris[:, 1]]
c = pts[tris[:, 2]]

area2d = 0.5*np.abs(
    (b[:,0]-a[:,0])*(c[:,1]-a[:,1])
    -
    (c[:,0]-a[:,0])*(b[:,1]-a[:,1])
)

# axisymmetric volume weighting
vol = 2.0*math.pi*r*area2d


# ============================================================
# INTERPOLATED STATE AT ARBITRARY TIME
# ============================================================

def state_at(target):

    # exact
    j = int(np.argmin(np.abs(times-target)))

    if abs(times[j]-target) < 1e-8:
        return read_o2(records[j][1]), times[j], times[j]

    # bracket
    k = bisect.bisect_left(times.tolist(), target)

    if k == 0 or k == len(records):
        raise RuntimeError(
            f"Cannot bracket target time {target}"
        )

    t0, f0 = records[k-1]
    t1, f1 = records[k]

    y0 = read_o2(f0)
    y1 = read_o2(f1)

    w = (target-t0)/(t1-t0)

    y = y0 + w*(y1-y0)

    return y, t0, t1


# ============================================================
# REPORT
# ============================================================

def plume_metrics(base, wet, current, mask):

    wet_excess = np.maximum(wet-base, 0.0)
    cur_excess = np.maximum(current-base, 0.0)

    wet_int = np.sum(
        wet_excess[mask]*vol[mask]
    )

    cur_int = np.sum(
        cur_excess[mask]*vol[mask]
    )

    residual = (
        cur_int/wet_int
        if wet_int > 0
        else np.nan
    )

    return wet_int, cur_int, residual


all_results = []


for C in CYCLES:

    name = C["name"]
    ts = C["start"]
    tw = C["wet_end"]
    te = C["end"]

    base, bs0, bs1 = state_at(ts)
    wet,  bw0, bw1 = state_at(tw)
    end,  be0, be1 = state_at(te)

    wet_excess = np.maximum(wet-base, 0.0)

    max_excess = np.max(wet_excess)

    plume_threshold = (
        PLUME_FRAC*max_excess
    )

    full_plume = (
        wet_excess >= plume_threshold
    )

    main_plume = (
        full_plume
        & (z <= MAIN_Z_MAX)
    )

    drainage_h = te-tw

    print()
    print("="*80)
    print(f"{name} O2 RECOVERY")
    print("="*80)

    print(
        f"cycle start       = {ts:.9f} h"
    )

    print(
        f"wet end           = {tw:.9f} h"
    )

    print(
        f"cycle end         = {te:.9f} h"
    )

    print(
        f"wetting duration  = {tw-ts:.6f} h"
        f" = {(tw-ts)/24:.6f} d"
    )

    print(
        f"drainage duration = {drainage_h:.6f} h"
        f" = {drainage_h/24:.6f} d"
    )

    print(
        f"wet VTU bracket   = "
        f"{bw0:.9f} -- {bw1:.9f} h"
    )

    print(
        f"max positive O2 enrichment = "
        f"{max_excess:.8e}"
    )

    print(
        f"10% plume threshold        = "
        f"{plume_threshold:.8e}"
    )


    for title, mask in [
        ("FULL POSITIVE PLUME", full_plume),
        ("MAIN DRYWELL PLUME (z <= 27 m)", main_plume),
    ]:

        wet_int, end_int, residual = plume_metrics(
            base,
            wet,
            end,
            mask
        )

        end_excess = np.maximum(
            end-base,
            0.0
        )

        local = (
            end_excess[mask]
            / max_excess
        )

        print()
        print("-"*80)
        print(title)
        print("-"*80)

        print(
            f"cells                 = "
            f"{np.count_nonzero(mask)}"
        )

        print(
            f"radial range          = "
            f"{np.min(r[mask]):.3f}"
            f" to "
            f"{np.max(r[mask]):.3f} m"
        )

        print(
            f"z range               = "
            f"{np.min(z[mask]):.3f}"
            f" to "
            f"{np.max(z[mask]):.3f} m"
        )

        print(
            f"end-wet integral      = "
            f"{wet_int:.8e}"
        )

        print(
            f"end-cycle integral    = "
            f"{end_int:.8e}"
        )

        print(
            f"POSITIVE RESIDUAL     = "
            f"{100*residual:.6f}%"
        )

        print(
            f"POSITIVE RECOVERY     = "
            f"{100*(1-residual):.6f}%"
        )

        for frac in [0.01, 0.05, 0.10]:

            pct = 100*np.mean(
                local > frac
            )

            print(
                f"cells still > "
                f"{100*frac:.0f}% max = "
                f"{pct:.6f}%"
            )

        below = 100*np.mean(
            end[mask] < base[mask]
        )

        print(
            f"cells below cycle-start background = "
            f"{below:.6f}%"
        )


    # --------------------------------------------------------
    # DAILY DRAINAGE TRAJECTORY
    # Use MAIN drywell plume
    # --------------------------------------------------------

    wet_int = np.sum(
        wet_excess[main_plume]
        * vol[main_plume]
    )

    print()
    print("-"*80)
    print(
        f"{name} MAIN-PLUME DRAINAGE TRAJECTORY"
    )
    print("-"*80)

    traj = []

    # t = wet-end
    check_times = [tw]

    # every 24 h after wetting
    t = tw + 24.0

    while t < te-1e-8:
        check_times.append(t)
        t += 24.0

    if abs(check_times[-1]-te) > 1e-8:
        check_times.append(te)

    for tt in check_times:

        state, _, _ = state_at(tt)

        excess = np.maximum(
            state-base,
            0.0
        )

        integ = np.sum(
            excess[main_plume]
            * vol[main_plume]
        )

        residual = (
            integ/wet_int
            if wet_int > 0
            else np.nan
        )

        drain_day = (tt-tw)/24.0

        traj.append(
            (drain_day, residual)
        )

        print(
            f"drain + {drain_day:7.3f} d"
            f" : residual = "
            f"{100*residual:10.6f}%"
        )


    final_res = traj[-1][1]

    all_results.append(
        (
            name,
            (tw-ts)/24.0,
            drainage_h/24.0,
            final_res
        )
    )


# ============================================================
# FINAL SUMMARY
# ============================================================

print()
print("="*80)
print("FINAL 60-DAY O2 RECOVERY SUMMARY")
print("="*80)

print(
    f"{'Cycle':<8}"
    f"{'Wet [d]':>12}"
    f"{'Drain [d]':>14}"
    f"{'Main plume residual':>24}"
    f"{'Recovery':>14}"
)

for name, wetd, draind, residual in all_results:

    print(
        f"{name:<8}"
        f"{wetd:12.6f}"
        f"{draind:14.6f}"
        f"{100*residual:23.6f}%"
        f"{100*(1-residual):13.6f}%"
    )

print()
print("Interpretation:")
print(
    "Residual here means positive O2 enrichment relative to the "
    "cycle-start O2 field, evaluated inside the plume footprint "
    "defined at end-wetting."
)
print(
    "It should be interpreted as cycle-specific positive-plume "
    "persistence, not automatically as whole-field return to the "
    "original 0-h background."
)
