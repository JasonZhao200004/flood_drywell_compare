#!/usr/bin/env python3

from pathlib import Path
import sys
import math
import xml.etree.ElementTree as ET

import numpy as np
import meshio


PVD = Path(sys.argv[1]).resolve()

C2_START_H = 480.0
C2_WET_END_H = 725.850601953
C2_END_H = 960.0

PLUME_FRAC = 0.10


def read_pvd(path):

    root = ET.parse(path).getroot()

    records = []

    for ds in root.findall(".//DataSet"):

        t = float(ds.attrib["timestep"])

        f = (
            path.parent /
            ds.attrib["file"]
        ).resolve()

        records.append((t, f))

    records.sort()

    if records[-1][0] > 2000:

        records = [
            (t/3600.0, f)
            for t, f in records
        ]

    return records


records = read_pvd(PVD)


def exact_or_before(target):

    exact = [
        (t, f)
        for t, f in records
        if abs(t-target) < 1e-8
    ]

    if exact:
        return exact[0]

    before = [
        (t, f)
        for t, f in records
        if t < target
    ]

    return before[-1]


tb, fb = exact_or_before(C2_START_H)
tw, fw = exact_or_before(C2_WET_END_H)
te, fe = exact_or_before(C2_END_H)


print("="*78)
print("C2 STATES")
print("="*78)

print(
    f"pre-C2 background = {tb:.9f} h"
)

print(
    f"C2 end-wet        = {tw:.9f} h "
    f"(target {C2_WET_END_H:.9f})"
)

print(
    f"C2 end-cycle      = {te:.9f} h "
    f"(target {C2_END_H:.9f})"
)

print(
    f"C2 drainage       = "
    f"{C2_END_H-C2_WET_END_H:.6f} h "
    f"= "
    f"{(C2_END_H-C2_WET_END_H)/24:.6f} d"
)


def read_state(path):

    mesh = meshio.read(path)

    pts = np.asarray(
        mesh.points,
        float
    )[:, :2]

    tris = []
    o2 = []

    for bi, block in enumerate(mesh.cells):

        if block.type != "triangle":
            continue

        conn = np.asarray(
            block.data,
            int
        )

        tris.append(conn)

        o2.append(
            np.asarray(
                mesh.cell_data[
                    "x^O2_liq"
                ][bi],
                float
            ).reshape(-1)
        )

    tris = np.vstack(tris)
    o2 = np.concatenate(o2)

    centers = pts[tris].mean(axis=1)

    a = pts[tris[:,0]]
    b = pts[tris[:,1]]
    c = pts[tris[:,2]]

    area2d = 0.5*np.abs(
        (b[:,0]-a[:,0])*(c[:,1]-a[:,1])
        -
        (c[:,0]-a[:,0])*(b[:,1]-a[:,1])
    )

    r = centers[:,0]

    vol = (
        2.0
        * math.pi
        * r
        * area2d
    )

    return (
        o2,
        centers[:,0],
        centers[:,1],
        vol
    )


CB, r, z, vol = read_state(fb)
CW, _, _, _ = read_state(fw)
CE, _, _, _ = read_state(fe)


wet_excess = np.maximum(
    CW-CB,
    0.0
)

end_excess = np.maximum(
    CE-CB,
    0.0
)

max_excess = np.max(wet_excess)

threshold = (
    PLUME_FRAC
    * max_excess
)

base_plume = (
    wet_excess >= threshold
)


def report(mask, title):

    wet = np.sum(
        wet_excess[mask]
        * vol[mask]
    )

    end = np.sum(
        end_excess[mask]
        * vol[mask]
    )

    residual = (
        end/wet
        if wet > 0
        else np.nan
    )

    local = (
        end_excess[mask]
        / max_excess
    )

    print()
    print("="*78)
    print(title)
    print("="*78)

    print(
        "cells                  =",
        np.count_nonzero(mask)
    )

    print(
        f"radial range           = "
        f"{np.min(r[mask]):.3f} "
        f"to "
        f"{np.max(r[mask]):.3f} m"
    )

    print(
        f"z range                = "
        f"{np.min(z[mask]):.3f} "
        f"to "
        f"{np.max(z[mask]):.3f} m"
    )

    print(
        f"end-wet integral       = "
        f"{wet:.8e}"
    )

    print(
        f"day-40 residual        = "
        f"{end:.8e}"
    )

    print(
        f"POSITIVE RESIDUAL      = "
        f"{100*residual:.6f}%"
    )

    print(
        f"POSITIVE RECOVERY      = "
        f"{100*(1-residual):.6f}%"
    )

    for frac in [
        0.01,
        0.05,
        0.10
    ]:

        rem = np.mean(
            local > frac
        )

        print(
            f"cells still > "
            f"{100*frac:.0f}% max = "
            f"{100*rem:.6f}%"
        )

    below = np.mean(
        CE[mask] < CB[mask]
    )

    print(
        f"cells below pre-C2 background = "
        f"{100*below:.6f}%"
    )


report(
    base_plume,
    "FULL C2 POSITIVE PLUME AT DAY 40"
)


main_plume = (
    base_plume
    & (z <= 27.0)
)

report(
    main_plume,
    "MAIN DRYWELL C2 PLUME (z <= 27 m) AT DAY 40"
)
