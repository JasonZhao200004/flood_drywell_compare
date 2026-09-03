#!/usr/bin/env python3

from pathlib import Path
import sys
import math
import csv
import xml.etree.ElementTree as ET

import numpy as np
import meshio


# ======================================================================
# LOCKED DESIGN
# ======================================================================

PVD = Path(sys.argv[1]).resolve()

TARGET_M3 = 3769.844

W1_H = 176.490541
W2_H = 354.547396832

C3_START_H = 480.0
EXPECTED_END_H = 960.0

HEAD_CM = 1200.0

P_ATM = 1.0e5
RHO_W = 1000.0
G = 9.81

DRYWELL_RADIUS = 0.61
DRYWELL_BOTTOM_Z = 14.0
DRYWELL_WALL_TOP_Z = 27.0

# Same reconstruction resolution used for C2.
DT_H = 0.1


# ======================================================================
# PVD
# ======================================================================

def read_pvd(path):

    root = ET.parse(path).getroot()

    records = []

    for ds in root.findall(".//DataSet"):

        t = float(ds.attrib["timestep"])

        f = (
            path.parent
            / ds.attrib["file"]
        ).resolve()

        records.append((t, f))

    records.sort(key=lambda x: x[0])

    if not records:
        raise RuntimeError(
            f"No DataSet records in {path}"
        )

    # DuMuX PVD may store seconds.
    if records[-1][0] > 2000.0:

        records = [
            (t/3600.0, f)
            for t, f in records
        ]

    return records


# ======================================================================
# TRIANGLE / FIELD HELPERS
# ======================================================================

def triangle_info(mesh):

    blocks = []
    conns = []

    for bi, block in enumerate(mesh.cells):

        if block.type.startswith("triangle"):

            arr = np.asarray(
                block.data,
                dtype=int
            )

            blocks.append(
                (bi, len(arr))
            )

            conns.append(arr)

    if not conns:
        raise RuntimeError(
            "No triangle cells"
        )

    return (
        np.vstack(conns),
        blocks
    )


def cell_field(mesh, name, blocks):

    if name not in mesh.cell_data:

        raise RuntimeError(
            f"Missing field {name}; "
            f"available="
            f"{sorted(mesh.cell_data.keys())}"
        )

    out = []

    for bi, n in blocks:

        arr = np.asarray(
            mesh.cell_data[name][bi],
            dtype=float
        ).reshape(-1)

        if len(arr) != n:

            raise RuntimeError(
                f"{name}: cell block mismatch"
            )

        out.append(arr)

    return np.concatenate(out)


# ======================================================================
# LAYERED PERMEABILITY
# ======================================================================

def permeability_from_z(z):

    z = np.asarray(
        z,
        dtype=float
    )

    return np.select(
        [
            z >= 22.0,
            z >= 18.5,
            z >= 12.5,
            z >= 8.5,
            z >= 4.0,
        ],
        [
            1.127e-11,
            3.709e-12,
            1.925e-12,
            1.161e-13,
            2.633e-13,
        ],
        default=1.812e-12,
    )


# ======================================================================
# IDENTIFY DRYWELL ROBIN BOUNDARY
# ======================================================================

def build_faces(vtu):

    mesh = meshio.read(vtu)

    pts = np.asarray(
        mesh.points,
        dtype=float
    )[:, :2]

    tris, _ = triangle_info(mesh)

    centers = (
        pts[tris].mean(axis=1)
    )

    owners = {}

    for ci, tri in enumerate(tris):

        for a, b in (
            (int(tri[0]), int(tri[1])),
            (int(tri[1]), int(tri[2])),
            (int(tri[2]), int(tri[0])),
        ):

            key = tuple(
                sorted((a, b))
            )

            owners.setdefault(
                key,
                []
            ).append(ci)

    exterior = [
        (edge, cells[0])
        for edge, cells
        in owners.items()
        if len(cells) == 1
    ]

    tol = 2.0e-4

    faces = []

    for (a, b), owner in exterior:

        p0 = pts[a]
        p1 = pts[b]

        mid = 0.5*(p0+p1)

        bottom = (
            abs(
                p0[1]
                - DRYWELL_BOTTOM_Z
            ) <= tol
            and
            abs(
                p1[1]
                - DRYWELL_BOTTOM_Z
            ) <= tol
            and
            min(
                p0[0],
                p1[0]
            ) >= -tol
            and
            max(
                p0[0],
                p1[0]
            )
            <= DRYWELL_RADIUS + tol
        )

        wall = (
            abs(
                p0[0]
                - DRYWELL_RADIUS
            ) <= tol
            and
            abs(
                p1[0]
                - DRYWELL_RADIUS
            ) <= tol
            and
            min(
                p0[1],
                p1[1]
            )
            >= DRYWELL_BOTTOM_Z - tol
            and
            max(
                p0[1],
                p1[1]
            )
            <= DRYWELL_WALL_TOP_Z + tol
        )

        if not (
            bottom or wall
        ):
            continue

        distance = float(
            np.linalg.norm(
                centers[owner] - mid
            )
        )

        if distance <= 0:
            continue

        if bottom:

            r0 = float(
                min(
                    p0[0],
                    p1[0]
                )
            )

            r1 = float(
                max(
                    p0[0],
                    p1[0]
                )
            )

            area = (
                math.pi
                * (
                    r1**2
                    - r0**2
                )
            )

            nz = 1.0
            kind = "bottom"

        else:

            R = 0.5*float(
                p0[0] + p1[0]
            )

            h = abs(
                float(
                    p1[1] - p0[1]
                )
            )

            area = (
                2.0
                * math.pi
                * R
                * h
            )

            nz = 0.0
            kind = "wall"

        faces.append({
            "owner":
                owner,

            "owner_z":
                float(
                    centers[
                        owner,
                        1
                    ]
                ),

            "z_mid":
                float(mid[1]),

            "distance":
                distance,

            "area":
                area,

            "nz":
                nz,

            "kind":
                kind,
        })

    if not faces:

        raise RuntimeError(
            "No drywell boundary faces found"
        )

    unique_owners = np.asarray(
        sorted(
            set(
                f["owner"]
                for f in faces
            )
        ),
        dtype=int
    )

    lookup = {
        owner: i
        for i, owner
        in enumerate(
            unique_owners
        )
    }

    for f in faces:

        f["state_i"] = (
            lookup[
                f["owner"]
            ]
        )

    return (
        faces,
        unique_owners
    )


# ======================================================================
# READ BOUNDARY-ADJACENT STATE
# ======================================================================

def read_state(vtu, owners):

    mesh = meshio.read(vtu)

    _, blocks = triangle_info(mesh)

    p = cell_field(
        mesh,
        "p_liq",
        blocks
    )[owners]

    rho = cell_field(
        mesh,
        "rho_liq",
        blocks
    )[owners]

    mob = cell_field(
        mesh,
        "mob_liq",
        blocks
    )[owners]

    return p, rho, mob


# ======================================================================
# RUN COMPLETENESS
# ======================================================================

records = read_pvd(PVD)

missing = [
    f
    for _, f in records
    if not f.exists()
]

end_offset_s = (
    records[-1][0]
    - EXPECTED_END_H
) * 3600.0


print("="*78)
print("RUN COMPLETENESS")
print("="*78)

print(
    "PVD              :",
    PVD
)

print(
    "outputs          :",
    len(records)
)

print(
    "first time [h]   :",
    records[0][0]
)

print(
    "last time [h]    :",
    records[-1][0]
)

print(
    "expected end [h] :",
    EXPECTED_END_H
)

print(
    "end offset [s]   :",
    f"{end_offset_s:+.6f}"
)

print(
    "missing VTUs     :",
    len(missing)
)

if missing:

    raise RuntimeError(
        f"Missing VTU: {missing[0]}"
    )

if (
    records[-1][0]
    < EXPECTED_END_H - 0.01
):

    raise RuntimeError(
        "C3 calibration ended "
        "materially before expected end"
    )


# ======================================================================
# BRACKET EXACT C3 START
# ======================================================================

before = [
    (t, f)
    for t, f in records
    if t <= C3_START_H
]

after = [
    (t, f)
    for t, f in records
    if t >= C3_START_H
]

if (
    not before
    or not after
):

    raise RuntimeError(
        "Cannot bracket exact C3 start"
    )

pre_record = before[-1]
post_record = after[0]


print()
print("="*78)
print("NEW CYCLE-2 START BRACKET")
print("="*78)

print(
    f"scheduled C3 start = "
    f"{C3_START_H:.9f} h"
)

print(
    f"VTU before         = "
    f"{pre_record[0]:.9f} h"
)

print(
    f"VTU after          = "
    f"{post_record[0]:.9f} h"
)

print(
    f"before offset      = "
    f"{(pre_record[0]-C3_START_H)*3600.0:+.6f} s"
)

print(
    f"after offset       = "
    f"{(post_record[0]-C3_START_H)*3600.0:+.6f} s"
)


# One record immediately before the event,
# plus all records after the event.
state_records = [
    pre_record
]

for t, f in records:

    if (
        t >=
        post_record[0] - 1e-12
    ):

        state_records.append(
            (t, f)
        )


# Deduplicate.
tmp = {}

for t, f in state_records:
    tmp[t] = f

state_records = sorted(
    tmp.items(),
    key=lambda x: x[0]
)


# ======================================================================
# GEOMETRY
# ======================================================================

faces, unique_owners = (
    build_faces(
        state_records[0][1]
    )
)


print()
print("="*78)
print("DRYWELL BOUNDARY")
print("="*78)

print(
    "faces        :",
    len(faces)
)

print(
    "bottom faces :",
    sum(
        f["kind"] == "bottom"
        for f in faces
    )
)

print(
    "wall faces   :",
    sum(
        f["kind"] == "wall"
        for f in faces
    )
)

print(
    "total area   :",
    sum(
        f["area"]
        for f in faces
    ),
    "m2"
)


# Compact face arrays.

state_i = np.asarray(
    [
        f["state_i"]
        for f in faces
    ],
    dtype=int
)

z_mid = np.asarray(
    [
        f["z_mid"]
        for f in faces
    ],
    dtype=float
)

distance = np.asarray(
    [
        f["distance"]
        for f in faces
    ],
    dtype=float
)

area = np.asarray(
    [
        f["area"]
        for f in faces
    ],
    dtype=float
)

nz = np.asarray(
    [
        f["nz"]
        for f in faces
    ],
    dtype=float
)

K = permeability_from_z(
    np.asarray(
        [
            f["owner_z"]
            for f in faces
        ],
        dtype=float
    )
)


# ======================================================================
# READ C3 BOUNDARY STATES
# ======================================================================

times = []
P = []
RHO = []
MOB = []


print()
print("="*78)
print("READING NEW CYCLE-2 BOUNDARY STATES")
print("="*78)


for j, (t, vtu) in enumerate(
    state_records
):

    p, rho, mob = (
        read_state(
            vtu,
            unique_owners
        )
    )

    times.append(t)
    P.append(p)
    RHO.append(rho)
    MOB.append(mob)

    if (
        j % 50 == 0
        or
        j == len(state_records)-1
    ):

        print(
            f"{j+1}/"
            f"{len(state_records)} "
            f"t={t:.6f} h"
        )


times = np.asarray(
    times,
    dtype=float
)

P = np.asarray(
    P,
    dtype=float
)

RHO = np.asarray(
    RHO,
    dtype=float
)

MOB = np.asarray(
    MOB,
    dtype=float
)


# ======================================================================
# INTERPOLATE BOUNDARY-CELL STATES
# ======================================================================

def interp_state(t):

    j = int(
        np.searchsorted(
            times,
            t
        )
    )

    if j == 0:

        return (
            P[0],
            RHO[0],
            MOB[0]
        )

    if j >= len(times):

        return (
            P[-1],
            RHO[-1],
            MOB[-1]
        )

    if abs(
        times[j] - t
    ) < 1e-12:

        return (
            P[j],
            RHO[j],
            MOB[j]
        )

    i = j - 1

    f = (
        (t-times[i])
        /
        (times[j]-times[i])
    )

    return (
        P[i]
        + f*(P[j]-P[i]),

        RHO[i]
        + f*(RHO[j]-RHO[i]),

        MOB[i]
        + f*(MOB[j]-MOB[i])
    )


# ======================================================================
# CURRENT DRYWELL ROBIN INFLOW
# ======================================================================

def gross_rate(t):

    p, rho, mob = (
        interp_state(t)
    )

    p_face = p[state_i]
    rho_face = rho[state_i]
    mob_face = mob[state_i]

    # C3 is an operational wet/dry switch.
    # No new startup ramp is imposed.
    local_head_cm = (
        HEAD_CM
        -
        (
            z_mid
            - DRYWELL_BOTTOM_Z
        )
        * 100.0
    )

    wet = (
        local_head_cm > 0.0
    )

    q_out = np.zeros(
        len(faces),
        dtype=float
    )

    if np.any(wet):

        p_boundary = (
            P_ATM
            +
            RHO_W
            * G
            * (
                local_head_cm[wet]
                / 100.0
            )
        )

        q_out[wet] = (
            K[wet]
            *
            mob_face[wet]
            *
            (
                (
                    p_face[wet]
                    -
                    p_boundary
                )
                /
                distance[wet]

                -

                rho_face[wet]
                * G
                * nz[wet]
            )
        )

    # DuMuX convention:
    # q_out < 0 = inward into soil.
    gross_in = float(
        np.sum(
            np.maximum(
                -q_out,
                0.0
            )
            * area
        )
    )

    return gross_in


# ======================================================================
# 0.1-H RECONSTRUCTION
# ======================================================================

t_end = float(
    times[-1]
)

grid = np.arange(
    C3_START_H,
    t_end + DT_H*0.25,
    DT_H
)

grid = grid[
    grid <= t_end + 1e-9
]

if (
    grid[-1]
    < t_end - 1e-8
):

    grid = np.append(
        grid,
        t_end
    )


q = np.empty(
    len(grid),
    dtype=float
)


print()
print("="*78)
print("RECONSTRUCTING NEW CYCLE-2 ROBIN INFLOW")
print("="*78)


for i, t in enumerate(grid):

    q[i] = gross_rate(t)

    if i % 500 == 0:

        print(
            f"{i+1}/"
            f"{len(grid)} "
            f"t={t:.3f} h, "
            f"Q="
            f"{q[i]*3600.0:.6f} "
            f"m3/h"
        )


# cumulative gross injection

cum = np.zeros(
    len(grid),
    dtype=float
)

dt_s = (
    np.diff(grid)
    * 3600.0
)

cum[1:] = np.cumsum(
    0.5
    * (
        q[:-1]
        + q[1:]
    )
    * dt_s
)


# ======================================================================
# SAVE C3 CURVE
# ======================================================================

csv_path = (
    PVD.parent
    /
    "cycle2_after20dC1_equalDose_volume.csv"
)

with csv_path.open(
    "w",
    newline=""
) as f:

    w = csv.writer(f)

    w.writerow([
        "absolute_time_h",
        "cycle3_elapsed_h",
        "gross_injection_rate_m3_h",
        "cycle3_cumulative_injection_m3",
    ])

    for t, qq, vv in zip(
        grid,
        q,
        cum
    ):

        w.writerow([
            f"{t:.9f}",
            f"{t-C3_START_H:.9f}",
            f"{qq*3600.0:.9f}",
            f"{vv:.9f}",
        ])


# ======================================================================
# TARGET CROSSING
# ======================================================================

idx = np.searchsorted(
    cum,
    TARGET_M3
)


print()
print("="*78)
print("NEW CYCLE 2 EQUAL-DOSE RESULT")
print("="*78)

print(
    f"Target                  = "
    f"{TARGET_M3:.6f} m3"
)

print(
    f"C3 start                = "
    f"{C3_START_H:.9f} h"
)

print(
    f"Volume at calibration end = "
    f"{cum[-1]:.6f} m3"
)


if idx >= len(cum):

    print()
    print("TARGET NOT REACHED")

    print(
        f"Shortfall               = "
        f"{TARGET_M3-cum[-1]:.6f} m3"
    )

    print(
        f"C3 wetting simulated    = "
        f"{grid[-1]-C3_START_H:.6f} h"
    )

    raise SystemExit(2)


i0 = idx - 1
i1 = idx

t0 = grid[i0]
t1 = grid[i1]

v0 = cum[i0]
v1 = cum[i1]

fraction = (
    TARGET_M3 - v0
) / (
    v1 - v0
)

cross_h = (
    t0
    +
    fraction
    * (t1-t0)
)

W3_H = (
    cross_h
    - C3_START_H
)

W3_D = (
    W3_H / 24.0
)


print()
print(
    f"Crossing bracket         = "
    f"{t0:.9f} to "
    f"{t1:.9f} h"
)

print(
    f"Volume bracket           = "
    f"{v0:.6f} to "
    f"{v1:.6f} m3"
)

print()
print(
    f"ABSOLUTE C3 END          = "
    f"{cross_h:.9f} h"
)

print(
    f"W3                       = "
    f"{W3_H:.9f} h"
)

print(
    f"W3                       = "
    f"{W3_D:.9f} days"
)

print()
print(
    f"W1                       = "
    f"{W1_H:.9f} h"
)

print(
    f"W2                       = "
    f"{W2_H:.9f} h"
)

print(
    f"W3 - W2                  = "
    f"{W3_H-W2_H:+.9f} h"
)

print(
    f"W3 / W2                  = "
    f"{W3_H/W2_H:.9f}"
)

print(
    f"W3 - W1                  = "
    f"{W3_H-W1_H:+.9f} h"
)

print(
    f"W3 / W1                  = "
    f"{W3_H/W1_H:.9f}"
)


# ======================================================================
# EARLY C3 TRANSIENT — SAME DIAGNOSTIC AS C2
# ======================================================================

print()
print("="*78)
print("EARLY NEW CYCLE-2 INFLOW CHECK")
print("="*78)

for target_h in [
    0,
    0.1,
    0.2,
    0.5,
    1,
    2,
    6,
    12,
    24,
    48
]:

    i = int(
        np.argmin(
            np.abs(
                (
                    grid
                    - C3_START_H
                )
                - target_h
            )
        )
    )

    print(
        f"C3 + "
        f"{grid[i]-C3_START_H:7.3f} h : "
        f"Q = "
        f"{q[i]*3600.0:12.6f} m3/h   "
        f"V = "
        f"{cum[i]:12.6f} m3"
    )


# ======================================================================
# SUMMARY
# ======================================================================

summary = (
    PVD.parent
    /
    "cycle2_after20dC1_equalDose_crossing.txt"
)

summary.write_text(
    "\n".join([
        "Cycle 3 equal-dose calibration",
        "==============================",
        "",
        f"target_m3 = {TARGET_M3:.9f}",
        f"cycle3_start_h = {C3_START_H:.9f}",
        "",
        f"crossing_bracket_h = {t0:.9f}, {t1:.9f}",
        f"volume_bracket_m3 = {v0:.9f}, {v1:.9f}",
        "",
        f"cycle3_end_h = {cross_h:.9f}",
        f"W3_h = {W3_H:.9f}",
        f"W3_days = {W3_D:.9f}",
        "",
        f"W1_h = {W1_H:.9f}",
        f"W2_h = {W2_H:.9f}",
        f"W3_minus_W2_h = {W3_H-W2_H:.9f}",
        f"W3_over_W2 = {W3_H/W2_H:.9f}",
        f"W3_minus_W1_h = {W3_H-W1_H:.9f}",
        f"W3_over_W1 = {W3_H/W1_H:.9f}",
        "",
    ])
)


print()
print(
    "CSV     :",
    csv_path
)

print(
    "Summary :",
    summary
)
