#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

PKG="/Users/zhaozhe/Documents/GitHub/floof_drywell_compare/DuMuX_floodmar_equalVolume_3cycle_20d_60d_GitHub"

echo "======================================================================"
echo "FINAL FLOOD MAR 60-DAY COMPLETION"
echo "======================================================================"
echo
echo "This WILL run DuMuX for ONLY the two missing segments:"
echo "  Part 1: 512.061111 h -> 955 h"
echo "  Part 2: 1030 h       -> 1440 h"
echo
echo "It will NOT rerun 0 -> 512 h."
echo "It will NOT delete existing verified results."
echo

# ----------------------------------------------------------------------
# Locate DuMuX project
# ----------------------------------------------------------------------

ROOT=""

for d in \
    "$HOME/dumux-work/dumux/dumux-floodmar" \
    "$HOME/Documents/GitHub/dumux-floodmar" \
    "$HOME/dumux-floodmar"
do
    if [[ -d "$d/build-cmake" ]]; then
        ROOT="$d"
        break
    fi
done

if [[ -z "$ROOT" ]]; then
    echo "ERROR: could not find dumux-floodmar build."
    echo
    echo "Checked:"
    echo "  $HOME/dumux-work/dumux/dumux-floodmar"
    echo "  $HOME/Documents/GitHub/dumux-floodmar"
    echo "  $HOME/dumux-floodmar"
    exit 1
fi

SRC="$ROOT/test/porousmediumflow/2p2c/floodmar_pond5m_head5cm_equalDose_3cycle_20d_oldformal"
BUILD="$ROOT/build-cmake"
RUN="$BUILD/test/porousmediumflow/2p2c/floodmar_pond5m_head5cm_equalDose_3cycle_20d_oldformal"

TARGET_NAME="floodmar_pond5m_head5cm_equalDose_3cycle_20d_restartwriter"
EXE="$RUN/$TARGET_NAME"

PARAMS="$PKG/model/params.input"
MESH="$PKG/model/floodmar_pond5m_depth10cm_nonames.msh"
IC="$PKG/model/oxygenIC_trueStagnantO2_2880h_2d_xO2.dat"

R512="$PKG/restart/floodmar_pond5m_head5cm_equalVolume_3cycle_60d_formal_v1_restart_from_00512.dat"
R1030="$PKG/restart/floodmar_C3continuous_restart_1030h.dat"

EXISTING="$PKG/existing_verified_6h/floodmar_existing_verified_6h.pvd"

NAME1="floodmar_FINAL_C2_512_to955_2h"
NAME2="floodmar_FINAL_C3_1030_to1440_2h"

RAW="$PKG/raw_mac_runs"
RAW1="$RAW/part1_512_to955_2h"
RAW2="$RAW/part2_1030_to1440_2h"

FINAL="$PKG/final_60d_6h"

mkdir -p "$RAW1" "$RAW2" "$FINAL"

# ----------------------------------------------------------------------
# Validate package
# ----------------------------------------------------------------------

echo "DuMuX root:"
echo "  $ROOT"
echo

for p in \
    "$PKG/model/source" \
    "$PARAMS" \
    "$MESH" \
    "$IC" \
    "$R512" \
    "$R1030" \
    "$EXISTING"
do
    if [[ ! -e "$p" ]]; then
        echo "ERROR: missing required package item:"
        echo "  $p"
        exit 1
    fi
done

NEXIST=$(find "$PKG/existing_verified_6h" -maxdepth 1 -name "*.vtu" | wc -l | tr -d ' ')

echo "Existing verified 6-h VTUs: $NEXIST"

if [[ "$NEXIST" != "109" ]]; then
    echo "ERROR: expected 109 existing verified VTUs."
    exit 1
fi

echo "Package validation: PASS"
echo

# ----------------------------------------------------------------------
# Sync packaged source into Mac DuMuX tree
# ----------------------------------------------------------------------

echo "======================================================================"
echo "SYNC FINAL MODEL SOURCE"
echo "======================================================================"

mkdir -p "$SRC"

rsync -a \
    "$PKG/model/source/" \
    "$SRC/"

echo "Source synced."
echo

# ----------------------------------------------------------------------
# Compile restart-capable executable
# ----------------------------------------------------------------------

echo "======================================================================"
echo "BUILD MODEL"
echo "======================================================================"

JOBS="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 4)"

cmake --build "$BUILD" \
    --target "$TARGET_NAME" \
    -j "$JOBS"

if [[ ! -x "$EXE" ]]; then
    echo "ERROR: executable not produced:"
    echo "  $EXE"
    exit 1
fi

echo
echo "Executable:"
echo "  $EXE"
echo

# ----------------------------------------------------------------------
# Protect against leftovers from an earlier aborted Mac attempt
# ----------------------------------------------------------------------

STAMP="$(date +%Y%m%d_%H%M%S)"
OLD="$PKG/mac_continue/old_attempt_$STAMP"
mkdir -p "$OLD"

for name in "$NAME1" "$NAME2"; do
    files=(
        "$RUN/${name}.pvd"
        "$RUN/${name}.log"
        "$RUN/${name}"*.vtu
    )

    if (( ${#files[@]} > 0 )); then
        for f in "${files[@]}"; do
            [[ -e "$f" ]] && mv "$f" "$OLD/"
        done
    fi
done

# ----------------------------------------------------------------------
# Common model options
# ----------------------------------------------------------------------

COMMON=(
    "$PARAMS"
    -Grid.File "$MESH"
    -Restart.Enable true
    -Restart.WriteAtEnd false

    -Problem.InitialOxygenProfileFile "$IC"
    -Problem.DrywellInjectionHead 5
    -Problem.DrywellStartupRampTime 1800
    -Problem.EnableGasPhaseOxygenAdvection true

    -Problem.Cycle1ShutdownStartHours 96.000000000000
    -Problem.Cycle2ShutdownStartHours 576.164044805955
    -Problem.Cycle3ShutdownStartHours 1060.953978421305

    -Dispersion.LongitudinalDispersivity 0.50
    -Dispersion.TransverseDispersivity 0.05
    -Dispersion.ApplyToGasPhase false

    -TimeLoop.DtInitial 1.08
    -TimeLoop.MaxTimeStepSize 300
    -TimeLoop.OutputInterval 7200

    -Newton.MaxRelativeShift 1e-4
    -Newton.MaxSteps 40
    -Newton.MaxTimeStepDivisions 15
    -Newton.TargetSteps 6
    -Newton.UseLineSearch false
    -Newton.EnableChop true
    -Newton.ChopMaxPressureChange 2e4
    -Newton.ChopMaxSaturationChange 0.10
    -Newton.ChopMaxLiquidMoleFractionChange 1e-5
    -Newton.ChopMaxGasMoleFractionChange 0.05
    -Newton.ChopSaturationBuffer 0

    -PrimaryVariableSwitch.GasPhaseHoldIterations 48
    -PrimaryVariableSwitch.GasDisappearanceTolerance 0.01
    -PrimaryVariableSwitch.GasAppearanceTolerance 0.02
    -PrimaryVariableSwitch.LiquidAppearanceTolerance 0.02
    -PrimaryVariableSwitch.AppearanceSaturation 1e-4
    -PrimaryVariableSwitch.Verbosity 0

    -Vtk.Precision Float64
)

cd "$RUN"

# ----------------------------------------------------------------------
# PART 1
# 512.061111 -> 955 h
# Generates final-C2 drainage states that were previously deleted.
# ----------------------------------------------------------------------

echo
echo "======================================================================"
echo "PART 1 / 2"
echo "FINAL C2: 512.061111 -> 955 h"
echo "======================================================================"
echo

"$EXE" "${COMMON[@]}" \
    -Restart.File "$R512" \
    -Problem.Name "$NAME1" \
    -TimeLoop.TEnd 3438000 \
    2>&1 | tee "$RUN/${NAME1}.log"

if [[ ! -f "$RUN/${NAME1}.pvd" ]]; then
    echo "ERROR: Part 1 PVD was not created."
    exit 1
fi

echo
echo "Part 1 complete."
echo

# ----------------------------------------------------------------------
# PART 2
# 1030 -> 1440 h
# Generates final C3 shutdown + complete final drainage.
# ----------------------------------------------------------------------

echo
echo "======================================================================"
echo "PART 2 / 2"
echo "FINAL C3: 1030 -> 1440 h"
echo "======================================================================"
echo

"$EXE" "${COMMON[@]}" \
    -Restart.File "$R1030" \
    -Problem.Name "$NAME2" \
    -TimeLoop.TEnd 5184000 \
    2>&1 | tee "$RUN/${NAME2}.log"

if [[ ! -f "$RUN/${NAME2}.pvd" ]]; then
    echo "ERROR: Part 2 PVD was not created."
    exit 1
fi

echo
echo "Part 2 complete."
echo

# ----------------------------------------------------------------------
# Move raw Mac runs into package.
# Nothing is deleted.
# ----------------------------------------------------------------------

echo "======================================================================"
echo "ARCHIVE RAW MAC RUNS"
echo "======================================================================"

rm -rf "$RAW1" "$RAW2"
mkdir -p "$RAW1" "$RAW2"

mv "$RUN/${NAME1}.pvd" "$RAW1/"
mv "$RUN/${NAME1}.log" "$RAW1/"
for f in "$RUN/${NAME1}"*.vtu; do
    mv "$f" "$RAW1/"
done

mv "$RUN/${NAME2}.pvd" "$RAW2/"
mv "$RUN/${NAME2}.log" "$RAW2/"
for f in "$RUN/${NAME2}"*.vtu; do
    mv "$f" "$RAW2/"
done

echo "Raw runs archived."
echo

# ----------------------------------------------------------------------
# Merge into FINAL 0-1440 h / nominal 6-h dataset.
#
# Existing:
#   0-576       = 97 states
#   960-1026    = 12 states
#
# Part 1:
#   582-954     = 63 states
#
# Part 2:
#   1032-1440   = 69 states
#
# Total = 241
# ----------------------------------------------------------------------

echo "======================================================================"
echo "BUILD FINAL 241-STATE DATASET"
echo "======================================================================"

python3 - "$PKG" "$NAME1" "$NAME2" <<'PY'
from pathlib import Path
import csv
import os
import shutil
import sys
import xml.etree.ElementTree as ET

PKG = Path(sys.argv[1])
NAME1 = sys.argv[2]
NAME2 = sys.argv[3]

EXISTING_PVD = (
    PKG /
    "existing_verified_6h" /
    "floodmar_existing_verified_6h.pvd"
)

P1_PVD = (
    PKG /
    "raw_mac_runs" /
    "part1_512_to955_2h" /
    f"{NAME1}.pvd"
)

P2_PVD = (
    PKG /
    "raw_mac_runs" /
    "part2_1030_to1440_2h" /
    f"{NAME2}.pvd"
)

FINAL = PKG / "final_60d_6h"
FINAL.mkdir(parents=True, exist_ok=True)

# Clean only a prior FINAL assembly.
# Raw simulations and original verified data are untouched.
for p in FINAL.glob("*.vtu"):
    p.unlink()

for p in [
    FINAL / "floodmar_equalVolume_60d_final_6h.pvd",
    FINAL / "OUTPUT_TIME_MAPPING.csv",
    FINAL / "README.txt",
]:
    if p.exists():
        p.unlink()

def read_pvd(path):
    if not path.exists():
        raise FileNotFoundError(path)

    root = ET.parse(path).getroot()

    raw = []
    for ds in root.findall(".//DataSet"):
        raw.append((
            float(ds.attrib["timestep"]),
            (path.parent / ds.attrib["file"]).resolve()
        ))

    if not raw:
        raise RuntimeError(f"No datasets in {path}")

    raw.sort(key=lambda x: x[0])

    # DuMuX PVD time in these models is seconds.
    # Also supports hour-valued PVDs if encountered.
    use_seconds = raw[-1][0] > 2000.0

    rows = [
        (
            t / 3600.0 if use_seconds else t,
            f
        )
        for t, f in raw
    ]

    missing = [f for _, f in rows if not f.exists()]
    if missing:
        raise RuntimeError(
            f"{path.name}: {len(missing)} referenced VTUs are missing"
        )

    print(
        f"{path.name}: "
        f"{len(rows)} records, "
        f"{rows[0][0]:.9f} -> {rows[-1][0]:.9f} h"
    )

    return rows

existing = read_pvd(EXISTING_PVD)
p1 = read_pvd(P1_PVD)
p2 = read_pvd(P2_PVD)

segments = [
    (
        "existing_0_576",
        list(range(0, 577, 6)),
        existing,
    ),
    (
        "mac_part1_final_C2",
        list(range(582, 955, 6)),
        p1,
    ),
    (
        "existing_960_1026",
        list(range(960, 1027, 6)),
        existing,
    ),
    (
        "mac_part2_final_C3",
        list(range(1032, 1441, 6)),
        p2,
    ),
]

expected_counts = [97, 63, 12, 69]

for (_, targets, _), expected in zip(segments, expected_counts):
    if len(targets) != expected:
        raise RuntimeError(
            f"Internal target-count error: "
            f"{len(targets)} != {expected}"
        )

vtk = ET.Element(
    "VTKFile",
    {
        "type": "Collection",
        "version": "0.1",
        "byte_order": "LittleEndian",
    },
)

collection = ET.SubElement(vtk, "Collection")

mapping = []
index = 0
previous_actual = -1e99

def link_or_copy(src, dst):
    try:
        os.link(src, dst)
        return "hardlink"
    except OSError:
        shutil.copy2(src, dst)
        return "copy"

for segment_name, targets, rows in segments:

    times = [x[0] for x in rows]

    for target in targets:

        actual, src = min(
            rows,
            key=lambda x: abs(x[0] - target)
        )

        diff_h = actual - target
        diff_min = diff_h * 60.0

        # Part 1 restart is at 512.061111 h, so nominal
        # integer-hour states are expected to be offset by
        # about +3.667 min.
        if abs(diff_min) > 10.0:
            raise RuntimeError(
                f"{segment_name}: target {target} h "
                f"nearest state={actual:.9f} h "
                f"({diff_min:+.3f} min)"
            )

        if actual <= previous_actual:
            raise RuntimeError(
                "Final time sequence is not strictly increasing"
            )

        previous_actual = actual

        outname = (
            f"floodmar_equalVolume_60d_"
            f"t{target:04d}h_"
            f"{index:04d}.vtu"
        )

        dst = FINAL / outname

        method = link_or_copy(src, dst)

        ET.SubElement(
            collection,
            "DataSet",
            timestep=f"{actual*3600.0:.12f}",
            group="",
            part="0",
            file=outname,
        )

        mapping.append([
            index,
            target,
            f"{actual:.12f}",
            f"{diff_min:+.6f}",
            segment_name,
            src.name,
            outname,
            method,
        ])

        index += 1

if index != 241:
    raise RuntimeError(
        f"Expected 241 final states; assembled {index}"
    )

PVD_OUT = FINAL / "floodmar_equalVolume_60d_final_6h.pvd"

ET.ElementTree(vtk).write(
    PVD_OUT,
    encoding="utf-8",
    xml_declaration=True,
)

with (FINAL / "OUTPUT_TIME_MAPPING.csv").open(
    "w",
    newline=""
) as f:
    w = csv.writer(f)
    w.writerow([
        "index",
        "nominal_time_h",
        "actual_simulation_time_h",
        "difference_min",
        "source_segment",
        "source_vtu",
        "final_vtu",
        "storage_method",
    ])
    w.writerows(mapping)

n = len(list(FINAL.glob("*.vtu")))

if n != 241:
    raise RuntimeError(
        f"Expected 241 VTUs in final folder; found {n}"
    )

actual_times = [
    float(r[2])
    for r in mapping
]

max_offset = max(
    abs(float(r[3]))
    for r in mapping
)

(FINAL / "README.txt").write_text(
f"""FINAL FLOOD MAR SAME-VOLUME DATASET

Simulation:
    0 to 1440 h = 60 days

Final schedule:
    C1 shutdown start = 96.000000000000 h
    C2 shutdown start = 576.164044805955 h
    C3 shutdown start = 1060.953978421305 h

Target injection:
    3956.009528 m3/cycle

Archive:
    nominal interval = 6 h
    states = 241

Sources:
    0-576 h:
        existing verified final-compatible trajectory

    582-954 h:
        Mac continuation from restart 512.061111111 h
        using final C2 schedule

    960-1026 h:
        existing trajectory generated after final C2
        and before final C3 shutdown

    1032-1440 h:
        Mac continuation from restart 1030 h
        using final C3 schedule

IMPORTANT:
    PVD timesteps preserve ACTUAL simulation times.
    The 512.061111-h restart causes the Part-1 states
    to differ from nominal integer-hour targets by about
    3.667 minutes.

    See OUTPUT_TIME_MAPPING.csv for exact times.

Maximum nominal/actual offset:
    {max_offset:.6f} min
"""
)

print()
print("="*78)
print("FINAL DATASET COMPLETE")
print("="*78)
print("Final folder :", FINAL)
print("VTUs         :", n)
print("PVD          :", PVD_OUT.name)
print("First state  :", f"{actual_times[0]:.9f} h")
print("Last state   :", f"{actual_times[-1]:.9f} h")
print("Max offset   :", f"{max_offset:.6f} min")
print()
print("Expected final VTUs = 241")
print("Actual final VTUs   =", n)
print()
print("PASS")
PY

echo
echo "======================================================================"
echo "ALL DONE"
echo "======================================================================"
echo
echo "FINAL DATASET:"
echo "  $FINAL"
echo
echo "Expected:"
echo "  241 VTUs"
echo "  1 final PVD"
echo "  OUTPUT_TIME_MAPPING.csv"
echo
echo "Raw Mac runs are also retained under:"
echo "  $RAW"
echo
echo "NO FINAL-TRAJECTORY DATA WAS DELETED."
