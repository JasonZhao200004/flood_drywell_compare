#!/usr/bin/env bash

set -Eeuo pipefail

# =============================================================================
# DRYWELL HYDROFIX — NEW FROM-SCRATCH 4-CYCLE EQUAL-VOLUME PIPELINE
#
# C1 =    0 ->  480 h
# C2 =  480 ->  960 h
# C3 =  960 -> 1440 h
# C4 = 1440 -> 1920 h
#
# Every cycle target:
# 3956.009528 m3 gross well -> soil
#
# Final acceptance tolerance:
# +/- 0.02 m3 per cycle
# =============================================================================

ROOT="$HOME/dumux-work/dumux/dumux-floodmar"

SRC="$ROOT/test/porousmediumflow/2p2c/drywell_head12m_reservoir_from_production"

BUILD="$ROOT/build-cmake"

PROBLEM="$SRC/problem_reservoir_equalVolume_3cycle_60d_HYDROFIX.hh"

PARAMS="$SRC/params.input"

MESH="$SRC/floodmar_nonames.msh"

IC="$SRC/oxygenIC_trueStagnantO2_2880h_2d_xO2.dat"

TARGET_NAME="drywell_head12m_HYDROFIX_cycle1_0to480h"

TARGET_M3="3956.009528"

TOL_M3="0.02"

OUT="/mnt/d/flood_drywell_compare/DuMuX_drywell_head12m_HYDROFIX_4cycle_equalVolume_80d_FROM0"

echo
echo "===================================================================================================="
echo "NEW DRYWELL HYDROFIX 4-CYCLE EQUAL-VOLUME MODEL"
echo "===================================================================================================="
echo "C1 :    0 ->  480 h"
echo "C2 :  480 ->  960 h"
echo "C3 :  960 -> 1440 h"
echo "C4 : 1440 -> 1920 h"
echo
echo "Target per cycle = $TARGET_M3 m3 gross well -> soil"
echo "Final tolerance  = +/- $TOL_M3 m3"
echo "Head             = 1200 cm"
echo "Startup ramp     = 1800 s"
echo

# -----------------------------------------------------------------------------
# SAFETY CHECKS
# -----------------------------------------------------------------------------

for f in "$PROBLEM" "$PARAMS" "$MESH" "$IC"; do
    if [ ! -e "$f" ]; then
        echo "ERROR: missing required file:"
        echo "  $f"
        exit 1
    fi
done

# Archive an old attempt instead of deleting it.
if [ -e "$OUT" ]; then
    OLD="${OUT}_backup_$(date +%Y%m%d_%H%M%S)"
    echo "Existing output found."
    echo "Moving to:"
    echo "  $OLD"
    mv "$OUT" "$OLD"
fi

mkdir -p \
    "$OUT/bin" \
    "$OUT/model_source" \
    "$OUT/calibration" \
    "$OUT/restarts" \
    "$OUT/scripts"

cp "$0" "$OUT/scripts/" 2>/dev/null || true

# -----------------------------------------------------------------------------
# BUILD A COMPLETELY SEPARATE 4-CYCLE EXECUTABLE
#
# Temporarily patch the trusted HYDROFIX header, compile, copy executable,
# then restore original header immediately.
# -----------------------------------------------------------------------------

BACKUP="/tmp/problem_HYDROFIX_original_4cycle_$$.hh"

cp -f "$PROBLEM" "$BACKUP"

restore_problem()
{
    if [ -f "$BACKUP" ]; then
        cp -f "$BACKUP" "$PROBLEM"
        rm -f "$BACKUP"
        echo
        echo "Trusted HYDROFIX problem source restored."
    fi
}

trap restore_problem EXIT INT TERM

echo
echo "===================================================================================================="
echo "1. CREATE 4-CYCLE PARAMETERIZED HYDROFIX SOURCE"
echo "===================================================================================================="

python3 - "$PROBLEM" "$OUT/model_source/problem_HYDROFIX_4cycle_equalVolume.hh" <<'PY'
from pathlib import Path
import re
import sys

p = Path(sys.argv[1])
archive = Path(sys.argv[2])

s = p.read_text()

MARK = "HYDROFIX_4CYCLE_EQUAL_VOLUME_FROM0"

if MARK in s:
    raise SystemExit("ERROR: unexpected existing 4-cycle patch in trusted source")


def bounds(text, signature):
    pos = text.find(signature)
    if pos < 0:
        raise RuntimeError(f"Cannot find function: {signature}")

    op = text.find("{", pos)
    if op < 0:
        raise RuntimeError(f"Cannot find opening brace: {signature}")

    depth = 0

    for i in range(op, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1

            if depth == 0:
                return pos, op, i

    raise RuntimeError(f"Cannot find closing brace: {signature}")


def replace_body(text, signature, body):
    pos, op, cl = bounds(text, signature)

    return (
        text[:op+1]
        + "\n"
        + body.rstrip()
        + "\n    "
        + text[cl:]
    )


# ============================================================================
# Helper: runtime-controlled externally supplied wetting windows
# ============================================================================

cycle_pos = s.find("    int cycleNumber(Scalar timeSeconds) const")

if cycle_pos < 0:
    raise RuntimeError("cycleNumber() not found")

helper = r'''
    // HYDROFIX_4CYCLE_EQUAL_VOLUME_FROM0
    //
    // Four 20-day cycles. Wetting cutoff times are runtime parameters
    // calibrated from model-native gross well->soil injection.
    bool externallySupplied4Cycle_(Scalar timeSeconds) const
    {
        constexpr Scalar eps = 1.0e-8;

        const Scalar c1End =
            getParam<Scalar>(
                "Problem.Cycle1WetEndHours",
                163.125247592924
            )*3600.0;

        const Scalar c2End =
            getParam<Scalar>(
                "Problem.Cycle2WetEndHours",
                737.084098745400
            )*3600.0;

        const Scalar c3End =
            getParam<Scalar>(
                "Problem.Cycle3WetEndHours",
                1233.636146520033
            )*3600.0;

        const Scalar c4End =
            getParam<Scalar>(
                "Problem.Cycle4WetEndHours",
                1750.0
            )*3600.0;

        if (timeSeconds <= c1End + eps)
            return true;

        if (
            timeSeconds > 480.0*3600.0 + eps
            && timeSeconds <= c2End + eps
        )
            return true;

        if (
            timeSeconds > 960.0*3600.0 + eps
            && timeSeconds <= c3End + eps
        )
            return true;

        if (
            timeSeconds > 1440.0*3600.0 + eps
            && timeSeconds <= c4End + eps
        )
            return true;

        return false;
    }

    Scalar cycleStartSeconds4Cycle_(Scalar timeSeconds) const
    {
        constexpr Scalar eps = 1.0e-8;

        if (timeSeconds > 1440.0*3600.0 + eps)
            return 1440.0*3600.0;

        if (timeSeconds > 960.0*3600.0 + eps)
            return 960.0*3600.0;

        if (timeSeconds > 480.0*3600.0 + eps)
            return 480.0*3600.0;

        return 0.0;
    }

'''

s = s[:cycle_pos] + helper + s[cycle_pos:]


# ============================================================================
# cycleNumber(): 1,2,3,4
# ============================================================================

s = replace_body(
    s,
    "int cycleNumber(Scalar timeSeconds) const",
r'''
        constexpr Scalar eps = 1.0e-8;

        if (timeSeconds <= 480.0*3600.0 + eps)
            return 1;

        if (timeSeconds <= 960.0*3600.0 + eps)
            return 2;

        if (timeSeconds <= 1440.0*3600.0 + eps)
            return 3;

        return 4;
'''
)


# ============================================================================
# appliedDrywellHeadCm(): SAME 1800-s ramp, reset at every cycle
# ============================================================================

s = replace_body(
    s,
    "Scalar appliedDrywellHeadCm(Scalar timeSeconds) const",
r'''
        const Scalar cycleStart =
            cycleStartSeconds4Cycle_(timeSeconds);

        const Scalar localWettingTime =
            std::max(
                Scalar(0.0),
                timeSeconds - cycleStart
            );

        const Scalar startupFactor =
            drywellStartupRampTime_ > 0.0
            ? std::clamp(
                  localWettingTime/drywellStartupRampTime_,
                  Scalar(0.0),
                  Scalar(1.0)
              )
            : Scalar(1.0);

        const Scalar scheduledHeadCm =
            drywellHeadCm_(timeSeconds);

        return scheduledHeadCm > 0.0
               ? startupFactor*scheduledHeadCm
               : scheduledHeadCm;
'''
)


# ============================================================================
# drywellHeadCm_(): runtime equal-volume schedule
# ============================================================================

s = replace_body(
    s,
    "Scalar drywellHeadCm_(Scalar timeSeconds) const",
r'''
        if (externallySupplied4Cycle_(timeSeconds))
            return drywellInjectionHeadCm_;

        // After external supply stops, the finite drywell reservoir controls
        // the remaining water head exactly as in the trusted HYDROFIX model.
        return reservoirHeadCm_;
'''
)


# ============================================================================
# nextBoundaryEventTime():
# synchronize solver with:
# - 4 cycle starts
# - end of every 1800-s startup ramp
# - all four calibrated supply cutoffs
# ============================================================================

s = replace_body(
    s,
    "Scalar nextBoundaryEventTime(Scalar currentTimeSeconds) const",
r'''
        constexpr Scalar eps = 1.0e-8;

        const Scalar c1 =
            getParam<Scalar>(
                "Problem.Cycle1WetEndHours",
                163.125247592924
            )*3600.0;

        const Scalar c2 =
            getParam<Scalar>(
                "Problem.Cycle2WetEndHours",
                737.084098745400
            )*3600.0;

        const Scalar c3 =
            getParam<Scalar>(
                "Problem.Cycle3WetEndHours",
                1233.636146520033
            )*3600.0;

        const Scalar c4 =
            getParam<Scalar>(
                "Problem.Cycle4WetEndHours",
                1750.0
            )*3600.0;

        const Scalar events[] = {
            drywellStartupRampTime_,
            c1,

            480.0*3600.0,
            480.0*3600.0 + drywellStartupRampTime_,
            c2,

            960.0*3600.0,
            960.0*3600.0 + drywellStartupRampTime_,
            c3,

            1440.0*3600.0,
            1440.0*3600.0 + drywellStartupRampTime_,
            c4,

            1920.0*3600.0
        };

        Scalar next = Scalar(1.0e100);

        for (const Scalar eventTime : events)
        {
            if (
                eventTime > currentTimeSeconds + eps
                && eventTime < next
            )
                next = eventTime;
        }

        return next;
'''
)


# ============================================================================
# Robin drywell BC:
# Replace ONLY the scheduled/ramped-head expression.
# This regex is deliberately narrow and cannot cross into another function.
# ============================================================================

if "void updateReservoirAfterAcceptedStep(" not in s:
    raise RuntimeError(
        "updateReservoirAfterAcceptedStep() is already missing BEFORE Robin patch"
    )

pat = re.compile(
    "const\\s+Scalar\\s+scheduledDrywellHeadCm\\s*=\\s*"
    "drywellHeadCm_\\(boundaryTime\\)\\s*;\\s*"
    "const\\s+Scalar\\s+rampedDrywellHeadCm\\s*=\\s*"
    "scheduledDrywellHeadCm\\s*>\\s*0\\.0\\s*"
    "\\?\\s*startupFactor\\s*\\*\\s*scheduledDrywellHeadCm\\s*"
    ":\\s*scheduledDrywellHeadCm\\s*;"
)

replacement = (
    "const Scalar rampedDrywellHeadCm =\\n"
    "                appliedDrywellHeadCm(boundaryTime);"
)

s2, n = pat.subn(replacement, s, count=1)

if n != 1:
    raise RuntimeError(
        f"Expected exactly one Robin scheduled-head block; found {n}"
    )

if "void updateReservoirAfterAcceptedStep(" not in s2:
    raise RuntimeError(
        "SAFETY STOP: Robin replacement removed updateReservoirAfterAcceptedStep()"
    )

s = s2

print("Robin replacement PASS")
print("updateReservoirAfterAcceptedStep() preserved")

# ============================================================================
# updateReservoirAfterAcceptedStep():
# replace old fixed C1/C2/C3 supply test with the new four-cycle runtime test.
# Keep the trusted finite-reservoir drainage calculation unchanged.
# ============================================================================

sig = "void updateReservoirAfterAcceptedStep("

pos, op, cl = bounds(s, sig)

body = s[op+1:cl]

marker = "if (reservoirHeadCm_ <= 0.0)"

m = body.find(marker)

if m < 0:
    raise RuntimeError(
        "Could not locate finite-reservoir drainage section"
    )

drainage_tail = body[m:]

new_prefix = r'''

        // HYDROFIX_4CYCLE_EQUAL_VOLUME_FROM0
        //
        // While external supply is active, the drywell is maintained by the
        // prescribed 12-m source. Once the calibrated gross-volume target has
        // been reached, finite-reservoir drainage resumes.
        if (externallySupplied4Cycle_(acceptedEndTimeSeconds))
        {
            reservoirHeadCm_ = drywellInjectionHeadCm_;
            lastReservoirNetFluxM3s_ = 0.0;
            return;
        }

'''

new_body = new_prefix + drainage_tail

s = s[:op+1] + new_body + s[cl:]


# Final source sanity checks.
required = [
    "HYDROFIX_4CYCLE_EQUAL_VOLUME_FROM0",
    "Problem.Cycle1WetEndHours",
    "Problem.Cycle2WetEndHours",
    "Problem.Cycle3WetEndHours",
    "Problem.Cycle4WetEndHours",
    "return 4;",
    "1440.0*3600.0 + drywellStartupRampTime_",
    "appliedDrywellHeadCm(boundaryTime)",
]

for x in required:
    if x not in s:
        raise RuntimeError(
            f"Final source verification failed: {x}"
        )

p.write_text(s)
archive.write_text(s)

print("SUCCESS: four-cycle HYDROFIX source created.")
print()
print("Physics retained:")
print("  Drywell head       = runtime 1200 cm")
print("  Startup ramp       = runtime 1800 s")
print("  Finite reservoir   = retained")
print("  HYDROFIX pressure  = retained")
print("  O2 / dispersion    = retained")
print()
print("Cycle starts:")
print("  C1 =    0 h")
print("  C2 =  480 h")
print("  C3 =  960 h")
print("  C4 = 1440 h")
PY

echo
echo "===================================================================================================="
echo "2. COMPILE NEW FOUR-CYCLE EXECUTABLE"
echo "===================================================================================================="

cmake --build "$BUILD" \
    --target "$TARGET_NAME" \
    -j2

EXE_BUILD="$(
    find "$BUILD" \
        -type f \
        -name "$TARGET_NAME" \
        -perm -u+x \
        2>/dev/null \
        | head -1
)"

if [ -z "$EXE_BUILD" ] || [ ! -x "$EXE_BUILD" ]; then
    echo "ERROR: compiled executable not found."
    exit 1
fi

EXE="$OUT/bin/drywell_HYDROFIX_4cycle_equalVolume"

cp -f "$EXE_BUILD" "$EXE"
chmod +x "$EXE"

echo
echo "New executable:"
echo "  $EXE"

# Restore the trusted source immediately.
restore_problem
trap - EXIT INT TERM


# =============================================================================
# COMMON RUNNER
# =============================================================================

hours_to_seconds()
{
    python3 - "$1" <<'PY'
import sys
print(f"{float(sys.argv[1])*3600.0:.12f}")
PY
}


run_case()
{
    local stage="$1"
    local name="$2"
    local tend_h="$3"
    local restart_file="$4"
    local restart_out="$5"
    local c1="$6"
    local c2="$7"
    local c3="$8"
    local c4="$9"

    mkdir -p "$stage"

    local tend_s
    tend_s="$(hours_to_seconds "$tend_h")"

    local args=(
        "$EXE"
        "$PARAMS"

        -Grid.File "$MESH"

        -Problem.Name "$name"
        -Problem.InitialOxygenProfileFile "$IC"
        -Problem.DrywellInjectionHead 1200
        -Problem.DrywellStartupRampTime 1800
        -Problem.EnableGasPhaseOxygenAdvection true

        -Problem.Cycle1WetEndHours "$c1"
        -Problem.Cycle2WetEndHours "$c2"
        -Problem.Cycle3WetEndHours "$c3"
        -Problem.Cycle4WetEndHours "$c4"

        -Dispersion.LongitudinalDispersivity 0.50
        -Dispersion.TransverseDispersivity 0.05
        -Dispersion.ApplyToGasPhase false
        -Dispersion.IncludeDrywellBoundaryVelocity false

        -TimeLoop.TEnd "$tend_s"
        -TimeLoop.DtInitial 1.08
        -TimeLoop.MaxTimeStepSize 360
        -TimeLoop.OutputInterval 21600

        -Newton.MaxRelativeShift 1e-4
        -Newton.MaxSteps 18
        -Newton.MaxTimeStepDivisions 15
        -Newton.TargetSteps 6
        -Newton.UseLineSearch false

        -Newton.EnableChop true
        -Newton.ChopMaxPressureChange 2e4
        -Newton.ChopMaxSaturationChange 0.10
        -Newton.ChopMaxLiquidMoleFractionChange 1e-5
        -Newton.ChopMaxGasMoleFractionChange 0.05
        -Newton.ChopSaturationBuffer 0

        -PrimaryVariableSwitch.GasAppearanceTolerance 0.02
        -PrimaryVariableSwitch.LiquidAppearanceTolerance 0.02
        -PrimaryVariableSwitch.AppearanceSaturation 1e-4
        -PrimaryVariableSwitch.Verbosity 0

        -LinearSolver.ResidualReduction 1e-8

        -Vtk.Precision Float64
    )

    if [ -n "$restart_file" ]; then
        args+=(
            -Restart.Enable true
            -Restart.File "$restart_file"
        )
    else
        args+=(
            -Restart.Enable false
        )
    fi

    if [ -n "$restart_out" ]; then
        args+=(
            -Restart.WriteAtEnd true
            -Restart.OutputFile "$restart_out"
        )
    else
        args+=(
            -Restart.WriteAtEnd false
        )
    fi

    echo
    echo "RUN:"
    echo "  $name"
    echo "  TEnd = $tend_h h"
    echo "  C1 cutoff = $c1"
    echo "  C2 cutoff = $c2"
    echo "  C3 cutoff = $c3"
    echo "  C4 cutoff = $c4"

    if ! (
        cd "$stage"
        "${args[@]}"
    ) > "$stage/run.log" 2>&1
    then
        echo
        echo "===================================================================================================="
        echo "MODEL FAILED: $name"
        echo "===================================================================================================="
        tail -100 "$stage/run.log"
        exit 1
    fi
}


# Remove calibration VTKs after their information has been extracted.
# Keep logs + injection_history CSV + restart snapshots.
cleanup_stage_vtk()
{
    local d="$1"

    rm -f \
        "$d"/*.vtu \
        "$d"/*.pvd \
        2>/dev/null || true
}


# =============================================================================
# PYTHON HELPERS FOR MODEL-NATIVE GROSS WELL->SOIL VOLUME
# =============================================================================

find_crossing()
{
    local csv="$1"
    local start_h="$2"

    python3 - "$csv" "$start_h" "$TARGET_M3" <<'PY'
import csv
import sys

path = sys.argv[1]
start = float(sys.argv[2])
target = float(sys.argv[3])

rows = []

with open(path, newline="") as f:
    r = csv.DictReader(f)

    if "time_h" not in (r.fieldnames or []):
        raise SystemExit("ERROR: time_h missing")

    if "step_well_to_soil_m3" not in (r.fieldnames or []):
        raise SystemExit(
            "ERROR: step_well_to_soil_m3 missing from model-native history"
        )

    for row in r:
        t = float(row["time_h"])

        if t <= start + 1e-9:
            continue

        step = max(
            0.0,
            float(row["step_well_to_soil_m3"])
        )

        rows.append((t, step))

cum = 0.0
prev_t = start
prev_v = 0.0

for t, step in rows:
    new_v = cum + step

    if new_v >= target:
        dv = new_v - cum

        if dv <= 0.0:
            tc = t
        else:
            f = (target - cum)/dv
            f = min(1.0, max(0.0, f))
            tc = prev_t + f*(t - prev_t)

        print(f"{tc:.12f}")
        raise SystemExit(0)

    prev_t = t
    prev_v = cum
    cum = new_v

raise SystemExit(3)
PY
}


measure_at_cutoff()
{
    local csv="$1"
    local start_h="$2"
    local cutoff_h="$3"

    python3 - "$csv" "$start_h" "$cutoff_h" <<'PY'
import csv
import sys

path = sys.argv[1]
start = float(sys.argv[2])
cut = float(sys.argv[3])

rows = []

with open(path, newline="") as f:
    r = csv.DictReader(f)

    fields = set(r.fieldnames or [])

    if "step_well_to_soil_m3" not in fields:
        raise SystemExit(
            "ERROR: missing step_well_to_soil_m3"
        )

    for row in r:
        t = float(row["time_h"])

        if not (
            t > start + 1e-9
            and t <= cut + 1e-7
        ):
            continue

        step = max(
            0.0,
            float(row["step_well_to_soil_m3"])
        )

        dt = float(row.get("dt_s", "0") or 0)

        rows.append((t, step, dt))

vol = sum(x[1] for x in rows)

# Estimate near-cutoff injection rate from the final ~1 h.
near = [
    x for x in rows
    if x[0] >= cut - 1.0
]

sv = sum(x[1] for x in near)
sd = sum(x[2] for x in near)

if sd > 0.0:
    rate_h = sv/sd*3600.0
else:
    rate_h = 0.0

print(f"{vol:.12f} {rate_h:.12f}")
PY
}


abs_le_tol()
{
    python3 - "$1" "$2" <<'PY'
import sys
print(
    "yes"
    if abs(float(sys.argv[1])) <= float(sys.argv[2])
    else "no"
)
PY
}


adjust_cutoff()
{
    python3 - "$1" "$2" "$3" "$4" <<'PY'
import sys

cut = float(sys.argv[1])
err = float(sys.argv[2])
rate = float(sys.argv[3])
cycle_end = float(sys.argv[4])

if rate <= 0.0:
    raise SystemExit("ERROR: non-positive near-cutoff injection rate")

new = cut + err/rate

new = max(0.0, min(cycle_end - 0.01, new))

print(f"{new:.12f}")
PY
}


# =============================================================================
# CALIBRATION STATE
# =============================================================================

C1="479.999"
C2="959.999"
C3="1439.999"
C4="1919.999"

STARTS=(0 480 960 1440)
ENDS=(480 960 1440 1920)

CURRENT_RESTART=""

echo
echo "===================================================================================================="
echo "3. SEQUENTIAL HYDROFIX CALIBRATION — C1, C2, C3"
echo "===================================================================================================="


for IDX in 0 1 2
do
    CYCLE=$((IDX + 1))
    START="${STARTS[$IDX]}"
    END="${ENDS[$IDX]}"

    echo
    echo "----------------------------------------------------------------------------------------------------"
    echo "CALIBRATING CYCLE $CYCLE"
    echo "Cycle start = $START h"
    echo "Target      = $TARGET_M3 m3"
    echo "----------------------------------------------------------------------------------------------------"

    CAL_END="$(
        python3 - "$START" <<'PY'
import sys
print(float(sys.argv[1]) + 360.0)
PY
    )"

    STAGE="$OUT/calibration/C${CYCLE}_continuous"

    NAME="drywell_HYDROFIX_C${CYCLE}_continuous_calibration"

    run_case \
        "$STAGE" \
        "$NAME" \
        "$CAL_END" \
        "$CURRENT_RESTART" \
        "" \
        "$C1" "$C2" "$C3" "$C4"

    CSV="$STAGE/${NAME}_injection_history.csv"

    set +e
    CUT="$(find_crossing "$CSV" "$START")"
    RC=$?
    set -e

    if [ "$RC" -ne 0 ]; then
        echo "Target not reached in first calibration window."
        echo "Extending calibration close to end of cycle..."

        CAL_END="$(
            python3 - "$START" <<'PY'
import sys
print(float(sys.argv[1]) + 470.0)
PY
        )"

        STAGE="$OUT/calibration/C${CYCLE}_continuous_extended"

        run_case \
            "$STAGE" \
            "$NAME" \
            "$CAL_END" \
            "$CURRENT_RESTART" \
            "" \
            "$C1" "$C2" "$C3" "$C4"

        CSV="$STAGE/${NAME}_injection_history.csv"

        CUT="$(find_crossing "$CSV" "$START")"
    fi

    cleanup_stage_vtk "$STAGE"

    echo
    echo "Initial C${CYCLE} cutoff estimate = $CUT h"

    case "$CYCLE" in
        1) C1="$CUT" ;;
        2) C2="$CUT" ;;
        3) C3="$CUT" ;;
    esac

    # -------------------------------------------------------------------------
    # Refine cutoff against a formal run using the exact event-synchronized
    # HYDROFIX model.
    # -------------------------------------------------------------------------

    for ITER in 1 2 3 4
    do
        VALID_END="$(
            python3 - "$CUT" "$END" <<'PY'
import sys
c = float(sys.argv[1])
e = float(sys.argv[2])
print(min(e, c + 0.25))
PY
        )"

        VSTAGE="$OUT/calibration/C${CYCLE}_validation_${ITER}"

        VNAME="drywell_HYDROFIX_C${CYCLE}_validation_${ITER}"

        run_case \
            "$VSTAGE" \
            "$VNAME" \
            "$VALID_END" \
            "$CURRENT_RESTART" \
            "" \
            "$C1" "$C2" "$C3" "$C4"

        VCSV="$VSTAGE/${VNAME}_injection_history.csv"

        read VOL RATE < <(
            measure_at_cutoff \
                "$VCSV" \
                "$START" \
                "$CUT"
        )

        ERR="$(
            python3 - "$TARGET_M3" "$VOL" <<'PY'
import sys
print(
    f"{float(sys.argv[1])-float(sys.argv[2]):.12f}"
)
PY
        )"

        echo
        echo "C${CYCLE} validation $ITER"
        echo "  cutoff = $CUT h"
        echo "  volume = $VOL m3"
        echo "  target = $TARGET_M3 m3"
        echo "  target-volume = $ERR m3"
        echo "  local rate = $RATE m3/h"

        cleanup_stage_vtk "$VSTAGE"

        if [ "$(abs_le_tol "$ERR" "$TOL_M3")" = "yes" ]; then
            echo "  PASS"
            break
        fi

        CUT="$(
            adjust_cutoff \
                "$CUT" \
                "$ERR" \
                "$RATE" \
                "$END"
        )"

        case "$CYCLE" in
            1) C1="$CUT" ;;
            2) C2="$CUT" ;;
            3) C3="$CUT" ;;
        esac

        echo "  new cutoff = $CUT h"

        if [ "$ITER" -eq 4 ]; then
            echo "ERROR: C${CYCLE} cutoff refinement did not converge."
            exit 1
        fi
    done

    # -------------------------------------------------------------------------
    # Run the calibrated cycle all the way to its 20-day boundary and write an
    # exact primary-variable restart for calibration of the next cycle.
    # -------------------------------------------------------------------------

    RESTART_OUT="$OUT/restarts/restart_${END}h.dat"

    FSTAGE="$OUT/calibration/C${CYCLE}_formal_to_${END}h"

    FNAME="drywell_HYDROFIX_C${CYCLE}_formal"

    run_case \
        "$FSTAGE" \
        "$FNAME" \
        "$END" \
        "$CURRENT_RESTART" \
        "$RESTART_OUT" \
        "$C1" "$C2" "$C3" "$C4"

    FCSV="$FSTAGE/${FNAME}_injection_history.csv"

    read FVOL FRATE < <(
        measure_at_cutoff \
            "$FCSV" \
            "$START" \
            "$CUT"
    )

    FERR="$(
        python3 - "$FVOL" "$TARGET_M3" <<'PY'
import sys
print(
    f"{float(sys.argv[1])-float(sys.argv[2]):+.12f}"
)
PY
    )"

    echo
    echo "FINAL C${CYCLE} FORMAL CHECK"
    echo "  cutoff = $CUT h"
    echo "  gross  = $FVOL m3"
    echo "  error  = $FERR m3"

    if [ ! -s "$RESTART_OUT" ]; then
        echo "ERROR: restart was not written:"
        echo "$RESTART_OUT"
        exit 1
    fi

    cleanup_stage_vtk "$FSTAGE"

    CURRENT_RESTART="$RESTART_OUT"
done


# =============================================================================
# C4 CALIBRATION
#
# IMPORTANT:
# Do NOT restart at 1440 h for the actual C4 calibration.
# Run from t=0 continuously through C1-C2-C3 into C4.
# This deliberately avoids today's restart/C4 numerical problem.
# =============================================================================

echo
echo "===================================================================================================="
echo "4. C4 CALIBRATION — CONTINUOUS FROM t=0"
echo "===================================================================================================="

C4="1919.999"

C4_STAGE="$OUT/calibration/C4_continuous_FROM0"

C4_NAME="drywell_HYDROFIX_C4_continuous_FROM0"

run_case \
    "$C4_STAGE" \
    "$C4_NAME" \
    "1800.0" \
    "" \
    "" \
    "$C1" "$C2" "$C3" "$C4"

C4_CSV="$C4_STAGE/${C4_NAME}_injection_history.csv"

set +e
C4_CUT="$(find_crossing "$C4_CSV" "1440.0")"
RC=$?
set -e

if [ "$RC" -ne 0 ]; then
    echo "C4 target not reached by 1800 h."
    echo "Extending continuous-from-zero calibration to 1910 h."

    cleanup_stage_vtk "$C4_STAGE"

    C4_STAGE="$OUT/calibration/C4_continuous_FROM0_extended"

    run_case \
        "$C4_STAGE" \
        "$C4_NAME" \
        "1910.0" \
        "" \
        "" \
        "$C1" "$C2" "$C3" "$C4"

    C4_CSV="$C4_STAGE/${C4_NAME}_injection_history.csv"

    C4_CUT="$(find_crossing "$C4_CSV" "1440.0")"
fi

C4="$C4_CUT"

cleanup_stage_vtk "$C4_STAGE"

echo
echo "Initial C4 cutoff from continuous 0-h trajectory:"
echo "  C4 cutoff = $C4 h"


# =============================================================================
# FINAL CONTINUOUS 80-DAY PRODUCTION RUN
#
# After every full run, calculate ACTUAL model-native gross well->soil volume
# in all four wetting periods.
#
# If any differs from target by >0.02 m3, adjust the cutoff and repeat.
# =============================================================================

echo
echo "===================================================================================================="
echo "5. FINAL CONTINUOUS 0 -> 1920 h PRODUCTION + FOUR-CYCLE QC"
echo "===================================================================================================="

ACCEPTED="no"
FINAL_STAGE=""

for GLOBAL_ITER in 1 2 3 4
do
    echo
    echo "----------------------------------------------------------------------------------------------------"
    echo "FULL 80-d PRODUCTION ITERATION $GLOBAL_ITER"
    echo "----------------------------------------------------------------------------------------------------"
    echo "C1 cutoff = $C1"
    echo "C2 cutoff = $C2"
    echo "C3 cutoff = $C3"
    echo "C4 cutoff = $C4"

    STAGE="$OUT/final_iter_${GLOBAL_ITER}"

    NAME="drywell_HYDROFIX_4cycle_equalVolume_80d_iter${GLOBAL_ITER}"

    run_case \
        "$STAGE" \
        "$NAME" \
        "1920.0" \
        "" \
        "$STAGE/restart_1920h.dat" \
        "$C1" "$C2" "$C3" "$C4"

    CSV="$STAGE/${NAME}_injection_history.csv"

    QC="$(
        python3 - \
            "$CSV" \
            "$TARGET_M3" \
            "$C1" "$C2" "$C3" "$C4" <<'PY'
import csv
import sys

path = sys.argv[1]
target = float(sys.argv[2])

cuts = [
    float(sys.argv[3]),
    float(sys.argv[4]),
    float(sys.argv[5]),
    float(sys.argv[6]),
]

starts = [0.0, 480.0, 960.0, 1440.0]

rows = []

with open(path, newline="") as f:
    r = csv.DictReader(f)

    fields = set(r.fieldnames or [])

    if "step_well_to_soil_m3" not in fields:
        raise SystemExit(
            "ERROR: step_well_to_soil_m3 missing"
        )

    for row in r:
        rows.append(row)

vals = []
rates = []

for start, cut in zip(starts, cuts):

    selected = []

    for row in rows:
        t = float(row["time_h"])

        if (
            t > start + 1e-9
            and t <= cut + 1e-7
        ):
            selected.append(row)

    vol = sum(
        max(
            0.0,
            float(r["step_well_to_soil_m3"])
        )
        for r in selected
    )

    near = [
        r for r in selected
        if float(r["time_h"]) >= cut - 1.0
    ]

    sv = sum(
        max(
            0.0,
            float(r["step_well_to_soil_m3"])
        )
        for r in near
    )

    sd = sum(
        float(r.get("dt_s", "0") or 0)
        for r in near
    )

    rate = sv/sd*3600.0 if sd > 0 else 0.0

    vals.append(vol)
    rates.append(rate)

for v, q in zip(vals, rates):
    print(f"{v:.12f}:{q:.12f}", end=" ")

PY
    )"

    read P1 P2 P3 P4 <<< "$QC"

    V1="${P1%%:*}"
    R1="${P1##*:}"

    V2="${P2%%:*}"
    R2="${P2##*:}"

    V3="${P3%%:*}"
    R3="${P3##*:}"

    V4="${P4%%:*}"
    R4="${P4##*:}"

    E1="$(
        python3 - "$TARGET_M3" "$V1" <<'PY'
import sys
print(float(sys.argv[1])-float(sys.argv[2]))
PY
    )"

    E2="$(
        python3 - "$TARGET_M3" "$V2" <<'PY'
import sys
print(float(sys.argv[1])-float(sys.argv[2]))
PY
    )"

    E3="$(
        python3 - "$TARGET_M3" "$V3" <<'PY'
import sys
print(float(sys.argv[1])-float(sys.argv[2]))
PY
    )"

    E4="$(
        python3 - "$TARGET_M3" "$V4" <<'PY'
import sys
print(float(sys.argv[1])-float(sys.argv[2]))
PY
    )"

    echo
    echo "===================================================================================================="
    echo "FOUR-CYCLE INJECTION QC — ITERATION $GLOBAL_ITER"
    echo "===================================================================================================="
    printf "C1 : %.9f m3   error = %+ .9f m3\n" "$V1" "$E1"
    printf "C2 : %.9f m3   error = %+ .9f m3\n" "$V2" "$E2"
    printf "C3 : %.9f m3   error = %+ .9f m3\n" "$V3" "$E3"
    printf "C4 : %.9f m3   error = %+ .9f m3\n" "$V4" "$E4"
    echo
    echo "Target = $TARGET_M3 m3"
    echo "Tolerance = +/- $TOL_M3 m3"

    PASS="$(
        python3 - \
            "$E1" "$E2" "$E3" "$E4" "$TOL_M3" <<'PY'
import sys

errs = list(map(float, sys.argv[1:5]))
tol = float(sys.argv[5])

print(
    "yes"
    if max(abs(x) for x in errs) <= tol
    else "no"
)
PY
    )"

    if [ "$PASS" = "yes" ]; then
        echo
        echo "===================================================================================================="
        echo "PASS — ALL FOUR CYCLES HAVE THE SAME INJECTION VOLUME"
        echo "===================================================================================================="

        ACCEPTED="yes"
        FINAL_STAGE="$STAGE"
        break
    fi

    if [ "$GLOBAL_ITER" -eq 4 ]; then
        echo
        echo "ERROR:"
        echo "Four-cycle volume QC did not converge to +/- $TOL_M3 m3."
        echo "No result will be labeled FINAL."
        exit 1
    fi

    echo
    echo "Adjusting all four cutoff times from the ACTUAL continuous 80-d result..."

    C1="$(adjust_cutoff "$C1" "$E1" "$R1" "480")"
    C2="$(adjust_cutoff "$C2" "$E2" "$R2" "960")"
    C3="$(adjust_cutoff "$C3" "$E3" "$R3" "1440")"
    C4="$(adjust_cutoff "$C4" "$E4" "$R4" "1920")"

    # Keep the history/log, discard large VTKs from a failed QC iteration.
    cleanup_stage_vtk "$STAGE"
done


# =============================================================================
# FINALIZE
# =============================================================================

if [ "$ACCEPTED" != "yes" ]; then
    echo "ERROR: no accepted final run."
    exit 1
fi

FINAL="$OUT/FINAL_80d_4cycle_equalVolume"

if [ -e "$FINAL" ]; then
    mv "$FINAL" "${FINAL}_old_$(date +%Y%m%d_%H%M%S)"
fi

mv "$FINAL_STAGE" "$FINAL"

cat > "$OUT/FINAL_SCHEDULE.csv" <<EOF
cycle,start_h,end_h,wetting_cutoff_h,target_gross_well_to_soil_m3
C1,0,480,$C1,$TARGET_M3
C2,480,960,$C2,$TARGET_M3
C3,960,1440,$C3,$TARGET_M3
C4,1440,1920,$C4,$TARGET_M3
EOF

cat > "$OUT/README_EQUAL_VOLUME.txt" <<EOF
DRYWELL HYDROFIX — 4-CYCLE EQUAL-VOLUME RUN
===========================================

Model starts from t = 0.
No 1440-h restart is used in the FINAL production run.

Cycle boundaries:
C1 =    0-480 h
C2 =  480-960 h
C3 =  960-1440 h
C4 = 1440-1920 h

Gross well-to-soil target per cycle:
$TARGET_M3 m3

Acceptance tolerance:
+/- $TOL_M3 m3 per cycle

Hydraulic head:
1200 cm

Startup ramp:
1800 s, restarted at each cycle start.

Finite-reservoir drainage:
retained from trusted HYDROFIX implementation.

Final calibrated cutoff times:
C1 = $C1 h
C2 = $C2 h
C3 = $C3 h
C4 = $C4 h

Final result directory:
$FINAL
EOF

echo
echo "===================================================================================================="
echo "DONE"
echo "===================================================================================================="
echo
echo "FINAL model:"
echo "  $FINAL"
echo
echo "Schedule:"
echo "  C1 cutoff = $C1 h"
echo "  C2 cutoff = $C2 h"
echo "  C3 cutoff = $C3 h"
echo "  C4 cutoff = $C4 h"
echo
echo "Target per cycle:"
echo "  $TARGET_M3 m3"
echo
echo "Final schedule file:"
echo "  $OUT/FINAL_SCHEDULE.csv"
echo
echo "This FINAL run is continuous from 0 to 1920 h."
echo "No C4 restart is used."
echo

