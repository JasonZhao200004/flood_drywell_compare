#!/usr/bin/env bash
set -euo pipefail

# Unzip the variant, then run this script. The archived original model and
# completed 60-day outputs are left untouched.
VARIANT="$(cd "$(dirname "$0")/.." && pwd)"
ORIGINAL="$HOME/Documents/GitHub/floof_drywell_compare/DuMuX_floodmar_equalVolume_3cycle_20d_60d_GitHub_FINAL"

DUMUX_ROOT=""
for candidate in \
    "$HOME/dumux-work/dumux/dumux-floodmar" \
    "$HOME/Documents/GitHub/dumux-floodmar" \
    "$HOME/dumux-floodmar"
do
    if [[ -d "$candidate/build-cmake" ]]; then
        DUMUX_ROOT="$candidate"
        break
    fi
done
if [[ -z "$DUMUX_ROOT" ]]; then
    echo "DuMuX build-cmake directory not found" >&2
    exit 1
fi

SRC="$DUMUX_ROOT/test/porousmediumflow/2p2c/floodmar_pond5m_head5cm_equalDose_3cycle_20d_oldformal"
RUN="$DUMUX_ROOT/build-cmake/test/porousmediumflow/2p2c/floodmar_pond5m_head5cm_equalDose_3cycle_20d_oldformal"
TARGET="floodmar_pond5m_head5cm_equalDose_3cycle_20d_restartwriter"
MODEL="$ORIGINAL/model"
RESULTS="$VARIANT/results"

for file in \
    "$SRC/problem_floodpond5m_depth10cm.hh" \
    "$VARIANT/model/source/problem_floodpond5m_depth10cm.hh" \
    "$MODEL/params.input" \
    "$MODEL/floodmar_pond5m_depth10cm_nonames.msh" \
    "$MODEL/oxygenIC_trueStagnantO2_2880h_2d_xO2.dat"
do
    [[ -f "$file" ]] || { echo "Missing: $file" >&2; exit 1; }
done

mkdir -p "$RESULTS"
# Save the active build-tree header in this variant before replacing it.
if [[ ! -f "$VARIANT/problem_before_atmospheric_variant.hh" ]]; then
    cp "$SRC/problem_floodpond5m_depth10cm.hh" \
       "$VARIANT/problem_before_atmospheric_variant.hh"
fi
cp "$VARIANT/model/source/problem_floodpond5m_depth10cm.hh" \
   "$SRC/problem_floodpond5m_depth10cm.hh"

cmake --build "$DUMUX_ROOT/build-cmake" --target "$TARGET" -j 4

# The boundary differs after C1 dryout; old intermediate restarts are invalid.
# Always start this variant at t=0 and keep its output in its own directory.
cd "$RESULTS"
"$RUN/$TARGET" "$MODEL/params.input" \
    -Grid.File "$MODEL/floodmar_pond5m_depth10cm_nonames.msh" \
    -Restart.Enable false \
    -Restart.WriteAtEnd false \
    -Problem.Name floodmar_equalVolume_60d_drain_atmospheric \
    -Problem.InitialOxygenProfileFile "$MODEL/oxygenIC_trueStagnantO2_2880h_2d_xO2.dat" \
    -Problem.DrywellInjectionHead 5 \
    -Problem.DrywellStartupRampTime 1800 \
    -Problem.EnableGasPhaseOxygenAdvection true \
    -Problem.Cycle1ShutdownStartHours 96.000000000000 \
    -Problem.Cycle2ShutdownStartHours 576.164044805955 \
    -Problem.Cycle3ShutdownStartHours 1060.953978421305 \
    -Dispersion.LongitudinalDispersivity 0.50 \
    -Dispersion.TransverseDispersivity 0.05 \
    -Dispersion.ApplyToGasPhase false \
    -TimeLoop.TEnd 5184000 \
    -TimeLoop.DtInitial 1.08 \
    -TimeLoop.MaxTimeStepSize 300 \
    -TimeLoop.OutputInterval 21600 \
    -Newton.MaxRelativeShift 1e-4 \
    -Newton.MaxSteps 40 \
    -Newton.MaxTimeStepDivisions 15 \
    -Newton.TargetSteps 6 \
    -Newton.UseLineSearch false \
    -Newton.EnableChop true \
    -Newton.ChopMaxPressureChange 2e4 \
    -Newton.ChopMaxSaturationChange 0.10 \
    -Newton.ChopMaxLiquidMoleFractionChange 1e-5 \
    -Newton.ChopMaxGasMoleFractionChange 0.05 \
    -Newton.ChopSaturationBuffer 0 \
    -PrimaryVariableSwitch.GasPhaseHoldIterations 48 \
    -PrimaryVariableSwitch.GasDisappearanceTolerance 0.01 \
    -PrimaryVariableSwitch.GasAppearanceTolerance 0.02 \
    -PrimaryVariableSwitch.LiquidAppearanceTolerance 0.02 \
    -PrimaryVariableSwitch.AppearanceSaturation 1e-4 \
    -PrimaryVariableSwitch.Verbosity 0 \
    -Vtk.Precision Float64 \
    2>&1 | tee floodmar_equalVolume_60d_drain_atmospheric.log
