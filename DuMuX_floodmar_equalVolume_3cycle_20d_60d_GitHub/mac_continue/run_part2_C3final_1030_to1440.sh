#!/usr/bin/env bash
set -euo pipefail

PKG="$(cd "$(dirname "$0")/.." && pwd)"
ROOT="$HOME/dumux-work/dumux/dumux-floodmar"
RUN="$ROOT/build-cmake/test/porousmediumflow/2p2c/floodmar_pond5m_head5cm_equalDose_3cycle_20d_oldformal"

EXE="$RUN/floodmar_pond5m_head5cm_equalDose_3cycle_20d_restartwriter"

"$EXE" "$PKG/model/params.input" \
  -Grid.File "$PKG/model/floodmar_pond5m_depth10cm_nonames.msh" \
  -Restart.Enable true \
  -Restart.File "$PKG/restart/floodmar_C3continuous_restart_1030h.dat" \
  -Restart.WriteAtEnd false \
  -Problem.Name floodmar_FINAL_C3_1030_to1440_2h \
  -Problem.InitialOxygenProfileFile "$PKG/model/oxygenIC_trueStagnantO2_2880h_2d_xO2.dat" \
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
  -TimeLoop.OutputInterval 7200 \
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
  -Vtk.Precision Float64
