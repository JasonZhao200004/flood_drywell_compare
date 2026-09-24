#!/usr/bin/env bash
set -euo pipefail

EXE="./drywell_head12m_reservoir_4d6d"
CFG="./model/config"

"$EXE" \
    "$CFG/params.input" \
    -Grid.File "$CFG/floodmar_nonames.msh" \
    -Restart.Enable false \
    -Problem.Name reservoir_production_240h \
    -Problem.InitialOxygenProfileFile "$CFG/oxygenIC_trueStagnantO2_2880h_2d_xO2.dat" \
    -Problem.DrywellInjectionHead 1200 \
    -Problem.DrywellStartupRampTime 1800 \
    -Problem.EnableGasPhaseOxygenAdvection true \
    -Dispersion.LongitudinalDispersivity 0.50 \
    -Dispersion.TransverseDispersivity 0.05 \
    -Dispersion.ApplyToGasPhase false \
    -Dispersion.IncludeDrywellBoundaryVelocity false \
    -TimeLoop.TEnd 864000 \
    -TimeLoop.DtInitial 1.08 \
    -TimeLoop.MaxTimeStepSize 360 \
    -TimeLoop.OutputInterval 3600 \
    -Newton.MaxRelativeShift 1e-4 \
    -Newton.MaxSteps 18 \
    -Newton.MaxTimeStepDivisions 15 \
    -Newton.TargetSteps 6 \
    -Newton.UseLineSearch false \
    -Newton.EnableChop true \
    -Newton.ChopMaxPressureChange 2e4 \
    -Newton.ChopMaxSaturationChange 0.10 \
    -Newton.ChopMaxLiquidMoleFractionChange 1e-5 \
    -Newton.ChopMaxGasMoleFractionChange 0.05 \
    -Newton.ChopSaturationBuffer 0 \
    -PrimaryVariableSwitch.GasAppearanceTolerance 0.02 \
    -PrimaryVariableSwitch.LiquidAppearanceTolerance 0.02 \
    -PrimaryVariableSwitch.AppearanceSaturation 1e-4 \
    -PrimaryVariableSwitch.Verbosity 0 \
    -LinearSolver.ResidualReduction 1e-8 \
    -Vtk.Precision Float64
