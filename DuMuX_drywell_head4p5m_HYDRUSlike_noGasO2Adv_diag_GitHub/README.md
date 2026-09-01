# DuMuX 4.5 m Drywell — HYDRUS-like No-Gas-O2-Advection Diagnostic

## Purpose

Diagnostic version of the validated 4.5 m-head DuMuX drywell model used to
isolate the effect of gas-phase O2 advection on the HYDRUS–DuMuX shallow
oxygen discrepancy.

This is not the formal production model. The formal production model retains
gas-phase O2 advection.

## Physics retained

- variably saturated liquid flow
- bulk gas Darcy flow and gas-pressure dynamics
- liquid-phase O2 advection
- liquid mechanical dispersion
- gas-phase O2 molecular diffusion
- equilibrium gas-liquid partitioning
- moisture-dependent O2 consumption
- PhaseGuard V2
- face-symmetric liquid dispersion

## Diagnostic processes disabled

- internal gas-phase O2 advective component transport
- surface pressure-driven advective O2 exchange
- gas-phase mechanical dispersion

## Implementation

Internal gas-phase O2 advection is suppressed using:

`floodmar_hydruslike_noGasO2Adv_localresidual.hh`

Surface advective O2 exchange is disabled at runtime using:

`Problem.EnableGasPhaseOxygenAdvection = false`

Bulk gas Darcy flow and gas-pressure dynamics remain active.

## Model configuration

- axisymmetric domain radius: 20 m
- domain depth: 30 m
- drywell radius: 0.61 m
- drywell depth: 16 m
- drywell bottom elevation: 14 m
- injection head: 4.5 m
- startup ramp: 1800 s
- longitudinal dispersivity: 0.50 m
- transverse dispersivity: 0.05 m
- gas mechanical dispersion: OFF

## Successful diagnostic simulation

Run:

`drywell_head4p5m_faceSymDisp_HYDRUSlike_allGasO2AdvOFF_diag240h_6h`

Duration:

- 240 h

Output interval:

- 6 h

Outputs:

- 41 VTU states including t = 0

## HYDRUS comparison

Near-well upper region:

- radius = 0–6 m
- depth = 0–6 m

At 96 h:

Formal gas-O2-advection-ON DuMuX:
- R2 = -0.4254
- RMSE = 2.1588e-6 g/cm3

HYDRUS-like diagnostic:
- R2 = 0.8204
- RMSE = 7.6624e-7 g/cm3

RMSE reduction = approximately 64.5%.

At 240 h:

Formal gas-O2-advection-ON DuMuX:
- R2 = -0.1592
- RMSE = 2.1415e-6 g/cm3

HYDRUS-like diagnostic:
- R2 = 0.8131
- RMSE = 8.5990e-7 g/cm3

RMSE reduction = approximately 59.9%.

These results indicate that gas-phase O2 advection is the dominant cause of
the localized shallow near-drywell HYDRUS–DuMuX oxygen discrepancy.

## Important

This configuration is a process-ablation diagnostic only.

Formal DuMuX production simulations retain gas-phase O2 advection.
