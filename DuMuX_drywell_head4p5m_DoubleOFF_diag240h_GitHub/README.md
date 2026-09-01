# DuMuX 4.5-m-head drywell — Double-OFF O2 diagnostic

## Purpose

Diagnostic comparison model derived from the validated 4.5-m-head
PhaseGuard V2 + face-symmetric dispersion drywell model.

This is the model referred to during development as the "Double-OFF"
case.

## Important physics definition

"Double-OFF" does NOT mean that all advection is disabled.

The diagnostic configuration is:

- liquid Darcy flow: ON
- liquid-phase O2 advection: ON
- liquid mechanical dispersion: ON
- bulk gas Darcy flow / gas-pressure dynamics: ON
- internal gas-phase O2 advection: OFF
- surface advective O2 exchange: OFF
- gas molecular O2 diffusion: ON
- gas mechanical dispersion: OFF

Thus the experiment isolates the influence of advective GAS-PHASE O2
transport while retaining hydraulic flow and dissolved-O2 transport.

## Drywell

- axisymmetric domain radius: 20 m
- drywell radius: 0.61 m
- drywell bottom elevation: 14 m
- injection head: 4.5 m
- production boundary: legacy whole-face midpoint
- initial O2 state: true-stagnant 2880-h O2 IC

## Completed diagnostic run

Result series:

drywell_head4p5m_phaseguardV2_faceSymDisp_gasO2AdvOFF_diag240h_6h

Simulation length:

240 h

The archived PVD and all referenced 6-h VTU outputs are stored in:

results/doubleOFF_240h/

## Scientific interpretation

This run was used as the gas-O2-advection-off control against the
formal 4.5-m DuMuX simulation with gas-phase O2 advection enabled.

It should not be described as a no-flow or no-advection hydraulic model.
