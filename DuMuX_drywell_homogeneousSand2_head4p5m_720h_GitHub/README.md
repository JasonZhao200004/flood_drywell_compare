# DuMuX Drywell MAR — Homogeneous Sand2

Homogeneous-texture sensitivity case for the
drywell MAR oxygen-transport model.

## Purpose

This case isolates the effect of subsurface
stratigraphy by replacing the original six-layer
profile with one homogeneous material.

The entire domain uses the hydraulic properties
of the original Layer 3 (`Sand2`).

The validated drywell geometry, operational
schedule, oxygen initial condition, oxygen
transport physics, and numerical formulation are
otherwise retained.

## Homogeneous material

Entire domain:

`Sand2`

Properties:

- intrinsic permeability = `1.925e-12 m2`
- porosity = `0.381022`

## Initial water condition

The whole domain starts at one uniform volumetric
water content:

`theta_w = 0.116926652 m3/m3`

Equivalent Sand2 pressure head:

`h = -29.770236345 cm`

The value was calculated as the mean initial water
content of the original Layer-3 / Sand2 interval.

## Initial oxygen condition

The O2 IC is intentionally NOT homogenized.

The case retains the same validated two-dimensional
TRUE-stagnant oxygen initial condition:

`oxygenIC_trueStagnantO2_2880h_2d_xO2.dat`

## Drywell operation

Injection head:

`450 cm = 4.5 m`

Schedule:

`4 d wet + 6 d drain`

Three cycles:

`720 h total`

## Transport physics

- gas-phase O2 advection: ON
- liquid mechanical dispersion: ON
- gas mechanical dispersion: OFF
- alphaL = 0.50 m
- alphaT = 0.05 m
- PhaseGuard V2
- face-symmetric CCTpfa dispersion

## Results

The full local simulation produced approximately
hourly output from 0 to 720 h.

For repository storage, this archive contains:

`0, 6, 12, ..., 720 h`

Total:

`121 VTUs`

The files are direct copies of the nearest real
DuMuX outputs.

No interpolation or modification of VTU fields was
performed.

Precise archived timestamps are recorded in:

`results/production_720h/TIME_MAPPING_6h.tsv`

Open the full archived time series with:

`DuMuX_homogeneousSand2_head4p5m_720h_every6h.pvd`

## Directory structure

### model/source

Final source/header implementation for this
homogeneous Sand2 target.

### model/config

- params.input
- floodmar_nonames.msh
- validated 2-D O2 IC
- homogeneous-case description
- exact runtime settings

### results/production_720h

- 121 six-hour VTUs
- PVD
- time mapping
- run summary

## Not included

The archive intentionally excludes:

- runtime log
- build files
- object files
- the complete approximately hourly 721-output series

The full-resolution run remains in the local DuMuX
build directory.
