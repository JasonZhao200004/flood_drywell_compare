# DuMuX Drywell MAR — 12 m Head

Final 720 h high-head drywell MAR sensitivity run.

## Geometry

2D axisymmetric model.

- radial domain: 0–20 m
- depth: 30 m
- ground surface: z = 30 m
- drywell radius: 0.61 m
- drywell diameter: 1.22 m
- drywell bottom: z = 14 m

The numerical horizontal coordinate represents
radial distance, so the physical axisymmetric
domain radius is 20 m.

## Hydraulic condition

DrywellInjectionHead:

`1200 cm = 12 m`

Head is referenced to the drywell bottom.

Therefore:

- drywell bottom = z 14 m
- maximum water elevation = z 26 m
- ground surface = z 30 m
- maximum water level is approximately 4 m below ground

## Initial condition

The production run was configured with:

`oxygenIC_trueStagnantO2_2880h_2d_xO2.dat`

This is the same TRUE-stagnant pre-equilibrated
oxygen IC used for the validated lower-head cases.

## Transport formulation

- gas O2 advection: ON
- liquid mechanical dispersion: ON
- gas mechanical dispersion: OFF
- alphaL = 0.50 m
- alphaT = 0.05 m
- PhaseGuard V2
- face-symmetric CCTpfa dispersion

## Production run

Duration:

`720 h`

Three 240 h cycles.

Original output:

`721 VTUs`

Final run check:

- 720 h reached successfully
- Newton failures = 12
- Linear failures = 0
- Fatal/abort/exception = 0

## Archived results

For repository storage, the original approximately
hourly outputs are reduced to the same temporal
resolution used for the HYDRUS reference:

`0–720 h every 6 h`

Total:

`121 VTU files`

Naming:

- dumux_head12m_000h.vtu
- dumux_head12m_006h.vtu
- ...
- dumux_head12m_720h.vtu

These files are copies of the nearest actual DuMuX
outputs. No field interpolation is applied.

Precise times are documented in:

`TIME_MAPPING_6h.tsv`

The PVD is:

`DuMuX_head12m_720h_every6h.pvd`

## Directory structure

### model/config

Contains:

- params.input
- floodmar_nonames.msh
- TRUE-stagnant oxygen IC
- exact 12 m runtime settings

### model/source

Contains the final source/header implementation.

### results/production_720h

Contains:

- 121 six-hour VTUs
- PVD
- time mapping
- run summary

## Intentionally excluded

- .log files
- build-system temporary files
- object files
- complete 721-file hourly output sequence

The complete hourly result remains in the local
DuMuX build directory.
