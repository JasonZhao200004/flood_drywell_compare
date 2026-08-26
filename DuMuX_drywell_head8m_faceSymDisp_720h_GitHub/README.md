# DuMuX Drywell MAR — 8 m Head

Final 720 h DuMuX drywell-MAR sensitivity simulation
with an 8 m injection head.

## Geometry

Axisymmetric domain:

- radial extent: r = 0–20 m
- depth: z = 0–30 m
- ground surface: z = 30 m
- drywell radius: 0.61 m
- drywell diameter: 1.22 m
- drywell bottom: z = 14 m

The numerical radial extent is 20 m. When revolved
around the symmetry axis, the modeled soil domain
therefore has a physical radius of 20 m.

## 8 m hydraulic condition

DrywellInjectionHead = 800 cm.

The head is referenced to the drywell bottom.

Therefore:

- drywell bottom elevation = z 14 m
- maximum water-surface elevation = z 22 m
- ground surface = z 30 m
- water level remains approximately 8 m below ground

## Final transport formulation

- TRUE-stagnant pre-equilibrated O2 initial condition
- gas O2 advection: ON
- PhaseGuard V2: ON
- face-symmetric CCTpfa dispersion: ON
- liquid mechanical dispersion: ON
- gas mechanical dispersion: OFF
- longitudinal dispersivity: 0.50 m
- transverse dispersivity: 0.05 m

## Simulation

- total duration: 720 h
- three 240 h cycles
- original output interval: approximately 1 h
- original outputs: 721 VTUs

Final solver check:

- Newton failures: 9
- Linear failures: 0
- Fatal-like messages: 0
- 720 h completed successfully

## Archived results

To match the temporal resolution of the HYDRUS
reference dataset, the GitHub archive contains:

**0–720 h every 6 h = 121 VTU files**

Filename convention:

- `dumux_head8m_000h.vtu`
- `dumux_head8m_006h.vtu`
- ...
- `dumux_head8m_720h.vtu`

The files are the nearest real DuMuX outputs to each
6 h target time.

No field interpolation was performed.

The exact original simulation times are listed in:

`TIME_MAPPING_6h.tsv`

and are preserved in:

`DuMuX_head8m_720h_every6h.pvd`

## Directory structure

### model/config

Contains:

- params.input
- floodmar_nonames.msh
- TRUE-stagnant O2 initial condition
- exact 8 m production runtime settings

### model/source

Contains the final custom DuMuX source/header files,
including the PhaseGuard V2 and face-symmetric
dispersion implementation.

### results/production_720h

Contains:

- 121 six-hour VTUs
- PVD
- time mapping table
- run summary

## Intentionally excluded

- simulation log files
- CMake build products
- object files
- temporary files
- complete 721-file hourly VTU sequence

The original hourly run remains in the local DuMuX
build directory.
