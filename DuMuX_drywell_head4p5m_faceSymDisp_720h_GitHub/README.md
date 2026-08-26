# DuMuX Drywell MAR O2 Model — 4.5 m Head

Final production and validation archive.

## Production model

- DuMuX CCTpfa
- 2D axisymmetric domain
- Radius: 20 m
- Depth: 30 m
- Drywell radius: 0.61 m
- Drywell bottom: z = 14 m
- Injection head: 450 cm
- Duration: 720 h
- Three 240 h cycles

## Final numerical formulation

- PhaseGuard V2
- Face-symmetric CCTpfa mechanical dispersion
- alphaL = 0.50 m
- alphaT = 0.05 m
- Liquid mechanical dispersion: ON
- Gas mechanical dispersion: OFF
- Gas O2 advection: ON
- TRUE-stagnant pre-equilibrated O2 initial condition

## Archived production output

The original production simulation wrote approximately hourly output.

For repository storage, the production result has been reduced to the
same temporal interval as the HYDRUS reference:

**0–720 h every 6 h = 121 VTU files**

Files are named:

- `dumux_000h.vtu`
- `dumux_006h.vtu`
- ...
- `dumux_720h.vtu`

The archived files are the nearest original DuMuX outputs to each 6 h
target. No VTU field interpolation was performed.

`TIME_MAPPING_6h.tsv` records:

- requested 6 h target time
- actual DuMuX output time
- time offset in minutes
- original VTU filename

`DuMuX_720h_every6h.pvd` contains the actual simulation timestamps.

## HYDRUS comparison

HYDRUS reference:

- 121 fields
- 0–720 h
- 6 h interval
- 13,480 nodes
- 26,347 triangles

Primary subsurface comparison region:

`r <= 6 m, z < 24 m`

Cycle-average RMSE/Cin:

- Cycle 1: 0.0573
- Cycle 2: 0.0523
- Cycle 3: 0.0485

At 720 h:

- RMSE/Cin = 0.0205
- bias/Cin = +0.0127
- Pearson r = 0.8893

## Numerical validation

- 720 h completed
- 721 original production outputs
- drywell-bottom maximum C/Cin = 0.9822
- phase=3 with Sgas <= 1e-8: 0 occurrences

Repeatable wet-to-drain transition maxima:

- Cycle 1: 1.1511 Cin at ~100.04 h
- Cycle 2: 1.1942 Cin at ~339.06 h
- Cycle 3: 1.1956 Cin at ~579.01 h

These small repeatable transition responses are distinct from the former
5–7 Cin numerical hotspot.

## Not included

- log files
- build files
- temporary files
- the complete set of 721 hourly VTUs

Only the 121 six-hour production outputs are archived.
