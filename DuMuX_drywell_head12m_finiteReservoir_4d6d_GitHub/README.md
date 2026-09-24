# DuMuX drywell MAR: 12 m head with finite-reservoir drainage

Archive date: 2026-09-24

## Purpose

This repository preserves a finite-reservoir drywell boundary trial derived
directly from the previously successful 720 h production drywell model.

The old prescribed dry-period head (`-25 cm`) is replaced by the finite volume
of water remaining inside the drywell after external supply is stopped.

## Boundary sequence

- `0-96 h`: original successful production wetting behavior.
- `96 h`: external supply OFF.
- `96-240 h`: finite 0-D reservoir with `V = pi r^2 H`.
- No artificial 96-101 h shutdown ramp.
- No prescribed `-25 cm` drainage head.

Reservoir radius is `0.61 m`; initial water depth is `12 m`, corresponding to
`14.0278395 m3`.

The reservoir state is updated once after each accepted DuMuX timestep:

`V_new = max(0, V_old + Q_well dt)`

where `Q_well < 0` is reservoir-to-soil flow and `Q_well > 0` is soil-to-well flow.

## Retained production physics

- gas-phase O2 advection: ON
- liquid mechanical dispersion: ON
- gas mechanical dispersion: OFF
- alpha_L = `0.50 m`
- alpha_T = `0.05 m`
- PhaseGuard V2
- face-symmetric CCTpfa dispersion
- true-stagnant 2880 h O2 initial condition
- six soil layers
- legacy whole-face midpoint drywell hydraulic boundary

## 240 h diagnostic

At `96 h`:

- head = `12.000 m`
- storage = `14.0278395 m3`

At `240 h`:

- head = `4.930447 m`
- storage = `5.7636260 m3`
- net reservoir flux = `-2.267160e-6 m3/s`

The reservoir does not empty by 240 h.

The first positive soil-to-well flux occurs at approximately `170.448 h`.
There are 203 positive-flux rows in the stored history.

Discrete reservoir water-balance residual:

- maximum absolute = `2.0358e-13 m3`
- mean absolute = `5.8238e-15 m3`

## Numerical limitation

This diagnostic uses `TimeLoop.MaxTimeStepSize = 360 s`.

The post-96 h falling-head trajectory is stable and mass-balanced, but timestep
sensitivity has not yet been completed. A recommended next check is post-shutoff
maximum timesteps of `360 s`, `60 s`, and `30 s`.

## Repository contents

- `model/source/`: baseline + finite-reservoir source
- `model/config/`: mesh, O2 IC, base params, run records
- `results/diagnostic_240h/`: reservoir history CSV and complete log
- `results/key_vtk/`: selected snapshots if available during packaging
- `scripts/analyze_reservoir_history.py`: stdlib-only history analysis
- `run_240h_example.sh`: example reproduction command

Full hourly VTUs are intentionally not archived to keep the repository compact.

## Analyze the archived reservoir history

```bash
python3 scripts/analyze_reservoir_history.py
```

No pandas installation is required.

## Production reference

`model/config/RUN_12m_720h.txt` preserves the runtime overrides from the
successful 720 h production model. Those runtime overrides, not `params.input`
alone, define the production numerical baseline.
