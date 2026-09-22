# DuMuX Flood-MAR equal-volume 60-day final dataset

## Final spatial dataset
- Duration: **1440 h (60 d)**
- 3 cycles, **480 h (20 d)** each
- **241 VTU states**
- Nominal archive interval: **6 h**
- Actual first/last state: **0.000000000 / 1440.000000000 h**
- Maximum absolute actual-vs-nominal time offset: **8.500000 min**

## PVD files
**Canonical quantitative PVD:**  
`outputs_6h/floodmar_equalVolume_60d_final_6h_ACTUAL_TIME.pvd`

This uses each selected state's actual DuMuX simulation time. Use this for
scientific analysis and any time integration.

**Aligned visualization PVD:**  
`outputs_6h/floodmar_equalVolume_60d_final_6h_NOMINAL_6H.pvd`

This labels the same 241 VTUs at exact 0, 6, 12, ..., 1440 h. Use it for
side-by-side visualization against another nominal 6-h dataset.

The VTU contents are identical in the two PVDs; only the timestep labels differ.

## Locked same-volume schedule
| Cycle | Start (h) | Ramp start (h) | Ramp end (h) | Target (m3) | Calibrated (m3) |
|---|---:|---:|---:|---:|---:|
| C1 | 0 | 96.000000000000 | 101.000000000000 | 3956.009528 | 3957.230719 |
| C2 | 480 | 576.164044805955 | 581.164044805955 | 3956.009528 | 3956.298762 |
| C3 | 960 | 1060.953978421305 | 1065.953978421305 | 3956.009528 | 3954.068732 |

## Final model settings
- Flood-MAR ponding head: **5 cm**
- startup/fill ramp parameter: **1800 s**
- gas-phase O2 advection: **ON**
- longitudinal dispersivity: **0.50 m**
- transverse dispersivity: **0.05 m**
- gas mechanical dispersion: **OFF**
- MaxTimeStepSize: **300 s**
- VTK precision: **Float64**

## Structure
- `outputs_6h/`: canonical 241-state final trajectory
- `model/`: archived source, parameters, mesh, initial O2 profile
- `calibration/`: C2/C3 equal-volume calibration histories
- `restarts/`: retained restart files
- `scripts/`: continuation/packaging scripts
- `logs/`: retained completion logs
- `provenance/`: time mapping and locked schedule

`raw_mac_runs/` is intentionally not duplicated here because it largely
duplicates the final trajectory. It remains in the original working package:
`/Users/zhaozhe/Documents/GitHub/floof_drywell_compare/DuMuX_floodmar_equalVolume_3cycle_20d_60d_GitHub/raw_mac_runs`

## GitHub
`.gitattributes` is configured for Git LFS for VTU, DAT and MSH files.
