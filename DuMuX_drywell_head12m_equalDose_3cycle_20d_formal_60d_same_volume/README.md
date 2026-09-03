# Drywell MAR – Final 60-Day Equal-Volume Three-Cycle Model

## Final experimental design

Three consecutive 20-day Drywell MAR cycles.

Target gross injection volume:

**3769.844 m3 per cycle**

### Cycle 1

- Wetting: 0 – 176.490541 h
- Wetting duration: 7.353773 d
- Drainage: 176.490541 – 480 h
- Drainage duration: 12.646227 d

### Cycle 2

- Wetting: 480 – 725.850601953 h
- Wetting duration: 10.243775 d
- Drainage: 725.850601953 – 960 h
- Drainage duration: 9.756225 d

### Cycle 3

- Wetting: 960 – 1210.359144226 h
- Wetting duration: 10.431631 d
- Drainage: 1210.359144226 – 1440 h
- Drainage duration: 9.568369 d

Total:

**1440 h = 60 days**

## Model setup

- Drywell injection head: 1200 cm
- Drainage head: -25 cm
- Startup ramp: 1800 s at global startup only
- No additional ramp at C2/C3 switching
- Gas-phase O2 advection: ON
- Liquid longitudinal dispersivity: 0.50 m
- Liquid transverse dispersivity: 0.05 m
- Gas mechanical dispersion: OFF
- Drywell boundary velocity excluded from mechanical dispersion
- Initial oxygen field: true-stagnant 2880 h pre-equilibrated state

## 6-hour visualization trajectory

The original formal run contains 1441 approximately hourly VTU outputs.

For repository storage and routine ParaView visualization, the trajectory
has been subsampled to one snapshot every 6 h:

**241 target time points: 0, 6, 12, ..., 1440 h**

Because some original DuMuX outputs near operational transitions are
slightly offset from integer hours, each 6-hour target uses the nearest
available original VTU.

The PVD stores the **actual simulation time** of each selected VTU.

Open:

`VTU_6h/drywell_head12m_equalDose_3cycle_20d_formal_60d_6h.pvd`

in ParaView.

`VTU_6h/VTU_6h_manifest.csv` records the target time, actual time, and
time offset of every selected output.

## Keyframes

The wetting cutoffs are not exact multiples of 6 h.

Therefore, separate visualization keyframes are retained for:

- t = 0
- C1 wet end
- C1 cycle end / C2 start
- C2 wet end
- C2 cycle end / C3 start
- C3 wet end
- C3 cycle end

Open:

`keyframes/final60d_keyframes.pvd`

for these event-focused snapshots.

For exact quantitative values at the non-output wetting cutoff times,
use the interpolation procedures in the analysis scripts rather than
the visualization snapshots.

## Oxygen recovery

Main-plume positive residual at the end of each 20-day cycle:

- C1: 0.000632%
- C2: 45.220106%
- C3: 44.440909%

C2 and C3 residual trajectories are non-monotonic during drainage.
The positive excess initially decreases substantially and later rises
again, so the cycle-end positive residual should not automatically be
interpreted as persistence of the original injection-derived oxygen
plume.

See:

`results/final60d_O2_recovery.txt`

for the complete analysis.

## Directory structure

- `model_source/` – final DuMuX model source
- `inputs/` – mesh, params, and oxygen initial condition
- `VTU_6h/` – 241 subsampled simulation outputs
- `keyframes/` – selected operational transition snapshots
- `analysis/` – postprocessing scripts
- `results/` – recovery and calibration results

The original 1441-output trajectory remains in the WSL run directory
and is not duplicated in this repository package.
