# Flood-MAR same-volume calibrated model archive

This archive was produced entirely from existing DuMuX output. **The model was not rerun.**

## Final equal-volume schedule

| Cycle | Ramp start (h) | Ramp end (h) | Volume (m3) | Error |
|---|---:|---:|---:|---:|
| C1 | 96.000000000000 | 101.000000000000 | 3957.230719 | +0.030869% |
| C2 | 576.164044805955 | 581.164044805955 | 3956.298762 | +0.007311% |
| C3 | 1060.953978421305 | 1065.953978421305 | 3954.068732 | -0.049059% |

Target volume: `3956.009528 m3/cycle`

## Output archive

Existing simulation states were sampled at approximately 6-hour intervals. The exact simulation time for every copied VTU is recorded in `OUTPUT_MANIFEST.csv`.

`outputs/floodmar_samevolume_existing_6h.pvd` references only the VTU files contained in the `outputs/` directory.

## Important limitation

Some original high-frequency VTU/PVD groups were deleted before this archive was created to recover disk space. Missing states were not reconstructed or interpolated. Therefore this archive contains only calibrated output files that still physically existed at packaging time.

## Included model files

- DuMuX source/header files
- `params.input`
- computational mesh
- 2880-h stagnant O2 initial condition
- CMake configuration
- C2/C3 calibration histories
- key restart files
- final schedule

## GitHub

VTU files are large binary scientific output. Git LFS is recommended.
