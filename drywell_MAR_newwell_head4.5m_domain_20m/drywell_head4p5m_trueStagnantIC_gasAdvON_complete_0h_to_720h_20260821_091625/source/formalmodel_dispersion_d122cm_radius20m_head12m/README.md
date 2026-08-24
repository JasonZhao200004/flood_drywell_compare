# Drywell: radius-20-m domain, 12-m head, dispersion and advection

This is a complete axisymmetric DuMuX project derived from the corrected
1.22-m-diameter drywell model.

## Model definition

- Radial domain: `0 <= r <= 20 m`; elevation: `0 <= z <= 30 m`.
- Drywell radius: `0.61 m`; bottom elevation: `14 m`.
- VariHead boundary: corrected lower drywell wall and bottom only.
- Injection head: `1200 cm` above the drywell bottom during injection periods.
- Drainage head: `-25 cm`, unchanged.
- Schedule: three cycles over 720 h, unchanged.
- Liquid dispersion: longitudinal `0.50 m`, transverse `0.05 m`.
- Dispersion is not applied to the gas phase, matching the previous model.
- Gas/liquid Darcy advection, molecular diffusion, gravity, moisture-dependent
  first-order O2 decay, and the 2880-h initial O2 profile remain active.
- Six soil layers and their hydraulic parameters are unchanged.

The corrected 12-m mesh is included as `floodmar_12m_source.msh`. The supplied
generator leaves the drywell and refined near-field mesh (`r <= 3 m`) exactly
unchanged and expands only the far field to 20 m. It preserves all physical
groups and verifies that no triangle becomes degenerate.

## Install and build

From the extracted project directory:

```bash
chmod +x install_and_build.sh
./install_and_build.sh
```

The installer copies the project into the DuMuX module when needed, registers
its CMake subdirectory, generates the 20-m mesh, compiles the executable, and
creates the required input-file links in the build run directory.

## Recommended one-hour test

```bash
cd ~/dumux-work/dumux/build-serial/dumux-floodmar/test/porousmediumflow/2p2c/formalmodel_dispersion_d122cm_radius20m_head12m

caffeinate -i ./floodmar_formal_dispersion_d122cm_radius20m_head12m params.input \
  -Problem.Name drywell_radius20m_head12m_test_1h \
  -TimeLoop.TEnd 3600 \
  -TimeLoop.DtInitial 1.08 \
  -TimeLoop.MaxTimeStepSize 60 \
  -TimeLoop.OutputInterval 1800 \
  -Newton.MaxRelativeShift 1e-7 \
  -Newton.MaxSteps 18 \
  2>&1 | tee drywell_radius20m_head12m_test_1h.log
```

## Full 720-hour run

```bash
cd ~/dumux-work/dumux/build-serial/dumux-floodmar/test/porousmediumflow/2p2c/formalmodel_dispersion_d122cm_radius20m_head12m

caffeinate -i ./floodmar_formal_dispersion_d122cm_radius20m_head12m params.input \
  -Problem.Name drywell_radius20m_head12m_dispersion_720h \
  2>&1 | tee drywell_radius20m_head12m_dispersion_720h.log
```

`caffeinate -i` prevents macOS idle sleep while the simulation is running.
Output is written every 6 h. The injection head can later be overridden without
recompiling, for example `-Problem.DrywellInjectionHead 900` for 9 m.

The startup ramp now raises the effective drywell water level from the bottom
to 12 m over 12 hours. It does not activate the complete 12 m wall at the
first time step. This reduces simultaneous phase switching along the well.

Because 12 m is a strong head, run the one-hour test first and inspect the time
step history before committing to 720 h.
