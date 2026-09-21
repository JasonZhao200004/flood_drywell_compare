# DuMuX 12-m drywell equal-volume 60-day run

Target gross injection volume:

**3956.009528 m3 per cycle**

Final schedule:

- C1: 0-480 h; wet end = 163.125247592924 h
- C2: 480-960 h; wet end = 737.084098745400 h
- C3: 960-1440 h; wet end = 1233.636146520033 h

Model configuration:

- Drywell injection head = 1200 cm (12 m)
- Initial startup ramp = 1800 s
- Cycle 2/3 refill transition = 0.5 h
- Drainage head after target volume = -25 cm
- Gas-phase O2 advection = ON
- alpha_L = 0.50 m
- alpha_T = 0.05 m
- Gas mechanical dispersion = OFF
- PhaseGuard V2
- Face-symmetric CCTpfa dispersion
- true-stagnant 2880-h O2 initial condition
- MaxTimeStepSize = 360 s
- Final archive output ~= 6 h

The three formal model segments are continuous through exact
primary-variable restart files at 480 h and 960 h.

Combined PVD:

outputs_6h/drywell_head12m_equalVolume_3cycle_20d_60d_6h.pvd

Calibration trajectories used 1-h VTU output. Gross drywell
injection was reconstructed from the model's Robin boundary
states on a 0.05-h integration grid. Calibration VTUs were
deleted after the cumulative-volume CSV files were saved.
