# Flood MAR: 5 m pond radius, 5 cm head, 20 m domain

Independent DuMuX 2p3c axisymmetric Flood-MAR case prepared for the existing
WSL installation. It is based on the previously developed 20 m HYDRUS-style
pond boundary model and uses the numerical components validated by the 12 m
drywell run.

## Model definition

- radial domain: 0--20 m; elevation: 0--30 m;
- excavated pond radius: 5 m; excavation depth: 0.5 m;
- recharge ponding head: +5 cm;
- drainage head: -25 cm;
- schedule: three 4 d recharge / 6 d drainage cycles over 720 h;
- transition ramps: 96--101 h, 240--245 h, 336--341 h, 480--485 h, and
  576--581 h; first recharge starts with the same 12 h startup ramp used by
  the current drywell comparison;
- six-layer hydraulic properties, 2880 h oxygen initial profile, influent O2
  concentration, moisture-dependent O2 decay, and atmospheric gas exchange
  retained from the established models;
- liquid-phase mechanical dispersion: alpha_L=0.50 m and alpha_T=0.05 m;
- UMFPACK, chopped Newton updates, and v4 phase-switch hysteresis enabled.

The pond area is pi*(5 m)^2 = 78.54 m2. Here `5 m` means **radius**, not
diameter.

The mesh deliberately omits the Gmsh `PhysicalNames` block while retaining
integer physical tags. This avoids the UGGrid mesh-reader crash previously
observed under WSL.

## Install and build under WSL

Copy or unzip this directory somewhere inside WSL, then run:

```bash
cd /path/to/floodmar_pond5m_head5cm_domain20m
chmod +x install_and_build_wsl.sh
./install_and_build_wsl.sh
```

## Required validation sequence

First run 1 h:

```bash
RUN_DIR=~/dumux-work/dumux/dumux-floodmar/build-cmake/test/porousmediumflow/2p2c/floodmar_pond5m_head5cm_domain20m
cd "$RUN_DIR"

./floodmar_flood_pond5m_head5cm_robust params.input \
  -Problem.Name floodmar_pond5m_head5cm_test_1h \
  -TimeLoop.TEnd 3600 \
  -TimeLoop.OutputInterval 1800 \
  2>&1 | tee floodmar_pond5m_head5cm_test_1h.log
```

If that completes, run the full 720 h case in the background:

```bash
LOG=floodmar_pond5m_head5cm_4d6d_720h.log

nohup stdbuf -oL -eL \
./floodmar_flood_pond5m_head5cm_robust params.input \
  -Problem.Name floodmar_pond5m_head5cm_4d6d_720h \
  > "$LOG" 2>&1 &

echo $! | tee floodmar_pond5m_head5cm_4d6d_720h.pid
```

Monitor it at any time:

```bash
PID=$(cat floodmar_pond5m_head5cm_4d6d_720h.pid)
LOG=floodmar_pond5m_head5cm_4d6d_720h.log
ps -p "$PID" -o pid,etime,%cpu,%mem,state,command
grep 'Time step .*done' "$LOG" | tail -n 1
echo -n 'Newton failures: '
grep -c 'Newton solver did not converge' "$LOG"
```
