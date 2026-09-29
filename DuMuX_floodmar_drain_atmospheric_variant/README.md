# Flood MAR 60-day drainage atmospheric-boundary variant

Source: `DuMuX_floodmar_equalVolume_3cycle_20d_60d_GitHub_FINAL`.

## Change

Only `model/source/problem_floodpond5m_depth10cm.hh` is modified. The preexisting exposed top surface keeps its atmospheric H2O/N2/O2 gas-flux law. The original pond infiltration faces (bottom, r ≤ 5 m at z = 29.90 m; lower side at r = 5 m, z = 29.90–29.95 m) use that same law whenever the local prescribed pond head is ≤ 0. The gas law imposes no liquid-water supply; it permits gas advection and O2 diffusion, with water vapor carried by the gas. While the pond head is positive, the original Robin liquid-flux law, including the 0.5 cm dryout smoothing, remains in place. The upper pond side wall remains as in the source (no flux).

The exact original schedule remains: 5 cm wetting head; 20-day cycles; drawdown starts at 96.000000000000, 576.164044805955, and 1060.953978421305 h, each with a 5-hour ramp. Changing drainage gas exchange can alter later water states and injection volumes, so the original equal-volume calibration is not guaranteed to remain exact.

## Run on the Mac

Unzip the variant, then execute:

```bash
bash /path/to/floodmar_drain_atmospheric_variant/scripts/run_from_initial_mac.sh
```

The script copies **only the modified problem header** into the existing DuMuX source tree, saves the previous active header as `problem_before_atmospheric_variant.hh` in this variant folder, builds the existing restart-writer target, and runs the full 0–1440 h simulation from the original initial condition. Its PVD, VTUs, and log go under `results/` in this variant folder. It does not use the original intermediate restart snapshots because their states were computed with the old drainage boundary.

The script reads the mesh, oxygen initial condition, and parameters from the original `..._GitHub_FINAL/model/` directory. It applies the original final-run command-line overrides because the archived `params.input` still contains older defaults. The script does not run automatically when the ZIP is extracted.

## Verification points

- In the wet period, the pond faces still use the 5 cm Robin liquid boundary and original 5-hour ramp.
- After each local head reaches zero, the pond faces use the same atmospheric gas/O2 law as the exposed top surface, with no liquid supply.
- Compare cumulative pond inflow and gas/O2 exchange near 101, 581.164, and 1065.954 h, and at 480, 960, 1440 h.
- Surface liquid exfiltration is blocked in the dry atmospheric period by this zero-liquid-supply gas law. A seepage-face or explicit surface-water-storage model would be a separate experiment.
