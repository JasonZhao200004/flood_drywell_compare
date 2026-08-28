

## Saved spatial outputs

To keep the repository compact, the complete 1 h VTU series is not
included.

Spatial model outputs are sampled approximately every 6 h from the
original 240 h simulation:

\[
t = 0,\ 6,\ 12,\ldots,240\ \mathrm{h}.
\]

This gives 41 retained VTU snapshots.

Because DuMuX adaptive time stepping may place individual PVD outputs
slightly away from the nominal output time, each retained VTU is the
available simulation output closest to the corresponding 6 h target.

See:

`outputs_6h/output_6h_time_mapping.csv`

for the exact simulated time associated with each snapshot.
