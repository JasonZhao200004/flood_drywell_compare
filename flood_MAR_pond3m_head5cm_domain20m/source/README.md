# Flood-MAR pond model

Independent axisymmetric DuMuX project derived from the validated CCTpfa
water–N2–O2 model.

## Geometry and boundaries

- Soil domain: radius 12 m, elevation 0–30 m.
- Pond: radius 3 m and excavation depth 0.5 m.
- VariHead: pond bottom (`y=29.5 m`, `r=0–3 m`) plus only the lowest
  5 cm of the vertical rim (`r=3 m`, `y=29.5–29.55 m`).
- Recharge head: +5 cm; drainage head: -25 cm.
- Recharge schedule: 0–96 h, 240–336 h, and 480–576 h, with the same
  5 h transitions used in the formal drywell model.
- Ground outside the pond is open to atmospheric gas; the outer radial
  side and the dry part of the pond rim are no-flow; the bottom is free
  drainage.
- Initial water state uses the validated piecewise pressure-head profile.
- Initial dissolved O2 uses the x=3 m profile at 2881.22 h (the closest
  stored output to 2880 h).

## Install into the DuMuX module

Copy this entire directory to:

```text
~/dumux-work/dumux/dumux-floodmar/test/porousmediumflow/2p2c/floodmodel
```

Add this line once to the parent `2p2c/CMakeLists.txt`:

```cmake
add_subdirectory(floodmodel)
```

Configure and build:

```bash
cmake -S ~/dumux-work/dumux/dumux-floodmar \
      -B ~/dumux-work/dumux/build-serial/dumux-floodmar

cmake --build ~/dumux-work/dumux/build-serial/dumux-floodmar \
  --target floodmar_flood -j 4
```

Run the supplied 1 h check first:

```bash
cd ~/dumux-work/dumux/build-serial/dumux-floodmar/test/porousmediumflow/2p2c/floodmodel

caffeinate -i ./floodmar_flood params.input \
  -Problem.Name floodmar_pond_test_1h \
  -TimeLoop.TEnd 3600 \
  -TimeLoop.OutputInterval 1800 \
  2>&1 | tee floodmar_pond_test_1h.log
```

Then run the complete 720 h model:

```bash
caffeinate -i ./floodmar_flood params.input \
  2>&1 | tee floodmar_pond3m_head5cm_720h.log
```

The mesh and oxygen profile are already included. To regenerate them:

```bash
python3 generate_pond_mesh.py
python3 extract_oxygen_profile.py
```

