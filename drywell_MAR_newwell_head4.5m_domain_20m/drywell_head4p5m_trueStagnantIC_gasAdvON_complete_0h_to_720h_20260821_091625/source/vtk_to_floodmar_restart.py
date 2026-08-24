#!/usr/bin/env python3
"""Convert a serial, ASCII DuMuX cell-data VTU file to FloodMAR restart v1."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import xml.etree.ElementTree as ET


REQUIRED_ARRAYS = (
    "phase presence",
    "p_liq",
    "S_gas",
    "x^N2_liq",
    "x^O2_liq",
    "x^H2O_gas",
    "x^O2_gas",
)


def scalar_cell_arrays(vtu_path: Path) -> tuple[int, dict[str, list[float]]]:
    root = ET.parse(vtu_path).getroot()
    piece = root.find(".//Piece")
    if piece is None:
        raise ValueError("VTU file has no Piece element")

    num_cells = int(piece.attrib["NumberOfCells"])
    cell_data = piece.find("CellData")
    if cell_data is None:
        raise ValueError("VTU file has no CellData")

    arrays: dict[str, list[float]] = {}
    for data_array in cell_data.findall("DataArray"):
        name = data_array.attrib.get("Name")
        components = int(data_array.attrib.get("NumberOfComponents", "1"))
        if name and components == 1:
            arrays[name] = [float(value) for value in (data_array.text or "").split()]

    missing = [name for name in REQUIRED_ARRAYS if name not in arrays]
    if missing:
        raise ValueError(f"VTU file is missing arrays: {', '.join(missing)}")

    for name in REQUIRED_ARRAYS:
        if len(arrays[name]) != num_cells:
            raise ValueError(
                f"Array {name!r} has {len(arrays[name])} values, expected {num_cells}"
            )

    return num_cells, arrays


def write_restart(vtu_path: Path, output_path: Path, restart_time: float) -> None:
    num_cells, arrays = scalar_cell_arrays(vtu_path)

    with output_path.open("w", encoding="utf-8", newline="\n") as output:
        output.write(f"FLOODMAR_RESTART_V1 {num_cells} {restart_time:.17g}\n")

        for dof in range(num_cells):
            state_value = arrays["phase presence"][dof]
            state = int(round(state_value))
            if not math.isclose(state_value, state, rel_tol=0.0, abs_tol=1.0e-10):
                raise ValueError(f"Non-integral phase state {state_value} at DOF {dof}")

            pressure = arrays["p_liq"][dof]

            # DuMuX 2pnc, p0s1 formulation:
            # state 1 = liquid only: xN2_liq, xO2_liq
            # state 2 = gas only:    xH2O_gas, xO2_gas
            # state 3 = both phases: S_gas, xO2_liq
            if state == 1:
                switch = arrays["x^N2_liq"][dof]
                oxygen = arrays["x^O2_liq"][dof]
            elif state == 2:
                switch = arrays["x^H2O_gas"][dof]
                oxygen = arrays["x^O2_gas"][dof]
            elif state == 3:
                switch = arrays["S_gas"][dof]
                oxygen = arrays["x^O2_liq"][dof]
            else:
                raise ValueError(f"Unsupported phase state {state} at DOF {dof}")

            values = (pressure, switch, oxygen)
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"Non-finite primary variable at DOF {dof}: {values}")

            output.write(
                f"{dof} {state} {pressure:.17g} {switch:.17g} {oxygen:.17g}\n"
            )

    print(f"Wrote {num_cells} DOFs to {output_path}")
    print(f"Restart time: {restart_time:.17g} s ({restart_time/3600.0:.6f} h)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("vtu", type=Path, help="serial ASCII .vtu state file")
    parser.add_argument("output", type=Path, help="output restart .dat file")
    parser.add_argument(
        "--time",
        type=float,
        required=True,
        help="absolute simulation time represented by the VTU, in seconds",
    )
    args = parser.parse_args()

    write_restart(args.vtu, args.output, args.time)


if __name__ == "__main__":
    main()
