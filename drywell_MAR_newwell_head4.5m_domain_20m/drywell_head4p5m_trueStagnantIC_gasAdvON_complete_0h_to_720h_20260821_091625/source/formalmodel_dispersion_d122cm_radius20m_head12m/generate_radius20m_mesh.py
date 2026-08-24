#!/usr/bin/env python3
"""Expand the corrected 12 m drywell mesh to a 20 m radial domain.

The drywell and the refined near field (r <= 3 m) are left unchanged.
Only the far field is stretched monotonically from r=3..12 m to r=3..20 m.
Physical boundary/material tags and element connectivity are preserved.
"""

from pathlib import Path

SOURCE = Path(__file__).with_name("floodmar_12m_source.msh")
TARGET = Path(__file__).with_name("floodmar.msh")
KEEP_RADIUS = 3.0
OLD_RADIUS = 12.0
NEW_RADIUS = 20.0


def mapped_radius(x: float) -> float:
    if x <= KEEP_RADIUS:
        return x
    return KEEP_RADIUS + (x - KEEP_RADIUS) * (
        (NEW_RADIUS - KEEP_RADIUS) / (OLD_RADIUS - KEEP_RADIUS)
    )


def main() -> None:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    start = lines.index("$Nodes")
    end = lines.index("$EndNodes")
    count = int(lines[start + 1])
    if end - start - 2 != count:
        raise RuntimeError("Unexpected Gmsh 2.2 node section")

    nodes = {}
    for i in range(start + 2, end):
        fields = lines[i].split()
        node_id = int(fields[0])
        x, y, z = map(float, fields[1:4])
        x_new = mapped_radius(x)
        nodes[node_id] = (x_new, y)
        lines[i] = f"{node_id} {x_new:.12g} {y:.12g} {z:.12g}"

    element_start = lines.index("$Elements")
    element_end = lines.index("$EndElements")
    minimum_area = float("inf")
    triangles = 0
    for line in lines[element_start + 2:element_end]:
        fields = line.split()
        if int(fields[1]) != 2:
            continue
        tag_count = int(fields[2])
        node_ids = list(map(int, fields[3 + tag_count:6 + tag_count]))
        (x1, y1), (x2, y2), (x3, y3) = (nodes[n] for n in node_ids)
        area = abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))/2.0
        minimum_area = min(minimum_area, area)
        triangles += 1

    if triangles == 0 or minimum_area <= 0.0:
        raise RuntimeError("Invalid or degenerate triangle found")
    if abs(max(x for x, _ in nodes.values()) - NEW_RADIUS) > 1e-10:
        raise RuntimeError("Outer radius was not mapped to 20 m")

    TARGET.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {TARGET}")
    print(f"nodes={count}, triangles={triangles}")
    print(f"radial range=0..{NEW_RADIUS} m; minimum triangle area={minimum_area:.6e} m2")


if __name__ == "__main__":
    main()
