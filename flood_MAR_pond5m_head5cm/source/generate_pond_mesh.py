#!/usr/bin/env python3
"""Generate the axisymmetric 20 m-radius x 30 m Flood-MAR pond mesh.

The soil domain has a 5 m-radius, 0.5 m-deep excavation. Physical line tag 2
marks the pond bottom and complete vertical rim. The problem code activates
only the submerged part of the rim, so its wetted height follows the positive
pond head. The output is Gmsh 2.2 ASCII without a PhysicalNames block because
the WSL/Linux UGGrid reader previously crashed on otherwise equivalent named
meshes. Integer physical tags are retained on every element.
"""

from pathlib import Path


OUT = Path(__file__).with_name("floodmar_pond5m_nonames.msh")
POND_RADIUS = 5.0
DOMAIN_RADIUS = 20.0


def segmented_grid(breaks, target_step):
    values = [float(breaks[0])]
    for a, b in zip(breaks[:-1], breaks[1:]):
        n = max(1, round((b - a) / target_step))
        values.extend(a + (b - a) * i / n for i in range(1, n + 1))
    return values


# Include the pond rim and every soil interface exactly.
x_all = segmented_grid([0.0, POND_RADIUS], 0.10)
x_all += segmented_grid([POND_RADIUS, DOMAIN_RADIUS], 0.25)[1:]
y_base = segmented_grid(
    [0.0, 4.0, 8.5, 12.5, 18.5, 22.0, 29.5], 0.20
)
# Local grading resolves 2.5 and 5 cm head tests without a hanging node.
y_upper = [29.5, 29.525, 29.55, 29.60, 29.70, 29.80, 29.90, 30.0]
x_upper = [x for x in x_all if x >= POND_RADIUS - 1.0e-12]

nodes = []
node_id = {}


def add_node(x, y):
    key = (round(x, 10), round(y, 10))
    if key not in node_id:
        node_id[key] = len(nodes) + 1
        nodes.append((float(x), float(y), 0.0))
    return node_id[key]


for y in y_base:
    for x in x_all:
        add_node(x, y)
for y in y_upper:
    for x in x_upper:
        add_node(x, y)


elements = []


def add_line(physical, a, b):
    elements.append((1, physical, [a, b]))


def material_tag(y_center):
    if y_center >= 22.0:
        return 101
    if y_center >= 18.5:
        return 102
    if y_center >= 12.5:
        return 103
    if y_center >= 8.5:
        return 104
    if y_center >= 4.0:
        return 105
    return 106


def add_quad(x0, x1, y0, y1, flip):
    a = add_node(x0, y0)
    b = add_node(x1, y0)
    c = add_node(x1, y1)
    d = add_node(x0, y1)
    tag = material_tag(0.5 * (y0 + y1))
    if flip:
        elements.append((2, tag, [a, b, d]))
        elements.append((2, tag, [b, c, d]))
    else:
        elements.append((2, tag, [a, b, c]))
        elements.append((2, tag, [a, c, d]))


# Boundary line elements.
for x0, x1 in zip(x_upper[:-1], x_upper[1:]):
    add_line(1, add_node(x0, 30.0), add_node(x1, 30.0))  # atmosphere
for x0, x1 in zip(x_all[:-1], x_all[1:]):
    if x1 <= POND_RADIUS + 1.0e-12:
        add_line(2, add_node(x0, 29.5), add_node(x1, 29.5))  # pond bottom
for y0, y1 in zip(y_upper[:-1], y_upper[1:]):
    add_line(2, add_node(POND_RADIUS, y0), add_node(POND_RADIUS, y1))
for y0, y1 in zip(y_base[:-1], y_base[1:]):
    add_line(3, add_node(0.0, y0), add_node(0.0, y1))  # axis
for x0, x1 in zip(x_all[:-1], x_all[1:]):
    add_line(4, add_node(x0, 0.0), add_node(x1, 0.0))  # free drainage
for y0, y1 in zip((y_base + y_upper[1:])[:-1], (y_base + y_upper[1:])[1:]):
    add_line(5, add_node(DOMAIN_RADIUS, y0), add_node(DOMAIN_RADIUS, y1))  # outer no-flow

# Soil triangles below the pond bottom and in the outside upper shoulder.
for j, (y0, y1) in enumerate(zip(y_base[:-1], y_base[1:])):
    for i, (x0, x1) in enumerate(zip(x_all[:-1], x_all[1:])):
        add_quad(x0, x1, y0, y1, (i + j) % 2 == 1)
for j, (y0, y1) in enumerate(zip(y_upper[:-1], y_upper[1:])):
    for i, (x0, x1) in enumerate(zip(x_upper[:-1], x_upper[1:])):
        add_quad(x0, x1, y0, y1, (i + j) % 2 == 1)

with OUT.open("w", encoding="ascii") as f:
    f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
    f.write(f"$Nodes\n{len(nodes)}\n")
    for i, (x, y, z) in enumerate(nodes, 1):
        f.write(f"{i} {x:.10g} {y:.10g} {z:.10g}\n")
    f.write("$EndNodes\n")
    f.write(f"$Elements\n{len(elements)}\n")
    for i, (kind, physical, conn) in enumerate(elements, 1):
        joined = " ".join(map(str, conn))
        f.write(f"{i} {kind} 2 {physical} {physical} {joined}\n")
    f.write("$EndElements\n")

print(f"Wrote {OUT}")
print(f"nodes={len(nodes)}, elements={len(elements)}")
