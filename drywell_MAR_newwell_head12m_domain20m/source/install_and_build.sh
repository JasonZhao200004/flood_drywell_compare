#!/usr/bin/env bash
set -euo pipefail

PROJECT_NAME="formalmodel_dispersion_d122cm_radius20m_head12m"
TARGET_NAME="floodmar_formal_dispersion_d122cm_radius20m_head12m"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DUMUX_WORK="${DUMUX_WORK:-$HOME/dumux-work/dumux}"
MODULE_ROOT="$DUMUX_WORK/dumux-floodmar"
SOURCE_PARENT="$MODULE_ROOT/test/porousmediumflow/2p2c"
TARGET_DIR="$SOURCE_PARENT/$PROJECT_NAME"
BUILD_ROOT="$DUMUX_WORK/build-serial/dumux-floodmar"
RUN_DIR="$BUILD_ROOT/test/porousmediumflow/2p2c/$PROJECT_NAME"

if [[ "$SCRIPT_DIR" != "$TARGET_DIR" ]]; then
    mkdir -p "$TARGET_DIR"
    rsync -a --exclude '*.log' --exclude '*.pvd' --exclude '*.vtu' \
        "$SCRIPT_DIR/" "$TARGET_DIR/"
fi

python3 "$TARGET_DIR/generate_radius20m_mesh.py"

PARENT_CMAKE="$SOURCE_PARENT/CMakeLists.txt"
ENTRY="add_subdirectory($PROJECT_NAME)"
if ! grep -Fqx "$ENTRY" "$PARENT_CMAKE"; then
    printf '\n%s\n' "$ENTRY" >> "$PARENT_CMAKE"
fi

cmake -S "$MODULE_ROOT" -B "$BUILD_ROOT"
cmake --build "$BUILD_ROOT" --target "$TARGET_NAME" -j 4

mkdir -p "$RUN_DIR"
ln -sf "$TARGET_DIR/params.input" "$RUN_DIR/params.input"
ln -sf "$TARGET_DIR/floodmar.msh" "$RUN_DIR/floodmar.msh"
ln -sf "$TARGET_DIR/oxygen_initial_2880.dat" "$RUN_DIR/oxygen_initial_2880.dat"

echo
echo "Built $TARGET_NAME"
echo "Run directory:"
echo "$RUN_DIR"
