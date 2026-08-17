#!/usr/bin/env bash
set -euo pipefail

project_name="floodmar_pond5m_head5cm_domain20m"
source_dir="$(cd "$(dirname "$0")" && pwd)"
dumux_root="${DUMUX_ROOT:-$HOME/dumux-work/dumux}"
module_root="$dumux_root/dumux-floodmar"
target_dir="$module_root/test/porousmediumflow/2p2c/$project_name"
parent_cmake="$module_root/test/porousmediumflow/2p2c/CMakeLists.txt"
build_root="$module_root/build-cmake"
run_dir="$build_root/test/porousmediumflow/2p2c/$project_name"
target="floodmar_flood_pond5m_head5cm_robust"

if [[ "$source_dir" != "$target_dir" ]]; then
    mkdir -p "$target_dir"
    cp -R "$source_dir"/. "$target_dir"/
fi

python3 "$target_dir/generate_pond_mesh.py"

if ! grep -Fxq "add_subdirectory($project_name)" "$parent_cmake"; then
    printf '\nadd_subdirectory(%s)\n' "$project_name" >> "$parent_cmake"
fi

cmake -S "$module_root" -B "$build_root"
cmake --build "$build_root" --target "$target" -j 8

mkdir -p "$run_dir"
for input in params.input floodmar_pond5m_nonames.msh oxygen_initial_2880.dat; do
    ln -sfn "$target_dir/$input" "$run_dir/$input"
done

printf '\nBuilt successfully.\nRun directory: %s\nTarget: %s\n' \
       "$run_dir" "$target"
