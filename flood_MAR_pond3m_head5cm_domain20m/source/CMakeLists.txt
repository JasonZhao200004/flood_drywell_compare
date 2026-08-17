# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links(FILES
    params.input
    floodmar_pond.msh
    oxygen_initial_2880.dat
)

dumux_add_test(
    NAME floodmar_flood
    LABELS porousmediumflow 2pnc floodmar
    SOURCES main.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMar
    COMMAND floodmar_flood
    CMD_ARGS params.input
)

target_link_libraries(
    floodmar_flood
    PRIVATE
    Dune::UGGrid
)

