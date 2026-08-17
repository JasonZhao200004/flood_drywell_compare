# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links(FILES
    params.input
    floodmar_pond5m_nonames.msh
    oxygen_initial_2880.dat
)

dumux_add_test(
    NAME floodmar_flood_pond5m_head5cm_robust
    LABELS porousmediumflow 2pnc floodmar dispersion umfpack chopped
    SOURCES main.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMarPond5m
    COMMAND floodmar_flood_pond5m_head5cm_robust
    CMD_ARGS params.input
)

target_link_libraries(
    floodmar_flood_pond5m_head5cm_robust
    PRIVATE
    Dune::UGGrid
)
