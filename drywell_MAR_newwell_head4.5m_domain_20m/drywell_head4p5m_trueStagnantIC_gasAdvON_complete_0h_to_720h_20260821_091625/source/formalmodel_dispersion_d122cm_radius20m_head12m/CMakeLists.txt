# SPDX-License-Identifier: GPL-3.0-or-later

add_input_file_links(FILES
    params.input
    floodmar.msh
    oxygen_initial_2880.dat
)

dumux_add_test(
    NAME floodmar_formal_dispersion_d122cm_radius20m_head12m
    LABELS porousmediumflow 2pnc floodmar dispersion
    SOURCES main.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMar
    COMMAND floodmar_formal_dispersion_d122cm_radius20m_head12m
    CMD_ARGS params.input
)

target_link_libraries(
    floodmar_formal_dispersion_d122cm_radius20m_head12m
    PRIVATE
    Dune::UGGrid
)
