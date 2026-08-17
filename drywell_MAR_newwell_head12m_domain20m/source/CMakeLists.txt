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

dumux_add_test(
    NAME floodmar_formal_dispersion_d122cm_radius20m_head12m_umfpack
    LABELS porousmediumflow 2pnc floodmar dispersion umfpack
    SOURCES main_umfpack.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMar
    COMMAND floodmar_formal_dispersion_d122cm_radius20m_head12m_umfpack
    CMD_ARGS params.input
)

target_link_libraries(
    floodmar_formal_dispersion_d122cm_radius20m_head12m_umfpack
    PRIVATE
    Dune::UGGrid
)

# BEGIN FLOODMAR PHASEGUARD V2 TARGET
# Independent experimental target.  The existing AMG and UMFPACK targets are
# intentionally left unchanged.
dumux_add_test(
    NAME floodmar_formal_dispersion_d122cm_radius20m_head12m_umfpack_phaseguard
    LABELS porousmediumflow 2pnc floodmar dispersion umfpack phaseguard
    SOURCES main_umfpack_phaseguard.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMarPhaseGuard
    COMMAND floodmar_formal_dispersion_d122cm_radius20m_head12m_umfpack_phaseguard
    CMD_ARGS params.input
)

target_link_libraries(
    floodmar_formal_dispersion_d122cm_radius20m_head12m_umfpack_phaseguard
    PRIVATE
    Dune::UGGrid
)
# END FLOODMAR PHASEGUARD V2 TARGET

# BEGIN FLOODMAR CHOPPED NEWTON V3 TARGET
# Independent experimental target. Existing targets remain unchanged.
dumux_add_test(
    NAME floodmar_formal_dispersion_d122cm_radius20m_head12m_umfpack_chopped
    LABELS porousmediumflow 2pnc floodmar dispersion umfpack chopped
    SOURCES main_umfpack_chopped.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMar
    COMMAND floodmar_formal_dispersion_d122cm_radius20m_head12m_umfpack_chopped
    CMD_ARGS params.input
)

target_link_libraries(
    floodmar_formal_dispersion_d122cm_radius20m_head12m_umfpack_chopped
    PRIVATE
    Dune::UGGrid
)
# END FLOODMAR CHOPPED NEWTON V3 TARGET

# BEGIN FLOODMAR RESTART PHASE HYSTERESIS V4 TARGET
dumux_add_test(
    NAME floodmar_formal_dispersion_d122cm_radius20m_head12m_restart_v4
    LABELS porousmediumflow 2pnc floodmar dispersion umfpack restart
    SOURCES main_umfpack_chopped_restart_v4.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMarRestartV4
    COMMAND floodmar_formal_dispersion_d122cm_radius20m_head12m_restart_v4
    CMD_ARGS params.input
)

target_link_libraries(
    floodmar_formal_dispersion_d122cm_radius20m_head12m_restart_v4
    PRIVATE
    Dune::UGGrid
)
# END FLOODMAR RESTART PHASE HYSTERESIS V4 TARGET
