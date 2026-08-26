# SPDX-License-Identifier: GPL-3.0-or-later
add_input_file_links(FILES
    params.input
    floodmar_nonames.msh
    oxygenIC_trueStagnantO2_2880h_2d_xO2.dat
)
dumux_add_test(
    NAME floodmar_head4p5m_hydrusmatch_trueStagnant
    LABELS porousmediumflow 2pnc floodmar dispersion hydrusmatch stagnant_o2
    SOURCES main_hydrusmatch.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMarHydrusMatch
    COMMAND floodmar_head4p5m_hydrusmatch_trueStagnant
    CMD_ARGS params.input
)
target_link_libraries(floodmar_head4p5m_hydrusmatch_trueStagnant PRIVATE Dune::UGGrid)

dumux_add_test(
    NAME floodmar_head4p5m_trueStagnantIC_gasAdvON
    LABELS porousmediumflow 2pnc floodmar dispersion gas_advection oxygenIC
    SOURCES main_gasAdvON.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMar
    COMMAND floodmar_head4p5m_trueStagnantIC_gasAdvON
    CMD_ARGS params.input
)

target_link_libraries(
    floodmar_head4p5m_trueStagnantIC_gasAdvON
    PRIVATE
    Dune::UGGrid
)


# Restart-capable robust 4.5 m drywell executable.
dumux_add_test(
    NAME floodmar_head4p5m_trueStagnantIC_gasAdvON_restart_v4
    SOURCES main_umfpack_chopped_restart_v4.cc
)

target_link_libraries(
    floodmar_head4p5m_trueStagnantIC_gasAdvON_restart_v4
    PRIVATE
    Dune::UGGrid
)

target_compile_definitions(
    floodmar_head4p5m_trueStagnantIC_gasAdvON_restart_v4
    PRIVATE
    TYPETAG=FloodMarRestartV4
)


# Controlled old-boundary comparison target.
dumux_add_test(
    NAME floodmar_head4p5m_oldboundary_gasAdvON_restart_v4
    SOURCES main_oldboundary_restart_v4.cc
)

target_link_libraries(
    floodmar_head4p5m_oldboundary_gasAdvON_restart_v4
    PRIVATE
    Dune::UGGrid
)

target_compile_definitions(
    floodmar_head4p5m_oldboundary_gasAdvON_restart_v4
    PRIVATE
    TYPETAG=FloodMarOldBoundaryRestartV4
)


# Box-discretization HYDRUS comparison.
# Keeps true internal gas-phase O2 advection.
dumux_add_test(
    NAME floodmar_head4p5m_oldboundary_box_gasAdvON
    SOURCES main_oldboundary_box.cc
)

target_link_libraries(
    floodmar_head4p5m_oldboundary_box_gasAdvON
    PRIVATE
    Dune::UGGrid
)

target_compile_definitions(
    floodmar_head4p5m_oldboundary_box_gasAdvON
    PRIVATE
    TYPETAG=FloodMarOldBoundaryBox
)


# ----------------------------------------------------------------------
# Old-boundary CCTpfa PhaseGuard V2
#
# Same physics / same chopped Newton solver as RestartV4.
# Only the primary-variable switch is replaced by the strict
# per-nonlinear-solve anti-chatter PhaseGuard V2.
# ----------------------------------------------------------------------

dumux_add_test(
    NAME floodmar_head4p5m_oldboundary_gasAdvON_phaseguard_v2
    SOURCES main_oldboundary_phaseguard_v2.cc
)

target_link_libraries(
    floodmar_head4p5m_oldboundary_gasAdvON_phaseguard_v2
    PRIVATE
    Dune::UGGrid
)

target_compile_definitions(
    floodmar_head4p5m_oldboundary_gasAdvON_phaseguard_v2
    PRIVATE
    TYPETAG=FloodMarOldBoundaryPhaseGuardV2
)


# ----------------------------------------------------------------------
# PhaseGuard V2 + face-symmetric CCTpfa dispersion
# ----------------------------------------------------------------------

dumux_add_test(
    NAME floodmar_head4p5m_oldboundary_gasAdvON_phaseguard_v2_facedisp
    SOURCES main_oldboundary_phaseguard_v2_facedisp.cc
)

target_link_libraries(
    floodmar_head4p5m_oldboundary_gasAdvON_phaseguard_v2_facedisp
    PRIVATE
    Dune::UGGrid
)

target_compile_definitions(
    floodmar_head4p5m_oldboundary_gasAdvON_phaseguard_v2_facedisp
    PRIVATE
    TYPETAG=FloodMarOldBoundaryPhaseGuardV2FaceDisp
)


# ============================================================
# Homogeneous Layer-3 Sand2 sensitivity case
# Uniform initial theta_w = 0.116926651515
# Equivalent Sand2 initial head = -29.770236344729 cm
# Baseline physics otherwise unchanged.
# ============================================================

dumux_add_test(
    NAME floodmar_head4p5m_homogeneousSand2_uniformIC_phaseguard_v2_facedisp
    SOURCES main_oldboundary_phaseguard_v2_facedisp_homogeneous_sand2.cc
)

target_link_libraries(
    floodmar_head4p5m_homogeneousSand2_uniformIC_phaseguard_v2_facedisp
    PRIVATE
    Dune::UGGrid
)

target_compile_definitions(
    floodmar_head4p5m_homogeneousSand2_uniformIC_phaseguard_v2_facedisp
    PRIVATE
    TYPETAG=FloodMarOldBoundaryPhaseGuardV2FaceDisp
)

