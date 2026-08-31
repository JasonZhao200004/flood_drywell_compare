# SPDX-License-Identifier: GPL-3.0-or-later

dumux_add_test(
    NAME floodmar_head4p5m_faceSymDisp_hydruslike_noGasO2Adv
    LABELS porousmediumflow 2pnc floodmar dispersion hydruslike no_gas_o2_advection
    SOURCES main_hydruslike_noGasO2Adv.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMarOldBoundaryPhaseGuardV2FaceDispNoGasO2Adv
    COMMAND floodmar_head4p5m_faceSymDisp_hydruslike_noGasO2Adv
)

target_link_libraries(
    floodmar_head4p5m_faceSymDisp_hydruslike_noGasO2Adv
    PRIVATE
    Dune::UGGrid
)
