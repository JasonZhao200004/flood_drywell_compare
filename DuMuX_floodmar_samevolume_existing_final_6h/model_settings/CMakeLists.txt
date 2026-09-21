add_input_file_links(FILES
    params.input
    floodmar_pond5m_depth10cm_nonames.msh
    oxygenIC_trueStagnantO2_2880h_2d_xO2.dat
)

dumux_add_test(
    NAME floodmar_pond5m_head5cm_equalDose_3cycle_20d_oldformal
    LABELS porousmediumflow 2pnc floodmar pond equalDose phaseguard facedisp
    SOURCES main_floodpond5m_depth10cm_phaseguard_v2_facedisp.cc
)

target_link_libraries(
    floodmar_pond5m_head5cm_equalDose_3cycle_20d_oldformal
    PRIVATE
    Dune::UGGrid
)

target_compile_definitions(
    floodmar_pond5m_head5cm_equalDose_3cycle_20d_oldformal
    PRIVATE
    TYPETAG=FloodMarOldBoundaryPhaseGuardV2FaceDisp
)


# Independent restart-writer target for C2/C3 calibration.
dumux_add_test(
    NAME floodmar_pond5m_head5cm_equalDose_3cycle_20d_restartwriter
    LABELS porousmediumflow 2p2c floodmar calibration restartwriter
    SOURCES main_floodpond5m_depth10cm_phaseguard_v2_facedisp_restartwriter.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMarOldBoundaryPhaseGuardV2FaceDisp
    COMMAND floodmar_pond5m_head5cm_equalDose_3cycle_20d_restartwriter
    CMD_ARGS params.input
)

target_link_libraries(
    floodmar_pond5m_head5cm_equalDose_3cycle_20d_restartwriter
    PRIVATE
    Dune::UGGrid
)
