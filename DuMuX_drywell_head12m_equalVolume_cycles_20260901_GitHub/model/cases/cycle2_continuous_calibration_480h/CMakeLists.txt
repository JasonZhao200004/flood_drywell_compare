add_input_file_links(FILES
    params.input
    floodmar_nonames.msh
    oxygenIC_trueStagnantO2_2880h_2d_xO2.dat
)

dumux_add_test(
    NAME floodmar_head12m_faceSymDisp_cycle2_calibration_480h
    LABELS porousmediumflow 2pnc floodmar drywell phaseguard facedisp continuouswet
    SOURCES main_oldboundary_phaseguard_v2_facedisp.cc
)

target_link_libraries(
    floodmar_head12m_faceSymDisp_cycle2_calibration_480h
    PRIVATE
    Dune::UGGrid
)

target_compile_definitions(
    floodmar_head12m_faceSymDisp_cycle2_calibration_480h
    PRIVATE
    TYPETAG=FloodMarOldBoundaryPhaseGuardV2FaceDisp
)
