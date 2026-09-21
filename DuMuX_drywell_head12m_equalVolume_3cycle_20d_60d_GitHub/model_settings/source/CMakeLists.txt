add_input_file_links(FILES
    params.input
    floodmar_nonames.msh
    oxygenIC_trueStagnantO2_2880h_2d_xO2.dat
)

dumux_add_test(
    NAME drywell_head12m_equalVolume_3cycle_20d_60d
    LABELS porousmediumflow 2pnc floodmar drywell phaseguard facedisp equalvolume
    SOURCES main_oldboundary_phaseguard_v2_facedisp.cc
)

target_link_libraries(
    drywell_head12m_equalVolume_3cycle_20d_60d
    PRIVATE
    Dune::UGGrid
)

target_compile_definitions(
    drywell_head12m_equalVolume_3cycle_20d_60d
    PRIVATE
    TYPETAG=FloodMarOldBoundaryPhaseGuardV2FaceDisp
)
