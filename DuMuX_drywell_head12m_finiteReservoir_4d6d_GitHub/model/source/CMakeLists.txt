add_input_file_links(FILES
    params.input
    floodmar_nonames.msh
    oxygenIC_trueStagnantO2_2880h_2d_xO2.dat
    oxygen_initial_2880.dat
)

dumux_add_test(
    NAME drywell_head12m_production_baseline
    SOURCES main_oldboundary_phaseguard_v2_facedisp.cc
)

target_link_libraries(
    drywell_head12m_production_baseline
    PRIVATE
    Dune::UGGrid
)

target_compile_definitions(
    drywell_head12m_production_baseline
    PRIVATE
    TYPETAG=FloodMarOldBoundaryPhaseGuardV2FaceDisp
)


dumux_add_test(
    NAME drywell_head12m_reservoir_4d6d
    SOURCES main_reservoir_4d6d.cc
)

target_link_libraries(
    drywell_head12m_reservoir_4d6d
    PRIVATE
    Dune::UGGrid
)

target_compile_definitions(
    drywell_head12m_reservoir_4d6d
    PRIVATE
    TYPETAG=FloodMarOldBoundaryPhaseGuardV2FaceDisp
)
