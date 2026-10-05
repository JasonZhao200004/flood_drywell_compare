dumux_add_test(
    NAME floodmar_HYDROFIX2_C4_sameVolume_from1440
    LABELS porousmediumflow 2pnc floodmar hydrofix2 c4 samevolume restart
    SOURCES main_floodpond5m_depth10cm_phaseguard_v2_facedisp.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMarOldBoundaryPhaseGuardV2FaceDisp
    COMMAND floodmar_HYDROFIX2_C4_sameVolume_from1440
    CMD_ARGS params.input
)

target_link_libraries(
    floodmar_HYDROFIX2_C4_sameVolume_from1440
    PRIVATE
    Dune::UGGrid
)
