dumux_add_test(
    NAME floodmar_TRUE_RIGHTSURFACEtag_restart101
    LABELS porousmediumflow 2pnc floodmar true-right-surface-tag restart
    SOURCES main_floodpond5m_depth10cm_phaseguard_v2_facedisp.cc
    COMPILE_DEFINITIONS TYPETAG=FloodMarOldBoundaryPhaseGuardV2FaceDisp
    COMMAND floodmar_TRUE_RIGHTSURFACEtag_restart101
    CMD_ARGS params.input
)
target_link_libraries(floodmar_TRUE_RIGHTSURFACEtag_restart101 PRIVATE Dune::UGGrid)
