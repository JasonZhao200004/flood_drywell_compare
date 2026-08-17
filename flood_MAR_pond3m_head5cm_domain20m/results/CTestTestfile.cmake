# CMake generated Testfile for 
# Source directory: /Users/zhaozhe/dumux-work/dumux/dumux-floodmar/test/porousmediumflow/2p2c/floodmodel
# Build directory: /Users/zhaozhe/dumux-work/dumux/build-serial/dumux-floodmar/test/porousmediumflow/2p2c/floodmodel
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(floodmar_flood "/Users/zhaozhe/dumux-work/dumux/build-serial/dumux-floodmar/test/porousmediumflow/2p2c/floodmodel/floodmar_flood" "params.input")
set_tests_properties(floodmar_flood PROPERTIES  LABELS "porousmediumflow;2pnc;floodmar" PROCESSORS "1" REQUIRED_FILES "/Users/zhaozhe/dumux-work/dumux/build-serial/dumux-floodmar/test/porousmediumflow/2p2c/floodmodel/floodmar_flood" SKIP_RETURN_CODE "77" TIMEOUT "300" WORKING_DIRECTORY "/Users/zhaozhe/dumux-work/dumux/build-serial/dumux-floodmar/test/porousmediumflow/2p2c/floodmodel" _BACKTRACE_TRIPLES "/Users/zhaozhe/dumux-work/dumux/dune-common/cmake/modules/DuneTestMacros.cmake;417;add_test;/Users/zhaozhe/dumux-work/dumux/dumux/cmake/modules/DumuxTestMacros.cmake;210;dune_add_test;/Users/zhaozhe/dumux-work/dumux/dumux-floodmar/test/porousmediumflow/2p2c/floodmodel/CMakeLists.txt;9;dumux_add_test;/Users/zhaozhe/dumux-work/dumux/dumux-floodmar/test/porousmediumflow/2p2c/floodmodel/CMakeLists.txt;0;")
