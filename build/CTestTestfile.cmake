# CMake generated Testfile for 
# Source directory: /home/alexsta1993/alexandros/Marmot
# Build directory: /home/alexsta1993/alexandros/Marmot/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(TestMarmotAutomaticDifferentiation "/home/alexsta1993/alexandros/Marmot/build/bin/TestMarmotAutomaticDifferentiation")
set_tests_properties(TestMarmotAutomaticDifferentiation PROPERTIES  WORKING_DIRECTORY "/home/alexsta1993/alexandros/Marmot/build" _BACKTRACE_TRIPLES "/home/alexsta1993/alexandros/Marmot/CMakeLists.txt;69;add_test;/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMathCore/test.cmake;5;add_marmot_test;/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMathCore/test.cmake;0;;/home/alexsta1993/alexandros/Marmot/CMakeLists.txt;297;include;/home/alexsta1993/alexandros/Marmot/CMakeLists.txt;0;")
add_test(TestMarmotNumericalIntegration "/home/alexsta1993/alexandros/Marmot/build/bin/TestMarmotNumericalIntegration")
set_tests_properties(TestMarmotNumericalIntegration PROPERTIES  WORKING_DIRECTORY "/home/alexsta1993/alexandros/Marmot/build" _BACKTRACE_TRIPLES "/home/alexsta1993/alexandros/Marmot/CMakeLists.txt;69;add_test;/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMathCore/test.cmake;14;add_marmot_test;/home/alexsta1993/alexandros/Marmot/modules/core/MarmotMathCore/test.cmake;0;;/home/alexsta1993/alexandros/Marmot/CMakeLists.txt;297;include;/home/alexsta1993/alexandros/Marmot/CMakeLists.txt;0;")
add_test(TestLinearElastic "/home/alexsta1993/alexandros/Marmot/build/bin/TestLinearElastic")
set_tests_properties(TestLinearElastic PROPERTIES  WORKING_DIRECTORY "/home/alexsta1993/alexandros/Marmot/build" _BACKTRACE_TRIPLES "/home/alexsta1993/alexandros/Marmot/CMakeLists.txt;69;add_test;/home/alexsta1993/alexandros/Marmot/modules/materials/LinearElastic/test.cmake;5;add_marmot_test;/home/alexsta1993/alexandros/Marmot/modules/materials/LinearElastic/test.cmake;0;;/home/alexsta1993/alexandros/Marmot/CMakeLists.txt;297;include;/home/alexsta1993/alexandros/Marmot/CMakeLists.txt;0;")
