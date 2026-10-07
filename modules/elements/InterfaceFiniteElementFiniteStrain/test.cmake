# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for InterfaceFiniteElementFiniteStrain
add_marmot_test("TestInterfaceFiniteElementFiniteStrain" "${CURR_TEST_SOURCE_DIR}/TestInterfaceFiniteElementFiniteStrain.cpp" REQUIRES CompressibleNeoHooke GradientEnhancedCompressibleNeoHookeDamage)
