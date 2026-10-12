# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for InterfaceFiniteElement
add_marmot_test("TestInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestInterfaceFiniteElement.cpp")

# Tests for YInterfaceFiniteElement
add_marmot_test("TestYInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestYInterfaceFiniteElement.cpp" REQUIRES LinearElastic VonMises)

# Tests for GradedInterfaceFiniteElement
add_marmot_test("TestGradedInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestGradedInterfaceFiniteElement.cpp" REQUIRES LinearElastic VonMises)

# Tests for GradedYInterfaceFiniteElement
add_marmot_test("TestGradedYInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestGradedYInterfaceFiniteElement.cpp" REQUIRES LinearElastic VonMises)
