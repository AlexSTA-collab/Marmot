# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for XInterfaceFiniteElement
add_marmot_test("TestXInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestXInterfaceFiniteElement.cpp")

# Tests for YInterfaceFiniteElement
add_marmot_test("TestYInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestYInterfaceFiniteElement.cpp")
