# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for XInterfaceFiniteElement
add_marmot_test("TestXInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestXInterfaceFiniteElement.cpp")

# Tests for YInterfaceFiniteElement
add_marmot_test("TestYInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestYInterfaceFiniteElement.cpp")

# Tests for GaussLobattoInterfaceFiniteElement
add_marmot_test("TestGaussLobattoInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestGaussLobattoInterfaceFiniteElement.cpp")

# Tests for GaussLobattoBBarInterfaceFiniteElement
add_marmot_test("TestGaussLobattoBBarInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestGaussLobattoBBarInterfaceFiniteElement.cpp")

# Tests for YNodalGradientInterfaceFiniteElement
add_marmot_test("TestYNodalGradientInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestYNodalGradientInterfaceFiniteElement.cpp")
