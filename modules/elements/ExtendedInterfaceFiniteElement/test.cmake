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

# Tests for WarpingInterfaceFiniteElement
add_marmot_test("TestWarpingInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestWarpingInterfaceFiniteElement.cpp")

# Q1 + MINI bubble, no pressure DOF (GLIQUAD4_BUBBLE)
add_marmot_test("TestBubbleOnlyInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestBubbleOnlyInterfaceFiniteElement.cpp")

# Diagnostic: which pressure mode the MINI bubble couples to
add_marmot_test("TestMiniBubblePressureModeCoupling" "${CURR_TEST_SOURCE_DIR}/TestMiniBubblePressureModeCoupling.cpp")

# Tests for the GLIQUAD4 SURFACE quadrature option (GLIQUAD4_SURFLOB)
add_marmot_test("TestGaussLobattoSurfaceQuadrature" "${CURR_TEST_SOURCE_DIR}/TestGaussLobattoSurfaceQuadrature.cpp")

# Tests for WarpingStabPressureMiniInterfaceFiniteElement
add_marmot_test("TestWarpingStabPressureMiniInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestWarpingStabPressureMiniInterfaceFiniteElement.cpp")

# Tests for ZInterfaceFiniteElement (ZIQUAD4)
add_marmot_test("TestZInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestZInterfaceFiniteElement.cpp")

# Tests for ZStabPressureInterfaceFiniteElement (ZIQUAD4_STABP)
add_marmot_test("TestZStabPressureInterfaceFiniteElement" "${CURR_TEST_SOURCE_DIR}/TestZStabPressureInterfaceFiniteElement.cpp")
