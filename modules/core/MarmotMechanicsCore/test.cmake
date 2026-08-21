# add current directory to the source for the tests
SET(CURR_TEST_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/test")

# Tests for HaighWestergaard
add_marmot_test("TestHaighWestergaard" "${CURR_TEST_SOURCE_DIR}/TestHaighWestergaard.cpp")

# Tests for MarmotElasticity
add_marmot_test("TestMarmotElasticity" "${CURR_TEST_SOURCE_DIR}/TestMarmotElasticity.cpp")

# Tests for MarmotKelvinChain
add_marmot_test("TestMarmotKelvinChain" "${CURR_TEST_SOURCE_DIR}/TestMarmotKelvinChain.cpp")

# Tests for MarmotWiechert
add_marmot_test("TestMarmotWiechert" "${CURR_TEST_SOURCE_DIR}/TestMarmotWiechert.cpp")

# Tests for MarmotKinematics
add_marmot_test("TestMarmotKinematics" "${CURR_TEST_SOURCE_DIR}/TestMarmotKinematics.cpp")

# Tests for MarmotLowerDimensionalStress
add_marmot_test("TestMarmotLowerDimensionalStress" "${CURR_TEST_SOURCE_DIR}/TestMarmotLowerDimensionalStress.cpp")

# Tests for MarmotLocalization
add_marmot_test("TestMarmotLocalization" "${CURR_TEST_SOURCE_DIR}/TestMarmotLocalization.cpp")

# Tests for MarmotPronySeries
add_marmot_test("TestMarmotPronySeries" "${CURR_TEST_SOURCE_DIR}/TestMarmotPronySeries.cpp")

# Tests for MarmotStateVarVectorManager
add_marmot_test("TestMarmotStateVarVectorManager" "${CURR_TEST_SOURCE_DIR}/TestMarmotStateVarVectorManager.cpp")

# Tests for MarmotViscoelasticity
add_marmot_test("TestMarmotViscoelasticity" "${CURR_TEST_SOURCE_DIR}/TestMarmotViscoelasticity.cpp")

# Tests for MarmotVoigt
add_marmot_test("TestMarmotVoigt" "${CURR_TEST_SOURCE_DIR}/TestMarmotVoigt.cpp")

# Tests for MenetreyWillam
add_marmot_test("TestMenetreyWillam" "${CURR_TEST_SOURCE_DIR}/TestMenetreyWillam.cpp")

# Tests for YieldSurfaceCombinatioNmanager
add_marmot_test("TestYieldSurfaceCombinationManager" "${CURR_TEST_SOURCE_DIR}/TestYieldSurfaceCombinationManager.cpp")

# Tests for NewmarkBetaIntegrator
add_marmot_test("TestNewmarkBetaIntegrator" "${CURR_TEST_SOURCE_DIR}/TestNewmarkBetaIntegrator.cpp")

# Tests for MarmotGeostaticStress
add_marmot_test("TestMarmotGeostaticStress" "${CURR_TEST_SOURCE_DIR}/TestMarmotGeostaticStress.cpp")

# Tests for MarmotInterfaceMaterialHypoElastic
add_marmot_test("TestMarmotInterfaceMaterialHypoElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotInterfaceMaterialHypoElastic.cpp")

# Tests for MarmotCorrectedInterfaceMaterialHypoElastic
add_marmot_test("TestMarmotCorrectedInterfaceMaterialHypoElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotCorrectedInterfaceMaterialHypoElastic.cpp")

# Tests for MarmotXInterfaceMaterialHypoElastic
add_marmot_test("TestMarmotXInterfaceMaterialHypoElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotXInterfaceMaterialHypoElastic.cpp")

# Tests for MarmotEquilibratedXInterfaceMaterialHypoElastic
add_marmot_test("TestMarmotEquilibratedXInterfaceMaterialHypoElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotEquilibratedXInterfaceMaterialHypoElastic.cpp")

# Tests for MarmotPartialMixedInterfaceMaterialHypoElastic
add_marmot_test("TestMarmotPartialMixedInterfaceMaterialHypoElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotPartialMixedInterfaceMaterialHypoElastic.cpp")

# Tests for MarmotGaussLobattoInterfaceMaterialHypoElastic
add_marmot_test("TestMarmotGaussLobattoInterfaceMaterialHypoElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotGaussLobattoInterfaceMaterialHypoElastic.cpp")

# Tests for MarmotGaussLobattoBBarInterfaceMaterialHypoElastic
add_marmot_test("TestMarmotGaussLobattoBBarInterfaceMaterialHypoElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotGaussLobattoBBarInterfaceMaterialHypoElastic.cpp")

# Tests for MarmotWarpingInterfaceMaterialHypoElastic
add_marmot_test("TestMarmotWarpingInterfaceMaterialHypoElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotWarpingInterfaceMaterialHypoElastic.cpp")

# Tests for MarmotWarpingStabPressureInterfaceMaterialHypoElastic
add_marmot_test("TestMarmotWarpingStabPressureInterfaceMaterialHypoElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotWarpingStabPressureInterfaceMaterialHypoElastic.cpp")

# Tests for MarmotZInterfaceMaterialHypoElastic
add_marmot_test("TestMarmotZInterfaceMaterialHypoElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotZInterfaceMaterialHypoElastic.cpp")

# Tests for MarmotZStabPressureInterfaceMaterialHypoElastic
add_marmot_test("TestMarmotZStabPressureInterfaceMaterialHypoElastic" "${CURR_TEST_SOURCE_DIR}/TestMarmotZStabPressureInterfaceMaterialHypoElastic.cpp")
