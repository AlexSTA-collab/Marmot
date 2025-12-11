#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotKelvinChain.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "Marmot/MarmotWiechertInterface.h"
#include "autodiff/forward/real.hpp"
#include <iomanip>

using namespace Marmot::Testing;
using namespace Marmot::Materials::WiechertInterface;
using namespace Marmot::ContinuumMechanics::Viscoelasticity;
using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;
using namespace Marmot;

void evaluateWIandUpdateStateVarsTestFunction()
{
  double factor = 1.35;

  // properties wiechert interface
  int    nMaxwell = 1;   // Number of Maxwell units in parallel
  double E        = 1e4; // Value of the elastic modulus
  double tau      = 10.; // Value of the relaxation time

  Properties elasticModuli   = initializeElasticModuli( nMaxwell, E );
  Properties relaxationTimes = initializeRelaxationTimes( nMaxwell, tau );

  // time increment
  int dT = 10;

  // arbitrary initial state vars
  StateVarMatrix_force_uu          stateVars_force_uu( 3, nMaxwell );
  StateVarMatrix_force_us          stateVars_force_us( 3, nMaxwell );
  StateVarMatrix_surface_stress_Z  stateVars_surface_stress_Z( 9, nMaxwell );
  StateVarMatrix_surface_stress_Y  stateVars_surface_stress_Y( 9, nMaxwell );
  StateVarMatrix_surface_stress_us stateVars_surface_stress_us( 9, nMaxwell );

  stateVars_force_uu << 0.01, 0.02, 0.03;
  stateVars_force_us << 0.01, 0.02, 0.03;
  stateVars_surface_stress_Z << 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09;
  stateVars_surface_stress_Y << 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09;
  stateVars_surface_stress_us << 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09;

  // initialize compliance and strain
  double   uniaxialStiffness = 0;
  Vector3d dforce_uu         = { 0., 0., 0. };
  Vector3d dforce_us         = { 0., 0., 0. };
  Vector9d dsurfaceStress_Z  = { 0., 0., 0., 0., 0., 0., 0., 0., 0. };
  Vector9d dsurfaceStress_Y  = { 0., 0., 0., 0., 0., 0., 0., 0., 0. };
  Vector9d dsurfaceStress_us = { 0., 0., 0., 0., 0., 0., 0., 0., 0. };

  double   corrStiffness         = 8533.6275441855286772;
  Vector3d corrDforce            = { 0.0085336275441855, 0.0170672550883711, 0.0256008826325566 };
  Vector9d corrDsurfaceStress_Z  = { 0.0085336275441855,
                                     0.0170672550883711,
                                     0.0256008826325566,
                                     0.0341345101767421,
                                     0.0426681377209276,
                                     0.0512017652651132,
                                     0.0597353928092987,
                                     0.0682690203534842,
                                     0.0768026478976698 };
  Vector9d corrDsurfaceStress_Y  = corrDsurfaceStress_Z; // Same expected values
  Vector9d corrDsurfaceStress_us = corrDsurfaceStress_Z; // Same expected values

  // Test evaluation of unified Wiechert interface model
  evaluateWiechert( dT,
                    elasticModuli,
                    relaxationTimes,
                    stateVars_force_uu,
                    stateVars_force_us,
                    stateVars_surface_stress_Z,
                    stateVars_surface_stress_Y,
                    stateVars_surface_stress_us,
                    uniaxialStiffness,
                    dforce_uu,
                    dforce_us,
                    dsurfaceStress_Z,
                    dsurfaceStress_Y,
                    dsurfaceStress_us,
                    factor );

  throwExceptionOnFailure( checkIfEqual( uniaxialStiffness, corrStiffness ),
                           "error in uniaxial stiffness of Wiechert model" );

  throwExceptionOnFailure( checkIfEqual< double >( dforce_uu, corrDforce ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in dforce_uu of Wiechert model" );

  throwExceptionOnFailure( checkIfEqual< double >( dforce_us, corrDforce ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in dforce_us of Wiechert model" );

  throwExceptionOnFailure( checkIfEqual< double >( dsurfaceStress_Z, corrDsurfaceStress_Z ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in dsurfaceStress_Z of Wiechert model" );

  throwExceptionOnFailure( checkIfEqual< double >( dsurfaceStress_Y, corrDsurfaceStress_Y ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in dsurfaceStress_Y of Wiechert model" );

  throwExceptionOnFailure( checkIfEqual< double >( dsurfaceStress_us, corrDsurfaceStress_us ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in dsurfaceStress_us of Wiechert model" );

  // Material properties and normal vector for interface material parameters
  double const                         nu_0      = 0.2;
  double const                         normal[3] = { 0., 0., 1. };
  Fastor::TensorMap< double const, 3 > normalFtensor( normal );

  auto [unitZ_ijkl,
        unitH_inv_ij,
        unitH_inv_nF_ijk,
        unitYn_H_inv_Fn_ijkl] = calculateInterfaceMaterialParameters( normalFtensor, nu_0 );

  // Convert to matrices for update functions
  Eigen::Matrix< double, 3, 3 > unitH_inv_voigt_full      = convert2ndOrderTensorToMatrix_3x3( unitH_inv_ij );
  Eigen::Matrix< double, 9, 9 > unitZ_voigt_full          = convert4thOrderTensorToMatrix_9x9( unitZ_ijkl );
  Eigen::Matrix< double, 3, 9 > unitH_inv_nF_ijk_3x9      = convert3rdOrderTensorToMatrix_3x9( unitH_inv_nF_ijk );
  Eigen::Matrix< double, 9, 3 > unitH_inv_nF_ijk_9x3      = convert3rdOrderTensorToMatrix_9x3( unitH_inv_nF_ijk );
  Eigen::Matrix< double, 9, 9 > unitYn_H_inv_Fn_ijkl_full = convert4thOrderTensorToMatrix_9x9( unitYn_H_inv_Fn_ijkl );

  // Test 1: updateStateVarMatrix_force_uu
  {
    Vector3d                djumpU_test = { 0.1, 0.2, 0.3 };
    StateVarMatrix_force_uu testStateVars_force_uu( 3, nMaxwell );
    testStateVars_force_uu << 0.01, 0.02, 0.03;

    updateStateVarMatrix_force_uu( dT,
                                   elasticModuli,
                                   relaxationTimes,
                                   testStateVars_force_uu,
                                   djumpU_test,
                                   unitH_inv_voigt_full );

    // Reference values from terminal output
    StateVarMatrix_force_uu expectedStateVars_force_uu( 3, nMaxwell );
    expectedStateVars_force_uu << 263.387244972977, 526.774489945955, 2107.07956581176;

    throwExceptionOnFailure( testStateVars_force_uu.isApprox( expectedStateVars_force_uu, 1e-10 ),
                             MakeString() << __PRETTY_FUNCTION__ << " error in updateStateVarMatrix_force_uu" );
  }

  // Test 2: updateStateVarMatrix_force_us
  {
    Vector9d                dsurfaceStrain_test = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
    StateVarMatrix_force_us testStateVars_force_us( 3, nMaxwell );
    testStateVars_force_us << 0.01, 0.02, 0.03;

    updateStateVarMatrix_force_us( dT,
                                   elasticModuli,
                                   relaxationTimes,
                                   testStateVars_force_us,
                                   dsurfaceStrain_test,
                                   unitH_inv_nF_ijk_3x9 );

    // Reference values from terminal output
    StateVarMatrix_force_us expectedStateVars_force_us( 3, nMaxwell );
    expectedStateVars_force_us << -3687.36624770551, -4214.12970126823, -12642.4001401879;

    throwExceptionOnFailure( testStateVars_force_us.isApprox( expectedStateVars_force_us, 1e-10 ),
                             MakeString() << __PRETTY_FUNCTION__ << " error in updateStateVarMatrix_force_us" );
  }

  // Test 3: updateStateVarMatrix_surface_stress_Z
  {
    Vector9d                        dsurfaceStrain_test = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
    StateVarMatrix_surface_stress_Z testStateVars_surface_stress_Z( 9, nMaxwell );
    testStateVars_surface_stress_Z << 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09;

    updateStateVarMatrix_surface_stress_Z( dT,
                                           elasticModuli,
                                           relaxationTimes,
                                           testStateVars_surface_stress_Z,
                                           dsurfaceStrain_test,
                                           unitZ_voigt_full );

    // Reference values from terminal output
    StateVarMatrix_surface_stress_Z expectedStateVars_surface_stress_Z( 9, nMaxwell );
    expectedStateVars_surface_stress_Z << 1316.92150968724, 1580.30875466022, 0.0110363832351433, 1580.31611224904,
      3424.00475429341, 0.0220727664702865, 0.025751560882001, 0.0294303552937154, 0.0331091497054298;

    throwExceptionOnFailure( testStateVars_surface_stress_Z.isApprox( expectedStateVars_surface_stress_Z, 1e-10 ),
                             MakeString() << __PRETTY_FUNCTION__ << " error in updateStateVarMatrix_surface_stress_Z" );
  }

  // Test 4: updateStateVarMatrix_surface_stress_Y
  {
    Vector9d                        dsurfaceStrain_test = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
    StateVarMatrix_surface_stress_Y testStateVars_surface_stress_Y( 9, nMaxwell );
    testStateVars_surface_stress_Y << 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09;

    updateStateVarMatrix_surface_stress_Y( dT,
                                           elasticModuli,
                                           relaxationTimes,
                                           testStateVars_surface_stress_Y,
                                           dsurfaceStrain_test,
                                           unitYn_H_inv_Fn_ijkl_full );

    // Reference values from terminal output
    StateVarMatrix_surface_stress_Y expectedStateVars_surface_stress_Y( 9, nMaxwell );
    expectedStateVars_surface_stress_Y << 6321.20926707999, 0.00735758882342885, 7374.75088938307, 0.0147151776468577,
      6321.22398225763, 8428.29619048057, 7374.76560456072, 8428.3035480694, 25284.855462292;

    throwExceptionOnFailure( testStateVars_surface_stress_Y.isApprox( expectedStateVars_surface_stress_Y, 1e-10 ),
                             MakeString() << __PRETTY_FUNCTION__ << " error in updateStateVarMatrix_surface_stress_Y" );
  }

  // Test 5: updateStateVarMatrix_surface_stress_us
  {
    Vector3d                         djumpU_test = { 0.1, 0.2, 0.3 };
    StateVarMatrix_surface_stress_us testStateVars_surface_stress_us( 9, nMaxwell );
    testStateVars_surface_stress_us << 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09;

    updateStateVarMatrix_surface_stress_us( dT,
                                            elasticModuli,
                                            relaxationTimes,
                                            testStateVars_surface_stress_us,
                                            djumpU_test,
                                            unitH_inv_nF_ijk_3x9 );

    // Reference values from terminal output
    StateVarMatrix_surface_stress_us expectedStateVars_surface_stress_us( 9, nMaxwell );
    expectedStateVars_surface_stress_us << 0.00367879441171442, 0.00735758882342885, -526.756095973896,
      0.0147151776468577, 0.0183939720585721, -1053.51219194779, 0.025751560882001, 0.0294303552937154,
      -4214.10394970735;

    throwExceptionOnFailure( testStateVars_surface_stress_us.isApprox( expectedStateVars_surface_stress_us, 1e-10 ),
                             MakeString()
                               << __PRETTY_FUNCTION__ << " error in updateStateVarMatrix_surface_stress_us" );
  }
}

void computeLambdaAndBetaTestFunction()
{
  double lambda, beta;
  double dT = 30;

  // case dT_tau >= 30.0
  double tau = 1 / ( 30 / dT );
  computeLambdaAndBeta( dT, tau, lambda, beta );

  throwExceptionOnFailure( checkIfEqual( beta, 0 ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in beta with dT_tau >= 30.0" );
  throwExceptionOnFailure( checkIfEqual( lambda, tau / dT ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in lambda with dT_tau >= 30.0" );

  // case dT_tau < 1e-6
  tau = 1 / ( 1e-7 / dT );
  computeLambdaAndBeta( dT, tau, lambda, beta );

  double corrLam = 1 - 0.5 * dT / tau + 1. / 6 * dT / tau * dT / tau;

  throwExceptionOnFailure( checkIfEqual( beta, 1 ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in beta with dT_tau < 1e-6" );
  throwExceptionOnFailure( checkIfEqual( lambda, corrLam ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in lambda with dT_tau < 1e-6" );

  // case else
  tau = 1 / ( 10 / dT );
  computeLambdaAndBeta( dT, tau, lambda, beta );

  double corrBeta = std::exp( -dT / tau );
  corrLam         = ( 1 - beta ) * ( tau / dT );

  throwExceptionOnFailure( checkIfEqual( beta, corrBeta ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in beta with 1e-6 <= dT_tau < 30.0" );
  throwExceptionOnFailure( checkIfEqual( lambda, corrLam ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in lambda with 1e-6 <= dT_tau < 30.0" );
}

int main()
{

  auto tests = std::vector< std::function< void() > >{ evaluateWIandUpdateStateVarsTestFunction,
                                                       computeLambdaAndBetaTestFunction };

  executeTestsAndCollectExceptions( tests );

  return 0;
}