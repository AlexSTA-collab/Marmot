#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotKelvinChain.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotViscoelasticity.h"
#include "Marmot/MarmotWiechertInterface.h"
#include "autodiff/forward/real.hpp"

using namespace Marmot::Testing;
using namespace Marmot::Materials::WiechertInterface;
using namespace Marmot::ContinuumMechanics::Viscoelasticity;
using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;
using namespace Marmot;

void evaluateWIandUpdateStateVarsTestFunction()
{
  double factor = 1.35;

  // properties wiechert interface
  int    nMaxwellRu = 1;   // Number of Maxwell units in parallel for the displacement jump
  int    nMaxwellRs = 1;   // Number of Maxwell units in parallel for the surface strain
  double ERu        = 1e4; // Value of the elastic modulus for the displacement jump
  double ERs        = 1e4; // Value of the elastic modulus for the surface strain
  double tauRu      = 10.; // Value of the relaxation time for the displacement jump
  double tauRs      = 10.; // Value of the relaxation time for the surface strain

  // double spacing = 5.; Not usied right now direct initialization of relaxation times

  Properties elasticModuliRu = initializeElasticModuliRu( nMaxwellRu, ERu );

  Properties elasticModuliRs = initializeElasticModuliRs( nMaxwellRs, ERs );

  Properties relaxationTimesRu = initializeRelaxationTimesRu( nMaxwellRu, tauRu );

  Properties relaxationTimesRs = initializeRelaxationTimesRs( nMaxwellRs, tauRs );

  // time increment
  int dT = 10;

  // arbitrary initial state vars
  StateVarMatrixRu stateVarsRu( 3, nMaxwellRu );
  StateVarMatrixRs stateVarsRs( 9, nMaxwellRs );
  stateVarsRu << 0.01, 0.02, 0.03;
  stateVarsRs << 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09;

  // initialize compliance and strain
  double   uniaxialStiffnessRu = 0;
  double   uniaxialStiffnessRs = 0;
  Vector3d dForce              = { 0., 0., 0. };
  Vector9d dSurfaceStress      = { 0., 0., 0., 0., 0., 0., 0., 0., 0. };

  // test retardation times
  // Properties retardationTimes = generateRetardationTimes( n, min, spacing );
  // for ( int i = 0; i < n; i++ )
  //  throwExceptionOnFailure( checkIfEqual( retardationTimes( i ), min * std::pow( spacing, i ) ),
  //                           MakeString() << __PRETTY_FUNCTION__ << " error in generation of retardation times" );
  double   corrStiffnessRu    = 8533.6275441855286772;
  double   corrStiffnessRs    = 8533.6275441855286772;
  Vector3d corrDforce         = { 0.0085336275441855, 0.0170672550883711, 0.0256008826325566 };
  Vector9d corrDsurfaceStress = { 0.0085336275441855,
                                  0.0170672550883711,
                                  0.0256008826325566,
                                  0.0341345101767421,
                                  0.0426681377209276,
                                  0.0512017652651132,
                                  0.0597353928092987,
                                  0.0682690203534842,
                                  0.0768026478976698 };

  // test evaluation of Wiechert interface model
  evaluateWiechertRu( dT, elasticModuliRu, relaxationTimesRu, stateVarsRu, uniaxialStiffnessRu, dForce, factor );
  evaluateWiechertRs( dT,
                      elasticModuliRs,
                      relaxationTimesRs,
                      stateVarsRs,
                      uniaxialStiffnessRs,
                      dSurfaceStress,
                      factor );

  throwExceptionOnFailure( checkIfEqual( uniaxialStiffnessRu, corrStiffnessRu ),
                           "error in uniaxial compliance of Wiechert model for the stiffness due to displacement "
                           "jump" );

  throwExceptionOnFailure( checkIfEqual( uniaxialStiffnessRs, corrStiffnessRs ),
                           "error in uniaxial compliance of Wiechert model for the stiffness due to surface strain" );

  throwExceptionOnFailure( checkIfEqual< double >( dForce, corrDforce ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in dForce of Wiechert model" );

  throwExceptionOnFailure( checkIfEqual< double >( dSurfaceStress, corrDsurfaceStress ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in dSurfaceStress of Wiechert model" );
  // set zeroth wiechert stiffness
  double zerothWiechertStiffnessRu = 1.0;
  double zerothWiechertStiffnessRs = 1.0;
  // material properties for the matrix and the inclusion and normal vector
  double const                         E_M       = 1.0;
  double const                         nu_M      = 0.2;
  double const                         E_I       = 1.0;
  double const                         nu_I      = 0.2;
  double const                         E_0       = 1e2;
  double const                         nu_0      = 0.2;
  double const                         normal[3] = { 0., 0., 1. };
  Fastor::TensorMap< double const, 3 > normalFtensor( normal );

  auto [H_inv_ij_effective,
        Z_ijkl_effective,
        unitH_inv_voigt_full,
        unitZ_voigt_full] = calculateEffectiveProperties( zerothWiechertStiffnessRu,
                                                          uniaxialStiffnessRu,
                                                          zerothWiechertStiffnessRs,
                                                          uniaxialStiffnessRs,
                                                          normalFtensor,
                                                          E_M,
                                                          nu_M,
                                                          E_I,
                                                          nu_I,
                                                          E_0,
                                                          nu_0 );

  // test state var update for the displacement jump
  Vector3d         dForceN = { 0.1, 0.2, 0.3 }; // arbitrary force increment
  StateVarMatrixRu corrStateVarsRuN( 3, nMaxwellRu );
  corrStateVarsRuN << 263.3872449729773848, 526.7744899459547696, 2107.0795658117604034;

  updateStateVarMatrixRu( dT, elasticModuliRu, relaxationTimesRu, stateVarsRu, dForceN, unitH_inv_voigt_full );

  throwExceptionOnFailure( checkIfEqual< double >( stateVarsRu, corrStateVarsRuN ),
                           MakeString() << __PRETTY_FUNCTION__
                                        << " error in state var update for the force due to the displacement jump" );

  // test state var update for the surface strain
  Vector9d dSurfaceStressN = { 0.1, 0.2, 0.3, 0.4, 0.4, 0.6, 0.7, 0.8, 0.9 }; // arbitrary surface strain increment
  StateVarMatrixRs corrStateVarsRsN( 9, nMaxwellRs );
  corrStateVarsRsN << 1185.2297265979575513, 1580.3087546602175735, 0.0110363832351433, 1580.3161122490409980,
    2765.5458388469987767, 0.0220727664702865, 0.0257515608820010, 0.0294303552937154, 0.0331091497054298;

  updateStateVarMatrixRs( dT, elasticModuliRs, relaxationTimesRs, stateVarsRs, dSurfaceStressN, unitZ_voigt_full );

  throwExceptionOnFailure( checkIfEqual< double >( stateVarsRs, corrStateVarsRsN ),
                           MakeString() << __PRETTY_FUNCTION__ << " error in state var update for the surface strain" );
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