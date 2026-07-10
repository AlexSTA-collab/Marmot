#include "Marmot/MarmotExtendedInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <functional>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Vector9d  = Eigen::Matrix< double, 9, 1 >;
  using Matrix21d = Eigen::Matrix< double, 21, 21, Eigen::RowMajor >;

  struct ExtendedEvaluation {
    Eigen::Matrix< double, 21, 1 > response;
    Matrix21d                      tangent;
  };

  void makeExtendedKinematics( const Eigen::Matrix< double, 21, 1 >& generalizedIncrement,
                               double*                               dU,
                               double*                               dSurfaceStrain )
  {
    Eigen::Map< Eigen::Matrix< double, 6, 1 > >  dUMap( dU );
    Eigen::Map< Eigen::Matrix< double, 18, 1 > > dSurfaceStrainMap( dSurfaceStrain );

    dUMap.setZero();
    dUMap.segment< 3 >( 0 ) = generalizedIncrement.segment< 3 >( 0 );

    const Vector9d averageSurfaceGradient = generalizedIncrement.segment< 9 >( 3 );
    const Vector9d surfaceGradientJump    = generalizedIncrement.segment< 9 >( 12 );

    dSurfaceStrainMap.segment< 9 >( 0 ) = averageSurfaceGradient + 0.5 * surfaceGradientJump;
    dSurfaceStrainMap.segment< 9 >( 9 ) = averageSurfaceGradient - 0.5 * surfaceGradientJump;
  }

  Matrix21d packExtendedTangent( double tangentBlocks[9][81] )
  {
    Matrix21d tangent;
    tangent.setZero();

    tangent.block< 3, 3 >( 0, 0 )   = Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( tangentBlocks[0] );
    tangent.block< 3, 9 >( 0, 3 )   = Eigen::Map< Eigen::Matrix< double, 3, 9, Eigen::RowMajor > >( tangentBlocks[1] );
    tangent.block< 3, 9 >( 0, 12 )  = Eigen::Map< Eigen::Matrix< double, 3, 9, Eigen::RowMajor > >( tangentBlocks[2] );
    tangent.block< 9, 3 >( 3, 0 )   = Eigen::Map< Eigen::Matrix< double, 9, 3, Eigen::RowMajor > >( tangentBlocks[3] );
    tangent.block< 9, 9 >( 3, 3 )   = Eigen::Map< Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >( tangentBlocks[4] );
    tangent.block< 9, 9 >( 3, 12 )  = Eigen::Map< Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >( tangentBlocks[5] );
    tangent.block< 9, 3 >( 12, 0 )  = Eigen::Map< Eigen::Matrix< double, 9, 3, Eigen::RowMajor > >( tangentBlocks[6] );
    tangent.block< 9, 9 >( 12, 3 )  = Eigen::Map< Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >( tangentBlocks[7] );
    tangent.block< 9, 9 >( 12, 12 ) = Eigen::Map< Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >( tangentBlocks[8] );

    return tangent;
  }

  ExtendedEvaluation evaluateExtendedMaterial( const std::string&                    materialName,
                                               const double*                         extendedProperties,
                                               int                                   nExtendedProperties,
                                               const Eigen::Matrix< double, 21, 1 >& generalizedIncrement )
  {
    const double normal[3] = { 0., 0., 1. };

    MarmotExtendedInterfaceMaterialHypoElastic material( materialName, extendedProperties, nExtendedProperties, 1 );
    Eigen::VectorXd                            stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    double dU[6]              = { 0. };
    double dSurfaceStrain[18] = { 0. };
    makeExtendedKinematics( generalizedIncrement, dU, dSurfaceStrain );

    Eigen::Vector3d force                = Eigen::Vector3d::Zero();
    Vector9d        averageSurfaceStress = Vector9d::Zero();
    Vector9d        jumpSurfaceStress    = Vector9d::Zero();

    double                                               tangentBlocks[9][81] = {};
    MarmotExtendedInterfaceMaterialHypoElastic::State    state{ force.data(),
                                                             averageSurfaceStress.data(),
                                                             jumpSurfaceStress.data(),
                                                             stateVars.data() };
    MarmotExtendedInterfaceMaterialHypoElastic::Tangents tangentBlockViews{
      tangentBlocks[0],
      tangentBlocks[1],
      tangentBlocks[2],
      tangentBlocks[3],
      tangentBlocks[4],
      tangentBlocks[5],
      tangentBlocks[6],
      tangentBlocks[7],
      tangentBlocks[8],
    };
    MarmotExtendedInterfaceMaterialHypoElastic::Deformation   deformation{ dU, dSurfaceStrain, normal };
    MarmotExtendedInterfaceMaterialHypoElastic::TimeIncrement time{ 0., 1. };

    material.computeStress( state, tangentBlockViews, deformation, time );

    ExtendedEvaluation evaluation;
    evaluation.response.segment< 3 >( 0 )  = force;
    evaluation.response.segment< 9 >( 3 )  = averageSurfaceStress;
    evaluation.response.segment< 9 >( 12 ) = jumpSurfaceStress;
    evaluation.tangent                     = packExtendedTangent( tangentBlocks );
    return evaluation;
  }

  ExtendedEvaluation evaluateExtendedLinearElastic( const double*                         extendedProperties,
                                                    const Eigen::Matrix< double, 21, 1 >& generalizedIncrement )
  {
    return evaluateExtendedMaterial( "LINEARELASTIC", extendedProperties, 7, generalizedIncrement );
  }

  void testExtendedMaterialReducesToStandardInterfaceForEqualSides()
  {
    const double standardProperties[3] = { 1e5, 0.3, 0.01 };
    const double extendedProperties[7] = { 0.01, 2., 1e5, 0.3, 2., 1e5, 0.3 };
    const double normal[3]             = { 0., 0., 1. };

    MarmotInterfaceMaterialHypoElastic         standardMaterial( "LINEARELASTIC", standardProperties, 3, 1 );
    MarmotExtendedInterfaceMaterialHypoElastic extendedMaterial( "LINEARELASTIC", extendedProperties, 7, 1 );

    Eigen::VectorXd standardStateVars( standardMaterial.getNumberOfRequiredStateVars() );
    Eigen::VectorXd extendedStateVars( extendedMaterial.getNumberOfRequiredStateVars() );
    standardMaterial.initializeYourself( standardStateVars.data(), standardStateVars.size() );
    extendedMaterial.initializeYourself( extendedStateVars.data(), extendedStateVars.size() );

    Eigen::Vector3d               forceStandard             = Eigen::Vector3d::Zero();
    Eigen::Vector3d               forceExtended             = Eigen::Vector3d::Zero();
    Eigen::Matrix< double, 9, 1 > surfaceStressStandard     = Eigen::Matrix< double, 9, 1 >::Zero();
    Eigen::Matrix< double, 9, 1 > surfaceStressExtended     = Eigen::Matrix< double, 9, 1 >::Zero();
    Eigen::Matrix< double, 9, 1 > surfaceStressJumpExtended = Eigen::Matrix< double, 9, 1 >::Zero();

    const double dU[6] = { 0., 1e-4, 0., 0., 0., 0. };
    const double dSurfaceStrain[18] =
      { 0., 2e-4, 0., 2e-4, 0., 0., 0., 0., 0., 0., 2e-4, 0., 2e-4, 0., 0., 0., 0., 0. };

    double standardQ[9]  = { 0. };
    double standardZ[81] = { 0. };
    double standardH[27] = { 0. };
    double standardY[81] = { 0. };

    MarmotInterfaceMaterialHypoElastic::State         standardState{ forceStandard.data(),
                                                             surfaceStressStandard.data(),
                                                             standardStateVars.data() };
    MarmotInterfaceMaterialHypoElastic::Tangents      standardTangents{ standardQ, standardZ, standardH, standardY };
    MarmotInterfaceMaterialHypoElastic::Deformation   standardDeformation{ dU, dSurfaceStrain, normal };
    MarmotInterfaceMaterialHypoElastic::TimeIncrement standardTime{ 0., 1. };

    standardMaterial.computeStress( standardState, standardTangents, standardDeformation, standardTime );

    double                                               extendedTangents[9][81] = {};
    MarmotExtendedInterfaceMaterialHypoElastic::State    extendedState{ forceExtended.data(),
                                                                     surfaceStressExtended.data(),
                                                                     surfaceStressJumpExtended.data(),
                                                                     extendedStateVars.data() };
    MarmotExtendedInterfaceMaterialHypoElastic::Tangents extendedTangentBlocks{
      extendedTangents[0],
      extendedTangents[1],
      extendedTangents[2],
      extendedTangents[3],
      extendedTangents[4],
      extendedTangents[5],
      extendedTangents[6],
      extendedTangents[7],
      extendedTangents[8],
    };
    MarmotExtendedInterfaceMaterialHypoElastic::Deformation   extendedDeformation{ dU, dSurfaceStrain, normal };
    MarmotExtendedInterfaceMaterialHypoElastic::TimeIncrement extendedTime{ 0., 1. };

    extendedMaterial.computeStress( extendedState, extendedTangentBlocks, extendedDeformation, extendedTime );

    throwExceptionOnFailure( checkIfEqual< double >( forceExtended, forceStandard, 1e-8 ),
                             "Extended interface force does not reduce to standard interface force." );
    throwExceptionOnFailure( checkIfEqual< double >( surfaceStressExtended, surfaceStressStandard, 1e-8 ),
                             "Extended interface average surface stress does not reduce to standard surface stress." );
    throwExceptionOnFailure( checkIfEqual< double >( surfaceStressJumpExtended,
                                                     Eigen::Matrix< double, 9, 1 >::Zero(),
                                                     1e-10 ),
                             "Extended interface surface stress jump is not zero for equal top and bottom sides." );
  }

  void testImplicitExtendedTangentMatchesFiniteDifference()
  {
    const double extendedProperties[7] = { 0.02, 2., 8e4, 0.22, 2., 1.5e5, 0.31 };

    Eigen::Matrix< double, 21, 1 > generalizedIncrement;
    generalizedIncrement.setZero();
    generalizedIncrement.segment< 3 >( 0 ) << 1.0e-5, -2.0e-5, 3.0e-5;
    generalizedIncrement.segment< 9 >( 3 ) << 2.0e-4, 1.0e-4, 0.0, -1.0e-4, 5.0e-5, 0.0, 0.0, 0.0, 0.0;
    generalizedIncrement.segment< 9 >( 12 ) << 7.0e-5, -3.0e-5, 0.0, 4.0e-5, -2.0e-5, 0.0, 0.0, 0.0, 0.0;

    const auto baseEvaluation = evaluateExtendedLinearElastic( extendedProperties, generalizedIncrement );

    Matrix21d finiteDifferenceTangent;
    finiteDifferenceTangent.setZero();
    for ( int i = 0; i < 21; ++i ) {
      Eigen::Matrix< double, 21, 1 > perturbedIncrement = generalizedIncrement;
      const double                   perturbation       = 1e-8 * std::max( 1.0, std::abs( generalizedIncrement[i] ) );
      perturbedIncrement[i] += perturbation;

      const auto perturbedEvaluation   = evaluateExtendedLinearElastic( extendedProperties, perturbedIncrement );
      finiteDifferenceTangent.col( i ) = ( perturbedEvaluation.response - baseEvaluation.response ) / perturbation;
    }

    const double error = ( baseEvaluation.tangent - finiteDifferenceTangent ).norm();
    const double scale = std::max( 1.0, finiteDifferenceTangent.norm() );
    throwExceptionOnFailure( error / scale < 1e-6,
                             "Implicit extended tangent does not match finite-difference tangent." );
  }

  void testSolvedNormalGradientJumpEnforcesTractionEquilibrium()
  {
    const double extendedProperties[7] = { 0.03, 2., 7.5e4, 0.18, 2., 1.9e5, 0.34 };
    const double normal[3]             = { 0., 0., 1. };

    MarmotExtendedInterfaceMaterialHypoElastic material( "LINEARELASTIC", extendedProperties, 7, 1 );

    Eigen::VectorXd stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    Eigen::Matrix< double, 21, 1 > generalizedIncrement;
    generalizedIncrement.setZero();
    generalizedIncrement.segment< 3 >( 0 ) << 2.5e-5, -1.5e-5, 1.0e-5;
    generalizedIncrement.segment< 9 >( 3 ) << 1.5e-4, 7.0e-5, 0.0, -4.0e-5, 9.0e-5, 0.0, 0.0, 0.0, 0.0;
    generalizedIncrement.segment< 9 >( 12 ) << -9.0e-5, 5.0e-5, 0.0, 6.0e-5, -3.0e-5, 0.0, 0.0, 0.0, 0.0;

    double dU[6]              = { 0. };
    double dSurfaceStrain[18] = { 0. };
    makeExtendedKinematics( generalizedIncrement, dU, dSurfaceStrain );

    Eigen::Vector3d force                = Eigen::Vector3d::Zero();
    Vector9d        averageSurfaceStress = Vector9d::Zero();
    Vector9d        jumpSurfaceStress    = Vector9d::Zero();
    double          tangentBlocks[9][81] = {};

    MarmotExtendedInterfaceMaterialHypoElastic::State    state{ force.data(),
                                                             averageSurfaceStress.data(),
                                                             jumpSurfaceStress.data(),
                                                             stateVars.data() };
    MarmotExtendedInterfaceMaterialHypoElastic::Tangents tangentBlockViews{
      tangentBlocks[0],
      tangentBlocks[1],
      tangentBlocks[2],
      tangentBlocks[3],
      tangentBlocks[4],
      tangentBlocks[5],
      tangentBlocks[6],
      tangentBlocks[7],
      tangentBlocks[8],
    };
    MarmotExtendedInterfaceMaterialHypoElastic::Deformation   deformation{ dU, dSurfaceStrain, normal };
    MarmotExtendedInterfaceMaterialHypoElastic::TimeIncrement time{ 0., 1. };

    material.computeStress( state, tangentBlockViews, deformation, time );

    const auto             topStressView    = material.getStateView( "topStress", stateVars.data() );
    const auto             bottomStressView = material.getStateView( "bottomStress", stateVars.data() );
    const Marmot::Vector6d topStress        = Eigen::Map< const Marmot::Vector6d >( topStressView.stateLocation );
    const Marmot::Vector6d bottomStress     = Eigen::Map< const Marmot::Vector6d >( bottomStressView.stateLocation );

    const Eigen::Matrix3d topStressTensor    = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( topStress );
    const Eigen::Matrix3d bottomStressTensor = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( bottomStress );
    const Eigen::Vector3d normalVector( normal[0], normal[1], normal[2] );

    const Eigen::Vector3d tractionJump = ( topStressTensor - bottomStressTensor ) * normalVector;
    throwExceptionOnFailure( tractionJump.norm() < 1e-8,
                             "Solved normal-gradient jump does not enforce top/bottom traction equilibrium." );
  }

  void testSolvedNormalGradientJumpEnforcesNonlinearTractionEquilibrium()
  {
    const double singleMaterialProperties[8] = { 210000., 0.3, 0.02, 120., 2100., 20., 20., 2400. };
    const double normal[3]                   = { 0., 0., 1. };

    MarmotExtendedInterfaceMaterialHypoElastic material( "VONMISES", singleMaterialProperties, 8, 1 );

    Eigen::VectorXd stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    Eigen::Matrix< double, 21, 1 > generalizedIncrement;
    generalizedIncrement.setZero();
    generalizedIncrement.segment< 3 >( 0 ) << 1.0e-4, -2.0e-4, 1.5e-4;
    generalizedIncrement.segment< 9 >( 3 ) << 4.0e-3, 1.5e-3, 0.0, 1.0e-3, -2.0e-3, 0.0, 0.0, 0.0, -2.0e-3;
    generalizedIncrement.segment< 9 >( 12 ) << 5.0e-3, -1.0e-3, 0.0, 2.0e-3, -3.0e-3, 0.0, 0.0, 0.0, 1.0e-3;

    double dU[6]              = { 0. };
    double dSurfaceStrain[18] = { 0. };
    makeExtendedKinematics( generalizedIncrement, dU, dSurfaceStrain );

    Eigen::Vector3d force                = Eigen::Vector3d::Zero();
    Vector9d        averageSurfaceStress = Vector9d::Zero();
    Vector9d        jumpSurfaceStress    = Vector9d::Zero();
    double          tangentBlocks[9][81] = {};

    MarmotExtendedInterfaceMaterialHypoElastic::State    state{ force.data(),
                                                             averageSurfaceStress.data(),
                                                             jumpSurfaceStress.data(),
                                                             stateVars.data() };
    MarmotExtendedInterfaceMaterialHypoElastic::Tangents tangentBlockViews{
      tangentBlocks[0],
      tangentBlocks[1],
      tangentBlocks[2],
      tangentBlocks[3],
      tangentBlocks[4],
      tangentBlocks[5],
      tangentBlocks[6],
      tangentBlocks[7],
      tangentBlocks[8],
    };
    MarmotExtendedInterfaceMaterialHypoElastic::Deformation   deformation{ dU, dSurfaceStrain, normal };
    MarmotExtendedInterfaceMaterialHypoElastic::TimeIncrement time{ 0., 1. };

    material.computeStress( state, tangentBlockViews, deformation, time );

    const auto             topStressView    = material.getStateView( "topStress", stateVars.data() );
    const auto             bottomStressView = material.getStateView( "bottomStress", stateVars.data() );
    const Marmot::Vector6d topStress        = Eigen::Map< const Marmot::Vector6d >( topStressView.stateLocation );
    const Marmot::Vector6d bottomStress     = Eigen::Map< const Marmot::Vector6d >( bottomStressView.stateLocation );

    const Eigen::Matrix3d topStressTensor    = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( topStress );
    const Eigen::Matrix3d bottomStressTensor = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( bottomStress );
    const Eigen::Vector3d normalVector( normal[0], normal[1], normal[2] );

    const Eigen::Vector3d tractionJump = ( topStressTensor - bottomStressTensor ) * normalVector;
    throwExceptionOnFailure( tractionJump.norm() < 1e-7,
                             "Solved normal-gradient jump does not enforce nonlinear top/bottom traction "
                             "equilibrium." );
  }

  void testImplicitExtendedTangentMatchesFiniteDifferenceForPlasticMaterial()
  {
    const double singleMaterialProperties[8] = { 210000., 0.3, 0.02, 120., 2100., 20., 20., 2400. };

    Eigen::Matrix< double, 21, 1 > generalizedIncrement;
    generalizedIncrement.setZero();
    generalizedIncrement.segment< 3 >( 0 ) << 1.0e-4, -2.0e-4, 1.5e-4;
    generalizedIncrement.segment< 9 >( 3 ) << 4.0e-3, 1.5e-3, 0.0, 1.0e-3, -2.0e-3, 0.0, 0.0, 0.0, -2.0e-3;
    generalizedIncrement.segment< 9 >( 12 ) << 5.0e-3, -1.0e-3, 0.0, 2.0e-3, -3.0e-3, 0.0, 0.0, 0.0, 1.0e-3;

    const auto baseEvaluation = evaluateExtendedMaterial( "VONMISES",
                                                          singleMaterialProperties,
                                                          8,
                                                          generalizedIncrement );

    Matrix21d finiteDifferenceTangent;
    finiteDifferenceTangent.setZero();
    for ( int i = 0; i < 21; ++i ) {
      Eigen::Matrix< double, 21, 1 > perturbedIncrement = generalizedIncrement;
      const double                   perturbation       = 1e-8 * std::max( 1.0, std::abs( generalizedIncrement[i] ) );
      perturbedIncrement[i] += perturbation;

      const auto perturbedEvaluation   = evaluateExtendedMaterial( "VONMISES",
                                                                 singleMaterialProperties,
                                                                 8,
                                                                 perturbedIncrement );
      finiteDifferenceTangent.col( i ) = ( perturbedEvaluation.response - baseEvaluation.response ) / perturbation;
    }

    const double error = ( baseEvaluation.tangent - finiteDifferenceTangent ).norm();
    const double scale = std::max( 1.0, finiteDifferenceTangent.norm() );
    throwExceptionOnFailure( error / scale < 5e-4,
                             "Plastic implicit extended tangent does not match finite-difference tangent." );
  }

  void testSingleMaterialInputKeepsIndependentTopAndBottomState()
  {
    const double singleMaterialProperties[8] = { 210000., 0.3, 0.01, 200., 2100., 20., 20., 2400. };
    const double normal[3]                   = { 0., 0., 1. };

    MarmotExtendedInterfaceMaterialHypoElastic material( "VONMISES", singleMaterialProperties, 8, 1 );

    Eigen::VectorXd stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    Eigen::Vector3d force                = Eigen::Vector3d::Zero();
    Vector9d        averageSurfaceStress = Vector9d::Zero();
    Vector9d        jumpSurfaceStress    = Vector9d::Zero();

    double dU[6]              = { 0. };
    double dSurfaceStrain[18] = { 0. };
    dSurfaceStrain[0]         = 0.02;
    dSurfaceStrain[4]         = -0.01;

    double                                               tangentBlocks[9][81] = {};
    MarmotExtendedInterfaceMaterialHypoElastic::State    state{ force.data(),
                                                             averageSurfaceStress.data(),
                                                             jumpSurfaceStress.data(),
                                                             stateVars.data() };
    MarmotExtendedInterfaceMaterialHypoElastic::Tangents tangentBlockViews{
      tangentBlocks[0],
      tangentBlocks[1],
      tangentBlocks[2],
      tangentBlocks[3],
      tangentBlocks[4],
      tangentBlocks[5],
      tangentBlocks[6],
      tangentBlocks[7],
      tangentBlocks[8],
    };
    MarmotExtendedInterfaceMaterialHypoElastic::Deformation   deformation{ dU, dSurfaceStrain, normal };
    MarmotExtendedInterfaceMaterialHypoElastic::TimeIncrement time{ 0., 1. };

    material.computeStress( state, tangentBlockViews, deformation, time );

    const auto bottomMaterialState = material.getStateView( "bottomMaterialStateVars", stateVars.data() );
    const auto topMaterialState    = material.getStateView( "topMaterialStateVars", stateVars.data() );

    throwExceptionOnFailure( bottomMaterialState.stateSize == 1 && topMaterialState.stateSize == 1,
                             "Von Mises top/bottom material state blocks should each contain kappa only." );
    throwExceptionOnFailure( std::abs( topMaterialState.stateLocation[0] - bottomMaterialState.stateLocation[0] ) >
                               1e-12,
                             "Single-material extended interface did not accumulate independent top/bottom kappa." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = {
    testExtendedMaterialReducesToStandardInterfaceForEqualSides,
    testImplicitExtendedTangentMatchesFiniteDifference,
    testSolvedNormalGradientJumpEnforcesTractionEquilibrium,
    testSolvedNormalGradientJumpEnforcesNonlinearTractionEquilibrium,
    testImplicitExtendedTangentMatchesFiniteDifferenceForPlasticMaterial,
    testSingleMaterialInputKeepsIndependentTopAndBottomState,
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
