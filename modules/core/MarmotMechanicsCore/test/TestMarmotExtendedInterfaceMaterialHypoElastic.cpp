#include "Marmot/MarmotExtendedInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Vector9d          = Eigen::Matrix< double, 9, 1 >;
  using Matrix21d         = Eigen::Matrix< double, 21, 21, Eigen::RowMajor >;
  using Matrix3dRowMajor  = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix3x9RowMajor = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;

  struct ExtendedEvaluation {
    Eigen::Matrix< double, 21, 1 > response;
    Matrix21d                      tangent;
  };

  struct TractionEquilibriumEvaluation {
    Eigen::Vector3d  residual;
    Matrix3dRowMajor jacobian;
  };

  Matrix3dRowMajor vectorToTensor( const Vector9d& vector )
  {
    return Eigen::Map< const Matrix3dRowMajor >( vector.data() );
  }

  Marmot::Vector6d strainToVoigt( const Matrix3dRowMajor& displacementGradient )
  {
    const Eigen::Matrix3d strain = 0.5 * ( displacementGradient + displacementGradient.transpose() );
    return Marmot::ContinuumMechanics::VoigtNotation::strainToVoigt( strain );
  }

  Matrix3dRowMajor stressToTensor( const Marmot::Vector6d& stress )
  {
    Matrix3dRowMajor stressTensor;
    stressTensor = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( stress );
    return stressTensor;
  }

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

  TractionEquilibriumEvaluation evaluateTractionEquilibriumEquation(
    const std::string&                         materialName,
    const double*                              extendedProperties,
    int                                        nExtendedProperties,
    const Eigen::Vector3d&                     averageNormalGradient,
    const Eigen::Vector3d&                     normalGradientJump,
    const Vector9d&                            averageSurfaceGradient,
    const Vector9d&                            surfaceGradientJump,
    const Eigen::Vector3d&                     normal,
    const MarmotMaterialHypoElastic::timeInfo& timeInfo )
  {
    MarmotExtendedInterfaceMaterialHypoElastic material( materialName, extendedProperties, nExtendedProperties, 1 );
    Eigen::VectorXd                            stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    const Matrix3dRowMajor averageSurfaceGradientTensor = vectorToTensor( averageSurfaceGradient );
    const Matrix3dRowMajor surfaceGradientJumpTensor    = vectorToTensor( surfaceGradientJump );

    Matrix3dRowMajor topDisplacementGradient    = averageSurfaceGradientTensor + 0.5 * surfaceGradientJumpTensor;
    Matrix3dRowMajor bottomDisplacementGradient = averageSurfaceGradientTensor - 0.5 * surfaceGradientJumpTensor;
    topDisplacementGradient += ( averageNormalGradient + 0.5 * normalGradientJump ) * normal.transpose();
    bottomDisplacementGradient += ( averageNormalGradient - 0.5 * normalGradientJump ) * normal.transpose();

    Marmot::Vector6d topStress     = Marmot::Vector6d::Zero();
    Marmot::Vector6d bottomStress  = Marmot::Vector6d::Zero();
    Marmot::Matrix6d topTangent    = Marmot::Matrix6d::Zero();
    Marmot::Matrix6d bottomTangent = Marmot::Matrix6d::Zero();

    auto topStateView    = material.getStateView( "topMaterialStateVars", stateVars.data() );
    auto bottomStateView = material.getStateView( "bottomMaterialStateVars", stateVars.data() );

    MarmotMaterialHypoElastic::state3D topState{ topStress, 0.0, 0.0, topStateView.stateLocation };
    MarmotMaterialHypoElastic::state3D bottomState{ bottomStress, 0.0, 0.0, bottomStateView.stateLocation };

    material.getTopMaterial().computeStress( topState, topTangent, strainToVoigt( topDisplacementGradient ), timeInfo );
    material.getBottomMaterial().computeStress( bottomState,
                                                bottomTangent,
                                                strainToVoigt( bottomDisplacementGradient ),
                                                timeInfo );

    const auto normalTensor                         = Marmot::FastorStandardTensors::Tensor3d( normal.data() );
    const auto [topZ, topQTensor, topHTensor, topY] = Marmot::Materials::InterfaceMaterialHelperFunctions::
      calculateInterfaceMaterialParameters( normalTensor, topTangent );
    const auto [bottomZ, bottomQTensor, bottomHTensor, bottomY] = Marmot::Materials::InterfaceMaterialHelperFunctions::
      calculateInterfaceMaterialParameters( normalTensor, bottomTangent );

    (void)topZ;
    (void)bottomZ;
    (void)topHTensor;
    (void)bottomHTensor;
    (void)topY;
    (void)bottomY;

    TractionEquilibriumEvaluation evaluation;
    evaluation.residual = ( stressToTensor( topState.stress ) - stressToTensor( bottomState.stress ) ) * normal;
    evaluation.jacobian = 0.5 * ( Eigen::Map< const Matrix3dRowMajor >( topQTensor.data() ) +
                                  Eigen::Map< const Matrix3dRowMajor >( bottomQTensor.data() ) );
    return evaluation;
  }

  Matrix3dRowMajor computeExplicitTractionEquilibriumJacobian(
    const std::function< TractionEquilibriumEvaluation( const Eigen::Vector3d& ) >& evaluator,
    const Eigen::Vector3d&                                                          normalGradientJump,
    double                                                                          relativePerturbation )
  {
    Matrix3dRowMajor perturbationJacobian;
    perturbationJacobian.setZero();

    for ( int i = 0; i < 3; ++i ) {
      Eigen::Vector3d plusNormalGradientJump  = normalGradientJump;
      Eigen::Vector3d minusNormalGradientJump = normalGradientJump;
      const double    perturbation = relativePerturbation * std::max( 1.0, std::abs( normalGradientJump[i] ) );

      plusNormalGradientJump[i] += perturbation;
      minusNormalGradientJump[i] -= perturbation;

      const auto plusEvaluation  = evaluator( plusNormalGradientJump );
      const auto minusEvaluation = evaluator( minusNormalGradientJump );

      perturbationJacobian.col( i ) = ( plusEvaluation.residual - minusEvaluation.residual ) / ( 2.0 * perturbation );
    }

    return perturbationJacobian;
  }

  Matrix21d computeExplicitPerturbationTangent(
    const std::function< ExtendedEvaluation( const Eigen::Matrix< double, 21, 1 >& ) >& evaluator,
    const Eigen::Matrix< double, 21, 1 >&                                               generalizedIncrement,
    double                                                                              relativePerturbation )
  {
    Matrix21d perturbationTangent;
    perturbationTangent.setZero();

    for ( int i = 0; i < 21; ++i ) {
      Eigen::Matrix< double, 21, 1 > plusIncrement  = generalizedIncrement;
      Eigen::Matrix< double, 21, 1 > minusIncrement = generalizedIncrement;
      const double perturbation = relativePerturbation * std::max( 1.0, std::abs( generalizedIncrement[i] ) );

      plusIncrement[i] += perturbation;
      minusIncrement[i] -= perturbation;

      const auto plusEvaluation  = evaluator( plusIncrement );
      const auto minusEvaluation = evaluator( minusIncrement );

      perturbationTangent.col( i ) = ( plusEvaluation.response - minusEvaluation.response ) / ( 2.0 * perturbation );
    }

    return perturbationTangent;
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

  void testImplicitExtendedTangentMatchesExplicitPerturbationTangent()
  {
    const double extendedProperties[7] = { 0.02, 2., 8e4, 0.22, 2., 1.5e5, 0.31 };

    Eigen::Matrix< double, 21, 1 > generalizedIncrement;
    generalizedIncrement.setZero();
    generalizedIncrement.segment< 3 >( 0 ) << 1.0e-5, -2.0e-5, 3.0e-5;
    generalizedIncrement.segment< 9 >( 3 ) << 2.0e-4, 1.0e-4, 0.0, -1.0e-4, 5.0e-5, 0.0, 0.0, 0.0, 0.0;
    generalizedIncrement.segment< 9 >( 12 ) << 7.0e-5, -3.0e-5, 0.0, 4.0e-5, -2.0e-5, 0.0, 0.0, 0.0, 0.0;

    const auto baseEvaluation = evaluateExtendedLinearElastic( extendedProperties, generalizedIncrement );

    const auto perturbationTangent = computeExplicitPerturbationTangent(
      [&]( const Eigen::Matrix< double, 21, 1 >& increment ) {
        return evaluateExtendedLinearElastic( extendedProperties, increment );
      },
      generalizedIncrement,
      1e-7 );

    const double error = ( baseEvaluation.tangent - perturbationTangent ).norm();
    const double scale = std::max( 1.0, perturbationTangent.norm() );
    throwExceptionOnFailure( error / scale < 1e-7,
                             "Implicit extended tangent does not match central explicit perturbation tangent." );
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

  void testTractionEquilibriumNewtonJacobianMatchesExplicitPerturbation()
  {
    const double extendedProperties[7] = { 0.03, 2., 7.5e4, 0.18, 2., 1.9e5, 0.34 };

    const Eigen::Vector3d normal = Eigen::Vector3d( 0.2, -0.3, 1.0 ).normalized();
    const Eigen::Vector3d averageNormalGradient( 8.0e-4, -4.0e-4, 2.0e-4 );
    const Eigen::Vector3d normalGradientJump( 3.0e-4, 2.0e-4, -1.0e-4 );

    Vector9d averageSurfaceGradient;
    averageSurfaceGradient << 2.0e-4, 1.0e-4, -3.0e-5, -1.0e-4, 5.0e-5, 2.0e-5, 4.0e-5, -2.0e-5, 1.0e-4;

    Vector9d surfaceGradientJump;
    surfaceGradientJump << 7.0e-5, -3.0e-5, 2.0e-5, 4.0e-5, -2.0e-5, -1.0e-5, 3.0e-5, 1.0e-5, -4.0e-5;

    const MarmotMaterialHypoElastic::timeInfo timeInfo{ 1.0, 1.0 };

    const auto baseEvaluation = evaluateTractionEquilibriumEquation( "LINEARELASTIC",
                                                                     extendedProperties,
                                                                     7,
                                                                     averageNormalGradient,
                                                                     normalGradientJump,
                                                                     averageSurfaceGradient,
                                                                     surfaceGradientJump,
                                                                     normal,
                                                                     timeInfo );

    const auto perturbationJacobian = computeExplicitTractionEquilibriumJacobian(
      [&]( const Eigen::Vector3d& perturbedNormalGradientJump ) {
        return evaluateTractionEquilibriumEquation( "LINEARELASTIC",
                                                    extendedProperties,
                                                    7,
                                                    averageNormalGradient,
                                                    perturbedNormalGradientJump,
                                                    averageSurfaceGradient,
                                                    surfaceGradientJump,
                                                    normal,
                                                    timeInfo );
      },
      normalGradientJump,
      1e-7 );

    const double error = ( baseEvaluation.jacobian - perturbationJacobian ).norm();
    const double scale = std::max( 1.0, perturbationJacobian.norm() );
    throwExceptionOnFailure( error / scale < 1e-7,
                             "Traction-equilibrium Newton Jacobian does not match explicit perturbation." );
  }

  void testPlasticTractionEquilibriumNewtonJacobianMatchesExplicitPerturbation()
  {
    const double singleMaterialProperties[8] = { 210000., 0.3, 0.02, 120., 2100., 20., 20., 2400. };

    const Eigen::Vector3d normal = Eigen::Vector3d( -0.25, 0.15, 1.0 ).normalized();
    const Eigen::Vector3d averageNormalGradient( 3.0e-3, -2.0e-3, 1.0e-3 );
    const Eigen::Vector3d normalGradientJump( 2.0e-3, 1.5e-3, -1.0e-3 );

    Vector9d averageSurfaceGradient;
    averageSurfaceGradient << 4.0e-3, 1.5e-3, -4.0e-4, 1.0e-3, -2.0e-3, 3.0e-4, 2.0e-4, -3.0e-4, -2.0e-3;

    Vector9d surfaceGradientJump;
    surfaceGradientJump << 5.0e-3, -1.0e-3, 5.0e-4, 2.0e-3, -3.0e-3, -4.0e-4, 3.0e-4, 2.0e-4, 1.0e-3;

    const MarmotMaterialHypoElastic::timeInfo timeInfo{ 1.0, 1.0 };

    const auto baseEvaluation = evaluateTractionEquilibriumEquation( "VONMISES",
                                                                     singleMaterialProperties,
                                                                     8,
                                                                     averageNormalGradient,
                                                                     normalGradientJump,
                                                                     averageSurfaceGradient,
                                                                     surfaceGradientJump,
                                                                     normal,
                                                                     timeInfo );

    const auto perturbationJacobian = computeExplicitTractionEquilibriumJacobian(
      [&]( const Eigen::Vector3d& perturbedNormalGradientJump ) {
        return evaluateTractionEquilibriumEquation( "VONMISES",
                                                    singleMaterialProperties,
                                                    8,
                                                    averageNormalGradient,
                                                    perturbedNormalGradientJump,
                                                    averageSurfaceGradient,
                                                    surfaceGradientJump,
                                                    normal,
                                                    timeInfo );
      },
      normalGradientJump,
      1e-7 );

    const double error = ( baseEvaluation.jacobian - perturbationJacobian ).norm();
    const double scale = std::max( 1.0, perturbationJacobian.norm() );
    throwExceptionOnFailure( error / scale < 5e-4,
                             "Plastic traction-equilibrium Newton Jacobian does not match explicit perturbation." );
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

  void testPlasticImplicitExtendedTangentMatchesExplicitPerturbationTangent()
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

    const auto perturbationTangent = computeExplicitPerturbationTangent(
      [&]( const Eigen::Matrix< double, 21, 1 >& increment ) {
        return evaluateExtendedMaterial( "VONMISES", singleMaterialProperties, 8, increment );
      },
      generalizedIncrement,
      1e-7 );

    const double error = ( baseEvaluation.tangent - perturbationTangent ).norm();
    const double scale = std::max( 1.0, perturbationTangent.norm() );
    throwExceptionOnFailure( error / scale < 5e-4,
                             "Plastic implicit extended tangent does not match central explicit perturbation "
                             "tangent." );
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
    testImplicitExtendedTangentMatchesExplicitPerturbationTangent,
    testSolvedNormalGradientJumpEnforcesTractionEquilibrium,
    testSolvedNormalGradientJumpEnforcesNonlinearTractionEquilibrium,
    testTractionEquilibriumNewtonJacobianMatchesExplicitPerturbation,
    testPlasticTractionEquilibriumNewtonJacobianMatchesExplicitPerturbation,
    testImplicitExtendedTangentMatchesFiniteDifferenceForPlasticMaterial,
    testPlasticImplicitExtendedTangentMatchesExplicitPerturbationTangent,
    testSingleMaterialInputKeepsIndependentTopAndBottomState,
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
