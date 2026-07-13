#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotExtendedInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
#include <memory>
#include <sstream>
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
    double                         alpha;
    bool                           alphaEvolutionActive;
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

  ExtendedEvaluation evaluateExtendedMaterialWithState( const std::string&                    materialName,
                                                        const double*                         extendedProperties,
                                                        int                                   nExtendedProperties,
                                                        const Eigen::Matrix< double, 21, 1 >& generalizedIncrement,
                                                        Eigen::VectorXd&                      stateVars,
                                                        double                                timeOld = 0.,
                                                        double                                dT      = 1.,
                                                        bool forceAlphaEvolutionActive                = false )
  {
    const double normal[3] = { 0., 0., 1. };

    MarmotExtendedInterfaceMaterialHypoElastic material( materialName, extendedProperties, nExtendedProperties, 1 );
    if ( stateVars.size() == 0 ) {
      stateVars.resize( material.getNumberOfRequiredStateVars() );
      material.initializeYourself( stateVars.data(), stateVars.size() );
    }
    throwExceptionOnFailure( stateVars.size() == material.getNumberOfRequiredStateVars(),
                             "Unexpected extended-material state size." );
    if ( forceAlphaEvolutionActive )
      material.getStateView( "alphaEvolutionActive", stateVars.data() ).stateLocation[0] = 1.0;

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
    MarmotExtendedInterfaceMaterialHypoElastic::TimeIncrement time{ timeOld, dT };

    material.computeStress( state, tangentBlockViews, deformation, time );

    ExtendedEvaluation evaluation;
    evaluation.response.segment< 3 >( 0 )  = force;
    evaluation.response.segment< 9 >( 3 )  = averageSurfaceStress;
    evaluation.response.segment< 9 >( 12 ) = jumpSurfaceStress;
    evaluation.tangent                     = packExtendedTangent( tangentBlocks );
    evaluation.alpha                       = material.getStateView( "alpha", stateVars.data() ).stateLocation[0];
    evaluation.alphaEvolutionActive        = material.getStateView( "alphaEvolutionActive", stateVars.data() )
                                        .stateLocation[0] > 0.5;
    return evaluation;
  }

  ExtendedEvaluation evaluateExtendedMaterial( const std::string&                    materialName,
                                               const double*                         extendedProperties,
                                               int                                   nExtendedProperties,
                                               const Eigen::Matrix< double, 21, 1 >& generalizedIncrement )
  {
    Eigen::VectorXd stateVars;
    return evaluateExtendedMaterialWithState( materialName,
                                              extendedProperties,
                                              nExtendedProperties,
                                              generalizedIncrement,
                                              stateVars );
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

  // Mirrors testExtendedMaterialReducesToStandardInterfaceForEqualSides, but
  // constructs the extended material through the *legacy* single-material
  // property layout [E, nu, h] (nMaterialProperties == 3, so the explicit
  // [h,nBottom,...] branch can never be selected) instead of the explicit
  // [h,nBottom,...,nTop,...] layout. MarmotInterfaceMaterialHypoElastic
  // parses the very same [E, nu, h, remaining...] layout in its own
  // constructor, so constructing it from the identical property array gives
  // an independent, already-verified reference for "did the legacy branch
  // extract h and reassemble {E, nu} correctly". This also closes the gap
  // that the legacy fallback was previously exercised only with VONMISES.
  void testLegacySingleMaterialLayoutReducesToStandardInterfaceForEqualSides()
  {
    const double standardProperties[3] = { 1e5, 0.3, 0.01 };
    const double legacyProperties[3]   = { 1e5, 0.3, 0.01 };
    const double normal[3]             = { 0., 0., 1. };

    MarmotInterfaceMaterialHypoElastic         standardMaterial( "LINEARELASTIC", standardProperties, 3, 1 );
    MarmotExtendedInterfaceMaterialHypoElastic extendedMaterial( "LINEARELASTIC", legacyProperties, 3, 1 );

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
                             "Legacy-layout extended interface force does not reduce to standard interface force." );
    throwExceptionOnFailure( checkIfEqual< double >( surfaceStressExtended, surfaceStressStandard, 1e-8 ),
                             "Legacy-layout extended interface average surface stress does not reduce to standard "
                             "surface stress." );
    throwExceptionOnFailure( checkIfEqual< double >( surfaceStressJumpExtended,
                                                     Eigen::Matrix< double, 9, 1 >::Zero(),
                                                     1e-10 ),
                             "Legacy-layout extended interface surface stress jump is not zero for equal top and "
                             "bottom sides." );
  }

  // Same idea as the test above, but with a base material (VONMISES) whose
  // legacy layout carries "remainingBaseMaterialProperties" beyond E and nu,
  // so the [materialProperties+3, materialProperties+nMaterialProperties)
  // slice actually has content to get wrong. Loading is kept well below the
  // yield stress (200) and perfectly symmetric between the two faces, so the
  // expected response is exactly the elastic response of the single
  // MarmotInterfaceMaterialHypoElastic parsing the identical property array.
  void testLegacyVonMisesLayoutReducesToStandardInterfaceInElasticRegime()
  {
    const double legacyVonMisesProperties[8] = { 210000., 0.3, 0.01, 200., 2100., 20., 20., 2400. };
    const double normal[3]                   = { 0., 0., 1. };

    MarmotInterfaceMaterialHypoElastic         standardMaterial( "VONMISES", legacyVonMisesProperties, 8, 1 );
    MarmotExtendedInterfaceMaterialHypoElastic extendedMaterial( "VONMISES", legacyVonMisesProperties, 8, 1 );

    Eigen::VectorXd standardStateVars( standardMaterial.getNumberOfRequiredStateVars() );
    Eigen::VectorXd extendedStateVars( extendedMaterial.getNumberOfRequiredStateVars() );
    standardMaterial.initializeYourself( standardStateVars.data(), standardStateVars.size() );
    extendedMaterial.initializeYourself( extendedStateVars.data(), extendedStateVars.size() );

    Eigen::Vector3d               forceStandard             = Eigen::Vector3d::Zero();
    Eigen::Vector3d               forceExtended             = Eigen::Vector3d::Zero();
    Eigen::Matrix< double, 9, 1 > surfaceStressStandard     = Eigen::Matrix< double, 9, 1 >::Zero();
    Eigen::Matrix< double, 9, 1 > surfaceStressExtended     = Eigen::Matrix< double, 9, 1 >::Zero();
    Eigen::Matrix< double, 9, 1 > surfaceStressJumpExtended = Eigen::Matrix< double, 9, 1 >::Zero();

    // Two orders of magnitude smaller than the LINEARELASTIC equal-sides
    // reference test: unlike that test, h = 0.01 here divides a *normal*
    // displacement jump into an average normal *gradient* (dU[1]/h), which
    // otherwise blows up the local shear strain enough to yield immediately
    // at a yield stress of only 200. These amplitudes keep every stress
    // component comfortably below yield (order E*strain ~ 0.4-8).
    const double dU[6] = { 0., 1e-6, 0., 0., 0., 0. };
    const double dSurfaceStrain[18] =
      { 0., 2e-6, 0., 2e-6, 0., 0., 0., 0., 0., 0., 2e-6, 0., 2e-6, 0., 0., 0., 0., 0. };

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
                             "Legacy VONMISES-layout extended interface force does not reduce to standard "
                             "interface force in the elastic regime." );
    throwExceptionOnFailure( checkIfEqual< double >( surfaceStressExtended, surfaceStressStandard, 1e-8 ),
                             "Legacy VONMISES-layout extended interface average surface stress does not reduce to "
                             "standard surface stress in the elastic regime." );
    throwExceptionOnFailure( checkIfEqual< double >( surfaceStressJumpExtended,
                                                     Eigen::Matrix< double, 9, 1 >::Zero(),
                                                     1e-10 ),
                             "Legacy VONMISES-layout extended interface surface stress jump is not zero for equal "
                             "top and bottom sides." );

    const auto bottomMaterialState = extendedMaterial.getStateView( "bottomMaterialStateVars",
                                                                    extendedStateVars.data() );
    const auto topMaterialState    = extendedMaterial.getStateView( "topMaterialStateVars", extendedStateVars.data() );
    throwExceptionOnFailure( std::abs( bottomMaterialState.stateLocation[0] ) < 1e-12 &&
                               std::abs( topMaterialState.stateLocation[0] ) < 1e-12,
                             "Elastic-regime legacy VONMISES loading unexpectedly accumulated plastic history." );
  }

  // Documents a genuine, currently-unresolved ambiguity in the layout
  // detection: hasExplicitTopBottomLayout is chosen whenever
  // nMaterialProperties >= 5 and materialProperties[1] happens to be
  // integer-valued. nu == 0.0 is a legitimate (if unusual) Poisson's ratio,
  // and any legacy-layout call whose base material needs 2 or more of its
  // own properties (e.g. VONMISES, which needs 5) reaches nMaterialProperties
  // >= 5. Such a call is therefore misclassified as the explicit
  // [h,nBottom,bottom...,nTop,top...] layout instead of the intended legacy
  // [E,nu,h,remaining...] layout. This does not silently produce a
  // wrong-but-plausible result: materialProperties[1] (nu=0.0) is read as
  // nBottom, rounds to 0, and the explicit-layout validation immediately
  // rejects nBottom <= 0. This test locks in that current (accepted)
  // behavior -- a clear std::invalid_argument at construction time -- as a
  // regression guard, and as documentation that nu == 0.0 is not usable via
  // the legacy layout once the base material contributes >= 2 properties of
  // its own. Callers who need nu == 0.0 with such a base material must use
  // the explicit [h,nBottom,...,nTop,...] layout instead.
  void testLegacyLayoutWithIntegerValuedNuIsMisclassifiedAndRejected()
  {
    const double ambiguousProperties[8] = { 210000., 0.0, 0.01, 200., 2100., 20., 20., 2400. };

    bool        threw = false;
    std::string exceptionMessage;
    try {
      MarmotExtendedInterfaceMaterialHypoElastic material( "VONMISES", ambiguousProperties, 8, 1 );
      (void)material;
    }
    catch ( const std::invalid_argument& e ) {
      threw            = true;
      exceptionMessage = e.what();
    }

    throwExceptionOnFailure( threw,
                             "A legacy-layout call with nu == 0.0 and >= 5 total properties was expected to be "
                             "misclassified as the explicit layout and rejected, but construction succeeded "
                             "instead; if the layout-detection heuristic changed, please update this regression "
                             "test to match the new documented behavior." );
    throwExceptionOnFailure( exceptionMessage.find( "Invalid extended interface material layout" ) != std::string::npos,
                             "Unexpected exception message for the nu==0.0 layout-ambiguity case: " +
                               exceptionMessage );
  }

  void testIndeterminateAlphaWithNonzeroNormalGradientJumpDoesNotFail()
  {
    const double singleMaterialProperties[8] = { 4e5, 0.3, 0.01, 5., 0.1, 0., 0., 0. };

    Eigen::Matrix< double, 21, 1 > generalizedIncrement;
    for ( int i = 0; i < generalizedIncrement.size(); ++i )
      generalizedIncrement[i] = 1e-10 * std::sin( 0.37 * ( i + 2 ) );

    const auto evaluation = evaluateExtendedMaterial( "VONMISES", singleMaterialProperties, 8, generalizedIncrement );
    throwExceptionOnFailure( evaluation.response.allFinite(), "Indeterminate-alpha response contains nan or inf." );
    throwExceptionOnFailure( evaluation.tangent.allFinite(), "Indeterminate-alpha tangent contains nan or inf." );
  }

  void testAlphaRemainsFixedDuringElasticLoading()
  {
    const double singleMaterialProperties[8] = { 4e5, 0.3, 0.01, 5., 0.1, 0., 0., 0. };

    Eigen::Matrix< double, 21, 1 > generalizedIncrement;
    for ( int i = 0; i < generalizedIncrement.size(); ++i )
      generalizedIncrement[i] = 1e-7 * std::sin( 0.41 * ( i + 1 ) );

    const auto evaluation = evaluateExtendedMaterial( "VONMISES", singleMaterialProperties, 8, generalizedIncrement );
    throwExceptionOnFailure( std::abs( evaluation.alpha - 0.5 ) < 1e-14,
                             "Alpha changed although no plastic history variable evolved." );
  }

  void testCommittedPlasticActivityControlsAlphaEvolution()
  {
    const double singleMaterialProperties[8] = { 210000., 0.3, 0.02, 120., 2100., 20., 20., 2400. };

    Eigen::Matrix< double, 21, 1 > plasticIncrement;
    plasticIncrement.setZero();
    plasticIncrement.segment< 3 >( 0 ) << 1.0e-4, -2.0e-4, 1.5e-4;
    plasticIncrement.segment< 9 >( 3 ) << 4.0e-3, 1.5e-3, 0.0, 1.0e-3, -2.0e-3, 0.0, 0.0, 0.0, -2.0e-3;
    plasticIncrement.segment< 9 >( 12 ) << 5.0e-3, -1.0e-3, 0.0, 2.0e-3, -3.0e-3, 0.0, 0.0, 0.0, 1.0e-3;

    Eigen::VectorXd stateVars;
    const auto      firstYielding = evaluateExtendedMaterialWithState( "VONMISES",
                                                                  singleMaterialProperties,
                                                                  8,
                                                                  plasticIncrement,
                                                                  stateVars,
                                                                  0.,
                                                                  1. );
    throwExceptionOnFailure( std::abs( firstYielding.alpha - 0.5 ) < 1e-14,
                             "Alpha moved during the first yielding increment." );
    throwExceptionOnFailure( firstYielding.alphaEvolutionActive,
                             "Plastic flow did not activate alpha evolution for the next increment." );

    const auto plasticContinuation = evaluateExtendedMaterialWithState( "VONMISES",
                                                                        singleMaterialProperties,
                                                                        8,
                                                                        plasticIncrement,
                                                                        stateVars,
                                                                        1.,
                                                                        1. );
    throwExceptionOnFailure( std::abs( plasticContinuation.alpha - firstYielding.alpha ) > 1e-8,
                             "Committed plastic activity did not enable alpha evolution." );
    throwExceptionOnFailure( plasticContinuation.alphaEvolutionActive,
                             "Continued plastic flow unexpectedly deactivated alpha evolution." );

    const Eigen::Matrix< double, 21, 1 > zeroIncrement        = Eigen::Matrix< double, 21, 1 >::Zero();
    const double                         alphaBeforeUnloading = plasticContinuation.alpha;
    const auto                           elasticUnloading     = evaluateExtendedMaterialWithState( "VONMISES",
                                                                     singleMaterialProperties,
                                                                     8,
                                                                     zeroIncrement,
                                                                     stateVars,
                                                                     2.,
                                                                     1. );
    throwExceptionOnFailure( !elasticUnloading.alphaEvolutionActive,
                             "Elastic unloading did not deactivate alpha evolution for the next increment." );

    const auto followingElasticIncrement = evaluateExtendedMaterialWithState( "VONMISES",
                                                                              singleMaterialProperties,
                                                                              8,
                                                                              zeroIncrement,
                                                                              stateVars,
                                                                              3.,
                                                                              1. );
    throwExceptionOnFailure( std::abs( followingElasticIncrement.alpha - elasticUnloading.alpha ) < 1e-14,
                             "Alpha moved after committed plastic activity was deactivated." );
    throwExceptionOnFailure( std::isfinite( alphaBeforeUnloading ), "Plastic alpha is not finite." );

    const Eigen::VectorXd committedState = stateVars;
    Eigen::VectorXd       repeatedStateA = committedState;
    Eigen::VectorXd       repeatedStateB = committedState;
    const auto            repeatedA      = evaluateExtendedMaterialWithState( "VONMISES",
                                                              singleMaterialProperties,
                                                              8,
                                                              zeroIncrement,
                                                              repeatedStateA,
                                                              4.,
                                                              1. );
    const auto            repeatedB      = evaluateExtendedMaterialWithState( "VONMISES",
                                                              singleMaterialProperties,
                                                              8,
                                                              zeroIncrement,
                                                              repeatedStateB,
                                                              4.,
                                                              1. );
    throwExceptionOnFailure( ( repeatedA.response - repeatedB.response ).norm() < 1e-14 &&
                               ( repeatedA.tangent - repeatedB.tangent ).norm() < 1e-12 &&
                               ( repeatedStateA - repeatedStateB ).norm() < 1e-14,
                             "Repeated trials from the same committed state are not deterministic." );
  }

  void testActiveAlphaSatisfiesLowerAndUpperBoundKKTConditions()
  {
    constexpr double alphaMinimum            = 1e-4;
    const double     lowerBoundProperties[7] = { 0.02, 2., 8e4, 0.25, 2., 1.6e5, 0.25 };
    const double     upperBoundProperties[7] = { 0.02, 2., 1.6e5, 0.25, 2., 8e4, 0.25 };

    Eigen::Matrix< double, 21, 1 > generalizedIncrement = Eigen::Matrix< double, 21, 1 >::Zero();
    generalizedIncrement[2]                             = 2e-5;

    Eigen::VectorXd lowerState;
    setenv( "MARMOT_EI_VALIDATE_LOCAL_DERIVATIVES", "1", 1 );
    const auto lower = evaluateExtendedMaterialWithState( "LINEARELASTIC",
                                                          lowerBoundProperties,
                                                          7,
                                                          generalizedIncrement,
                                                          lowerState,
                                                          0.,
                                                          1.,
                                                          true );
    unsetenv( "MARMOT_EI_VALIDATE_LOCAL_DERIVATIVES" );
    throwExceptionOnFailure( std::abs( lower.alpha - alphaMinimum ) < 1e-10,
                             "Active alpha did not lock at the lower KKT bound." );
    throwExceptionOnFailure( lower.response.allFinite() && lower.tangent.allFinite(),
                             "Lower-bound active-set response or tangent contains nan or inf." );

    Eigen::VectorXd upperState;
    const auto      upper = evaluateExtendedMaterialWithState( "LINEARELASTIC",
                                                          upperBoundProperties,
                                                          7,
                                                          generalizedIncrement,
                                                          upperState,
                                                          0.,
                                                          1.,
                                                          true );
    throwExceptionOnFailure( std::abs( upper.alpha - ( 1.0 - alphaMinimum ) ) < 1e-10,
                             "Active alpha did not lock at the upper KKT bound." );
    throwExceptionOnFailure( upper.response.allFinite() && upper.tangent.allFinite(),
                             "Upper-bound active-set response or tangent contains nan or inf." );
  }

  // A perfectly-plastic Von Mises sublayer (yield stress essentially zero,
  // no hardening/softening) makes the elastic-branch acoustic Jacobian
  // (1-alpha)*Qtop + alpha*Qbottom singular for *any* nonzero kinematic
  // increment: the very first trial already sits fully inside the singular
  // return-mapped tangent, independent of step size. This is therefore not a
  // "too large a step" failure that a cutback could ever repair, but it is a
  // genuine, naturally reachable Marmot::StressUpdateFailed raised from
  // solveNormalGradientJumpForAlpha (i.e. reached before the commit=true
  // MaterialTrial evaluation in computeStress runs). It is used here purely
  // to exercise the "state must not be mutated by a failed trial" invariant
  // documented in evaluateSideTrial/computeMaterialTrial (only commit==true
  // callers write into the real bottomStress/topStress/*MaterialStateVars
  // pointers; every trial evaluation used while iterating operates on a
  // local stateCopy and is discarded on failure).
  void testNaturalStressUpdateFailurePreservesCommittedState()
  {
    const double singleMaterialProperties[8] = { 210000., 0.3, 1e-8, 0., 0., 0., 0., 2400. };

    MarmotExtendedInterfaceMaterialHypoElastic material( "VONMISES", singleMaterialProperties, 8, 1 );

    Eigen::VectorXd stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );
    const Eigen::VectorXd stateBeforeFailedAttempt = stateVars;

    Eigen::Matrix< double, 21, 1 > generalizedIncrement;
    generalizedIncrement.setZero();
    generalizedIncrement.segment< 3 >( 0 ) << 1.0e-4, -2.0e-4, 1.5e-4;
    generalizedIncrement.segment< 9 >( 3 ) << 4.0e-3, 1.5e-3, 0.0, 1.0e-3, -2.0e-3, 0.0, 0.0, 0.0, -2.0e-3;
    generalizedIncrement.segment< 9 >( 12 ) << 5.0e-3, -1.0e-3, 0.0, 2.0e-3, -3.0e-3, 0.0, 0.0, 0.0, 1.0e-3;

    double dU[6]              = { 0. };
    double dSurfaceStrain[18] = { 0. };
    makeExtendedKinematics( generalizedIncrement, dU, dSurfaceStrain );
    const double normal[3] = { 0., 0., 1. };

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

    bool threw = false;
    try {
      material.computeStress( state, tangentBlockViews, deformation, time );
    }
    catch ( const Marmot::StressUpdateFailed& ) {
      threw = true;
    }

    throwExceptionOnFailure( threw,
                             "Perfectly-plastic single-material extended interface did not raise "
                             "StressUpdateFailed for a nonzero increment as expected." );
    throwExceptionOnFailure( ( stateVars - stateBeforeFailedAttempt ).lpNorm< Eigen::Infinity >() == 0.0,
                             "Committed state vars were mutated by a trial that raised StressUpdateFailed; "
                             "a failed increment must leave the committed state byte-identical so that a "
                             "time-step cutback can safely retry from the same starting point." );
  }

  // Models the actual cutback protocol an external FE driver implements
  // around MarmotElement::computeYourself / MarmotExtendedInterfaceMaterialHypoElastic:
  // an increment attempt is evaluated on a private copy of the committed
  // state; if it is rejected (here simulated by discarding the mutated copy
  // instead of persisting it, exactly as a driver would do after catching
  // Marmot::StressUpdateFailed and requesting pNewDT<1), the *real* committed
  // stateVars snapshot is untouched and the retried increment is computed
  // fresh from that same snapshot. This verifies the retried increment is
  // deterministic and depends only on the pre-attempt committed snapshot,
  // not on whatever the rejected trial happened to compute.
  void testCutbackRetryFromPreservedSnapshotIsDeterministicAfterDiscardedAttempt()
  {
    const double singleMaterialProperties[8] = { 210000., 0.3, 0.02, 120., 2100., 20., 20., 2400. };

    Eigen::Matrix< double, 21, 1 > plasticIncrement;
    plasticIncrement.setZero();
    plasticIncrement.segment< 3 >( 0 ) << 1.0e-4, -2.0e-4, 1.5e-4;
    plasticIncrement.segment< 9 >( 3 ) << 4.0e-3, 1.5e-3, 0.0, 1.0e-3, -2.0e-3, 0.0, 0.0, 0.0, -2.0e-3;
    plasticIncrement.segment< 9 >( 12 ) << 5.0e-3, -1.0e-3, 0.0, 2.0e-3, -3.0e-3, 0.0, 0.0, 0.0, 1.0e-3;

    // Build up some real plastic history so the committed snapshot we cut
    // back from is nontrivial (nonzero stresses, evolved kappa, alpha
    // evolution active) rather than the pristine initial state.
    Eigen::VectorXd stateVars;
    const auto      firstYielding = evaluateExtendedMaterialWithState( "VONMISES",
                                                                  singleMaterialProperties,
                                                                  8,
                                                                  plasticIncrement,
                                                                  stateVars,
                                                                  0.,
                                                                  1. );
    throwExceptionOnFailure( firstYielding.alphaEvolutionActive,
                             "Precondition failed: expected plastic flow to activate alpha evolution." );

    const Eigen::VectorXd preAttemptCommittedSnapshot = stateVars;

    // "Attempt" the next increment at the full, dT=1 load and then discard
    // the result entirely -- as an external driver would after receiving a
    // pNewDT<1 cutback request -- instead of persisting it back into the
    // committed buffer.
    Eigen::VectorXd discardedAttemptState = preAttemptCommittedSnapshot;
    const auto      discardedAttempt      = evaluateExtendedMaterialWithState( "VONMISES",
                                                                     singleMaterialProperties,
                                                                     8,
                                                                     plasticIncrement,
                                                                     discardedAttemptState,
                                                                     1.,
                                                                     1. );
    (void)discardedAttempt;

    // Retry with a smaller dT (cutback), starting again from the untouched
    // pre-attempt committed snapshot, exactly as if the discarded attempt
    // above had never happened.
    const double    cutbackDT   = 0.25;
    Eigen::VectorXd retryStateA = preAttemptCommittedSnapshot;
    const auto      retryA      = evaluateExtendedMaterialWithState( "VONMISES",
                                                           singleMaterialProperties,
                                                           8,
                                                           plasticIncrement,
                                                           retryStateA,
                                                           1.,
                                                           cutbackDT );

    // An independent second retry from a fresh copy of the same pre-attempt
    // snapshot must reproduce the first retry bit-for-bit: the outcome must
    // depend only on the preserved committed snapshot and the retried
    // increment/dT, never on the discarded attempt or on evaluation order.
    Eigen::VectorXd retryStateB = preAttemptCommittedSnapshot;
    const auto      retryB      = evaluateExtendedMaterialWithState( "VONMISES",
                                                           singleMaterialProperties,
                                                           8,
                                                           plasticIncrement,
                                                           retryStateB,
                                                           1.,
                                                           cutbackDT );

    throwExceptionOnFailure( retryA.response.allFinite() && retryA.tangent.allFinite(),
                             "Cutback retry response or tangent contains nan or inf." );
    throwExceptionOnFailure( ( retryA.response - retryB.response ).norm() < 1e-14 &&
                               ( retryA.tangent - retryB.tangent ).norm() < 1e-12 &&
                               ( retryStateA - retryStateB ).norm() < 1e-14,
                             "Cutback retries from an identical preserved committed snapshot are not "
                             "deterministic, or leaked state from the discarded failed-dT=1 attempt." );

    // The preserved snapshot itself must still be exactly what it was before
    // either the discarded attempt or the retries ran, since evaluateExtendedMaterialWithState
    // always operates on the copy it is given rather than mutating shared state.
    throwExceptionOnFailure( ( preAttemptCommittedSnapshot - stateVars ).lpNorm< Eigen::Infinity >() == 0.0,
                             "Pre-attempt committed snapshot was unexpectedly modified." );
  }

  // ===========================================================================
  // Independent closed-form reference for a *moving* interior kink position.
  //
  // Background / why the existing alpha tests are not enough (see task
  // background): testAlphaRemainsFixedDuringElasticLoading and
  // testCommittedPlasticActivityControlsAlphaEvolution only check that alpha
  // stays put or moves *somewhere*, and
  // testActiveAlphaSatisfiesLowerAndUpperBoundKKTConditions only checks that
  // alpha saturates at one of the KKT bounds alphaMinimum / 1-alphaMinimum.
  // None of them pins alpha to a specific interior number derived by a route
  // that does not itself re-run the condensation code.
  //
  // Derivation. Restrict to a *pure tangential* (mode-II) interface slip:
  // only the x-component of the displacement jump dU is nonzero, everything
  // else (surface gradients) is zero, normal = (0,0,1). Then, by isotropy,
  // only the x-component of the local unknowns is active, and the whole
  // condensed problem collapses to a single scalar "two nonlinear springs in
  // series" problem:
  //
  //   gammaTop    = ubar + (1-alpha)*g
  //   gammaBottom = ubar - alpha*g
  //
  // (ubar = averageNormalGradient_x = dU_x/h, g = normalGradientJump_x),
  // which by construction always satisfies the *identity*
  //   alpha*gammaTop + (1-alpha)*gammaBottom = ubar                      (K)
  // for *every* alpha, g -- this is pure kinematics, not an equation to
  // solve; it comes directly from MarmotExtendedInterfaceMaterialHypoElastic.cpp's
  // reconstructFaceGradients().
  //
  // Each sublayer is an elastic-linear-hardening Von Mises material loaded in
  // simple shear. For an isotropic material with only one nonzero
  // displacement-gradient component, the deviatoric stress is exactly the
  // shear stress tau = G*gamma_elastic (elastic law, always), and the
  // classical J2 conversion sigma_eq = sqrt(3)*tau together with the
  // standard equivalent-plastic-strain identity dKappa = dGammaPlastic/sqrt(3)
  // for pure-shear associative flow (verified against VonMisesModel's own
  // return map below) turns the isotropic linear hardening law
  // fy(kappa) = sigmaY + HLin*kappa into the *plastic* hardening law
  // tau = tauY + Hs*gammaPlastic with tauY = sigmaY/sqrt(3), Hs = HLin/3.
  // Combined with the *elastic* law gammaElastic = tau/G and the additive
  // split gamma = gammaElastic + gammaPlastic, eliminating gammaPlastic gives
  // the *apparent* (tangent) traction-slip law actually observed in total
  // strain gamma:
  //
  //   tau(gamma) = G*gamma                              , gamma <= gammaY
  //   tau(gamma) = tauY + Htan*(gamma-gammaY)            , gamma >  gammaY
  //
  // with gammaY = tauY/G and Htan = G*Hs/(G+Hs) -- the elastic and plastic
  // branches combine like two compliances *in series* (1/Htan = 1/G + 1/Hs),
  // NOT Htan = Hs. (This subtlety was caught and fixed by directly probing
  // VonMisesModel's stress-vs-kappa response for a pure-shear increment
  // outside this repo: the naive Htan = Hs guess reproduced the elastic
  // branch exactly but was off by a factor of Hs/(G+Hs) on the hardening
  // branch -- e.g. ~25% low for Hs = HLin/3 comparable to G -- which is
  // exactly explained by gammaElastic continuing to grow with tau past first
  // yield instead of freezing at gammaY.)
  //
  // Because the whole load is applied in one single monotonic increment
  // starting from the virgin (zero stress/zero kappa) state, VonMisesModel's
  // own incrementalPotential for that increment is exactly the area under
  // this tau(gamma) curve from 0 to gamma (elastic triangle + hardening
  // trapezoid):
  //   Psi(gamma) = 0.5*G*gammaY^2 + tauY*(gamma-gammaY) + 0.5*Htan*(gamma-gammaY)^2
  //   (for gamma > gammaY; Psi(gamma)=0.5*G*gamma^2 elastically).
  //
  // computeMaterialTrial()/solveCoupledLocalProblem() in
  // MarmotExtendedInterfaceMaterialHypoElastic.cpp define (see the g- and
  // alpha-stationarity residuals assembled in computeLocalResidual /
  // MaterialTrial::alphaResidual):
  //   dPi/dg     = alpha*(1-alpha)*(tauTop - tauBottom)                  (traction continuity)
  //   dPi/dalpha = PsiTop - PsiBottom - (alpha*tauTop+(1-alpha)*tauBottom)*g
  // where Pi = alpha*PsiTop + (1-alpha)*PsiBottom is exactly
  // MaterialTrial::reducedPotential. At a converged, non-active-bound local
  // solution both vanish. Traction continuity forces tauTop=tauBottom=:t, and
  // then the weighted traction in the alpha-residual collapses to plain t,
  // so:  dPi/dalpha = PsiTop(t) - PsiBottom(t) - t*g(t).
  //
  // For *both* sublayers on their hardening branch (gamma_i > gammaY_i),
  // write gamma_i(t) = gammaY_i + (t-tauY_i)/Htan_i =: a_i + b_i*t with
  // b_i=1/Htan_i, a_i=tauY_i*(1/G-b_i), and (elementary algebra, expand the
  // trapezoid formula in t):
  //   Psi_i(t) = 0.5*b_i*t^2 + d_i ,  d_i := 0.5*tauY_i^2*(1/G-b_i).
  // g(t) = gammaTop(t)-gammaBottom(t) = (a_top-a_bottom) + (b_top-b_bottom)*t.
  // Substituting into dPi/dalpha=0 and collecting in t gives a *materials-only*
  // quadratic (independent of alpha and of the applied ubar!):
  //   -0.5*deltaB*t^2 - deltaA*t + deltaD = 0                            (S)
  // with deltaB=b_top-b_bottom, deltaA=a_top-a_bottom, deltaD=d_top-d_bottom.
  // (Sanity check performed while deriving this: for tauY_top=tauY_bottom
  // equation (S) collapses to -0.5*deltaB*(t-tauY)^2=0, a repeated root --
  // exactly reproducing the *monotonic*, boundary-seeking behaviour already
  // covered by testActiveAlphaSatisfiesLowerAndUpperBoundKKTConditions. A
  // genuine sign change -- and hence an interior stationary point -- only
  // appears once tauY_top != tauY_bottom, i.e. once the two sublayers have
  // different yield stresses in addition to different hardening moduli.)
  //
  // Any root t* of (S) with t* > max(tauY_top,tauY_bottom) is a traction at
  // which both branches are simultaneously admissible; the corresponding
  // strains gammaTop* = gammaTop(t*), gammaBottom* = gammaBottom(t*) are then
  // fixed, *alpha-independent* numbers, and by (K) the interior alpha that
  // reaches t* under a given ubar is simply
  //   alpha* = (ubar - gammaBottom*) / (gammaTop* - gammaBottom*).        (A)
  // Conversely, picking any target alpha* in (0,1) and setting
  //   ubar = alpha* * gammaTop* + (1-alpha*) * gammaBottom*
  // guarantees the coupled local problem's stationary point sits exactly at
  // that alpha*. This whole derivation is plain scalar algebra on the
  // sublayers' own bilinear traction laws -- it never calls
  // MarmotExtendedInterfaceMaterialHypoElastic's g/alpha Newton solve, so it
  // is a genuinely independent route to the same number.
  //
  // This was checked independently offline (not part of this repo) by
  // directly root-finding the *original* piecewise system (bisecting the
  // traction t for each alpha off a grid, no use of equation (S)) and
  // confirming the resulting profile Pi(alpha) has zero central-difference
  // derivative and strictly positive curvature exactly at the alpha* used
  // below, for every case tested -- i.e. equation (S)/(A) was cross-checked
  // against a brute-force reference before being encoded here.
  struct BilinearKinkPrediction {
    double alphaTarget;
    double tangentialJump;
    double criticalTraction;
    double gammaTopAtCritical;
    double gammaBottomAtCritical;
  };

  // Implements only the closed-form arithmetic of the derivation above.
  // Deliberately does not touch MarmotExtendedInterfaceMaterialHypoElastic,
  // MarmotMaterialHypoElastic, or VonMisesModel in any way.
  BilinearKinkPrediction predictBilinearKinkAlpha( double E,
                                                   double nu,
                                                   double sigmaYTop,
                                                   double HLinTop,
                                                   double sigmaYBottom,
                                                   double HLinBottom,
                                                   double h,
                                                   double alphaTarget )
  {
    const double sqrt3 = std::sqrt( 3.0 );
    const double G     = E / ( 2.0 * ( 1.0 + nu ) );

    const double tauYTop    = sigmaYTop / sqrt3;
    const double tauYBottom = sigmaYBottom / sqrt3;
    const double HsTop      = HLinTop / 3.0;
    const double HsBottom   = HLinBottom / 3.0;
    // Apparent (tangent) hardening modulus in total-strain space: elastic
    // and plastic compliances add in series, Htan = 1/(1/G+1/Hs).
    const double HtanTop    = G * HsTop / ( G + HsTop );
    const double HtanBottom = G * HsBottom / ( G + HsBottom );

    const double bTop    = 1.0 / HtanTop;
    const double bBottom = 1.0 / HtanBottom;
    const double aTop    = tauYTop * ( 1.0 / G - bTop );
    const double aBottom = tauYBottom * ( 1.0 / G - bBottom );
    const double dTop    = 0.5 * tauYTop * tauYTop * ( 1.0 / G - bTop );
    const double dBottom = 0.5 * tauYBottom * tauYBottom * ( 1.0 / G - bBottom );

    const double deltaB = bTop - bBottom;
    const double deltaA = aTop - aBottom;
    const double deltaD = dTop - dBottom;

    // Equation (S): -0.5*deltaB*t^2 - deltaA*t + deltaD = 0.
    const double quadraticA   = -0.5 * deltaB;
    const double quadraticB   = -deltaA;
    const double quadraticC   = deltaD;
    const double discriminant = quadraticB * quadraticB - 4.0 * quadraticA * quadraticC;
    throwExceptionOnFailure( discriminant >= 0.0,
                             "predictBilinearKinkAlpha: equation (S) has no real root for these material "
                             "parameters." );

    const double root1   = ( -quadraticB + std::sqrt( discriminant ) ) / ( 2.0 * quadraticA );
    const double root2   = ( -quadraticB - std::sqrt( discriminant ) ) / ( 2.0 * quadraticA );
    const double maxTauY = std::max( tauYTop, tauYBottom );

    double tStar = std::numeric_limits< double >::quiet_NaN();
    if ( root1 > maxTauY )
      tStar = root1;
    else if ( root2 > maxTauY )
      tStar = root2;
    else
      throwExceptionOnFailure( false,
                               "predictBilinearKinkAlpha: neither root of equation (S) exceeds "
                               "max(tauYTop,tauYBottom); no admissible both-plastic stationary traction." );

    const double gammaTopStar    = tauYTop / G + ( tStar - tauYTop ) / HtanTop;
    const double gammaBottomStar = tauYBottom / G + ( tStar - tauYBottom ) / HtanBottom;

    const double ubar = alphaTarget * gammaTopStar + ( 1.0 - alphaTarget ) * gammaBottomStar;
    const double dux  = ubar * h;

    return { alphaTarget, dux, tStar, gammaTopStar, gammaBottomStar };
  }

  // Builds the explicit [h,nBottom,bottom...,nTop,top...] VONMISES property
  // layout (7 base properties per side: E, nu, yieldStress, HLin,
  // deltaYieldStress, delta, density) for a given (sigmaY, HLin) pair per
  // side, sharing the same E, nu, h.
  std::vector< double > makeBilinearKinkVonMisesProperties( double E,
                                                            double nu,
                                                            double sigmaYBottom,
                                                            double HLinBottom,
                                                            double sigmaYTop,
                                                            double HLinTop,
                                                            double h )
  {
    return std::vector<
      double >{ h, 7., E, nu, sigmaYBottom, HLinBottom, 0., 0., 0., 7., E, nu, sigmaYTop, HLinTop, 0., 0., 0. };
  }

  // Runs one bilinear-kink case end to end: builds the properties, imposes
  // the closed-form-predicted tangential jump from a virgin state with alpha
  // evolution forced active (mirroring
  // testActiveAlphaSatisfiesLowerAndUpperBoundKKTConditions's use of the
  // forceAlphaEvolutionActive test hook), and checks the material's own
  // converged alpha against the independently-derived alphaTarget.
  void checkBilinearKinkCaseConvergesToClosedFormAlpha( const std::string& caseLabel,
                                                        double             E,
                                                        double             nu,
                                                        double             sigmaYTop,
                                                        double             HLinTop,
                                                        double             sigmaYBottom,
                                                        double             HLinBottom,
                                                        double             h,
                                                        double             alphaTarget,
                                                        double             alphaTolerance )
  {
    const BilinearKinkPrediction
      prediction = predictBilinearKinkAlpha( E, nu, sigmaYTop, HLinTop, sigmaYBottom, HLinBottom, h, alphaTarget );

    const std::vector< double >
      properties = makeBilinearKinkVonMisesProperties( E, nu, sigmaYBottom, HLinBottom, sigmaYTop, HLinTop, h );

    Eigen::Matrix< double, 21, 1 > generalizedIncrement = Eigen::Matrix< double, 21, 1 >::Zero();
    generalizedIncrement[0]                             = prediction.tangentialJump;

    Eigen::VectorXd stateVars;
    const auto      evaluation = evaluateExtendedMaterialWithState( "VONMISES",
                                                               properties.data(),
                                                               static_cast< int >( properties.size() ),
                                                               generalizedIncrement,
                                                               stateVars,
                                                               0.,
                                                               1.,
                                                               true );

    std::ostringstream message;
    message << "Bilinear-kink case '" << caseLabel << "': converged alpha=" << evaluation.alpha
            << " does not match the independently-derived closed-form alpha*=" << alphaTarget
            << " (criticalTraction=" << prediction.criticalTraction << ", gammaTop*=" << prediction.gammaTopAtCritical
            << ", gammaBottom*=" << prediction.gammaBottomAtCritical << ").";
    throwExceptionOnFailure( std::abs( evaluation.alpha - alphaTarget ) < alphaTolerance, message.str() );

    throwExceptionOnFailure( evaluation.response.allFinite() && evaluation.tangent.allFinite(),
                             "Bilinear-kink case '" + caseLabel + "': response or tangent contains nan/inf." );

    // Consistency check on the physical premise of the derivation: both
    // sublayers must actually have yielded (kappa > 0) for equation (S) to
    // apply. Query the layout through a throwaway material instance built
    // from the same properties (mirrors how
    // testSingleMaterialInputKeepsIndependentTopAndBottomState reads back
    // per-side state).
    MarmotExtendedInterfaceMaterialHypoElastic layoutProbe( "VONMISES",
                                                            properties.data(),
                                                            static_cast< int >( properties.size() ),
                                                            1 );
    const double topKappa    = layoutProbe.getStateView( "topMaterialStateVars", stateVars.data() ).stateLocation[0];
    const double bottomKappa = layoutProbe.getStateView( "bottomMaterialStateVars", stateVars.data() ).stateLocation[0];
    throwExceptionOnFailure( topKappa > 1e-10 && bottomKappa > 1e-10,
                             "Bilinear-kink case '" + caseLabel +
                               "': premise of the derivation (both sublayers plastically yielded) does not hold "
                               "(topKappa=" +
                               std::to_string( topKappa ) + ", bottomKappa=" + std::to_string( bottomKappa ) + ")." );
  }

  // testActiveAlphaSatisfiesLowerAndUpperBoundKKTConditions already shows
  // that, for two *purely elastic* (or equal-yield-stress plastic, see the
  // derivation above) mismatched sublayers, the alpha-stationarity condition
  // is monotonic and alpha is always driven to a KKT bound. This test uses
  // two Von Mises sublayers with *different* yield stresses AND different
  // (positive) linear hardening moduli, loaded in pure tangential slip past
  // both yield points, to reach a genuine, non-trivial *interior* stationary
  // point whose location is pinned down by the closed-form equations (S) and
  // (A) derived above -- a case the existing tests never exercise.
  void testAlphaConvergesToClosedFormBilinearKinkValue()
  {
    // Three independent (sigmaYTop/sigmaYBottom, HLinTop/HLinBottom) ratios,
    // each with its own target alpha well inside (alphaMinimum, 1-alphaMinimum),
    // to demonstrate the closed-form match is not a coincidence of one
    // particular parameter choice.
    checkBilinearKinkCaseConvergesToClosedFormAlpha( "A (yield ratio 2, hardening ratio 1/16)",
                                                     210000.,
                                                     0.3,
                                                     300.,
                                                     5000.,
                                                     150.,
                                                     80000.,
                                                     0.02,
                                                     0.3,
                                                     1e-6 );
    checkBilinearKinkCaseConvergesToClosedFormAlpha( "B (yield ratio 2, hardening ratio 1/30)",
                                                     210000.,
                                                     0.3,
                                                     400.,
                                                     2000.,
                                                     200.,
                                                     60000.,
                                                     0.02,
                                                     0.4,
                                                     1e-6 );
    checkBilinearKinkCaseConvergesToClosedFormAlpha( "C (yield ratio 5/3, hardening ratio 1/20)",
                                                     210000.,
                                                     0.3,
                                                     250.,
                                                     3000.,
                                                     150.,
                                                     60000.,
                                                     0.015,
                                                     0.6,
                                                     1e-6 );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = {
    testNaturalStressUpdateFailurePreservesCommittedState,
    testCutbackRetryFromPreservedSnapshotIsDeterministicAfterDiscardedAttempt,
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
    testLegacySingleMaterialLayoutReducesToStandardInterfaceForEqualSides,
    testLegacyVonMisesLayoutReducesToStandardInterfaceInElasticRegime,
    testLegacyLayoutWithIntegerValuedNuIsMisclassifiedAndRejected,
    testIndeterminateAlphaWithNonzeroNormalGradientJumpDoesNotFail,
    testAlphaRemainsFixedDuringElasticLoading,
    testCommittedPlasticActivityControlsAlphaEvolution,
    testActiveAlphaSatisfiesLowerAndUpperBoundKKTConditions,
    testAlphaConvergesToClosedFormBilinearKinkValue,
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
