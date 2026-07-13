#include "Marmot/MarmotExtendedInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Matrix3dRowMajor   = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
using Matrix9dRowMajor   = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
using Matrix3x9RowMajor  = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;
using Matrix9x3RowMajor  = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;
using Matrix9x21RowMajor = Eigen::Matrix< double, 9, 21, Eigen::RowMajor >;
using Matrix3x21RowMajor = Eigen::Matrix< double, 3, 21, Eigen::RowMajor >;
using Matrix4x21RowMajor = Eigen::Matrix< double, 4, 21, Eigen::RowMajor >;
using Matrix21dRowMajor  = Eigen::Matrix< double, 21, 21, Eigen::RowMajor >;
using Matrix4dRowMajor   = Eigen::Matrix< double, 4, 4, Eigen::RowMajor >;
using RowVector21d       = Eigen::Matrix< double, 1, 21 >;
using Vector9d           = Eigen::Matrix< double, 9, 1 >;

namespace {

  constexpr double alphaMinimum        = 1e-4;
  constexpr double alphaMaximum        = 1.0 - alphaMinimum;
  constexpr double alphaBoundTolerance = 1e-7;

  struct InterfaceKinematics {
    Eigen::Vector3d displacementJump;
    Eigen::Vector3d averageNormalGradient;
    Vector9d        averageSurfaceGradient;
    Vector9d        surfaceGradientJump;
  };

  struct FaceDisplacementGradients {
    Matrix3dRowMajor top;
    Matrix3dRowMajor bottom;
  };

  struct SideTrial {
    Marmot::Vector6d  stress               = Marmot::Vector6d::Zero();
    Marmot::Matrix6d  tangent              = Marmot::Matrix6d::Zero();
    Matrix3dRowMajor  stressTensor         = Matrix3dRowMajor::Zero();
    Matrix3dRowMajor  Q                    = Matrix3dRowMajor::Zero();
    Matrix3x9RowMajor H                    = Matrix3x9RowMajor::Zero();
    double            incrementalPotential = 0.0;
    bool              internalStateChanged = false;
  };

  struct MaterialTrial {
    SideTrial top;
    SideTrial bottom;

    Matrix3dRowMajor weightedStressTensor        = Matrix3dRowMajor::Zero();
    Matrix3dRowMajor weightedJumpResultantTensor = Matrix3dRowMajor::Zero();

    Eigen::Vector3d topTraction    = Eigen::Vector3d::Zero();
    Eigen::Vector3d bottomTraction = Eigen::Vector3d::Zero();
    Eigen::Vector3d tractionJump   = Eigen::Vector3d::Zero();

    double reducedPotential = 0.0;
    double alphaResidual    = 0.0;
  };

  struct LocalSolution {
    Eigen::Vector3d normalGradientJump   = Eigen::Vector3d::Zero();
    double          alpha                = 0.5;
    bool            alphaIsActive        = false;
    bool            alphaIsIndeterminate = false;
    MaterialTrial   trial;
  };

  Matrix4dRowMajor computeLocalHessian( const LocalSolution& solution );

  struct ExtendedResponse {
    Eigen::Vector3d force                = Eigen::Vector3d::Zero();
    Vector9d        averageSurfaceStress = Vector9d::Zero();
    Vector9d        jumpSurfaceStress    = Vector9d::Zero();
  };

  struct FixedLocalDerivatives {
    Matrix9x21RowMajor dTopStressDZ    = Matrix9x21RowMajor::Zero();
    Matrix9x21RowMajor dBottomStressDZ = Matrix9x21RowMajor::Zero();

    Matrix3x21RowMajor dTopTractionDZ    = Matrix3x21RowMajor::Zero();
    Matrix3x21RowMajor dBottomTractionDZ = Matrix3x21RowMajor::Zero();

    RowVector21d dTopPotentialDZ    = RowVector21d::Zero();
    RowVector21d dBottomPotentialDZ = RowVector21d::Zero();

    Matrix9x3RowMajor dTopStressDG    = Matrix9x3RowMajor::Zero();
    Matrix9x3RowMajor dBottomStressDG = Matrix9x3RowMajor::Zero();

    Vector9d dTopStressDAlpha    = Vector9d::Zero();
    Vector9d dBottomStressDAlpha = Vector9d::Zero();
  };

  bool isIntegerProperty( double value )
  {
    return std::abs( value - std::round( value ) ) < 1e-12;
  }

  Matrix3dRowMajor vectorToTensor( const Vector9d& vector )
  {
    return Eigen::Map< const Matrix3dRowMajor >( vector.data() );
  }

  Vector9d tensorToVector( const Matrix3dRowMajor& tensor )
  {
    return Eigen::Map< const Vector9d >( tensor.data() );
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

  template < class DeformationType >
  InterfaceKinematics extractInterfaceKinematics( const DeformationType& deformation, double h )
  {
    const Eigen::Map< const Eigen::Matrix< double, 6, 1 > >  dU( deformation.dU.data() );
    const Eigen::Map< const Eigen::Matrix< double, 18, 1 > > dSurfaceStrain( deformation.dSurfaceStrain.data() );

    InterfaceKinematics kinematics;
    kinematics.displacementJump       = dU.segment< 3 >( 0 ) - dU.segment< 3 >( 3 );
    kinematics.averageNormalGradient  = kinematics.displacementJump / h;
    kinematics.averageSurfaceGradient = 0.5 * ( dSurfaceStrain.segment< 9 >( 0 ) + dSurfaceStrain.segment< 9 >( 9 ) );
    kinematics.surfaceGradientJump    = dSurfaceStrain.segment< 9 >( 0 ) - dSurfaceStrain.segment< 9 >( 9 );
    return kinematics;
  }

  FaceDisplacementGradients reconstructFaceGradients( const Eigen::Vector3d& averageNormalGradient,
                                                      const Eigen::Vector3d& normalGradientJump,
                                                      const Vector9d&        averageSurfaceGradient,
                                                      const Vector9d&        surfaceGradientJump,
                                                      const Eigen::Vector3d& normal,
                                                      double                 alpha )
  {
    FaceDisplacementGradients gradients;
    gradients.top    = vectorToTensor( averageSurfaceGradient ) + 0.5 * vectorToTensor( surfaceGradientJump );
    gradients.bottom = vectorToTensor( averageSurfaceGradient ) - 0.5 * vectorToTensor( surfaceGradientJump );

    const Eigen::Vector3d topNormalGradient    = averageNormalGradient + ( 1.0 - alpha ) * normalGradientJump;
    const Eigen::Vector3d bottomNormalGradient = averageNormalGradient - alpha * normalGradientJump;

    gradients.top += topNormalGradient * normal.transpose();
    gradients.bottom += bottomNormalGradient * normal.transpose();
    return gradients;
  }

  Vector9d stressGradientColumn( const Marmot::Matrix6d& tangent, const Matrix3dRowMajor& dDisplacementGradient )
  {
    const Marmot::Vector6d dStress = tangent * strainToVoigt( dDisplacementGradient );
    return tensorToVector( stressToTensor( dStress ) );
  }

  Matrix9x3RowMajor computeNormalGradientStressMap( const Marmot::Matrix6d& tangent, const Eigen::Vector3d& normal )
  {
    Matrix9x3RowMajor map = Matrix9x3RowMajor::Zero();
    for ( int column = 0; column < 3; ++column ) {
      Matrix3dRowMajor dGradient = Matrix3dRowMajor::Zero();
      dGradient.row( column )    = normal.transpose();
      map.col( column )          = stressGradientColumn( tangent, dGradient );
    }
    return map;
  }

  Matrix3x9RowMajor stressToTractionMap( const Eigen::Vector3d& normal )
  {
    Matrix3x9RowMajor map = Matrix3x9RowMajor::Zero();
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        map( i, 3 * i + j ) = normal[j];
    return map;
  }

  double evaluateStressPowerAtPathParameter( MarmotMaterialHypoElastic&                 material,
                                             const Marmot::Vector6d&                    oldStress,
                                             const double*                              oldStateVars,
                                             int                                        nStateVars,
                                             const Marmot::Vector6d&                    strainIncrement,
                                             double                                     lambda,
                                             const MarmotMaterialHypoElastic::timeInfo& timeInfo )
  {
    std::vector< double > stateCopy( static_cast< size_t >( std::max( 0, nStateVars ) ) );
    if ( nStateVars > 0 )
      std::copy( oldStateVars, oldStateVars + nStateVars, stateCopy.begin() );

    Marmot::Vector6d                   stress  = oldStress;
    Marmot::Matrix6d                   tangent = Marmot::Matrix6d::Zero();
    MarmotMaterialHypoElastic::state3D trialState{ stress, 0.0, 0.0, nStateVars > 0 ? stateCopy.data() : nullptr };

    material.computeStress( trialState, tangent, lambda * strainIncrement, timeInfo );
    return trialState.stress.dot( strainIncrement );
  }

  double adaptiveSimpsonRecursive( const std::function< double( double ) >& f,
                                   double                                   a,
                                   double                                   b,
                                   double                                   fa,
                                   double                                   fm,
                                   double                                   fb,
                                   double                                   whole,
                                   double                                   tolerance,
                                   int                                      depth )
  {
    const double m   = 0.5 * ( a + b );
    const double lm  = 0.5 * ( a + m );
    const double rm  = 0.5 * ( m + b );
    const double flm = f( lm );
    const double frm = f( rm );

    const double left    = ( m - a ) * ( fa + 4.0 * flm + fm ) / 6.0;
    const double right   = ( b - m ) * ( fm + 4.0 * frm + fb ) / 6.0;
    const double refined = left + right;
    const double error   = refined - whole;

    if ( depth <= 0 || std::abs( error ) <= 15.0 * tolerance )
      return refined + error / 15.0;

    return adaptiveSimpsonRecursive( f, a, m, fa, flm, fm, left, 0.5 * tolerance, depth - 1 ) +
           adaptiveSimpsonRecursive( f, m, b, fm, frm, fb, right, 0.5 * tolerance, depth - 1 );
  }

  double computeCondensedIncrementalPotential( MarmotMaterialHypoElastic&                 material,
                                               const Marmot::Vector6d&                    oldStress,
                                               const double*                              oldStateVars,
                                               int                                        nStateVars,
                                               const Matrix3dRowMajor&                    displacementGradientIncrement,
                                               const MarmotMaterialHypoElastic::timeInfo& timeInfo )
  {
    const Marmot::Vector6d strainIncrement = strainToVoigt( displacementGradientIncrement );
    if ( strainIncrement.norm() <= 1e-18 )
      return 0.0;

    const auto integrand = [&]( double lambda ) {
      return evaluateStressPowerAtPathParameter( material,
                                                 oldStress,
                                                 oldStateVars,
                                                 nStateVars,
                                                 strainIncrement,
                                                 lambda,
                                                 timeInfo );
    };

    const double fa        = integrand( 0.0 );
    const double fm        = integrand( 0.5 );
    const double fb        = integrand( 1.0 );
    const double whole     = ( fa + 4.0 * fm + fb ) / 6.0;
    const double tolerance = 1e-10 * std::max( 1.0, std::abs( whole ) );
    return adaptiveSimpsonRecursive( integrand, 0.0, 1.0, fa, fm, fb, whole, tolerance, 12 );
  }

} // namespace

MarmotExtendedInterfaceMaterialHypoElastic::MarmotExtendedInterfaceMaterialHypoElastic(
  const std::string& materialName_,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : materialProperties( matProperties_ ),
    nMaterialProperties( nMaterialProperties_ ),
    materialName( materialName_ ),
    materialNumber( materialNumber_ )
{
  if ( nMaterialProperties < 3 )
    throw std::invalid_argument( "MarmotExtendedInterfaceMaterialHypoElastic requires material properties." );

  const bool hasExplicitTopBottomLayout = nMaterialProperties >= 5 && isIntegerProperty( materialProperties[1] );

  if ( hasExplicitTopBottomLayout ) {
    h                         = materialProperties[0];
    const int nBottom         = static_cast< int >( std::round( materialProperties[1] ) );
    const int topSizePosition = 2 + nBottom;

    if ( nBottom <= 0 || topSizePosition >= nMaterialProperties ||
         !isIntegerProperty( materialProperties[topSizePosition] ) )
      throw std::invalid_argument(
        "Invalid extended interface material layout. Expected [h,nBottom,bottom...,nTop,top...]." );

    const int nTop = static_cast< int >( std::round( materialProperties[topSizePosition] ) );
    if ( nTop <= 0 || topSizePosition + 1 + nTop != nMaterialProperties )
      throw std::invalid_argument(
        "Invalid extended interface material layout. Expected [h,nBottom,bottom...,nTop,top...]." );

    bottomMaterialProperties.assign( materialProperties + 2, materialProperties + 2 + nBottom );
    topMaterialProperties.assign( materialProperties + topSizePosition + 1,
                                  materialProperties + topSizePosition + 1 + nTop );
  }
  else {
    h = materialProperties[2];
    bottomMaterialProperties.reserve( nMaterialProperties - 1 );
    bottomMaterialProperties.push_back( materialProperties[0] );
    bottomMaterialProperties.push_back( materialProperties[1] );
    bottomMaterialProperties.insert( bottomMaterialProperties.end(),
                                     materialProperties + 3,
                                     materialProperties + nMaterialProperties );
    topMaterialProperties = bottomMaterialProperties;
  }

  if ( h <= 0.0 )
    throw std::invalid_argument( "MarmotExtendedInterfaceMaterialHypoElastic requires h > 0." );

  characteristicElementLength = 0.0;

  bottomMaterial.reset(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     bottomMaterialProperties.data(),
                                                                     static_cast< int >(
                                                                       bottomMaterialProperties.size() ),
                                                                     materialNumber ) );

  topMaterial.reset(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     topMaterialProperties.data(),
                                                                     static_cast< int >( topMaterialProperties.size() ),
                                                                     materialNumber ) );

  if ( !bottomMaterial || !topMaterial )
    throw std::invalid_argument( "Unknown base material for MarmotExtendedInterfaceMaterialHypoElastic: " +
                                 materialName );

  stateLayout.add( "bottomStress", 6 );
  stateLayout.add( "topStress", 6 );
  stateLayout.add( "normalGradientJump", 3 );
  stateLayout.add( "alpha", 1 );
  stateLayout.add( "alphaEvolutionActive", 1 );
  stateLayout.add( "bottomMaterialStateVars", bottomMaterial->getNumberOfRequiredStateVars() );
  stateLayout.add( "topMaterialStateVars", topMaterial->getNumberOfRequiredStateVars() );
  stateLayout.finalize();
}

void MarmotExtendedInterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  bottomMaterial->setCharacteristicElementLength( length );
  topMaterial->setCharacteristicElementLength( length );
}

namespace {

  SideTrial evaluateSideTrial( MarmotMaterialHypoElastic&                 material,
                               const Marmot::Vector6d&                    oldStress,
                               const double*                              oldStateVars,
                               int                                        nStateVars,
                               const Matrix3dRowMajor&                    displacementGradient,
                               const Eigen::Vector3d&                     normal,
                               const MarmotMaterialHypoElastic::timeInfo& timeInfo,
                               bool                                       computePotential,
                               bool                                       commit,
                               double*                                    committedStress,
                               double*                                    committedStateVars )
  {
    std::vector< double > stateCopy( static_cast< size_t >( std::max( 0, nStateVars ) ) );
    if ( nStateVars > 0 )
      std::copy( oldStateVars, oldStateVars + nStateVars, stateCopy.begin() );

    Marmot::Vector6d                   stress  = oldStress;
    Marmot::Matrix6d                   tangent = Marmot::Matrix6d::Zero();
    MarmotMaterialHypoElastic::state3D sideState{ stress, 0.0, 0.0, nStateVars > 0 ? stateCopy.data() : nullptr };

    const Marmot::Vector6d strainIncrement              = strainToVoigt( displacementGradient );
    double                 incrementalPotential         = 0.0;
    bool                   hasExactIncrementalPotential = false;
    if ( computePotential ) {
      hasExactIncrementalPotential = material.computeStressAndIncrementalPotential( sideState,
                                                                                    tangent,
                                                                                    strainIncrement,
                                                                                    timeInfo,
                                                                                    incrementalPotential );
    }
    else {
      material.computeStress( sideState, tangent, strainIncrement, timeInfo );
    }

    if ( computePotential && !hasExactIncrementalPotential )
      incrementalPotential = computeCondensedIncrementalPotential( material,
                                                                   oldStress,
                                                                   oldStateVars,
                                                                   nStateVars,
                                                                   displacementGradient,
                                                                   timeInfo );

    bool internalStateChanged = false;
    for ( int i = 0; i < nStateVars; ++i ) {
      const double scale = std::max( { 1.0, std::abs( oldStateVars[i] ), std::abs( stateCopy[i] ) } );
      if ( std::abs( stateCopy[i] - oldStateVars[i] ) > 1e-12 * scale ) {
        internalStateChanged = true;
        break;
      }
    }

    if ( commit ) {
      Eigen::Map< Marmot::Vector6d > committedStressMap( committedStress );
      committedStressMap = sideState.stress;
      if ( nStateVars > 0 )
        std::copy( stateCopy.begin(), stateCopy.end(), committedStateVars );
    }

    const auto normalTensor             = Marmot::FastorStandardTensors::Tensor3d( normal.data() );
    const auto [Z, QTensor, HTensor, Y] = Marmot::Materials::InterfaceMaterialHelperFunctions::
      calculateInterfaceMaterialParameters( normalTensor, tangent );

    SideTrial result;
    result.stress               = sideState.stress;
    result.tangent              = tangent;
    result.stressTensor         = stressToTensor( sideState.stress );
    result.Q                    = Eigen::Map< const Matrix3dRowMajor >( QTensor.data() );
    result.H                    = Eigen::Map< const Matrix3x9RowMajor >( HTensor.data() );
    result.incrementalPotential = incrementalPotential;
    result.internalStateChanged = internalStateChanged;

    (void)Z;
    (void)Y;
    return result;
  }

  MaterialTrial computeMaterialTrial( MarmotExtendedInterfaceMaterialHypoElastic& material,
                                      MarmotStateLayoutDynamic&                   stateLayout,
                                      double*                                     stateVars,
                                      const Eigen::Vector3d&                      averageNormalGradient,
                                      const Eigen::Vector3d&                      normalGradientJump,
                                      const Vector9d&                             averageSurfaceGradient,
                                      const Vector9d&                             surfaceGradientJump,
                                      const Eigen::Vector3d&                      normal,
                                      const MarmotMaterialHypoElastic::timeInfo&  timeInfo,
                                      double                                      alpha,
                                      bool                                        computePotential,
                                      bool                                        commit )
  {
    double* bottomStressPtr = stateLayout.getPtr( stateVars, "bottomStress" );
    double* topStressPtr    = stateLayout.getPtr( stateVars, "topStress" );
    double* bottomStatePtr  = stateLayout.getPtr( stateVars, "bottomMaterialStateVars" );
    double* topStatePtr     = stateLayout.getPtr( stateVars, "topMaterialStateVars" );

    const int nBottomStateVars = material.getStateView( "bottomMaterialStateVars", stateVars ).stateSize;
    const int nTopStateVars    = material.getStateView( "topMaterialStateVars", stateVars ).stateSize;

    const Marmot::Vector6d oldBottomStress = Eigen::Map< const Marmot::Vector6d >( bottomStressPtr );
    const Marmot::Vector6d oldTopStress    = Eigen::Map< const Marmot::Vector6d >( topStressPtr );

    const auto gradients = reconstructFaceGradients( averageNormalGradient,
                                                     normalGradientJump,
                                                     averageSurfaceGradient,
                                                     surfaceGradientJump,
                                                     normal,
                                                     alpha );

    MaterialTrial trial;
    trial.top = evaluateSideTrial( material.getTopMaterial(),
                                   oldTopStress,
                                   topStatePtr,
                                   nTopStateVars,
                                   gradients.top,
                                   normal,
                                   timeInfo,
                                   computePotential,
                                   commit,
                                   topStressPtr,
                                   topStatePtr );

    trial.bottom = evaluateSideTrial( material.getBottomMaterial(),
                                      oldBottomStress,
                                      bottomStatePtr,
                                      nBottomStateVars,
                                      gradients.bottom,
                                      normal,
                                      timeInfo,
                                      computePotential,
                                      commit,
                                      bottomStressPtr,
                                      bottomStatePtr );

    trial.topTraction    = trial.top.stressTensor * normal;
    trial.bottomTraction = trial.bottom.stressTensor * normal;
    trial.tractionJump   = trial.topTraction - trial.bottomTraction;

    trial.weightedStressTensor = alpha * trial.top.stressTensor + ( 1.0 - alpha ) * trial.bottom.stressTensor;

    trial.weightedJumpResultantTensor = 0.5 * ( alpha * trial.top.stressTensor -
                                                ( 1.0 - alpha ) * trial.bottom.stressTensor );

    if ( computePotential ) {
      trial.reducedPotential = alpha * trial.top.incrementalPotential +
                               ( 1.0 - alpha ) * trial.bottom.incrementalPotential;

      const Eigen::Vector3d weightedTraction = alpha * trial.topTraction + ( 1.0 - alpha ) * trial.bottomTraction;

      trial.alphaResidual = trial.top.incrementalPotential - trial.bottom.incrementalPotential -
                            weightedTraction.dot( normalGradientJump );
    }

    return trial;
  }

  double tractionResidualTolerance( const MaterialTrial& trial )
  {
    const double scale = std::max( { 1.0, trial.topTraction.norm(), trial.bottomTraction.norm() } );
    return 1e-9 + 1e-8 * scale;
  }

  Eigen::Vector3d solveNormalGradientJumpForAlpha( MarmotExtendedInterfaceMaterialHypoElastic& material,
                                                   MarmotStateLayoutDynamic&                   stateLayout,
                                                   double*                                     stateVars,
                                                   const Eigen::Vector3d&                      averageNormalGradient,
                                                   const Vector9d&                             averageSurfaceGradient,
                                                   const Vector9d&                             surfaceGradientJump,
                                                   const Eigen::Vector3d&                      normal,
                                                   const MarmotMaterialHypoElastic::timeInfo&  timeInfo,
                                                   double                                      alpha,
                                                   const Eigen::Vector3d&                      initialGuess )
  {
    constexpr int maxIterations           = 50;
    constexpr int maxLineSearchIterations = 25;

    Eigen::Vector3d g                = initialGuess;
    double          lastResidualNorm = std::numeric_limits< double >::infinity();
    double          lastTolerance    = 0.0;

    for ( int iteration = 0; iteration < maxIterations; ++iteration ) {
      const auto trial = computeMaterialTrial( material,
                                               stateLayout,
                                               stateVars,
                                               averageNormalGradient,
                                               g,
                                               averageSurfaceGradient,
                                               surfaceGradientJump,
                                               normal,
                                               timeInfo,
                                               alpha,
                                               false,
                                               false );

      const Eigen::Vector3d residual  = trial.tractionJump;
      const double          tolerance = tractionResidualTolerance( trial );
      lastResidualNorm                = residual.norm();
      lastTolerance                   = tolerance;
      if ( residual.norm() <= tolerance )
        return g;

      const Matrix3dRowMajor               jacobian = ( 1.0 - alpha ) * trial.top.Q + alpha * trial.bottom.Q;
      Eigen::FullPivLU< Matrix3dRowMajor > lu( jacobian );
      if ( !lu.isInvertible() )
        throw Marmot::StressUpdateFailed(
          "MarmotExtendedInterfaceMaterialHypoElastic: local acoustic Jacobian is singular." );

      const Eigen::Vector3d direction = lu.solve( residual );
      if ( !direction.allFinite() )
        throw Marmot::StressUpdateFailed(
          "MarmotExtendedInterfaceMaterialHypoElastic: non-finite local Newton direction." );

      const double phi0       = 0.5 * residual.squaredNorm();
      double       stepLength = 1.0;
      bool         accepted   = false;

      for ( int lineSearchIteration = 0; lineSearchIteration < maxLineSearchIterations; ++lineSearchIteration ) {
        const Eigen::Vector3d candidate      = g - stepLength * direction;
        const auto            candidateTrial = computeMaterialTrial( material,
                                                          stateLayout,
                                                          stateVars,
                                                          averageNormalGradient,
                                                          candidate,
                                                          averageSurfaceGradient,
                                                          surfaceGradientJump,
                                                          normal,
                                                          timeInfo,
                                                          alpha,
                                                          false,
                                                          false );

        const double candidatePhi = 0.5 * candidateTrial.tractionJump.squaredNorm();
        if ( candidateTrial.tractionJump.norm() <= tractionResidualTolerance( candidateTrial ) ||
             candidatePhi <= phi0 - 1e-4 * stepLength * residual.squaredNorm() ) {
          g        = candidate;
          accepted = true;
          break;
        }
        stepLength *= 0.5;
      }

      if ( !accepted )
        throw Marmot::StressUpdateFailed(
          "MarmotExtendedInterfaceMaterialHypoElastic: traction-equilibrium line search failed." );
    }

    std::ostringstream message;
    message << "MarmotExtendedInterfaceMaterialHypoElastic: traction-equilibrium Newton iteration failed"
            << " (residual=" << lastResidualNorm << ", tolerance=" << lastTolerance << ", alpha=" << alpha
            << ", g=" << g.transpose() << ").";
    throw Marmot::StressUpdateFailed( message.str() );
  }

  double alphaResidualTolerance( const LocalSolution& solution )
  {
    const double scale = std::max( { 1.0,
                                     std::abs( solution.trial.top.incrementalPotential ),
                                     std::abs( solution.trial.bottom.incrementalPotential ),
                                     std::abs( solution.trial.alphaResidual ) } );
    return 1e-8 * scale + 5e-8;
  }

  Eigen::Vector4d computeLocalResidual( const LocalSolution& solution )
  {
    Eigen::Vector4d residual;
    const double    A    = solution.alpha * ( 1.0 - solution.alpha );
    residual.head< 3 >() = A * solution.trial.tractionJump;
    residual[3]          = solution.trial.alphaResidual;
    return residual;
  }

  double projectedAlphaResidual( const LocalSolution& solution )
  {
    if ( solution.alpha <= alphaMinimum + alphaBoundTolerance )
      return std::min( 0.0, solution.trial.alphaResidual );
    if ( solution.alpha >= alphaMaximum - alphaBoundTolerance )
      return std::max( 0.0, solution.trial.alphaResidual );
    return solution.trial.alphaResidual;
  }

  LocalSolution solveCoupledLocalProblem( MarmotExtendedInterfaceMaterialHypoElastic& material,
                                          MarmotStateLayoutDynamic&                   stateLayout,
                                          double*                                     stateVars,
                                          const Eigen::Vector3d&                      averageNormalGradient,
                                          const Vector9d&                             averageSurfaceGradient,
                                          const Vector9d&                             surfaceGradientJump,
                                          const Eigen::Vector3d&                      normal,
                                          const MarmotMaterialHypoElastic::timeInfo&  timeInfo,
                                          double                                      initialAlpha,
                                          const Eigen::Vector3d&                      initialG,
                                          bool                                        alphaEvolutionActive )
  {
    constexpr int    maxIterations           = 100;
    constexpr int    maxLineSearchIterations = 25;
    constexpr double armijoParameter         = 1e-4;

    const auto evaluate = [&]( const Eigen::Vector3d& g, double alpha ) {
      LocalSolution solution;
      solution.normalGradientJump = g;
      solution.alpha              = std::clamp( alpha, alphaMinimum, alphaMaximum );
      solution.trial              = computeMaterialTrial( material,
                                             stateLayout,
                                             stateVars,
                                             averageNormalGradient,
                                             solution.normalGradientJump,
                                             averageSurfaceGradient,
                                             surfaceGradientJump,
                                             normal,
                                             timeInfo,
                                             solution.alpha,
                                             true,
                                             false );
      return solution;
    };

    LocalSolution current = evaluate( initialG, std::clamp( initialAlpha, alphaMinimum, alphaMaximum ) );

    // Opt-in diagnostic used by the regression test and available when a new
    // constitutive model is connected to this potential-based formulation.
    // It is intentionally evaluated away from the converged state so all
    // off-equilibrium terms in the coupled Hessian are exercised.
    if ( alphaEvolutionActive && std::getenv( "MARMOT_EI_VALIDATE_LOCAL_DERIVATIVES" ) != nullptr ) {
      const Eigen::Vector4d x = ( Eigen::Vector4d() << current.normalGradientJump, current.alpha ).finished();
      Eigen::Vector4d       finiteDifferenceGradient = Eigen::Vector4d::Zero();
      Eigen::Matrix4d       finiteDifferenceHessian  = Eigen::Matrix4d::Zero();
      for ( int column = 0; column < 4; ++column ) {
        const double    perturbation = 1e-7 * std::max( 1.0, std::abs( x[column] ) );
        Eigen::Vector4d plusX        = x;
        Eigen::Vector4d minusX       = x;
        plusX[column] += perturbation;
        minusX[column] -= perturbation;
        plusX[3]                         = std::clamp( plusX[3], alphaMinimum, alphaMaximum );
        minusX[3]                        = std::clamp( minusX[3], alphaMinimum, alphaMaximum );
        const LocalSolution plus         = evaluate( plusX.head< 3 >(), plusX[3] );
        const LocalSolution minus        = evaluate( minusX.head< 3 >(), minusX[3] );
        const double        denominator  = plusX[column] - minusX[column];
        finiteDifferenceGradient[column] = ( plus.trial.reducedPotential - minus.trial.reducedPotential ) / denominator;
        finiteDifferenceHessian.col( column ) = ( computeLocalResidual( plus ) - computeLocalResidual( minus ) ) /
                                                denominator;
      }

      const Eigen::Vector4d residual = computeLocalResidual( current );
      const Eigen::Matrix4d hessian  = computeLocalHessian( current );
      const double gradientError = ( finiteDifferenceGradient - residual ).norm() / std::max( 1.0, residual.norm() );
      const double hessianError  = ( finiteDifferenceHessian - hessian ).norm() / std::max( 1.0, hessian.norm() );
      const double symmetryError = ( hessian - hessian.transpose() ).norm() / std::max( 1.0, hessian.norm() );
      if ( gradientError > 2e-6 || hessianError > 2e-6 || symmetryError > 2e-10 ) {
        std::ostringstream message;
        message << "MarmotExtendedInterfaceMaterialHypoElastic: local derivative validation failed"
                << " (gradientError=" << gradientError << ", hessianError=" << hessianError
                << ", symmetryError=" << symmetryError << ").";
        throw Marmot::StressUpdateFailed( message.str() );
      }
    }

    // Keep the constitutive branch fixed throughout the increment.  The flag
    // comes from the previously committed increment and is updated only in
    // the trial state written below, so Newton retries and cutbacks see the
    // same alpha equation.
    if ( !alphaEvolutionActive ) {
      current.normalGradientJump   = solveNormalGradientJumpForAlpha( material,
                                                                    stateLayout,
                                                                    stateVars,
                                                                    averageNormalGradient,
                                                                    averageSurfaceGradient,
                                                                    surfaceGradientJump,
                                                                    normal,
                                                                    timeInfo,
                                                                    current.alpha,
                                                                    current.normalGradientJump );
      current                      = evaluate( current.normalGradientJump, current.alpha );
      current.alphaIsIndeterminate = true;
      return current;
    }

    if ( averageNormalGradient.norm() <= 1e-18 && averageSurfaceGradient.norm() <= 1e-18 &&
         surfaceGradientJump.norm() <= 1e-18 ) {
      current                      = evaluate( Eigen::Vector3d::Zero(), current.alpha );
      current.alphaIsIndeterminate = true;
      return current;
    }

    auto normalizedMerit = []( const LocalSolution& solution ) {
      const double tractionTolerance = tractionResidualTolerance( solution.trial );
      const double kktTolerance      = alphaResidualTolerance( solution );
      const double tractionMeasure   = solution.trial.tractionJump.norm() / tractionTolerance;
      const double alphaMeasure      = projectedAlphaResidual( solution ) / kktTolerance;
      return 0.5 * ( tractionMeasure * tractionMeasure + alphaMeasure * alphaMeasure );
    };

    const auto plasticActiveSetChanged = []( const LocalSolution& first, const LocalSolution& second ) {
      return first.trial.top.internalStateChanged != second.trial.top.internalStateChanged ||
             first.trial.bottom.internalStateChanged != second.trial.bottom.internalStateChanged;
    };

    for ( int iteration = 0; iteration < maxIterations; ++iteration ) {
      if ( current.alpha <= alphaMinimum + alphaBoundTolerance && current.alpha != alphaMinimum )
        current = evaluate( current.normalGradientJump, alphaMinimum );
      else if ( current.alpha >= alphaMaximum - alphaBoundTolerance && current.alpha != alphaMaximum )
        current = evaluate( current.normalGradientJump, alphaMaximum );

      const double tractionTolerance = tractionResidualTolerance( current.trial );
      const double kktTolerance      = alphaResidualTolerance( current );
      const bool   tractionConverged = current.trial.tractionJump.norm() <= tractionTolerance;
      const bool   alphaConverged    = std::abs( projectedAlphaResidual( current ) ) <= kktTolerance;

      // When the normal-gradient jump vanishes, alpha is a local gauge: both
      // reconstructed face gradients are insensitive to moving the kink to
      // first relevant order.  At very small cutbacks the path-integrated
      // potentials retain more quadrature/return-map noise than this effect.
      if ( tractionConverged && current.normalGradientJump.norm() <= 1e-7 ) {
        current.alphaIsIndeterminate = true;
        return current;
      }

      if ( tractionConverged && alphaConverged ) {
        const bool lowerActive = current.alpha <= alphaMinimum + alphaBoundTolerance &&
                                 current.trial.alphaResidual >= -kktTolerance;
        const bool upperActive = current.alpha >= alphaMaximum - alphaBoundTolerance &&
                                 current.trial.alphaResidual <= kktTolerance;
        current.alphaIsActive = lowerActive || upperActive;
        if ( !current.alphaIsActive && current.normalGradientJump.norm() <= 1e-14 )
          current.alphaIsIndeterminate = true;
        return current;
      }

      const bool lowerActive = current.alpha <= alphaMinimum + alphaBoundTolerance &&
                               current.trial.alphaResidual >= -kktTolerance;
      const bool upperActive = current.alpha >= alphaMaximum - alphaBoundTolerance &&
                               current.trial.alphaResidual <= kktTolerance;
      const bool alphaActive = lowerActive || upperActive;

      const Eigen::Vector4d  residual = computeLocalResidual( current );
      const Matrix4dRowMajor hessian  = computeLocalHessian( current );
      Eigen::Vector4d        step     = Eigen::Vector4d::Zero();

      if ( alphaActive ) {
        const Matrix3dRowMajor K = ( 1.0 - current.alpha ) * current.trial.top.Q +
                                   current.alpha * current.trial.bottom.Q;
        Eigen::FullPivLU< Matrix3dRowMajor > lu( K );
        if ( !lu.isInvertible() )
          throw Marmot::StressUpdateFailed(
            "MarmotExtendedInterfaceMaterialHypoElastic: active-bound acoustic Jacobian is singular." );
        step.head< 3 >() = -lu.solve( current.trial.tractionJump );
        // Plastic return maps can leave a stress-level residual a few ulps
        // above the strict traction tolerance even though the corresponding
        // kinematic Newton correction is at machine precision.
        if ( step.head< 3 >().norm() <= 1e-12 * std::max( 1.0, current.normalGradientJump.norm() ) ) {
          current.alphaIsActive = true;
          return current;
        }
      }
      else {
        Eigen::FullPivLU< Matrix4dRowMajor > lu( hessian );
        if ( lu.isInvertible() )
          step = -lu.solve( residual );

        // Newton is used only when it is a potential-descent direction.
        // Otherwise modify the symmetric Hessian eigenvalues to obtain one.
        if ( !step.allFinite() || residual.dot( step ) >= 0.0 ) {
          const Eigen::Matrix4d                            symmetricHessian = 0.5 * ( Eigen::Matrix4d( hessian ) +
                                                           Eigen::Matrix4d( hessian ).transpose() );
          Eigen::SelfAdjointEigenSolver< Eigen::Matrix4d > eigenSolver( symmetricHessian );
          if ( eigenSolver.info() != Eigen::Success )
            step.setZero();
          else {
            const double          eigenvalueScale    = std::max( 1.0, eigenSolver.eigenvalues().cwiseAbs().maxCoeff() );
            const double          eigenvalueFloor    = 1e-8 * eigenvalueScale;
            const Eigen::Vector4d inverseEigenvalues = eigenSolver.eigenvalues().unaryExpr(
              [=]( double value ) { return 1.0 / std::max( value, eigenvalueFloor ); } );
            step = -eigenSolver.eigenvectors() * inverseEigenvalues.asDiagonal() *
                   eigenSolver.eigenvectors().transpose() * residual;
          }
        }
      }

      if ( !step.allFinite() )
        step.setZero();

      const Eigen::Vector4d x          = ( Eigen::Vector4d() << current.normalGradientJump, current.alpha ).finished();
      bool                  accepted   = false;
      double                stepLength = 1.0;
      const double          merit0     = normalizedMerit( current );

      // First globalize with the incremental potential. Projection is applied
      // before testing Armijo, so the directional derivative uses the actual
      // feasible displacement.
      for ( int lineSearchIteration = 0; lineSearchIteration < maxLineSearchIterations; ++lineSearchIteration ) {
        Eigen::Vector4d candidateX       = x + stepLength * step;
        candidateX[3]                    = std::clamp( candidateX[3], alphaMinimum, alphaMaximum );
        const Eigen::Vector4d actualStep = candidateX - x;
        try {
          LocalSolution candidate             = evaluate( candidateX.head< 3 >(), candidateX[3] );
          const double  directionalDerivative = residual.dot( actualStep );
          const double  tractionPhi           = 0.5 * current.trial.tractionJump.squaredNorm();
          const double  candidateTractionPhi  = 0.5 * candidate.trial.tractionJump.squaredNorm();
          const bool    activeBoundAccepted   = alphaActive &&
                                           ( candidate.trial.tractionJump.norm() <=
                                               tractionResidualTolerance( candidate.trial ) ||
                                             candidateTractionPhi <=
                                               tractionPhi - armijoParameter * stepLength *
                                                               current.trial.tractionJump.squaredNorm() );
          const bool interiorAccepted = !alphaActive && directionalDerivative < 0.0 &&
                                        candidate.trial.reducedPotential <=
                                          current.trial.reducedPotential + armijoParameter * directionalDerivative;
          const bool branchChangeAccepted = !plasticActiveSetChanged( current, candidate ) ||
                                            normalizedMerit( candidate ) < merit0;
          if ( ( activeBoundAccepted || interiorAccepted ) && branchChangeAccepted ) {
            current  = std::move( candidate );
            accepted = true;
            break;
          }
        }
        catch ( const Marmot::StressUpdateFailed& ) {
        }
        stepLength *= 0.5;
      }

      if ( accepted )
        continue;

      // Fallback: scaled Levenberg-Marquardt steps minimize the traction and
      // projected KKT residuals without allowing the A=a(1-a) factor to mask
      // traction disequilibrium near a bound.
      const double    A                = current.alpha * ( 1.0 - current.alpha );
      Eigen::Matrix4d scaledJacobian   = hessian;
      Eigen::Vector4d scaledResidual   = residual;
      const double    tractionRowScale = 1.0 / ( std::max( A, 1e-12 ) * tractionTolerance );
      scaledJacobian.topRows< 3 >() *= tractionRowScale;
      scaledResidual.head< 3 >() *= tractionRowScale;
      scaledJacobian.row( 3 ) /= kktTolerance;
      scaledResidual[3] = projectedAlphaResidual( current ) / kktTolerance;

      if ( alphaActive ) {
        scaledJacobian.row( 3 ).setZero();
        scaledJacobian.col( 3 ).setZero();
        scaledJacobian( 3, 3 ) = 1.0;
        scaledResidual[3]      = 0.0;
      }

      const Eigen::Matrix4d normalMatrix        = scaledJacobian.transpose() * scaledJacobian;
      const Eigen::Vector4d normalResidual      = scaledJacobian.transpose() * scaledResidual;
      const double          regularizationScale = std::max( 1.0, normalMatrix.diagonal().cwiseAbs().maxCoeff() );

      for ( int lmIteration = 0; lmIteration < 12 && !accepted; ++lmIteration ) {
        const double    lambda      = std::pow( 10.0, lmIteration - 8 ) * regularizationScale;
        Eigen::Matrix4d regularized = normalMatrix;
        regularized.diagonal().array() += lambda;
        Eigen::LDLT< Eigen::Matrix4d > ldlt( regularized );
        if ( ldlt.info() != Eigen::Success )
          continue;
        const Eigen::Vector4d lmStep = -ldlt.solve( normalResidual );
        if ( !lmStep.allFinite() )
          continue;

        double lmStepLength = 1.0;
        for ( int lineSearchIteration = 0; lineSearchIteration < maxLineSearchIterations; ++lineSearchIteration ) {
          Eigen::Vector4d candidateX = x + lmStepLength * lmStep;
          candidateX[3]              = std::clamp( candidateX[3], alphaMinimum, alphaMaximum );
          try {
            LocalSolution candidate = evaluate( candidateX.head< 3 >(), candidateX[3] );
            if ( normalizedMerit( candidate ) < merit0 ) {
              current  = std::move( candidate );
              accepted = true;
              break;
            }
          }
          catch ( const Marmot::StressUpdateFailed& ) {
          }
          lmStepLength *= 0.5;
        }
      }

      if ( !accepted ) {
        std::ostringstream message;
        message << "MarmotExtendedInterfaceMaterialHypoElastic: coupled local globalization failed"
                << " (tractionResidual=" << current.trial.tractionJump.norm()
                << ", tractionTolerance=" << tractionTolerance << ", alphaResidual=" << current.trial.alphaResidual
                << ", alphaTolerance=" << kktTolerance << ", alpha=" << current.alpha
                << ", g=" << current.normalGradientJump.transpose() << ").";
        throw Marmot::StressUpdateFailed( message.str() );
      }
    }

    const double finalTractionTolerance = tractionResidualTolerance( current.trial );
    const double finalKktTolerance      = alphaResidualTolerance( current );
    if ( current.trial.tractionJump.norm() <= finalTractionTolerance &&
         std::abs( projectedAlphaResidual( current ) ) <= finalKktTolerance ) {
      const bool lowerActive = current.alpha <= alphaMinimum + alphaBoundTolerance &&
                               current.trial.alphaResidual >= -finalKktTolerance;
      const bool upperActive = current.alpha >= alphaMaximum - alphaBoundTolerance &&
                               current.trial.alphaResidual <= finalKktTolerance;
      current.alphaIsActive = lowerActive || upperActive;
      if ( !current.alphaIsActive && current.normalGradientJump.norm() <= 1e-14 )
        current.alphaIsIndeterminate = true;
      return current;
    }

    std::ostringstream message;
    message << "MarmotExtendedInterfaceMaterialHypoElastic: coupled local iteration failed"
            << " (tractionResidual=" << current.trial.tractionJump.norm()
            << ", alphaResidual=" << current.trial.alphaResidual << ", alpha=" << current.alpha
            << ", g=" << current.normalGradientJump.transpose() << ").";
    throw Marmot::StressUpdateFailed( message.str() );
  }

  Matrix9x21RowMajor computeFixedStressDerivativeDZ( const Marmot::Matrix6d& tangent,
                                                     const Eigen::Vector3d&  normal,
                                                     double                  h,
                                                     double                  jumpSurfaceSign )
  {
    Matrix9x21RowMajor derivative = Matrix9x21RowMajor::Zero();

    for ( int column = 0; column < 3; ++column ) {
      Matrix3dRowMajor dGradient = Matrix3dRowMajor::Zero();
      dGradient.row( column )    = normal.transpose() / h;
      derivative.col( column )   = stressGradientColumn( tangent, dGradient );
    }

    for ( int column = 0; column < 9; ++column ) {
      Matrix3dRowMajor dGradient          = Matrix3dRowMajor::Zero();
      dGradient( column / 3, column % 3 ) = 1.0;
      derivative.col( 3 + column )        = stressGradientColumn( tangent, dGradient );

      dGradient.setZero();
      dGradient( column / 3, column % 3 ) = jumpSurfaceSign;
      derivative.col( 12 + column )       = stressGradientColumn( tangent, dGradient );
    }
    return derivative;
  }

  RowVector21d computeFixedPotentialDerivativeDZ( const Matrix3dRowMajor& stressTensor,
                                                  const Eigen::Vector3d&  traction,
                                                  double                  h,
                                                  double                  jumpSurfaceSign )
  {
    RowVector21d derivative       = RowVector21d::Zero();
    derivative.segment< 3 >( 0 )  = traction.transpose() / h;
    const Vector9d stressVector   = tensorToVector( stressTensor );
    derivative.segment< 9 >( 3 )  = stressVector.transpose();
    derivative.segment< 9 >( 12 ) = jumpSurfaceSign * stressVector.transpose();
    return derivative;
  }

  FixedLocalDerivatives computeFixedLocalDerivatives( const MaterialTrial&   trial,
                                                      const Eigen::Vector3d& normal,
                                                      const Eigen::Vector3d& g,
                                                      double                 h,
                                                      double                 alpha )
  {
    FixedLocalDerivatives derivatives;
    derivatives.dTopStressDZ    = computeFixedStressDerivativeDZ( trial.top.tangent, normal, h, 0.5 );
    derivatives.dBottomStressDZ = computeFixedStressDerivativeDZ( trial.bottom.tangent, normal, h, -0.5 );

    const Matrix3x9RowMajor tractionMap = stressToTractionMap( normal );
    derivatives.dTopTractionDZ          = tractionMap * derivatives.dTopStressDZ;
    derivatives.dBottomTractionDZ       = tractionMap * derivatives.dBottomStressDZ;

    derivatives.dTopPotentialDZ    = computeFixedPotentialDerivativeDZ( trial.top.stressTensor,
                                                                     trial.topTraction,
                                                                     h,
                                                                     0.5 );
    derivatives.dBottomPotentialDZ = computeFixedPotentialDerivativeDZ( trial.bottom.stressTensor,
                                                                        trial.bottomTraction,
                                                                        h,
                                                                        -0.5 );

    const Matrix9x3RowMajor topNormalMap    = computeNormalGradientStressMap( trial.top.tangent, normal );
    const Matrix9x3RowMajor bottomNormalMap = computeNormalGradientStressMap( trial.bottom.tangent, normal );

    derivatives.dTopStressDG        = ( 1.0 - alpha ) * topNormalMap;
    derivatives.dBottomStressDG     = -alpha * bottomNormalMap;
    derivatives.dTopStressDAlpha    = -topNormalMap * g;
    derivatives.dBottomStressDAlpha = -bottomNormalMap * g;
    return derivatives;
  }

  Matrix4dRowMajor computeLocalHessian( const LocalSolution& solution )
  {
    const double           a     = solution.alpha;
    const double           A     = a * ( 1.0 - a );
    const auto&            trial = solution.trial;
    const Eigen::Vector3d& g     = solution.normalGradientJump;
    const Eigen::Vector3d  r     = trial.tractionJump;

    const Matrix3dRowMajor K      = ( 1.0 - a ) * trial.top.Q + a * trial.bottom.Q;
    const Matrix3dRowMajor deltaQ = trial.top.Q - trial.bottom.Q;
    const Matrix3dRowMajor qHat   = a * trial.top.Q + ( 1.0 - a ) * trial.bottom.Q;

    Matrix4dRowMajor hessian      = Matrix4dRowMajor::Zero();
    hessian.block< 3, 3 >( 0, 0 ) = A * K;
    hessian.block< 3, 1 >( 0, 3 ) = ( 1.0 - 2.0 * a ) * r - A * deltaQ * g;
    hessian.block< 1, 3 >( 3, 0 ) = ( 1.0 - 2.0 * a ) * r.transpose() - A * g.transpose() * deltaQ;
    hessian( 3, 3 )               = -2.0 * r.dot( g ) + g.dot( qHat * g );
    return hessian;
  }

  Matrix4x21RowMajor computeLocalExternalDerivative( const LocalSolution&         solution,
                                                     const FixedLocalDerivatives& derivatives )
  {
    const double           a = solution.alpha;
    const double           A = a * ( 1.0 - a );
    const Eigen::Vector3d& g = solution.normalGradientJump;

    Matrix4x21RowMajor derivative     = Matrix4x21RowMajor::Zero();
    derivative.block< 3, 21 >( 0, 0 ) = A * ( derivatives.dTopTractionDZ - derivatives.dBottomTractionDZ );

    derivative.row( 3 ) = derivatives.dTopPotentialDZ - derivatives.dBottomPotentialDZ -
                          g.transpose() *
                            ( a * derivatives.dTopTractionDZ + ( 1.0 - a ) * derivatives.dBottomTractionDZ );
    return derivative;
  }

  Matrix21dRowMajor computeCondensedTangent( LocalSolution& solution, const Eigen::Vector3d& normal, double h )
  {
    const double           a     = solution.alpha;
    const auto&            trial = solution.trial;
    const Eigen::Vector3d& g     = solution.normalGradientJump;

    const FixedLocalDerivatives fixed = computeFixedLocalDerivatives( trial, normal, g, h, a );

    Matrix4x21RowMajor dLocalDZ = Matrix4x21RowMajor::Zero();

    const auto condenseWithFixedAlpha = [&]() {
      const Matrix3dRowMajor               K = ( 1.0 - a ) * trial.top.Q + a * trial.bottom.Q;
      Eigen::FullPivLU< Matrix3dRowMajor > lu( K );
      if ( !lu.isInvertible() )
        throw Marmot::StressUpdateFailed(
          "MarmotExtendedInterfaceMaterialHypoElastic: condensed g-Jacobian is singular." );

      dLocalDZ.block< 3, 21 >( 0, 0 ) = -lu.solve( fixed.dTopTractionDZ - fixed.dBottomTractionDZ );
    };

    if ( solution.alphaIsActive || solution.alphaIsIndeterminate ) {
      condenseWithFixedAlpha();
    }
    else {
      const Matrix4dRowMajor               hessian = computeLocalHessian( solution );
      Eigen::FullPivLU< Matrix4dRowMajor > lu( hessian );
      if ( lu.isInvertible() ) {
        dLocalDZ = -lu.solve( computeLocalExternalDerivative( solution, fixed ) );
      }
      else {
        // Alpha is a free local gauge when the reduced potential has no
        // curvature in the kink-position direction.  This occurs, for
        // example, for identical sublayers at vanishingly small strain.
        // Keep the converged alpha fixed and condense only the regular
        // traction-equilibrium block instead of rejecting every cutback.
        solution.alphaIsIndeterminate = true;
        condenseWithFixedAlpha();
      }
    }

    const Eigen::Matrix< double, 3, 21, Eigen::RowMajor > dgDZ = dLocalDZ.block< 3, 21 >( 0, 0 );
    const RowVector21d                                    daDZ = dLocalDZ.row( 3 );

    Matrix9x21RowMajor dTopStress    = fixed.dTopStressDZ + fixed.dTopStressDG * dgDZ + fixed.dTopStressDAlpha * daDZ;
    Matrix9x21RowMajor dBottomStress = fixed.dBottomStressDZ + fixed.dBottomStressDG * dgDZ +
                                       fixed.dBottomStressDAlpha * daDZ;

    const Vector9d topStressVector    = tensorToVector( trial.top.stressTensor );
    const Vector9d bottomStressVector = tensorToVector( trial.bottom.stressTensor );

    Matrix9x21RowMajor dWeightedStress = a * dTopStress + ( 1.0 - a ) * dBottomStress +
                                         ( topStressVector - bottomStressVector ) * daDZ;

    Matrix9x21RowMajor dWeightedJumpResultant = 0.5 * ( a * dTopStress - ( 1.0 - a ) * dBottomStress +
                                                        ( topStressVector + bottomStressVector ) * daDZ );

    const Matrix3x9RowMajor tractionMap = stressToTractionMap( normal );
    Matrix21dRowMajor       tangent     = Matrix21dRowMajor::Zero();
    tangent.block< 3, 21 >( 0, 0 )      = tractionMap * dWeightedStress;
    tangent.block< 9, 21 >( 3, 0 )      = h * dWeightedStress;
    tangent.block< 9, 21 >( 12, 0 )     = h * dWeightedJumpResultant;
    return tangent;
  }

  ExtendedResponse makeExtendedResponse( const MaterialTrial& trial, const Eigen::Vector3d& normal, double h )
  {
    ExtendedResponse response;
    response.force                = trial.weightedStressTensor * normal;
    response.averageSurfaceStress = h * tensorToVector( trial.weightedStressTensor );
    response.jumpSurfaceStress    = h * tensorToVector( trial.weightedJumpResultantTensor );
    return response;
  }

  template < class StateType >
  void writeResponseToState( StateType& state, const ExtendedResponse& response )
  {
    state.force                = Marmot::FastorStandardTensors::Tensor3d( response.force.data() );
    state.averageSurfaceStress = Marmot::FastorStandardTensors::Tensor33d(
      vectorToTensor( response.averageSurfaceStress ).data() );
    state.jumpSurfaceStress = Marmot::FastorStandardTensors::Tensor33d(
      vectorToTensor( response.jumpSurfaceStress ).data() );
  }

  template < class TangentsType >
  void writeTangentBlocks( TangentsType& tangents, const Matrix21dRowMajor& tangent )
  {
    Eigen::Map< Matrix3dRowMajor >( tangents.forceJumpU )                   = tangent.block< 3, 3 >( 0, 0 );
    Eigen::Map< Matrix3x9RowMajor >( tangents.forceAverageSurfaceGradient ) = tangent.block< 3, 9 >( 0, 3 );
    Eigen::Map< Matrix3x9RowMajor >( tangents.forceJumpSurfaceGradient )    = tangent.block< 3, 9 >( 0, 12 );

    Eigen::Map< Matrix9x3RowMajor >( tangents.averageSurfaceStressJumpU ) = tangent.block< 9, 3 >( 3, 0 );
    Eigen::Map< Matrix9dRowMajor >( tangents.averageSurfaceStressAverageSurfaceGradient ) = tangent.block< 9, 9 >( 3,
                                                                                                                   3 );
    Eigen::Map< Matrix9dRowMajor >( tangents.averageSurfaceStressJumpSurfaceGradient ) = tangent.block< 9, 9 >( 3, 12 );

    Eigen::Map< Matrix9x3RowMajor >( tangents.jumpSurfaceStressJumpU )                 = tangent.block< 9, 3 >( 12, 0 );
    Eigen::Map< Matrix9dRowMajor >( tangents.jumpSurfaceStressAverageSurfaceGradient ) = tangent.block< 9, 9 >( 12, 3 );
    Eigen::Map< Matrix9dRowMajor >( tangents.jumpSurfaceStressJumpSurfaceGradient ) = tangent.block< 9, 9 >( 12, 12 );
  }

} // namespace

void MarmotExtendedInterfaceMaterialHypoElastic::computeStress( State&               state,
                                                                Tangents&            tangents,
                                                                const Deformation&   deformation,
                                                                const TimeIncrement& timeIncrement )
{
  const Eigen::Map< const Eigen::Vector3d > normal( deformation.normal.data() );
  const InterfaceKinematics                 kinematics = extractInterfaceKinematics( deformation, h );

  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld + timeIncrement.dT, timeIncrement.dT };

  double alphaOld = *stateLayout.getPtr( state.stateVars, "alpha" );
  if ( !std::isfinite( alphaOld ) || alphaOld <= 0.0 || alphaOld >= 1.0 )
    alphaOld = 0.5;

  Eigen::Vector3d gOld = Eigen::Map< const Eigen::Vector3d >(
    stateLayout.getPtr( state.stateVars, "normalGradientJump" ) );
  if ( !gOld.allFinite() )
    gOld.setZero();

  const bool alphaEvolutionActive = *stateLayout.getPtr( state.stateVars, "alphaEvolutionActive" ) > 0.5;

  LocalSolution solution = solveCoupledLocalProblem( *this,
                                                     stateLayout,
                                                     state.stateVars,
                                                     kinematics.averageNormalGradient,
                                                     kinematics.averageSurfaceGradient,
                                                     kinematics.surfaceGradientJump,
                                                     normal,
                                                     timeInfo,
                                                     alphaOld,
                                                     gOld,
                                                     alphaEvolutionActive );

  MaterialTrial committedTrial = computeMaterialTrial( *this,
                                                       stateLayout,
                                                       state.stateVars,
                                                       kinematics.averageNormalGradient,
                                                       solution.normalGradientJump,
                                                       kinematics.averageSurfaceGradient,
                                                       kinematics.surfaceGradientJump,
                                                       normal,
                                                       timeInfo,
                                                       solution.alpha,
                                                       true,
                                                       true );

  solution.trial = committedTrial;

  Eigen::Map< Eigen::Vector3d > storedNormalGradientJump( stateLayout.getPtr( state.stateVars, "normalGradientJump" ) );
  storedNormalGradientJump                                       = solution.normalGradientJump;
  *stateLayout.getPtr( state.stateVars, "alpha" )                = solution.alpha;
  *stateLayout.getPtr( state.stateVars, "alphaEvolutionActive" ) = ( committedTrial.top.internalStateChanged ||
                                                                     committedTrial.bottom.internalStateChanged )
                                                                     ? 1.0
                                                                     : 0.0;

  const Matrix21dRowMajor tangent  = computeCondensedTangent( solution, normal, h );
  const ExtendedResponse  response = makeExtendedResponse( committedTrial, normal, h );

  writeResponseToState( state, response );
  writeTangentBlocks( tangents, tangent );
}

void MarmotExtendedInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int )
{
  Eigen::Map< Marmot::Vector6d > bottomStress( stateLayout.getPtr( stateVars, "bottomStress" ) );
  Eigen::Map< Marmot::Vector6d > topStress( stateLayout.getPtr( stateVars, "topStress" ) );
  Eigen::Map< Eigen::Vector3d >  normalGradientJump( stateLayout.getPtr( stateVars, "normalGradientJump" ) );

  bottomStress.setZero();
  topStress.setZero();
  normalGradientJump.setZero();
  *stateLayout.getPtr( stateVars, "alpha" )                = 0.5;
  *stateLayout.getPtr( stateVars, "alphaEvolutionActive" ) = 0.0;

  bottomMaterial->initializeYourself( stateLayout.getPtr( stateVars, "bottomMaterialStateVars" ),
                                      bottomMaterial->getNumberOfRequiredStateVars() );

  topMaterial->initializeYourself( stateLayout.getPtr( stateVars, "topMaterialStateVars" ),
                                   topMaterial->getNumberOfRequiredStateVars() );
}

double MarmotExtendedInterfaceMaterialHypoElastic::getDensity()
{
  return 0.5 * ( bottomMaterial->getDensity( nullptr ) + topMaterial->getDensity( nullptr ) );
}
