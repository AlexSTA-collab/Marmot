#include "Marmot/MarmotExtendedInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

using Matrix3dRowMajor  = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
using Matrix9dRowMajor  = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
using Matrix3x9RowMajor = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;
using Matrix9x3RowMajor = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;

namespace {

  using Vector9d = Eigen::Matrix< double, 9, 1 >;

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

  struct MaterialTrial {
    Marmot::Vector6d  topStress;
    Marmot::Vector6d  bottomStress;
    Marmot::Matrix6d  topTangent;
    Marmot::Matrix6d  bottomTangent;
    Matrix3dRowMajor  topStressTensor;
    Matrix3dRowMajor  bottomStressTensor;
    Matrix3dRowMajor  averageStressTensor;
    Matrix3dRowMajor  jumpStressTensor;
    Matrix3dRowMajor  averageQ;
    Matrix3dRowMajor  jumpQ;
    Matrix3x9RowMajor averageH;
    Matrix3x9RowMajor jumpH;
  };

  struct ExtendedResponse {
    Eigen::Vector3d force;
    Vector9d        averageSurfaceStress;
    Vector9d        jumpSurfaceStress;
  };

  struct NormalGradientJumpTangents {
    Eigen::Matrix3d   dNormalGradientJumpDJumpU;
    Matrix3x9RowMajor dNormalGradientJumpDAverageSurfaceGradient;
    Matrix3x9RowMajor dNormalGradientJumpDJumpSurfaceGradient;
  };

  using Matrix9x21RowMajor = Eigen::Matrix< double, 9, 21, Eigen::RowMajor >;
  using Matrix21dRowMajor  = Eigen::Matrix< double, 21, 21, Eigen::RowMajor >;

} // namespace

MarmotExtendedInterfaceMaterialHypoElastic::MarmotExtendedInterfaceMaterialHypoElastic( const std::string& materialName,
                                                                                        const double* matProperties_,
                                                                                        int nMaterialProperties_,
                                                                                        int materialNumber_ )
  : materialProperties( matProperties_ ), nMaterialProperties( nMaterialProperties_ ), materialNumber( materialNumber_ )
{
  if ( nMaterialProperties < 3 ) {
    throw std::invalid_argument( "MarmotExtendedInterfaceMaterialHypoElastic requires material properties." );
  }

  const bool hasExplicitTopBottomLayout = nMaterialProperties >= 5 && isIntegerProperty( materialProperties[1] );

  if ( hasExplicitTopBottomLayout ) {
    h                         = materialProperties[0];
    const int nBottom         = static_cast< int >( std::round( materialProperties[1] ) );
    const int topSizePosition = 2 + nBottom;

    if ( nBottom <= 0 || topSizePosition >= nMaterialProperties ||
         !isIntegerProperty( materialProperties[topSizePosition] ) ) {
      throw std::invalid_argument( "Invalid extended interface material layout. Expected [h, nBottom, "
                                   "bottomProperties..., nTop, topProperties...]." );
    }

    const int nTop = static_cast< int >( std::round( materialProperties[topSizePosition] ) );

    if ( nTop <= 0 || topSizePosition + 1 + nTop != nMaterialProperties ) {
      throw std::invalid_argument( "Invalid extended interface material layout. Expected [h, nBottom, "
                                   "bottomProperties..., nTop, topProperties...]." );
    }

    bottomMaterialProperties.assign( materialProperties + 2, materialProperties + 2 + nBottom );
    topMaterialProperties.assign( materialProperties + topSizePosition + 1,
                                  materialProperties + topSizePosition + 1 + nTop );
  }
  else {
    if ( nMaterialProperties < 3 ) {
      throw std::invalid_argument(
        "Legacy MarmotExtendedInterfaceMaterialHypoElastic layout requires at least E, nu, and h." );
    }

    h = materialProperties[2];
    bottomMaterialProperties.reserve( nMaterialProperties - 1 );
    bottomMaterialProperties.push_back( materialProperties[0] );
    bottomMaterialProperties.push_back( materialProperties[1] );
    bottomMaterialProperties.insert( bottomMaterialProperties.end(),
                                     materialProperties + 3,
                                     materialProperties + nMaterialProperties );
    topMaterialProperties = bottomMaterialProperties;
  }

  bottomMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     bottomMaterialProperties.data(),
                                                                     static_cast< int >(
                                                                       bottomMaterialProperties.size() ),
                                                                     materialNumber ) );
  topMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     topMaterialProperties.data(),
                                                                     static_cast< int >( topMaterialProperties.size() ),
                                                                     materialNumber ) );

  if ( !bottomMaterial || !topMaterial ) {
    throw std::invalid_argument( "Unknown base material for MarmotExtendedInterfaceMaterialHypoElastic: " +
                                 materialName );
  }

  stateLayout.add( "bottomStress", 6 );
  stateLayout.add( "topStress", 6 );
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

  MaterialTrial computeMaterialTrial( MarmotExtendedInterfaceMaterialHypoElastic& material,
                                      MarmotStateLayoutDynamic&                   stateLayout,
                                      double*                                     stateVars,
                                      const Eigen::Vector3d&                      averageNormalGradient,
                                      const Eigen::Vector3d&                      normalGradientJump,
                                      const Vector9d&                             averageSurfaceGradient,
                                      const Vector9d&                             surfaceGradientJump,
                                      const Eigen::Vector3d&                      normal,
                                      const MarmotMaterialHypoElastic::timeInfo&  timeInfo,
                                      bool                                        commit )
  {
    auto* bottomStressPtr = stateLayout.getPtr( stateVars, "bottomStress" );
    auto* topStressPtr    = stateLayout.getPtr( stateVars, "topStress" );
    auto* bottomStatePtr  = stateLayout.getPtr( stateVars, "bottomMaterialStateVars" );
    auto* topStatePtr     = stateLayout.getPtr( stateVars, "topMaterialStateVars" );

    const int nBottomStateVars = material.getStateView( "bottomMaterialStateVars", stateVars ).stateSize;
    const int nTopStateVars    = material.getStateView( "topMaterialStateVars", stateVars ).stateSize;

    std::vector< double > bottomStateVarsCopy;
    std::vector< double > topStateVarsCopy;

    if ( !commit ) {
      bottomStateVarsCopy.assign( bottomStatePtr, bottomStatePtr + nBottomStateVars );
      topStateVarsCopy.assign( topStatePtr, topStatePtr + nTopStateVars );
      bottomStatePtr = bottomStateVarsCopy.data();
      topStatePtr    = topStateVarsCopy.data();
    }

    const Matrix3dRowMajor averageSurfaceGradientTensor = vectorToTensor( averageSurfaceGradient );
    const Matrix3dRowMajor surfaceGradientJumpTensor    = vectorToTensor( surfaceGradientJump );

    Matrix3dRowMajor topDisplacementGradient    = averageSurfaceGradientTensor + 0.5 * surfaceGradientJumpTensor;
    Matrix3dRowMajor bottomDisplacementGradient = averageSurfaceGradientTensor - 0.5 * surfaceGradientJumpTensor;

    topDisplacementGradient += ( averageNormalGradient + 0.5 * normalGradientJump ) * normal.transpose();
    bottomDisplacementGradient += ( averageNormalGradient - 0.5 * normalGradientJump ) * normal.transpose();

    Marmot::Vector6d topStress    = Eigen::Map< const Marmot::Vector6d >( topStressPtr );
    Marmot::Vector6d bottomStress = Eigen::Map< const Marmot::Vector6d >( bottomStressPtr );

    Marmot::Matrix6d topTangent    = Marmot::Matrix6d::Zero();
    Marmot::Matrix6d bottomTangent = Marmot::Matrix6d::Zero();

    MarmotMaterialHypoElastic::state3D topState{ topStress, 0.0, 0.0, topStatePtr };
    MarmotMaterialHypoElastic::state3D bottomState{ bottomStress, 0.0, 0.0, bottomStatePtr };

    material.getTopMaterial().computeStress( topState, topTangent, strainToVoigt( topDisplacementGradient ), timeInfo );
    material.getBottomMaterial().computeStress( bottomState,
                                                bottomTangent,
                                                strainToVoigt( bottomDisplacementGradient ),
                                                timeInfo );

    if ( commit ) {
      Eigen::Map< Marmot::Vector6d > topStressMap( topStressPtr );
      Eigen::Map< Marmot::Vector6d > bottomStressMap( bottomStressPtr );
      topStressMap    = topState.stress;
      bottomStressMap = bottomState.stress;
    }

    const auto normalTensor                         = Marmot::FastorStandardTensors::Tensor3d( normal.data() );
    const auto [topZ, topQTensor, topHTensor, topY] = Marmot::Materials::InterfaceMaterialHelperFunctions::
      calculateInterfaceMaterialParameters( normalTensor, topTangent );
    const auto [bottomZ, bottomQTensor, bottomHTensor, bottomY] = Marmot::Materials::InterfaceMaterialHelperFunctions::
      calculateInterfaceMaterialParameters( normalTensor, bottomTangent );

    MaterialTrial trial;
    trial.topStress           = topState.stress;
    trial.bottomStress        = bottomState.stress;
    trial.topTangent          = topTangent;
    trial.bottomTangent       = bottomTangent;
    trial.topStressTensor     = stressToTensor( topState.stress );
    trial.bottomStressTensor  = stressToTensor( bottomState.stress );
    trial.averageStressTensor = 0.5 * ( trial.topStressTensor + trial.bottomStressTensor );
    trial.jumpStressTensor    = trial.topStressTensor - trial.bottomStressTensor;
    trial.averageQ            = 0.5 * ( Eigen::Map< const Matrix3dRowMajor >( topQTensor.data() ) +
                             Eigen::Map< const Matrix3dRowMajor >( bottomQTensor.data() ) );
    trial.jumpQ               = Eigen::Map< const Matrix3dRowMajor >( topQTensor.data() ) -
                  Eigen::Map< const Matrix3dRowMajor >( bottomQTensor.data() );
    trial.averageH = 0.5 * ( Eigen::Map< const Matrix3x9RowMajor >( topHTensor.data() ) +
                             Eigen::Map< const Matrix3x9RowMajor >( bottomHTensor.data() ) );
    trial.jumpH    = Eigen::Map< const Matrix3x9RowMajor >( topHTensor.data() ) -
                  Eigen::Map< const Matrix3x9RowMajor >( bottomHTensor.data() );

    (void)topZ;
    (void)bottomZ;
    (void)topY;
    (void)bottomY;

    return trial;
  }

  Eigen::Vector3d computeTractionJumpResidual( const MaterialTrial& trial, const Eigen::Vector3d& normal )
  {
    return ( trial.topStressTensor - trial.bottomStressTensor ) * normal;
  }

  Eigen::Vector3d solveNormalGradientJump( MarmotExtendedInterfaceMaterialHypoElastic& material,
                                           MarmotStateLayoutDynamic&                   stateLayout,
                                           double*                                     stateVars,
                                           const Eigen::Vector3d&                      averageNormalGradient,
                                           const Vector9d&                             averageSurfaceGradient,
                                           const Vector9d&                             surfaceGradientJump,
                                           const Eigen::Vector3d&                      normal,
                                           const MarmotMaterialHypoElastic::timeInfo&  timeInfo )
  {
    Eigen::Vector3d normalGradientJump = Eigen::Vector3d::Zero();

    for ( int iteration = 0; iteration < 12; ++iteration ) {
      const auto trial = computeMaterialTrial( material,
                                               stateLayout,
                                               stateVars,
                                               averageNormalGradient,
                                               normalGradientJump,
                                               averageSurfaceGradient,
                                               surfaceGradientJump,
                                               normal,
                                               timeInfo,
                                               false );

      const Eigen::Vector3d residual = computeTractionJumpResidual( trial, normal );

      if ( residual.norm() < 1e-11 * std::max( 1.0, normalGradientJump.norm() ) ) {
        return normalGradientJump;
      }

      const Eigen::Vector3d correction = trial.averageQ.fullPivLu().solve( residual );
      normalGradientJump -= correction;

      if ( correction.norm() < 1e-11 * std::max( 1.0, normalGradientJump.norm() ) ) {
        return normalGradientJump;
      }
    }

    throw std::runtime_error( "MarmotExtendedInterfaceMaterialHypoElastic local Newton iteration failed." );
  }

  NormalGradientJumpTangents computeNormalGradientJumpTangents( const MaterialTrial& trial, double h )
  {
    NormalGradientJumpTangents tangents;
    tangents.dNormalGradientJumpDJumpU                  = -( 1. / h ) * trial.averageQ.fullPivLu().solve( trial.jumpQ );
    tangents.dNormalGradientJumpDAverageSurfaceGradient = -trial.averageQ.fullPivLu().solve( trial.jumpH );
    tangents.dNormalGradientJumpDJumpSurfaceGradient    = -trial.averageQ.fullPivLu().solve( trial.averageH );
    return tangents;
  }

  Vector9d stressVoigtToTensorVector( const Marmot::Vector6d& stressVoigt )
  {
    return tensorToVector( stressToTensor( stressVoigt ) );
  }

  Vector9d stressGradientColumn( const Marmot::Matrix6d& tangent, const Matrix3dRowMajor& dDisplacementGradient )
  {
    return stressVoigtToTensorVector( tangent * strainToVoigt( dDisplacementGradient ) );
  }

  Matrix9x21RowMajor computeStressTensorDerivative( const Marmot::Matrix6d&  tangent,
                                                    const Eigen::Matrix3d&   dNormalGradientDJumpU,
                                                    const Matrix3x9RowMajor& dNormalGradientDAverageSurfaceGradient,
                                                    const Matrix3x9RowMajor& dNormalGradientDJumpSurfaceGradient,
                                                    const Eigen::Vector3d&   normal,
                                                    double                   averageSurfaceGradientSign,
                                                    double                   jumpSurfaceGradientSign )
  {
    Matrix9x21RowMajor derivative = Matrix9x21RowMajor::Zero();

    for ( int column = 0; column < 3; ++column ) {
      Matrix3dRowMajor dDisplacementGradient = Matrix3dRowMajor::Zero();
      dDisplacementGradient += dNormalGradientDJumpU.col( column ) * normal.transpose();
      derivative.col( column ) = stressGradientColumn( tangent, dDisplacementGradient );
    }

    for ( int column = 0; column < 9; ++column ) {
      Matrix3dRowMajor dDisplacementGradient = Matrix3dRowMajor::Zero();
      dDisplacementGradient( column / 3, column % 3 ) += averageSurfaceGradientSign;
      dDisplacementGradient += dNormalGradientDAverageSurfaceGradient.col( column ) * normal.transpose();
      derivative.col( 3 + column ) = stressGradientColumn( tangent, dDisplacementGradient );
    }

    for ( int column = 0; column < 9; ++column ) {
      Matrix3dRowMajor dDisplacementGradient = Matrix3dRowMajor::Zero();
      dDisplacementGradient( column / 3, column % 3 ) += jumpSurfaceGradientSign;
      dDisplacementGradient += dNormalGradientDJumpSurfaceGradient.col( column ) * normal.transpose();
      derivative.col( 12 + column ) = stressGradientColumn( tangent, dDisplacementGradient );
    }

    return derivative;
  }

  Matrix21dRowMajor computeImplicitTangent( const MaterialTrial&              trial,
                                            const NormalGradientJumpTangents& normalGradientJumpTangents,
                                            const Eigen::Vector3d&            normal,
                                            double                            h )
  {
    Eigen::Matrix3d dTopNormalGradientDJumpU = ( 1. / h ) * Eigen::Matrix3d::Identity();
    dTopNormalGradientDJumpU += 0.5 * normalGradientJumpTangents.dNormalGradientJumpDJumpU;

    Eigen::Matrix3d dBottomNormalGradientDJumpU = ( 1. / h ) * Eigen::Matrix3d::Identity();
    dBottomNormalGradientDJumpU -= 0.5 * normalGradientJumpTangents.dNormalGradientJumpDJumpU;

    Matrix3x9RowMajor dTopNormalGradientDAverageSurfaceGradient = 0.5 * normalGradientJumpTangents
                                                                          .dNormalGradientJumpDAverageSurfaceGradient;
    Matrix3x9RowMajor dBottomNormalGradientDAverageSurfaceGradient = -0.5 *
                                                                     normalGradientJumpTangents
                                                                       .dNormalGradientJumpDAverageSurfaceGradient;
    Matrix3x9RowMajor dTopNormalGradientDJumpSurfaceGradient = 0.5 * normalGradientJumpTangents
                                                                       .dNormalGradientJumpDJumpSurfaceGradient;
    Matrix3x9RowMajor dBottomNormalGradientDJumpSurfaceGradient = -0.5 * normalGradientJumpTangents
                                                                           .dNormalGradientJumpDJumpSurfaceGradient;

    const Matrix9x21RowMajor dTopStress = computeStressTensorDerivative( trial.topTangent,
                                                                         dTopNormalGradientDJumpU,
                                                                         dTopNormalGradientDAverageSurfaceGradient,
                                                                         dTopNormalGradientDJumpSurfaceGradient,
                                                                         normal,
                                                                         1.0,
                                                                         0.5 );
    const Matrix9x21RowMajor
      dBottomStress = computeStressTensorDerivative( trial.bottomTangent,
                                                     dBottomNormalGradientDJumpU,
                                                     dBottomNormalGradientDAverageSurfaceGradient,
                                                     dBottomNormalGradientDJumpSurfaceGradient,
                                                     normal,
                                                     1.0,
                                                     -0.5 );

    const Matrix9x21RowMajor dAverageStress = 0.5 * ( dTopStress + dBottomStress );
    const Matrix9x21RowMajor dJumpStress    = dTopStress - dBottomStress;

    Eigen::Matrix< double, 3, 9, Eigen::RowMajor > dForceDAverageStress;
    dForceDAverageStress.setZero();
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 3; ++j ) {
        dForceDAverageStress( i, i * 3 + j ) = normal[j];
      }
    }

    Matrix21dRowMajor tangent       = Matrix21dRowMajor::Zero();
    tangent.block< 3, 21 >( 0, 0 )  = dForceDAverageStress * dAverageStress;
    tangent.block< 9, 21 >( 3, 0 )  = h * dAverageStress;
    tangent.block< 9, 21 >( 12, 0 ) = 0.25 * h * dJumpStress;

    return tangent;
  }

  ExtendedResponse evaluateExtendedResponse( MarmotExtendedInterfaceMaterialHypoElastic& material,
                                             MarmotStateLayoutDynamic&                   stateLayout,
                                             double*                                     stateVars,
                                             double                                      h,
                                             const Eigen::Vector3d&                      displacementJump,
                                             const Vector9d&                             averageSurfaceGradient,
                                             const Vector9d&                             surfaceGradientJump,
                                             const Eigen::Vector3d&                      normal,
                                             const MarmotMaterialHypoElastic::timeInfo&  timeInfo,
                                             bool                                        commit )
  {
    const Eigen::Vector3d averageNormalGradient = ( 1. / h ) * displacementJump;
    const Eigen::Vector3d normalGradientJump    = solveNormalGradientJump( material,
                                                                        stateLayout,
                                                                        stateVars,
                                                                        averageNormalGradient,
                                                                        averageSurfaceGradient,
                                                                        surfaceGradientJump,
                                                                        normal,
                                                                        timeInfo );

    const auto trial = computeMaterialTrial( material,
                                             stateLayout,
                                             stateVars,
                                             averageNormalGradient,
                                             normalGradientJump,
                                             averageSurfaceGradient,
                                             surfaceGradientJump,
                                             normal,
                                             timeInfo,
                                             commit );

    ExtendedResponse response;
    response.force                = trial.averageStressTensor * normal;
    response.averageSurfaceStress = h * tensorToVector( trial.averageStressTensor );
    response.jumpSurfaceStress    = 0.25 * h * tensorToVector( trial.jumpStressTensor );
    return response;
  }

} // namespace

void MarmotExtendedInterfaceMaterialHypoElastic::computeStress( State&               state,
                                                                Tangents&            tangents,
                                                                const Deformation&   deformation,
                                                                const TimeIncrement& timeIncrement )
{
  const Eigen::Map< const Eigen::Vector3d >                normal( deformation.normal.data() );
  const Eigen::Map< const Eigen::Matrix< double, 6, 1 > >  dU( deformation.dU.data() );
  const Eigen::Map< const Eigen::Matrix< double, 18, 1 > > dSurfaceStrain( deformation.dSurfaceStrain.data() );

  const Eigen::Vector3d displacementJump = dU.segment< 3 >( 0 ) - dU.segment< 3 >( 3 );
  const Vector9d averageSurfaceGradient = 0.5 * ( dSurfaceStrain.segment< 9 >( 0 ) + dSurfaceStrain.segment< 9 >( 9 ) );
  const Vector9d surfaceGradientJump    = dSurfaceStrain.segment< 9 >( 0 ) - dSurfaceStrain.segment< 9 >( 9 );

  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld + timeIncrement.dT, timeIncrement.dT };

  const Eigen::Vector3d averageNormalGradient      = ( 1. / h ) * displacementJump;
  const Eigen::Vector3d normalGradientJump         = solveNormalGradientJump( *this,
                                                                      stateLayout,
                                                                      state.stateVars,
                                                                      averageNormalGradient,
                                                                      averageSurfaceGradient,
                                                                      surfaceGradientJump,
                                                                      normal,
                                                                      timeInfo );
  const auto            trial                      = computeMaterialTrial( *this,
                                           stateLayout,
                                           state.stateVars,
                                           averageNormalGradient,
                                           normalGradientJump,
                                           averageSurfaceGradient,
                                           surfaceGradientJump,
                                           normal,
                                           timeInfo,
                                           false );
  const auto            normalGradientJumpTangents = computeNormalGradientJumpTangents( trial, h );
  const auto            tangent = computeImplicitTangent( trial, normalGradientJumpTangents, normal, h );

  const auto committedResponse = evaluateExtendedResponse( *this,
                                                           stateLayout,
                                                           state.stateVars,
                                                           h,
                                                           displacementJump,
                                                           averageSurfaceGradient,
                                                           surfaceGradientJump,
                                                           normal,
                                                           timeInfo,
                                                           true );

  state.force                = Marmot::FastorStandardTensors::Tensor3d( committedResponse.force.data() );
  state.averageSurfaceStress = Marmot::FastorStandardTensors::Tensor33d(
    vectorToTensor( committedResponse.averageSurfaceStress ).data() );
  state.jumpSurfaceStress = Marmot::FastorStandardTensors::Tensor33d(
    vectorToTensor( committedResponse.jumpSurfaceStress ).data() );

  Eigen::Map< Matrix3dRowMajor >( tangents.forceJumpU )                   = tangent.block< 3, 3 >( 0, 0 );
  Eigen::Map< Matrix3x9RowMajor >( tangents.forceAverageSurfaceGradient ) = tangent.block< 3, 9 >( 0, 3 );
  Eigen::Map< Matrix3x9RowMajor >( tangents.forceJumpSurfaceGradient )    = tangent.block< 3, 9 >( 0, 12 );

  Eigen::Map< Matrix9x3RowMajor >( tangents.averageSurfaceStressJumpU )                 = tangent.block< 9, 3 >( 3, 0 );
  Eigen::Map< Matrix9dRowMajor >( tangents.averageSurfaceStressAverageSurfaceGradient ) = tangent.block< 9, 9 >( 3, 3 );
  Eigen::Map< Matrix9dRowMajor >( tangents.averageSurfaceStressJumpSurfaceGradient ) = tangent.block< 9, 9 >( 3, 12 );

  Eigen::Map< Matrix9x3RowMajor >( tangents.jumpSurfaceStressJumpU )                 = tangent.block< 9, 3 >( 12, 0 );
  Eigen::Map< Matrix9dRowMajor >( tangents.jumpSurfaceStressAverageSurfaceGradient ) = tangent.block< 9, 9 >( 12, 3 );
  Eigen::Map< Matrix9dRowMajor >( tangents.jumpSurfaceStressJumpSurfaceGradient )    = tangent.block< 9, 9 >( 12, 12 );
}

void MarmotExtendedInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int )
{
  Eigen::Map< Marmot::Vector6d >( stateLayout.getPtr( stateVars, "bottomStress" ) ).setZero();
  Eigen::Map< Marmot::Vector6d >( stateLayout.getPtr( stateVars, "topStress" ) ).setZero();

  bottomMaterial->initializeYourself( stateLayout.getPtr( stateVars, "bottomMaterialStateVars" ),
                                      bottomMaterial->getNumberOfRequiredStateVars() );
  topMaterial->initializeYourself( stateLayout.getPtr( stateVars, "topMaterialStateVars" ),
                                   topMaterial->getNumberOfRequiredStateVars() );
}

double MarmotExtendedInterfaceMaterialHypoElastic::getDensity()
{
  return 0.5 * ( bottomMaterial->getDensity( nullptr ) + topMaterial->getDensity( nullptr ) );
}
