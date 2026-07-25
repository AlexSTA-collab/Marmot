#include "Marmot/MarmotXInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

using namespace Marmot;
using namespace Marmot::FastorStandardTensors;

namespace {

  using Vector3d          = Eigen::Matrix< double, 3, 1 >;
  using Vector9d          = Eigen::Matrix< double, 9, 1 >;
  using Matrix3dRowMajor  = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor  = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;
  using Matrix3x9RowMajor = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;

  constexpr int flatIndex( int i, int j )
  {
    return 3 * i + j;
  }

  /**
   * Convert the 6x6 algorithmic tangent returned by a hypoelastic material
   * into the 9x9 tangent mapping a full displacement-gradient increment to
   * the symmetric Cauchy-stress increment. Identical construction to
   * MarmotCorrectedInterfaceMaterialHypoElastic's helper of the same name.
   */
  Matrix9dRowMajor fullGradientTangent( const Matrix6d& tangentVoigt )
  {
    Matrix9dRowMajor tangentFull = Matrix9dRowMajor::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        Matrix3dRowMajor dGradient = Matrix3dRowMajor::Zero();
        dGradient( k, l )          = 1.0;

        const Matrix3dRowMajor dStrain      = 0.5 * ( dGradient + dGradient.transpose() );
        const Vector6d         dStrainVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt( dStrain );
        const Vector6d         dStressVoigt = tangentVoigt * dStrainVoigt;
        const Eigen::Matrix3d  dStress      = ContinuumMechanics::VoigtNotation::voigtToStress( dStressVoigt );

        const int column = flatIndex( k, l );
        for ( int i = 0; i < 3; ++i ) {
          for ( int j = 0; j < 3; ++j ) {
            tangentFull( flatIndex( i, j ), column ) = dStress( i, j );
          }
        }
      }
    }

    return tangentFull;
  }

  struct InterfaceGeometry {
    double   normalSeparation;
    Vector3d tangentialSeparation;
  };

  InterfaceGeometry evaluateInterfaceGeometry( const Vector3d& normal,
                                               const Vector3d& separationVector,
                                               double          constitutiveThickness )
  {
    constexpr double tolerance = 1.0e-12;

    if ( constitutiveThickness <= 0.0 ) {
      throw std::invalid_argument( "MarmotXInterfaceMaterialHypoElastic: interface thickness h must be positive." );
    }

    if ( separationVector.norm() <= tolerance ) {
      return { constitutiveThickness, Vector3d::Zero() };
    }

    const double normalSeparation = separationVector.dot( normal );
    if ( normalSeparation <= tolerance ) {
      throw std::invalid_argument(
        "MarmotXInterfaceMaterialHypoElastic: the top-bottom connector must have a positive normal component." );
    }

    return { normalSeparation, separationVector - normalSeparation * normal };
  }

  /** Maps the raw mesh jump increment to the reconstructed full gradient. */
  Matrix9x3RowMajor jumpToGradient( const Vector3d& normal, double normalSeparation )
  {
    Matrix9x3RowMajor map = Matrix9x3RowMajor::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        map( flatIndex( k, l ), k ) = normal( l ) / normalSeparation;
      }
    }

    return map;
  }

  /**
   * Maps one side's surface-gradient increment A^s to the reconstructed full
   * gradient
   *
   *   G^s = A^s + 1/ell ( [u] - A^s d_tau ) \otimes n.
   */
  Matrix9dRowMajor surfaceToGradient( const Vector3d& normal,
                                      const Vector3d& tangentialSeparation,
                                      double          normalSeparation )
  {
    Matrix9dRowMajor map = Matrix9dRowMajor::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        const int row = flatIndex( k, l );

        for ( int a = 0; a < 3; ++a ) {
          for ( int b = 0; b < 3; ++b ) {
            if ( k != a )
              continue;

            const double identityPart     = l == b ? 1.0 : 0.0;
            map( row, flatIndex( a, b ) ) = identityPart - tangentialSeparation( b ) * normal( l ) / normalSeparation;
          }
        }
      }
    }

    return map;
  }

  /** Maps one side's Cauchy stress to its generalized force conjugate to [u]. */
  Matrix3x9RowMajor stressToForce( const Vector3d& normal, double sideThickness, double normalSeparation )
  {
    Matrix3x9RowMajor map   = Matrix3x9RowMajor::Zero();
    const double      scale = sideThickness / normalSeparation;

    for ( int i = 0; i < 3; ++i ) {
      for ( int b = 0; b < 3; ++b ) {
        map( i, flatIndex( i, b ) ) = scale * normal( b );
      }
    }

    return map;
  }

  /**
   * Maps one side's Cauchy stress to its generalized surface resultant
   * conjugate to A^s:
   *
   *   S^s = (h/2) sigma^s - force^s \otimes d_tau.
   */
  Matrix9dRowMajor stressToSurfaceResultant( const Vector3d& normal,
                                             const Vector3d& tangentialSeparation,
                                             double          sideThickness,
                                             double          normalSeparation )
  {
    Matrix9dRowMajor map        = Matrix9dRowMajor::Zero();
    const double     forceScale = sideThickness / normalSeparation;

    for ( int i = 0; i < 3; ++i ) {
      for ( int p = 0; p < 3; ++p ) {
        const int row = flatIndex( i, p );

        for ( int a = 0; a < 3; ++a ) {
          for ( int b = 0; b < 3; ++b ) {
            if ( i != a )
              continue;

            const double directPart       = p == b ? sideThickness : 0.0;
            const double correction       = forceScale * tangentialSeparation( p ) * normal( b );
            map( row, flatIndex( a, b ) ) = directPart - correction;
          }
        }
      }
    }

    return map;
  }

  /** All generalized outputs of one side's (+ or -) constitutive update. */
  struct OneSideResult {
    Vector3d          force;
    Matrix3dRowMajor  surfaceStress;
    Matrix3dRowMajor  Q;
    Matrix3x9RowMajor H;
    Matrix9x3RowMajor K;
    Matrix9dRowMajor  Z;
  };

  /**
   * Evaluate one side's (+ or -) constitutive response, integrated over
   * half the constitutive thickness. Mirrors
   * MarmotCorrectedInterfaceMaterialHypoElastic::computeStress exactly,
   * with the averaged surface gradient replaced by this side's own A^s and
   * the full thickness h replaced by the half-thickness sideThickness.
   */
  OneSideResult updateOneSide( MarmotMaterialHypoElastic&                 sideMaterial,
                               double*                                    sideMaterialStateVars,
                               const Vector3d&                            normal,
                               const Vector3d&                            dTangential,
                               double                                     ell,
                               double                                     sideThickness,
                               const Matrix3dRowMajor&                    A,
                               const Vector3d&                            jumpIncrement,
                               const Vector3d&                            forceCurrent,
                               const Matrix3dRowMajor&                    surfaceStressCurrent,
                               const MarmotMaterialHypoElastic::timeInfo& timeInfo )
  {
    const Vector3d         q                    = jumpIncrement - A * dTangential;
    const Matrix3dRowMajor G                    = A + ( q / ell ) * normal.transpose();
    const Matrix3dRowMajor strainIncrement      = 0.5 * ( G + G.transpose() );
    const Vector6d         strainIncrementVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt( strainIncrement );

    // Recover this side's current Cauchy stress from
    //   surfaceStressCurrent = sideThickness * sigma - forceCurrent \otimes d_tau.
    Matrix3dRowMajor stressCurrent = ( surfaceStressCurrent + forceCurrent * dTangential.transpose() ) / sideThickness;
    stressCurrent                  = 0.5 * ( stressCurrent + stressCurrent.transpose() );
    const Eigen::Matrix3d stressCurrentSym( stressCurrent );
    const Vector6d        stressVoigt = ContinuumMechanics::VoigtNotation::stressToVoigt( stressCurrentSym );

    Matrix6d                           tangentVoigt = Matrix6d::Zero();
    MarmotMaterialHypoElastic::state3D baseState{ stressVoigt, 0.0, 0.0, sideMaterialStateVars };

    sideMaterial.computeStress( baseState, tangentVoigt, strainIncrementVoigt, timeInfo );

    const Eigen::Matrix3d  stressUpdatedEigen = ContinuumMechanics::VoigtNotation::voigtToStress( baseState.stress );
    const Matrix3dRowMajor stressUpdated      = stressUpdatedEigen;
    const Eigen::Map< const Vector9d > stressVector( stressUpdated.data() );

    const Matrix9x3RowMajor BJump    = jumpToGradient( normal, ell );
    const Matrix9dRowMajor  BSurface = surfaceToGradient( normal, dTangential, ell );
    const Matrix3x9RowMajor RForce   = stressToForce( normal, sideThickness, ell );
    const Matrix9dRowMajor  RSurface = stressToSurfaceResultant( normal, dTangential, sideThickness, ell );
    const Matrix9dRowMajor  CFull    = fullGradientTangent( tangentVoigt );

    OneSideResult result;
    result.force                             = RForce * stressVector;
    const Vector9d surfaceStressVectorResult = RSurface * stressVector;
    result.surfaceStress                     = Eigen::Map< const Matrix3dRowMajor >( surfaceStressVectorResult.data() );
    result.Q                                 = RForce * CFull * BJump;
    result.H                                 = RForce * CFull * BSurface;
    result.K                                 = RSurface * CFull * BJump;
    result.Z                                 = RSurface * CFull * BSurface;
    return result;
  }

} // namespace

MarmotXInterfaceMaterialHypoElastic::MarmotXInterfaceMaterialHypoElastic( const std::string& materialName,
                                                                          const double*      matProperties_,
                                                                          int                nMaterialProperties_,
                                                                          int                materialNumber_ )
  : materialProperties( matProperties_ ), nMaterialProperties( nMaterialProperties_ ), materialNumber( materialNumber_ )
{
  if ( nMaterialProperties < 3 ) {
    throw std::invalid_argument(
      "MarmotXInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h." );
  }

  h = materialProperties[2];
  if ( h <= 0.0 ) {
    throw std::invalid_argument( "MarmotXInterfaceMaterialHypoElastic requires h > 0." );
  }

  baseMaterialProperties.reserve( nMaterialProperties - 1 );
  baseMaterialProperties.push_back( materialProperties[0] );
  baseMaterialProperties.push_back( materialProperties[1] );
  baseMaterialProperties.insert( baseMaterialProperties.end(),
                                 materialProperties + 3,
                                 materialProperties + nMaterialProperties );

  // Top and bottom sides share the same material name and properties but
  // are otherwise fully independent instances with their own history.
  topMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     baseMaterialProperties.data(),
                                                                     static_cast< int >(
                                                                       baseMaterialProperties.size() ),
                                                                     materialNumber ) );
  bottomMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     baseMaterialProperties.data(),
                                                                     static_cast< int >(
                                                                       baseMaterialProperties.size() ),
                                                                     materialNumber ) );

  stateLayout.add( "topMaterialStateVars", topMaterial->getNumberOfRequiredStateVars() );
  stateLayout.add( "bottomMaterialStateVars", bottomMaterial->getNumberOfRequiredStateVars() );
  stateLayout.finalize();
}

void MarmotXInterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  if ( topMaterial ) {
    topMaterial->setCharacteristicElementLength( length );
  }
  if ( bottomMaterial ) {
    bottomMaterial->setCharacteristicElementLength( length );
  }
}

void MarmotXInterfaceMaterialHypoElastic::computeStress( State&               state,
                                                         Tangents&            tangents,
                                                         const Deformation&   deformation,
                                                         const TimeIncrement& timeIncrement )
{
  if ( !topMaterial || !bottomMaterial ) {
    throw std::logic_error( "MarmotXInterfaceMaterialHypoElastic has no base material." );
  }

  const Eigen::Map< const Vector3d > normalMap( deformation.normal.data() );
  const Eigen::Map< const Vector3d > separationMap( deformation.separationVector.data() );

  Vector3d     normal     = normalMap;
  const double normalNorm = normal.norm();
  if ( normalNorm <= 1.0e-12 ) {
    throw std::invalid_argument( "MarmotXInterfaceMaterialHypoElastic: interface normal is zero." );
  }
  normal /= normalNorm;

  const InterfaceGeometry geometry      = evaluateInterfaceGeometry( normal, separationMap, h );
  const double            ell           = geometry.normalSeparation;
  const Vector3d&         dTangential   = geometry.tangentialSeparation;
  const double            sideThickness = 0.5 * h;

  const Eigen::Map< const Eigen::Matrix< double, 6, 1 > >  dU( deformation.dU.data() );
  const Eigen::Map< const Eigen::Matrix< double, 18, 1 > > dSurfaceGradient( deformation.dSurfaceStrain.data() );

  const Vector3d jumpIncrement = dU.template segment< 3 >( 0 ) - dU.template segment< 3 >( 3 );

  const Eigen::Map< const Matrix3dRowMajor > APlus( dSurfaceGradient.data() );
  const Eigen::Map< const Matrix3dRowMajor > AMinus( dSurfaceGradient.data() + 9 );

  const Eigen::Map< const Vector3d >         forcePlusCurrent( state.forcePlus.data() );
  const Eigen::Map< const Vector3d >         forceMinusCurrent( state.forceMinus.data() );
  const Eigen::Map< const Matrix3dRowMajor > surfaceStressPlusCurrent( state.surfaceStressPlus.data() );
  const Eigen::Map< const Matrix3dRowMajor > surfaceStressMinusCurrent( state.surfaceStressMinus.data() );

  double* topStateVars    = stateLayout.getPtr( state.stateVars, "topMaterialStateVars" );
  double* bottomStateVars = stateLayout.getPtr( state.stateVars, "bottomMaterialStateVars" );

  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld, timeIncrement.dT };

  const OneSideResult plusResult = updateOneSide( *topMaterial,
                                                  topStateVars,
                                                  normal,
                                                  dTangential,
                                                  ell,
                                                  sideThickness,
                                                  APlus,
                                                  jumpIncrement,
                                                  forcePlusCurrent,
                                                  surfaceStressPlusCurrent,
                                                  timeInfo );

  const OneSideResult minusResult = updateOneSide( *bottomMaterial,
                                                   bottomStateVars,
                                                   normal,
                                                   dTangential,
                                                   ell,
                                                   sideThickness,
                                                   AMinus,
                                                   jumpIncrement,
                                                   forceMinusCurrent,
                                                   surfaceStressMinusCurrent,
                                                   timeInfo );

  Eigen::Map< Vector3d >         forcePlus( state.forcePlus.data() );
  Eigen::Map< Vector3d >         forceMinus( state.forceMinus.data() );
  Eigen::Map< Matrix3dRowMajor > surfaceStressPlus( state.surfaceStressPlus.data() );
  Eigen::Map< Matrix3dRowMajor > surfaceStressMinus( state.surfaceStressMinus.data() );

  forcePlus          = plusResult.force;
  forceMinus         = minusResult.force;
  surfaceStressPlus  = plusResult.surfaceStress;
  surfaceStressMinus = minusResult.surfaceStress;

  Eigen::Map< Matrix3dRowMajor >  Qplus( tangents.Q_plus.data() );
  Eigen::Map< Matrix3dRowMajor >  Qminus( tangents.Q_minus.data() );
  Eigen::Map< Matrix3x9RowMajor > Hplus( tangents.H_plus.data() );
  Eigen::Map< Matrix3x9RowMajor > Hminus( tangents.H_minus.data() );
  Eigen::Map< Matrix9x3RowMajor > Kplus( tangents.K_plus.data() );
  Eigen::Map< Matrix9x3RowMajor > Kminus( tangents.K_minus.data() );
  Eigen::Map< Matrix9dRowMajor >  Zplus( tangents.Z_plus.data() );
  Eigen::Map< Matrix9dRowMajor >  Zminus( tangents.Z_minus.data() );

  Qplus  = plusResult.Q;
  Hplus  = plusResult.H;
  Kplus  = plusResult.K;
  Zplus  = plusResult.Z;
  Qminus = minusResult.Q;
  Hminus = minusResult.H;
  Kminus = minusResult.K;
  Zminus = minusResult.Z;
}

void MarmotXInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int nStateVars )
{
  if ( !topMaterial || !bottomMaterial ) {
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }
    return;
  }

  topMaterial->initializeYourself( stateLayout.getPtr( stateVars, "topMaterialStateVars" ),
                                   topMaterial->getNumberOfRequiredStateVars() );
  bottomMaterial->initializeYourself( stateLayout.getPtr( stateVars, "bottomMaterialStateVars" ),
                                      bottomMaterial->getNumberOfRequiredStateVars() );
}

double MarmotXInterfaceMaterialHypoElastic::getDensity()
{
  if ( !topMaterial ) {
    return -1;
  }

  return topMaterial->getDensity( nullptr );
}
