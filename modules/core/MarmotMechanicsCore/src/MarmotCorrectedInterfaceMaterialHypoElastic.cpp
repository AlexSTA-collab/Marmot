#include "Marmot/MarmotCorrectedInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

using namespace Marmot;
using namespace Marmot::FastorStandardTensors;

namespace {

  using Vector3d          = Eigen::Matrix< double, 3, 1 >;
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
   * the symmetric Cauchy-stress increment.
   *
   * This construction uses Marmot's own strain/stress Voigt conversions, so
   * the shear convention is inherited correctly. No major symmetry of the
   * 6x6 algorithmic tangent is assumed.
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
      throw std::invalid_argument(
        "MarmotCorrectedInterfaceMaterialHypoElastic: interface thickness h must be positive." );
    }

    if ( separationVector.norm() <= tolerance ) {
      return { constitutiveThickness, Vector3d::Zero() };
    }

    const double normalSeparation = separationVector.dot( normal );
    if ( normalSeparation <= tolerance ) {
      throw std::invalid_argument(
        "MarmotCorrectedInterfaceMaterialHypoElastic: the top-bottom connector must have a positive normal "
        "component." );
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
   * Maps the average surface-gradient increment A to the reconstructed full
   * gradient
   *
   *   G = A + 1/ell ( [u] - A d_tau ) \otimes n.
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

  /** Maps the Cauchy stress to the generalized force conjugate to [u]. */
  Matrix3x9RowMajor stressToForce( const Vector3d& normal, double constitutiveThickness, double normalSeparation )
  {
    Matrix3x9RowMajor map   = Matrix3x9RowMajor::Zero();
    const double      scale = constitutiveThickness / normalSeparation;

    for ( int i = 0; i < 3; ++i ) {
      for ( int b = 0; b < 3; ++b ) {
        map( i, flatIndex( i, b ) ) = scale * normal( b );
      }
    }

    return map;
  }

  /**
   * Maps the Cauchy stress to the generalized surface resultant conjugate to A:
   *
   *   S = h sigma - force \otimes d_tau.
   */
  Matrix9dRowMajor stressToSurfaceResultant( const Vector3d& normal,
                                             const Vector3d& tangentialSeparation,
                                             double          constitutiveThickness,
                                             double          normalSeparation )
  {
    Matrix9dRowMajor map        = Matrix9dRowMajor::Zero();
    const double     forceScale = constitutiveThickness / normalSeparation;

    for ( int i = 0; i < 3; ++i ) {
      for ( int p = 0; p < 3; ++p ) {
        const int row = flatIndex( i, p );

        for ( int a = 0; a < 3; ++a ) {
          for ( int b = 0; b < 3; ++b ) {
            if ( i != a )
              continue;

            const double directPart       = p == b ? constitutiveThickness : 0.0;
            const double correction       = forceScale * tangentialSeparation( p ) * normal( b );
            map( row, flatIndex( a, b ) ) = directPart - correction;
          }
        }
      }
    }

    return map;
  }

} // namespace

MarmotCorrectedInterfaceMaterialHypoElastic::MarmotCorrectedInterfaceMaterialHypoElastic(
  const std::string& materialName,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : materialProperties( matProperties_ ), nMaterialProperties( nMaterialProperties_ ), materialNumber( materialNumber_ )
{
  if ( nMaterialProperties < 3 ) {
    throw std::invalid_argument(
      "MarmotCorrectedInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h." );
  }

  h = materialProperties[2];
  if ( h <= 0.0 ) {
    throw std::invalid_argument( "MarmotCorrectedInterfaceMaterialHypoElastic requires h > 0." );
  }

  baseMaterialProperties.reserve( nMaterialProperties - 1 );
  baseMaterialProperties.push_back( materialProperties[0] );
  baseMaterialProperties.push_back( materialProperties[1] );
  baseMaterialProperties.insert( baseMaterialProperties.end(),
                                 materialProperties + 3,
                                 materialProperties + nMaterialProperties );

  baseMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     baseMaterialProperties.data(),
                                                                     static_cast< int >(
                                                                       baseMaterialProperties.size() ),
                                                                     materialNumber ) );

  stateLayout.add( "baseMaterialStateVars", baseMaterial->getNumberOfRequiredStateVars() );
  stateLayout.finalize();
}

void MarmotCorrectedInterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  if ( baseMaterial ) {
    baseMaterial->setCharacteristicElementLength( length );
  }
}

void MarmotCorrectedInterfaceMaterialHypoElastic::computeStress( State&               state,
                                                                 Tangents&            tangents,
                                                                 const Deformation&   deformation,
                                                                 const TimeIncrement& timeIncrement )
{
  if ( !baseMaterial ) {
    throw std::logic_error( "MarmotCorrectedInterfaceMaterialHypoElastic has no base material." );
  }

  Eigen::Map< Eigen::Vector3d >  force( state.force.data() );
  Eigen::Map< Matrix3dRowMajor > surfaceResultant( state.surfaceStress.data() );

  Eigen::Map< Matrix3dRowMajor >  Q( tangents.Q_ij.data() );
  Eigen::Map< Matrix9dRowMajor >  Z( tangents.Z_ijkl.data() );
  Eigen::Map< Matrix3x9RowMajor > H( tangents.H_ijk.data() );
  Eigen::Map< Matrix9x3RowMajor > K( tangents.K_ijk.data() );

  const Eigen::Map< const Vector3d > normalMap( deformation.normal.data() );
  const Eigen::Map< const Vector3d > separationMap( deformation.separationVector.data() );

  Vector3d     normal     = normalMap;
  const double normalNorm = normal.norm();
  if ( normalNorm <= 1.0e-12 ) {
    throw std::invalid_argument( "MarmotCorrectedInterfaceMaterialHypoElastic: interface normal is zero." );
  }
  normal /= normalNorm;

  const InterfaceGeometry geometry = evaluateInterfaceGeometry( normal, separationMap, h );
  const double            ell      = geometry.normalSeparation;

  // Diagnostic hook: MARMOT_CORRECTED_IFACE_DTAU_SCALE scales the tangential
  // connector d_tau in ALL of its appearances (kinematic reconstruction,
  // generalized-stress maps, tangent blocks, and state recovery), i.e. a
  // variationally consistent lambda-scaling of the geometric coupling.
  // lambda = 1 (default) is the physical model; 0 disables the coupling;
  // -1 reverses d_tau.
  static const double dTauScale = []() {
    const char* s = std::getenv( "MARMOT_CORRECTED_IFACE_DTAU_SCALE" );
    return s ? std::atof( s ) : 1.0;
  }();
  const Vector3d dTangential = dTauScale * geometry.tangentialSeparation;

  const Eigen::Map< const Eigen::Matrix< double, 6, 1 > >  dU( deformation.dU.data() );
  const Eigen::Map< const Eigen::Matrix< double, 18, 1 > > dSurfaceGradient( deformation.dSurfaceStrain.data() );

  const Vector3d jumpIncrement = dU.template segment< 3 >( 0 ) - dU.template segment< 3 >( 3 );

  const Eigen::Map< const Matrix3dRowMajor > topSurfaceGradient( dSurfaceGradient.data() );
  const Eigen::Map< const Matrix3dRowMajor > bottomSurfaceGradient( dSurfaceGradient.data() + 9 );
  const Matrix3dRowMajor averageSurfaceGradient = 0.5 * ( topSurfaceGradient + bottomSurfaceGradient );

  const Vector3d         correctedNormalJump  = jumpIncrement - averageSurfaceGradient * dTangential;
  const Matrix3dRowMajor displacementGradient = averageSurfaceGradient +
                                                ( correctedNormalJump / ell ) * normal.transpose();
  const Matrix3dRowMajor strainIncrement = 0.5 * ( displacementGradient + displacementGradient.transpose() );

  const Vector6d strainIncrementVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt( strainIncrement );

  // The stored interface state is generalized. Recover the Cauchy stress from
  //   surfaceResultant = h sigma - force \otimes d_tau.
  Matrix3dRowMajor stressCurrent = ( surfaceResultant + force * dTangential.transpose() ) / h;
  stressCurrent                  = 0.5 * ( stressCurrent + stressCurrent.transpose() );
  const Eigen::Matrix3d stressCurrentSym( stressCurrent );
  const Vector6d        stressVoigt = ContinuumMechanics::VoigtNotation::stressToVoigt( stressCurrentSym );

  Matrix6d tangentVoigt          = Matrix6d::Zero();
  double*  baseMaterialStateVars = stateLayout.getPtr( state.stateVars, "baseMaterialStateVars" );
  MarmotMaterialHypoElastic::state3D        baseState{ stressVoigt, 0.0, 0.0, baseMaterialStateVars };
  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld, timeIncrement.dT };

  baseMaterial->computeStress( baseState, tangentVoigt, strainIncrementVoigt, timeInfo );

  const Eigen::Matrix3d  stressUpdatedEigen = ContinuumMechanics::VoigtNotation::voigtToStress( baseState.stress );
  const Matrix3dRowMajor stressUpdated      = stressUpdatedEigen;
  const Eigen::Map< const Eigen::Matrix< double, 9, 1 > > stressVector( stressUpdated.data() );

  const Matrix9x3RowMajor BJump    = jumpToGradient( normal, ell );
  const Matrix9dRowMajor  BSurface = surfaceToGradient( normal, dTangential, ell );
  const Matrix3x9RowMajor RForce   = stressToForce( normal, h, ell );
  const Matrix9dRowMajor  RSurface = stressToSurfaceResultant( normal, dTangential, h, ell );
  const Matrix9dRowMajor  CFull    = fullGradientTangent( tangentVoigt );

  force                                                      = RForce * stressVector;
  const Eigen::Matrix< double, 9, 1 > surfaceResultantVector = RSurface * stressVector;
  surfaceResultant = Eigen::Map< const Matrix3dRowMajor >( surfaceResultantVector.data() );

  // Four independent tangent blocks. No major symmetry is used.
  Q = RForce * CFull * BJump;
  H = RForce * CFull * BSurface;
  K = RSurface * CFull * BJump;
  Z = RSurface * CFull * BSurface;
}

void MarmotCorrectedInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int nStateVars )
{
  if ( !baseMaterial ) {
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }
    return;
  }

  baseMaterial->initializeYourself( stateLayout.getPtr( stateVars, "baseMaterialStateVars" ),
                                    baseMaterial->getNumberOfRequiredStateVars() );
}

double MarmotCorrectedInterfaceMaterialHypoElastic::getDensity()
{
  if ( !baseMaterial ) {
    return -1;
  }

  return baseMaterial->getDensity( nullptr );
}
