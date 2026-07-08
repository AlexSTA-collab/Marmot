#include "Marmot/MarmotInterfaceMaterialFiniteStrain.h"

#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

#include <Eigen/Dense>

#include <stdexcept>

using namespace Marmot;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

MarmotInterfaceMaterialFiniteStrain::MarmotInterfaceMaterialFiniteStrain( const std::string& materialName,
                                                                          const double*      materialProperties_,
                                                                          int                nMaterialProperties_,
                                                                          int                materialNumber_ )
  : materialProperties( materialProperties_ ),
    nMaterialProperties( nMaterialProperties_ ),
    materialNumber( materialNumber_ )
{
  if ( nMaterialProperties < 2 ) {
    throw std::invalid_argument(
      "MarmotInterfaceMaterialFiniteStrain requires interface thickness h followed by base material properties." );
  }

  h = materialProperties[0];
  if ( h <= 0.0 ) {
    throw std::invalid_argument( "MarmotInterfaceMaterialFiniteStrain requires positive interface thickness h." );
  }

  baseMaterialProperties.assign( materialProperties + 1, materialProperties + nMaterialProperties );

  baseMaterial = std::unique_ptr< MarmotMaterialFiniteStrain >(
    MarmotLibrary::MarmotMaterialFiniteStrainFactory::createMaterial( materialName,
                                                                      baseMaterialProperties.data(),
                                                                      static_cast< int >(
                                                                        baseMaterialProperties.size() ),
                                                                      materialNumber ) );

  stateLayout.add( "baseMaterialStateVars", baseMaterial->getNumberOfRequiredStateVars() );
  stateLayout.finalize();
}

void MarmotInterfaceMaterialFiniteStrain::computeStress( State&               state,
                                                         Tangents&            tangents,
                                                         const Deformation&   deformation,
                                                         const TimeIncrement& timeIncrement )
{
  using namespace Fastor;

  Tensor33d averageSurfaceGradient( 0.0 );
  Tensor3d  jumpU( 0.0 );
  Tensor3d  normal( deformation.normal );

  for ( int i = 0; i < 3; ++i ) {
    jumpU( i ) = deformation.dU[i] - deformation.dU[3 + i];

    for ( int J = 0; J < 3; ++J ) {
      const int index                = i * 3 + J;
      averageSurfaceGradient( i, J ) = 0.5 * ( deformation.dSurfaceGradient[index] +
                                               deformation.dSurfaceGradient[9 + index] );
    }
  }

  Tensor33d F( 0.0 );
  F.eye();
  F += averageSurfaceGradient + ( 1.0 / h ) * einsum< i, j, to_ij >( jumpU, normal );

  double* baseMaterialStateVars = stateLayout.getPtr( state.stateVars, "baseMaterialStateVars" );

  MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, baseMaterialStateVars );
  MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 >    materialTangents;
  MarmotMaterialFiniteStrain::Deformation< 3 >          materialDeformation{ F };
  const MarmotMaterialFiniteStrain::TimeIncrement       materialTimeIncrement{ timeIncrement.timeOld + timeIncrement.dT,
                                                                         timeIncrement.dT };

  baseMaterial->computeStress( response, materialTangents, materialDeformation, materialTimeIncrement );

  const Tensor33d FInvT = transpose( inverse( F ) );
  Tensor33d       P( 0.0 );
  for ( int i = 0; i < 3; ++i )
    for ( int J = 0; J < 3; ++J )
      for ( int j = 0; j < 3; ++j )
        P( i, J ) += response.tau( i, j ) * FInvT( j, J );

  Tensor33d T( 0.0 );
  T.eye();
  T -= einsum< i, j, to_ij >( normal, normal );

  Tensor3333d A( 0.0 );
  for ( int i = 0; i < 3; ++i ) {
    for ( int J = 0; J < 3; ++J ) {
      for ( int k = 0; k < 3; ++k ) {
        for ( int L = 0; L < 3; ++L ) {
          double value = 0.0;

          for ( int j = 0; j < 3; ++j ) {
            value += materialTangents.dTau_dF( i, j, k, L ) * FInvT( j, J );
            value -= response.tau( i, j ) * FInvT( k, J ) * FInvT( j, L );
          }

          A( i, J, k, L ) = value;
        }
      }
    }
  }

  Eigen::Map< Eigen::Matrix< double, 3, 1 > >                  force( state.force );
  Eigen::Map< Eigen::Matrix< double, 9, 1 > >                  surfaceStress( state.surfaceStress );
  Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > Q( tangents.Q_ik );
  Eigen::Map< Eigen::Matrix< double, 3, 9, Eigen::RowMajor > > HJumpAverage( tangents.H_ikR );
  Eigen::Map< Eigen::Matrix< double, 9, 3, Eigen::RowMajor > > HAverageJump( tangents.H_iRk );
  Eigen::Map< Eigen::Matrix< double, 9, 9, Eigen::RowMajor > > AAverage( tangents.A_iRkP );

  force.setZero();
  surfaceStress.setZero();
  Q.setZero();
  HJumpAverage.setZero();
  HAverageJump.setZero();
  AAverage.setZero();

  for ( int i = 0; i < 3; ++i ) {
    for ( int J = 0; J < 3; ++J ) {
      force( i ) += P( i, J ) * normal( J );

      for ( int R = 0; R < 3; ++R )
        surfaceStress( i * 3 + R ) += h * P( i, J ) * T( R, J );

      for ( int k = 0; k < 3; ++k ) {
        for ( int L = 0; L < 3; ++L ) {
          Q( i, k ) += ( 1.0 / h ) * A( i, J, k, L ) * normal( J ) * normal( L );

          for ( int R = 0; R < 3; ++R ) {
            HJumpAverage( i, k * 3 + R ) += A( i, J, k, L ) * normal( J ) * T( R, L );
            HAverageJump( i * 3 + R, k ) += A( i, J, k, L ) * T( R, J ) * normal( L );

            for ( int Pdir = 0; Pdir < 3; ++Pdir ) {
              AAverage( i * 3 + R, k * 3 + Pdir ) += h * A( i, J, k, L ) * T( R, J ) * T( Pdir, L );
            }
          }
        }
      }
    }
  }
}

void MarmotInterfaceMaterialFiniteStrain::initializeYourself( double* stateVars, int nStateVars )
{
  baseMaterial->initializeYourself( stateLayout.getPtr( stateVars, "baseMaterialStateVars" ),
                                    baseMaterial->getNumberOfRequiredStateVars() );
}

double MarmotInterfaceMaterialFiniteStrain::getDensity() const
{
  return baseMaterial->getDensity( nullptr );
}
