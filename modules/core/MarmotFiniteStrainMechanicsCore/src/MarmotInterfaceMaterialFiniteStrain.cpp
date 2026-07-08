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

  const Tensor3d  topDisplacement( deformation.dU );
  const Tensor3d  bottomDisplacement( deformation.dU + 3 );
  const Tensor33d topSurfaceGradient( deformation.dSurfaceGradient );
  const Tensor33d bottomSurfaceGradient( deformation.dSurfaceGradient + 9 );
  const Tensor3d  normal( deformation.normal );

  const Tensor3d  jumpU                  = topDisplacement - bottomDisplacement;
  const Tensor33d averageSurfaceGradient = 0.5 * ( topSurfaceGradient + bottomSurfaceGradient );

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
  const Tensor33d P     = einsum< ij, jk, to_ik >( response.tau, FInvT );

  Tensor33d T( 0.0 );
  T.eye();
  T -= einsum< i, j, to_ij >( normal, normal );

  const Tensor3333d A = einsum< imkl, mj, to_ijkl >( materialTangents.dTau_dF, FInvT ) -
                        einsum< im, kj, ml, to_ijkl >( response.tau, FInvT, FInvT );

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

  const Tensor3d    forceTensor         = einsum< ij, j, to_i >( P, normal );
  const Tensor33d   surfaceStressTensor = h * einsum< ij, kj, to_ik >( P, T );
  const Tensor33d   QTensor             = ( 1.0 / h ) * einsum< ijkl, j, l, to_ik >( A, normal, normal );
  const Tensor333d  HJumpAverageTensor  = einsum< ijkl, j, ml, to_ikm >( A, normal, T );
  const Tensor333d  HAverageJumpTensor  = einsum< ijkl, mj, l, to_imk >( A, T, normal );
  const Tensor3333d AAverageTensor      = h * einsum< ijkl, mj, nl, to_imkn >( A, T, T );

  force         = Eigen::Map< const Eigen::Matrix< double, 3, 1 > >( forceTensor.data() );
  surfaceStress = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( surfaceStressTensor.data() );
  Q             = Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( QTensor.data() );
  HJumpAverage  = Eigen::Map< const Eigen::Matrix< double, 3, 9, Eigen::RowMajor > >( HJumpAverageTensor.data() );
  HAverageJump  = Eigen::Map< const Eigen::Matrix< double, 9, 3, Eigen::RowMajor > >( HAverageJumpTensor.data() );
  AAverage      = Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >( AAverageTensor.data() );
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

MarmotInterfaceMaterialFiniteStrain* MarmotLibrary::MarmotInterfaceMaterialFiniteStrainFactory::createMaterial(
  const std::string& materialName,
  const double*      materialProperties,
  int                nMaterialProperties,
  int                materialNumber )
{
  auto& map = materialFactoryFunctionByName();
  auto  it  = map.find( materialName );

  if ( it != map.end() ) {
    return it->second( materialProperties, nMaterialProperties, materialNumber );
  }

  std::string baseMaterialName = materialName;

  return new MarmotInterfaceMaterialFiniteStrain( baseMaterialName,
                                                  materialProperties,
                                                  nMaterialProperties,
                                                  materialNumber );
}
