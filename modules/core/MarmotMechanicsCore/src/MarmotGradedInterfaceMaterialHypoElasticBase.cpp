#include "Marmot/MarmotGradedInterfaceMaterialHypoElasticBase.h"

#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotThicknessProfile.h"
#include "Marmot/MarmotUtils.h"
#include "Marmot/MarmotVoigt.h"

#include <algorithm>
#include <stdexcept>

using namespace Marmot;

MarmotGradedInterfaceMaterialHypoElasticBase::MarmotGradedInterfaceMaterialHypoElasticBase(
  const std::string& materialName,
  const double*      matProperties,
  int                nMaterialProperties,
  int                materialNumber_ )
  : materialNumber( materialNumber_ ), h( checkedMaterialProperty( matProperties, nMaterialProperties, 0 ) )
{
  auto profile = Marmot::Materials::ThicknessProfile::fromMaterialProperties( matProperties, nMaterialProperties );

  tangentialGradient = profile.tangentialGradient;

  cumulativeWeights.resize( profile.points.size(), profile.points.size() );
  for ( std::size_t q = 0; q < profile.points.size(); q++ )
    for ( std::size_t p = 0; p < profile.points.size(); p++ )
      cumulativeWeights( q, p ) = profile.cumulativeWeights[q][p];

  std::size_t offset = 0;
  for ( auto& profilePoint : profile.points ) {
    ThicknessPoint point;
    point.weight     = profilePoint.weight;
    point.position   = profilePoint.position;
    point.properties = std::move( profilePoint.properties );

    point.material = std::unique_ptr< MarmotMaterialHypoElastic >(
      MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                       point.properties.data(),
                                                                       static_cast< int >( point.properties.size() ),
                                                                       materialNumber ) );
    if ( !point.material )
      throw std::invalid_argument( "Unknown base material for a graded interface material: " + materialName );

    point.nBaseStateVars = point.material->getNumberOfRequiredStateVars();
    point.offset         = offset;
    offset += 6 + point.nBaseStateVars;
    points.push_back( std::move( point ) );
  }

  for ( std::size_t q = 0; q < points.size(); q++ ) {
    stateLayout.add( "stress" + std::to_string( q ), 6 );
    stateLayout.add( "baseMaterialStateVars" + std::to_string( q ), points[q].nBaseStateVars );
  }
  stateLayout.finalize();
}

Eigen::Matrix< double, 6, 3 > MarmotGradedInterfaceMaterialHypoElasticBase::normalStrainOperator(
  const Eigen::Vector3d& n )
{
  Eigen::Matrix< double, 6, 3 > S;
  for ( int j = 0; j < 3; j++ ) {
    Eigen::Vector3d e            = Eigen::Vector3d::Zero();
    e( j )                       = 1.0;
    const Eigen::Matrix3d strain = 0.5 * ( e * n.transpose() + n * e.transpose() );
    S.col( j )                   = ContinuumMechanics::VoigtNotation::strainToVoigt( strain );
  }
  return S;
}

Eigen::Matrix< double, 6, 9 > MarmotGradedInterfaceMaterialHypoElasticBase::gradientToVoigtStrain()
{
  Eigen::Matrix< double, 6, 9 > E;
  for ( int k = 0; k < 3; k++ )
    for ( int l = 0; l < 3; l++ ) {
      Eigen::Matrix3d unit         = Eigen::Matrix3d::Zero();
      unit( k, l )                 = 1.0;
      const Eigen::Matrix3d strain = 0.5 * ( unit + unit.transpose() );
      E.col( 3 * k + l )           = ContinuumMechanics::VoigtNotation::strainToVoigt( strain );
    }
  return E;
}

Eigen::Matrix< double, 9, 6 > MarmotGradedInterfaceMaterialHypoElasticBase::voigtStressToVector()
{
  Eigen::Matrix< double, 9, 6 > T;
  for ( int j = 0; j < 6; j++ ) {
    Vector6d unit                = Vector6d::Zero();
    unit( j )                    = 1.0;
    const Eigen::Matrix3d stress = ContinuumMechanics::VoigtNotation::voigtToStress( unit );
    for ( int k = 0; k < 3; k++ )
      for ( int l = 0; l < 3; l++ )
        T( 3 * k + l, j ) = stress( k, l );
  }
  return T;
}

void MarmotGradedInterfaceMaterialHypoElasticBase::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  for ( auto& point : points )
    point.material->setCharacteristicElementLength( length );
}

MarmotGradedInterfaceMaterialHypoElasticBase::BulkResponse MarmotGradedInterfaceMaterialHypoElasticBase::evaluateBulk(
  int                                        q,
  const double*                              stateVars,
  const Vector6d&                            strainIncrement,
  const MarmotMaterialHypoElastic::timeInfo& timeInfo ) const
{
  const auto& point = points[q];

  BulkResponse response;
  response.state.assign( stateVars + point.offset, stateVars + point.offset + 6 + point.nBaseStateVars );

  MarmotMaterialHypoElastic::state3D baseState;
  baseState.stress    = Eigen::Map< const Vector6d >( response.state.data() );
  baseState.stateVars = response.state.data() + 6;

  response.D.setZero();
  point.material->computeStress( baseState, response.D, strainIncrement, timeInfo );

  Eigen::Map< Vector6d > storedStress( response.state.data() );
  storedStress   = baseState.stress;
  response.sigma = baseState.stress;
  return response;
}

void MarmotGradedInterfaceMaterialHypoElasticBase::initializeYourself( double* stateVars, int )
{
  for ( auto& point : points ) {
    std::fill( stateVars + point.offset, stateVars + point.offset + 6, 0.0 );
    point.material->initializeYourself( stateVars + point.offset + 6, static_cast< int >( point.nBaseStateVars ) );
  }
}

double MarmotGradedInterfaceMaterialHypoElasticBase::getDensity()
{
  double density = 0.;
  for ( auto& point : points )
    density += point.weight * point.material->getDensity( nullptr );
  return density;
}
