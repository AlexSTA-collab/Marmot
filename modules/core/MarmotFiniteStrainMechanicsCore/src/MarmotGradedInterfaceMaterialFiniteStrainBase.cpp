#include "Marmot/MarmotGradedInterfaceMaterialFiniteStrainBase.h"

#include "Marmot/MarmotInterfaceFiniteStrainHelperFunctions.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotThicknessProfile.h"
#include "Marmot/MarmotUtils.h"

#include <algorithm>
#include <stdexcept>

using namespace Marmot;
using namespace Marmot::FastorStandardTensors;

MarmotGradedInterfaceMaterialFiniteStrainBase::MarmotGradedInterfaceMaterialFiniteStrainBase(
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

    point.material = std::unique_ptr< MarmotMaterialFiniteStrain >(
      MarmotLibrary::MarmotMaterialFiniteStrainFactory::createMaterial( materialName,
                                                                        point.properties.data(),
                                                                        static_cast< int >( point.properties.size() ),
                                                                        materialNumber ) );
    if ( !point.material )
      throw std::invalid_argument( "Unknown base material for a graded finite-strain interface material: " +
                                   materialName );

    point.nBaseStateVars = point.material->getNumberOfRequiredStateVars();
    point.offset         = offset;
    offset += 3 + point.nBaseStateVars;
    points.push_back( std::move( point ) );
  }

  for ( std::size_t q = 0; q < points.size(); q++ ) {
    stateLayout.add( "correction" + std::to_string( q ), 3 );
    stateLayout.add( "baseMaterialStateVars" + std::to_string( q ), points[q].nBaseStateVars );
  }
  stateLayout.finalize();
}

void MarmotGradedInterfaceMaterialFiniteStrainBase::setCharacteristicElementLength( double length )
{
  for ( auto& point : points )
    point.material->setCharacteristicElementLength( length );
}

MarmotGradedInterfaceMaterialFiniteStrainBase::BulkResponse MarmotGradedInterfaceMaterialFiniteStrainBase::evaluateBulk(
  int                                              q,
  const double*                                    stateVars,
  const Tensor33d&                                 F,
  const MarmotMaterialFiniteStrain::TimeIncrement& timeIncrement ) const
{
  const auto& point = points[q];

  BulkResponse bulk;
  bulk.state.assign( stateVars + point.offset, stateVars + point.offset + 3 + point.nBaseStateVars );

  MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, bulk.state.data() + 3 );
  MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 >    moduli;
  const MarmotMaterialFiniteStrain::Deformation< 3 >    deformation{ F };

  point.material->computeStress( response, moduli, deformation, timeIncrement );

  Marmot::Materials::InterfaceFiniteStrainHelperFunctions::firstPiolaKirchhoff( response.tau,
                                                                                moduli.dTau_dF,
                                                                                F,
                                                                                bulk.P,
                                                                                bulk.A );
  return bulk;
}

void MarmotGradedInterfaceMaterialFiniteStrainBase::initializeYourself( double* stateVars, int )
{
  for ( auto& point : points ) {
    std::fill( stateVars + point.offset, stateVars + point.offset + 3, 0.0 );
    point.material->initializeYourself( stateVars + point.offset + 3, static_cast< int >( point.nBaseStateVars ) );
  }
}

double MarmotGradedInterfaceMaterialFiniteStrainBase::getDensity() const
{
  double density = 0.;
  for ( const auto& point : points )
    density += point.weight * point.material->getDensity( nullptr );
  return density;
}
