#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotMaterialFiniteStrainFromGradientEnhanced.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotTesting.h"

#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  // bulk modulus, shear modulus, damage threshold, softening parameter, nonlocal radius, density
  const std::vector< double > damageProperties   = { 3500.0, 1500.0, 1.0e-3, 5.0e-3, 0.1, 1.0 };
  const std::vector< double > neoHookeProperties = { 3500.0, 1500.0 };

  Fastor::Tensor< double, 3, 3 > deformationGradient()
  {
    Fastor::Tensor< double, 3, 3 > F( 0.0 );
    F.eye();
    F( 0, 0 ) += 0.02;
    F( 0, 1 ) += 0.01;
    F( 1, 2 ) -= 0.015;
    F( 2, 2 ) -= 0.01;
    return F;
  }

  /// The wrapped material below its damage threshold is the neo-Hookean material of its bulk.
  void testWrappedMaterialMatchesTheLocalMaterial()
  {
    const std::unique_ptr< MarmotMaterialFiniteStrain > wrapped = MarmotLibrary::
      createFiniteStrainMaterialAllowingGradientEnhanced( "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE",
                                                          damageProperties.data(),
                                                          static_cast< int >( damageProperties.size() ),
                                                          1 );
    const std::unique_ptr< MarmotMaterialFiniteStrain > local(
      MarmotLibrary::MarmotMaterialFiniteStrainFactory::createMaterial( "COMPRESSIBLENEOHOOKE",
                                                                        neoHookeProperties.data(),
                                                                        static_cast< int >( neoHookeProperties.size() ),
                                                                        1 ) );

    std::vector< double > wrappedState( wrapped->getNumberOfRequiredStateVars(), 1.0 );
    wrapped->initializeYourself( wrappedState.data(), static_cast< int >( wrappedState.size() ) );
    throwExceptionOnFailure( wrapped->getNumberOfRequiredStateVars() > 0,
                             "The wrapped material keeps its state variables." );
    for ( double value : wrappedState )
      throwExceptionOnFailure( value == 0.0, "The state variables of the wrapped material are initialized." );

    std::vector< double > localState( local->getNumberOfRequiredStateVars(), 0.0 );

    MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > wrappedResponse( Fastor::Tensor< double, 3, 3 >( 0.0 ),
                                                                           0.0,
                                                                           0.0,
                                                                           wrappedState.data() );
    MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > localResponse( Fastor::Tensor< double, 3, 3 >( 0.0 ),
                                                                         0.0,
                                                                         0.0,
                                                                         localState.data() );
    MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 >    wrappedTangents, localTangents;

    const MarmotMaterialFiniteStrain::Deformation< 3 > deformation{ deformationGradient() };
    const MarmotMaterialFiniteStrain::TimeIncrement    time{ 0.0, 1.0 };

    wrapped->setCharacteristicElementLength( 0.5 );
    wrapped->computeStress( wrappedResponse, wrappedTangents, deformation, time );
    local->computeStress( localResponse, localTangents, deformation, time );

    double stressDifference = 0.0, tangentDifference = 0.0, stressScale = 0.0, tangentScale = 0.0;
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j ) {
        stressDifference = std::max( stressDifference,
                                     std::abs( wrappedResponse.tau( i, j ) - localResponse.tau( i, j ) ) );
        stressScale      = std::max( stressScale, std::abs( localResponse.tau( i, j ) ) );
        for ( int k = 0; k < 3; ++k )
          for ( int l = 0; l < 3; ++l ) {
            tangentDifference = std::max( tangentDifference,
                                          std::abs( wrappedTangents.dTau_dF( i, j, k, l ) -
                                                    localTangents.dTau_dF( i, j, k, l ) ) );
            tangentScale      = std::max( tangentScale, std::abs( localTangents.dTau_dF( i, j, k, l ) ) );
          }
      }

    throwExceptionOnFailure( stressScale > 0.0 && stressDifference < 1e-12 * stressScale,
                             "The stress differs from the local material." );
    throwExceptionOnFailure( tangentDifference < 1e-12 * tangentScale, "The tangent differs from the local material." );
    throwExceptionOnFailure( std::abs( wrapped->getDensity( wrappedState.data() ) - 1.0 ) < 1e-14,
                             "The density is the one of the wrapped material." );
  }

  /// A nonlocal field above the damage threshold makes the stress depend on it, which the adapter must refuse.
  void testMaterialWhoseStressDependsOnTheNonlocalFieldIsRefused()
  {
    auto gradientEnhanced = std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain >(
      MarmotLibrary::MarmotMaterialGradientEnhancedFiniteStrainFactory::
        createMaterial( "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE",
                        damageProperties.data(),
                        static_cast< int >( damageProperties.size() ),
                        1 ) );

    MarmotMaterialFiniteStrainFromGradientEnhanced adapter( std::move( gradientEnhanced ),
                                                            damageProperties.data(),
                                                            static_cast< int >( damageProperties.size() ),
                                                            1,
                                                            3.0e-3 );

    std::vector< double > state( adapter.getNumberOfRequiredStateVars(), 0.0 );
    adapter.initializeYourself( state.data(), static_cast< int >( state.size() ) );

    MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > response( Fastor::Tensor< double, 3, 3 >( 0.0 ),
                                                                    0.0,
                                                                    0.0,
                                                                    state.data() );
    MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 >    tangents;

    bool refused = false;
    try {
      adapter.computeStress( response,
                             tangents,
                             MarmotMaterialFiniteStrain::Deformation< 3 >{ deformationGradient() },
                             { 0.0, 1.0 } );
    }
    catch ( const std::invalid_argument& ) {
      refused = true;
    }
    throwExceptionOnFailure( refused, "A material whose stress depends on the nonlocal field must be refused." );
  }

  void testInvalidMaterialsAreRefused()
  {
    bool noMaterial = false;
    try {
      MarmotMaterialFiniteStrainFromGradientEnhanced( nullptr, neoHookeProperties.data(), 2, 1 );
    }
    catch ( const std::invalid_argument& ) {
      noMaterial = true;
    }
    throwExceptionOnFailure( noMaterial, "An adapter without a material is refused." );

    bool unknown = false;
    try {
      MarmotLibrary::createFiniteStrainMaterialAllowingGradientEnhanced( "NOSUCHMATERIAL",
                                                                         neoHookeProperties.data(),
                                                                         2,
                                                                         1 );
    }
    catch ( const std::invalid_argument& error ) {
      unknown = std::string( error.what() ).find( "NOSUCHMATERIAL" ) != std::string::npos;
    }
    throwExceptionOnFailure( unknown, "An unknown material is reported by name." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testWrappedMaterialMatchesTheLocalMaterial,
                                                   testMaterialWhoseStressDependsOnTheNonlocalFieldIsRefused,
                                                   testInvalidMaterialsAreRefused };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
