#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotStabPressureXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Material   = MarmotStabPressureXInterfaceMaterialHypoElastic;
  using EqMaterial = MarmotEquilibratedXInterfaceMaterialHypoElastic;

  using Vector21d = Eigen::Matrix< double, 21, 1 >;
  using Vector22d = Eigen::Matrix< double, 22, 1 >;
  using Matrix22d = Eigen::Matrix< double, 22, 22 >;

  const std::vector< double > elasticProperties = { 2e5, 0.3, 0.1 };
  const std::vector< double > plasticProperties = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };

  const double normalZ[3]        = { 0., 0., 1. };
  const double normalTilt[3]     = { std::sin( 0.17 ), 0., std::cos( 0.17 ) };
  const double noSeparation[3]   = { 0., 0., 0. };
  const double tiltSeparation[3] = { 0.1 * std::sin( 0.17 ) + 0.02 * std::cos( 0.17 ),
                                     0.,
                                     0.1 * std::cos( 0.17 ) - 0.02 * std::sin( 0.17 ) };

  /// Generalized stress and volumetric residual of one update, flattened as (f, S+, S-, r_p).
  struct Result {
    Vector22d response = Vector22d::Zero();
    Matrix22d tangent  = Matrix22d::Zero(); ///< rows: (f, S+, S-, r_p), columns: (w, A+, A-, p)
  };

  Result update( Material&              material,
                 std::vector< double >& stateVars,
                 const Vector21d&       x,
                 double                 dPressure,
                 const double*          normal,
                 const double*          separation,
                 double                 timeOld = 0.0,
                 double                 dT      = 1.0 )
  {
    Material::State state{ Material::Tensor3d( 0.0 ),
                           Material::Tensor33d( 0.0 ),
                           Material::Tensor33d( 0.0 ),
                           0.0,
                           stateVars.data() };

    Material::Tangents tangents{ Material::Tensor33d( 0.0 ),
                                 Material::Tensor333d( 0.0 ),
                                 Material::Tensor333d( 0.0 ),
                                 Material::Tensor333d( 0.0 ),
                                 Material::Tensor3333d( 0.0 ),
                                 Material::Tensor3333d( 0.0 ),
                                 Material::Tensor333d( 0.0 ),
                                 Material::Tensor3333d( 0.0 ),
                                 Material::Tensor3333d( 0.0 ),
                                 Material::Tensor3d( 0.0 ),
                                 Material::Tensor33d( 0.0 ),
                                 Material::Tensor33d( 0.0 ),
                                 Material::Tensor3d( 0.0 ),
                                 Material::Tensor33d( 0.0 ),
                                 Material::Tensor33d( 0.0 ),
                                 0.0 };

    // x = (w, A+, A-): the displacement increments of the two sides differ by w
    double dU[6] = { x( 0 ), x( 1 ), x( 2 ), 0., 0., 0. };

    Material::Deformation deformation{ Material::Tensor6d( dU ),
                                       Material::Tensor18d( x.data() + 3 ),
                                       Material::Tensor3d( normal ),
                                       Material::Tensor3d( separation ),
                                       dPressure };

    material.computeStress( state, tangents, deformation, Material::TimeIncrement{ timeOld, dT } );

    Result result;
    result.response.segment< 3 >( 0 ) = Eigen::Map< const Eigen::Vector3d >( state.generalizedForce.data() );
    result.response.segment< 9 >( 3 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >(
      state.surfaceStressPlus.data() );
    result.response.segment< 9 >( 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >(
      state.surfaceStressMinus.data() );
    result.response( 21 ) = state.volumetricResidual;

    constexpr int RM          = Eigen::RowMajor;
    auto&         T           = result.tangent;
    T.block< 3, 3 >( 0, 0 )   = Eigen::Map< const Eigen::Matrix< double, 3, 3, RM > >( tangents.Q_ww.data() );
    T.block< 3, 9 >( 0, 3 )   = Eigen::Map< const Eigen::Matrix< double, 3, 9, RM > >( tangents.Q_wAp.data() );
    T.block< 3, 9 >( 0, 12 )  = Eigen::Map< const Eigen::Matrix< double, 3, 9, RM > >( tangents.Q_wAm.data() );
    T.block< 9, 3 >( 3, 0 )   = Eigen::Map< const Eigen::Matrix< double, 9, 3, RM > >( tangents.Q_Apw.data() );
    T.block< 9, 9 >( 3, 3 )   = Eigen::Map< const Eigen::Matrix< double, 9, 9, RM > >( tangents.Q_ApAp.data() );
    T.block< 9, 9 >( 3, 12 )  = Eigen::Map< const Eigen::Matrix< double, 9, 9, RM > >( tangents.Q_ApAm.data() );
    T.block< 9, 3 >( 12, 0 )  = Eigen::Map< const Eigen::Matrix< double, 9, 3, RM > >( tangents.Q_Amw.data() );
    T.block< 9, 9 >( 12, 3 )  = Eigen::Map< const Eigen::Matrix< double, 9, 9, RM > >( tangents.Q_AmAp.data() );
    T.block< 9, 9 >( 12, 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 9, RM > >( tangents.Q_AmAm.data() );
    T.block< 3, 1 >( 0, 21 )  = Eigen::Map< const Eigen::Vector3d >( tangents.Q_wp.data() );
    T.block< 9, 1 >( 3, 21 )  = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( tangents.Q_App.data() );
    T.block< 9, 1 >( 12, 21 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( tangents.Q_Amp.data() );
    T.block< 1, 3 >( 21, 0 )  = Eigen::Map< const Eigen::Matrix< double, 1, 3 > >( tangents.Q_pw.data() );
    T.block< 1, 9 >( 21, 3 )  = Eigen::Map< const Eigen::Matrix< double, 1, 9 > >( tangents.Q_pAp.data() );
    T.block< 1, 9 >( 21, 12 ) = Eigen::Map< const Eigen::Matrix< double, 1, 9 > >( tangents.Q_pAm.data() );
    T( 21, 21 )               = tangents.Q_pp;
    return result;
  }

  Vector21d makeStrain( double scale, int pattern )
  {
    Vector21d x;
    for ( int i = 0; i < 21; ++i )
      x( i ) = scale * 1.0e-3 * ( ( ( i * 7 + pattern * 5 ) % 13 ) - 6 ) / 6.0;
    return x;
  }

  std::vector< double > initialState( Material& material )
  {
    std::vector< double > stateVars( material.getNumberOfRequiredStateVars(), 0.0 );
    material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );
    return stateVars;
  }

  /// The tangent, including the columns of the pressure, is the derivative of the response.
  void assertTangentMatchesFiniteDifference( const std::string&           materialName,
                                             const std::vector< double >& properties,
                                             const double*                normal,
                                             const double*                separation,
                                             double                       tolerance,
                                             const std::string&           context )
  {
    Material material( materialName, properties.data(), static_cast< int >( properties.size() ), 1 );
    material.setCharacteristicElementLength( 1.0 );

    // a plastic history first, so that the two faces differ
    std::vector< double > history = initialState( material );
    update( material, history, makeStrain( 1.0, 2 ), 1.0e-1, normal, separation );

    const Vector21d x0 = makeStrain( 0.25, 3 );
    const double    p0 = 2.0e-2;

    const auto evaluateAt = [&]( const Vector22d& xp ) {
      std::vector< double > state = history;
      return update( material, state, xp.head< 21 >(), xp( 21 ), normal, separation );
    };

    Vector22d xp0;
    xp0 << x0, p0;
    const Result base = evaluateAt( xp0 );
    throwExceptionOnFailure( base.response.allFinite() && base.tangent.allFinite(), context + ": nan or inf." );

    Matrix22d finiteDifference;
    for ( int i = 0; i < 22; ++i ) {
      const double perturbation = 1e-7;
      Vector22d    plus = xp0, minus = xp0;
      plus( i ) += perturbation;
      minus( i ) -= perturbation;
      finiteDifference.col( i ) = ( evaluateAt( plus ).response - evaluateAt( minus ).response ) /
                                  ( 2.0 * perturbation );
    }

    const double error = ( finiteDifference - base.tangent ).norm();
    const double scale = std::max( 1.0, base.tangent.norm() );
    throwExceptionOnFailure( error / scale < tolerance,
                             context + ": tangent differs from the finite difference, relative error " +
                               std::to_string( error / scale ) );
  }

  void testTangentIsTheDerivative()
  {
    assertTangentMatchesFiniteDifference( "VONMISES",
                                          plasticProperties,
                                          normalZ,
                                          noSeparation,
                                          1e-5,
                                          "plastic, coincident" );
    assertTangentMatchesFiniteDifference( "VONMISES",
                                          plasticProperties,
                                          normalTilt,
                                          tiltSeparation,
                                          1e-5,
                                          "plastic, tilted" );
    assertTangentMatchesFiniteDifference( "LINEARELASTIC",
                                          elasticProperties,
                                          normalTilt,
                                          tiltSeparation,
                                          1e-6,
                                          "elastic, tilted" );
  }

  /// The pressure column is purely geometric and the volumetric residual does not depend on the local solve:
  /// with the pressure of the constitutive law the residual vanishes, and the response is the elastic one.
  void testConstitutivePressureGivesAnElasticResponse()
  {
    Material   mixed( "LINEARELASTIC", elasticProperties.data(), 3, 1 );
    EqMaterial plain( "LINEARELASTIC", elasticProperties.data(), 3, 1 );
    mixed.setCharacteristicElementLength( 1.0 );
    plain.setCharacteristicElementLength( 1.0 );

    // equal surface gradients of the two faces and a jump along the normal: a homogeneous strain, z = 0
    Vector21d x = Vector21d::Zero();
    x( 2 )      = 3.0e-4; // normal jump
    for ( int face = 0; face < 2; ++face ) {
      x( 3 + 9 * face + 0 ) = 2.0e-3;  // A_xx
      x( 3 + 9 * face + 4 ) = -1.0e-3; // A_yy
    }

    const double E = elasticProperties[0], nu = elasticProperties[1], h = elasticProperties[2];
    const double K                = E / ( 3.0 * ( 1.0 - 2.0 * nu ) );
    const double volumetricStrain = 2.0e-3 - 1.0e-3 + 3.0e-4 / h;

    std::vector< double > mixedState  = initialState( mixed );
    const Result          mixedResult = update( mixed, mixedState, x, -K * volumetricStrain, normalZ, noSeparation );

    std::vector< double > plainState( plain.getNumberOfRequiredStateVars(), 0.0 );
    plain.initializeYourself( plainState.data(), static_cast< int >( plainState.size() ) );
    EqMaterial::State       state{ EqMaterial::Tensor3d( 0.0 ),
                             EqMaterial::Tensor33d( 0.0 ),
                             EqMaterial::Tensor33d( 0.0 ),
                             plainState.data(),
                             nullptr };
    EqMaterial::Tangents    tangents{ EqMaterial::Tensor33d( 0.0 ),
                                   EqMaterial::Tensor333d( 0.0 ),
                                   EqMaterial::Tensor333d( 0.0 ),
                                   EqMaterial::Tensor333d( 0.0 ),
                                   EqMaterial::Tensor3333d( 0.0 ),
                                   EqMaterial::Tensor3333d( 0.0 ),
                                   EqMaterial::Tensor333d( 0.0 ),
                                   EqMaterial::Tensor3333d( 0.0 ),
                                   EqMaterial::Tensor3333d( 0.0 ) };
    double                  dU[6] = { 0., 0., x( 2 ), 0., 0., 0. };
    EqMaterial::Deformation deformation{ EqMaterial::Tensor6d( dU ),
                                         EqMaterial::Tensor18d( x.data() + 3 ),
                                         EqMaterial::Tensor3d( normalZ ),
                                         EqMaterial::Tensor3d( noSeparation ) };
    plain.computeStress( state, tangents, deformation, EqMaterial::TimeIncrement{ 0.0, 1.0 } );

    const Eigen::Map< const Eigen::Vector3d > plainForce( state.generalizedForce.data() );
    const double                              scale = std::max( 1.0, plainForce.norm() );
    throwExceptionOnFailure( ( mixedResult.response.head< 3 >() - plainForce ).norm() / scale < 1e-9,
                             "The force of the mixed material differs from the equilibrated one." );
    throwExceptionOnFailure( std::abs( mixedResult.response( 21 ) ) < 1e-12,
                             "The volumetric residual of the constitutive pressure has to vanish." );
  }

  void testRejectsIncompressibleAndNegativeThickness()
  {
    const std::vector< double > incompressible = { 2e5, 0.5, 0.1 };
    const std::vector< double > noThickness    = { 2e5, 0.3, 0.0 };

    bool thrown = false;
    try {
      Material material( "LINEARELASTIC", incompressible.data(), 3, 1 );
    }
    catch ( const std::invalid_argument& ) {
      thrown = true;
    }
    throwExceptionOnFailure( thrown, "nu = 0.5 has no finite bulk modulus and must be rejected." );

    thrown = false;
    try {
      Material material( "LINEARELASTIC", noThickness.data(), 3, 1 );
    }
    catch ( const std::invalid_argument& ) {
      thrown = true;
    }
    throwExceptionOnFailure( thrown, "h = 0 must be rejected." );
  }

  /// A zero increment after a plastic one leaves the stress where it is. The jump z enters the face strains as an
  /// increment, so the local solve must not start from the z of the previous increment: applied to a zero increment it
  /// would move the faces, and at a load reversal it has the wrong sign.
  void testZeroIncrementAfterPlasticIncrementIsInert()
  {
    const std::vector< double > softPlastic = { 210000., 0.3, 0.2, 120., 0.1, 0.0, 0.0, 0 };
    Material                    material( "VONMISES", softPlastic.data(), static_cast< int >( softPlastic.size() ), 1 );
    material.setCharacteristicElementLength( 1.0 );

    std::vector< double > stateVars = initialState( material );
    const Result          loaded    = update( material, stateVars, makeStrain( 5.0, 2 ), 1.0e-1, normalTilt, tiltSeparation );
    const Result          inert     = update(
      material, stateVars, Vector21d::Zero(), 0.0, normalTilt, tiltSeparation, 1.0, 0.0 );

    throwExceptionOnFailure( loaded.response.allFinite() && inert.response.allFinite(), "nan or inf." );
    const double scale = std::max( 1.0, loaded.response.head< 21 >().norm() );
    throwExceptionOnFailure( ( inert.response.head< 21 >() - loaded.response.head< 21 >() ).norm() / scale < 1e-8,
                             "A zero increment changed the generalized stress." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testTangentIsTheDerivative,
                                                   testConstitutivePressureGivesAnElasticResponse,
                                                   testRejectsIncompressibleAndNegativeThickness,
                                                   testZeroIncrementAfterPlasticIncrementIsInert };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
