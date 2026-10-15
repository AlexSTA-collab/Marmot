#include "Marmot/MarmotEquilibratedInterfaceMaterialFiniteStrain.h"
#include "Marmot/MarmotGradedEquilibratedInterfaceMaterialFiniteStrain.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Graded       = MarmotGradedEquilibratedInterfaceMaterialFiniteStrain;
  using Equilibrated = MarmotEquilibratedInterfaceMaterialFiniteStrain;

  constexpr double h = 0.1;

  const std::vector< double > neoHooke = { 175000., 80769.2, 1.0 };

  /// Generalized stress p = ( force, S+, S- ) and the tangent over the generalized strain ( w, A+, A- ).
  struct Response {
    Eigen::Matrix< double, 21, 1 >  p         = Eigen::Matrix< double, 21, 1 >::Zero();
    Eigen::Matrix< double, 21, 1 >  resultant = Eigen::Matrix< double, 21, 1 >::Zero(); ///< moments of the stress
    Eigen::Matrix< double, 21, 21 > K         = Eigen::Matrix< double, 21, 21 >::Zero();
  };

  /// Total kinematics of an interface point.
  struct Kinematics {
    Eigen::Matrix< double, 6, 1 >  U         = Eigen::Matrix< double, 6, 1 >::Zero();
    Eigen::Matrix< double, 18, 1 > gradient  = Eigen::Matrix< double, 18, 1 >::Zero();
    Eigen::Vector3d                normal    = Eigen::Vector3d( 0., 0., 1. );
    Eigen::Vector3d                connector = Eigen::Vector3d( 0., 0., h );
  };

  /// Material properties of the graded material: [h, mode, nGauss, K, m, (position, properties) x K, flag].
  std::vector< double > card( int                                                              mode,
                              int                                                              nGauss,
                              const std::vector< std::pair< double, std::vector< double > > >& groups,
                              int                                                              tangentialGradient = -1 )
  {
    std::vector< double > properties = { h,
                                         double( mode ),
                                         double( nGauss ),
                                         double( groups.size() ),
                                         double( groups[0].second.size() ) };
    for ( const auto& [position, bulk] : groups ) {
      properties.push_back( position );
      properties.insert( properties.end(), bulk.begin(), bulk.end() );
    }
    if ( tangentialGradient >= 0 )
      properties.push_back( double( tangentialGradient ) );
    return properties;
  }

  template < class InterfaceMaterial >
  Response evaluate( InterfaceMaterial& material, double* stateVars, const Kinematics& kinematics, double timeOld = 0. )
  {
    typename InterfaceMaterial::State state{ typename InterfaceMaterial::Tensor3d( 0.0 ),
                                             typename InterfaceMaterial::Tensor33d( 0.0 ),
                                             typename InterfaceMaterial::Tensor33d( 0.0 ),
                                             stateVars };

    typename InterfaceMaterial::Tangents tangents{ typename InterfaceMaterial::Tensor33d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor333d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor333d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor333d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor3333d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor3333d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor333d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor3333d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor3333d( 0.0 ) };

    typename InterfaceMaterial::Deformation
                                              deformation{ typename InterfaceMaterial::Tensor6d( kinematics.U.data() ),
                   typename InterfaceMaterial::Tensor18d( kinematics.gradient.data() ),
                   typename InterfaceMaterial::Tensor3d( kinematics.normal.data() ),
                   typename InterfaceMaterial::Tensor3d( kinematics.connector.data() ) };
    typename InterfaceMaterial::TimeIncrement timeIncrement{ timeOld, 1.0 };

    material.computeStress( state, tangents, deformation, timeIncrement );

    Response response;
    if constexpr ( requires { material.stressResultant(); } )
      response.resultant = material.stressResultant();
    response.p.segment< 3 >( 0 )  = Eigen::Map< const Eigen::Vector3d >( state.generalizedForce.data() );
    response.p.segment< 9 >( 3 )  = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( state.surfaceStressPlus.data() );
    response.p.segment< 9 >( 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >(
      state.surfaceStressMinus.data() );

    using RowMajor   = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >;
    const auto block = [&]( const double* data, int rows, int cols ) {
      return Eigen::Map< const RowMajor >( data, rows, cols );
    };
    response.K.block< 3, 3 >( 0, 0 )   = block( tangents.Q_ww.data(), 3, 3 );
    response.K.block< 3, 9 >( 0, 3 )   = block( tangents.Q_wAp.data(), 3, 9 );
    response.K.block< 3, 9 >( 0, 12 )  = block( tangents.Q_wAm.data(), 3, 9 );
    response.K.block< 9, 3 >( 3, 0 )   = block( tangents.Q_Apw.data(), 9, 3 );
    response.K.block< 9, 9 >( 3, 3 )   = block( tangents.Q_ApAp.data(), 9, 9 );
    response.K.block< 9, 9 >( 3, 12 )  = block( tangents.Q_ApAm.data(), 9, 9 );
    response.K.block< 9, 3 >( 12, 0 )  = block( tangents.Q_Amw.data(), 9, 3 );
    response.K.block< 9, 9 >( 12, 3 )  = block( tangents.Q_AmAp.data(), 9, 9 );
    response.K.block< 9, 9 >( 12, 12 ) = block( tangents.Q_AmAm.data(), 9, 9 );
    return response;
  }

  struct GradedPoint {
    std::unique_ptr< Graded > material;
    std::vector< double >     stateVars;

    explicit GradedPoint( const std::vector< double >& properties )
    {
      material = std::make_unique< Graded >( "COMPRESSIBLENEOHOOKE",
                                             properties.data(),
                                             static_cast< int >( properties.size() ),
                                             1 );
      stateVars.assign( material->getNumberOfRequiredStateVars(), 0.0 );
      material->initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );
    }
  };

  /// A total kinematics with every component nonzero, whose surface gradients are tangential.
  Kinematics makeKinematics( double scale, bool inclined )
  {
    Kinematics kinematics;
    if ( inclined ) {
      kinematics.normal    = Eigen::Vector3d( 1., 2., 2. ) / 3.;
      kinematics.connector = 1.05 * h * kinematics.normal + 0.1 * h * Eigen::Vector3d( 2., -2., 1. ) / 3.;
    }
    for ( int i = 0; i < 3; i++ ) {
      kinematics.U( i )     = scale * h * std::sin( 0.7 * i + 0.2 );
      kinematics.U( 3 + i ) = -0.3 * scale * h * std::cos( 0.5 * i + 0.1 );
    }
    for ( int i = 0; i < 18; i++ )
      kinematics.gradient( i ) = scale * std::sin( 1.1 * i + 0.4 );
    for ( int block = 0; block < 2; block++ ) {
      Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > A( kinematics.gradient.data() + 9 * block );
      A = A - ( A * kinematics.normal ) * kinematics.normal.transpose();
    }
    return kinematics;
  }

  /// With two points at the faces and without the tangential gradient, the material is the equilibrated one.
  void testTwoFacesAreEquilibratedMaterial()
  {
    const std::vector< double > equilibratedCard = { neoHooke[0], neoHooke[1], h, neoHooke[2] };

    for ( const bool inclined : { false, true } )
      for ( const double scale : { 1e-3, 0.05, 0.2 } ) {
        Equilibrated          equilibrated( "COMPRESSIBLENEOHOOKE",
                                   equilibratedCard.data(),
                                   static_cast< int >( equilibratedCard.size() ),
                                   1 );
        std::vector< double > equilibratedState( equilibrated.getNumberOfRequiredStateVars(), 0.0 );
        equilibrated.initializeYourself( equilibratedState.data(), static_cast< int >( equilibratedState.size() ) );

        GradedPoint graded( card( 2, 1, { { 0.0, neoHooke }, { 1.0, neoHooke } }, 0 ) );

        const Kinematics kinematics = makeKinematics( scale, inclined );
        const Response   expected   = evaluate( equilibrated, equilibratedState.data(), kinematics );
        const Response   actual     = evaluate( *graded.material, graded.stateVars.data(), kinematics );

        throwExceptionOnFailure( checkIfEqual< double >( actual.p, expected.p, 1e-8 * ( 1. + expected.p.norm() ) ),
                                 "The generalized stress is not the one of the equilibrated material at the scale " +
                                   std::to_string( scale ) );
        throwExceptionOnFailure( checkIfEqual< double >( Eigen::VectorXd(
                                                           Eigen::Map< const Eigen::VectorXd >( actual.K.data(),
                                                                                                441 ) ),
                                                         Eigen::VectorXd(
                                                           Eigen::Map< const Eigen::VectorXd >( expected.K.data(),
                                                                                                441 ) ),
                                                         1e-8 * ( 1. + expected.K.norm() ) ),
                                 "The tangent is not the one of the equilibrated material at the scale " +
                                   std::to_string( scale ) );
      }
  }

  /**
   * The tangent is the derivative of the generalized stress through a layered stack. Without the closure of the
   * tangential gradient it is exact. With it, the tangent neglects the dependence of the derivative of the
   * correction on the deformation, which is of the order of the strain times the closure terms.
   */
  void testTangentIsConsistent()
  {
    const std::vector< double > stiff = { 400000., 200000., 1.0 };

    for ( const bool closure : { false, true } ) {
      const auto layers = card( 0, 2, { { 0.3, neoHooke }, { 0.7, stiff } }, closure ? 1 : 0 );

      for ( const bool inclined : { false, true } )
        for ( const double scale : { 1e-3, 0.05 } ) {
          const Kinematics kinematics = makeKinematics( scale, inclined );

          GradedPoint    reference( layers );
          const Response response = evaluate( *reference.material, reference.stateVars.data(), kinematics );

          Eigen::Matrix< double, 21, 21 > numerical;
          for ( int j = 0; j < 21; j++ ) {
            const double step = 1e-7 * ( j < 3 ? h : 1.0 );
            Response     sides[2];
            for ( int side = 0; side < 2; side++ ) {
              GradedPoint  point( layers );
              Kinematics   perturbed = kinematics;
              const double delta     = side == 0 ? -step : step;
              if ( j < 3 )
                perturbed.U( j ) += delta;
              else
                perturbed.gradient( j - 3 ) += delta;
              sides[side] = evaluate( *point.material, point.stateVars.data(), perturbed );
            }
            numerical.col( j ) = ( sides[1].p - sides[0].p ) / ( 2. * step );
          }

          const double tolerance = ( closure ? scale : 1e-6 ) * response.K.cwiseAbs().maxCoeff();
          for ( int i = 0; i < 21; i++ )
            for ( int j = 0; j < 21; j++ )
              throwExceptionOnFailure( std::abs( numerical( i, j ) - response.K( i, j ) ) < tolerance,
                                       "The tangent of the finite-strain stack is not the derivative of the stress "
                                       "at (" +
                                         std::to_string( i ) + "," + std::to_string( j ) + ")" );
        }
    }
  }

  /// A rigid rotation of the whole layer is stress free, also for a stack with a stiffness contrast.
  void testRigidRotationIsStressFree()
  {
    const std::vector< double > stiff  = { 400000., 200000., 1.0 };
    const auto                  layers = card( 0, 2, { { 0.3, neoHooke }, { 0.7, stiff } } );

    const double          angle = 0.4;
    const Eigen::Vector3d axis  = Eigen::Vector3d( 1., 2., 2. ).normalized();
    const Eigen::Matrix3d R     = Eigen::AngleAxisd( angle, axis ).toRotationMatrix();

    Kinematics kinematics;
    kinematics.normal    = Eigen::Vector3d( 0., 0., 1. );
    kinematics.connector = Eigen::Vector3d( 0., 0., h );

    // u = ( R - I ) X, the bottom face at X = 0 and the top face at the connector
    const Eigen::Vector3d topDisplacement = ( R - Eigen::Matrix3d::Identity() ) * kinematics.connector;
    kinematics.U.segment< 3 >( 0 )        = topDisplacement;

    // the surface gradient is ( R - I ) restricted to the tangent plane
    Eigen::Matrix3d gradient = ( R - Eigen::Matrix3d::Identity() );
    gradient.col( 2 ).setZero();
    for ( int block = 0; block < 2; block++ )
      Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( kinematics.gradient.data() + 9 * block ) = gradient;

    GradedPoint    point( layers );
    const Response response = evaluate( *point.material, point.stateVars.data(), kinematics );

    throwExceptionOnFailure( response.p.cwiseAbs().maxCoeff() < 1e-8,
                             "A rigid rotation of a graded layer is not stress free." );
  }

  /**
   * At small strains the layered stack is the solution of plane strain elasticity, as in the test of the
   * small-strain material: the mean of the stress and the moments h int y P dy and h int (1-y) P dy at x = 0 of a
   * stack of the lower layer K = 175000, G = 80769.2 (40 % of the thickness) and the upper layer K = 400000,
   * G = 200000, with the jump z0 = 2e-8 of the normal displacement and the gradient g = 1e-7 of the displacement of
   * the top face along the interface. The reference was derived with the computer algebra system sympy.
   */
  void testLayeredStackAgainstExactPlaneStrainSolution()
  {
    const std::vector< double > stiff  = { 400000., 200000., 1.0 };
    const auto                  layers = card( 0, 3, { { 0.4, neoHooke }, { 0.6, stiff } } );

    GradedPoint point( layers );

    Kinematics kinematics;
    kinematics.U( 2 )                = 2e-8;
    kinematics.gradient( 9 + 0 )     = 0.;   // A-
    kinematics.gradient( 3 * 2 + 0 ) = 1e-7; // A+_zx

    const Response response = evaluate( *point.material, point.stateVars.data(), kinematics );

    const auto S = [&]( int offset, int i, int j ) { return response.resultant( offset + 3 * i + j ); };

    struct Reference {
      const char* name;
      double      value;
      double      actual;
    };
    const Reference references[] = {
      { "force x", 0.00563602843723195, response.resultant( 0 ) },
      { "force z", 0.08639435025596491, response.resultant( 2 ) },
      { "S+ xx", 0.0017476343783169167, S( 3, 0, 0 ) },
      { "S+ yy", 0.0017476343783169167, S( 3, 1, 1 ) },
      { "S+ zz", 0.0043197175127982455, S( 3, 2, 2 ) },
      { "S+ xz", 0.00013419048634875927, S( 3, 0, 2 ) },
      { "S- xx", 0.001806876497909772, S( 12, 0, 0 ) },
      { "S- zz", 0.0043197175127982455, S( 12, 2, 2 ) },
      { "S- xz", 0.00042941235737443576, S( 12, 0, 2 ) },
    };
    for ( const auto& reference : references )
      throwExceptionOnFailure( std::abs( reference.actual - reference.value ) < 1e-5 * std::abs( reference.value ),
                               std::string( "The finite-strain stack differs from the plane strain solution in " ) +
                                 reference.name );
  }

  void testInvalidPropertiesAreRefused()
  {
    bool thrown = false;
    try {
      const auto properties = card( 0, 1, { { 1.0, neoHooke } } );
      Graded     material( "NOSUCHMATERIAL", properties.data(), static_cast< int >( properties.size() ), 1 );
    }
    catch ( const std::invalid_argument& ) {
      thrown = true;
    }
    throwExceptionOnFailure( thrown, "An unknown bulk material is not refused." );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{ testTwoFacesAreEquilibratedMaterial,
                                                       testTangentIsConsistent,
                                                       testRigidRotationIsStressFree,
                                                       testLayeredStackAgainstExactPlaneStrainSolution,
                                                       testInvalidPropertiesAreRefused };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
