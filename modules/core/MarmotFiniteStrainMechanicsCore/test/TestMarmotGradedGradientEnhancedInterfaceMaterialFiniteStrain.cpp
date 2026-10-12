#include "Marmot/MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Material = MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain;

  constexpr double h = 0.1;

  struct Geometry {
    Eigen::Vector3d normal;
    Eigen::Vector3d connector;
  };

  Geometry inclinedGeometry( double alpha )
  {
    return { Eigen::Vector3d( -std::sin( alpha ), 0.0, std::cos( alpha ) ), Eigen::Vector3d( 0.0, 0.0, h ) };
  }

  /// Kinematics: the displacements of both sides, their (equal) surface gradients and the nonlocal fields.
  struct Kinematics {
    std::array< double, 6 >  U{};
    std::array< double, 18 > surfaceGradient{};
    Eigen::VectorXd          nonlocal;
  };

  Kinematics arbitraryKinematics( double scale, int nStations )
  {
    Kinematics kinematics;
    for ( int i = 0; i < 3; i++ ) {
      kinematics.U[i]     = scale * h * std::sin( 0.7 * i + 0.2 );
      kinematics.U[3 + i] = -0.3 * scale * h * std::cos( 0.5 * i + 0.1 );
    }
    for ( int i = 0; i < 9; i++ ) {
      const double entry                = 0.6 * scale * std::sin( 1.1 * i + 0.4 );
      kinematics.surfaceGradient[i]     = entry;
      kinematics.surfaceGradient[9 + i] = entry;
    }
    kinematics.nonlocal.resize( nStations );
    for ( int q = 0; q < nStations; q++ )
      kinematics.nonlocal( q ) = 0.002 * ( 1.0 + 0.7 * q );
    return kinematics;
  }

  struct Response {
    Eigen::Vector3d                                force;
    Eigen::Matrix< double, 9, 1 >                  surfaceStress;
    Eigen::VectorXd                                L;
    Eigen::VectorXd                                radius;
    Eigen::Matrix3d                                Q;
    Eigen::Matrix< double, 3, 9, Eigen::RowMajor > H_ikR;
    Eigen::Matrix< double, 9, 3, Eigen::RowMajor > H_iRk;
    Eigen::Matrix< double, 9, 9, Eigen::RowMajor > A;
    Eigen::MatrixXd                                dForce_dN, dSurface_dN, dL_dJump, dL_dAverage, dL_dN;
  };

  std::vector< double > gradedCard( const std::vector< double >& positions, const std::vector< double >& bulk )
  {
    std::vector< double > properties = { h, 2.0, 1.0, double( positions.size() ), double( bulk.size() ) };
    for ( const double position : positions ) {
      properties.push_back( position );
      properties.insert( properties.end(), bulk.begin(), bulk.end() );
    }
    return properties;
  }

  Response evaluate( const std::string&           name,
                     const std::vector< double >& card,
                     const Geometry&              geometry,
                     const Kinematics&            kinematics,
                     std::vector< double >*       stateVarsOut = nullptr )
  {
    Material              material( name, card.data(), card.size(), 1 );
    std::vector< double > stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    Material::State       state{ Material::Tensor3d( 0.0 ),
                           Material::Tensor33d( 0.0 ),
                           stateVars.data(),
                           Eigen::VectorXd(),
                           Eigen::VectorXd() };
    Material::Tangents    tangents{ Material::Tensor33d( 0.0 ),
                                 Material::Tensor333d( 0.0 ),
                                 Material::Tensor333d( 0.0 ),
                                 Material::Tensor3333d( 0.0 ),
                                    {},
                                    {},
                                    {},
                                    {},
                                    {} };
    Material::Deformation deformation{ Material::Tensor6d( kinematics.U.data() ),
                                       Material::Tensor18d( kinematics.surfaceGradient.data() ),
                                       Material::Tensor3d( geometry.normal.data() ),
                                       Material::Tensor3d( geometry.connector.data() ),
                                       kinematics.nonlocal };

    material.computeStress( state, tangents, deformation, { 0.0, 1.0 } );

    Response response;
    response.force         = Eigen::Map< const Eigen::Vector3d >( state.force.data() );
    response.surfaceStress = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( state.surfaceStress.data() );
    response.L             = state.localFields;
    response.radius        = state.nonlocalRadii;
    response.Q             = Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( tangents.Q_ik.data() );
    response.H_ikR       = Eigen::Map< const Eigen::Matrix< double, 3, 9, Eigen::RowMajor > >( tangents.H_ikR.data() );
    response.H_iRk       = Eigen::Map< const Eigen::Matrix< double, 9, 3, Eigen::RowMajor > >( tangents.H_iRk.data() );
    response.A           = Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >( tangents.A_iRkP.data() );
    response.dForce_dN   = tangents.dForce_dNonlocal;
    response.dSurface_dN = tangents.dSurfaceStress_dNonlocal;
    response.dL_dJump    = tangents.dLocal_dJump;
    response.dL_dAverage = tangents.dLocal_dAverage;
    response.dL_dN       = tangents.dLocal_dNonlocal;

    if ( stateVarsOut )
      *stateVarsOut = stateVars;
    return response;
  }

  /// All tangents are the derivatives of the response, including the ones of the coupling to the nonlocal fields.
  void checkTangents( const std::string&           name,
                      const std::vector< double >& bulk,
                      const std::vector< double >& positions,
                      double                       scale,
                      double                       alpha,
                      const std::string&           context )
  {
    const auto card       = gradedCard( positions, bulk );
    const int  nStations  = static_cast< int >( positions.size() );
    const auto geometry   = inclinedGeometry( alpha );
    const auto kinematics = arbitraryKinematics( scale, nStations );

    const Response base = evaluate( name, card, geometry, kinematics );

    throwExceptionOnFailure( base.L.maxCoeff() > 0.0,
                             context + ": the test does not drive the local field, so it does not test its tangents." );

    const double eps = 1e-7;

    const auto close =
      [&]( const Eigen::VectorXd& numerical, const Eigen::VectorXd& analytical, const std::string& what ) {
        const double tolerance = 2e-4 * std::max( analytical.norm(), 1e-4 * std::max( 1.0, numerical.norm() ) ) + 1e-8;
        throwExceptionOnFailure( ( numerical - analytical ).norm() <= tolerance,
                                 context + ": " + what + " does not match the finite difference, error " +
                                   std::to_string( ( numerical - analytical ).norm() ) + " for a norm of " +
                                   std::to_string( analytical.norm() ) + "." );
      };

    // jump
    for ( int k = 0; k < 3; k++ ) {
      Kinematics perturbed = kinematics;
      perturbed.U[k] += eps;
      const Response r = evaluate( name, card, geometry, perturbed );
      close( ( r.force - base.force ) / eps, base.Q.col( k ), "dForce/dJump" );
      close( ( r.surfaceStress - base.surfaceStress ) / eps, base.H_iRk.col( k ), "dSurfaceStress/dJump" );
      close( ( r.L - base.L ) / eps, base.dL_dJump.col( k ), "dLocal/dJump" );
    }

    // average surface gradient, along the tangential part e_k (x) T e_R
    const Eigen::Matrix3d T = Eigen::Matrix3d::Identity() - geometry.normal * geometry.normal.transpose();
    for ( int k = 0; k < 3; k++ )
      for ( int R = 0; R < 3; R++ ) {
        Kinematics perturbed = kinematics;
        for ( int P = 0; P < 3; P++ ) {
          perturbed.surfaceGradient[3 * k + P] += eps * T( P, R );
          perturbed.surfaceGradient[9 + 3 * k + P] += eps * T( P, R );
        }
        const Response r = evaluate( name, card, geometry, perturbed );
        close( ( r.force - base.force ) / eps, base.H_ikR.col( 3 * k + R ), "dForce/dSurfaceGradient" );
        close( ( r.surfaceStress - base.surfaceStress ) / eps,
               base.A.col( 3 * k + R ),
               "dSurfaceStress/dSurfaceGradient" );
        close( ( r.L - base.L ) / eps, base.dL_dAverage.col( 3 * k + R ), "dLocal/dSurfaceGradient" );
      }

    // nonlocal field of every station
    for ( int f = 0; f < nStations; f++ ) {
      Kinematics perturbed = kinematics;
      perturbed.nonlocal( f ) += eps;
      const Response r = evaluate( name, card, geometry, perturbed );
      close( ( r.force - base.force ) / eps, base.dForce_dN.col( f ), "dForce/dNonlocal" );
      close( ( r.surfaceStress - base.surfaceStress ) / eps, base.dSurface_dN.col( f ), "dSurfaceStress/dNonlocal" );
      close( ( r.L - base.L ) / eps, base.dL_dN.col( f ), "dLocal/dNonlocal" );
    }
  }

  const std::vector< double > neoHookeDamage = { 3500.0, 1500.0, 1.0e-3, 5.0e-3, 0.2, 1.0 };
  const std::string           neoHookeName   = "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE";

  const std::vector< double > dufour = { 1793.9, 828.0, 12.0,   21.6, 4051.0, 15.8, 552.0, 243.0, 0.0, 0.0,
                                         15.4,   1e-6,  0.0435, 0.3,  0.5,    1.75, 0.99,  0.12,  1.0, 1.0 };

  void testTangentsOfTheNeoHookeDamageModel()
  {
    checkTangents( neoHookeName, neoHookeDamage, { 0.0, 0.5, 1.0 }, 0.15, 0.0, "damage, three stations" );
    checkTangents( neoHookeName, neoHookeDamage, { 0.0, 0.3, 1.0 }, 0.15, 0.4, "damage, inclined, uneven stations" );
    checkTangents( neoHookeName, neoHookeDamage, { 0.0, 1.0 }, 0.15, 0.4, "damage, two stations" );
  }

  void testTangentsOfTheDufourModel()
  {
    checkTangents( "DUFOURMODEL", dufour, { 0.0, 0.5, 1.0 }, 0.15, 0.0, "Dufour, three stations" );
    checkTangents( "DUFOURMODEL", dufour, { 0.0, 1.0 }, 0.15, 0.4, "Dufour, two stations, inclined" );
  }

  /// With equal stations and fields, the layer is homogeneous: the response does not depend on the number of stations.
  void testHomogeneousLayerDoesNotDependOnTheStations()
  {
    const auto geometry = inclinedGeometry( 0.3 );

    const auto single = [&]( int n ) {
      std::vector< double > positions;
      for ( int q = 0; q < n; q++ )
        positions.push_back( double( q ) / ( n - 1 ) );
      auto kinematics     = arbitraryKinematics( 0.1, n );
      kinematics.nonlocal = Eigen::VectorXd::Constant( n, 0.003 );
      return evaluate( neoHookeName, gradedCard( positions, neoHookeDamage ), geometry, kinematics );
    };

    const Response two   = single( 2 );
    const Response three = single( 3 );
    const Response five  = single( 5 );

    throwExceptionOnFailure( ( two.force - three.force ).norm() < 1e-9 * std::max( 1.0, two.force.norm() ) &&
                               ( two.force - five.force ).norm() < 1e-9 * std::max( 1.0, two.force.norm() ),
                             "The force of a homogeneous layer depends on the number of stations." );
    throwExceptionOnFailure( ( two.surfaceStress - five.surfaceStress ).norm() <
                               1e-9 * std::max( 1.0, two.surfaceStress.norm() ),
                             "The surface stress of a homogeneous layer depends on the number of stations." );
    throwExceptionOnFailure( std::abs( two.L( 0 ) - five.L( 0 ) ) < 1e-9 * std::max( 1.0, std::abs( two.L( 0 ) ) ),
                             "The local field of a homogeneous layer depends on the number of stations." );
  }

  template < class Callable >
  bool throwsInvalidArgument( Callable&& callable )
  {
    try {
      callable();
    }
    catch ( const std::invalid_argument& ) {
      return true;
    }
    return false;
  }

  void testInvalidInput()
  {
    const auto valid = gradedCard( { 0.0, 0.5, 1.0 }, neoHookeDamage );

    throwExceptionOnFailure( !throwsInvalidArgument(
                               [&] { Material m( neoHookeName, valid.data(), valid.size(), 1 ); } ),
                             "A valid card was rejected." );
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { Material m( "COMPRESSIBLENEOHOOKE", valid.data(), valid.size(), 1 ); } ),
                             "A bulk material without a nonlocal field was accepted." );

    auto layers = valid;
    layers[1]   = 0.0;
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { Material m( neoHookeName, layers.data(), layers.size(), 1 ); } ),
                             "A card that is not a profile of stations was accepted." );

    Material              material( neoHookeName, valid.data(), valid.size(), 1 );
    std::vector< double > stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    throwExceptionOnFailure( material.numberOfStations() == 3, "Three stations were expected." );
    throwExceptionOnFailure( std::abs( material.stationWeight( 0 ) - 0.25 ) < 1e-14 &&
                               std::abs( material.stationWeight( 1 ) - 0.5 ) < 1e-14 &&
                               std::abs( material.stationWeight( 2 ) - 0.25 ) < 1e-14,
                             "The weights of the stations are not the ones of the trapezoidal rule." );
    throwExceptionOnFailure( material.stationPosition( 1 ) == 0.5 && material.thickness() == h,
                             "The geometry of the stations is wrong." );
    throwExceptionOnFailure( std::isfinite( material.getDensity( stateVars.data() ) ), "The density is not finite." );
    material.setCharacteristicElementLength( 0.1 );

    // a wrong number of fields
    Kinematics kinematics = arbitraryKinematics( 0.1, 2 );
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { evaluate( neoHookeName, valid, inclinedGeometry( 0. ), kinematics ); } ),
                             "A wrong number of nonlocal fields was accepted." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testTangentsOfTheNeoHookeDamageModel,
                                                   testTangentsOfTheDufourModel,
                                                   testHomogeneousLayerDoesNotDependOnTheStations,
                                                   testInvalidInput };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
