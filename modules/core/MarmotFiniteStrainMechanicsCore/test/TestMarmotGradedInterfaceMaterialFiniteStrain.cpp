#include "Marmot/MarmotGradedInterfaceMaterialFiniteStrain.h"
#include "Marmot/MarmotInterfaceMaterialFiniteStrain.h"
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

  constexpr double h = 0.1;

  using RowMajor3 = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;

  /// Geometry of an interface inclined by alpha with its faces separated vertically.
  struct Geometry {
    Eigen::Vector3d normal;
    Eigen::Vector3d connector;
  };

  Geometry inclinedGeometry( double alpha )
  {
    return { Eigen::Vector3d( -std::sin( alpha ), 0.0, std::cos( alpha ) ), Eigen::Vector3d( 0.0, 0.0, h ) };
  }

  /// Kinematics: the displacements of both sides and their (equal) surface gradients.
  struct Kinematics {
    std::array< double, 6 >  U{};
    std::array< double, 18 > surfaceGradient{};
  };

  Kinematics arbitraryKinematics( double scale )
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
    return kinematics;
  }

  struct Response {
    Eigen::Vector3d                                force         = Eigen::Vector3d::Zero();
    RowMajor3                                      surfaceStress = RowMajor3::Zero();
    RowMajor3                                      Q             = RowMajor3::Zero();
    Eigen::Matrix< double, 3, 9, Eigen::RowMajor > H_ikR;
    Eigen::Matrix< double, 9, 3, Eigen::RowMajor > H_iRk;
    Eigen::Matrix< double, 9, 9, Eigen::RowMajor > A;
  };

  template < class InterfaceMaterial >
  Response evaluate( InterfaceMaterial& material,
                     double*            stateVars,
                     const Geometry&    geometry,
                     const Kinematics&  kinematics,
                     double             timeOld = 0.0,
                     double             dT      = 1.0 )
  {
    typename InterfaceMaterial::State         state{ typename InterfaceMaterial::Tensor3d( 0.0 ),
                                             typename InterfaceMaterial::Tensor33d( 0.0 ),
                                             stateVars };
    typename InterfaceMaterial::Tangents      tangents{ typename InterfaceMaterial::Tensor33d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor333d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor333d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor3333d( 0.0 ) };
    typename InterfaceMaterial::Deformation   deformation{ typename InterfaceMaterial::Tensor6d( kinematics.U.data() ),
                                                         typename InterfaceMaterial::Tensor18d(
                                                           kinematics.surfaceGradient.data() ),
                                                         typename InterfaceMaterial::Tensor3d( geometry.normal.data() ),
                                                         typename InterfaceMaterial::Tensor3d(
                                                           geometry.connector.data() ) };
    typename InterfaceMaterial::TimeIncrement timeIncrement{ timeOld, dT };

    material.computeStress( state, tangents, deformation, timeIncrement );

    Response response;
    response.force         = Eigen::Map< const Eigen::Vector3d >( state.force.data() );
    response.surfaceStress = Eigen::Map< const RowMajor3 >( state.surfaceStress.data() );
    response.Q             = Eigen::Map< const RowMajor3 >( tangents.Q_ik.data() );
    response.H_ikR = Eigen::Map< const Eigen::Matrix< double, 3, 9, Eigen::RowMajor > >( tangents.H_ikR.data() );
    response.H_iRk = Eigen::Map< const Eigen::Matrix< double, 9, 3, Eigen::RowMajor > >( tangents.H_iRk.data() );
    response.A     = Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >( tangents.A_iRkP.data() );
    return response;
  }

  /// Material properties of the graded layer: [h, mode, nGauss, K, m, (position, bulk properties) x K].
  std::vector< double > gradedCard( const std::vector< std::pair< double, std::vector< double > > >& layers,
                                    int                                                              mode   = 0,
                                    int                                                              nGauss = 2 )
  {
    std::vector< double > properties = { h,
                                         double( mode ),
                                         double( nGauss ),
                                         double( layers.size() ),
                                         double( layers[0].second.size() ) };
    for ( const auto& [position, bulk] : layers ) {
      properties.push_back( position );
      properties.insert( properties.end(), bulk.begin(), bulk.end() );
    }
    return properties;
  }

  using Graded = MarmotGradedInterfaceMaterialFiniteStrain;
  using Plain  = MarmotInterfaceMaterialFiniteStrain;

  void checkResponsesEqual( const Response& a, const Response& b, const std::string& what, double tolerance = 1e-8 )
  {
    throwExceptionOnFailure( checkIfEqual< double >( a.force, b.force, tolerance ), what + ": force differs." );
    throwExceptionOnFailure( checkIfEqual< double >( a.surfaceStress, b.surfaceStress, tolerance ),
                             what + ": surface stress differs." );
    throwExceptionOnFailure( checkIfEqual< double >( a.Q, b.Q, tolerance ), what + ": Q differs." );
    throwExceptionOnFailure( checkIfEqual< double >( a.H_ikR, b.H_ikR, tolerance ), what + ": H_ikR differs." );
    throwExceptionOnFailure( checkIfEqual< double >( a.H_iRk, b.H_iRk, tolerance ), what + ": H_iRk differs." );
    throwExceptionOnFailure( checkIfEqual< double >( a.A, b.A, tolerance ), what + ": A differs." );
  }

  const std::vector< double > neoHooke = { 175000., 80769.2, 1.0 };

  /// A homogeneous layer, however it is described, is the interface material of the plain element.
  void testHomogeneousLayerEqualsPlainInterfaceMaterial()
  {
    std::vector< double > plainCard = { h };
    plainCard.insert( plainCard.end(), neoHooke.begin(), neoHooke.end() );

    const std::vector< std::vector< double > > cards = {
      gradedCard( { { 1.0, neoHooke } }, 0, 1 ),
      gradedCard( { { 0.3, neoHooke }, { 0.7, neoHooke } }, 0, 3 ),
      gradedCard( { { 0.0, neoHooke }, { 0.4, neoHooke }, { 1.0, neoHooke } }, 1, 2 ),
    };

    for ( const double alpha : { 0.0, 0.35 } ) {
      const Geometry geometry = inclinedGeometry( alpha );

      for ( const double scale : { 1e-3, 0.05, 0.2 } ) {
        const Kinematics kinematics = arbitraryKinematics( scale );

        Plain                 plain( "COMPRESSIBLENEOHOOKE", plainCard.data(), plainCard.size(), 1 );
        std::vector< double > plainState( plain.getNumberOfRequiredStateVars() );
        plain.initializeYourself( plainState.data(), plainState.size() );
        const Response expected = evaluate( plain, plainState.data(), geometry, kinematics );

        for ( const auto& card : cards ) {
          Graded                graded( "COMPRESSIBLENEOHOOKE", card.data(), card.size(), 1 );
          std::vector< double > gradedState( graded.getNumberOfRequiredStateVars() );
          graded.initializeYourself( gradedState.data(), gradedState.size() );
          const Response actual = evaluate( graded, gradedState.data(), geometry, kinematics );
          checkResponsesEqual( actual, expected, "A homogeneous layer at scale " + std::to_string( scale ) );
        }
      }
    }
  }

  /// In the small-strain limit the layers act in series for a normal jump.
  void testLayersActInSeriesForSmallStrains()
  {
    const std::vector< double > soft = { 50000., 20000., 1.0 };
    const std::vector< double > hard = { 400000., 200000., 1.0 };

    const std::array< double, 2 > fractions = { 0.3, 0.7 };
    const auto                    card      = gradedCard( { { fractions[0], soft }, { fractions[1], hard } }, 0, 2 );

    Graded                material( "COMPRESSIBLENEOHOOKE", card.data(), card.size(), 1 );
    std::vector< double > stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    const Geometry geometry = { Eigen::Vector3d( 0., 0., 1. ), Eigen::Vector3d( 0., 0., h ) };

    const double jump = 1e-9;
    Kinematics   kinematics;
    kinematics.U[2] = jump;

    const Response response = evaluate( material, stateVars.data(), geometry, kinematics );

    const auto   constrainedModulus = []( const std::vector< double >& p ) { return p[0] + 4. / 3. * p[1]; };
    const double compliance = fractions[0] / constrainedModulus( soft ) + fractions[1] / constrainedModulus( hard );

    throwExceptionOnFailure( checkIfEqual( response.force( 2 ), jump / h / compliance, 1e-6 ),
                             "The normal traction of the stack is not the one of the series model." );
    throwExceptionOnFailure( checkIfEqual( response.Q( 2, 2 ), 1. / ( h * compliance ), 1e-6 ),
                             "The normal stiffness of the stack is not the one of the series model." );
  }

  /// The tangents are the derivatives of the response at large deformation, through a stiffness contrast.
  void testTangentIsConsistent()
  {
    const std::vector< double > soft = { 50000., 20000., 1.0 };
    const std::vector< double > hard = { 400000., 200000., 1.0 };

    const auto card = gradedCard( { { 0.4, soft }, { 0.6, hard } }, 0, 2 );

    for ( const double alpha : { 0.0, 0.4 } ) {
      const Geometry   geometry   = inclinedGeometry( alpha );
      const Kinematics kinematics = arbitraryKinematics( 0.15 );

      const auto run = [&]( const Kinematics& k ) {
        Graded                material( "COMPRESSIBLENEOHOOKE", card.data(), card.size(), 1 );
        std::vector< double > stateVars( material.getNumberOfRequiredStateVars() );
        material.initializeYourself( stateVars.data(), stateVars.size() );
        return evaluate( material, stateVars.data(), geometry, k );
      };

      const Response base = run( kinematics );

      const double eps = 1e-7;

      // jump: the first side carries it
      for ( int k = 0; k < 3; k++ ) {
        Kinematics perturbed = kinematics;
        perturbed.U[k] += eps;
        const Response r = run( perturbed );

        const Eigen::Vector3d dForce = ( r.force - base.force ) / eps;
        throwExceptionOnFailure( checkIfEqual< double >( dForce, base.Q.col( k ), 1e-4 * base.Q.col( k ).norm() ),
                                 "dForce/dJump is not Q." );

        const Eigen::Matrix3d               dSurface = ( r.surfaceStress - base.surfaceStress ) / eps;
        const Eigen::Matrix< double, 9, 1 > expected = base.H_iRk.col( k );
        throwExceptionOnFailure( checkIfEqual< double >( Eigen::Map< const Eigen::Matrix< double, 9, 1 > >(
                                                           RowMajor3( dSurface ).data() ),
                                                         expected,
                                                         1e-4 * expected.norm() + 1e-9 ),
                                 "dSurfaceStress/dJump is not H_iRk." );
      }

      // average surface gradient: both sides move together; the tangents are the derivatives with respect to the
      // tangential part, so the perturbation of the component (k, R) is e_k (x) T e_R
      const Eigen::Matrix3d T = Eigen::Matrix3d::Identity() - geometry.normal * geometry.normal.transpose();

      for ( int k = 0; k < 3; k++ ) {
        for ( int R = 0; R < 3; R++ ) {
          Kinematics perturbed = kinematics;
          for ( int P = 0; P < 3; P++ ) {
            perturbed.surfaceGradient[3 * k + P] += eps * T( P, R );
            perturbed.surfaceGradient[9 + 3 * k + P] += eps * T( P, R );
          }
          const Response r      = run( perturbed );
          const int      column = 3 * k + R;

          const Eigen::Vector3d dForce = ( r.force - base.force ) / eps;
          throwExceptionOnFailure( checkIfEqual< double >( dForce,
                                                           base.H_ikR.col( column ),
                                                           1e-4 * base.H_ikR.col( column ).norm() + 1e-9 ),
                                   "dForce/dSurfaceGradient is not H_ikR." );

          const RowMajor3                     dSurface = ( r.surfaceStress - base.surfaceStress ) / eps;
          const Eigen::Matrix< double, 9, 1 > expected = base.A.col( column );
          throwExceptionOnFailure( checkIfEqual< double >( Eigen::Map< const Eigen::Matrix< double, 9, 1 > >(
                                                             dSurface.data() ),
                                                           expected,
                                                           1e-4 * expected.norm() + 1e-9 ),
                                   "dSurfaceStress/dSurfaceGradient is not A_iRkP." );
        }
      }
    }
  }

  /// A rigid rotation of a heterogeneous stack causes no force and no surface stress.
  void testRigidRotationGivesNoStress()
  {
    const std::vector< double > soft = { 50000., 20000., 1.0 };
    const std::vector< double > hard = { 400000., 200000., 1.0 };
    const auto                  card = gradedCard( { { 0.4, soft }, { 0.6, hard } }, 0, 2 );

    Graded                material( "COMPRESSIBLENEOHOOKE", card.data(), card.size(), 1 );
    std::vector< double > stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    const Geometry geometry = inclinedGeometry( 0.4 );

    // x = R X with a rotation about the y axis: u = ( R - I ) X, and the surface gradient of u is ( R - I )
    // restricted to the tangential directions
    const double          angle = 0.6;
    Eigen::Matrix3d       R     = Eigen::AngleAxisd( angle, Eigen::Vector3d::UnitY() ).toRotationMatrix();
    const Eigen::Matrix3d RmI   = R - Eigen::Matrix3d::Identity();

    const Eigen::Vector3d topPoint( 0.2, 0.1, h ), bottomPoint( 0.2, 0.1, 0.0 );
    const Eigen::Vector3d uTop = RmI * topPoint, uBottom = RmI * bottomPoint;

    Eigen::Matrix3d T = Eigen::Matrix3d::Identity() - geometry.normal * geometry.normal.transpose();
    // surface gradient of u: grad_s u = ( R - I ) T, row k column R
    const Eigen::Matrix3d surfaceGradientTensor = RmI * T;

    Kinematics kinematics;
    for ( int i = 0; i < 3; i++ ) {
      kinematics.U[i]     = uTop( i );
      kinematics.U[3 + i] = uBottom( i );
      for ( int j = 0; j < 3; j++ ) {
        kinematics.surfaceGradient[3 * i + j]     = surfaceGradientTensor( i, j );
        kinematics.surfaceGradient[9 + 3 * i + j] = surfaceGradientTensor( i, j );
      }
    }

    const Response response = evaluate( material, stateVars.data(), geometry, kinematics );

    throwExceptionOnFailure( response.force.norm() < 1e-8 * neoHooke[0],
                             "A rigid rotation of the stack causes a force of norm " +
                               std::to_string( response.force.norm() ) + "." );
    throwExceptionOnFailure( response.surfaceStress.norm() < 1e-8 * neoHooke[0] * h,
                             "A rigid rotation of the stack causes a surface stress." );
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

  void testInvalidCardsAndState()
  {
    const auto valid = gradedCard( { { 0.5, neoHooke }, { 0.5, neoHooke } }, 0, 2 );

    throwExceptionOnFailure( !throwsInvalidArgument(
                               [&] { Graded material( "COMPRESSIBLENEOHOOKE", valid.data(), valid.size(), 1 ); } ),
                             "A valid card was rejected." );
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { Graded material( "NOT_A_MATERIAL", valid.data(), valid.size(), 1 ); } ),
                             "An unregistered bulk material was accepted." );

    auto thickness = valid;
    thickness[0]   = -1.;
    throwExceptionOnFailure( throwsInvalidArgument( [&] {
                               Graded material( "COMPRESSIBLENEOHOOKE", thickness.data(), thickness.size(), 1 );
                             } ),
                             "A negative thickness was accepted." );

    Graded material( "COMPRESSIBLENEOHOOKE", valid.data(), valid.size(), 1 );
    throwExceptionOnFailure( material.numberOfThicknessPoints() == 4,
                             "Two layers with two Gauss points make four points." );
    throwExceptionOnFailure( checkIfEqual( material.getDensity(), neoHooke[2], 1e-12 ),
                             "The density is not the mean of the densities of the layers." );

    std::vector< double > stateVars( material.getNumberOfRequiredStateVars(), 7. );
    material.initializeYourself( stateVars.data(), stateVars.size() );
    const StateView view = material.getStateView( "correction0", stateVars.data() );
    throwExceptionOnFailure( view.stateSize == 3 && view.stateLocation[0] == 0. && view.stateLocation[2] == 0.,
                             "The corrections are not initialized to zero." );

    material.setCharacteristicElementLength( 0.1 );

    // a coincident-face mesh has no normal fibre
    const Geometry coincident = { Eigen::Vector3d( 0., 0., 1. ), Eigen::Vector3d( 0., 0., 0. ) };
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { evaluate( material, stateVars.data(), coincident, Kinematics() ); } ),
                             "Coincident faces were accepted." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testHomogeneousLayerEqualsPlainInterfaceMaterial,
                                                   testLayersActInSeriesForSmallStrains,
                                                   testTangentIsConsistent,
                                                   testRigidRotationGivesNoStress,
                                                   testInvalidCardsAndState };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
