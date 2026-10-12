#include "Marmot/MarmotGradedInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Material = MarmotGradedInterfaceMaterialHypoElastic;

  constexpr double h = 0.02;

  struct Response {
    Eigen::Vector3d force = Eigen::Vector3d::Zero();
    Eigen::Matrix< double, 3, 3, Eigen::RowMajor >
                    surfaceStress = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >::Zero();
    Eigen::VectorXd Q             = Eigen::VectorXd::Zero( 9 );
    Eigen::VectorXd Z             = Eigen::VectorXd::Zero( 81 );
    Eigen::VectorXd H             = Eigen::VectorXd::Zero( 27 );
    Eigen::VectorXd Y             = Eigen::VectorXd::Zero( 81 );
  };

  /// Material properties of the graded material: [h, mode, nGauss, K, m, (position, properties) x K]
  std::vector< double > card( int                                                              mode,
                              int                                                              nGauss,
                              const std::vector< std::pair< double, std::vector< double > > >& groups )
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
    return properties;
  }

  /// Evaluate any interface material of the family at the given increment; the state variables are updated.
  template < class InterfaceMaterial >
  Response evaluate( InterfaceMaterial& material,
                     double*            stateVars,
                     const double*      dU,
                     const double*      dSurfaceStrain,
                     const double*      normal,
                     double             timeOld  = 0.0,
                     double             dT       = 1.0,
                     const Response&    previous = Response() )
  {
    Response response;

    typename InterfaceMaterial::State         state{ typename InterfaceMaterial::Tensor3d( previous.force.data() ),
                                             typename InterfaceMaterial::Tensor33d( previous.surfaceStress.data() ),
                                             stateVars };
    typename InterfaceMaterial::Tangents      tangents{ typename InterfaceMaterial::Tensor33d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor3333d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor333d( 0.0 ),
                                                   typename InterfaceMaterial::Tensor3333d( 0.0 ) };
    typename InterfaceMaterial::Deformation   deformation{ typename InterfaceMaterial::Tensor6d( dU ),
                                                         typename InterfaceMaterial::Tensor18d( dSurfaceStrain ),
                                                         typename InterfaceMaterial::Tensor3d( normal ) };
    typename InterfaceMaterial::TimeIncrement timeIncrement{ timeOld, dT };

    material.computeStress( state, tangents, deformation, timeIncrement );

    response.force         = Marmot::mapEigenToFastor( state.force );
    response.surfaceStress = Marmot::mapEigenToFastor( state.surfaceStress );
    response.Q             = Eigen::Map< const Eigen::VectorXd >( tangents.Q_ij.data(), 9 );
    response.Z             = Eigen::Map< const Eigen::VectorXd >( tangents.Z_ijkl.data(), 81 );
    response.H             = Eigen::Map< const Eigen::VectorXd >( tangents.H_ijk.data(), 27 );
    response.Y             = Eigen::Map< const Eigen::VectorXd >( tangents.Y_ijkl.data(), 81 );
    return response;
  }

  /// Normal-normal stress of the integration point through the thickness, whose stress starts at the given entry.
  double normalStressAt( const Eigen::VectorXd& stateVars, int entry )
  {
    const Marmot::Vector6d stress = stateVars.segment< 6 >( entry );
    return Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( stress )( 2, 2 );
  }

  /// Row-major flattening of a surface stress.
  Eigen::VectorXd flatten( const Eigen::Matrix< double, 3, 3, Eigen::RowMajor >& matrix )
  {
    return Eigen::Map< const Eigen::VectorXd >( matrix.data(), 9 );
  }

  const double obliqueNormal[3] = { 1. / 3., 2. / 3., 2. / 3. };

  /// A set of arbitrary increments of the jump and of the surface gradients of the two sides.
  struct Increment {
    double dU[6];
    double dSurfaceStrain[18];
  };

  Increment makeIncrement( double scale, double tangentialShare = 1.0 )
  {
    Increment    increment;
    const double jump[6] = { 0.8, -0.3, 0.5, -0.1, 0.2, 0.05 };
    for ( int i = 0; i < 6; i++ )
      increment.dU[i] = scale * jump[i] * h;
    for ( int i = 0; i < 18; i++ )
      increment.dSurfaceStrain[i] = tangentialShare * scale * std::sin( 1.3 * i + 0.4 ) * 0.7;
    return increment;
  }

  void checkResponsesEqual( const Response& a, const Response& b, const std::string& what, double tolerance = 1e-8 )
  {
    throwExceptionOnFailure( checkIfEqual< double >( a.force, b.force, tolerance ), what + ": force differs." );
    throwExceptionOnFailure( checkIfEqual< double >( a.surfaceStress, b.surfaceStress, tolerance ),
                             what + ": surface stress differs." );
    throwExceptionOnFailure( checkIfEqual< double >( a.Q, b.Q, tolerance ), what + ": Q differs." );
    throwExceptionOnFailure( checkIfEqual< double >( a.Z, b.Z, tolerance ), what + ": Z differs." );
    throwExceptionOnFailure( checkIfEqual< double >( a.H, b.H, tolerance ), what + ": H differs." );
    throwExceptionOnFailure( checkIfEqual< double >( a.Y, b.Y, tolerance ), what + ": Y differs." );
  }

  /// A homogeneous layer, however it is described, is the interface material of the plain element.
  void testHomogeneousLayerEqualsPlainInterfaceMaterial( const std::string&           bulkName,
                                                         const std::vector< double >& bulk,
                                                         double                       relativeScale )
  {
    std::vector< double > plainCard = { bulk[0], bulk[1], h };
    plainCard.insert( plainCard.end(), bulk.begin() + 2, bulk.end() );
    MarmotInterfaceMaterialHypoElastic plain( bulkName, plainCard.data(), plainCard.size(), 1 );

    const std::vector< std::vector< double > > cards = {
      card( 0, 1, { { 1.0, bulk } } ),
      card( 0, 3, { { 0.3, bulk }, { 0.7, bulk } } ),
      card( 1, 2, { { 0.0, bulk }, { 0.4, bulk }, { 1.0, bulk } } ),
    };

    for ( const auto& properties : cards ) {
      Material graded( bulkName, properties.data(), properties.size(), 1 );

      Eigen::VectorXd plainState( plain.getNumberOfRequiredStateVars() );
      Eigen::VectorXd gradedState( graded.getNumberOfRequiredStateVars() );
      plain.initializeYourself( plainState.data(), plainState.size() );
      graded.initializeYourself( gradedState.data(), gradedState.size() );

      double   timeOld = 0.0;
      Response expected, actual; // the plain material takes its previous stress from the response
      for ( const double scale : { 1.0, 2.0, -3.0 } ) {
        const Increment increment = makeIncrement( relativeScale * scale );
        expected                  = evaluate( plain,
                             plainState.data(),
                             increment.dU,
                             increment.dSurfaceStrain,
                             obliqueNormal,
                             timeOld,
                             1.0,
                             expected );
        actual                    = evaluate( graded,
                           gradedState.data(),
                           increment.dU,
                           increment.dSurfaceStrain,
                           obliqueNormal,
                           timeOld,
                           1.0,
                           actual );
        checkResponsesEqual( actual, expected, bulkName + ": a homogeneous layer" );
        timeOld += 1.0;
      }
    }
  }

  void testHomogeneousLinearElasticLayer()
  {
    testHomogeneousLayerEqualsPlainInterfaceMaterial( "LINEARELASTIC", { 1e5, 0.3 }, 1e-4 );
  }

  void testHomogeneousVonMisesLayer()
  {
    testHomogeneousLayerEqualsPlainInterfaceMaterial( "VONMISES", { 1e5, 0.3, 100., 10., 0., 1., 2400. }, 1e-3 );
  }

  double constrainedModulus( double E, double nu )
  {
    return E * ( 1. - nu ) / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
  }

  /// Normal stiffness of a stack of layers loaded by a pure normal jump: the layers act in series.
  void testLayersActInSeries()
  {
    const double nu           = 0.3;
    const double fractions[2] = { 0.3, 0.7 }, E[2] = { 1e5, 4e5 };

    const auto properties = card( 0, 2, { { fractions[0], { E[0], nu } }, { fractions[1], { E[1], nu } } } );
    Material   material( "LINEARELASTIC", properties.data(), properties.size(), 1 );

    Eigen::VectorXd stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    const double normal[3] = { 0., 0., 1. };
    const double jump      = 1e-4;
    const double dU[6]     = { 0., 0., jump, 0., 0., 0. };
    const double dS[18]    = {};

    const Response response = evaluate( material, stateVars.data(), dU, dS, normal );

    const double compliance = fractions[0] / constrainedModulus( E[0], nu ) +
                              fractions[1] / constrainedModulus( E[1], nu );
    throwExceptionOnFailure( checkIfEqual( response.force( 2 ), jump / h / compliance, 1e-9 ),
                             "The normal traction of a stack of layers is not the one of the series model." );
    throwExceptionOnFailure( checkIfEqual( response.Q( 8 ), 1. / ( h * compliance ), 1e-9 ),
                             "The normal stiffness of a stack of layers is not the one of the series model." );

    // The two layers do not see the same strain, but they carry the same traction.
    const double traction0 = normalStressAt( stateVars, 0 );
    const double traction1 = normalStressAt( stateVars, 6 );
    throwExceptionOnFailure( checkIfEqual( traction0, traction1, 1e-9 ), "The layers do not carry the same traction." );
    throwExceptionOnFailure( checkIfEqual( traction0, response.force( 2 ), 1e-9 ),
                             "The traction of the layers is not the force." );

    // A stack without a contrast in the stiffness is a homogeneous layer
    const auto      homogeneous = card( 0, 2, { { fractions[0], { E[1], nu } }, { fractions[1], { E[1], nu } } } );
    Material        homogeneousMaterial( "LINEARELASTIC", homogeneous.data(), homogeneous.size(), 1 );
    Eigen::VectorXd homogeneousState( homogeneousMaterial.getNumberOfRequiredStateVars() );
    homogeneousMaterial.initializeYourself( homogeneousState.data(), homogeneousState.size() );
    const Response homogeneousResponse = evaluate( homogeneousMaterial, homogeneousState.data(), dU, dS, normal );
    throwExceptionOnFailure( checkIfEqual( homogeneousResponse.force( 2 ),
                                           constrainedModulus( E[1], nu ) * jump / h,
                                           1e-9 ),
                             "A stack without contrast is not the homogeneous layer." );
  }

  /// The tangential stiffness of a stack of layers is the one of the parallel model.
  void testLayersActInParallelInThePlane()
  {
    const double nu           = 0.0;
    const double fractions[2] = { 0.25, 0.75 }, E[2] = { 1e5, 3e5 };

    const auto properties = card( 0, 1, { { fractions[0], { E[0], nu } }, { fractions[1], { E[1], nu } } } );
    Material   material( "LINEARELASTIC", properties.data(), properties.size(), 1 );

    Eigen::VectorXd stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    // tangential strain eps_xx = 1e-3 through the average surface gradient, normal (0,0,1)
    const double normal[3] = { 0., 0., 1. };
    const double dU[6]     = {};
    double       dS[18]    = {};
    dS[0]                  = 1e-3;
    dS[9]                  = 1e-3;

    const Response response = evaluate( material, stateVars.data(), dU, dS, normal );
    const double   expected = ( fractions[0] * E[0] + fractions[1] * E[1] ) * 1e-3;
    throwExceptionOnFailure( checkIfEqual( response.surfaceStress( 0, 0 ), h * expected, 1e-9 ),
                             "The in-plane stress of a stack of layers is not the one of the parallel model." );
  }

  /// A linear grading of the modulus: the normal compliance is the integral of 1 / M(y) through the thickness.
  void testLinearGradingOfTheModulus()
  {
    const double nu = 0.25, E0 = 1e5, E1 = 3e5;

    const auto properties = card( 1, 16, { { 0.0, { E0, nu } }, { 1.0, { E1, nu } } } );
    Material   material( "LINEARELASTIC", properties.data(), properties.size(), 1 );

    Eigen::VectorXd stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    const double normal[3] = { 0., 0., 1. };
    const double jump      = 2e-4;
    const double dU[6]     = { 0., 0., jump, 0., 0., 0. };
    const double dS[18]    = {};

    const Response response = evaluate( material, stateVars.data(), dU, dS, normal );

    const double compliance = std::log( E1 / E0 ) / ( constrainedModulus( 1., nu ) * ( E1 - E0 ) );
    throwExceptionOnFailure( checkIfEqual( response.force( 2 ), jump / h / compliance, 1e-9 ),
                             "The normal traction of a graded layer is not the one of the integral of the "
                             "compliance." );

    // all points carry the traction of the layer, although only the first and the last are compared
    const double tractionBottom = normalStressAt( stateVars, 0 );
    const double tractionTop    = normalStressAt( stateVars, 6 * ( material.numberOfThicknessPoints() - 1 ) );
    throwExceptionOnFailure( checkIfEqual( tractionBottom, tractionTop, 1e-9 ),
                             "The graded layer does not carry a constant traction." );
  }

  /// The tangents are the derivatives of the response, also where the layers yield differently.
  void testTangentIsConsistent()
  {
    const std::vector< double > soft = { 1e5, 0.3, 40., 10., 0., 1., 2400. };
    const std::vector< double > hard = { 2e5, 0.3, 150., 10., 0., 1., 2400. };

    const auto properties = card( 0, 2, { { 0.4, soft }, { 0.6, hard } } );
    Material   material( "VONMISES", properties.data(), properties.size(), 1 );

    Eigen::VectorXd stateVars( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), stateVars.size() );

    // load into the plastic range of the soft layer, then check the tangent of the next step
    const Increment first = makeIncrement( 2e-3, 0.05 );
    evaluate( material, stateVars.data(), first.dU, first.dSurfaceStrain, obliqueNormal, 0.0, 1.0 );

    // the soft layers have hardened, the hard layers have not yet yielded: the layers have different histories
    const double kappaSoft = material.getStateView( "baseMaterialStateVars0", stateVars.data() ).stateLocation[0];
    const double kappaHard = material.getStateView( "baseMaterialStateVars2", stateVars.data() ).stateLocation[0];
    throwExceptionOnFailure( kappaSoft > 0. && kappaHard == 0.,
                             "The test does not load the soft layer into the plastic range alone." );

    const Eigen::VectorXd committed = stateVars;

    const Increment second = makeIncrement( 5e-4, 0.05 );
    Eigen::VectorXd trial  = committed;
    const Response base = evaluate( material, trial.data(), second.dU, second.dSurfaceStrain, obliqueNormal, 1.0, 1.0 );

    const double perturbation = 1e-7;

    const auto perturbed = [&]( const Increment& increment ) {
      Eigen::VectorXd state = committed;
      return evaluate( material, state.data(), increment.dU, increment.dSurfaceStrain, obliqueNormal, 1.0, 1.0 );
    };

    // the jump enters through the first side, the average surface gradient through both sides
    for ( int j = 0; j < 3; j++ ) {
      Increment plus = second;
      plus.dU[j] += perturbation;

      const Response        response = perturbed( plus );
      const Eigen::Vector3d dForce   = ( response.force - base.force ) / perturbation;
      const Eigen::Vector3d
        expectedForce = Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( base.Q.data() ).col( j );
      throwExceptionOnFailure( checkIfEqual< double >( dForce, expectedForce, 1e-4 * expectedForce.norm() ),
                               "dForce/dJump is not the tangent Q." );

      const Eigen::VectorXd dSurface = flatten( ( response.surfaceStress - base.surfaceStress ) / perturbation );
      const Eigen::Matrix< double, 9, 3, Eigen::RowMajor >
        HT = Eigen::Map< const Eigen::Matrix< double, 3, 9, Eigen::RowMajor > >( base.H.data() ).transpose();
      throwExceptionOnFailure( checkIfEqual< double >( dSurface, HT.col( j ), 1e-4 * HT.col( j ).norm() + 1e-9 ),
                               "dSurfaceStress/dJump is not the transpose of the tangent H." );
    }

    for ( int k = 0; k < 9; k++ ) {
      Increment plus = second;
      plus.dSurfaceStrain[k] += perturbation;
      plus.dSurfaceStrain[9 + k] += perturbation;

      const Response response = perturbed( plus );

      const Eigen::Vector3d dForce = ( response.force - base.force ) / perturbation;
      const Eigen::Vector3d
        expectedForce = Eigen::Map< const Eigen::Matrix< double, 3, 9, Eigen::RowMajor > >( base.H.data() ).col( k );
      throwExceptionOnFailure( checkIfEqual< double >( dForce, expectedForce, 1e-4 * expectedForce.norm() + 1e-9 ),
                               "dForce/dSurfaceGradient is not the tangent H." );

      const Eigen::VectorXd dSurface        = flatten( ( response.surfaceStress - base.surfaceStress ) / perturbation );
      const Eigen::VectorXd expectedSurface = ( Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >(
                                                  base.Z.data() ) +
                                                Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >(
                                                  base.Y.data() ) )
                                                .col( k );
      throwExceptionOnFailure( checkIfEqual< double >( dSurface,
                                                       expectedSurface,
                                                       1e-4 * expectedSurface.norm() + 1e-9 ),
                               "dSurfaceStress/dSurfaceGradient is not the tangent Z + Y." );
    }
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

  void testInvalidCardsAreRejected()
  {
    const std::vector< double > bulk = { 1e5, 0.3 };

    const auto rejected = [&]( std::vector< double > properties, const std::string& materialName = "LINEARELASTIC" ) {
      return throwsInvalidArgument(
        [&] { Material material( materialName, properties.data(), properties.size(), 1 ); } );
    };

    const auto valid = card( 1, 2, { { 0.0, bulk }, { 1.0, bulk } } );
    throwExceptionOnFailure( !rejected( valid ), "A valid card was rejected." );

    auto thickness = valid;
    thickness[0]   = 0.;
    throwExceptionOnFailure( rejected( thickness ), "A vanishing thickness was accepted." );

    auto mode = valid;
    mode[1]   = 2.;
    throwExceptionOnFailure( rejected( mode ), "An unknown mode was accepted." );

    auto gauss = valid;
    gauss[2]   = 0.;
    throwExceptionOnFailure( rejected( gauss ), "Without a Gauss point was accepted." );

    throwExceptionOnFailure( rejected( card( 1, 2, { { 0.0, bulk } } ) ), "A single station was accepted." );
    throwExceptionOnFailure( rejected( card( 1, 2, { { 0.1, bulk }, { 1.0, bulk } } ) ),
                             "A first station away from 0 was accepted." );
    throwExceptionOnFailure( rejected( card( 1, 2, { { 0.0, bulk }, { 0.9, bulk } } ) ),
                             "A last station away from 1 was accepted." );
    throwExceptionOnFailure( rejected( card( 1, 2, { { 0.0, bulk }, { 0.5, bulk }, { 0.4, bulk }, { 1.0, bulk } } ) ),
                             "Stations that do not increase were accepted." );
    throwExceptionOnFailure( rejected( card( 0, 2, { { 0.0, bulk }, { 1.0, bulk } } ) ),
                             "A layer without thickness was accepted." );

    auto truncated = valid;
    truncated.pop_back();
    throwExceptionOnFailure( rejected( truncated ), "A card with a missing property was accepted." );
    throwExceptionOnFailure( rejected( { h, 0., 1. } ), "A card that is too short was accepted." );
    throwExceptionOnFailure( rejected( valid, "NOT_A_REGISTERED_MATERIAL" ),
                             "An unregistered bulk material was accepted." );
  }

  void testStateAndDensity()
  {
    const std::vector< double > bulk       = { 2e5, 0.2, 0.5, 0.1, 10., 0.0001, 1., 2400. };
    const auto                  properties = card( 0, 2, { { 1.0, bulk }, { 1.0, bulk } } );
    Material                    material( "LINEARVISCOELASTICWIECHERT", properties.data(), properties.size(), 1 );

    throwExceptionOnFailure( material.numberOfThicknessPoints() == 4,
                             "Two layers with two Gauss points make four points." );

    Eigen::VectorXd stateVars = Eigen::VectorXd::Constant( material.getNumberOfRequiredStateVars(), 7. );
    material.initializeYourself( stateVars.data(), stateVars.size() );
    throwExceptionOnFailure( stateVars.head( 6 ).isZero(), "The stress is not initialized to zero." );

    const StateView view = material.getStateView( "stress0", stateVars.data() );
    throwExceptionOnFailure( view.stateLocation == stateVars.data() && view.stateSize == 6,
                             "The stress is not a state view." );

    material.setCharacteristicElementLength( 0.1 );
    throwExceptionOnFailure( material.characteristicElementLength == 0.1, "The characteristic length is not stored." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testHomogeneousLinearElasticLayer,
                                                   testHomogeneousVonMisesLayer,
                                                   testLayersActInSeries,
                                                   testLayersActInParallelInThePlane,
                                                   testLinearGradingOfTheModulus,
                                                   testTangentIsConsistent,
                                                   testInvalidCardsAreRejected,
                                                   testStateAndDensity };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
