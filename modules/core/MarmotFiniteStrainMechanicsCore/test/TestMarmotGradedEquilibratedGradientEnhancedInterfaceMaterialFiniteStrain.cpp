#include "Marmot/MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain.h"
#include "Marmot/MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain.h"
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

  using Graded       = MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain;
  using Equilibrated = MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain;

  constexpr double h = 0.1;

  const std::string           bulkName = "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE";
  const std::vector< double > damaged  = { 3500.0, 1500.0, 1.0e-3, 5.0e-3, 0.2, 1.0 };

  /// Response of the layer: the generalized stress, the fields of the stations and the derivatives.
  struct Response {
    Eigen::Matrix< double, 21, 1 >  p = Eigen::Matrix< double, 21, 1 >::Zero();
    Eigen::Matrix< double, 21, 21 > K = Eigen::Matrix< double, 21, 21 >::Zero();
    Eigen::MatrixXd                 dStress_dNonlocal;
    Eigen::MatrixXd                 dLocal_dStrain;
    Eigen::MatrixXd                 dLocal_dNonlocal;
    Eigen::VectorXd                 local;
    Eigen::VectorXd                 radius;
  };

  /// Total kinematics of an interface point and the nonlocal fields of the stations.
  struct Kinematics {
    Eigen::Matrix< double, 6, 1 >  U         = Eigen::Matrix< double, 6, 1 >::Zero();
    Eigen::Matrix< double, 18, 1 > gradient  = Eigen::Matrix< double, 18, 1 >::Zero();
    Eigen::Vector3d                normal    = Eigen::Vector3d( 0., 0., 1. );
    Eigen::Vector3d                connector = Eigen::Vector3d( 0., 0., h );
    Eigen::VectorXd                nonlocal;
  };

  /// Properties of the graded material: [h, 2, 1, K, m, (position, properties) x K, flag].
  std::vector< double > card( const std::vector< double >& positions, int tangentialGradient )
  {
    std::vector< double > properties = { h, 2., 1., double( positions.size() ), double( damaged.size() ) };
    for ( const double position : positions ) {
      properties.push_back( position );
      properties.insert( properties.end(), damaged.begin(), damaged.end() );
    }
    if ( tangentialGradient >= 0 )
      properties.push_back( double( tangentialGradient ) );
    return properties;
  }

  Kinematics makeKinematics( double scale, bool inclined, int nStations )
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
    kinematics.nonlocal.resize( nStations );
    for ( int q = 0; q < nStations; q++ )
      kinematics.nonlocal( q ) = 0.002 * ( 1.0 + 0.7 * q );
    return kinematics;
  }

  void extractBlocks( Response&     response,
                      const double* ww,
                      const double* wAp,
                      const double* wAm,
                      const double* Apw,
                      const double* ApAp,
                      const double* ApAm,
                      const double* Amw,
                      const double* AmAp,
                      const double* AmAm )
  {
    using RowMajor   = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >;
    const auto block = [&]( const double* data, int rows, int cols ) {
      return Eigen::Map< const RowMajor >( data, rows, cols );
    };
    response.K.block< 3, 3 >( 0, 0 )   = block( ww, 3, 3 );
    response.K.block< 3, 9 >( 0, 3 )   = block( wAp, 3, 9 );
    response.K.block< 3, 9 >( 0, 12 )  = block( wAm, 3, 9 );
    response.K.block< 9, 3 >( 3, 0 )   = block( Apw, 9, 3 );
    response.K.block< 9, 9 >( 3, 3 )   = block( ApAp, 9, 9 );
    response.K.block< 9, 9 >( 3, 12 )  = block( ApAm, 9, 9 );
    response.K.block< 9, 3 >( 12, 0 )  = block( Amw, 9, 3 );
    response.K.block< 9, 9 >( 12, 3 )  = block( AmAp, 9, 9 );
    response.K.block< 9, 9 >( 12, 12 ) = block( AmAm, 9, 9 );
  }

  struct GradedPoint {
    std::unique_ptr< Graded > material;
    std::vector< double >     stateVars;

    explicit GradedPoint( const std::vector< double >& properties )
    {
      material = std::make_unique< Graded >( bulkName, properties.data(), static_cast< int >( properties.size() ), 1 );
      stateVars.assign( material->getNumberOfRequiredStateVars(), 0.0 );
      material->initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );
    }
  };

  Response evaluate( Graded& material, double* stateVars, const Kinematics& kinematics, double timeOld = 0. )
  {
    Graded::State state{ Graded::Tensor3d( 0.0 ),
                         Graded::Tensor33d( 0.0 ),
                         Graded::Tensor33d( 0.0 ),
                         stateVars,
                         Eigen::VectorXd(),
                         Eigen::VectorXd() };

    Graded::Tangents tangents{ Graded::Tensor33d( 0.0 ),
                               Graded::Tensor333d( 0.0 ),
                               Graded::Tensor333d( 0.0 ),
                               Graded::Tensor333d( 0.0 ),
                               Graded::Tensor3333d( 0.0 ),
                               Graded::Tensor3333d( 0.0 ),
                               Graded::Tensor333d( 0.0 ),
                               Graded::Tensor3333d( 0.0 ),
                               Graded::Tensor3333d( 0.0 ),
                               {},
                               {},
                               {} };

    Graded::Deformation deformation{ Graded::Tensor6d( kinematics.U.data() ),
                                     Graded::Tensor18d( kinematics.gradient.data() ),
                                     Graded::Tensor3d( kinematics.normal.data() ),
                                     Graded::Tensor3d( kinematics.connector.data() ),
                                     kinematics.nonlocal };

    material.computeStress( state, tangents, deformation, { timeOld, 1.0 } );

    Response response;
    response.p.segment< 3 >( 0 )  = Eigen::Map< const Eigen::Vector3d >( state.generalizedForce.data() );
    response.p.segment< 9 >( 3 )  = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( state.surfaceStressPlus.data() );
    response.p.segment< 9 >( 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >(
      state.surfaceStressMinus.data() );
    extractBlocks( response,
                   tangents.Q_ww.data(),
                   tangents.Q_wAp.data(),
                   tangents.Q_wAm.data(),
                   tangents.Q_Apw.data(),
                   tangents.Q_ApAp.data(),
                   tangents.Q_ApAm.data(),
                   tangents.Q_Amw.data(),
                   tangents.Q_AmAp.data(),
                   tangents.Q_AmAm.data() );
    response.dStress_dNonlocal = tangents.dStress_dNonlocal;
    response.dLocal_dStrain    = tangents.dLocal_dStrain;
    response.dLocal_dNonlocal  = tangents.dLocal_dNonlocal;
    response.local             = state.localFields;
    response.radius            = state.nonlocalRadii;
    return response;
  }

  /// With two stations at the faces and without the tangential gradient, the material is the face-fields material.
  void testTwoStationsAreFaceFieldsMaterial()
  {
    const std::vector< double > equilibratedCard =
      { damaged[0], damaged[1], h, damaged[2], damaged[3], damaged[4], damaged[5] };

    for ( const bool inclined : { false, true } )
      for ( const double scale : { 1e-3, 0.05 } ) {
        Equilibrated          equilibrated( bulkName,
                                   equilibratedCard.data(),
                                   static_cast< int >( equilibratedCard.size() ),
                                   1 );
        std::vector< double > equilibratedState( equilibrated.getNumberOfRequiredStateVars(), 0.0 );
        equilibrated.initializeYourself( equilibratedState.data(), static_cast< int >( equilibratedState.size() ) );

        GradedPoint graded( card( { 0.0, 1.0 }, 0 ) );

        const Kinematics kinematics = makeKinematics( scale, inclined, 2 );

        Equilibrated::FaceFieldsState state;
        state.generalizedForce   = Equilibrated::Tensor3d( 0.0 );
        state.surfaceStressPlus  = Equilibrated::Tensor33d( 0.0 );
        state.surfaceStressMinus = Equilibrated::Tensor33d( 0.0 );
        state.stateVars          = equilibratedState.data();

        Equilibrated::FaceFieldsTangents    tangents;
        Equilibrated::FaceFieldsDeformation deformation;
        deformation.U               = Equilibrated::Tensor6d( kinematics.U.data() );
        deformation.surfaceGradient = Equilibrated::Tensor18d( kinematics.gradient.data() );
        deformation.normal          = Equilibrated::Tensor3d( kinematics.normal.data() );
        deformation.connector       = Equilibrated::Tensor3d( kinematics.connector.data() );
        // the first station is the bottom face and the second one the top face
        deformation.nonlocalFieldBottom = kinematics.nonlocal( 0 );
        deformation.nonlocalFieldTop    = kinematics.nonlocal( 1 );

        equilibrated.computeStress( state, tangents, deformation, { 0.0, 1.0 } );

        const Response actual = evaluate( *graded.material, graded.stateVars.data(), kinematics );

        Eigen::Matrix< double, 21, 1 > expectedP;
        expectedP.segment< 3 >( 0 ) = Eigen::Map< const Eigen::Vector3d >( state.generalizedForce.data() );
        expectedP.segment< 9 >( 3 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >(
          state.surfaceStressPlus.data() );
        expectedP.segment< 9 >( 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >(
          state.surfaceStressMinus.data() );

        const std::string context = std::string( inclined ? "inclined" : "straight" ) + " at the scale " +
                                    std::to_string( scale );

        throwExceptionOnFailure( checkIfEqual< double >( actual.p, expectedP, 1e-8 * ( 1. + expectedP.norm() ) ),
                                 "The generalized stress is not the one of the face-fields material, " + context );

        // station 1 is the top face and station 0 the bottom face
        const Eigen::Matrix< double, 21, 1 > dTop = Eigen::Map< const Eigen::Matrix< double, 21, 1 > >(
          tangents.dGeneralizedStress_dNonlocalTop.data() );
        const Eigen::Matrix< double, 21, 1 > dBottom = Eigen::Map< const Eigen::Matrix< double, 21, 1 > >(
          tangents.dGeneralizedStress_dNonlocalBottom.data() );
        throwExceptionOnFailure( checkIfEqual< double >( Eigen::VectorXd( actual.dStress_dNonlocal.col( 1 ) ),
                                                         Eigen::VectorXd( dTop ),
                                                         1e-8 * ( 1. + dTop.norm() ) ),
                                 "The derivative of the stress with respect to the top field differs, " + context );
        throwExceptionOnFailure( checkIfEqual< double >( Eigen::VectorXd( actual.dStress_dNonlocal.col( 0 ) ),
                                                         Eigen::VectorXd( dBottom ),
                                                         1e-8 * ( 1. + dBottom.norm() ) ),
                                 "The derivative of the stress with respect to the bottom field differs, " + context );

        const double Ltop = state.localFieldTop, Lbottom = state.localFieldBottom;
        throwExceptionOnFailure( std::abs( actual.local( 1 ) - Ltop ) < 1e-10 * ( 1. + std::abs( Ltop ) ) &&
                                   std::abs( actual.local( 0 ) - Lbottom ) < 1e-10 * ( 1. + std::abs( Lbottom ) ),
                                 "The local fields differ from the face-fields material, " + context );

        const Eigen::Matrix< double, 1, 21 > dLTop = Eigen::Map< const Eigen::Matrix< double, 1, 21 > >(
          tangents.dLocalFieldTop_dGeneralizedStrain.data() );
        throwExceptionOnFailure( ( actual.dLocal_dStrain.row( 1 ) - dLTop ).norm() < 1e-8 * ( 1. + dLTop.norm() ),
                                 "The derivative of the local field of the top face differs, " + context );
        throwExceptionOnFailure( std::abs( actual.dLocal_dNonlocal( 1, 1 ) - tangents.dLocalFieldTop_dNonlocalTop ) <
                                     1e-8 &&
                                   std::abs( actual.dLocal_dNonlocal( 0, 1 ) -
                                             tangents.dLocalFieldBottom_dNonlocalTop ) < 1e-8,
                                 "The derivative of the local field with respect to the fields differs, " + context );
      }
  }

  /// Derivatives of the response by central differences, with respect to the strain and to the fields.
  void checkDerivatives( const std::vector< double >& properties,
                         const Kinematics&            kinematics,
                         double                       tolerance,
                         const std::string&           what )
  {
    GradedPoint    reference( properties );
    const Response response = evaluate( *reference.material, reference.stateVars.data(), kinematics );
    const int      nq       = static_cast< int >( kinematics.nonlocal.size() );

    const auto evaluateAt = [&]( const Kinematics& perturbed ) {
      GradedPoint point( properties );
      return evaluate( *point.material, point.stateVars.data(), perturbed );
    };

    const double stressScale = response.K.cwiseAbs().maxCoeff();

    // with respect to the generalized strain
    for ( int j = 0; j < 21; j++ ) {
      const double step = 1e-7 * ( j < 3 ? h : 1.0 );
      Response     sides[2];
      for ( int side = 0; side < 2; side++ ) {
        Kinematics   perturbed = kinematics;
        const double delta     = side == 0 ? -step : step;
        if ( j < 3 )
          perturbed.U( j ) += delta;
        else
          perturbed.gradient( j - 3 ) += delta;
        sides[side] = evaluateAt( perturbed );
      }
      const Eigen::VectorXd dStress = ( sides[1].p - sides[0].p ) / ( 2. * step );
      const Eigen::VectorXd dLocal  = ( sides[1].local - sides[0].local ) / ( 2. * step );

      throwExceptionOnFailure( ( dStress - response.K.col( j ) ).cwiseAbs().maxCoeff() < tolerance * stressScale,
                               what +
                                 ": the tangent is not the derivative of the stress with respect to the "
                                 "strain component " +
                                 std::to_string( j ) );
      const double localScale = response.dLocal_dStrain.cwiseAbs().maxCoeff() + 1e-12;
      throwExceptionOnFailure( ( dLocal - response.dLocal_dStrain.col( j ) ).cwiseAbs().maxCoeff() <
                                 tolerance * localScale,
                               what + ": the derivative of the local fields is wrong for the strain component " +
                                 std::to_string( j ) );
    }

    // with respect to the nonlocal fields
    for ( int f = 0; f < nq; f++ ) {
      const double step = 1e-7;
      Response     sides[2];
      for ( int side = 0; side < 2; side++ ) {
        Kinematics perturbed = kinematics;
        perturbed.nonlocal( f ) += side == 0 ? -step : step;
        sides[side] = evaluateAt( perturbed );
      }
      const Eigen::VectorXd dStress = ( sides[1].p - sides[0].p ) / ( 2. * step );
      const Eigen::VectorXd dLocal  = ( sides[1].local - sides[0].local ) / ( 2. * step );

      const double scale = response.dStress_dNonlocal.cwiseAbs().maxCoeff() + 1e-12;
      throwExceptionOnFailure( ( dStress - response.dStress_dNonlocal.col( f ) ).cwiseAbs().maxCoeff() <
                                 tolerance * scale,
                               what + ": the derivative of the stress with respect to the field " +
                                 std::to_string( f ) + " is wrong" );
      const double localScale = response.dLocal_dNonlocal.cwiseAbs().maxCoeff() + 1e-12;
      throwExceptionOnFailure( ( dLocal - response.dLocal_dNonlocal.col( f ) ).cwiseAbs().maxCoeff() <
                                 tolerance * localScale,
                               what + ": the derivative of the local fields with respect to the field " +
                                 std::to_string( f ) + " is wrong" );
    }
  }

  void testDerivativesAreConsistent()
  {
    for ( const bool inclined : { false, true } ) {
      // without the closure of the tangential gradient the derivatives are exact
      checkDerivatives( card( { 0.0, 0.5, 1.0 }, 0 ), makeKinematics( 0.05, inclined, 3 ), 1e-5, "without closure" );
      // with it, they neglect the dependence of the derivative of the correction on the state and on the fields,
      // which is the second derivative of the damage
      checkDerivatives( card( { 0.0, 0.4, 1.0 }, 1 ), makeKinematics( 1e-3, inclined, 3 ), 2e-1, "with closure" );
    }
  }

  /// Damage differs through the thickness when the fields of the stations differ.
  void testDamageVariesThroughTheThickness()
  {
    GradedPoint point( card( { 0.0, 0.5, 1.0 }, 1 ) );

    Kinematics kinematics = makeKinematics( 0.05, false, 3 );
    kinematics.nonlocal << 0.0, 0.002, 0.004;

    const Response response = evaluate( *point.material, point.stateVars.data(), kinematics );

    throwExceptionOnFailure( response.local.size() == 3, "A local field is required for every station." );
    throwExceptionOnFailure( response.radius.minCoeff() > 0., "A nonlocal radius is required for every station." );
  }

  void testWrongNumberOfFieldsIsRefused()
  {
    GradedPoint point( card( { 0.0, 0.5, 1.0 }, 1 ) );
    Kinematics  kinematics = makeKinematics( 0.05, false, 2 );

    bool refused = false;
    try {
      evaluate( *point.material, point.stateVars.data(), kinematics );
    }
    catch ( const std::invalid_argument& ) {
      refused = true;
    }
    throwExceptionOnFailure( refused, "A wrong number of nonlocal fields is not refused." );
  }

  void testModeMustPlaceThePointsAtTheStations()
  {
    std::vector< double > properties   = card( { 0.0, 1.0 }, -1 );
    properties[1]                      = 0.0; // layers
    properties[5]                      = 0.5;
    properties[5 + 1 + damaged.size()] = 0.5;

    bool refused = false;
    try {
      Graded material( bulkName, properties.data(), static_cast< int >( properties.size() ), 1 );
    }
    catch ( const std::invalid_argument& ) {
      refused = true;
    }
    throwExceptionOnFailure( refused, "A mode other than 2 is not refused." );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{ testTwoStationsAreFaceFieldsMaterial,
                                                       testDerivativesAreConsistent,
                                                       testDamageVariesThroughTheThickness,
                                                       testWrongNumberOfFieldsIsRefused,
                                                       testModeMustPlaceThePointsAtTheStations };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
