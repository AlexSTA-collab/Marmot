#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotGradedXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Graded       = MarmotGradedXInterfaceMaterialHypoElastic;
  using Equilibrated = MarmotEquilibratedXInterfaceMaterialHypoElastic;

  constexpr double h = 0.02;

  /// Generalized stress p = ( force, S+, S- ) and the tangent over the generalized strain ( w, A+, A- ).
  struct Response {
    Eigen::Matrix< double, 21, 1 >  p         = Eigen::Matrix< double, 21, 1 >::Zero();
    Eigen::Matrix< double, 21, 1 >  resultant = Eigen::Matrix< double, 21, 1 >::Zero(); ///< moments of the stress
    Eigen::Matrix< double, 21, 21 > K         = Eigen::Matrix< double, 21, 21 >::Zero();
  };

  /// Generalized strain increment ( w, A+, A- ) and the geometry of an interface point.
  struct Kinematics {
    Eigen::Matrix< double, 21, 1 > x         = Eigen::Matrix< double, 21, 1 >::Zero();
    Eigen::Vector3d                normal    = Eigen::Vector3d( 0., 0., 1. );
    Eigen::Vector3d                connector = Eigen::Vector3d::Zero();
  };

  /// Material properties of the graded material: [h, mode, nGauss, K, m, (position, properties) x K]
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

  /// Evaluate a material with the contract of the equilibrated material for the increment; the state is updated.
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

    // the displacement of the top side is the jump, the one of the bottom side vanishes
    double dU[6] = { kinematics.x( 0 ), kinematics.x( 1 ), kinematics.x( 2 ), 0., 0., 0. };
    double dA[18];
    for ( int i = 0; i < 18; i++ )
      dA[i] = kinematics.x( 3 + i );

    typename InterfaceMaterial::Deformation   deformation{ typename InterfaceMaterial::Tensor6d( dU ),
                                                         typename InterfaceMaterial::Tensor18d( dA ),
                                                         typename InterfaceMaterial::Tensor3d(
                                                           kinematics.normal.data() ),
                                                         typename InterfaceMaterial::Tensor3d(
                                                           kinematics.connector.data() ) };
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

  /// A new graded material and its initialized state.
  struct GradedPoint {
    std::unique_ptr< Graded > material;
    std::vector< double >     stateVars;

    GradedPoint( const std::string& bulkName, const std::vector< double >& properties )
    {
      material = std::make_unique< Graded >( bulkName, properties.data(), static_cast< int >( properties.size() ), 1 );
      stateVars.assign( material->getNumberOfRequiredStateVars(), 0.0 );
      material->initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );
    }
  };

  /// An increment with every component of the generalized strain nonzero, scaled.
  Kinematics makeKinematics( double scale, bool withConnector )
  {
    Kinematics kinematics;
    kinematics.normal = Eigen::Vector3d( 1., 2., 2. ) / 3.;
    for ( int i = 0; i < 3; i++ )
      kinematics.x( i ) = scale * h * std::sin( 0.7 * i + 0.3 );
    for ( int i = 3; i < 21; i++ )
      kinematics.x( i ) = scale * std::sin( 1.1 * i + 0.4 );

    // surface gradients are tangential: remove the normal part of every row
    for ( int block = 0; block < 2; block++ ) {
      Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > A( kinematics.x.data() + 3 + 9 * block );
      A = A - ( A * kinematics.normal ) * kinematics.normal.transpose();
    }
    if ( withConnector ) {
      const Eigen::Vector3d tangent = Eigen::Vector3d( 2., -2., 1. ) / 3.;
      kinematics.connector          = 1.1 * h * kinematics.normal + 0.15 * h * tangent;
    }
    return kinematics;
  }

  /// With two points at the faces, the material is the equilibrated two-sided material.
  void testTwoFacesAreEquilibratedMaterial( const std::string&           bulkName,
                                            const std::vector< double >& bulk,
                                            double                       scale,
                                            bool                         withConnector )
  {
    // [E, nu, h, remaining properties] of the equilibrated material
    std::vector< double > equilibratedCard = { bulk[0], bulk[1], h };
    equilibratedCard.insert( equilibratedCard.end(), bulk.begin() + 2, bulk.end() );

    Equilibrated equilibrated( bulkName, equilibratedCard.data(), static_cast< int >( equilibratedCard.size() ), 1 );
    std::vector< double > equilibratedState( equilibrated.getNumberOfRequiredStateVars(), 0.0 );
    equilibrated.initializeYourself( equilibratedState.data(), static_cast< int >( equilibratedState.size() ) );

    GradedPoint graded( bulkName, card( 2, 1, { { 0.0, bulk }, { 1.0, bulk } }, 0 ) );

    for ( int step = 0; step < 3; step++ ) {
      const Kinematics kinematics = makeKinematics( scale * ( step + 1 ), withConnector );

      const Response expected = evaluate( equilibrated, equilibratedState.data(), kinematics, step );
      const Response actual   = evaluate( *graded.material, graded.stateVars.data(), kinematics, step );

      if ( std::getenv( "DEBUG_TWO_FACES" ) )
        std::cout << bulkName << " connector " << withConnector << " step " << step << "\nactual   "
                  << actual.p.transpose() << "\nexpected " << expected.p.transpose() << std::endl;
      throwExceptionOnFailure( checkIfEqual< double >( actual.p, expected.p, 1e-8 * ( 1. + expected.p.norm() ) ),
                               bulkName + ": the generalized stress is not the one of the equilibrated material." );
      throwExceptionOnFailure( checkIfEqual< double >( Eigen::VectorXd(
                                                         Eigen::Map< const Eigen::VectorXd >( actual.K.data(), 441 ) ),
                                                       Eigen::VectorXd(
                                                         Eigen::Map< const Eigen::VectorXd >( expected.K.data(),
                                                                                              441 ) ),
                                                       1e-8 * ( 1. + expected.K.norm() ) ),
                               bulkName + ": the tangent is not the one of the equilibrated material." );
    }
  }

  void testElasticTwoFacesAreEquilibratedMaterial()
  {
    testTwoFacesAreEquilibratedMaterial( "LINEARELASTIC", { 4000., 0.3 }, 1e-3, false );
    testTwoFacesAreEquilibratedMaterial( "LINEARELASTIC", { 4000., 0.3 }, 1e-3, true );
  }

  void testPlasticTwoFacesAreEquilibratedMaterial()
  {
    const std::vector< double > vonMises = { 1e4, 0.3, 40., 10., 0., 1., 2400. };
    testTwoFacesAreEquilibratedMaterial( "VONMISES", vonMises, 2e-3, false );
    testTwoFacesAreEquilibratedMaterial( "VONMISES", vonMises, 2e-3, true );
  }

  /// Tangent of a response by central differences of the generalized stress.
  Eigen::Matrix< double, 21, 21 > numericalTangent( const std::string&           bulkName,
                                                    const std::vector< double >& properties,
                                                    const Kinematics&            kinematics,
                                                    const GradedPoint&           committed )
  {
    Eigen::Matrix< double, 21, 21 > K = Eigen::Matrix< double, 21, 21 >::Zero();
    for ( int j = 0; j < 21; j++ ) {
      const double step = ( j < 3 ? h : 1.0 ) * 1e-7;

      Response sides[2];
      for ( int side = 0; side < 2; side++ ) {
        GradedPoint point( bulkName, properties );
        point.stateVars      = committed.stateVars;
        Kinematics perturbed = kinematics;
        perturbed.x( j ) += ( side == 0 ? -step : step );
        sides[side] = evaluate( *point.material, point.stateVars.data(), perturbed );
      }
      K.col( j ) = ( sides[1].p - sides[0].p ) / ( 2. * step );
    }
    return K;
  }

  void checkTangent( const std::string&           bulkName,
                     const std::vector< double >& properties,
                     const Kinematics&            kinematics,
                     double                       tolerance,
                     const std::string&           what )
  {
    GradedPoint point( bulkName, properties );
    const auto  committed = point.stateVars;
    GradedPoint reference( bulkName, properties );
    reference.stateVars     = committed;
    const Response response = evaluate( *reference.material, reference.stateVars.data(), kinematics );

    const auto   numerical = numericalTangent( bulkName, properties, kinematics, point );
    const double scale     = response.K.cwiseAbs().maxCoeff();

    for ( int i = 0; i < 21; i++ )
      for ( int j = 0; j < 21; j++ )
        throwExceptionOnFailure( std::abs( numerical( i, j ) - response.K( i, j ) ) < tolerance * scale,
                                 what + ": the tangent is not the derivative of the generalized stress at (" +
                                   std::to_string( i ) + "," + std::to_string( j ) + ")" );
  }

  void testTangentOfElasticStackIsConsistent()
  {
    // the tangential gradient is active: the surface gradients of the two faces differ
    const auto layers = card( 0, 2, { { 0.3, { 4000., 0.3 } }, { 0.7, { 20000., 0.2 } } } );
    checkTangent( "LINEARELASTIC", layers, makeKinematics( 1e-3, false ), 1e-6, "elastic stack" );
    checkTangent( "LINEARELASTIC", layers, makeKinematics( 1e-3, true ), 1e-6, "elastic stack, faces apart" );

    const auto graded = card( 1, 3, { { 0.0, { 4000., 0.3 } }, { 1.0, { 16000., 0.3 } } } );
    checkTangent( "LINEARELASTIC", graded, makeKinematics( 1e-3, true ), 1e-6, "continuously graded layer" );
  }

  void testTangentOfPlasticStackIsConsistent()
  {
    // the soft layer yields in the increment, the hard one does not
    const auto layers = card( 0,
                              2,
                              { { 0.4, { 1e4, 0.3, 40., 10., 0., 1., 2400. } },
                                { 0.6, { 2e4, 0.3, 4000., 10., 0., 1., 2400. } } } );
    // the tangent neglects the derivative of the tangential gradient terms, which is small for small differences of the
    // gradients
    Kinematics kinematics          = makeKinematics( 2e-3, false );
    kinematics.x.segment< 9 >( 3 ) = kinematics.x.segment< 9 >( 12 ) * 1.0001;
    checkTangent( "VONMISES", layers, kinematics, 1e-5, "plastic stack" );
  }

  /// Two elastic layers without lateral contraction: the tangential gradient changes the shear traction by a known
  /// amount.
  void testTangentialGradientOfLayeredStackAgainstAnalyticalValue()
  {
    const double E1 = 4000., E2 = 16000.;
    const auto   layers = card( 0, 1, { { 0.5, { E1, 0. } }, { 0.5, { E2, 0. } } } );
    GradedPoint  point( "LINEARELASTIC", layers );

    // normal jump j, and a variation of the jump along x: the surface gradient of the top is du_z/dx = g
    const double j = 1e-4 * h, g = 2e-3;

    Kinematics kinematics;
    kinematics.normal             = Eigen::Vector3d( 0., 0., 1. );
    kinematics.x( 2 )             = j;
    kinematics.x( 3 + 3 * 2 + 0 ) = g; // A+_zx

    const Response response = evaluate( *point.material, point.stateVars.data(), kinematics );

    const double Eh             = 1. / ( 0.5 / E1 + 0.5 / E2 );
    const double Gh             = 0.5 * Eh;
    const double c0             = Eh / E1 - 1.;
    const double expectedShear  = Gh * g * ( 0.5 + 0.25 * c0 );
    const double expectedNormal = Eh * j / h;

    throwExceptionOnFailure( std::abs( response.p( 2 ) - expectedNormal ) < 1e-9 * std::abs( expectedNormal ),
                             "The normal force of the layered stack is not the series value." );
    throwExceptionOnFailure( std::abs( response.p( 0 ) - expectedShear ) < 1e-9 * std::abs( expectedShear ),
                             "The shear force of the layered stack does not contain the tangential gradient of the "
                             "correction." );
  }

  /**
   * A stack of two elastic layers, whose faces are displaced such that the jump of the normal displacement varies
   * linearly along the interface: u_z = z0 at the top and 0 at the bottom, with du_z/dx = g at the top. The exact
   * solution of plane strain for this case is polynomial, linear in x, in every layer, with continuous
   * displacements and tractions at the interface of the layers. It was derived with the computer algebra system
   * sympy for the thickness 0.02, z0 = 3e-5, g = 2e-3, the lower layer E = 4000, nu = 0.3 of 40 % of the
   * thickness and the upper layer E = 20000, nu = 0.2. The mean of the stress and the moments h int y sigma dy and
   * h int (1-y) sigma dy at x = 0 are the generalized stresses of the material. The material reproduces them,
   * which includes the linear variation of the shear stress through the thickness of every layer, that is the
   * tangential divergence of the stress, and the effect of the tangential gradient of the correction.
   */
  void testLayeredStackAgainstExactPlaneStrainSolution()
  {
    const auto  layers = card( 0, 3, { { 0.4, { 4000., 0.3 } }, { 0.6, { 20000., 0.2 } } } );
    GradedPoint point( "LINEARELASTIC", layers );

    Kinematics kinematics;
    kinematics.normal             = Eigen::Vector3d( 0., 0., 1. );
    kinematics.x( 2 )             = 3e-5;
    kinematics.x( 3 + 3 * 2 + 0 ) = 2e-3; // A+_zx

    const Response response = evaluate( *point.material, point.stateVars.data(), kinematics );

    // row-major index of the entry (i, j) of the surface stress, which starts at 3 for S+ and at 12 for S-
    const auto S = [&]( int offset, int i, int j ) { return response.resultant( offset + 3 * i + j ); };

    struct Reference {
      const char* name;
      double      value;
      double      actual;
    };
    const Reference references[] = {
      { "force x", 2.8004826074396316, response.resultant( 0 ) },
      { "force z", 14.809590973201692, response.resultant( 2 ) },
      { "S+ xx", 0.041255289139633285, S( 3, 0, 0 ) },
      { "S+ yy", 0.041255289139633285, S( 3, 1, 1 ) },
      { "S+ zz", 0.14809590973201692, S( 3, 2, 2 ) },
      { "S+ xz", 0.01770863425493228, S( 3, 0, 2 ) },
      { "S- xx", 0.05394922425952045, S( 12, 0, 0 ) },
      { "S- zz", 0.14809590973201692, S( 12, 2, 2 ) },
      { "S- xz", 0.03830101789386035, S( 12, 0, 2 ) },
    };
    for ( const auto& reference : references )
      throwExceptionOnFailure( std::abs( reference.actual - reference.value ) < 1e-8 * std::abs( reference.value ),
                               std::string( "The stack differs from the exact plane strain solution in " ) +
                                 reference.name );
  }

  /**
   * The second bending mode of the same stack: the top face stretches along the interface, u_x = 2e-3 x at the top and
   * zero at the bottom, so that the axial strain of the layer varies linearly through the thickness. The exact
   * solution of plane strain was derived with sympy as above. The mean of the stress and the moments at x = 0 are
   * the generalized stresses.
   */
  void testLayeredStackInDifferentialStretchAgainstExactPlaneStrainSolution()
  {
    const auto  layers = card( 0, 3, { { 0.4, { 4000., 0.3 } }, { 0.6, { 20000., 0.2 } } } );
    GradedPoint point( "LINEARELASTIC", layers );

    Kinematics kinematics;
    kinematics.normal = Eigen::Vector3d( 0., 0., 1. );
    kinematics.x( 3 ) = 2e-3; // A+_xx

    const Response response = evaluate( *point.material, point.stateVars.data(), kinematics );

    const auto S = [&]( int offset, int i, int j ) { return response.resultant( offset + 3 * i + j ); };

    struct Reference {
      const char* name;
      double      value;
      double      actual;
    };
    const Reference references[] = {
      { "S+ xx", 0.33095616396676775, S( 3, 0, 0 ) },
      { "S- xx", 0.16398835317089433, S( 12, 0, 0 ) },
      { "S+ yy", 0.07183969810331392, S( 3, 1, 1 ) },
      { "S- yy", 0.047201204576516824, S( 12, 1, 1 ) },
      { "S+ zz", 0.019579021303833103, S( 3, 2, 2 ) },
      { "S- zz", 0.03965934258897367, S( 12, 2, 2 ) },
      { "force z", 2.9619181946403383, response.resultant( 2 ) },
    };
    for ( const auto& reference : references )
      throwExceptionOnFailure( std::abs( reference.actual - reference.value ) < 1e-8 * std::abs( reference.value ),
                               std::string( "The stack in differential stretch differs from the exact plane strain "
                                            "solution in " ) +
                                 reference.name );
  }

  /**
   * The generalized stress is the derivative of the energy of the layer: its test function is the derivative of the
   * strain of the points, including the response of the microproblem and the tangential gradient of the correction.
   * For an elastic stack the tangent is therefore symmetric, with and without the tangential gradient, and the
   * generalized stress is the tangent applied to the generalized strain.
   */
  void testGeneralizedStressIsTheDerivativeOfAnEnergy()
  {
    for ( int withTangentialGradient = 0; withTangentialGradient <= 1; withTangentialGradient++ ) {
      GradedPoint point( "LINEARELASTIC",
                         card( 0, 3, { { 0.4, { 4000., 0.3 } }, { 0.6, { 20000., 0.2 } } }, withTangentialGradient ) );

      const Kinematics kinematics = makeKinematics( 1e-3, true );
      const Response   response   = evaluate( *point.material, point.stateVars.data(), kinematics );

      const Eigen::Matrix< double, 21, 21 > asymmetry = response.K - response.K.transpose();
      throwExceptionOnFailure( asymmetry.norm() < 1e-12 * response.K.norm(),
                               "The tangent of the elastic stack is not symmetric, tangential gradient flag " +
                                 std::to_string( withTangentialGradient ) + "." );

      // the increment of the generalized strain has the jump scaled by the thickness in makeKinematics
      const Eigen::Matrix< double, 21, 1 > expected = response.K * kinematics.x;
      throwExceptionOnFailure( checkIfEqual< double >( response.p, expected, 1e-9 * ( 1. + expected.norm() ) ),
                               "The generalized stress of the elastic stack is not the tangent applied to the "
                               "strain, tangential gradient flag " +
                                 std::to_string( withTangentialGradient ) + "." );
    }
  }

  /// A homogeneous layer does not depend on its decomposition into layers.
  void testHomogeneousLayerDoesNotDependOnItsDecomposition()
  {
    const std::vector< double > bulk = { 5000., 0.3 };
    GradedPoint                 single( "LINEARELASTIC", card( 0, 3, { { 1.0, bulk } } ) );
    GradedPoint                 stack( "LINEARELASTIC", card( 0, 3, { { 0.2, bulk }, { 0.5, bulk }, { 0.3, bulk } } ) );

    const Kinematics kinematics = makeKinematics( 1e-3, true );
    const Response   a          = evaluate( *single.material, single.stateVars.data(), kinematics );
    const Response   b          = evaluate( *stack.material, stack.stateVars.data(), kinematics );

    throwExceptionOnFailure( checkIfEqual< double >( a.p, b.p, 1e-9 * ( 1. + a.p.norm() ) ),
                             "A homogeneous layer depends on its decomposition into layers." );
  }

  void testInvalidPropertiesAreRefused()
  {
    bool thrown = false;
    try {
      const auto properties = card( 0, 1, { { 1.0, { 5000., 0.3 } } } );
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
  auto tests = std::vector<
    std::function< void() > >{ testElasticTwoFacesAreEquilibratedMaterial,
                               testPlasticTwoFacesAreEquilibratedMaterial,
                               testTangentOfElasticStackIsConsistent,
                               testTangentOfPlasticStackIsConsistent,
                               testTangentialGradientOfLayeredStackAgainstAnalyticalValue,
                               testLayeredStackAgainstExactPlaneStrainSolution,
                               testLayeredStackInDifferentialStretchAgainstExactPlaneStrainSolution,
                               testGeneralizedStressIsTheDerivativeOfAnEnergy,
                               testHomogeneousLayerDoesNotDependOnItsDecomposition,
                               testInvalidPropertiesAreRefused };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
