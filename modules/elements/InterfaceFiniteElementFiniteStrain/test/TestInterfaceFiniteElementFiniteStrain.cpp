#include "Marmot/GYInterfaceFiniteElementFiniteStrain.h"
#include "Marmot/InterfaceFiniteElementFiniteStrain.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Elements;

namespace {

  constexpr double pi = 3.14159265358979323846;

  template < int nDim, int nNodes >
  using DofVector = Eigen::Matrix< double, nDim * nNodes, 1 >;

  template < int nDim, int nNodes >
  using Coordinates = std::array< double, nDim * nNodes >;

  /// the interface material properties: the thickness, the bulk modulus and the shear modulus
  std::vector< double > neoHookeProperties( double thickness )
  {
    return { thickness, 3500.0, 1500.0 };
  }

  /// the properties of the equilibrated interface material: the bulk modulus, the shear modulus and the thickness
  std::vector< double > equilibratedNeoHookeProperties( double thickness )
  {
    return { 3500.0, 1500.0, thickness };
  }

  /// A bilinear interface whose bottom side is inclined by the angle about the y axis, and whose top side is
  /// offset VERTICALLY by the gap, so that the connector is not parallel to the normal.
  Coordinates< 3, 8 > inclinedQuadrilateral( double angle, double gap )
  {
    const double x[4] = { -0.5, 0.5, 0.5, -0.5 };
    const double y[4] = { -0.5, -0.5, 0.5, 0.5 };

    Coordinates< 3, 8 > coordinates{};
    for ( int side = 0; side < 2; ++side )
      for ( int node = 0; node < 4; ++node ) {
        const int n            = 4 * side + node;
        coordinates[3 * n + 0] = x[node];
        coordinates[3 * n + 1] = y[node];
        coordinates[3 * n + 2] = std::tan( angle ) * x[node] + side * gap;
      }
    return coordinates;
  }

  /// A serendipity interface: corners first, then the midsides, bottom side first, with the top offset by the gap.
  Coordinates< 3, 16 > serendipityQuadrilateral( double gap )
  {
    const double xi[8][2] = { { -0.5, -0.5 },
                              { 0.5, -0.5 },
                              { 0.5, 0.5 },
                              { -0.5, 0.5 },
                              { 0.0, -0.5 },
                              { 0.5, 0.0 },
                              { 0.0, 0.5 },
                              { -0.5, 0.0 } };

    Coordinates< 3, 16 > coordinates{};
    for ( int side = 0; side < 2; ++side )
      for ( int node = 0; node < 8; ++node ) {
        const int n            = 8 * side + node;
        coordinates[3 * n + 0] = xi[node][0];
        coordinates[3 * n + 1] = xi[node][1];
        coordinates[3 * n + 2] = 0.03 * xi[node][0] * xi[node][1] + side * gap;
      }
    return coordinates;
  }

  /// The inclined line of the interface studies: inclined by the angle, with the top offset vertically by the gap.
  Coordinates< 2, 4 > inclinedLine( double angle, double gap )
  {
    return { 0.0, 0.0, std::cos( angle ), std::sin( angle ), 0.0, gap, std::cos( angle ), std::sin( angle ) + gap };
  }

  struct Evaluation {
    Eigen::VectorXd residual;
    Eigen::MatrixXd tangent;
  };

  /// Evaluate the element for the total displacements, starting from the initial state.
  template < class ElementType, int nDim, int nNodes >
  Evaluation evaluate( const Coordinates< nDim, nNodes >& coordinates,
                       const std::vector< double >&       materialProperties,
                       const Eigen::VectorXd&             Q,
                       const std::string&                 materialName = "COMPRESSIBLENEOHOOKE" )
  {
    constexpr int nDofs = ElementType::sizeLoadVector;

    ElementType element( 3,
                         FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                         ElementType::SectionType::Interface );
    element.assignNodeCoordinates( coordinates.data() );

    const std::array< double, 1 > thickness = { 1.0 };
    ElementProperties             properties( thickness.data(), 1 );
    element.assignProperty( properties );
    element.assignMaterial( materialName, materialProperties.data(), static_cast< int >( materialProperties.size() ) );

    std::vector< double > stateVars( element.getNumberOfRequiredStateVars(), 0.0 );
    element.assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element.initializeYourself();
    element.setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

    std::array< double, nDofs >         Qarray{};
    std::array< double, nDofs >         Pe{};
    std::array< double, nDofs * nDofs > Ke{};
    std::copy( Q.data(), Q.data() + nDofs, Qarray.begin() );

    // the material has no history here, so the increment is irrelevant to it
    element.computeKernels( Qarray.data(), Qarray.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

    Evaluation evaluation;
    evaluation.residual = Eigen::Map< Eigen::VectorXd >( Pe.data(), nDofs );
    // the element tangent is stored column major, like the one of the element base
    evaluation.tangent = Eigen::Map< Eigen::Matrix< double, nDofs, nDofs > >( Ke.data() );
    return evaluation;
  }

  template < int nDofs >
  Eigen::Matrix< double, nDofs, 1 > makeDisplacements( double scale, int pattern )
  {
    Eigen::Matrix< double, nDofs, 1 > Q;
    for ( int i = 0; i < nDofs; ++i )
      Q( i ) = scale * ( ( ( i * 7 + pattern * 5 ) % 13 ) - 6 ) / 6.0;
    return Q;
  }

  /// The tangent has to be the derivative of the residual with respect to the total displacements.
  template < class ElementType, int nDim, int nNodes >
  void assertTangentMatchesFiniteDifference( const Coordinates< nDim, nNodes >& coordinates,
                                             const std::vector< double >&       materialProperties,
                                             const Eigen::VectorXd&             Q,
                                             const std::string&                 context,
                                             const std::string&                 materialName = "COMPRESSIBLENEOHOOKE" )
  {
    constexpr int nDofs = ElementType::sizeLoadVector;

    const auto base = evaluate< ElementType, nDim, nNodes >( coordinates, materialProperties, Q, materialName );
    throwExceptionOnFailure( base.residual.allFinite() && base.tangent.allFinite(), context + ": nan or inf." );
    throwExceptionOnFailure( base.residual.norm() > 0.0, context + ": the displacements must produce a residual." );

    Eigen::MatrixXd finiteDifference( nDofs, nDofs );
    for ( int i = 0; i < nDofs; ++i ) {
      const double eps = 1e-6;

      Eigen::VectorXd plus  = Q;
      Eigen::VectorXd minus = Q;
      plus( i ) += eps;
      minus( i ) -= eps;

      finiteDifference.col(
        i ) = ( evaluate< ElementType, nDim, nNodes >( coordinates, materialProperties, plus, materialName ).residual -
                evaluate< ElementType, nDim, nNodes >( coordinates, materialProperties, minus, materialName )
                  .residual ) /
              ( 2.0 * eps );
    }

    const double error = ( finiteDifference - base.tangent ).norm() / std::max( 1.0, base.tangent.norm() );
    throwExceptionOnFailure( error < 1e-6,
                             context +
                               ": the tangent does not match the finite difference of the residual, "
                               "relative error " +
                               std::to_string( error ) );
  }

  void testTangentMatchesFiniteDifference()
  {
    const double angle = 10.0 * pi / 180.0;
    const double gap   = 0.1;

    assertTangentMatchesFiniteDifference< InterfaceFiniteElementFiniteStrain< 3, 8 >,
                                          3,
                                          8 >( inclinedQuadrilateral( angle, gap ),
                                               neoHookeProperties( gap * std::cos( angle ) ),
                                               makeDisplacements< 24 >( 3e-3, 1 ),
                                               "3D, inclined" );

    assertTangentMatchesFiniteDifference< InterfaceFiniteElementFiniteStrain< 3, 16 >,
                                          3,
                                          16 >( serendipityQuadrilateral( gap ),
                                                neoHookeProperties( gap ),
                                                makeDisplacements< 48 >( 3e-3, 2 ),
                                                "3D, serendipity" );

    assertTangentMatchesFiniteDifference< InterfaceFiniteElementFiniteStrain< 2, 4 >,
                                          2,
                                          4 >( inclinedLine( angle, gap ),
                                               neoHookeProperties( gap * std::cos( angle ) ),
                                               makeDisplacements< 8 >( 5e-3, 3 ),
                                               "2D, inclined" );
  }

  /**
   * A rigid rotation of the whole element must not cause stresses. This is the property of the connector-based
   * reduction: on an inclined interface the connector is not parallel to the normal, and the jump alone would make
   * the layer stressed under a rotation.
   */
  void testRigidRotationIsStressFree()
  {
    const double angle = 10.0 * pi / 180.0;
    const double gap   = 0.1;
    const double theta = 0.3;

    const auto coordinates = inclinedQuadrilateral( angle, gap );

    Eigen::Matrix3d R;
    R << std::cos( theta ), 0.0, std::sin( theta ), 0.0, 1.0, 0.0, -std::sin( theta ), 0.0, std::cos( theta );

    DofVector< 3, 8 > Q;
    for ( int node = 0; node < 8; ++node ) {
      const Eigen::Vector3d X( coordinates[3 * node], coordinates[3 * node + 1], coordinates[3 * node + 2] );
      Q.segment< 3 >( 3 * node ) = ( R - Eigen::Matrix3d::Identity() ) * X;
    }

    const auto rotated = evaluate< InterfaceFiniteElementFiniteStrain< 3, 8 >, 3, 8 >( coordinates,
                                                                                       neoHookeProperties(
                                                                                         gap * std::cos( angle ) ),
                                                                                       Q );

    // the scale of the forces of a deformation of the order of the rotation, for comparison
    const auto stretched = evaluate< InterfaceFiniteElementFiniteStrain< 3, 8 >, 3, 8 >( coordinates,
                                                                                         neoHookeProperties(
                                                                                           gap * std::cos( angle ) ),
                                                                                         makeDisplacements< 24 >( 1e-2,
                                                                                                                  1 ) );

    throwExceptionOnFailure( rotated.residual.norm() < 1e-9 * stretched.residual.norm(),
                             "A rigid rotation of an inclined interface must be stress free, residual " +
                               std::to_string( rotated.residual.norm() ) );
  }

  void testRigidTranslationIsStressFree()
  {
    const auto coordinates = serendipityQuadrilateral( 0.1 );

    DofVector< 3, 16 > Q;
    for ( int node = 0; node < 16; ++node )
      Q.segment< 3 >( 3 * node ) = Eigen::Vector3d( 0.3, -0.2, 0.1 );

    const auto translated = evaluate< InterfaceFiniteElementFiniteStrain< 3, 16 >, 3, 16 >( coordinates,
                                                                                            neoHookeProperties( 0.1 ),
                                                                                            Q );
    throwExceptionOnFailure( translated.residual.norm() < 1e-9, "A rigid translation must not cause forces." );
  }

  void testCoincidentFacesAreRejected()
  {
    const auto coordinates = inclinedQuadrilateral( 0.0, 0.0 );

    bool rejected = false;
    try {
      evaluate< InterfaceFiniteElementFiniteStrain< 3, 8 >, 3, 8 >( coordinates,
                                                                    neoHookeProperties( 0.1 ),
                                                                    makeDisplacements< 24 >( 1e-3, 1 ) );
    }
    catch ( const std::invalid_argument& ) {
      rejected = true;
    }
    throwExceptionOnFailure( rejected, "Coincident faces carry no normal fibre and must be rejected." );
  }

  void testInvalidMaterialsAreRejected()
  {
    const auto build = []( const std::string& name, const std::vector< double >& properties ) {
      MarmotInterfaceMaterialFiniteStrain material( name,
                                                    properties.data(),
                                                    static_cast< int >( properties.size() ),
                                                    1 );
    };

    const auto rejects = [&]( const std::string& name, const std::vector< double >& properties ) {
      try {
        build( name, properties );
      }
      catch ( const std::invalid_argument& ) {
        return true;
      }
      return false;
    };

    throwExceptionOnFailure( rejects( "COMPRESSIBLENEOHOOKE", { 0.0, 3500.0, 1500.0 } ),
                             "A zero thickness is rejected." );
    throwExceptionOnFailure( rejects( "COMPRESSIBLENEOHOOKE", { -0.1, 3500.0, 1500.0 } ),
                             "A negative thickness is rejected." );
    throwExceptionOnFailure( rejects( "COMPRESSIBLENEOHOOKE", {} ), "Missing properties are rejected." );
    throwExceptionOnFailure( rejects( "NOSUCHMATERIAL", { 0.1, 3500.0, 1500.0 } ),
                             "An unknown base material is rejected." );
  }

  /// The equilibrated element: tangent and objectivity.
  void testEquilibratedElement()
  {
    const double angle = 10.0 * pi / 180.0;
    const double gap   = 0.1;
    const double hn    = gap * std::cos( angle );

    assertTangentMatchesFiniteDifference< YInterfaceFiniteElementFiniteStrain< 3, 8 >,
                                          3,
                                          8 >( inclinedQuadrilateral( angle, gap ),
                                               equilibratedNeoHookeProperties( hn ),
                                               makeDisplacements< 24 >( 3e-3, 1 ),
                                               "Y 3D, inclined" );

    assertTangentMatchesFiniteDifference< YInterfaceFiniteElementFiniteStrain< 3, 16 >,
                                          3,
                                          16 >( serendipityQuadrilateral( gap ),
                                                equilibratedNeoHookeProperties( gap ),
                                                makeDisplacements< 48 >( 3e-3, 2 ),
                                                "Y 3D, serendipity" );

    assertTangentMatchesFiniteDifference< YInterfaceFiniteElementFiniteStrain< 2, 4 >,
                                          2,
                                          4 >( inclinedLine( angle, gap ),
                                               equilibratedNeoHookeProperties( hn ),
                                               makeDisplacements< 8 >( 5e-3, 3 ),
                                               "Y 2D, inclined" );

    // a rigid rotation is stress free
    const double theta       = 0.3;
    const auto   coordinates = inclinedQuadrilateral( angle, gap );

    Eigen::Matrix3d R;
    R << std::cos( theta ), 0.0, std::sin( theta ), 0.0, 1.0, 0.0, -std::sin( theta ), 0.0, std::cos( theta );

    DofVector< 3, 8 > Q;
    for ( int node = 0; node < 8; ++node ) {
      const Eigen::Vector3d X( coordinates[3 * node], coordinates[3 * node + 1], coordinates[3 * node + 2] );
      Q.segment< 3 >( 3 * node ) = ( R - Eigen::Matrix3d::Identity() ) * X;
    }

    const auto rotated = evaluate< YInterfaceFiniteElementFiniteStrain< 3, 8 >, 3, 8 >( coordinates,
                                                                                        equilibratedNeoHookeProperties(
                                                                                          hn ),
                                                                                        Q );
    const auto
      stretched = evaluate< YInterfaceFiniteElementFiniteStrain< 3, 8 >, 3, 8 >( coordinates,
                                                                                 equilibratedNeoHookeProperties( hn ),
                                                                                 makeDisplacements< 24 >( 1e-2, 1 ) );

    throwExceptionOnFailure( rotated.residual.norm() < 1e-9 * stretched.residual.norm(),
                             "A rigid rotation must be stress free, residual " +
                               std::to_string( rotated.residual.norm() ) );
  }

  /// With a uniform shift of the top face the two faces carry the same stress, and the equilibrated element gives the
  /// force of the plain element.
  void testEquilibratedElementAgreesWithThePlainElementForAUniformShift()
  {
    const auto coordinates = inclinedQuadrilateral( 0.0, 0.1 );

    DofVector< 3, 8 > shift = DofVector< 3, 8 >::Zero();
    for ( int node = 4; node < 8; ++node )
      shift.segment< 3 >( 3 * node ) = Eigen::Vector3d( 1.5e-3, -0.8e-3, 2.0e-3 );

    const auto plain = evaluate< InterfaceFiniteElementFiniteStrain< 3, 8 >, 3, 8 >( coordinates,
                                                                                     neoHookeProperties( 0.1 ),
                                                                                     shift );
    const auto
      equilibrated = evaluate< YInterfaceFiniteElementFiniteStrain< 3, 8 >, 3, 8 >( coordinates,
                                                                                    equilibratedNeoHookeProperties(
                                                                                      0.1 ),
                                                                                    shift );

    throwExceptionOnFailure( ( plain.residual - equilibrated.residual ).norm() < 1e-9 * plain.residual.norm(),
                             "The equilibrated element must give the residual of the plain element for a uniform "
                             "shift." );
  }

  // bulk modulus, shear modulus, thickness, damage threshold, softening parameter, nonlocal radius, density
  std::vector< double > gradientEnhancedProperties( double thickness )
  {
    return { 3500.0, 1500.0, thickness, 1.0e-3, 5.0e-3, 0.2, 1.0 };
  }

  const std::string gradientEnhancedMaterial = "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE";

  /// The displacements and the nonlocal field of an element with a nonlocal field.
  template < class ElementType >
  Eigen::VectorXd makeDisplacementsAndNonlocalField( double scale, int pattern, double nonlocal )
  {
    Eigen::VectorXd Q = makeDisplacements< ElementType::sizeLoadVector >( scale, pattern );
    Q.tail( ElementType::nExtraDofs ).setConstant( nonlocal );
    for ( int i = 0; i < ElementType::nExtraDofs; ++i )
      Q( ElementType::nDofDisplacement + i ) += 1.0e-4 * i;
    return Q;
  }

  /// The element with the nonlocal field: the tangent over the displacements and the nonlocal field, with the damage
  /// active.
  void testGradientEnhancedElement()
  {
    const double angle = 10.0 * pi / 180.0;
    const double gap   = 0.1;
    const double hn    = gap * std::cos( angle );

    using Element3d8  = GYInterfaceFiniteElementFiniteStrain< 3, 8 >;
    using Element3d16 = GYInterfaceFiniteElementFiniteStrain< 3, 16 >;
    using Element2d   = GYInterfaceFiniteElementFiniteStrain< 2, 4 >;

    assertTangentMatchesFiniteDifference< Element3d8, 3, 8 >( inclinedQuadrilateral( angle, gap ),
                                                              gradientEnhancedProperties( hn ),
                                                              makeDisplacementsAndNonlocalField< Element3d8 >( 3e-3,
                                                                                                               1,
                                                                                                               3e-3 ),
                                                              "GY 3D, inclined",
                                                              gradientEnhancedMaterial );

    assertTangentMatchesFiniteDifference< Element3d16, 3, 16 >( serendipityQuadrilateral( gap ),
                                                                gradientEnhancedProperties( gap ),
                                                                makeDisplacementsAndNonlocalField<
                                                                  Element3d16 >( 3e-3, 2, 3e-3 ),
                                                                "GY 3D, serendipity",
                                                                gradientEnhancedMaterial );

    assertTangentMatchesFiniteDifference< Element2d, 2, 4 >( inclinedLine( angle, gap ),
                                                             gradientEnhancedProperties( hn ),
                                                             makeDisplacementsAndNonlocalField< Element2d >( 5e-3,
                                                                                                             3,
                                                                                                             3e-3 ),
                                                             "GY 2D, inclined",
                                                             gradientEnhancedMaterial );
  }

  /// A rigid rotation of the element with no nonlocal field is stress free and has no local field to drive it.
  void testGradientEnhancedElementIsObjective()
  {
    const double angle       = 10.0 * pi / 180.0;
    const double gap         = 0.1;
    const double theta       = 0.3;
    const auto   coordinates = inclinedQuadrilateral( angle, gap );

    using Element = GYInterfaceFiniteElementFiniteStrain< 3, 8 >;

    Eigen::Matrix3d R;
    R << std::cos( theta ), 0.0, std::sin( theta ), 0.0, 1.0, 0.0, -std::sin( theta ), 0.0, std::cos( theta );

    Eigen::VectorXd Q = Eigen::VectorXd::Zero( Element::sizeLoadVector );
    for ( int node = 0; node < 8; ++node ) {
      const Eigen::Vector3d X( coordinates[3 * node], coordinates[3 * node + 1], coordinates[3 * node + 2] );
      Q.segment< 3 >( 3 * node ) = ( R - Eigen::Matrix3d::Identity() ) * X;
    }

    const auto rotated   = evaluate< Element, 3, 8 >( coordinates,
                                                    gradientEnhancedProperties( gap * std::cos( angle ) ),
                                                    Q,
                                                    gradientEnhancedMaterial );
    const auto stretched = evaluate< Element, 3, 8 >( coordinates,
                                                      gradientEnhancedProperties( gap * std::cos( angle ) ),
                                                      makeDisplacementsAndNonlocalField< Element >( 1e-2, 1, 0.0 ),
                                                      gradientEnhancedMaterial );

    throwExceptionOnFailure( rotated.residual.norm() < 1e-9 * stretched.residual.norm(),
                             "A rigid rotation must be stress free, residual " +
                               std::to_string( rotated.residual.norm() ) );
  }

  /// The degrees of freedom of the nonlocal field sit on the nodes of the bottom side.
  void testGradientEnhancedElementDofs()
  {
    using Element = GYInterfaceFiniteElementFiniteStrain< 3, 8 >;

    Element element( 1, FiniteElement::Quadrature::IntegrationTypes::FullIntegration, Element::SectionType::Interface );

    throwExceptionOnFailure( element.getNDofPerElement() == 8 * 3 + 4,
                             "The element has the displacements and four nonlocal dofs." );

    const auto fields = element.getNodeFields();
    throwExceptionOnFailure( fields.size() == 8, "One list of fields per node." );
    for ( size_t node = 0; node < fields.size(); ++node )
      throwExceptionOnFailure( fields[node].size() == ( node < 4 ? 2u : 1u ),
                               "Only the bottom nodes carry the nonlocal field." );
    throwExceptionOnFailure( fields[0][1] == "nonlocal damage",
                             "The nonlocal field is named like the one of the bulk element." );

    // node by node: the four bottom nodes with [u, abar], then the four top nodes with [u]
    auto permutation = element.getDofIndicesPermutationPattern();
    throwExceptionOnFailure( permutation.size() == 28, "A place for every dof." );
    throwExceptionOnFailure( permutation[0] == 0 && permutation[2] == 2 && permutation[3] == 4 &&
                               permutation[12] == 16 && permutation[24] == 3 && permutation[27] == 15,
                             "The permutation places the degrees of freedom node by node." );
    std::sort( permutation.begin(), permutation.end() );
    for ( int i = 0; i < 28; ++i )
      throwExceptionOnFailure( permutation[i] == i, "The permutation is a bijection." );
  }

  /// A state of uniform strain and uniform nonlocal field balances in the nonlocal equation.
  void testGradientEnhancedElementNamedState()
  {
    using Element = GYInterfaceFiniteElementFiniteStrain< 3, 8 >;

    const auto coordinates = inclinedQuadrilateral( 0.0, 0.1 );

    Element element( 1, FiniteElement::Quadrature::IntegrationTypes::FullIntegration, Element::SectionType::Interface );
    element.assignNodeCoordinates( coordinates.data() );
    const auto properties = gradientEnhancedProperties( 0.1 );
    element.assignMaterial( gradientEnhancedMaterial, properties.data(), static_cast< int >( properties.size() ) );

    std::vector< double > stateVars( element.getNumberOfRequiredStateVars(), 0.0 );
    element.assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element.initializeYourself();
    element.setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

    throwExceptionOnFailure( element.getStateView( "localField", 0 ).stateSize == 1,
                             "The local field is a named state." );
    throwExceptionOnFailure( element.getStateView( "nonlocalField", 0 ).stateSize == 1,
                             "The nonlocal field is a named state." );
    throwExceptionOnFailure( element.getStateView( "topMaterialStateVars", 0 ).stateSize > 0,
                             "The state of the material is reachable." );

    // A uniform shift of the top side and a uniform nonlocal field: the nonlocal residual is the integral of
    // h ( abar - Lbar ) over the midsurface, which has the area 1 here, and the gradient term does not contribute.
    const double nonlocal = 2.0e-3;

    Eigen::VectorXd Q = Eigen::VectorXd::Zero( Element::sizeLoadVector );
    for ( int node = 4; node < 8; ++node )
      Q.segment< 3 >( 3 * node ) = Eigen::Vector3d( 1.5e-3, -0.8e-3, 2.0e-3 );
    Q.tail( Element::nExtraDofs ).setConstant( nonlocal );

    std::vector< double > Pe( Element::sizeLoadVector, 0.0 ),
      Ke( Element::sizeLoadVector * Element::sizeLoadVector, 0.0 );
    element.computeKernels( Q.data(), Q.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

    const double localField = element.getStateView( "localField", 0 ).stateLocation[0];
    throwExceptionOnFailure( localField > 0.0, "The shift strains the layer, which gives a local field." );

    double nonlocalResidual = 0.0;
    for ( int i = 0; i < Element::nExtraDofs; ++i )
      nonlocalResidual += Pe[Element::nDofDisplacement + i];
    throwExceptionOnFailure( std::abs( nonlocalResidual - 0.1 * ( nonlocal - localField ) ) < 1e-12,
                             "The nonlocal residual is not the integral of h ( abar - Lbar )." );
  }

  void testRegistration()
  {
    for ( const std::string name : { "ILINE2UL",
                                     "IQUAD4UL",
                                     "IQUAD8UL",
                                     "YILINE2UL",
                                     "YIQUAD4UL",
                                     "YIQUAD8UL",
                                     "GYILINE2UL",
                                     "GYIQUAD4UL",
                                     "GYIQUAD8UL" } )
      throwExceptionOnFailure( std::unique_ptr< MarmotElement >(
                                 MarmotLibrary::MarmotElementFactory::createElement( name, 1 ) ) != nullptr,
                               name + " is not registered with the element factory." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testTangentMatchesFiniteDifference,
                                                   testRigidRotationIsStressFree,
                                                   testGradientEnhancedElement,
                                                   testGradientEnhancedElementIsObjective,
                                                   testGradientEnhancedElementDofs,
                                                   testGradientEnhancedElementNamedState,
                                                   testEquilibratedElement,
                                                   testEquilibratedElementAgreesWithThePlainElementForAUniformShift,
                                                   testRigidTranslationIsStressFree,
                                                   testCoincidentFacesAreRejected,
                                                   testInvalidMaterialsAreRejected,
                                                   testRegistration };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
