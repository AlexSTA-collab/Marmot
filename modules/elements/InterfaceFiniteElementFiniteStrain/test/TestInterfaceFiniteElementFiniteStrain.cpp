#include "Marmot/GYInterfaceFiniteElementFiniteStrain.h"
#include "Marmot/GYSInterfaceFiniteElementFiniteStrain.h"
#include "Marmot/InterfaceFiniteElementFiniteStrain.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
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

  /// the card of the interface element with the Dufour model: K, G, the thickness, and the other properties of the
  /// model
  std::vector< double > dufourProperties( double thickness )
  {
    return { 1793.9, 828.0, thickness, 12.0, 21.6, 4051.0, 15.8, 552.0, 243.0, 0.0, 0.0,
             15.4,   1e-6,  0.0435,    0.3,  0.5,  1.75,   0.99, 0.12,  1.0,   1.0 };
  }

  /// The element with the viscoplastic damage model: the return map and the damage are active, and the tangent over the
  /// displacements and the nonlocal field has to match the finite differences of the residual.
  void testGradientEnhancedElementWithTheDufourModel()
  {
    using Element = GYInterfaceFiniteElementFiniteStrain< 3, 8 >;

    const double angle = 10.0 * pi / 180.0;
    const double gap   = 0.1;

    const Eigen::VectorXd Q = makeDisplacementsAndNonlocalField< Element >( 6e-3, 1, 0.05 );

    // the state must be plastic, otherwise the test would not cover the return map
    const auto coordinates = inclinedQuadrilateral( angle, gap );
    {
      Element element( 1,
                       FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                       Element::SectionType::Interface );
      element.assignNodeCoordinates( coordinates.data() );
      const auto properties = dufourProperties( gap * std::cos( angle ) );
      element.assignMaterial( "DUFOURMODEL", properties.data(), static_cast< int >( properties.size() ) );

      std::vector< double > stateVars( element.getNumberOfRequiredStateVars(), 0.0 );
      element.assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
      element.initializeYourself();
      element.setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

      std::vector< double > Pe( Element::sizeLoadVector, 0.0 ),
        Ke( Element::sizeLoadVector * Element::sizeLoadVector, 0.0 );
      element.computeKernels( Q.data(), Q.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

      const double* topState = element.getStateView( "topMaterialStateVars", 0 ).stateLocation;
      // alphaP is the tenth entry of the state of the model
      throwExceptionOnFailure( topState[9] > 0.0, "The Dufour model should be plastic in this test." );
    }

    assertTangentMatchesFiniteDifference< Element, 3, 8 >( coordinates,
                                                           dufourProperties( gap * std::cos( angle ) ),
                                                           Q,
                                                           "GY with the Dufour model",
                                                           "DUFOURMODEL" );
  }

  void testMaterialAccessors()
  {
    // thickness, bulk modulus, shear modulus and density
    const std::vector< double >         properties = { 0.1, 3500.0, 1500.0, 2.0 };
    MarmotInterfaceMaterialFiniteStrain material( "COMPRESSIBLENEOHOOKE", properties.data(), 4, 1 );
    material.setCharacteristicElementLength( 0.5 );

    std::vector< double > stateVars( material.getNumberOfRequiredStateVars(), 0.0 );
    material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );

    throwExceptionOnFailure( material.getStateView( "baseMaterialStateVars", stateVars.data() ).stateSize ==
                               material.getNumberOfRequiredStateVars(),
                             "The state of the bulk material is a named state." );

    const double                                  bulkProperties[] = { 3500.0, 1500.0, 2.0 };
    std::unique_ptr< MarmotMaterialFiniteStrain > bulk(
      MarmotLibrary::MarmotMaterialFiniteStrainFactory::createMaterial( "COMPRESSIBLENEOHOOKE",
                                                                        bulkProperties,
                                                                        3,
                                                                        1 ) );
    throwExceptionOnFailure( material.getDensity() == bulk->getDensity( nullptr ),
                             "The density is the one of the bulk material." );
  }

  /// The displacements, and a nonlocal field for each face: the bottom nodes carry the field of the bottom face.
  template < class ElementType >
  Eigen::VectorXd makeDisplacementsAndFaceFields( double scale, int pattern, double nonlocalBottom, double nonlocalTop )
  {
    constexpr int nPerFace = ElementType::nExtraDofs / 2;

    Eigen::VectorXd Q = makeDisplacements< ElementType::sizeLoadVector >( scale, pattern );
    for ( int i = 0; i < nPerFace; ++i ) {
      Q( ElementType::nDofDisplacement + i )            = nonlocalBottom + 1.0e-4 * i;
      Q( ElementType::nDofDisplacement + nPerFace + i ) = nonlocalTop + 1.0e-4 * i;
    }
    return Q;
  }

  /// The element with a nonlocal field for each face: the tangent over the displacements and the fields of both faces,
  /// with different fields on the faces and the damage active.
  void testFaceFieldsElement()
  {
    const double angle = 10.0 * pi / 180.0;
    const double gap   = 0.1;
    const double hn    = gap * std::cos( angle );

    using Element3d8  = GY2InterfaceFiniteElementFiniteStrain< 3, 8 >;
    using Element3d16 = GY2InterfaceFiniteElementFiniteStrain< 3, 16 >;
    using Element2d   = GY2InterfaceFiniteElementFiniteStrain< 2, 4 >;

    assertTangentMatchesFiniteDifference< Element3d8, 3, 8 >( inclinedQuadrilateral( angle, gap ),
                                                              gradientEnhancedProperties( hn ),
                                                              makeDisplacementsAndFaceFields< Element3d8 >( 3e-3,
                                                                                                            1,
                                                                                                            3e-3,
                                                                                                            5e-3 ),
                                                              "GY2 3D, inclined",
                                                              gradientEnhancedMaterial );

    assertTangentMatchesFiniteDifference< Element3d16, 3, 16 >( serendipityQuadrilateral( gap ),
                                                                gradientEnhancedProperties( gap ),
                                                                makeDisplacementsAndFaceFields< Element3d16 >( 3e-3,
                                                                                                               2,
                                                                                                               4e-3,
                                                                                                               2e-3 ),
                                                                "GY2 3D, serendipity",
                                                                gradientEnhancedMaterial );

    assertTangentMatchesFiniteDifference< Element2d, 2, 4 >( inclinedLine( angle, gap ),
                                                             gradientEnhancedProperties( hn ),
                                                             makeDisplacementsAndFaceFields< Element2d >( 5e-3,
                                                                                                          3,
                                                                                                          3e-3,
                                                                                                          6e-3 ),
                                                             "GY2 2D, inclined",
                                                             gradientEnhancedMaterial );
  }

  /// A rigid rotation without nonlocal fields is stress free.
  void testFaceFieldsElementIsObjective()
  {
    const double angle       = 10.0 * pi / 180.0;
    const double gap         = 0.1;
    const double theta       = 0.3;
    const auto   coordinates = inclinedQuadrilateral( angle, gap );

    using Element = GY2InterfaceFiniteElementFiniteStrain< 3, 8 >;

    Eigen::Matrix3d R;
    R << std::cos( theta ), 0.0, std::sin( theta ), 0.0, 1.0, 0.0, -std::sin( theta ), 0.0, std::cos( theta );

    Eigen::VectorXd Q = Eigen::VectorXd::Zero( Element::sizeLoadVector );
    for ( int node = 0; node < 8; ++node ) {
      const Eigen::Vector3d X( coordinates[3 * node], coordinates[3 * node + 1], coordinates[3 * node + 2] );
      Q.segment< 3 >( 3 * node ) = ( R - Eigen::Matrix3d::Identity() ) * X;
    }

    const auto properties = gradientEnhancedProperties( gap * std::cos( angle ) );
    const auto rotated    = evaluate< Element, 3, 8 >( coordinates, properties, Q, gradientEnhancedMaterial );
    const auto stretched  = evaluate< Element, 3, 8 >( coordinates,
                                                      properties,
                                                      makeDisplacementsAndFaceFields< Element >( 1e-2, 1, 0.0, 0.0 ),
                                                      gradientEnhancedMaterial );

    throwExceptionOnFailure( rotated.residual.norm() < 1e-9 * stretched.residual.norm(),
                             "A rigid rotation must be stress free, residual " +
                               std::to_string( rotated.residual.norm() ) );
  }

  /// A nonlocal field on every node, the faces damage independently, and they are coupled through the thickness.
  void testFaceFieldsAreIndependent()
  {
    using Element = GY2InterfaceFiniteElementFiniteStrain< 3, 8 >;

    Element element( 1, FiniteElement::Quadrature::IntegrationTypes::FullIntegration, Element::SectionType::Interface );

    throwExceptionOnFailure( element.getNDofPerElement() == 8 * 3 + 8,
                             "The displacements and a nonlocal field on every node." );

    const auto fields = element.getNodeFields();
    throwExceptionOnFailure( fields.size() == 8, "One list of fields per node." );
    for ( const auto& nodeFields : fields )
      throwExceptionOnFailure( nodeFields.size() == 2 && nodeFields[1] == "nonlocal damage",
                               "Every node carries the field." );

    // node by node: every node has [u, N], the bottom nodes before the top nodes
    auto permutation = element.getDofIndicesPermutationPattern();
    throwExceptionOnFailure( permutation.size() == 32, "A place for every dof." );
    throwExceptionOnFailure( permutation[0] == 0 && permutation[3] == 4 && permutation[12] == 16 &&
                               permutation[24] == 3 && permutation[28] == 19 && permutation[31] == 31,
                             "The permutation places the degrees of freedom node by node." );
    std::sort( permutation.begin(), permutation.end() );
    for ( int i = 0; i < 32; ++i )
      throwExceptionOnFailure( permutation[i] == i, "The permutation is a bijection." );

    // The field of the top face stays below the damage threshold kappa0 = 1e-3 while the one of the bottom face is
    // above it: only the bottom face is damaged.
    const auto coordinates = inclinedQuadrilateral( 0.0, 0.1 );
    element.assignNodeCoordinates( coordinates.data() );
    const auto properties = gradientEnhancedProperties( 0.1 );
    element.assignMaterial( gradientEnhancedMaterial, properties.data(), static_cast< int >( properties.size() ) );

    std::vector< double > stateVars( element.getNumberOfRequiredStateVars(), 0.0 );
    element.assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element.initializeYourself();
    element.setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

    Eigen::VectorXd Q = Eigen::VectorXd::Zero( Element::sizeLoadVector );
    for ( int node = 4; node < 8; ++node )
      Q.segment< 3 >( 3 * node ) = Eigen::Vector3d( 1.5e-3, -0.8e-3, 2.0e-3 );
    const double nonlocalBottom = 3.0e-3;
    const double nonlocalTop    = 0.5e-3;
    Q.segment< 4 >( Element::nDofDisplacement ).setConstant( nonlocalBottom );
    Q.segment< 4 >( Element::nDofDisplacement + 4 ).setConstant( nonlocalTop );

    std::vector< double > Pe( Element::sizeLoadVector, 0.0 ),
      Ke( Element::sizeLoadVector * Element::sizeLoadVector, 0.0 );
    element.computeKernels( Q.data(), Q.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

    const double kappaTop    = element.getStateView( "topMaterialStateVars", 0 ).stateLocation[0];
    const double kappaBottom = element.getStateView( "bottomMaterialStateVars", 0 ).stateLocation[0];
    throwExceptionOnFailure( kappaTop <= 1.0e-3 + 1e-14 && std::abs( kappaBottom - nonlocalBottom ) < 1e-14,
                             "Only the face whose field exceeds the threshold is damaged." );

    // The residual of the nonlocal field of a face is the integral over the midsurface of area 1 of h/2 ( N - L ) and
    // of the coupling through the thickness, l^2 / h ( N+ - N- ) with the sign + on the top and - on the bottom face.
    // The radius is 0.2 and the thickness 0.1.
    const double localTop    = element.getStateView( "localFieldPlus", 0 ).stateLocation[0];
    const double localBottom = element.getStateView( "localFieldMinus", 0 ).stateLocation[0];
    const double coupling    = 0.2 * 0.2 / 0.1 * ( nonlocalTop - nonlocalBottom );
    double       residualTop = 0.0, residualBottom = 0.0;
    for ( int i = 0; i < 4; ++i ) {
      residualBottom += Pe[Element::nDofDisplacement + i];
      residualTop += Pe[Element::nDofDisplacement + 4 + i];
    }
    throwExceptionOnFailure( std::abs( residualTop - ( 0.05 * ( nonlocalTop - localTop ) + coupling ) ) < 1e-12 &&
                               std::abs( residualBottom - ( 0.05 * ( nonlocalBottom - localBottom ) - coupling ) ) <
                                 1e-12,
                             "The nonlocal residual of a face is the integral of h/2 ( N - L ) plus the coupling "
                             "through the thickness." );
  }

  /// For the same field on both faces, the coupling through the thickness cancels in the sum over the faces, and the
  /// nonlocal residual of the element with a field for each face sums to the one of the element with a shared field.
  void testFaceFieldsElementAgreesWithTheSharedFieldInTheSum()
  {
    const auto coordinates = inclinedQuadrilateral( 0.0, 0.1 );
    const auto properties  = gradientEnhancedProperties( 0.1 );

    using Shared = GYInterfaceFiniteElementFiniteStrain< 3, 8 >;
    using Faces  = GY2InterfaceFiniteElementFiniteStrain< 3, 8 >;

    const double nonlocal = 2.0e-3;

    Eigen::VectorXd QShared = Eigen::VectorXd::Zero( Shared::sizeLoadVector );
    Eigen::VectorXd QFaces  = Eigen::VectorXd::Zero( Faces::sizeLoadVector );
    for ( int node = 4; node < 8; ++node ) {
      QShared.segment< 3 >( 3 * node ) = Eigen::Vector3d( 1.5e-3, -0.8e-3, 2.0e-3 );
      QFaces.segment< 3 >( 3 * node )  = Eigen::Vector3d( 1.5e-3, -0.8e-3, 2.0e-3 );
    }
    QShared.tail( Shared::nExtraDofs ).setConstant( nonlocal );
    QFaces.tail( Faces::nExtraDofs ).setConstant( nonlocal );

    const auto shared = evaluate< Shared, 3, 8 >( coordinates, properties, QShared, gradientEnhancedMaterial );
    const auto faces  = evaluate< Faces, 3, 8 >( coordinates, properties, QFaces, gradientEnhancedMaterial );

    throwExceptionOnFailure( std::abs( shared.residual.tail( Shared::nExtraDofs ).sum() -
                                       faces.residual.tail( Faces::nExtraDofs ).sum() ) < 1e-12,
                             "The nonlocal residuals must sum to the same value." );
    throwExceptionOnFailure( ( shared.residual.head( 24 ) - faces.residual.head( 24 ) ).norm() <
                               1e-12 * shared.residual.head( 24 ).norm(),
                             "The momentum residuals agree for the same fields." );
  }

  /// The element with a field for each face and the viscoplastic damage model.
  void testFaceFieldsElementWithTheDufourModel()
  {
    using Element = GY2InterfaceFiniteElementFiniteStrain< 3, 8 >;

    const double angle = 10.0 * pi / 180.0;
    const double gap   = 0.1;

    assertTangentMatchesFiniteDifference< Element, 3, 8 >( inclinedQuadrilateral( angle, gap ),
                                                           dufourProperties( gap * std::cos( angle ) ),
                                                           makeDisplacementsAndFaceFields< Element >( 6e-3,
                                                                                                      1,
                                                                                                      0.05,
                                                                                                      0.02 ),
                                                           "GY2 with the Dufour model",
                                                           "DUFOURMODEL" );
  }

  /// the properties of the graded interface material: the thickness, two layers of the given relative thickness
  std::vector< double > gradedNeoHookeProperties( double thickness, double softFraction, double contrast )
  {
    return { thickness,
             0.0,
             2.0,
             2.0,
             2.0,
             softFraction,
             3500.0,
             1500.0,
             1.0 - softFraction,
             contrast * 3500.0,
             contrast * 1500.0 };
  }

  /// A homogeneous graded layer, however it is divided, is the plain element.
  void testGradedElementEqualsThePlainElementForAHomogeneousLayer()
  {
    const double angle = 10.0 * pi / 180.0;
    const double gap   = 0.1;

    const auto compare = []( const Evaluation& expected, const Evaluation& actual, const std::string& context ) {
      throwExceptionOnFailure( ( expected.residual - actual.residual ).norm() <=
                                 1e-9 * std::max( 1.0, expected.residual.norm() ),
                               context + ": the residual of a homogeneous graded layer is not the plain one." );
      throwExceptionOnFailure( ( expected.tangent - actual.tangent ).norm() <=
                                 1e-9 * std::max( 1.0, expected.tangent.norm() ),
                               context + ": the tangent of a homogeneous graded layer is not the plain one." );
    };

    {
      const double thickness = gap * std::cos( angle );
      const auto   Q         = makeDisplacements< 24 >( 3e-2, 1 );
      compare( evaluate< InterfaceFiniteElementFiniteStrain< 3, 8 >, 3, 8 >( inclinedQuadrilateral( angle, gap ),
                                                                             neoHookeProperties( thickness ),
                                                                             Q ),
               evaluate< GradedInterfaceFiniteElementFiniteStrain< 3, 8 >, 3, 8 >( inclinedQuadrilateral( angle, gap ),
                                                                                   gradedNeoHookeProperties( thickness,
                                                                                                             0.4,
                                                                                                             1.0 ),
                                                                                   Q ),
               "3D, inclined" );
    }
    {
      const auto Q = makeDisplacements< 8 >( 5e-2, 3 );
      compare( evaluate< InterfaceFiniteElementFiniteStrain< 2, 4 >, 2, 4 >( inclinedLine( angle, gap ),
                                                                             neoHookeProperties( gap *
                                                                                                 std::cos( angle ) ),
                                                                             Q ),
               evaluate< GradedInterfaceFiniteElementFiniteStrain< 2, 4 >,
                         2,
                         4 >( inclinedLine( angle, gap ),
                              gradedNeoHookeProperties( gap * std::cos( angle ), 0.5, 1.0 ),
                              Q ),
               "2D, inclined" );
    }
  }

  /// The tangent of the graded element is the derivative of its residual, through a stiffness contrast.
  void testGradedElementTangentMatchesFiniteDifference()
  {
    const double angle = 10.0 * pi / 180.0;
    const double gap   = 0.1;

    assertTangentMatchesFiniteDifference< GradedInterfaceFiniteElementFiniteStrain< 3, 8 >,
                                          3,
                                          8 >( inclinedQuadrilateral( angle, gap ),
                                               gradedNeoHookeProperties( gap * std::cos( angle ), 0.4, 6.0 ),
                                               makeDisplacements< 24 >( 3e-2, 1 ),
                                               "graded 3D, inclined" );

    assertTangentMatchesFiniteDifference< GradedInterfaceFiniteElementFiniteStrain< 3, 16 >,
                                          3,
                                          16 >( serendipityQuadrilateral( gap ),
                                                gradedNeoHookeProperties( gap, 0.3, 6.0 ),
                                                makeDisplacements< 48 >( 3e-2, 2 ),
                                                "graded 3D, serendipity" );

    assertTangentMatchesFiniteDifference< GradedInterfaceFiniteElementFiniteStrain< 2, 4 >,
                                          2,
                                          4 >( inclinedLine( angle, gap ),
                                               gradedNeoHookeProperties( gap * std::cos( angle ), 0.5, 6.0 ),
                                               makeDisplacements< 8 >( 5e-2, 3 ),
                                               "graded 2D, inclined" );
  }

  /// the card of the graded gradient-enhanced layer: the thickness, the mode 2, the stations and their bulk properties
  std::vector< double > stationCard( double                       thickness,
                                     const std::vector< double >& positions,
                                     const std::vector< double >& bulk )
  {
    std::vector< double > card = { thickness, 2.0, 1.0, double( positions.size() ), double( bulk.size() ) };
    for ( const double position : positions ) {
      card.push_back( position );
      card.insert( card.end(), bulk.begin(), bulk.end() );
    }
    return card;
  }

  const std::vector< double > neoHookeDamageBulk = { 3500.0, 1500.0, 1.0e-3, 5.0e-3, 0.2, 1.0 };
  const std::vector< double > dufourBulk         = { 1793.9, 828.0, 12.0,   21.6, 4051.0, 15.8, 552.0, 243.0, 0.0, 0.0,
                                                     15.4,   1e-6,  0.0435, 0.3,  0.5,    1.75, 0.99,  0.12,  1.0, 1.0 };

  /// The displacements of the element, and a nonlocal field at every station with a value of its own.
  template < class ElementType, int nStations >
  Eigen::VectorXd makeDisplacementsAndStationFields( double scale, int pattern, double nonlocal )
  {
    constexpr int nPerStation = ElementType::nExtraDofs / nStations;

    Eigen::VectorXd Q = makeDisplacements< ElementType::sizeLoadVector >( scale, pattern );
    for ( int i = 0; i < ElementType::nExtraDofs; ++i ) {
      const int station                      = i / nPerStation;
      Q( ElementType::nDofDisplacement + i ) = nonlocal * ( 1.0 + 0.6 * station ) + 1.0e-4 * ( i % nPerStation );
    }
    return Q;
  }

  /// The element with a nonlocal field at every station: the tangent over the displacements and all fields, with the
  /// damage active.
  void testStationFieldsElement()
  {
    const double angle = 10.0 * pi / 180.0;
    const double gap   = 0.1;
    const double hn    = gap * std::cos( angle );

    using Element3d8  = GYSInterfaceFiniteElementFiniteStrain< 3, 8, 3 >;
    using Element3d16 = GYSInterfaceFiniteElementFiniteStrain< 3, 16, 3 >;
    using Element2d   = GYSInterfaceFiniteElementFiniteStrain< 2, 4, 4 >;

    assertTangentMatchesFiniteDifference< Element3d8, 3, 8 >( inclinedQuadrilateral( angle, gap ),
                                                              stationCard( hn, { 0.0, 0.5, 1.0 }, neoHookeDamageBulk ),
                                                              makeDisplacementsAndStationFields< Element3d8,
                                                                                                 3 >( 3e-3, 1, 3e-3 ),
                                                              "GYS3 3D, inclined",
                                                              gradientEnhancedMaterial );

    assertTangentMatchesFiniteDifference< Element3d16, 3, 16 >( serendipityQuadrilateral( gap ),
                                                                stationCard( gap,
                                                                             { 0.0, 0.3, 1.0 },
                                                                             neoHookeDamageBulk ),
                                                                makeDisplacementsAndStationFields< Element3d16,
                                                                                                   3 >( 3e-3, 2, 4e-3 ),
                                                                "GYS3 3D, serendipity",
                                                                gradientEnhancedMaterial );

    assertTangentMatchesFiniteDifference< Element2d, 2, 4 >( inclinedLine( angle, gap ),
                                                             stationCard( hn,
                                                                          { 0.0, 0.25, 0.7, 1.0 },
                                                                          neoHookeDamageBulk ),
                                                             makeDisplacementsAndStationFields< Element2d, 4 >( 5e-3,
                                                                                                                3,
                                                                                                                3e-3 ),
                                                             "GYS4 2D, inclined",
                                                             gradientEnhancedMaterial );

    assertTangentMatchesFiniteDifference< Element3d8, 3, 8 >( inclinedQuadrilateral( angle, gap ),
                                                              stationCard( hn, { 0.0, 0.5, 1.0 }, dufourBulk ),
                                                              makeDisplacementsAndStationFields< Element3d8,
                                                                                                 3 >( 6e-3, 1, 0.05 ),
                                                              "GYS3 with the Dufour model",
                                                              "DUFOURMODEL" );
  }

  /**
   * With two stations at the faces, the element is the element with a nonlocal field for each face, whenever the two
   * faces have the same surface gradient. The element of this test has one tangential gradient for the layer, so
   * that a virtual displacement moves the surface gradient of both faces together, and the surface stress of the
   * faces works on their sum. The nonlocal residuals are the ones of the other element, and so is the virtual work of
   * every virtual displacement with a common surface gradient of the faces (the top displacement is the bottom one
   * plus a constant). The residuals of the displacements differ in general where the faces carry different surface
   * stresses, which is the modelling difference between the two elements.
   */
  void testTwoStationsAreTheFaceFieldsElement()
  {
    const double angle = 10.0 * pi / 180.0;
    const double gap   = 0.1;
    const double hn    = gap * std::cos( angle );

    using Stations = GYSInterfaceFiniteElementFiniteStrain< 3, 8, 2 >;
    using Faces    = GY2InterfaceFiniteElementFiniteStrain< 3, 8 >;

    const auto coordinates = inclinedQuadrilateral( angle, gap );

    // a displacement of the faces with a common surface gradient: the top is the bottom plus a constant
    const auto commonSurfaceGradient = []( Eigen::VectorXd& Q, const Eigen::Vector3d& offset ) {
      for ( int A = 0; A < 4; ++A )
        Q.segment< 3 >( 3 * ( 4 + A ) ) = Q.segment< 3 >( 3 * A ) + offset;
    };

    for ( const std::string name : std::vector< std::string >{ "DUFOURMODEL", gradientEnhancedMaterial } ) {
      const bool                  dufour         = name == "DUFOURMODEL";
      const std::vector< double > faceProperties = dufour ? dufourProperties( hn ) : gradientEnhancedProperties( hn );
      const std::vector< double > stationProperties = stationCard( hn,
                                                                   { 0.0, 1.0 },
                                                                   dufour ? dufourBulk : neoHookeDamageBulk );

      Eigen::VectorXd Q = makeDisplacementsAndFaceFields< Faces >( dufour ? 6e-3 : 3e-3,
                                                                   1,
                                                                   dufour ? 0.05 : 3e-3,
                                                                   dufour ? 0.02 : 5e-3 );
      commonSurfaceGradient( Q, Eigen::Vector3d( 0.002, -0.001, 0.004 ) );

      const auto faces    = evaluate< Faces, 3, 8 >( coordinates, faceProperties, Q, name );
      const auto stations = evaluate< Stations, 3, 8 >( coordinates, stationProperties, Q, name );

      constexpr int nDofD = 24;
      constexpr int nDofA = 4;

      throwExceptionOnFailure( faces.residual.head( nDofD ).norm() > 0.0, name + ": the test needs a residual." );

      // the nonlocal equations of both faces are the same
      throwExceptionOnFailure( ( faces.residual.tail( 2 * nDofA ) - stations.residual.tail( 2 * nDofA ) ).norm() <=
                                 1e-9 * std::max( 1e-3, faces.residual.tail( 2 * nDofA ).norm() ),
                               name + ": the nonlocal residuals of two stations are not the ones of the faces." );

      // the virtual work of displacements with a common surface gradient of the faces is the same
      for ( int pattern = 0; pattern < 4; ++pattern ) {
        Eigen::VectorXd virtualDisplacement = makeDisplacements< Faces::sizeLoadVector >( 1.0, pattern + 3 );
        virtualDisplacement.tail( 2 * nDofA ).setZero();
        commonSurfaceGradient( virtualDisplacement, Eigen::Vector3d( 0.3, -0.2, 0.5 * pattern ) );

        const double workFaces    = faces.residual.dot( virtualDisplacement );
        const double workStations = stations.residual.dot( virtualDisplacement );

        throwExceptionOnFailure( std::abs( workFaces - workStations ) <= 1e-8 * std::max( 1.0, std::abs( workFaces ) ),
                                 name + ": the virtual work of two stations is not the one of the faces, " +
                                   std::to_string( workStations ) + " against " + std::to_string( workFaces ) + "." );
      }

      // the two elements are different elements: without a common surface gradient the residuals of the
      // displacements are not the same where the faces carry different surface stresses
    }
  }

  /// A material with a different number of stations than the element has fields is refused.
  void testStationCountMismatchIsRefused()
  {
    using Element = GYSInterfaceFiniteElementFiniteStrain< 3, 8, 3 >;

    bool refused = false;
    try {
      evaluate< Element, 3, 8 >( inclinedQuadrilateral( 0.0, 0.1 ),
                                 stationCard( 0.1, { 0.0, 1.0 }, neoHookeDamageBulk ),
                                 makeDisplacementsAndStationFields< Element, 3 >( 3e-3, 1, 3e-3 ),
                                 gradientEnhancedMaterial );
    }
    catch ( const std::invalid_argument& ) {
      refused = true;
    }
    throwExceptionOnFailure( refused, "A material with two stations was accepted by an element with three fields." );
  }

  void testRegistration()
  {
    for ( const std::string name :
          { "ILINE2UL",     "IQUAD4UL",     "IQUAD8UL",     "YILINE2UL",    "YIQUAD4UL",    "YIQUAD8UL",
            "GYILINE2UL",   "GYIQUAD4UL",   "GYIQUAD8UL",   "GY2ILINE2UL",  "GY2IQUAD4UL",  "GY2IQUAD8UL",
            "HILINE2UL",    "HIQUAD4UL",    "HIQUAD8UL",    "GYS2ILINE2UL", "GYS2IQUAD4UL", "GYS2IQUAD8UL",
            "GYS3ILINE2UL", "GYS3IQUAD4UL", "GYS3IQUAD8UL", "GYS4ILINE2UL", "GYS4IQUAD4UL", "GYS4IQUAD8UL" } )
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
                                                   testGradientEnhancedElementWithTheDufourModel,
                                                   testFaceFieldsElement,
                                                   testFaceFieldsElementIsObjective,
                                                   testFaceFieldsAreIndependent,
                                                   testFaceFieldsElementAgreesWithTheSharedFieldInTheSum,
                                                   testFaceFieldsElementWithTheDufourModel,
                                                   testEquilibratedElement,
                                                   testEquilibratedElementAgreesWithThePlainElementForAUniformShift,
                                                   testRigidTranslationIsStressFree,
                                                   testCoincidentFacesAreRejected,
                                                   testInvalidMaterialsAreRejected,
                                                   testMaterialAccessors,
                                                   testGradedElementEqualsThePlainElementForAHomogeneousLayer,
                                                   testGradedElementTangentMatchesFiniteDifference,
                                                   testStationFieldsElement,
                                                   testTwoStationsAreTheFaceFieldsElement,
                                                   testStationCountMismatchIsRefused,
                                                   testRegistration };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
