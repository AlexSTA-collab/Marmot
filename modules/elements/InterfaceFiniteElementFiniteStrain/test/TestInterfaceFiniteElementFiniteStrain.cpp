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
                       const DofVector< nDim, nNodes >&   Q )
  {
    constexpr int nDofs = nDim * nNodes;

    ElementType element( 3,
                         FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                         ElementType::SectionType::Interface );
    element.assignNodeCoordinates( coordinates.data() );

    const std::array< double, 1 > thickness = { 1.0 };
    ElementProperties             properties( thickness.data(), 1 );
    element.assignProperty( properties );
    element.assignMaterial( "COMPRESSIBLENEOHOOKE",
                            materialProperties.data(),
                            static_cast< int >( materialProperties.size() ) );

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
    evaluation.tangent  = Eigen::Map< Eigen::Matrix< double, nDofs, nDofs, Eigen::RowMajor > >( Ke.data() );
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
                                             const DofVector< nDim, nNodes >&   Q,
                                             const std::string&                 context )
  {
    constexpr int nDofs = nDim * nNodes;

    const auto base = evaluate< ElementType, nDim, nNodes >( coordinates, materialProperties, Q );
    throwExceptionOnFailure( base.residual.allFinite() && base.tangent.allFinite(), context + ": nan or inf." );
    throwExceptionOnFailure( base.residual.norm() > 0.0, context + ": the displacements must produce a residual." );

    Eigen::MatrixXd finiteDifference( nDofs, nDofs );
    for ( int i = 0; i < nDofs; ++i ) {
      const double eps = 1e-6;

      DofVector< nDim, nNodes > plus  = Q;
      DofVector< nDim, nNodes > minus = Q;
      plus( i ) += eps;
      minus( i ) -= eps;

      finiteDifference.col(
        i ) = ( evaluate< ElementType, nDim, nNodes >( coordinates, materialProperties, plus ).residual -
                evaluate< ElementType, nDim, nNodes >( coordinates, materialProperties, minus ).residual ) /
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

  void testRegistration()
  {
    for ( const std::string name : { "ILINE2UL", "IQUAD4UL", "IQUAD8UL", "YILINE2UL", "YIQUAD4UL", "YIQUAD8UL" } )
      throwExceptionOnFailure( std::unique_ptr< MarmotElement >(
                                 MarmotLibrary::MarmotElementFactory::createElement( name, 1 ) ) != nullptr,
                               name + " is not registered with the element factory." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testTangentMatchesFiniteDifference,
                                                   testRigidRotationIsStressFree,
                                                   testEquilibratedElement,
                                                   testEquilibratedElementAgreesWithThePlainElementForAUniformShift,
                                                   testRigidTranslationIsStressFree,
                                                   testCoincidentFacesAreRejected,
                                                   testInvalidMaterialsAreRejected,
                                                   testRegistration };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
