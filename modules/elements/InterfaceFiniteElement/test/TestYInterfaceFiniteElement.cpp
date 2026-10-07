#include "Marmot/InterfaceFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/YInterfaceFiniteElement.h"

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

  /// bottom nodes first, then the top nodes at the gap of 0.1
  const std::array< double, 24 > flatCoordinates = { -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0,
                                                     -0.5, -0.5, 0.1, 0.5, -0.5, 0.1, 0.5, 0.5, 0.1, -0.5, 0.5, 0.1 };

  const std::array< double, 24 > skewedCoordinates = { -0.5, -0.5, 0.088163,  0.5,  -0.5, -0.088163,
                                                       0.5,  0.5,  -0.088163, -0.5, 0.5,  0.088163,
                                                       -0.5, -0.5, 0.188163,  0.5,  -0.5, 0.011837,
                                                       0.5,  0.5,  0.011837,  -0.5, 0.5,  0.188163 };

  /// the top nodes lie below the bottom nodes
  const std::array< double, 24 > invertedCoordinates = { -0.5, -0.5, 0.1, 0.5, -0.5, 0.1,  0.5, 0.5,
                                                         0.1,  -0.5, 0.5, 0.1, -0.5, -0.5, 0.0, 0.5,
                                                         -0.5, 0.0,  0.5, 0.5, 0.0,  -0.5, 0.5, 0.0 };

  /// the top nodes are shifted along the interface only, so that the connector has no normal component
  const std::array< double, 24 > slidingCoordinates = { -0.5, -0.5, 0.0, 0.5, -0.5, 0.0,  0.5, 0.5,
                                                        0.0,  -0.5, 0.5, 0.0, -0.4, -0.5, 0.0, 0.6,
                                                        -0.5, 0.0,  0.6, 0.5, 0.0,  -0.4, 0.5, 0.0 };

  /// a line in 2D, bottom nodes first
  const std::array< double, 8 > lineCoordinates = { 0.0, 0.0, 1.0, 0.05, 0.0, 0.11, 1.0, 0.1 };

  const std::vector< double > elasticProperties = { 2e5, 0.3, 0.1 };
  const std::vector< double > plasticProperties = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };

  template < int nDim, int nNodes >
  using ElementDofVector = Eigen::Matrix< double, nDim * nNodes, 1 >;

  /// Residual and tangent of one evaluation of an element.
  template < int nDim, int nNodes >
  struct Evaluation {
    ElementDofVector< nDim, nNodes >                                       residual;
    Eigen::Matrix< double, nDim * nNodes, nDim * nNodes, Eigen::RowMajor > tangent;
  };

  /**
   * Create an element of the given type on the given nodes, with a material and the state variables initialized.
   * The coordinates have to outlive the element, which only views them.
   */
  template < class ElementType, int nDim, int nNodes >
  std::unique_ptr< ElementType > makeElement( const std::array< double, nDim * nNodes >& coordinates,
                                              const std::string&                         materialName,
                                              const std::vector< double >&               materialProperties,
                                              std::vector< double >&                     stateVars )
  {
    auto element = std::make_unique< ElementType >( 7,
                                                    FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                    ElementType::SectionType::Interface );
    element->assignNodeCoordinates( coordinates.data() );

    static const std::array< double, 1 > thickness = { 1.0 };
    ElementProperties                    properties( thickness.data(), static_cast< int >( thickness.size() ) );
    element->assignProperty( properties );

    element->assignMaterial( materialName, materialProperties.data(), static_cast< int >( materialProperties.size() ) );

    stateVars.assign( element->getNumberOfRequiredStateVars(), 0.0 );
    element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element->initializeYourself();
    element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
    return element;
  }

  /**
   * Evaluate an element for an increment of the displacements, starting from the given committed state, which
   * is updated in place. Without a committed state the element starts from its initial state.
   */
  template < class ElementType, int nDim, int nNodes >
  Evaluation< nDim, nNodes > evaluate( const std::array< double, nDim * nNodes >& coordinates,
                                       const std::string&                         materialName,
                                       const std::vector< double >&               materialProperties,
                                       const ElementDofVector< nDim, nNodes >&    dQ,
                                       std::vector< double >*                     committedState = nullptr )
  {
    constexpr int nDofs = nDim * nNodes;

    std::vector< double > stateVars;
    auto element = makeElement< ElementType, nDim, nNodes >( coordinates, materialName, materialProperties, stateVars );
    if ( committedState ) {
      stateVars = *committedState;
      element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    }

    std::array< double, nDofs >         Q{};
    std::array< double, nDofs >         dQArray{};
    std::array< double, nDofs >         Pe{};
    std::array< double, nDofs * nDofs > Ke{};
    std::copy( dQ.data(), dQ.data() + nDofs, dQArray.begin() );

    element->computeKernels( Q.data(), dQArray.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

    if ( committedState )
      *committedState = stateVars;

    Evaluation< nDim, nNodes > evaluation;
    evaluation.residual = Eigen::Map< ElementDofVector< nDim, nNodes > >( Pe.data() );
    evaluation.tangent  = Eigen::Map< Eigen::Matrix< double, nDofs, nDofs, Eigen::RowMajor > >( Ke.data() );
    return evaluation;
  }

  /// A smooth increment of the displacements, scaled.
  template < int nDofs >
  Eigen::Matrix< double, nDofs, 1 > makeIncrement( double scale, int pattern )
  {
    Eigen::Matrix< double, nDofs, 1 > dQ;
    for ( int i = 0; i < nDofs; ++i )
      dQ( i ) = scale * 1.0e-3 * ( ( ( i * 7 + pattern * 5 ) % 13 ) - 6 ) / 6.0;
    return dQ;
  }

  /**
   * The tangent has to be the derivative of the residual: compare it with the central finite difference of the
   * residual of the same element and the same committed state.
   */
  template < class ElementType, int nDim, int nNodes >
  void assertTangentMatchesFiniteDifference( const std::array< double, nDim * nNodes >& coordinates,
                                             const std::string&                         materialName,
                                             const std::vector< double >&               materialProperties,
                                             const ElementDofVector< nDim, nNodes >&    dQ,
                                             const std::vector< double >*               committedState,
                                             double                                     tolerance,
                                             const std::string&                         context )
  {
    constexpr int nDofs = nDim * nNodes;

    const auto evaluateFrom = [&]( const ElementDofVector< nDim, nNodes >& increment ) {
      std::vector< double > state = committedState ? *committedState : std::vector< double >();
      return evaluate< ElementType, nDim, nNodes >( coordinates,
                                                    materialName,
                                                    materialProperties,
                                                    increment,
                                                    committedState ? &state : nullptr );
    };

    const auto base = evaluateFrom( dQ );
    throwExceptionOnFailure( base.residual.allFinite() && base.tangent.allFinite(),
                             context + ": residual or tangent contains nan or inf." );
    throwExceptionOnFailure( base.residual.norm() > 0.0, context + ": the increment must produce a residual." );

    Eigen::Matrix< double, nDofs, nDofs, Eigen::RowMajor > finiteDifference;
    for ( int i = 0; i < nDofs; ++i ) {
      const double perturbation = 1e-7 * std::max( 1.0, std::abs( dQ( i ) ) );

      ElementDofVector< nDim, nNodes > plus  = dQ;
      ElementDofVector< nDim, nNodes > minus = dQ;
      plus( i ) += perturbation;
      minus( i ) -= perturbation;

      finiteDifference.col( i ) = ( evaluateFrom( plus ).residual - evaluateFrom( minus ).residual ) /
                                  ( 2.0 * perturbation );
    }

    const double error = ( finiteDifference - base.tangent ).norm();
    const double scale = std::max( 1.0, base.tangent.norm() );
    throwExceptionOnFailure( error / scale < tolerance,
                             context +
                               ": the tangent does not match the finite difference of the residual, "
                               "relative error " +
                               std::to_string( error / scale ) );
  }

  void testTangentMatchesResidualFiniteDifferenceElastic()
  {
    assertTangentMatchesFiniteDifference< YInterfaceFiniteElement< 3, 8 >, 3, 8 >( flatCoordinates,
                                                                                   "LINEARELASTIC",
                                                                                   elasticProperties,
                                                                                   makeIncrement< 24 >( 0.1, 1 ),
                                                                                   nullptr,
                                                                                   1e-6,
                                                                                   "Elastic, flat" );
  }

  /// A plastic history on the faces, so that the tangent of the two faces differs and the cross blocks matter.
  template < int nDim, int nNodes >
  std::vector< double > makePlasticHistory( const std::array< double, nDim * nNodes >& coordinates )
  {
    std::vector< double > state;
    {
      std::vector< double > initial;
      auto                  element = makeElement< YInterfaceFiniteElement< nDim, nNodes >, nDim, nNodes >( coordinates,
                                                                                           "VONMISES",
                                                                                           plasticProperties,
                                                                                           initial );
      state                         = initial;
    }
    evaluate< YInterfaceFiniteElement< nDim, nNodes >, nDim, nNodes >( coordinates,
                                                                       "VONMISES",
                                                                       plasticProperties,
                                                                       makeIncrement< nDim * nNodes >( 1.0, 2 ),
                                                                       &state );
    return state;
  }

  void testTangentMatchesResidualFiniteDifferencePlasticAfterHistory()
  {
    const std::vector< double > history = makePlasticHistory< 3, 8 >( flatCoordinates );

    assertTangentMatchesFiniteDifference< YInterfaceFiniteElement< 3, 8 >, 3, 8 >( flatCoordinates,
                                                                                   "VONMISES",
                                                                                   plasticProperties,
                                                                                   makeIncrement< 24 >( 0.25, 3 ),
                                                                                   &history,
                                                                                   1e-5,
                                                                                   "Plastic, flat" );
  }

  void testTangentMatchesResidualFiniteDifferencePlasticSkewed()
  {
    const std::vector< double > history = makePlasticHistory< 3, 8 >( skewedCoordinates );

    assertTangentMatchesFiniteDifference< YInterfaceFiniteElement< 3, 8 >, 3, 8 >( skewedCoordinates,
                                                                                   "VONMISES",
                                                                                   plasticProperties,
                                                                                   makeIncrement< 24 >( 0.25, 3 ),
                                                                                   &history,
                                                                                   1e-5,
                                                                                   "Plastic, skewed" );
  }

  void testTangentMatchesResidualFiniteDifferenceTwoDimensional()
  {
    assertTangentMatchesFiniteDifference< YInterfaceFiniteElement< 2, 4 >, 2, 4 >( lineCoordinates,
                                                                                   "LINEARELASTIC",
                                                                                   elasticProperties,
                                                                                   makeIncrement< 8 >( 0.1, 1 ),
                                                                                   nullptr,
                                                                                   1e-6,
                                                                                   "Elastic, 2D" );

    const std::vector< double > history = makePlasticHistory< 2, 4 >( lineCoordinates );
    assertTangentMatchesFiniteDifference< YInterfaceFiniteElement< 2, 4 >, 2, 4 >( lineCoordinates,
                                                                                   "VONMISES",
                                                                                   plasticProperties,
                                                                                   makeIncrement< 8 >( 0.25, 3 ),
                                                                                   &history,
                                                                                   1e-5,
                                                                                   "Plastic, 2D" );
  }

  /**
   * For a uniform shift of the top nodes, the surface gradients vanish and the jump is the shift. With the thickness
   * equal to the gap of the faces, both faces see the same normal gradient and the response is the one of the plain
   * element, which averages the two sides.
   */
  void testAgreesWithThePlainElementForAUniformShift()
  {
    constexpr int nDofs = 24;

    Eigen::Matrix< double, nDofs, 1 > shift = Eigen::Matrix< double, nDofs, 1 >::Zero();
    for ( int node = 4; node < 8; ++node ) {
      shift( 3 * node + 0 ) = 1.5e-4;
      shift( 3 * node + 1 ) = -0.8e-4;
      shift( 3 * node + 2 ) = 2.0e-4;
    }

    const auto equilibrated = evaluate< YInterfaceFiniteElement< 3, 8 >, 3, 8 >( flatCoordinates,
                                                                                 "LINEARELASTIC",
                                                                                 elasticProperties,
                                                                                 shift );
    const auto plain        = evaluate< InterfaceFiniteElement< 3, 8 >, 3, 8 >( flatCoordinates,
                                                                         "LINEARELASTIC",
                                                                         elasticProperties,
                                                                         shift );

    throwExceptionOnFailure( ( equilibrated.residual - plain.residual ).norm() < 1e-9 * plain.residual.norm(),
                             "The equilibrated element must give the residual of the plain element for a "
                             "uniform shift of the top face." );

    // the shift is a direction in which both residuals agree, so their tangents agree in that direction
    throwExceptionOnFailure( ( ( equilibrated.tangent - plain.tangent ) * shift ).norm() <
                               1e-9 * ( plain.tangent * shift ).norm(),
                             "The tangents must agree in the direction of a uniform shift of the top face." );
  }

  void testNormalPointsFromTheBottomToTheTopFace()
  {
    for ( const auto* coordinates : { &flatCoordinates, &skewedCoordinates, &invertedCoordinates } ) {
      std::vector< double > stateVars;
      auto                  element = makeElement< YInterfaceFiniteElement< 3, 8 >, 3, 8 >( *coordinates,
                                                                           "LINEARELASTIC",
                                                                           elasticProperties,
                                                                           stateVars );

      for ( const auto& qp : element->qps ) {
        const auto& data = qp.formulationData;
        throwExceptionOnFailure( data.separationVector.norm() > 0.0, "The faces are apart at every quadrature point." );
        throwExceptionOnFailure( data.separationVector.dot( qp.normal ) > 0.0,
                                 "The normal must point from the bottom to the top face, also for inverted nodes." );
        throwExceptionOnFailure( std::abs( data.normalSeparation - data.separationVector.dot( qp.normal ) ) < 1e-14,
                                 "The normal separation is the component of the connector along the normal." );
        throwExceptionOnFailure( std::abs( data.tangentialSeparation.dot( qp.normal ) ) < 1e-12,
                                 "The tangential separation has no component along the normal." );
      }
    }
  }

  void testConnectorWithoutNormalComponentThrows()
  {
    std::vector< double > stateVars;

    bool rejected = false;
    try {
      makeElement< YInterfaceFiniteElement< 3, 8 >, 3, 8 >( slidingCoordinates,
                                                            "LINEARELASTIC",
                                                            elasticProperties,
                                                            stateVars );
    }
    catch ( const std::invalid_argument& ) {
      rejected = true;
    }
    throwExceptionOnFailure( rejected, "Faces without a separation along the normal must be rejected." );
  }

  void testStateVariablesAreNamed()
  {
    std::vector< double > stateVars;
    auto                  element = makeElement< YInterfaceFiniteElement< 3, 8 >, 3, 8 >( flatCoordinates,
                                                                         "VONMISES",
                                                                         plasticProperties,
                                                                         stateVars );

    throwExceptionOnFailure( element->getStateView( "generalizedForce", 0 ).stateSize == 3,
                             "generalizedForce has three entries." );
    throwExceptionOnFailure( element->getStateView( "surfaceStressPlus", 0 ).stateSize == 9,
                             "surfaceStressPlus has nine entries." );
    throwExceptionOnFailure( element->getStateView( "surfaceStressMinus", 0 ).stateSize == 9,
                             "surfaceStressMinus has nine entries." );
    throwExceptionOnFailure( element->getStateView( "displacement", 0 ).stateSize == 6,
                             "displacement holds both sides." );
    throwExceptionOnFailure( element->getStateView( "surfaceStrain", 0 ).stateSize == 18,
                             "surfaceStrain holds both sides." );

    // the state variables of the material of one quadrature point follow the 48 entries of the quadrature point
    const int nMaterial = element->material->getNumberOfRequiredStateVars();
    throwExceptionOnFailure( element->getNumberOfRequiredStateVars() == 4 * ( 48 + nMaterial ),
                             "Four quadrature points with 48 own entries and the state of the material each." );
    throwExceptionOnFailure( element->getStateView( "topMaterialStateVars", 0 ).stateSize > 0,
                             "The state of the material is reachable through the element." );
  }

  /// the nodes of a serendipity quadrilateral interface: corners, then midsides, bottom side first; the top is skewed
  std::array< double, 48 > serendipityCoordinates( double skew )
  {
    const std::array< std::array< double, 2 >, 8 > parametric = { { { -0.5, -0.5 },
                                                                    { 0.5, -0.5 },
                                                                    { 0.5, 0.5 },
                                                                    { -0.5, 0.5 },
                                                                    { 0.0, -0.5 },
                                                                    { 0.5, 0.0 },
                                                                    { 0.0, 0.5 },
                                                                    { -0.5, 0.0 } } };

    std::array< double, 48 > x{};
    for ( int side = 0; side < 2; ++side )
      for ( int node = 0; node < 8; ++node ) {
        const double sign              = side == 0 ? -1.0 : 1.0;
        x[3 * ( 8 * side + node )]     = parametric[node][0];
        x[3 * ( 8 * side + node ) + 1] = parametric[node][1];
        x[3 * ( 8 * side + node ) + 2] = 0.05 * ( side == 0 ? 0.0 : 2.0 ) +
                                         sign * skew * parametric[node][0] * parametric[node][1] + 0.05;
      }
    return x;
  }

  void testSerendipityElement()
  {
    const auto flat   = serendipityCoordinates( 0.0 );
    const auto warped = serendipityCoordinates( 0.04 );

    {
      std::vector< double > stateVars;
      auto                  element = makeElement< YInterfaceFiniteElement< 3, 16 >, 3, 16 >( flat,
                                                                             "LINEARELASTIC",
                                                                             elasticProperties,
                                                                             stateVars );
      throwExceptionOnFailure( element->getNumberOfQuadraturePoints() == 9,
                               "Full integration of the Quad8 interface has 9 points." );
    }

    assertTangentMatchesFiniteDifference< YInterfaceFiniteElement< 3, 16 >, 3, 16 >( warped,
                                                                                     "LINEARELASTIC",
                                                                                     elasticProperties,
                                                                                     makeIncrement< 48 >( 0.1, 1 ),
                                                                                     nullptr,
                                                                                     1e-6,
                                                                                     "Elastic, serendipity" );

    const std::vector< double > history = makePlasticHistory< 3, 16 >( warped );
    assertTangentMatchesFiniteDifference< YInterfaceFiniteElement< 3, 16 >, 3, 16 >( warped,
                                                                                     "VONMISES",
                                                                                     plasticProperties,
                                                                                     makeIncrement< 48 >( 0.25, 3 ),
                                                                                     &history,
                                                                                     1e-5,
                                                                                     "Plastic, serendipity" );

    // a rigid translation of the whole element leaves it stress free
    Eigen::Matrix< double, 48, 1 > translation;
    for ( int node = 0; node < 16; ++node )
      translation.segment< 3 >( 3 * node ) = Eigen::Vector3d( 0.3, -0.2, 0.1 );
    const auto rigid = evaluate< YInterfaceFiniteElement< 3, 16 >, 3, 16 >( warped,
                                                                            "LINEARELASTIC",
                                                                            elasticProperties,
                                                                            translation );
    throwExceptionOnFailure( rigid.residual.norm() < 1e-9 * 2e5, "A rigid translation must not cause forces." );
  }

  void testRegistration()
  {
    for ( const std::string name : { "IQUAD8", "YIQUAD8" } )
      throwExceptionOnFailure( std::unique_ptr< MarmotElement >(
                                 MarmotLibrary::MarmotElementFactory::createElement( name, 1 ) ) != nullptr,
                               name + " is not registered with the element factory." );

    const std::pair< std::string, std::string > registered[] = { { "YIQUAD4", "hexa8" }, { "YILINE2", "bar2" } };
    for ( const auto& [name, shape] : registered ) {
      std::unique_ptr< MarmotElement > element( MarmotLibrary::MarmotElementFactory::createElement( name, 1 ) );
      throwExceptionOnFailure( element != nullptr, name + " is not registered with the element factory." );
      throwExceptionOnFailure( element->getElementShape() == shape,
                               name + " reports the wrong result-file shape: " + element->getElementShape() );
    }
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testTangentMatchesResidualFiniteDifferenceElastic,
                                                   testTangentMatchesResidualFiniteDifferencePlasticAfterHistory,
                                                   testTangentMatchesResidualFiniteDifferencePlasticSkewed,
                                                   testTangentMatchesResidualFiniteDifferenceTwoDimensional,
                                                   testAgreesWithThePlainElementForAUniformShift,
                                                   testNormalPointsFromTheBottomToTheTopFace,
                                                   testConnectorWithoutNormalComponentThrows,
                                                   testStateVariablesAreNamed,
                                                   testSerendipityElement,
                                                   testRegistration };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
