#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/YInterfaceFiniteElement.h"
#include "Marmot/YStabPressureMiniInterfaceFiniteElement.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Elements;

namespace {

  using MiniElement  = YStabPressureMiniInterfaceFiniteElement< 3, 8 >;
  using PlainElement = YInterfaceFiniteElement< 3, 8 >;

  constexpr int nDofU = 24;
  constexpr int nDofs = 28;

  using DofVector = Eigen::Matrix< double, nDofs, 1 >;
  using KeMatrix  = Eigen::Matrix< double, nDofs, nDofs >;

  /// bottom nodes first, then the top nodes at the gap of 0.1
  const std::array< double, 24 > flatCoordinates = { -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0,
                                                     -0.5, -0.5, 0.1, 0.5, -0.5, 0.1, 0.5, 0.5, 0.1, -0.5, 0.5, 0.1 };

  const std::array< double, 24 > skewedCoordinates = { -0.5, -0.5, 0.088163,  0.5,  -0.5, -0.088163,
                                                       0.5,  0.5,  -0.088163, -0.5, 0.5,  0.088163,
                                                       -0.5, -0.5, 0.188163,  0.5,  -0.5, 0.011837,
                                                       0.5,  0.5,  0.011837,  -0.5, 0.5,  0.188163 };

  const std::vector< double > elasticProperties = { 2e5, 0.3, 0.1 };
  const std::vector< double > plasticProperties = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };

  /// Create an element with a material, the state variables initialized, and the thickness 1 and a stabilisation of
  /// 0.2.
  std::unique_ptr< MiniElement > makeMiniElement( const std::array< double, 24 >& coordinates,
                                                  const std::string&              materialName,
                                                  const std::vector< double >&    materialProperties,
                                                  std::vector< double >&          stateVars )
  {
    auto element = std::make_unique< MiniElement >( 7,
                                                    FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                    MiniElement::SectionType::Interface );
    element->assignNodeCoordinates( coordinates.data() );

    static const std::array< double, 2 > properties = { 1.0, 0.2 };
    ElementProperties elementProperties( properties.data(), static_cast< int >( properties.size() ) );
    element->assignProperty( elementProperties );
    element->assignMaterial( materialName, materialProperties.data(), static_cast< int >( materialProperties.size() ) );

    stateVars.assign( element->getNumberOfRequiredStateVars(), 0.0 );
    element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element->initializeYourself();
    element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
    return element;
  }

  struct Evaluation {
    DofVector residual;
    KeMatrix  tangent;
  };

  /// One evaluation, from the committed state if given, which is updated. The total values are those of the increment.
  Evaluation evaluateMini( const std::array< double, 24 >& coordinates,
                           const std::string&              materialName,
                           const std::vector< double >&    materialProperties,
                           const DofVector&                dQ,
                           const DofVector&                QStart,
                           std::vector< double >*          committedState = nullptr )
  {
    std::vector< double > stateVars;
    auto                  element = makeMiniElement( coordinates, materialName, materialProperties, stateVars );
    if ( committedState ) {
      stateVars = *committedState;
      element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    }

    const DofVector                     QTotal = QStart + dQ;
    std::array< double, nDofs >         Pe{};
    std::array< double, nDofs * nDofs > Ke{};
    element->computeKernels( QTotal.data(), dQ.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

    if ( committedState )
      *committedState = stateVars;

    Evaluation evaluation;
    evaluation.residual = Eigen::Map< DofVector >( Pe.data() );
    evaluation.tangent  = Eigen::Map< KeMatrix >( Ke.data() );
    return evaluation;
  }

  DofVector makeIncrement( double scale, int pattern )
  {
    DofVector dQ;
    for ( int i = 0; i < nDofs; ++i )
      dQ( i ) = scale * 1.0e-3 * ( ( ( i * 7 + pattern * 5 ) % 13 ) - 6 ) / 6.0;
    return dQ;
  }

  /// The tangent has to be the derivative of the condensed residual, bubble and pressure included.
  void assertTangentMatchesFiniteDifference( const std::array< double, 24 >& coordinates,
                                             const std::string&              materialName,
                                             const std::vector< double >&    materialProperties,
                                             const DofVector&                dQ,
                                             const DofVector&                QStart,
                                             const std::vector< double >*    committedState,
                                             double                          tolerance,
                                             const std::string&              context )
  {
    const auto evaluateFrom = [&]( const DofVector& increment ) {
      std::vector< double > state = committedState ? *committedState : std::vector< double >();
      return evaluateMini( coordinates,
                           materialName,
                           materialProperties,
                           increment,
                           QStart,
                           committedState ? &state : nullptr );
    };

    const auto base = evaluateFrom( dQ );
    throwExceptionOnFailure( base.residual.allFinite() && base.tangent.allFinite(),
                             context + ": residual or tangent contains nan or inf." );
    throwExceptionOnFailure( base.residual.norm() > 0.0, context + ": the increment must produce a residual." );

    KeMatrix finiteDifference;
    for ( int i = 0; i < nDofs; ++i ) {
      const double perturbation = 1e-7 * std::max( 1.0, std::abs( dQ( i ) ) );

      DofVector plus  = dQ;
      DofVector minus = dQ;
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
    assertTangentMatchesFiniteDifference( flatCoordinates,
                                          "LINEARELASTIC",
                                          elasticProperties,
                                          makeIncrement( 0.1, 1 ),
                                          DofVector::Zero(),
                                          nullptr,
                                          1e-6,
                                          "Elastic, flat" );
    assertTangentMatchesFiniteDifference( skewedCoordinates,
                                          "LINEARELASTIC",
                                          elasticProperties,
                                          makeIncrement( 0.1, 4 ),
                                          DofVector::Zero(),
                                          nullptr,
                                          1e-6,
                                          "Elastic, skewed" );
  }

  /// A plastic history on the faces, so that the tangent of the two faces differs and the cross blocks matter.
  std::vector< double > makePlasticHistory( const std::array< double, 24 >& coordinates, DofVector& QAfterHistory )
  {
    std::vector< double > state;
    {
      std::vector< double > initial;
      auto                  element = makeMiniElement( coordinates, "VONMISES", plasticProperties, initial );
      state                         = initial;
    }
    const DofVector dQ = makeIncrement( 1.0, 2 );
    evaluateMini( coordinates, "VONMISES", plasticProperties, dQ, DofVector::Zero(), &state );
    QAfterHistory = dQ;
    return state;
  }

  void testTangentMatchesResidualFiniteDifferencePlastic()
  {
    for ( const auto* coordinates : { &flatCoordinates, &skewedCoordinates } ) {
      DofVector  QHistory;
      const auto history = makePlasticHistory( *coordinates, QHistory );
      assertTangentMatchesFiniteDifference( *coordinates,
                                            "VONMISES",
                                            plasticProperties,
                                            makeIncrement( 0.25, 3 ),
                                            QHistory,
                                            &history,
                                            1e-5,
                                            coordinates == &flatCoordinates ? "Plastic, flat" : "Plastic, skewed" );
    }
  }

  /// A uniform strain with the pressure of the constitutive law is stressed like the plain element, and the
  /// volumetric residual and the residual of the bubble vanish.
  void testAgreesWithThePlainElementForAUniformStrain()
  {
    const double a = 2.0e-3;
    const double b = -1.0e-3;
    const double w = 3.0e-4;

    DofVector dQ = DofVector::Zero();
    for ( int node = 0; node < 8; ++node ) {
      const double x     = flatCoordinates[3 * node];
      const double y     = flatCoordinates[3 * node + 1];
      dQ( 3 * node )     = a * x;
      dQ( 3 * node + 1 ) = b * y;
      dQ( 3 * node + 2 ) = node >= 4 ? w : 0.0;
    }

    const double E = elasticProperties[0], nu = elasticProperties[1], h = elasticProperties[2];
    const double K                = E / ( 3.0 * ( 1.0 - 2.0 * nu ) );
    const double volumetricStrain = a + b + w / h;
    for ( int node = 0; node < 4; ++node )
      dQ( nDofU + node ) = -K * volumetricStrain;

    const auto mini = evaluateMini( flatCoordinates, "LINEARELASTIC", elasticProperties, dQ, DofVector::Zero() );

    std::vector< double > plainState;
    auto                  plain = std::make_unique< PlainElement >( 7,
                                                   FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                   PlainElement::SectionType::Interface );
    plain->assignNodeCoordinates( flatCoordinates.data() );
    static const std::array< double, 1 > thickness = { 1.0 };
    ElementProperties                    properties( thickness.data(), 1 );
    plain->assignProperty( properties );
    plain->assignMaterial( "LINEARELASTIC", elasticProperties.data(), static_cast< int >( elasticProperties.size() ) );
    plainState.assign( plain->getNumberOfRequiredStateVars(), 0.0 );
    plain->assignStateVars( plainState.data(), static_cast< int >( plainState.size() ) );
    plain->initializeYourself();
    plain->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

    std::array< double, nDofU >         Q{};
    std::array< double, nDofU >         Pe{};
    std::array< double, nDofU * nDofU > Ke{};
    std::array< double, nDofU >         dU{};
    std::copy( dQ.data(), dQ.data() + nDofU, dU.begin() );
    plain->computeKernels( Q.data(), dU.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

    const Eigen::Map< Eigen::Matrix< double, nDofU, 1 > > plainResidual( Pe.data() );
    const double                                          scale = std::max( 1.0, plainResidual.norm() );
    throwExceptionOnFailure( ( mini.residual.head< nDofU >() - plainResidual ).norm() / scale < 1e-9,
                             "The displacement residual of the uniform strain differs from the plain element by " +
                               std::to_string( ( mini.residual.head< nDofU >() - plainResidual ).norm() / scale ) );
    throwExceptionOnFailure( mini.residual.tail< 4 >().norm() < 1e-9 * K,
                             "The pressure residual of the uniform strain does not vanish: " +
                               std::to_string( mini.residual.tail< 4 >().norm() ) );
  }

  void testRigidTranslationIsStressFree()
  {
    DofVector translation = DofVector::Zero();
    for ( int node = 0; node < 8; ++node )
      translation.segment< 3 >( 3 * node ) = Eigen::Vector3d( 0.3, -0.2, 0.1 );

    const auto rigid = evaluateMini( skewedCoordinates,
                                     "LINEARELASTIC",
                                     elasticProperties,
                                     translation,
                                     DofVector::Zero() );
    throwExceptionOnFailure( rigid.residual.norm() < 1e-9 * 2e5, "A rigid translation must not cause forces." );
  }

  void testDegreesOfFreedomAndNodalFields()
  {
    std::unique_ptr< MarmotElement > element(
      MarmotLibrary::MarmotElementFactory::createElement( "YIQUAD4_STABP_MINI", 1 ) );
    throwExceptionOnFailure( element != nullptr, "YIQUAD4_STABP_MINI is not registered with the element factory." );
    throwExceptionOnFailure( element->getElementShape() == "hexa8", "Wrong result-file shape." );
    throwExceptionOnFailure( element->getNDofPerElement() == nDofs,
                             "The element has 24 displacements and 4 pressures." );
    throwExceptionOnFailure( element->getNNodes() == 8, "The element has 8 nodes." );

    const auto fields = element->getNodeFields();
    throwExceptionOnFailure( fields.size() == 8, "One list of fields per node." );
    for ( int node = 0; node < 8; ++node ) {
      const bool hasPressure = std::find( fields[node].begin(), fields[node].end(), "interfacePressure" ) !=
                               fields[node].end();
      throwExceptionOnFailure( hasPressure == ( node < 4 ),
                               "The pressure lives on the nodes of the bottom side only." );
    }

    // 4 bottom nodes with [u(3), p] and 4 top nodes with u(3): the pattern is a permutation of 0..27
    auto permutation = element->getDofIndicesPermutationPattern();
    throwExceptionOnFailure( permutation.size() == nDofs, "The pattern places every degree of freedom." );
    std::sort( permutation.begin(), permutation.end() );
    std::vector< int > expected( nDofs );
    std::iota( expected.begin(), expected.end(), 0 );
    throwExceptionOnFailure( permutation == expected, "The pattern has to be a permutation." );
  }

  void testStateVariablesAreNamed()
  {
    std::vector< double > stateVars;
    auto                  element = makeMiniElement( flatCoordinates, "LINEARELASTIC", elasticProperties, stateVars );
    for ( const std::string name : { "generalizedForce",
                                     "surfaceStressPlus",
                                     "surfaceStressMinus",
                                     "volumetricResidual",
                                     "bubbleAlpha",
                                     "pressure" } ) {
      const StateView view = element->getStateView( name, 0 );
      throwExceptionOnFailure( view.stateSize > 0, name + " is not a state variable of the element." );
    }
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testTangentMatchesResidualFiniteDifferenceElastic,
                                                   testTangentMatchesResidualFiniteDifferencePlastic,
                                                   testAgreesWithThePlainElementForAUniformStrain,
                                                   testRigidTranslationIsStressFree,
                                                   testDegreesOfFreedomAndNodalFields,
                                                   testStateVariablesAreNamed };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
