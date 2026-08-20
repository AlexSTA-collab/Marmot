/* GLIQUAD4_BUBBLE: Q1 + MINI displacement bubble, condensed, no pressure DOF.
 *
 * The tests pin the two claims the experiment rests on:
 *   1. with the bubble suppressed the element IS GLIQUAD4, bit-for-bit -- so
 *      the only difference between the Q1 and BUBBLE_ONLY production runs is
 *      the bubble;
 *   2. with the bubble active the condensed tangent is consistent.
 */
#include "Marmot/BubbleOnlyInterfaceFiniteElement.h"
#include "Marmot/GaussLobattoInterfaceFiniteElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Elements;

namespace {

  using BubbleElement = BubbleOnlyInterfaceFiniteElement;
  using GLElement     = GaussLobattoInterfaceFiniteElement< 3, 8 >;

  constexpr int nDofs = 24;

  std::array< double, 24 > flatCoordinates = {
    -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0,
    -0.5, -0.5, 0.1, 0.5, -0.5, 0.1, 0.5, 0.5, 0.1, -0.5, 0.5, 0.1,
  };

  const std::vector< double > elasticProperties = { 1.0e5, 0.3, 0.01 };
  const std::vector< double > plasticProperties = { 1.0e5, 0.3, 0.01, 5.0, 1.0e3, 0.0, 0.0, 0.0 };

  bool suppressBubble = false;

  template < class ElementT >
  std::unique_ptr< ElementT > makeElement( const std::string& materialName, const std::vector< double >& properties )
  {
    std::unique_ptr< ElementT > element;
    if constexpr ( std::is_same_v< ElementT, BubbleOnlyInterfaceFiniteElement > ) {
      element = std::make_unique< ElementT >( 7,
                                              FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                              ElementT::SectionType::Interface,
                                              suppressBubble );
    }
    else {
      element = std::make_unique< ElementT >( 7,
                                              FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                              ElementT::SectionType::Interface );
    }
    element->assignNodeCoordinates( flatCoordinates.data() );
    static std::array< double, 1 > elementPropertyValues = { 1.0 };
    ElementProperties              elementProperties( elementPropertyValues.data(), 1 );
    element->assignProperty( elementProperties );
    element->assignMaterial( materialName, properties.data(), static_cast< int >( properties.size() ) );
    return element;
  }

  struct Evaluation {
    Eigen::Matrix< double, nDofs, 1 >     residual;
    Eigen::Matrix< double, nDofs, nDofs > tangent;
    double                                bubbleMagnitude = 0.0;
  };

  template < class ElementT >
  Evaluation evaluate( const std::string&                       materialName,
                       const std::vector< double >&             properties,
                       const Eigen::Matrix< double, nDofs, 1 >& dU,
                       bool                                     readBubble = false )
  {
    auto element = makeElement< ElementT >( materialName, properties );

    std::vector< double > stateVars( element->getNumberOfRequiredStateVars(), 0.0 );
    element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element->initializeYourself();
    element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

    std::array< double, nDofs >         total{};
    std::array< double, nDofs >         increment{};
    std::array< double, nDofs >         Pe{};
    std::array< double, nDofs * nDofs > Ke{};
    std::copy( dU.data(), dU.data() + dU.size(), increment.begin() );
    element->computeKernels( total.data(), increment.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

    Evaluation evaluation;
    evaluation.residual = Eigen::Map< Eigen::Matrix< double, nDofs, 1 > >( Pe.data() );
    evaluation.tangent  = Eigen::Map< Eigen::Matrix< double, nDofs, nDofs > >( Ke.data() ); // column-major
    if ( readBubble ) {
      const auto view = element->getStateView( "bubbleAlpha", 0 );
      for ( int i = 0; i < view.stateSize; ++i ) {
        evaluation.bubbleMagnitude = std::max( evaluation.bubbleMagnitude, std::abs( view.stateLocation[i] ) );
      }
    }
    return evaluation;
  }

  Eigen::Matrix< double, nDofs, 1 > makeGeneralIncrement()
  {
    constexpr int                    half = nDofs / 2;
    const std::array< double, half > bottomPattern =
      { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4, -0.5e-4, 0.8e-4, 1.0e-4, 0.2e-4, -0.3e-4, 0.4e-4 };
    const std::array< double, 3 > topOffset = { 2.0e-5, -1.0e-5, 3.0e-5 };

    Eigen::Matrix< double, nDofs, 1 > dU;
    for ( int i = 0; i < half; ++i ) {
      dU( i )        = bottomPattern[i];
      dU( half + i ) = bottomPattern[i] + topOffset[i % 3];
    }
    return dU;
  }

  // ------------------------------------------------------------------
  // 1. Suppressing the bubble must reproduce GLIQUAD4 EXACTLY. This is what
  //    licenses reading the production Q1 run and the BUBBLE_ONLY run as a
  //    controlled pair: it proves the material reduction (warping slots zero ->
  //    Gauss-Lobatto) and that the bubble is the only difference.
  // ------------------------------------------------------------------
  void TestWithoutBubbleItIsExactlyGaussLobatto()
  {
    suppressBubble = true;
    const auto dU  = makeGeneralIncrement();

    for ( const auto& [name, properties] :
          std::vector< std::pair< std::string, std::vector< double > > >{ { "LINEARELASTIC", elasticProperties },
                                                                          { "VONMISES", plasticProperties } } ) {
      const auto bubble = evaluate< BubbleElement >( name, properties, dU );
      const auto gl     = evaluate< GLElement >( name, properties, dU );

      const double residualError = ( bubble.residual - gl.residual ).lpNorm< Eigen::Infinity >() /
                                   std::max( 1e-30, gl.residual.lpNorm< Eigen::Infinity >() );
      const double tangentError = ( bubble.tangent - gl.tangent ).lpNorm< Eigen::Infinity >() /
                                  std::max( 1e-30, gl.tangent.lpNorm< Eigen::Infinity >() );

      std::cout << "bubble suppressed vs GLIQUAD4 (" << name << "): residual " << residualError << ", tangent "
                << tangentError << "\n";

      throwExceptionOnFailure( residualError < 1e-12 && tangentError < 1e-12,
                               "with the bubble suppressed the element must reproduce GLIQUAD4 exactly (" + name +
                                 "): residual " + std::to_string( residualError ) + ", tangent " +
                                 std::to_string( tangentError ) );
    }
    suppressBubble = false;
  }

  // ------------------------------------------------------------------
  // 2. With the bubble active: consistent condensed tangent, and the bubble
  //    genuinely activates (otherwise the production run would silently be Q1).
  // ------------------------------------------------------------------
  void TestBubbleTangentIsConsistentAndBubbleActivates()
  {
    const auto dU = makeGeneralIncrement();

    for ( const auto& [name, properties] :
          std::vector< std::pair< std::string, std::vector< double > > >{ { "LINEARELASTIC", elasticProperties },
                                                                          { "VONMISES", plasticProperties } } ) {
      const auto reference = evaluate< BubbleElement >( name, properties, dU, true );

      Eigen::Matrix< double, nDofs, nDofs > finiteDifference;
      const double                          perturbation = 1.0e-9;
      for ( int column = 0; column < nDofs; ++column ) {
        auto forward = dU, backward = dU;
        forward( column ) += perturbation;
        backward( column ) -= perturbation;
        const auto plus                = evaluate< BubbleElement >( name, properties, forward );
        const auto minus               = evaluate< BubbleElement >( name, properties, backward );
        finiteDifference.col( column ) = -( plus.residual - minus.residual ) / ( 2.0 * perturbation );
      }

      const double relativeError = ( reference.tangent - finiteDifference ).lpNorm< Eigen::Infinity >() /
                                   std::max( 1.0, finiteDifference.lpNorm< Eigen::Infinity >() );

      std::cout << "GLIQUAD4_BUBBLE tangent vs FD (" << name << "): " << relativeError << ", bubble amplitude "
                << reference.bubbleMagnitude << "\n";

      throwExceptionOnFailure( relativeError < 1e-6,
                               "condensed tangent does not match central differences (" + name +
                                 "): " + std::to_string( relativeError ) );
      throwExceptionOnFailure( reference.bubbleMagnitude > 0.0,
                               "the bubble amplitude stayed exactly zero -- the enrichment is inactive (" + name +
                                 ")." );
    }
  }

  // ------------------------------------------------------------------
  // 3. Condensing an internal field that is relaxed to equilibrium can only
  //    soften the element.
  // ------------------------------------------------------------------
  void TestBubbleCondensationSoftens()
  {
    const auto dU = makeGeneralIncrement();

    const auto bubble = evaluate< BubbleElement >( "LINEARELASTIC", elasticProperties, dU );
    const auto gl     = evaluate< GLElement >( "LINEARELASTIC", elasticProperties, dU );

    const double bubbleEnergy = dU.transpose() * bubble.tangent * dU;
    const double glEnergy     = dU.transpose() * gl.tangent * dU;

    std::cout << "energy: GLIQUAD4_BUBBLE = " << bubbleEnergy << ", GLIQUAD4 = " << glEnergy << "\n";

    throwExceptionOnFailure( bubbleEnergy > 0.0, "the condensed element must retain positive energy." );
    throwExceptionOnFailure( bubbleEnergy <= glEnergy * ( 1.0 + 1e-10 ),
                             "condensing the bubble must not stiffen the element: " + std::to_string( bubbleEnergy ) +
                               " > " + std::to_string( glEnergy ) );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{ TestWithoutBubbleItIsExactlyGaussLobatto,
                                                       TestBubbleTangentIsConsistentAndBubbleActivates,
                                                       TestBubbleCondensationSoftens };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
