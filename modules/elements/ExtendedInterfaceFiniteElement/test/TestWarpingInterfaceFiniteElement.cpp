#include "Marmot/GaussLobattoInterfaceFiniteElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/WarpingInterfaceFiniteElement.h"

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Elements;

namespace {

  using WElement  = WarpingInterfaceFiniteElement< 3, 8 >;
  using GLElement = GaussLobattoInterfaceFiniteElement< 3, 8 >;

  constexpr int nElementDofs = 24;

  std::array< double, 24 > flatCoordinates = {
    -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0,
    -0.5, -0.5, 0.1, 0.5, -0.5, 0.1, 0.5, 0.5, 0.1, -0.5, 0.5, 0.1,
  };

  std::array< double, 24 > skewedCoordinates = {
    -0.500000, -0.500000, 0.088163,  0.500000,  -0.500000, -0.088163,
    0.500000,  0.500000,  -0.088163, -0.500000, 0.500000,  0.088163,

    -0.500000, -0.500000, 0.188163,  0.500000,  -0.500000, 0.011837,
    0.500000,  0.500000,  0.011837,  -0.500000, 0.500000,  0.188163,
  };

  const std::vector< double > elasticProperties = { 1.0e5, 0.3, 0.01 };
  //                                      E,     nu,  h,    fy,   H,     Aexp, Hexp, id
  const std::vector< double > plasticProperties = { 1.0e5, 0.3, 0.01, 5.0, 1.0e3, 0.0, 0.0, 0.0 };

  template < class ElementT >
  std::unique_ptr< ElementT > makeElement( const std::string&           materialName,
                                           const std::vector< double >& materialProperties,
                                           bool                         skewedGeometry )
  {
    auto element = std::make_unique< ElementT >( 7,
                                                 FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                 ElementT::SectionType::Interface );
    element->assignNodeCoordinates( skewedGeometry ? skewedCoordinates.data() : flatCoordinates.data() );

    static std::array< double, 1 > elementPropertyValues = { 1.0 };
    ElementProperties              elementProperties( elementPropertyValues.data(),
                                         static_cast< int >( elementPropertyValues.size() ) );
    element->assignProperty( elementProperties );
    element->assignMaterial( materialName, materialProperties.data(), static_cast< int >( materialProperties.size() ) );
    return element;
  }

  struct ElementEvaluation {
    Eigen::Matrix< double, nElementDofs, 1 >                             residual;
    Eigen::Matrix< double, nElementDofs, nElementDofs, Eigen::RowMajor > tangent;
    std::vector< double >                                                stateOut;
  };

  template < class ElementT >
  ElementEvaluation evaluate( const std::string&                              materialName,
                              const std::vector< double >&                    materialProperties,
                              const Eigen::Matrix< double, nElementDofs, 1 >& dU,
                              const std::vector< double >*                    initialStateVars = nullptr,
                              bool                                            skewedGeometry   = false )
  {
    auto element = makeElement< ElementT >( materialName, materialProperties, skewedGeometry );

    std::vector< double > stateVars;
    if ( initialStateVars ) {
      stateVars = *initialStateVars;
      element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
      element->initializeYourself();
    }
    else {
      stateVars.assign( element->getNumberOfRequiredStateVars(), 0.0 );
      element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
      element->initializeYourself();
      element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
    }

    std::array< double, nElementDofs >                total{};
    std::array< double, nElementDofs >                increment{};
    std::array< double, nElementDofs >                Pe{};
    std::array< double, nElementDofs * nElementDofs > Ke{};
    std::copy( dU.data(), dU.data() + dU.size(), increment.begin() );

    element->computeKernels( total.data(), increment.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

    ElementEvaluation evaluation;
    evaluation.residual = Eigen::Map< Eigen::Matrix< double, nElementDofs, 1 > >( Pe.data() );
    evaluation.tangent  = Eigen::Map< Eigen::Matrix< double, nElementDofs, nElementDofs, Eigen::RowMajor > >(
      Ke.data() );
    evaluation.stateOut = stateVars;
    return evaluation;
  }

  Eigen::Matrix< double, nElementDofs, 1 > makeGeneralIncrement()
  {
    constexpr int                    half = nElementDofs / 2;
    const std::array< double, half > bottomPattern =
      { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4, -0.5e-4, 0.8e-4, 1.0e-4, 0.2e-4, -0.3e-4, 0.4e-4 };
    const std::array< double, 3 > topOffset = { 2.0e-5, -1.0e-5, 3.0e-5 };

    Eigen::Matrix< double, nElementDofs, 1 > dU;
    for ( int i = 0; i < half; ++i ) {
      dU( i )        = bottomPattern[i];
      dU( half + i ) = bottomPattern[i] + topOffset[i % 3];
    }
    return dU;
  }

  /** Uniform through-thickness loading: both faces get the SAME in-plane field. */
  Eigen::Matrix< double, nElementDofs, 1 > makeUniformIncrement()
  {
    constexpr int                            half = nElementDofs / 2;
    Eigen::Matrix< double, nElementDofs, 1 > dU   = Eigen::Matrix< double, nElementDofs, 1 >::Zero();

    // rigid in-plane shear of both faces by the same amount + a normal opening
    for ( int node = 0; node < 4; ++node ) {
      const double x = flatCoordinates[3 * node + 0];
      const double y = flatCoordinates[3 * node + 1];

      const double ux = 1.0e-4 * x;
      const double uy = 0.5e-4 * y;

      dU( 3 * node + 0 )        = ux;
      dU( 3 * node + 1 )        = uy;
      dU( 3 * node + 2 )        = 0.0;
      dU( half + 3 * node + 0 ) = ux;
      dU( half + 3 * node + 1 ) = uy;
      dU( half + 3 * node + 2 ) = 2.0e-5;
    }
    return dU;
  }

  // ------------------------------------------------------------------
  // 1. Condensed element tangent vs central differences of the condensed
  //    residual. This is the test that actually validates the static
  //    condensation, the inner warping Newton and the pseudo-inverse
  //    together: if the internal problem were not solved to convergence, or
  //    if the Schur complement used the wrong operator, the residual and the
  //    tangent would disagree.
  // ------------------------------------------------------------------
  void checkTangentAgainstFiniteDifference( const std::string&           materialName,
                                            const std::vector< double >& materialProperties,
                                            bool                         skewedGeometry,
                                            const std::string&           label )
  {
    const auto dU = makeGeneralIncrement();

    const auto reference = evaluate< WElement >( materialName, materialProperties, dU, nullptr, skewedGeometry );

    Eigen::Matrix< double, nElementDofs, nElementDofs > finiteDifference;
    const double                                        perturbation = 1.0e-9;
    for ( int column = 0; column < nElementDofs; ++column ) {
      auto forward = dU, backward = dU;
      forward( column ) += perturbation;
      backward( column ) -= perturbation;

      const auto plus  = evaluate< WElement >( materialName, materialProperties, forward, nullptr, skewedGeometry );
      const auto minus = evaluate< WElement >( materialName, materialProperties, backward, nullptr, skewedGeometry );

      // Pe is MINUS the internal residual, so d(Pe)/d(u) = -K.
      finiteDifference.col( column ) = -( plus.residual - minus.residual ) / ( 2.0 * perturbation );
    }

    const double error         = ( reference.tangent - finiteDifference ).lpNorm< Eigen::Infinity >();
    const double scale         = std::max( 1.0, finiteDifference.lpNorm< Eigen::Infinity >() );
    const double relativeError = error / scale;

    std::cout << "condensed element tangent vs FD (" << label << "): relative error = " << relativeError << "\n";

    throwExceptionOnFailure( relativeError < 1e-6,
                             "WarpingInterfaceFiniteElement condensed tangent does not match central differences (" +
                               label + "): relative error " + std::to_string( relativeError ) );
  }

  void TestElasticTangentMatchesFiniteDifference()
  {
    checkTangentAgainstFiniteDifference( "LINEARELASTIC", elasticProperties, false, "elastic, flat" );
  }

  void TestPlasticTangentMatchesFiniteDifference()
  {
    checkTangentAgainstFiniteDifference( "VONMISES", plasticProperties, false, "plastic, flat" );
  }

  void TestSkewedPlasticTangentMatchesFiniteDifference()
  {
    checkTangentAgainstFiniteDifference( "VONMISES", plasticProperties, true, "plastic, skewed" );
  }

  // ------------------------------------------------------------------
  // 2. The internal warping problem must actually be solved: the converged
  //    amplitudes must be nonzero for a general increment (otherwise the
  //    element would silently degenerate to GLIQUAD4), and the element must
  //    remain well-posed despite the documented null space of the internal
  //    block.
  // ------------------------------------------------------------------
  void TestWarpingAmplitudesAreActivatedAndStored()
  {
    const auto dU      = makeGeneralIncrement();
    auto       element = makeElement< WElement >( "LINEARELASTIC", elasticProperties, false );

    std::vector< double > stateVars( element->getNumberOfRequiredStateVars(), 0.0 );
    element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element->initializeYourself();
    element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

    std::array< double, nElementDofs >                total{};
    std::array< double, nElementDofs >                increment{};
    std::array< double, nElementDofs >                Pe{};
    std::array< double, nElementDofs * nElementDofs > Ke{};
    std::copy( dU.data(), dU.data() + dU.size(), increment.begin() );
    element->computeKernels( total.data(), increment.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

    const auto view = element->getStateView( "warpingAmplitudes", 0 );
    double     norm = 0.0;
    for ( int i = 0; i < view.stateSize; ++i ) {
      norm = std::max( norm, std::abs( view.stateLocation[i] ) );
    }

    std::cout << "converged warping amplitude magnitude = " << norm << " over " << view.stateSize << " internal DOF\n";

    throwExceptionOnFailure( view.stateSize == WElement::nWarpDof,
                             "unexpected number of stored warping amplitudes: " + std::to_string( view.stateSize ) );
    throwExceptionOnFailure( norm > 0.0,
                             "the internal warping amplitudes stayed exactly zero -- the enrichment is inactive." );
    throwExceptionOnFailure( std::isfinite( norm ),
                             "the internal warping amplitudes are not finite -- the condensation is ill-posed." );

    // Every entry of the condensed operators must be finite: this is what a
    // failed pseudo-inverse would break first.
    for ( int i = 0; i < nElementDofs; ++i ) {
      throwExceptionOnFailure( std::isfinite( Pe[i] ), "condensed residual has a non-finite entry." );
      for ( int j = 0; j < nElementDofs; ++j ) {
        throwExceptionOnFailure( std::isfinite( Ke[i * nElementDofs + j] ),
                                 "condensed tangent has a non-finite entry." );
      }
    }
  }

  // ------------------------------------------------------------------
  // 3. THE PATCH TEST.
  //
  // Under a through-thickness-UNIFORM in-plane field the exact solution has no
  // warping, so a consistent enrichment must leave the state alone: the
  // converged amplitudes must be zero and the RESIDUAL must equal the
  // Gauss-Lobatto element's exactly. This is what fails if the warping modes do
  // not vanish on the element boundary (non-bubble modes {xi, eta, xi*eta} give
  // a 21% residual error here, because int_Ae grad_s M dA no longer vanishes and
  // a uniform stress state drives spurious warping).
  //
  // The TANGENT is deliberately NOT required to match: even where the converged
  // amplitudes are zero, the Schur complement -K_db K_bb^+ K_bd is not, because
  // the warping can relax PERTURBATIONS about the uniform state. A softer
  // tangent is the enrichment doing its job; only a stiffer one would be wrong.
  // ------------------------------------------------------------------
  void TestUniformLoadingSatisfiesThePatchTest()
  {
    const auto dU = makeUniformIncrement();

    const auto warping      = evaluate< WElement >( "LINEARELASTIC", elasticProperties, dU );
    const auto gaussLobatto = evaluate< GLElement >( "LINEARELASTIC", elasticProperties, dU );

    const double residualError = ( warping.residual - gaussLobatto.residual ).lpNorm< Eigen::Infinity >();
    const double residualScale = std::max( 1.0, gaussLobatto.residual.lpNorm< Eigen::Infinity >() );

    const double warpingEnergy      = dU.transpose() * warping.tangent * dU;
    const double gaussLobattoEnergy = dU.transpose() * gaussLobatto.tangent * dU;

    std::cout << "patch test (uniform loading): residual rel. error = " << residualError / residualScale
              << ", energy WIQUAD4/GLIQUAD4 = " << warpingEnergy / gaussLobattoEnergy << "\n";

    throwExceptionOnFailure( residualError / residualScale < 1e-12,
                             "PATCH TEST FAILED: under through-thickness-uniform loading the warping element must "
                             "reproduce the Gauss-Lobatto residual exactly, relative error " +
                               std::to_string( residualError / residualScale ) );
    throwExceptionOnFailure( warpingEnergy <= gaussLobattoEnergy * ( 1.0 + 1e-10 ),
                             "the enrichment must not stiffen the uniform state." );
  }

  // ------------------------------------------------------------------
  // 4. The enrichment must SOFTEN the element: adding internal degrees of
  //    freedom that are relaxed to equilibrium can only lower the energy for
  //    a given boundary displacement, so the condensed stiffness cannot
  //    exceed the unenriched one in the energy norm. A violation would mean
  //    the Schur complement carries the wrong sign somewhere.
  // ------------------------------------------------------------------
  void TestCondensationSoftensRelativeToGaussLobatto()
  {
    const auto dU = makeGeneralIncrement();

    const auto warping      = evaluate< WElement >( "LINEARELASTIC", elasticProperties, dU );
    const auto gaussLobatto = evaluate< GLElement >( "LINEARELASTIC", elasticProperties, dU );

    const double warpingEnergy      = dU.transpose() * warping.tangent * dU;
    const double gaussLobattoEnergy = dU.transpose() * gaussLobatto.tangent * dU;

    std::cout << "energy in the increment direction: WIQUAD4 = " << warpingEnergy
              << ", GLIQUAD4 = " << gaussLobattoEnergy << "\n";

    throwExceptionOnFailure( warpingEnergy <= gaussLobattoEnergy * ( 1.0 + 1e-10 ),
                             "static condensation of the warping fields must not stiffen the element: " +
                               std::to_string( warpingEnergy ) + " > " + std::to_string( gaussLobattoEnergy ) );
    throwExceptionOnFailure( warpingEnergy > 0.0, "the condensed element must retain positive energy." );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{ TestElasticTangentMatchesFiniteDifference,
                                                       TestPlasticTangentMatchesFiniteDifference,
                                                       TestSkewedPlasticTangentMatchesFiniteDifference,
                                                       TestWarpingAmplitudesAreActivatedAndStored,
                                                       TestUniformLoadingSatisfiesThePatchTest,
                                                       TestCondensationSoftensRelativeToGaussLobatto };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
