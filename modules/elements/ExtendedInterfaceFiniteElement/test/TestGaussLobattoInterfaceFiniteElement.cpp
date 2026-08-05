#include "Marmot/GaussLobattoInterfaceFiniteElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/YInterfaceFiniteElement.h"

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

  static std::array< double, 24 > flatCoordinates = {
    -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0,
    -0.5, -0.5, 0.1, 0.5, -0.5, 0.1, 0.5, 0.5, 0.1, -0.5, 0.5, 0.1,
  };

  static std::array< double, 24 > skewedCoordinates = {
    -0.500000, -0.500000, 0.088163,  0.500000,  -0.500000, -0.088163,
    0.500000,  0.500000,  -0.088163, -0.500000, 0.500000,  0.088163,

    -0.500000, -0.500000, 0.188163,  0.500000,  -0.500000, 0.011837,
    0.500000,  0.500000,  0.011837,  -0.500000, 0.500000,  0.188163,
  };

  template < class ElementT >
  std::unique_ptr< ElementT > makeInterfaceElementWithMaterial( const std::string& materialName,
                                                                const double*      materialProperties,
                                                                int                nMaterialProperties,
                                                                bool               skewedGeometry = false )
  {
    const int  elId    = 7;
    const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
    const auto secType = ElementT::SectionType::Interface;

    auto element = std::make_unique< ElementT >( elId, intType, secType );
    element->assignNodeCoordinates( skewedGeometry ? skewedCoordinates.data() : flatCoordinates.data() );

    static std::array< double, 1 > elPropsVec = { 1.0 };
    ElementProperties              elProps( elPropsVec.data(), static_cast< int >( elPropsVec.size() ) );
    element->assignProperty( elProps );

    element->assignMaterial( materialName, materialProperties, nMaterialProperties );
    return element;
  }

  template < class ElementT >
  void initializeStateAndMaterial( ElementT& element, std::vector< double >& stateVars )
  {
    stateVars.assign( element.getNumberOfRequiredStateVars(), 0.0 );
    element.assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element.initializeYourself();
    element.setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
  }

  struct ElementEvaluation {
    Eigen::Matrix< double, 24, 1 >                   residual;
    Eigen::Matrix< double, 24, 24, Eigen::RowMajor > tangent;
  };

  template < class ElementT >
  ElementEvaluation evaluateElementResponse( const std::string&                    materialName,
                                             const std::vector< double >&          materialProperties,
                                             const Eigen::Matrix< double, 24, 1 >& dU,
                                             const std::vector< double >*          initialStateVars = nullptr,
                                             bool                                  skewedGeometry   = false )
  {
    constexpr int nElementDofs = 24;

    auto element = makeInterfaceElementWithMaterial< ElementT >( materialName,
                                                                 materialProperties.data(),
                                                                 static_cast< int >( materialProperties.size() ),
                                                                 skewedGeometry );

    std::vector< double > stateVars;
    if ( initialStateVars ) {
      stateVars = *initialStateVars;
      element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
      element->initializeYourself();
    }
    else {
      initializeStateAndMaterial( *element, stateVars );
    }

    std::array< double, nElementDofs >                U{};
    std::array< double, nElementDofs >                dQ{};
    std::array< double, nElementDofs >                Pe{};
    std::array< double, nElementDofs * nElementDofs > Ke{};
    std::copy( dU.data(), dU.data() + dU.size(), dQ.begin() );

    element->computeKernels( U.data(), dQ.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

    ElementEvaluation evaluation;
    evaluation.residual = Eigen::Map< Eigen::Matrix< double, nElementDofs, 1 > >( Pe.data() );
    evaluation.tangent  = Eigen::Map< Eigen::Matrix< double, nElementDofs, nElementDofs, Eigen::RowMajor > >(
      Ke.data() );
    return evaluation;
  }

  template < class ElementT >
  std::vector< double > makeElementStateAfterIncrement( const std::string&                    materialName,
                                                        const std::vector< double >&          materialProperties,
                                                        const Eigen::Matrix< double, 24, 1 >& dU,
                                                        bool                                  skewedGeometry = false )
  {
    constexpr int nElementDofs = 24;

    auto element = makeInterfaceElementWithMaterial< ElementT >( materialName,
                                                                 materialProperties.data(),
                                                                 static_cast< int >( materialProperties.size() ),
                                                                 skewedGeometry );

    std::vector< double > stateVars;
    initializeStateAndMaterial( *element, stateVars );

    std::array< double, nElementDofs >                U{};
    std::array< double, nElementDofs >                dQ{};
    std::array< double, nElementDofs >                Pe{};
    std::array< double, nElementDofs * nElementDofs > Ke{};
    std::copy( dU.data(), dU.data() + dU.size(), dQ.begin() );

    element->computeKernels( U.data(), dQ.data(), Pe.data(), Ke.data(), 0.0, 1.0 );
    return stateVars;
  }

  Eigen::Matrix< double, 24, 1 > makeHistoryIncrement( double scale )
  {
    Eigen::Matrix< double, 24, 1 > dU;
    dU.setZero();

    dU( 0 )  = scale * 0.0e-3;
    dU( 1 )  = scale * 0.0e-3;
    dU( 2 )  = scale * 0.0e-3;
    dU( 3 )  = scale * -2.0e-3;
    dU( 4 )  = scale * 0.4e-3;
    dU( 5 )  = scale * 0.2e-3;
    dU( 6 )  = scale * -1.0e-3;
    dU( 7 )  = scale * -1.7e-3;
    dU( 8 )  = scale * 0.4e-3;
    dU( 9 )  = scale * 0.8e-3;
    dU( 10 ) = scale * -0.3e-3;
    dU( 11 ) = scale * -0.1e-3;

    dU( 12 ) = scale * 0.3e-3;
    dU( 13 ) = scale * -0.2e-3;
    dU( 14 ) = scale * 0.6e-3;
    dU( 15 ) = scale * 2.2e-3;
    dU( 16 ) = scale * -0.7e-3;
    dU( 17 ) = scale * -0.3e-3;
    dU( 18 ) = scale * 1.2e-3;
    dU( 19 ) = scale * 1.9e-3;
    dU( 20 ) = scale * -0.6e-3;
    dU( 21 ) = scale * -0.9e-3;
    dU( 22 ) = scale * 0.5e-3;
    dU( 23 ) = scale * 0.2e-3;

    return dU;
  }

  Eigen::Matrix< double, 24, 24, Eigen::RowMajor > computeCentralDifferenceElementResidualJacobian(
    const std::function< ElementEvaluation( const Eigen::Matrix< double, 24, 1 >& ) >& evaluator,
    const Eigen::Matrix< double, 24, 1 >&                                              dU,
    double                                                                             relativePerturbation )
  {
    constexpr int nElementDofs = 24;

    Eigen::Matrix< double, nElementDofs, nElementDofs, Eigen::RowMajor > perturbationJacobian;
    perturbationJacobian.setZero();

    for ( int i = 0; i < nElementDofs; ++i ) {
      Eigen::Matrix< double, nElementDofs, 1 > plusDU  = dU;
      Eigen::Matrix< double, nElementDofs, 1 > minusDU = dU;
      const double perturbation                        = relativePerturbation * std::max( 1.0, std::abs( dU( i ) ) );

      plusDU( i ) += perturbation;
      minusDU( i ) -= perturbation;

      const auto plusEvaluation  = evaluator( plusDU );
      const auto minusEvaluation = evaluator( minusDU );

      perturbationJacobian.col( i ) = ( plusEvaluation.residual - minusEvaluation.residual ) / ( 2.0 * perturbation );
    }

    return perturbationJacobian;
  }

} // namespace

void TestGaussLobattoInterfaceElementTangentMatchesResidualFiniteDifference()
{
  std::cout << "\n--- TestGaussLobattoInterfaceElementTangentMatchesResidualFiniteDifference ---\n";

  using ElementT = GaussLobattoInterfaceFiniteElement< 3, 8 >;

  constexpr int nElementDofs = 24;

  Eigen::Matrix< double, nElementDofs, 1 > dU;
  dU.setZero();
  dU( 0 )  = -1.0e-5;
  dU( 1 )  = 2.0e-5;
  dU( 3 )  = 1.5e-5;
  dU( 7 )  = -1.0e-5;
  dU( 12 ) = 2.5e-5;
  dU( 13 ) = -1.5e-5;
  dU( 16 ) = 1.0e-5;
  dU( 20 ) = 2.0e-5;

  const std::vector< double > materialProperties = { 0.02, 2., 8.0e4, 0.22, 2., 1.5e5, 0.31 };

  const auto baseEvaluation = evaluateElementResponse< ElementT >( "LINEARELASTIC", materialProperties, dU );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    [&]( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateElementResponse< ElementT >( "LINEARELASTIC", materialProperties, perturbedDU );
    },
    dU,
    1e-7 );

  const auto   expectedResidualTangent = -baseEvaluation.tangent;
  const double error                   = ( perturbationResidualJacobian - expectedResidualTangent ).norm();
  const double scale                   = std::max( 1.0, expectedResidualTangent.norm() );

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(), "GL interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(), "GL interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen GL interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 1e-6,
                           "GL interface element tangent does not match residual finite difference." );
}

void TestGaussLobattoInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory()
{
  std::cout << "\n--- TestGaussLobattoInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory ---\n";

  using ElementT = GaussLobattoInterfaceFiniteElement< 3, 8 >;

  constexpr int nElementDofs = 24;

  const std::vector< double > materialProperties = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };

  const Eigen::Matrix< double, nElementDofs, 1 > historyDU = makeHistoryIncrement( 1.0 );
  const std::vector< double > historyStateVars             = makeElementStateAfterIncrement< ElementT >( "VONMISES",
                                                                                             materialProperties,
                                                                                             historyDU );

  const Eigen::Matrix< double, nElementDofs, 1 > dU = makeHistoryIncrement( 0.25 );

  const auto baseEvaluation = evaluateElementResponse< ElementT >( "VONMISES",
                                                                   materialProperties,
                                                                   dU,
                                                                   &historyStateVars );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    [&]( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateElementResponse< ElementT >( "VONMISES", materialProperties, perturbedDU, &historyStateVars );
    },
    dU,
    1e-7 );

  const auto   expectedResidualTangent = -baseEvaluation.tangent;
  const double error                   = ( perturbationResidualJacobian - expectedResidualTangent ).norm();
  const double scale                   = std::max( 1.0, expectedResidualTangent.norm() );
  std::cout << "plastic flat GLIQUAD4 relative tangent error = " << error / scale << "\n";

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(), "Plastic GL interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(), "Plastic GL interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen plastic GL interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 2e-3,
                           "Plastic GL interface element tangent does not match residual finite difference after "
                           "committed history." );
}

void TestSkewedGaussLobattoInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory()
{
  std::cout << "\n--- TestSkewedGaussLobattoInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory "
               "---\n";

  using ElementT = GaussLobattoInterfaceFiniteElement< 3, 8 >;

  constexpr int nElementDofs = 24;

  const std::vector< double > materialProperties = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };

  const Eigen::Matrix< double, nElementDofs, 1 > historyDU = makeHistoryIncrement( 1.0 );
  const std::vector< double > historyStateVars             = makeElementStateAfterIncrement< ElementT >( "VONMISES",
                                                                                             materialProperties,
                                                                                             historyDU,
                                                                                             true );

  const Eigen::Matrix< double, nElementDofs, 1 > dU = makeHistoryIncrement( 0.25 );

  const auto baseEvaluation = evaluateElementResponse< ElementT >( "VONMISES",
                                                                   materialProperties,
                                                                   dU,
                                                                   &historyStateVars,
                                                                   true );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    [&]( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateElementResponse< ElementT >( "VONMISES",
                                                  materialProperties,
                                                  perturbedDU,
                                                  &historyStateVars,
                                                  true );
    },
    dU,
    1e-7 );

  const auto   expectedResidualTangent = -baseEvaluation.tangent;
  const double error                   = ( perturbationResidualJacobian - expectedResidualTangent ).norm();
  const double scale                   = std::max( 1.0, expectedResidualTangent.norm() );
  std::cout << "plastic skewed GLIQUAD4 relative tangent error = " << error / scale << "\n";

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(),
                           "Skewed plastic GL interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(),
                           "Skewed plastic GL interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen skewed plastic GL interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 2e-3,
                           "Skewed plastic GL interface element tangent does not match residual finite difference "
                           "after committed history." );
}

/**
 * Residual parity with YIQUAD4 under UNIFORM through-thickness loading: when
 * A+ == A- at every quadrature point (achieved here by an identical top and
 * bottom in-plane DOF pattern, offset only by a small rigid normal-direction
 * gap), all five Gauss-Lobatto stations see the SAME tangential gradient and
 * the SAME normal-gradient-compatibility target, so the whole formulation
 * degenerates to the ordinary (uniform through-thickness) response -- which
 * must match YIQUAD4's (the validated two-station equilibrated element) to
 * tight tolerance. The TANGENT is not compared (analogous to the X-vs-
 * Corrected uniform-loading test): probing bottom/top-asymmetric DOF
 * directions exercises the bending response that is the entire point of the
 * multi-station formulation and the two elements are expected to respond
 * differently there.
 */
void TestGaussLobattoInterfaceElementUniformLoadingResidualMatchesYElement()
{
  std::cout << "\n--- TestGaussLobattoInterfaceElementUniformLoadingResidualMatchesYElement ---\n";

  using GLElementT = GaussLobattoInterfaceFiniteElement< 3, 8 >;
  using YElementT  = YInterfaceFiniteElement< 3, 8 >;

  constexpr int nElementDofs = 24;
  constexpr int halfNDof     = 12;

  const std::array< double, halfNDof > bottomPattern =
    { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4, -0.5e-4, 0.8e-4, 1.0e-4, 0.2e-4, -0.3e-4, 0.4e-4 };
  const std::array< double, 3 > topOffset = { 2.0e-5, -1.0e-5, 3.0e-5 };

  Eigen::Matrix< double, nElementDofs, 1 > dU;
  for ( int i = 0; i < halfNDof; ++i ) {
    dU( i )            = bottomPattern[i];
    dU( halfNDof + i ) = bottomPattern[i] + topOffset[i % 3];
  }

  const std::vector< double > materialProperties = { 1e5, 0.3, 0.01 };

  const auto glEvaluation = evaluateElementResponse< GLElementT >( "LINEARELASTIC", materialProperties, dU );
  const auto yEvaluation  = evaluateElementResponse< YElementT >( "LINEARELASTIC", materialProperties, dU );

  const double error = ( glEvaluation.residual - yEvaluation.residual ).norm();
  const double scale = std::max( 1.0, yEvaluation.residual.norm() );

  std::cout << "uniform-loading GL vs Y residual relative error = " << error / scale << "\n";

  throwExceptionOnFailure( error / scale < 1e-8,
                           "GaussLobattoInterfaceFiniteElement residual under uniform through-thickness loading "
                           "does not match YInterfaceFiniteElement." );
}

int main()
{
  auto tests = std::vector< std::function<
    void() > >{ TestGaussLobattoInterfaceElementTangentMatchesResidualFiniteDifference,
                TestGaussLobattoInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory,
                TestSkewedGaussLobattoInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory,
                TestGaussLobattoInterfaceElementUniformLoadingResidualMatchesYElement };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
