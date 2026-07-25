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

  std::unique_ptr< YInterfaceFiniteElement< 3, 8 > > makeYInterfaceElementWithMaterial(
    const std::string& materialName,
    const double*      materialProperties,
    int                nMaterialProperties,
    bool               skewedGeometry = false )
  {
    constexpr int nDim   = 3;
    constexpr int nNodes = 8;

    const int  elId    = 7;
    const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
    const auto secType = YInterfaceFiniteElement< nDim, nNodes >::SectionType::Interface;

    auto element = std::make_unique< YInterfaceFiniteElement< nDim, nNodes > >( elId, intType, secType );
    element->assignNodeCoordinates( skewedGeometry ? skewedCoordinates.data() : flatCoordinates.data() );

    static std::array< double, 1 > elPropsVec = { 1.0 };
    ElementProperties              elProps( elPropsVec.data(), static_cast< int >( elPropsVec.size() ) );
    element->assignProperty( elProps );

    element->assignMaterial( materialName, materialProperties, nMaterialProperties );
    return element;
  }

  template < int nDim, int nNodes >
  void initializeStateAndMaterial( YInterfaceFiniteElement< nDim, nNodes >& element, std::vector< double >& stateVars )
  {
    stateVars.assign( element.getNumberOfRequiredStateVars(), 0.0 );
    element.assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element.initializeYourself();
    element.setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
  }

  struct YElementEvaluation {
    Eigen::Matrix< double, 24, 1 >                   residual;
    Eigen::Matrix< double, 24, 24, Eigen::RowMajor > tangent;
  };

  YElementEvaluation evaluateYElementResponse( const std::string&                    materialName,
                                               const std::vector< double >&          materialProperties,
                                               const Eigen::Matrix< double, 24, 1 >& dU,
                                               const std::vector< double >*          initialStateVars = nullptr,
                                               bool                                  skewedGeometry   = false )
  {
    constexpr int nElementDofs = 24;

    auto element = makeYInterfaceElementWithMaterial( materialName,
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

    YElementEvaluation evaluation;
    evaluation.residual = Eigen::Map< Eigen::Matrix< double, nElementDofs, 1 > >( Pe.data() );
    evaluation.tangent  = Eigen::Map< Eigen::Matrix< double, nElementDofs, nElementDofs, Eigen::RowMajor > >(
      Ke.data() );
    return evaluation;
  }

  std::vector< double > makeYElementStateAfterIncrement( const std::string&                    materialName,
                                                         const std::vector< double >&          materialProperties,
                                                         const Eigen::Matrix< double, 24, 1 >& dU,
                                                         bool                                  skewedGeometry = false )
  {
    constexpr int nElementDofs = 24;

    auto element = makeYInterfaceElementWithMaterial( materialName,
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

  Eigen::Matrix< double, 24, 1 > makeYHistoryIncrement( double scale )
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
    const std::function< YElementEvaluation( const Eigen::Matrix< double, 24, 1 >& ) >& evaluator,
    const Eigen::Matrix< double, 24, 1 >&                                               dU,
    double                                                                              relativePerturbation )
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

void TestYInterfaceElementTangentMatchesResidualFiniteDifference()
{
  std::cout << "\n--- TestYInterfaceElementTangentMatchesResidualFiniteDifference ---\n";

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

  const auto baseEvaluation = evaluateYElementResponse( "LINEARELASTIC", materialProperties, dU );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    [&]( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateYElementResponse( "LINEARELASTIC", materialProperties, perturbedDU );
    },
    dU,
    1e-7 );

  const auto   expectedResidualTangent = -baseEvaluation.tangent;
  const double error                   = ( perturbationResidualJacobian - expectedResidualTangent ).norm();
  const double scale                   = std::max( 1.0, expectedResidualTangent.norm() );

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(), "Y interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(), "Y interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen Y interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 1e-6,
                           "Y interface element tangent does not match residual finite difference." );
}

void TestYInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory()
{
  std::cout << "\n--- TestYInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory ---\n";

  constexpr int nElementDofs = 24;

  const std::vector< double > materialProperties = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };

  const Eigen::Matrix< double, nElementDofs, 1 > historyDU        = makeYHistoryIncrement( 1.0 );
  const std::vector< double >                    historyStateVars = makeYElementStateAfterIncrement( "VONMISES",
                                                                                  materialProperties,
                                                                                  historyDU );

  const Eigen::Matrix< double, nElementDofs, 1 > dU = makeYHistoryIncrement( 0.25 );

  const auto baseEvaluation = evaluateYElementResponse( "VONMISES", materialProperties, dU, &historyStateVars );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    [&]( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateYElementResponse( "VONMISES", materialProperties, perturbedDU, &historyStateVars );
    },
    dU,
    1e-7 );

  const auto   expectedResidualTangent = -baseEvaluation.tangent;
  const double error                   = ( perturbationResidualJacobian - expectedResidualTangent ).norm();
  const double scale                   = std::max( 1.0, expectedResidualTangent.norm() );
  std::cout << "plastic flat YIQUAD relative tangent error = " << error / scale << "\n";

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(), "Plastic Y interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(), "Plastic Y interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen plastic Y interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 2e-3,
                           "Plastic Y interface element tangent does not match residual finite difference after "
                           "committed history." );
}

void TestSkewedYInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory()
{
  std::cout << "\n--- TestSkewedYInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory ---\n";

  constexpr int nElementDofs = 24;

  const std::vector< double > materialProperties = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };

  const Eigen::Matrix< double, nElementDofs, 1 > historyDU        = makeYHistoryIncrement( 1.0 );
  const std::vector< double >                    historyStateVars = makeYElementStateAfterIncrement( "VONMISES",
                                                                                  materialProperties,
                                                                                  historyDU,
                                                                                  true );

  const Eigen::Matrix< double, nElementDofs, 1 > dU = makeYHistoryIncrement( 0.25 );

  const auto baseEvaluation = evaluateYElementResponse( "VONMISES", materialProperties, dU, &historyStateVars, true );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    [&]( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateYElementResponse( "VONMISES", materialProperties, perturbedDU, &historyStateVars, true );
    },
    dU,
    1e-7 );

  const auto   expectedResidualTangent = -baseEvaluation.tangent;
  const double error                   = ( perturbationResidualJacobian - expectedResidualTangent ).norm();
  const double scale                   = std::max( 1.0, expectedResidualTangent.norm() );
  std::cout << "plastic skewed YIQUAD relative tangent error = " << error / scale << "\n";

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(),
                           "Skewed plastic Y interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(),
                           "Skewed plastic Y interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen skewed plastic Y interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 2e-3,
                           "Skewed plastic Y interface element tangent does not match residual finite difference "
                           "after committed history." );
}

int main()
{
  auto tests = std::vector<
    std::function< void() > >{ TestYInterfaceElementTangentMatchesResidualFiniteDifference,
                               TestYInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory,
                               TestSkewedYInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
