#include "Marmot/CorrectedInterfaceFiniteElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/XInterfaceFiniteElement.h"

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

  template < typename DerivedA, typename DerivedB >
  void assertMatrixNear( const Eigen::MatrixBase< DerivedA >& actual,
                         const Eigen::MatrixBase< DerivedB >& expected,
                         double                               tol,
                         const std::string&                   message )
  {
    throwExceptionOnFailure( actual.rows() == expected.rows() && actual.cols() == expected.cols(),
                             message + ": matrix shape mismatch." );

    const double err = ( actual - expected ).template lpNorm< Eigen::Infinity >();
    throwExceptionOnFailure( err < tol, message + ": max error = " + std::to_string( err ) );
  }

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

  std::unique_ptr< XInterfaceFiniteElement< 3, 8 > > makeXInterfaceElementWithMaterial(
    const std::string& materialName,
    const double*      materialProperties,
    int                nMaterialProperties,
    bool               skewedGeometry = false )
  {
    constexpr int nDim   = 3;
    constexpr int nNodes = 8;

    const int  elId    = 5;
    const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
    const auto secType = XInterfaceFiniteElement< nDim, nNodes >::SectionType::Interface;

    auto element = std::make_unique< XInterfaceFiniteElement< nDim, nNodes > >( elId, intType, secType );
    element->assignNodeCoordinates( skewedGeometry ? skewedCoordinates.data() : flatCoordinates.data() );

    static std::array< double, 1 > elPropsVec = { 1.0 };
    ElementProperties              elProps( elPropsVec.data(), static_cast< int >( elPropsVec.size() ) );
    element->assignProperty( elProps );

    element->assignMaterial( materialName, materialProperties, nMaterialProperties );
    return element;
  }

  std::unique_ptr< CorrectedInterfaceFiniteElement< 3, 8 > > makeCorrectedInterfaceElementWithMaterial(
    const std::string& materialName,
    const double*      materialProperties,
    int                nMaterialProperties,
    bool               skewedGeometry = false )
  {
    constexpr int nDim   = 3;
    constexpr int nNodes = 8;

    const int  elId    = 6;
    const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
    const auto secType = CorrectedInterfaceFiniteElement< nDim, nNodes >::SectionType::Interface;

    auto element = std::make_unique< CorrectedInterfaceFiniteElement< nDim, nNodes > >( elId, intType, secType );
    element->assignNodeCoordinates( skewedGeometry ? skewedCoordinates.data() : flatCoordinates.data() );

    static std::array< double, 1 > elPropsVec = { 1.0 };
    ElementProperties              elProps( elPropsVec.data(), static_cast< int >( elPropsVec.size() ) );
    element->assignProperty( elProps );

    element->assignMaterial( materialName, materialProperties, nMaterialProperties );
    return element;
  }

  template < int nDim, int nNodes >
  void initializeStateAndMaterial( XInterfaceFiniteElement< nDim, nNodes >& element, std::vector< double >& stateVars )
  {
    stateVars.assign( element.getNumberOfRequiredStateVars(), 0.0 );
    element.assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element.initializeYourself();
    element.setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
  }

  template < int nDim, int nNodes >
  void initializeStateAndMaterial( CorrectedInterfaceFiniteElement< nDim, nNodes >& element,
                                   std::vector< double >&                           stateVars )
  {
    stateVars.assign( element.getNumberOfRequiredStateVars(), 0.0 );
    element.assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element.initializeYourself();
    element.setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
  }

  struct XElementEvaluation {
    Eigen::Matrix< double, 24, 1 >                   residual;
    Eigen::Matrix< double, 24, 24, Eigen::RowMajor > tangent;
  };

  XElementEvaluation evaluateXElementResponse( const std::string&                    materialName,
                                               const std::vector< double >&          materialProperties,
                                               const Eigen::Matrix< double, 24, 1 >& dU,
                                               const std::vector< double >*          initialStateVars = nullptr,
                                               bool                                  skewedGeometry   = false )
  {
    constexpr int nElementDofs = 24;

    auto element = makeXInterfaceElementWithMaterial( materialName,
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

    XElementEvaluation evaluation;
    evaluation.residual = Eigen::Map< Eigen::Matrix< double, nElementDofs, 1 > >( Pe.data() );
    evaluation.tangent  = Eigen::Map< Eigen::Matrix< double, nElementDofs, nElementDofs, Eigen::RowMajor > >(
      Ke.data() );
    return evaluation;
  }

  std::vector< double > makeXElementStateAfterIncrement( const std::string&                    materialName,
                                                         const std::vector< double >&          materialProperties,
                                                         const Eigen::Matrix< double, 24, 1 >& dU,
                                                         bool                                  skewedGeometry = false )
  {
    constexpr int nElementDofs = 24;

    auto element = makeXInterfaceElementWithMaterial( materialName,
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

  Eigen::Matrix< double, 24, 1 > makeXHistoryIncrement( double scale )
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
    const std::function< XElementEvaluation( const Eigen::Matrix< double, 24, 1 >& ) >& evaluator,
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

void TestXInterfaceElementTangentMatchesResidualFiniteDifference()
{
  std::cout << "\n--- TestXInterfaceElementTangentMatchesResidualFiniteDifference ---\n";

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

  const auto baseEvaluation = evaluateXElementResponse( "LINEARELASTIC", materialProperties, dU );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    [&]( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateXElementResponse( "LINEARELASTIC", materialProperties, perturbedDU );
    },
    dU,
    1e-7 );

  const auto   expectedResidualTangent = -baseEvaluation.tangent;
  const double error                   = ( perturbationResidualJacobian - expectedResidualTangent ).norm();
  const double scale                   = std::max( 1.0, expectedResidualTangent.norm() );

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(), "X interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(), "X interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen X interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 1e-6,
                           "X interface element tangent does not match residual finite difference." );
}

void TestXInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory()
{
  std::cout << "\n--- TestXInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory ---\n";

  constexpr int nElementDofs = 24;

  const std::vector< double > materialProperties = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };

  const Eigen::Matrix< double, nElementDofs, 1 > historyDU        = makeXHistoryIncrement( 1.0 );
  const std::vector< double >                    historyStateVars = makeXElementStateAfterIncrement( "VONMISES",
                                                                                  materialProperties,
                                                                                  historyDU );

  const Eigen::Matrix< double, nElementDofs, 1 > dU = makeXHistoryIncrement( 0.25 );

  const auto baseEvaluation = evaluateXElementResponse( "VONMISES", materialProperties, dU, &historyStateVars );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    [&]( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateXElementResponse( "VONMISES", materialProperties, perturbedDU, &historyStateVars );
    },
    dU,
    1e-7 );

  const auto   expectedResidualTangent = -baseEvaluation.tangent;
  const double error                   = ( perturbationResidualJacobian - expectedResidualTangent ).norm();
  const double scale                   = std::max( 1.0, expectedResidualTangent.norm() );
  std::cout << "plastic flat XIQUAD relative tangent error = " << error / scale << "\n";

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(), "Plastic X interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(), "Plastic X interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen plastic X interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 2e-3,
                           "Plastic X interface element tangent does not match residual finite difference after "
                           "committed history." );
}

void TestSkewedXInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory()
{
  std::cout << "\n--- TestSkewedXInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory ---\n";

  constexpr int nElementDofs = 24;

  const std::vector< double > materialProperties = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };

  const Eigen::Matrix< double, nElementDofs, 1 > historyDU        = makeXHistoryIncrement( 1.0 );
  const std::vector< double >                    historyStateVars = makeXElementStateAfterIncrement( "VONMISES",
                                                                                  materialProperties,
                                                                                  historyDU,
                                                                                  true );

  const Eigen::Matrix< double, nElementDofs, 1 > dU = makeXHistoryIncrement( 0.25 );

  const auto baseEvaluation = evaluateXElementResponse( "VONMISES", materialProperties, dU, &historyStateVars, true );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    [&]( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateXElementResponse( "VONMISES", materialProperties, perturbedDU, &historyStateVars, true );
    },
    dU,
    1e-7 );

  const auto   expectedResidualTangent = -baseEvaluation.tangent;
  const double error                   = ( perturbationResidualJacobian - expectedResidualTangent ).norm();
  const double scale                   = std::max( 1.0, expectedResidualTangent.norm() );
  std::cout << "plastic skewed XIQUAD relative tangent error = " << error / scale << "\n";

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(),
                           "Skewed plastic X interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(),
                           "Skewed plastic X interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen skewed plastic X interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 2e-3,
                           "Skewed plastic X interface element tangent does not match residual finite difference "
                           "after committed history." );
}

/**
 * Sanity/parity check (element level): if the SAME in-plane nodal pattern is
 * applied to both the top and bottom side nodes (differing only by a
 * rigid, uniform offset on the top side, which BmatSide annihilates -- see
 * TestSingleInputFileElementGeometryMatrices in the plain-interface test
 * suite), there is no through-thickness asymmetry: A+ == A- at every
 * quadrature point. In that case the X element's RESIDUAL at this state must
 * reduce EXACTLY to the Corrected element's residual (their generalized
 * states differ only by the h/2 vs h integration split, which cancels once
 * summed back through BmatAverage = 0.5*(BPlusFull + BMinusFull)).
 *
 * The TANGENT is intentionally NOT compared here: Ke probes sensitivity to
 * ALL 24 dof directions, including bottom/top-asymmetric perturbations that
 * break A+ == A-. Along exactly those directions the two elements are
 * supposed to respond differently -- giving the layer stiffness against the
 * rotation/bending mode is the entire point of the X formulation, so a
 * divergent tangent at this state is the expected, correct outcome, not a
 * defect.
 */
void TestXInterfaceElementUniformLoadingResidualMatchesCorrectedElement()
{
  std::cout << "\n--- TestXInterfaceElementUniformLoadingResidualMatchesCorrectedElement ---\n";

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

  auto xElement         = makeXInterfaceElementWithMaterial( "LINEARELASTIC",
                                                     materialProperties.data(),
                                                     static_cast< int >( materialProperties.size() ) );
  auto correctedElement = makeCorrectedInterfaceElementWithMaterial( "LINEARELASTIC",
                                                                     materialProperties.data(),
                                                                     static_cast< int >( materialProperties.size() ) );

  std::vector< double > xStateVars, correctedStateVars;
  initializeStateAndMaterial( *xElement, xStateVars );
  initializeStateAndMaterial( *correctedElement, correctedStateVars );

  std::array< double, nElementDofs >                U{};
  std::array< double, nElementDofs >                dQ{};
  std::array< double, nElementDofs >                xPe{}, correctedPe{};
  std::array< double, nElementDofs * nElementDofs > xKe{}, correctedKe{};
  std::copy( dU.data(), dU.data() + dU.size(), dQ.begin() );

  xElement->computeKernels( U.data(), dQ.data(), xPe.data(), xKe.data(), 0.0, 1.0 );
  correctedElement->computeKernels( U.data(), dQ.data(), correctedPe.data(), correctedKe.data(), 0.0, 1.0 );

  Eigen::Map< Eigen::Matrix< double, nElementDofs, 1 > > xPeMap( xPe.data() );
  Eigen::Map< Eigen::Matrix< double, nElementDofs, 1 > > correctedPeMap( correctedPe.data() );
  (void)xKe;
  (void)correctedKe;

  throwExceptionOnFailure( xPeMap.norm() > 0.0, "Uniform-loading X element residual should be nonzero." );

  assertMatrixNear( xPeMap,
                    correctedPeMap,
                    1e-10,
                    "Uniform-loading (A+ == A- at every QP) X element residual differs from Corrected element." );

  for ( size_t q = 0; q < xElement->qps.size(); ++q ) {
    const Eigen::Vector3d combinedForce = xElement->qps[q].managedStateVars->forcePlus +
                                          xElement->qps[q].managedStateVars->forceMinus;
    const Eigen::Matrix< double, 9, 1 > combinedSurfaceStress = xElement->qps[q].managedStateVars->surfaceStressPlus +
                                                                xElement->qps[q].managedStateVars->surfaceStressMinus;

    assertMatrixNear( combinedForce,
                      correctedElement->qps[q].managedStateVars->force,
                      1e-10,
                      "Uniform-loading combined force differs from Corrected element at qp " + std::to_string( q ) );
    assertMatrixNear( combinedSurfaceStress,
                      correctedElement->qps[q].managedStateVars->surfaceStress,
                      1e-10,
                      "Uniform-loading combined surface stress differs from Corrected element at qp " +
                        std::to_string( q ) );
  }
}

int main()
{
  auto tests = std::vector<
    std::function< void() > >{ TestXInterfaceElementTangentMatchesResidualFiniteDifference,
                               TestXInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory,
                               TestSkewedXInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory,
                               TestXInterfaceElementUniformLoadingResidualMatchesCorrectedElement };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
