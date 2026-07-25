#include "Marmot/CorrectedInterfaceFiniteElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"

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

  std::unique_ptr< CorrectedInterfaceFiniteElement< 3, 8 > > makeCorrectedInterfaceElementWithMaterial(
    const std::string& materialName,
    const double*      materialProperties,
    int                nMaterialProperties,
    bool               skewedGeometry = false )
  {
    constexpr int nDim   = 3;
    constexpr int nNodes = 8;

    const int  elId    = 4;
    const auto intType = FiniteElement::Quadrature::IntegrationTypes::FullIntegration;
    const auto secType = CorrectedInterfaceFiniteElement< nDim, nNodes >::SectionType::Interface;

    auto element = std::make_unique< CorrectedInterfaceFiniteElement< nDim, nNodes > >( elId, intType, secType );

    static std::array< double, nDim* nNodes > flatCoordinates = {
      -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0,
      -0.5, -0.5, 0.1, 0.5, -0.5, 0.1, 0.5, 0.5, 0.1, -0.5, 0.5, 0.1,
    };

    static std::array< double, nDim* nNodes > skewedCoordinates = {
      -0.500000, -0.500000, 0.088163,  0.500000,  -0.500000, -0.088163,
      0.500000,  0.500000,  -0.088163, -0.500000, 0.500000,  0.088163,

      -0.500000, -0.500000, 0.188163,  0.500000,  -0.500000, 0.011837,
      0.500000,  0.500000,  0.011837,  -0.500000, 0.500000,  0.188163,
    };

    element->assignNodeCoordinates( skewedGeometry ? skewedCoordinates.data() : flatCoordinates.data() );

    static std::array< double, 1 > elPropsVec = { 1.0 };
    ElementProperties              elProps( elPropsVec.data(), static_cast< int >( elPropsVec.size() ) );
    element->assignProperty( elProps );

    element->assignMaterial( materialName, materialProperties, nMaterialProperties );
    return element;
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

  struct CorrectedElementEvaluation {
    Eigen::Matrix< double, 24, 1 >                   residual;
    Eigen::Matrix< double, 24, 24, Eigen::RowMajor > tangent;
  };

  CorrectedElementEvaluation evaluateCorrectedElementResponse( const std::string&                    materialName,
                                                               const std::vector< double >&          materialProperties,
                                                               const Eigen::Matrix< double, 24, 1 >& dU,
                                                               const std::vector< double >* initialStateVars = nullptr,
                                                               bool                         skewedGeometry   = false )
  {
    constexpr int nElementDofs = 24;

    auto element = makeCorrectedInterfaceElementWithMaterial( materialName,
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

    CorrectedElementEvaluation evaluation;
    evaluation.residual = Eigen::Map< Eigen::Matrix< double, nElementDofs, 1 > >( Pe.data() );
    evaluation.tangent  = Eigen::Map< Eigen::Matrix< double, nElementDofs, nElementDofs, Eigen::RowMajor > >(
      Ke.data() );
    return evaluation;
  }

  CorrectedElementEvaluation evaluateCorrectedElementResponse( const Eigen::Matrix< double, 24, 1 >& dU )
  {
    const std::vector< double > materialProperties = { 0.02, 2., 8.0e4, 0.22, 2., 1.5e5, 0.31 };
    return evaluateCorrectedElementResponse( "LINEARELASTIC", materialProperties, dU );
  }

  std::vector< double > makeCorrectedElementStateAfterIncrement( const std::string&           materialName,
                                                                 const std::vector< double >& materialProperties,
                                                                 const Eigen::Matrix< double, 24, 1 >& dU,
                                                                 bool skewedGeometry = false )
  {
    constexpr int nElementDofs = 24;

    auto element = makeCorrectedInterfaceElementWithMaterial( materialName,
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

  Eigen::Matrix< double, 24, 1 > makeCorrectedHistoryIncrement( double scale )
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
    const std::function< CorrectedElementEvaluation( const Eigen::Matrix< double, 24, 1 >& ) >& evaluator,
    const Eigen::Matrix< double, 24, 1 >&                                                       dU,
    double                                                                                      relativePerturbation )
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

void TestCorrectedInterfaceElementTangentMatchesResidualFiniteDifference()
{
  std::cout << "\n--- TestCorrectedInterfaceElementTangentMatchesResidualFiniteDifference ---\n";

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

  const auto baseEvaluation = evaluateCorrectedElementResponse( dU );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    []( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateCorrectedElementResponse( perturbedDU );
    },
    dU,
    1e-7 );

  const auto   expectedResidualTangent = -baseEvaluation.tangent;
  const double error                   = ( perturbationResidualJacobian - expectedResidualTangent ).norm();
  const double scale                   = std::max( 1.0, expectedResidualTangent.norm() );

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(), "Corrected interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(), "Corrected interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen corrected interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 1e-6,
                           "Corrected interface element tangent does not match residual finite difference." );
}

void TestCorrectedInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory()
{
  std::cout << "\n--- TestCorrectedInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory ---\n";

  constexpr int nElementDofs = 24;

  const std::vector< double > materialProperties = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };

  const Eigen::Matrix< double, nElementDofs, 1 > historyDU        = makeCorrectedHistoryIncrement( 1.0 );
  const std::vector< double >                    historyStateVars = makeCorrectedElementStateAfterIncrement( "VONMISES",
                                                                                          materialProperties,
                                                                                          historyDU );

  const Eigen::Matrix< double, nElementDofs, 1 > dU = makeCorrectedHistoryIncrement( 0.25 );

  const auto baseEvaluation = evaluateCorrectedElementResponse( "VONMISES", materialProperties, dU, &historyStateVars );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    [&]( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateCorrectedElementResponse( "VONMISES", materialProperties, perturbedDU, &historyStateVars );
    },
    dU,
    1e-7 );

  const auto   expectedResidualTangent = -baseEvaluation.tangent;
  const double error                   = ( perturbationResidualJacobian - expectedResidualTangent ).norm();
  const double scale                   = std::max( 1.0, expectedResidualTangent.norm() );
  std::cout << "plastic flat EIQUAD relative tangent error = " << error / scale << "\n";

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(),
                           "Plastic corrected interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(),
                           "Plastic corrected interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen plastic corrected interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 2e-3,
                           "Plastic corrected interface element tangent does not match residual finite difference "
                           "after committed history." );
}

void TestSkewedCorrectedInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory()
{
  std::cout << "\n--- TestSkewedCorrectedInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory "
               "---\n";

  constexpr int nElementDofs = 24;

  const std::vector< double > materialProperties = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };

  const Eigen::Matrix< double, nElementDofs, 1 > historyDU        = makeCorrectedHistoryIncrement( 1.0 );
  const std::vector< double >                    historyStateVars = makeCorrectedElementStateAfterIncrement( "VONMISES",
                                                                                          materialProperties,
                                                                                          historyDU,
                                                                                          true );

  const Eigen::Matrix< double, nElementDofs, 1 > dU = makeCorrectedHistoryIncrement( 0.25 );

  const auto baseEvaluation = evaluateCorrectedElementResponse( "VONMISES",
                                                                materialProperties,
                                                                dU,
                                                                &historyStateVars,
                                                                true );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    [&]( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateCorrectedElementResponse( "VONMISES", materialProperties, perturbedDU, &historyStateVars, true );
    },
    dU,
    1e-7 );

  const auto   expectedResidualTangent = -baseEvaluation.tangent;
  const double error                   = ( perturbationResidualJacobian - expectedResidualTangent ).norm();
  const double scale                   = std::max( 1.0, expectedResidualTangent.norm() );
  std::cout << "plastic skewed EIQUAD relative tangent error = " << error / scale << "\n";

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(),
                           "Skewed plastic corrected interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(),
                           "Skewed plastic corrected interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen skewed plastic corrected interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 2e-3,
                           "Skewed plastic corrected interface element tangent does not match residual finite "
                           "difference after committed history." );
}

void TestSkewedCorrectedInterfaceElementWiechertTangentMatchesResidualFiniteDifferenceAfterHistory()
{
  std::cout << "\n--- TestSkewedCorrectedInterfaceElementWiechertTangentMatchesResidualFiniteDifferenceAfterHistory "
               "---\n";

  constexpr int nElementDofs = 24;

  const std::vector< double > materialProperties = { 2e5, 0.2, 0.2, 0.5, 0.1, 10., 1e-4, 1., 2400. };

  const Eigen::Matrix< double, nElementDofs, 1 > historyDU = makeCorrectedHistoryIncrement( 0.6 );
  const std::vector< double > historyStateVars = makeCorrectedElementStateAfterIncrement( "LINEARVISCOELASTICWIECHERT",
                                                                                          materialProperties,
                                                                                          historyDU,
                                                                                          true );

  const Eigen::Matrix< double, nElementDofs, 1 > dU = makeCorrectedHistoryIncrement( 0.15 );

  const auto baseEvaluation = evaluateCorrectedElementResponse( "LINEARVISCOELASTICWIECHERT",
                                                                materialProperties,
                                                                dU,
                                                                &historyStateVars,
                                                                true );

  const auto perturbationResidualJacobian = computeCentralDifferenceElementResidualJacobian(
    [&]( const Eigen::Matrix< double, nElementDofs, 1 >& perturbedDU ) {
      return evaluateCorrectedElementResponse( "LINEARVISCOELASTICWIECHERT",
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
  std::cout << "Wiechert skewed EIQUAD relative tangent error = " << error / scale << "\n";

  throwExceptionOnFailure( baseEvaluation.residual.allFinite(),
                           "Skewed Wiechert corrected interface residual contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.tangent.allFinite(),
                           "Skewed Wiechert corrected interface tangent contains nan or inf." );
  throwExceptionOnFailure( baseEvaluation.residual.norm() > 0.0,
                           "Chosen skewed Wiechert corrected interface increment should produce a nonzero residual." );
  throwExceptionOnFailure( error / scale < 1e-6,
                           "Skewed Wiechert corrected interface element tangent does not match residual finite "
                           "difference after committed history." );
}

// End-to-end natural-failure cutback test for the *real* corrected-interface
// material (as opposed to a test double). A Von Mises base material with a
// tiny yield stress and a steep exponential-softening branch (deltaFy < 0
// with a very large delta) makes the scalar return-mapping Newton in
// VonMisesModel::computeStress overshoot and hit its iteration cap for the
// very first plastically-loading quadrature point, so
// CorrectedInterfaceFiniteElement::computeKernels genuinely raises
// Marmot::StressUpdateFailed here. This verifies (a) the exception really
// propagates out of computeKernels, (b) the persisted quadrature-point state
// is left byte-identical by the failed attempt, and (c)
// MarmotElement::computeYourself still translates it into a pNewDT<1
// cutback request instead of rethrowing.
void TestCorrectedInterfaceElementNaturalStressUpdateFailurePreservesStateAndRequestsCutback()
{
  std::cout << "\n--- TestCorrectedInterfaceElementNaturalStressUpdateFailurePreservesStateAndRequestsCutback ---\n";

  constexpr int               nElementDofs       = 24;
  const std::vector< double > materialProperties = { 210000., 0.3, 1e-8, 1e-3, 0., -20., 1e5, 2400. };

  auto                  element = makeCorrectedInterfaceElementWithMaterial( "VONMISES",
                                                            materialProperties.data(),
                                                            static_cast< int >( materialProperties.size() ) );
  std::vector< double > stateVars;
  initializeStateAndMaterial( *element, stateVars );
  const std::vector< double > stateBeforeFailedAttempt = stateVars;

  const Eigen::Matrix< double, nElementDofs, 1 > dU = makeCorrectedHistoryIncrement( 1.0 );

  std::array< double, nElementDofs >                U{};
  std::array< double, nElementDofs >                dQ{};
  std::array< double, nElementDofs >                Pe{};
  std::array< double, nElementDofs * nElementDofs > Ke{};
  std::copy( dU.data(), dU.data() + dU.size(), dQ.begin() );

  bool stressUpdateFailed = false;
  try {
    element->computeKernels( U.data(), dQ.data(), Pe.data(), Ke.data(), 0.0, 1.0 );
  }
  catch ( const Marmot::StressUpdateFailed& ) {
    stressUpdateFailed = true;
  }

  throwExceptionOnFailure( stressUpdateFailed,
                           "Perfectly-plastic single-material corrected interface element did not raise "
                           "StressUpdateFailed for a nonzero increment as expected." );

  double stateMaxAbsDiff = 0.0;
  for ( size_t i = 0; i < stateVars.size(); ++i )
    stateMaxAbsDiff = std::max( stateMaxAbsDiff, std::abs( stateVars[i] - stateBeforeFailedAttempt[i] ) );

  throwExceptionOnFailure( stateMaxAbsDiff == 0.0,
                           "Corrected interface element's persisted quadrature-point state was mutated by a "
                           "failed (StressUpdateFailed) increment attempt; a rejected trial must leave the "
                           "committed state untouched so a time-step cutback can safely retry." );

  double                  pNewDT               = 1e36;
  std::array< double, 2 > time                 = { 0.0, 0.0 };
  bool                    computeYourselfThrew = false;
  try {
    element->computeYourself( U.data(), dQ.data(), Pe.data(), Ke.data(), time.data(), 1.0, pNewDT );
  }
  catch ( const std::exception& ) {
    computeYourselfThrew = true;
  }

  throwExceptionOnFailure( !computeYourselfThrew,
                           "MarmotElement::computeYourself should translate the corrected-interface material's "
                           "StressUpdateFailed into pNewDT rather than rethrowing." );
  throwExceptionOnFailure( pNewDT < 1.0,
                           "MarmotElement::computeYourself did not request a cutback after the corrected-interface "
                           "material's StressUpdateFailed." );
}

int main()
{
  auto tests = std::vector< std::function<
    void() > >{ TestCorrectedInterfaceElementNaturalStressUpdateFailurePreservesStateAndRequestsCutback,
                TestCorrectedInterfaceElementTangentMatchesResidualFiniteDifference,
                TestCorrectedInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory,
                TestSkewedCorrectedInterfaceElementPlasticTangentMatchesResidualFiniteDifferenceAfterHistory,
                TestSkewedCorrectedInterfaceElementWiechertTangentMatchesResidualFiniteDifferenceAfterHistory };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
