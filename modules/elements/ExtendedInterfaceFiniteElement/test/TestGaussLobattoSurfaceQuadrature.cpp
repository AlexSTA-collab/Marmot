/* Controlled test of the SURFACE quadrature option of GLIQUAD4.
 *
 * The element integrates over two nested domains:
 *   * the 2D parametric SURFACE (xi, eta)      -- the rule under test here,
 *   * the through-thickness direction zeta     -- LobattoRule<NStations>,
 *                                                 inside the material.
 * Only the first is being changed: xi,eta = +/-1/sqrt(3)  ->  +/-1. These tests
 * exist to prove that nothing ELSE moved with it.
 */
#include "Marmot/GaussLobattoInterfaceFiniteElement.h"
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

  using Element = GaussLobattoInterfaceFiniteElement< 3, 8 >;
  using Scheme  = Element::SurfaceIntegrationScheme;

  constexpr int nDofs = 24;

  std::array< double, 24 > flatCoordinates = {
    -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0,
    -0.5, -0.5, 0.1, 0.5, -0.5, 0.1, 0.5, 0.5, 0.1, -0.5, 0.5, 0.1,
  };

  const std::vector< double > elasticProperties = { 1.0e5, 0.3, 0.01 };
  //                                      E,     nu,  h,    fy,  H,     Aexp, Hexp, id
  const std::vector< double > plasticProperties = { 1.0e5, 0.3, 0.01, 5.0, 1.0e3, 0.0, 0.0, 0.0 };

  std::unique_ptr< Element > makeElement( const std::string&           materialName,
                                          const std::vector< double >& materialProperties,
                                          Scheme                       scheme )
  {
    auto element = std::make_unique< Element >( 7,
                                                FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                Element::SectionType::Interface,
                                                scheme );
    element->assignNodeCoordinates( flatCoordinates.data() );

    static std::array< double, 1 > elementPropertyValues = { 1.0 };
    ElementProperties              elementProperties( elementPropertyValues.data(), 1 );
    element->assignProperty( elementProperties );
    element->assignMaterial( materialName, materialProperties.data(), static_cast< int >( materialProperties.size() ) );
    return element;
  }

  struct Evaluation {
    Eigen::Matrix< double, nDofs, 1 >     residual;
    Eigen::Matrix< double, nDofs, nDofs > tangent;
  };

  Evaluation evaluate( const std::string&                       materialName,
                       const std::vector< double >&             materialProperties,
                       Scheme                                   scheme,
                       const Eigen::Matrix< double, nDofs, 1 >& dU )
  {
    auto element = makeElement( materialName, materialProperties, scheme );

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
    // column-major: the convention of every Marmot element and of the host assembler
    evaluation.tangent = Eigen::Map< Eigen::Matrix< double, nDofs, nDofs > >( Ke.data() );
    return evaluation;
  }

  /** Through-thickness-uniform in-plane field: the integrand is CONSTANT over
   *  the surface, so any rule that is exact for constants must agree. */
  Eigen::Matrix< double, nDofs, 1 > makeConstantGradientIncrement()
  {
    constexpr int                     half = nDofs / 2;
    Eigen::Matrix< double, nDofs, 1 > dU   = Eigen::Matrix< double, nDofs, 1 >::Zero();
    for ( int node = 0; node < 4; ++node ) {
      const double x            = flatCoordinates[3 * node + 0];
      const double y            = flatCoordinates[3 * node + 1];
      const double ux           = 1.0e-4 * x;
      const double uy           = 0.5e-4 * y;
      dU( 3 * node + 0 )        = ux;
      dU( 3 * node + 1 )        = uy;
      dU( half + 3 * node + 0 ) = ux;
      dU( half + 3 * node + 1 ) = uy;
      dU( half + 3 * node + 2 ) = 2.0e-5;
    }
    return dU;
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
  // 1. The points really are where they are supposed to be, and only they
  //    moved: same count, same weights, same total measure.
  // ------------------------------------------------------------------
  void TestSurfacePointsAreAtTheNodesAndWeightsAreUnchanged()
  {
    const auto gauss   = Element::surfaceQuadrature( Marmot::FiniteElement::ElementShapes::Quad4,
                                                   FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                   Scheme::Gauss2x2 );
    const auto lobatto = Element::surfaceQuadrature( Marmot::FiniteElement::ElementShapes::Quad4,
                                                     FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                     Scheme::Lobatto2x2 );

    throwExceptionOnFailure( gauss.size() == 4 && lobatto.size() == 4,
                             "both surface rules must have four points, got " + std::to_string( gauss.size() ) +
                               " and " + std::to_string( lobatto.size() ) );

    double gaussWeight = 0.0, lobattoWeight = 0.0;
    for ( size_t i = 0; i < gauss.size(); ++i ) {
      gaussWeight += gauss[i].weight;
      lobattoWeight += lobatto[i].weight;

      for ( int d = 0; d < 2; ++d ) {
        throwExceptionOnFailure( std::abs( std::abs( gauss[i].xi( d ) ) - 1.0 / std::sqrt( 3.0 ) ) < 1e-14,
                                 "Gauss surface point " + std::to_string( i ) + " is not at +/-1/sqrt(3)." );
        throwExceptionOnFailure( std::abs( std::abs( lobatto[i].xi( d ) ) - 1.0 ) < 1e-14,
                                 "Lobatto surface point " + std::to_string( i ) + " is not at +/-1." );
      }
      throwExceptionOnFailure( std::abs( lobatto[i].weight - 1.0 ) < 1e-14,
                               "every 2x2 Lobatto surface weight must be 1." );
    }

    std::cout << "surface rules: 4 points each, total weight Gauss = " << gaussWeight << ", Lobatto = " << lobattoWeight
              << "\n";

    throwExceptionOnFailure( std::abs( gaussWeight - lobattoWeight ) < 1e-14,
                             "the two surface rules must carry the same total measure: " +
                               std::to_string( gaussWeight ) + " vs " + std::to_string( lobattoWeight ) );
  }

  // ------------------------------------------------------------------
  // 2. Nothing but the surface rule changed: the through-thickness discretisation
  //    (station count, and therefore the whole material state block) must be
  //    bit-identical between the two elements.
  // ------------------------------------------------------------------
  void TestThroughThicknessDiscretisationIsUntouched()
  {
    auto gauss   = makeElement( "VONMISES", plasticProperties, Scheme::Gauss2x2 );
    auto lobatto = makeElement( "VONMISES", plasticProperties, Scheme::Lobatto2x2 );

    const int nGauss   = gauss->getNumberOfRequiredStateVars();
    const int nLobatto = lobatto->getNumberOfRequiredStateVars();

    std::cout << "state vars: Gauss = " << nGauss << ", Lobatto = " << nLobatto
              << " (stations = " << Element::Material::nStations << ")\n";

    throwExceptionOnFailure( nGauss == nLobatto,
                             "the two surface rules must require identical state storage -- if they do not, the "
                             "through-thickness discretisation moved too: " +
                               std::to_string( nGauss ) + " vs " + std::to_string( nLobatto ) );
    throwExceptionOnFailure( gauss->getNumberOfQuadraturePoints() == lobatto->getNumberOfQuadraturePoints(),
                             "the two surface rules must have the same number of surface points." );
    throwExceptionOnFailure( Element::Material::nStations == 5,
                             "the production through-thickness rule must remain five Lobatto stations." );
  }

  // ------------------------------------------------------------------
  // 3. THE CONTROL, and the accuracy price of the new rule.
  //
  // Under a through-thickness-uniform, linear in-plane field the stress is
  // CONSTANT over the surface, so the residual integrand B^T sigma is LINEAR in
  // (xi, eta) -- grad(N) of a Q1 quad is linear, not constant. Both rules are
  // exact to degree 1 per direction, so the RESIDUAL must agree to round-off.
  //
  // The TANGENT must NOT be expected to agree. B^T C B is QUADRATIC in
  // (xi, eta), and there the two rules part company: 2x2 Gauss is exact to
  // degree 3 per direction, 2-point Lobatto only to degree 1. The Lobatto rule
  // therefore UNDER-INTEGRATES the stiffness -- measured here at ~125%
  // difference. That is not a bug, it is the defining property of the rule, and
  // it is the reason the zero-energy-mode check below is mandatory rather than
  // optional.
  // ------------------------------------------------------------------
  void TestConstantIntegrandGivesIdenticalResultsForBothRules()
  {
    const auto dU = makeConstantGradientIncrement();

    const auto gauss   = evaluate( "LINEARELASTIC", elasticProperties, Scheme::Gauss2x2, dU );
    const auto lobatto = evaluate( "LINEARELASTIC", elasticProperties, Scheme::Lobatto2x2, dU );

    const double residualError = ( gauss.residual - lobatto.residual ).lpNorm< Eigen::Infinity >() /
                                 std::max( 1.0e-30, gauss.residual.lpNorm< Eigen::Infinity >() );
    const double tangentError = ( gauss.tangent - lobatto.tangent ).lpNorm< Eigen::Infinity >() /
                                std::max( 1.0e-30, gauss.tangent.lpNorm< Eigen::Infinity >() );

    std::cout << "constant stress state: residual rel. diff = " << residualError
              << " (must vanish), tangent rel. diff = " << tangentError << " (under-integration, expected)\n";

    throwExceptionOnFailure( residualError < 1e-12,
                             "for a constant stress state both surface rules must give the same residual -- the "
                             "integrand is linear and both rules are exact to degree 1: relative difference " +
                               std::to_string( residualError ) );
    throwExceptionOnFailure( tangentError > 1e-3,
                             "the 2-point Lobatto rule is only exact to degree 1 per direction, so it MUST "
                             "under-integrate the quadratic stiffness integrand; a vanishing difference would mean "
                             "the rule was not actually applied. Got " +
                               std::to_string( tangentError ) );
  }

  // ------------------------------------------------------------------
  // 3b. THE RISK THE UNDER-INTEGRATION CREATES: spurious zero-energy modes.
  //
  // A quadrature that is exact only to degree 1 can fail to see deformation
  // modes that the exact stiffness resists -- classic hourglassing. Compare the
  // eigenvalue spectra of the two elastic element tangents and count the modes
  // with negligible energy. The interface element's legitimate null space is
  // the same for both rules (it comes from the kinematics, not the quadrature),
  // so ANY extra soft mode under Lobatto is spurious and would invalidate the
  // experiment before it starts.
  // ------------------------------------------------------------------
  void TestLobattoSurfaceIntroducesNoExtraZeroEnergyModes()
  {
    const auto dU = Eigen::Matrix< double, nDofs, 1 >::Zero().eval();

    const auto gauss   = evaluate( "LINEARELASTIC", elasticProperties, Scheme::Gauss2x2, dU );
    const auto lobatto = evaluate( "LINEARELASTIC", elasticProperties, Scheme::Lobatto2x2, dU );

    auto softModeCount = [&]( const Eigen::Matrix< double, nDofs, nDofs >& K ) {
      const Eigen::Matrix< double, nDofs, nDofs >      symmetric = 0.5 * ( K + K.transpose() );
      Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > solver( symmetric );
      const Eigen::VectorXd                            eigenvalues = solver.eigenvalues();
      const double                                     largest     = eigenvalues.cwiseAbs().maxCoeff();
      int                                              count       = 0;
      for ( int i = 0; i < eigenvalues.size(); ++i ) {
        if ( std::abs( eigenvalues( i ) ) < 1.0e-9 * largest ) {
          ++count;
        }
      }
      return std::make_pair( count, largest );
    };

    const auto [gaussSoft, gaussScale]     = softModeCount( gauss.tangent );
    const auto [lobattoSoft, lobattoScale] = softModeCount( lobatto.tangent );

    std::cout << "zero-energy modes (of " << nDofs << "): Gauss = " << gaussSoft << ", Lobatto = " << lobattoSoft
              << "  [max eigenvalue " << gaussScale << " vs " << lobattoScale << "]\n";

    throwExceptionOnFailure( lobattoSoft <= gaussSoft,
                             "HOURGLASSING: the 2-point Lobatto surface rule under-integrates the stiffness and has "
                             "introduced " +
                               std::to_string( lobattoSoft - gaussSoft ) +
                               " extra zero-energy mode(s) relative to 2x2 Gauss (" + std::to_string( lobattoSoft ) +
                               " vs " + std::to_string( gaussSoft ) +
                               "). The comparison would be measuring a rank-deficient element." );
  }

  // ------------------------------------------------------------------
  // 4. The Lobatto-surface element must still be a correct element: consistent
  //    tangent, elastic and plastic.
  // ------------------------------------------------------------------
  void checkTangent( const std::string&           materialName,
                     const std::vector< double >& properties,
                     const std::string&           label )
  {
    const auto dU = makeGeneralIncrement();

    const auto reference = evaluate( materialName, properties, Scheme::Lobatto2x2, dU );

    Eigen::Matrix< double, nDofs, nDofs > finiteDifference;
    const double                          perturbation = 1.0e-9;
    for ( int column = 0; column < nDofs; ++column ) {
      auto forward = dU, backward = dU;
      forward( column ) += perturbation;
      backward( column ) -= perturbation;
      const auto plus                = evaluate( materialName, properties, Scheme::Lobatto2x2, forward );
      const auto minus               = evaluate( materialName, properties, Scheme::Lobatto2x2, backward );
      finiteDifference.col( column ) = -( plus.residual - minus.residual ) / ( 2.0 * perturbation );
    }

    const double relativeError = ( reference.tangent - finiteDifference ).lpNorm< Eigen::Infinity >() /
                                 std::max( 1.0, finiteDifference.lpNorm< Eigen::Infinity >() );

    std::cout << "Lobatto-surface tangent vs FD (" << label << "): relative error = " << relativeError << "\n";

    throwExceptionOnFailure( relativeError < 1e-6,
                             "GLIQUAD4_SURFLOB tangent does not match central differences (" + label +
                               "): relative error " + std::to_string( relativeError ) );
  }

  void TestLobattoSurfaceTangentIsConsistent()
  {
    checkTangent( "LINEARELASTIC", elasticProperties, "elastic" );
    checkTangent( "VONMISES", plasticProperties, "plastic" );
  }

  // ------------------------------------------------------------------
  // 5. The rules must NOT agree on a general (non-constant) integrand --
  //    otherwise the experiment would be vacuous and any observed change in the
  //    oscillation could not be attributed to the quadrature.
  // ------------------------------------------------------------------
  void TestTheTwoRulesActuallyDifferOnANonConstantIntegrand()
  {
    const auto dU = makeGeneralIncrement();

    const auto gauss   = evaluate( "VONMISES", plasticProperties, Scheme::Gauss2x2, dU );
    const auto lobatto = evaluate( "VONMISES", plasticProperties, Scheme::Lobatto2x2, dU );

    const double difference = ( gauss.residual - lobatto.residual ).lpNorm< Eigen::Infinity >() /
                              std::max( 1.0e-30, gauss.residual.lpNorm< Eigen::Infinity >() );

    std::cout << "non-constant integrand: residual rel. difference between the rules = " << difference << "\n";

    throwExceptionOnFailure( difference > 1e-6,
                             "the two surface rules must differ on a non-constant integrand, otherwise the experiment "
                             "measures nothing; got " +
                               std::to_string( difference ) );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{ TestSurfacePointsAreAtTheNodesAndWeightsAreUnchanged,
                                                       TestThroughThicknessDiscretisationIsUntouched,
                                                       TestConstantIntegrandGivesIdenticalResultsForBothRules,
                                                       TestLobattoSurfaceIntroducesNoExtraZeroEnergyModes,
                                                       TestLobattoSurfaceTangentIsConsistent,
                                                       TestTheTwoRulesActuallyDifferOnANonConstantIntegrand };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
