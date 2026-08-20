/* Diagnostic: WHICH pressure mode does the MINI bubble actually rescue?
 *
 * Hypothesis tested: the Q1 displacement space is weakly coupled to the
 * alternating (checkerboard) pressure mode, and the internal MINI displacement
 * bubble supplies the missing coupling.
 *
 * MEASURED RESULT -- the second half of that hypothesis is FALSE.
 *
 *   pressure mode      |B_u^T q|     |B_a^T q|    B_u rel. to constant
 *   constant           4.322e-04     6.8e-21              1.000
 *   linear xi          1.386e-04     1.481e-04            0.321
 *   linear eta         1.386e-04     1.481e-04            0.321
 *   alternating        4.431e-05     6.8e-21              0.103
 *
 * The Q1 part holds: the checkerboard is the weakest of the four modes, ~10x
 * below the constant mode. But the bubble does NOT rescue it -- its coupling to
 * the alternating mode is 6.8e-21, i.e. EXACTLY zero to machine precision, and
 * likewise exactly zero to the constant mode. The bubble couples to the two
 * LINEAR pressure modes and to nothing else.
 *
 * Both zeros are structural, not numerical:
 *   * constant mode: the bubble vanishes on the element boundary, so
 *     int grad(b) dA = 0 and a constant pressure exerts no force on it;
 *   * alternating mode: on a regular element sum_qp (N.c) grad(b) J0xW cancels
 *     by parity, c being the checkerboard nodal pattern.
 *
 * The coupling is also EXACTLY geometric -- zero drift from an elastic state
 * (kappa = 0) to a fully plastic one (kappa = 0.117) -- which confirms in
 * measurement what the material only documents in a comment (d r_p / d gamma = 0).
 *
 * RESOLUTION -- the bubble IS the mechanism; the mode was misidentified.
 *
 * The ablation (MARMOT_MINI_BUBBLE=OFF, everything else byte-identical) makes
 * the oscillation return: interior-row alternation 0.0034 -> 0.1346, a factor
 * of 39. So the bubble is decisive despite being blind to q_alt.
 *
 * A spatial spectrum of the interior row explains why. The spurious mode that
 * the bubble suppresses is NOT the two-element checkerboard: with the bubble
 * off, the fluctuation is dominated by 3-5 ELEMENT wavelengths (amplitudes
 * 0.038 at 4.0 elements, 0.028 at 3.4, 0.025 at 4.8), while the true
 * two-element Nyquist content is only 0.5% of the total. Restricted to ONE
 * element, a 3-5 element wavelength is essentially the LINEAR pressure mode --
 * exactly and only the mode the bubble couples to.
 *
 * So the local table and the ablation agree once the naming is fixed:
 *   * element-local q_alt (per-element checkerboard): bubble blind to it, and
 *     it carries almost none of the actual oscillation;
 *   * element-local LINEAR modes: bubble couples strongly, and they are what
 *     the global short-wavelength oscillation is built from.
 *
 * The element's own local coupling matrices are used, NOT a generic solid
 * element expression:
 *
 *     B_u = d(R_p)/d(u)      4 x 24   pressure  <-> Q1 displacement
 *     B_a = d(R_p)/d(beta)   4 x  3   pressure  <-> MINI bubble
 *
 * with R_p the element's actual volumetric residual, which carries the full
 * interface kinematics -- the jump term, both surface-gradient terms and the
 * through-thickness reconstruction -- integrated on the real nested
 * surface(2x2 Gauss)/thickness quadrature with the real interface measure.
 * Neither quadrature is changed here.
 *
 * Geometry and material are those of a REAL interior element of the production
 * benchmark (stiff, angle 10 deg, h = 0.01), taken from the converged mesh at
 * the centre of the interface, where the checkerboard is strong.
 */
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/YStabPressureMiniInterfaceFiniteElement.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Elements;

namespace {

  using Element = YStabPressureMiniInterfaceFiniteElement;

  constexpr int nDofs = Element::sizeLoadVector; // 28
  constexpr int nDofU = Element::nDofU;          // 24
  constexpr int nDofP = Element::nDofP;          // 4

  /** Interior element 434 of the production interface mesh, centroid
   *  (0.016665, 0.016665, 0.00206). Bottom face first, counter-clockwise. */
  std::array< double, 24 > interiorElementCoordinates = {
    -0.000000000, 0.000000000,  0.000000000, 0.033330001, 0.000000000,  -0.005880000, 0.033330001, 0.033330001,
    -0.005880000, -0.000000000, 0.033330001, 0.000000000, -0.000000000, 0.000000000,  0.010000000, 0.033330001,
    0.000000000,  0.004120000,  0.033330001, 0.033330001, 0.004120000,  -0.000000000, 0.033330001, 0.010000000,
  };

  /** The production interface material: E_0, nu_0, h, fy, H, Aexp, Hexp, id. */
  const std::vector< double > interfaceProperties = { 4.0e5, 0.3, 0.00984807753012208, 5.0, 0.1, 0.0, 0.0, 0.0 };

  std::unique_ptr< Element > makeElement()
  {
    auto element = std::make_unique< Element >( 434,
                                                FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                                Element::Interface );
    element->assignNodeCoordinates( interiorElementCoordinates.data() );

    static std::array< double, 1 > elementPropertyValues = { 1.0 };
    ElementProperties              elementProperties( elementPropertyValues.data(), 1 );
    element->assignProperty( elementProperties );
    element->assignMaterial( "VONMISES", interfaceProperties.data(), static_cast< int >( interfaceProperties.size() ) );
    return element;
  }

  /** The four Q1 pressure modes on the midsurface, in the element's node order
   *  (bottom face, counter-clockwise). Each has unit norm. */
  struct PressureMode {
    std::string                       name;
    Eigen::Matrix< double, nDofP, 1 > q;
  };

  std::vector< PressureMode > pressureModes()
  {
    // NOTE the explicit return type: `return 0.5 * q;` with a deduced type
    // returns an Eigen EXPRESSION referencing the local q, which dangles.
    auto make = []( double a, double b, double c, double d ) -> Eigen::Matrix< double, nDofP, 1 > {
      Eigen::Matrix< double, nDofP, 1 > q;
      q << a, b, c, d;
      return 0.5 * q;
    };
    return { { "constant      q0", make( 1, 1, 1, 1 ) },
             { "linear xi     qx", make( -1, 1, 1, -1 ) },
             { "linear eta    qe", make( -1, -1, 1, 1 ) },
             { "alternating qxe", make( 1, -1, 1, -1 ) } };
  }

  /** Drive the element through nSteps equal increments of a shear-dominated
   *  field and report the coupling matrices at the converged state. */
  struct State {
    Eigen::Matrix< double, nDofP, nDofU > Bu;
    Eigen::Matrix< double, nDofP, 3 >     Ba;
    double                                kappa; // equivalent plastic strain, top station
  };

  State driveTo( int nSteps, double amplitude )
  {
    auto element = makeElement();

    std::vector< double > stateVars( element->getNumberOfRequiredStateVars(), 0.0 );
    element->assignStateVars( stateVars.data(), static_cast< int >( stateVars.size() ) );
    element->initializeYourself();
    element->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

    // in-plane shear of the top face relative to the bottom, the deformation
    // the production benchmark actually applies at 10 degrees
    Eigen::Matrix< double, nDofs, 1 > increment = Eigen::Matrix< double, nDofs, 1 >::Zero();
    for ( int node = 0; node < 4; ++node ) {
      increment( 12 + 3 * node + 0 ) = amplitude * std::cos( 10.0 * M_PI / 180.0 );
      increment( 12 + 3 * node + 2 ) = -amplitude * std::sin( 10.0 * M_PI / 180.0 );
    }

    std::array< double, nDofs >         total{};
    std::array< double, nDofs >         Pe{};
    std::array< double, nDofs * nDofs > Ke{};

    for ( int step = 0; step < nSteps; ++step ) {
      for ( int i = 0; i < nDofs; ++i ) {
        total[i] += increment( i );
      }
      element->computeKernels( total.data(), increment.data(), Pe.data(), Ke.data(), step * 1.0, 1.0 );
    }

    State state;
    state.Bu = element->pressureDisplacementCoupling;
    state.Ba = element->pressureBubbleCoupling;
    // equivalent plastic strain at the top station: the material forwards the
    // "Top"/"Bottom" suffix to its embedded bulk materials
    state.kappa = element->getStateView( "kappaTop", 0 ).stateLocation[0];
    return state;
  }

  void reportTable( const std::string& label, const State& state )
  {
    std::printf( "\n--- %s   (equivalent plastic strain kappa = %.6e)\n", label.c_str(), state.kappa );
    std::printf( "%-18s %16s %16s %14s\n", "pressure mode", "|B_u^T q|", "|B_a^T q|", "ratio a/u" );

    double constantCoupling = 0.0;
    for ( const auto& mode : pressureModes() ) {
      const double cQ1     = ( state.Bu.transpose() * mode.q ).norm();
      const double cBubble = ( state.Ba.transpose() * mode.q ).norm();
      if ( constantCoupling == 0.0 ) {
        constantCoupling = cQ1;
      }
      std::printf( "%-18s %16.6e %16.6e %14.4f\n",
                   mode.name.c_str(),
                   cQ1,
                   cBubble,
                   cQ1 > 0.0 ? cBubble / cQ1 : INFINITY );
    }
    std::printf( "%-18s %16s %16s\n", "  (relative to the constant mode)", "", "" );
    for ( const auto& mode : pressureModes() ) {
      const double cQ1 = ( state.Bu.transpose() * mode.q ).norm();
      std::printf( "%-18s %16.6f\n", mode.name.c_str(), cQ1 / constantCoupling );
    }
  }

  void TestPressureModeCoupling()
  {
    // (1) early, still elastic  (2) later, plasticity developed.
    // The later state uses the ACTUAL converged element state and algorithmic
    // tangent -- the material is never reset to elastic.
    const State elastic = driveTo( 1, 1.0e-8 );
    const State plastic = driveTo( 40, 5.0e-5 );

    std::cout << "\n================================================================================\n";
    std::cout << "MINI: coupling of each Q1 pressure mode to the displacement spaces\n";
    std::cout << "real interior element, real nested quadrature, real interface kinematics\n";
    std::cout << "================================================================================";
    reportTable( "EARLY (elastic)", elastic );
    reportTable( "LATER (plastic)", plastic );

    // Report, do not presume: print the decisive ratios explicitly.
    auto modes    = pressureModes();
    auto coupling = []( const State& s, const Eigen::Matrix< double, nDofP, 1 >& q ) {
      return std::make_pair( ( s.Bu.transpose() * q ).norm(), ( s.Ba.transpose() * q ).norm() );
    };
    const auto [uConst, aConst] = coupling( plastic, modes[0].q );
    const auto [uAlt, aAlt]     = coupling( plastic, modes[3].q );

    std::printf( "\nDECISIVE RATIOS (plastic state):\n" );
    std::printf( "  |B_u^T q_alt| / |B_u^T q_const| = %.6e   (small => Q1 barely sees the checkerboard)\n",
                 uAlt / uConst );
    std::printf( "  |B_a^T q_alt| / |B_u^T q_alt|   = %.6e   (MEASURED ~1e-16: the bubble is blind to it)\n",
                 uAlt > 0.0 ? aAlt / uAlt : INFINITY );
    std::printf( "  |B_a^T q_lin| / |B_u^T q_lin|   = %.6e   (the bubble couples to the LINEAR modes instead)\n",
                 ( plastic.Bu.transpose() * modes[1].q ).norm() > 0.0
                   ? ( plastic.Ba.transpose() * modes[1].q ).norm() / ( plastic.Bu.transpose() * modes[1].q ).norm()
                   : INFINITY );

    const auto [uLin, aLin] = coupling( plastic, modes[1].q );

    throwExceptionOnFailure( std::isfinite( uConst ) && std::isfinite( uAlt ) && std::isfinite( aAlt ),
                             "coupling norms must be finite." );
    throwExceptionOnFailure( uConst > 0.0, "the constant pressure mode must couple to the Q1 displacement space." );

    // The Q1 half of the hypothesis: the checkerboard IS the weakest mode the
    // Q1 displacement space sees -- but it is not invisible.
    throwExceptionOnFailure( uAlt < 0.2 * uConst && uAlt > 0.0,
                             "the alternating pressure mode should be the WEAKEST but non-zero Q1 coupling; measured "
                             "ratio " +
                               std::to_string( uAlt / uConst ) );

    // The bubble half: REFUTED, and pinned here so a future change that alters
    // it is visible. The bubble sees the linear modes and nothing else.
    throwExceptionOnFailure( aLin > 0.5 * uLin,
                             "the MINI bubble must couple appreciably to the LINEAR pressure modes; measured " +
                               std::to_string( aLin ) + " vs Q1 " + std::to_string( uLin ) );
    throwExceptionOnFailure( aAlt < 1.0e-12 * aLin,
                             "MECHANISM CHANGED: the MINI bubble is structurally blind to the alternating pressure "
                             "mode (parity cancellation of sum_qp (N.c) grad(b)), but measured " +
                               std::to_string( aAlt ) + " against " + std::to_string( aLin ) + " for a linear mode." );
    throwExceptionOnFailure( aConst < 1.0e-12 * aLin,
                             "the MINI bubble must be blind to a CONSTANT pressure (int grad(b) dA = 0), measured " +
                               std::to_string( aConst ) );
  }

  /** The coupling blocks of THIS material are documented as purely geometric
   *  (d r_p / d gamma = 0), so they must not drift as plasticity develops. This
   *  turns that documented claim into a measured one. */
  void TestCouplingIsGeometricAndDoesNotDriftWithPlasticity()
  {
    const State elastic = driveTo( 1, 1.0e-8 );
    const State plastic = driveTo( 40, 5.0e-5 );

    const double drift = ( plastic.Bu - elastic.Bu ).lpNorm< Eigen::Infinity >() /
                         std::max( 1e-30, elastic.Bu.lpNorm< Eigen::Infinity >() );
    const double driftBubble = ( plastic.Ba - elastic.Ba ).lpNorm< Eigen::Infinity >() /
                               std::max( 1e-30, elastic.Ba.lpNorm< Eigen::Infinity >() );

    std::printf( "\ncoupling drift elastic -> plastic:  B_u %.3e,  B_a %.3e\n", drift, driftBubble );
    std::printf( "plasticity check: kappa elastic = %.3e, kappa plastic = %.3e\n", elastic.kappa, plastic.kappa );

    throwExceptionOnFailure( elastic.kappa == 0.0 && plastic.kappa > 1.0e-4,
                             "the two states must genuinely bracket yield (kappa " + std::to_string( elastic.kappa ) +
                               " -> " + std::to_string( plastic.kappa ) +
                               "), otherwise the comparison proves "
                               "nothing." );
    throwExceptionOnFailure( drift < 1e-10 && driftBubble < 1e-10,
                             "the pressure coupling blocks are documented as purely geometric, so they must not "
                             "change as plasticity develops; measured drift " +
                               std::to_string( drift ) + " / " + std::to_string( driftBubble ) );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{ TestPressureModeCoupling,
                                                       TestCouplingIsGeometricAndDoesNotDriftWithPlasticity };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
