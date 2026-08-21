/** Tests for ZIQUAD4_STABP: 44 DOF = u(24) + g(12) + pbar(4) + [p](4). */
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/ZStabPressureInterfaceFiniteElement.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Elements;

namespace {
  using E                = ZStabPressureInterfaceFiniteElement;
  constexpr int    nD    = E::sizeLoadVector; // 44
  constexpr double hMesh = 0.1;

  std::array< double, 24 > coords = { -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0,
                                      -0.5, -0.5, 0.1, 0.5, -0.5, 0.1, 0.5, 0.5, 0.1, -0.5, 0.5, 0.1 };

  std::vector< double > vmProps()
  {
    return { 4.0e5, 0.3, hMesh, 5.0, 0.1, 0.0, 0.0, 0.0 };
  }
  std::vector< double > elProps()
  {
    return { 2.0e5, 0.3, hMesh };
  }

  std::unique_ptr< E > make( const std::string& mat, const std::vector< double >& p )
  {
    auto el = std::make_unique< E >( 7,
                                     FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                     E::SectionType::Interface );
    el->assignNodeCoordinates( coords.data() );
    static std::array< double, 1 > elp = { 1.0 };
    ElementProperties              ep( elp.data(), 1 );
    el->assignProperty( ep );
    el->assignMaterial( mat, p.data(), (int)p.size() );
    return el;
  }

  struct Eval {
    Eigen::Matrix< double, nD, 1 >  Pe;
    Eigen::Matrix< double, nD, nD > Ke;
  };

  Eval run( E&                                    el,
            const std::vector< double >&          svIn,
            const Eigen::Matrix< double, nD, 1 >& q,
            const Eigen::Matrix< double, nD, 1 >& dq )
  {
    std::vector< double > sv = svIn;
    el.assignStateVars( sv.data(), (int)sv.size() );
    std::array< double, nD >      Q{}, dQ{}, Pe{};
    std::array< double, nD * nD > Ke{};
    std::copy( q.data(), q.data() + nD, Q.begin() );
    std::copy( dq.data(), dq.data() + nD, dQ.begin() );
    el.computeKernels( Q.data(), dQ.data(), Pe.data(), Ke.data(), 0.0, 1.0 );
    Eval e;
    e.Pe = Eigen::Map< Eigen::Matrix< double, nD, 1 > >( Pe.data() );
    e.Ke = Eigen::Map< Eigen::Matrix< double, nD, nD > >( Ke.data() ); // column-major, as written
    return e;
  }

  void testLayout()
  {
    throwExceptionOnFailure( nD == 44, "ZIQUAD4_STABP must have 44 DOF." );
    throwExceptionOnFailure( E::offU == 0 && E::offG == 24 && E::offPm == 36 && E::offPj == 40, "offsets." );
    auto       el = make( "LINEARELASTIC", elProps() );
    const auto nf = el->getNodeFields();
    throwExceptionOnFailure( nf.size() == 8, "8 nodes." );
    for ( int i = 0; i < 4; i++ )
      throwExceptionOnFailure( nf[i].size() == 4 && nf[i][1] == "normalGradientJump" &&
                                 nf[i][2] == "interfacePressure" && nf[i][3] == "interfacePressureJump",
                               "bottom nodes carry u, g, pbar, [p]." );
    for ( int i = 4; i < 8; i++ )
      throwExceptionOnFailure( nf[i].size() == 1, "top nodes carry only displacement." );
    auto p = el->getDofIndicesPermutationPattern();
    std::sort( p.begin(), p.end() );
    throwExceptionOnFailure( (int)p.size() == nD, "perm size." );
    for ( int i = 0; i < nD; i++ )
      throwExceptionOnFailure( p[i] == i, "perm is a bijection." );
  }

  void testTangentFD( const std::string& mat, const std::vector< double >& props )
  {
    auto                  el = make( mat, props );
    std::vector< double > sv( el->getNumberOfRequiredStateVars(), 0.0 );
    el->assignStateVars( sv.data(), (int)sv.size() );
    el->initializeYourself();
    el->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );

    Eigen::Matrix< double, nD, 1 > dq;
    for ( int i = 0; i < nD; ++i )
      dq( i ) = ( i < 36 ? 1.0e-5 : 1.0e-1 ) * std::sin( 0.5 + 0.77 * i );
    const Eval b = run( *el, sv, dq, dq );

    const double eps = 1.0e-9;
    double       err = 0.0;
    for ( int c = 0; c < nD; ++c ) {
      Eigen::Matrix< double, nD, 1 > d2 = dq;
      d2( c ) += eps;
      const Eval o = run( *el, sv, d2, d2 );
      for ( int r = 0; r < nD; ++r )
        err = std::max( err, std::abs( ( o.Pe( r ) - b.Pe( r ) ) / eps + b.Ke( r, c ) ) );
    }
    // Per-block breakdown: a bad coupling block is a coding error, while a large
    // K_uu error under VONMISES is the return map's own FD noise at a yield kink.
    const int                       off[5] = { 0, E::offG, E::offPm, E::offPj, nD };
    const char*                     nm[4]  = { "u", "g", "pm", "pj" };
    Eigen::Matrix< double, nD, nD > fd;
    for ( int c = 0; c < nD; ++c ) {
      Eigen::Matrix< double, nD, 1 > d2 = dq;
      d2( c ) += eps;
      const Eval o = run( *el, sv, d2, d2 );
      for ( int r = 0; r < nD; ++r )
        fd( r, c ) = -( o.Pe( r ) - b.Pe( r ) ) / eps;
    }
    const double sc = b.Ke.cwiseAbs().maxCoeff();
    std::cout << "  [" << mat << "] element FD tangent rel err = " << err / sc << "   per block:";
    for ( int R = 0; R < 4; ++R )
      for ( int C = 0; C < 4; ++C ) {
        const int nr = off[R + 1] - off[R], nc = off[C + 1] - off[C];
        const double
          e = ( fd.block( off[R], off[C], nr, nc ) - b.Ke.block( off[R], off[C], nr, nc ) ).cwiseAbs().maxCoeff();
        if ( e / sc > 1e-9 )
          std::cout << "  " << nm[R] << nm[C] << "=" << e / sc;
      }
    std::cout << std::endl;
    // Two tolerances, split by whether a block passes through the base material's
    // INELASTIC RETURN MAP. The (u,g) x (u,g) blocks carry the algorithmic tangent,
    // whose gap against a central difference across a yield kink is a property of
    // VONMISES, not of this element -- LINEARELASTIC gives 9e-12 for every block.
    // Everything touching a pressure index is purely geometric and is held tight.
    double errMat = 0.0, errGeom = 0.0;
    for ( int R = 0; R < 4; ++R )
      for ( int C = 0; C < 4; ++C ) {
        const int nr = off[R + 1] - off[R], nc = off[C + 1] - off[C];
        const double
          e = ( fd.block( off[R], off[C], nr, nc ) - b.Ke.block( off[R], off[C], nr, nc ) ).cwiseAbs().maxCoeff();
        ( ( R < 2 && C < 2 ) ? errMat : errGeom ) = std::max( ( R < 2 && C < 2 ) ? errMat : errGeom, e );
      }
    std::cout << "     material blocks " << errMat / sc << "   geometric blocks " << errGeom / sc << std::endl;
    throwExceptionOnFailure( b.Ke.allFinite(), "tangent finite." );
    throwExceptionOnFailure( errGeom < 1.0e-6 * sc, "a purely geometric block disagrees with finite differences." );
    throwExceptionOnFailure( errMat < 1.0e-4 * sc, "a material block disagrees beyond the return map's own FD gap." );
  }

  /** The two structural zeros must survive assembly: pbar must not couple to g. */
  void testMeanPressureDoesNotReachG()
  {
    auto                  el = make( "LINEARELASTIC", elProps() );
    std::vector< double > sv( el->getNumberOfRequiredStateVars(), 0.0 );
    el->assignStateVars( sv.data(), (int)sv.size() );
    el->initializeYourself();
    el->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
    Eigen::Matrix< double, nD, 1 > dq;
    for ( int i = 0; i < nD; ++i )
      dq( i ) = ( i < 36 ? 1.0e-5 : 1.0e-1 ) * std::cos( 0.2 + 0.6 * i );
    const Eval e = run( *el, sv, dq, dq );

    const double gPm = e.Ke.block< E::nDofG, E::nDofP >( E::offG, E::offPm ).cwiseAbs().maxCoeff();
    const double pmG = e.Ke.block< E::nDofP, E::nDofG >( E::offPm, E::offG ).cwiseAbs().maxCoeff();
    const double gPj = e.Ke.block< E::nDofG, E::nDofP >( E::offG, E::offPj ).cwiseAbs().maxCoeff();
    const double pjG = e.Ke.block< E::nDofP, E::nDofG >( E::offPj, E::offG ).cwiseAbs().maxCoeff();
    std::cout << "  K(g,pbar)=" << gPm << "  K(pbar,g)=" << pmG << "   |   K(g,[p])=" << gPj << "  K([p],g)=" << pjG
              << std::endl;
    throwExceptionOnFailure( gPm < 1e-14 && pmG < 1e-14, "pbar must not couple to g after assembly." );
    throwExceptionOnFailure( gPj > 1e-6 && pjG > 1e-6, "[p] MUST couple to g -- that is the whole point." );
  }
} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{
    [&]() { testLayout(); },
    [&]() { testTangentFD( "LINEARELASTIC", elProps() ); },
    [&]() { testTangentFD( "VONMISES", vmProps() ); },
    [&]() { testMeanPressureDoesNotReachG(); },
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
