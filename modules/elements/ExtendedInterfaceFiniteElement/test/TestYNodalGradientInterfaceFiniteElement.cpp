#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/YInterfaceFiniteElement.h"
#include "Marmot/YNodalGradientInterfaceFiniteElement.h"

#include <Eigen/Dense>
#include <algorithm>
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

  using E                       = YNodalGradientInterfaceFiniteElement;
  constexpr int    nElementDofs = E::sizeLoadVector; // 48
  constexpr double hMesh        = 0.1;               // top-face z = hMesh, = material h

  static std::array< double, 24 > flatCoordinates = {
    -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0,
    -0.5, -0.5, 0.1, 0.5, -0.5, 0.1, 0.5, 0.5, 0.1, -0.5, 0.5, 0.1,
  };
  static const std::array< std::array< double, 2 >, 4 > flatXY = {
    { { -0.5, -0.5 }, { 0.5, -0.5 }, { 0.5, 0.5 }, { -0.5, 0.5 } } };

  std::unique_ptr< E > makeElement( const std::string&           mat,
                                    const std::vector< double >& props,
                                    const double*                coords = flatCoordinates.data() )
  {
    auto el = std::make_unique< E >( 11,
                                     FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                     E::SectionType::Interface );
    el->assignNodeCoordinates( coords );
    static std::array< double, 1 > elp = { hMesh };
    ElementProperties              ep( elp.data(), 1 );
    el->assignProperty( ep );
    el->assignMaterial( mat, props.data(), static_cast< int >( props.size() ) );
    return el;
  }

  void initState( E& el, std::vector< double >& sv )
  {
    sv.assign( el.getNumberOfRequiredStateVars(), 0.0 );
    el.assignStateVars( sv.data(), static_cast< int >( sv.size() ) );
    el.initializeYourself();
    el.setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
  }

  struct Eval {
    Eigen::Matrix< double, nElementDofs, 1 >            Pe;
    Eigen::Matrix< double, nElementDofs, nElementDofs > Ke; // col-major, matching element output
  };

  Eval evaluate( const std::string&                              mat,
                 const std::vector< double >&                    props,
                 const Eigen::Matrix< double, nElementDofs, 1 >& QTotalBase,
                 const Eigen::Matrix< double, nElementDofs, 1 >& dQ,
                 const std::vector< double >*                    initialSv = nullptr )
  {
    auto                  el = makeElement( mat, props );
    std::vector< double > sv;
    if ( initialSv ) {
      sv = *initialSv;
      el->assignStateVars( sv.data(), static_cast< int >( sv.size() ) );
      el->initializeYourself();
    }
    else
      initState( *el, sv );

    const Eigen::Matrix< double, nElementDofs, 1 >    QTotal = QTotalBase + dQ;
    std::array< double, nElementDofs >                Qa{}, dQa{}, Pe{};
    std::array< double, nElementDofs * nElementDofs > Ke{};
    std::copy( QTotal.data(), QTotal.data() + nElementDofs, Qa.begin() );
    std::copy( dQ.data(), dQ.data() + nElementDofs, dQa.begin() );
    el->computeKernels( Qa.data(), dQa.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

    Eval e;
    e.Pe = Eigen::Map< Eigen::Matrix< double, nElementDofs, 1 > >( Pe.data() );
    e.Ke = Eigen::Map< Eigen::Matrix< double, nElementDofs, nElementDofs > >( Ke.data() ); // col-major
    return e;
  }

  // Compatible affine state: u_i = c_i + F_ij x_j on both faces (top at z=hMesh),
  // g_A = F_.3 (= [u]/h), t_A = given traction. tau1=x, tau2=y, n=z.
  Eigen::Matrix< double, nElementDofs, 1 > affineState( const Eigen::Vector3d& c,
                                                        const Eigen::Matrix3d& F,
                                                        const Eigen::Vector3d& tNodal )
  {
    Eigen::Matrix< double, nElementDofs, 1 > q;
    q.setZero();
    const Eigen::Vector3d gVal = F.col( 2 ); // dU/dz
    for ( int A = 0; A < 4; A++ ) {
      const Eigen::Vector3d xB( flatXY[A][0], flatXY[A][1], 0.0 );
      const Eigen::Vector3d xT( flatXY[A][0], flatXY[A][1], hMesh );
      q.segment< 3 >( E::offD + 3 * A )                = c + F * xB; // u- bottom
      q.segment< 3 >( E::offD + E::nSideDofU + 3 * A ) = c + F * xT; // u+ top
      q.segment< 3 >( E::offG + 3 * A )                = gVal;       // g
      q.segment< 3 >( E::offT + 3 * A )                = tNodal;     // t
    }
    return q;
  }

} // namespace

void Test_DofLayout()
{
  std::cout << "\n--- Test_DofLayout ---\n";
  throwExceptionOnFailure( E::sizeLoadVector == 48, "must be 48 DOF." );
  throwExceptionOnFailure( E::offD == 0 && E::offG == 24 && E::offT == 36, "offsets d/g/t = 0/24/36." );

  auto el = makeElement( "LINEARELASTIC", { 8.0e4, 0.3, hMesh } );
  throwExceptionOnFailure( el->getNDofPerElement() == 48, "getNDofPerElement 48." );
  const auto nf = el->getNodeFields();
  throwExceptionOnFailure( nf.size() == 8, "8 nodes." );
  for ( int i = 0; i < 4; i++ )
    throwExceptionOnFailure( nf[i].size() == 3 && nf[i][0] == "displacement" && nf[i][1] == "normalGradientAverage" &&
                               nf[i][2] == "commonTraction",
                             "bottom nodes carry [displacement, normalGradientAverage, commonTraction]." );
  for ( int i = 4; i < 8; i++ )
    throwExceptionOnFailure( nf[i].size() == 1 && nf[i][0] == "displacement", "top nodes carry [displacement]." );

  auto perm = el->getDofIndicesPermutationPattern();
  throwExceptionOnFailure( perm.size() == 48, "perm size 48." );
  std::vector< int > s = perm;
  std::sort( s.begin(), s.end() );
  for ( int i = 0; i < 48; i++ )
    throwExceptionOnFailure( s[i] == i, "perm is a bijection onto [0,48)." );
}

void Test_CompatibleAffinePatch()
{
  std::cout << "\n--- Test_CompatibleAffinePatch ---\n";
  const std::vector< double > props = { 8.0e4, 0.3, hMesh };

  Eigen::Matrix3d F;
  F << 0.5e-3, 0.1e-3, -0.2e-3, 0.15e-3, 0.4e-3, 0.05e-3, -0.1e-3, 0.2e-3, 0.6e-3;
  const Eigen::Vector3d c( 0.1e-3, -0.2e-3, 0.05e-3 );

  // first pass with zero traction to read t^mat (constant for an affine state)
  auto                  el = makeElement( "LINEARELASTIC", props );
  std::vector< double > sv;
  initState( *el, sv );
  Eigen::Matrix< double, nElementDofs, 1 >          q0 = affineState( c, F, Eigen::Vector3d::Zero() );
  std::array< double, nElementDofs >                Qa{}, dQa{}, Pe{};
  std::array< double, nElementDofs * nElementDofs > Ke{};
  std::copy( q0.data(), q0.data() + nElementDofs, Qa.begin() );
  std::copy( q0.data(), q0.data() + nElementDofs, dQa.begin() );
  el->computeKernels( Qa.data(), dQa.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

  const Eigen::Vector3d tMat0 = el->qps[0].managedStateVars->commonTraction;
  const Eigen::Vector3d gQp   = el->qps[0].managedStateVars->normalGradientAtQp;

  // g at QP must equal [u]/h = F.col(2)
  const double gDefect = ( gQp - F.col( 2 ) ).norm();
  std::cout << "affine: |g - [u]/h| = " << gDefect << ", t^mat = " << tMat0.transpose() << "\n";
  throwExceptionOnFailure( gDefect < 1e-12, "g must equal [u]/h for a compatible affine state." );

  // t+ = t- from the local kernel: check via generalizedForce equality of sigma+ n and sigma- n
  {
    auto sp = el->qps[0].material->getStateView( "stressPlus", el->qps[0].managedStateVars->materialStateVars.data() );
    auto sm = el->qps[0].material->getStateView( "stressMinus", el->qps[0].managedStateVars->materialStateVars.data() );
    Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > SP( sp.stateLocation );
    Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > SM( sm.stateLocation );
    const Eigen::Vector3d                                        n( 0, 0, 1 );
    const double teq = ( SP * n - SM * n ).norm() / std::max( 1.0, ( SP * n ).norm() );
    std::cout << "affine: |t+ - t-|/|t+| = " << teq << "\n";
    throwExceptionOnFailure( teq < 1e-8, "local traction equilibrium t+=t-." );
  }

  // now set nodal t = t^mat, re-evaluate; R_q and R_p must vanish
  const auto   full = evaluate( "LINEARELASTIC",
                              props,
                              Eigen::Matrix< double, nElementDofs, 1 >::Zero(),
                              affineState( c, F, tMat0 ) );
  const double Rq   = full.Pe.segment< E::nDofG >( E::offG ).norm();
  const double Rp   = full.Pe.segment< E::nDofT >( E::offT ).norm();
  std::cout << "affine: |R_q| = " << Rq << ", |R_p| = " << Rp << "\n";
  throwExceptionOnFailure( Rq < 1e-9, "R_q (t = t^mat) must vanish." );
  throwExceptionOnFailure( Rp < 1e-12, "R_p (g = [u]/h) must vanish." );
}

Eigen::Matrix< double, nElementDofs, nElementDofs > fdTangent( const std::string&                              mat,
                                                               const std::vector< double >&                    props,
                                                               const Eigen::Matrix< double, nElementDofs, 1 >& base,
                                                               const Eigen::Matrix< double, nElementDofs, 1 >& dQ,
                                                               const std::vector< double >*                    sv,
                                                               double                                          rel )
{
  Eigen::Matrix< double, nElementDofs, nElementDofs > Jfd;
  for ( int j = 0; j < nElementDofs; j++ ) {
    const double                             eps = rel * std::max( 1.0, std::abs( dQ( j ) ) );
    Eigen::Matrix< double, nElementDofs, 1 > dp = dQ, dm = dQ;
    dp( j ) += eps;
    dm( j ) -= eps;
    const auto ep = evaluate( mat, props, base, dp, sv );
    const auto em = evaluate( mat, props, base, dm, sv );
    Jfd.col( j )  = ( ep.Pe - em.Pe ) / ( 2.0 * eps );
  }
  return Jfd;
}

void reportTangent( const char*                                     label,
                    const std::string&                              mat,
                    const std::vector< double >&                    props,
                    const Eigen::Matrix< double, nElementDofs, 1 >& base,
                    const Eigen::Matrix< double, nElementDofs, 1 >& dQ,
                    const std::vector< double >*                    sv,
                    double                                          rel,
                    double                                          tol )
{
  const auto   base_e = evaluate( mat, props, base, dQ, sv );
  const auto   Jfd    = fdTangent( mat, props, base, dQ, sv, rel );
  const auto   diff   = Jfd - ( -base_e.Ke );
  const double scale  = std::max( 1.0, base_e.Ke.norm() );
  const double eD     = diff.middleCols( E::offD, E::nDofD ).norm() / scale;
  const double eG     = diff.middleCols( E::offG, E::nDofG ).norm() / scale;
  const double eT     = diff.middleCols( E::offT, E::nDofT ).norm() / scale;
  std::cout << label << ": disp-col err=" << eD << " g-col err=" << eG << " t-col err=" << eT << "\n";
  throwExceptionOnFailure( base_e.Ke.allFinite(), "tangent finite." );
  throwExceptionOnFailure( eD < tol && eG < tol && eT < tol, "tangent-FD error too large." );
}

void Test_TangentFD_Elastic()
{
  std::cout << "\n--- Test_TangentFD_Elastic ---\n";
  const std::vector< double > props = { 8.0e4, 0.3, hMesh };
  Eigen::Matrix3d             F;
  F << 0.3e-3, 0.1e-3, -0.1e-3, 0.1e-3, 0.2e-3, 0.05e-3, -0.05e-3, 0.15e-3, 0.4e-3;
  Eigen::Matrix< double, nElementDofs, 1 > dQ = affineState( Eigen::Vector3d( 0.05e-3, 0, 0 ),
                                                             F,
                                                             Eigen::Vector3d( 1.0, -0.5, 2.0 ) );
  for ( int j = 0; j < nElementDofs; j++ )
    dQ( j ) += 1.0e-6 * std::sin( 0.7 * j + 0.3 );
  reportTangent( "elastic",
                 "LINEARELASTIC",
                 props,
                 Eigen::Matrix< double, nElementDofs, 1 >::Zero(),
                 dQ,
                 nullptr,
                 1e-6,
                 1e-4 );
}

void Test_TangentFD_Plastic()
{
  std::cout << "\n--- Test_TangentFD_Plastic ---\n";
  const std::vector< double > props = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };
  Eigen::Matrix3d             Fh;
  Fh << -2.0e-3, 0.4e-3, 0.2e-3, -1.0e-3, -1.7e-3, 0.4e-3, 0.8e-3, -0.3e-3, 2.2e-3;
  const Eigen::Vector3d                    histT( 5.0, -2.0, 30.0 );
  Eigen::Matrix< double, nElementDofs, 1 > hist = affineState( Eigen::Vector3d::Zero(), Fh, histT );

  auto                  el = makeElement( "VONMISES", props );
  std::vector< double > sv;
  initState( *el, sv );
  {
    std::array< double, nElementDofs >                Qa{}, dQa{}, Pe{};
    std::array< double, nElementDofs * nElementDofs > Ke{};
    std::copy( hist.data(), hist.data() + nElementDofs, Qa.begin() );
    std::copy( hist.data(), hist.data() + nElementDofs, dQa.begin() );
    el->computeKernels( Qa.data(), dQa.data(), Pe.data(), Ke.data(), 0.0, 1.0 );
  }
  std::vector< double > histSv = sv;

  Eigen::Matrix< double, nElementDofs, 1 > dQ = 0.05 * hist;
  reportTangent( "plastic", "VONMISES", props, hist, dQ, &histSv, 1e-6, 5e-3 );
}

void Test_EquivalenceCommonTraction()
{
  std::cout << "\n--- Test_EquivalenceCommonTraction ---\n";
  // Under a uniform compatible state, the mixed element's common traction t^mat
  // must equal YIQUAD4's generalized force (both = <sigma> n).
  const std::vector< double > props = { 8.0e4, 0.3, hMesh };
  Eigen::Matrix3d             F;
  F << 0.5e-3, 0.1e-3, -0.2e-3, 0.15e-3, 0.4e-3, 0.05e-3, -0.1e-3, 0.2e-3, 0.6e-3;
  const Eigen::Vector3d c( 0, 0, 0 );

  auto                  el = makeElement( "LINEARELASTIC", props );
  std::vector< double > sv;
  initState( *el, sv );
  Eigen::Matrix< double, nElementDofs, 1 > q = affineState( c, F, Eigen::Vector3d::Zero() );
  {
    std::array< double, nElementDofs >                Qa{}, dQa{}, Pe{};
    std::array< double, nElementDofs * nElementDofs > Ke{};
    std::copy( q.data(), q.data() + nElementDofs, Qa.begin() );
    std::copy( q.data(), q.data() + nElementDofs, dQa.begin() );
    el->computeKernels( Qa.data(), dQa.data(), Pe.data(), Ke.data(), 0.0, 1.0 );
  }
  const Eigen::Vector3d tMat = el->qps[0].managedStateVars->commonTraction;

  // YIQUAD4 with the same displacement field
  Eigen::Matrix< double, 24, 1 > dUY;
  for ( int A = 0; A < 4; A++ ) {
    const Eigen::Vector3d xB( flatXY[A][0], flatXY[A][1], 0.0 );
    const Eigen::Vector3d xT( flatXY[A][0], flatXY[A][1], hMesh );
    dUY.segment< 3 >( 3 * A )      = c + F * xB;
    dUY.segment< 3 >( 12 + 3 * A ) = c + F * xT;
  }
  auto y = std::make_unique<
    YInterfaceFiniteElement< 3, 8 > >( 11,
                                       FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                       YInterfaceFiniteElement< 3, 8 >::SectionType::Interface );
  y->assignNodeCoordinates( flatCoordinates.data() );
  static std::array< double, 1 > elp = { hMesh };
  ElementProperties              ep( elp.data(), 1 );
  y->assignProperty( ep );
  y->assignMaterial( "LINEARELASTIC", props.data(), static_cast< int >( props.size() ) );
  std::vector< double > ysv( y->getNumberOfRequiredStateVars(), 0.0 );
  y->assignStateVars( ysv.data(), static_cast< int >( ysv.size() ) );
  y->initializeYourself();
  y->setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
  {
    std::array< double, 24 >      U{}, dq{}, Pe{};
    std::array< double, 24 * 24 > Ke{};
    std::copy( dUY.data(), dUY.data() + 24, dq.begin() );
    y->computeKernels( U.data(), dq.data(), Pe.data(), Ke.data(), 0.0, 1.0 );
  }
  const Eigen::Vector3d fY = y->qps[0].managedStateVars->generalizedForce;

  const double defect = ( tMat - fY ).norm() / std::max( 1.0, fY.norm() );
  std::cout << "t^mat = " << tMat.transpose() << ", YIQUAD4 f = " << fY.transpose() << ", rel defect = " << defect
            << "\n";
  throwExceptionOnFailure( defect < 1e-10, "common traction must equal YIQUAD4 generalized force." );
}

void Test_TwoElementSharedNodeAndRank()
{
  std::cout << "\n--- Test_TwoElementSharedNodeAndRank ---\n";
  const std::vector< double > props = { 8.0e4, 0.3, hMesh };

  std::array< double, 24 > cA = { -1.0, -0.5, 0.0, 0.0, -0.5, 0.0, 0.0, 0.5, 0.0, -1.0, 0.5, 0.0,
                                  -1.0, -0.5, 0.1, 0.0, -0.5, 0.1, 0.0, 0.5, 0.1, -1.0, 0.5, 0.1 };
  std::array< double, 24 > cB = { 0.0, -0.5, 0.0, 1.0, -0.5, 0.0, 1.0, 0.5, 0.0, 0.0, 0.5, 0.0,
                                  0.0, -0.5, 0.1, 1.0, -0.5, 0.1, 1.0, 0.5, 0.1, 0.0, 0.5, 0.1 };
  // global midsurface/displacement bottom nodes: 0=(-1,-.5) 1=(0,-.5) 2=(1,-.5) 3=(-1,.5) 4=(0,.5) 5=(1,.5)
  // top nodes 6..11 same xy. g,t live on bottom nodes (midsurface).
  const std::array< int, 4 > bA = { 0, 1, 4, 3 }, tA = { 6, 7, 10, 9 };
  const std::array< int, 4 > bB = { 1, 2, 5, 4 }, tB = { 7, 8, 11, 10 };
  constexpr int              nGN = 12,
                nGD              = 12 *
                      6; // 6 slots per node (3 disp + 3 g + 3 t on bottom; top only disp but keep 6 for simplicity)
  // canonical global dof: node k -> 6k + [disp(3), g(3)... ] but g/t only on bottom. Use per-node 9 for bottom, 3 for
  // top. Simpler: assign global slots: bottom node k(0..5): 9 slots (disp3,g3,t3); top node k(6..11): 3 slots (disp3).
  auto      gdofBottomDisp = [&]( int gn, int c ) { return 9 * gn + c; };
  auto      gdofBottomG    = [&]( int gn, int c ) { return 9 * gn + 3 + c; };
  auto      gdofBottomT    = [&]( int gn, int c ) { return 9 * gn + 6 + c; };
  auto      gdofTopDisp    = [&]( int gn, int c ) { return 9 * 6 + 3 * ( gn - 6 ) + c; };
  const int nGlobal        = 9 * 6 + 3 * 6; // 54 + 18 = 72

  auto localToGlobal = [&]( const std::array< int, 4 >& bmap, const std::array< int, 4 >& tmap ) {
    std::array< int, nElementDofs > m{};
    for ( int A = 0; A < 4; A++ )
      for ( int c = 0; c < 3; c++ ) {
        m[E::offD + 3 * A + c]                = gdofBottomDisp( bmap[A], c );
        m[E::offD + E::nSideDofU + 3 * A + c] = gdofTopDisp( tmap[A], c );
        m[E::offG + 3 * A + c]                = gdofBottomG( bmap[A], c );
        m[E::offT + 3 * A + c]                = gdofBottomT( bmap[A], c );
      }
    return m;
  };
  const auto mA = localToGlobal( bA, tA );
  const auto mB = localToGlobal( bB, tB );

  // shared bottom node (0,-0.5)=global 1: A local node 1, B local node 0
  throwExceptionOnFailure( mA[E::offT + 3 * 1] == mB[E::offT + 3 * 0],
                           "shared bottom node must map t to the same global dof." );
  throwExceptionOnFailure( mA[E::offG + 3 * 1] == mB[E::offG + 3 * 0],
                           "shared bottom node must map g to the same global dof." );

  auto                  elA = makeElement( "LINEARELASTIC", props, cA.data() );
  auto                  elB = makeElement( "LINEARELASTIC", props, cB.data() );
  std::vector< double > svA, svB;
  initState( *elA, svA );
  initState( *elB, svB );

  Eigen::Matrix< double, nElementDofs, 1 > qA = Eigen::Matrix< double, nElementDofs, 1 >::Zero();
  Eigen::Matrix< double, nElementDofs, 1 > qB = Eigen::Matrix< double, nElementDofs, 1 >::Zero();
  for ( int i = 0; i < nElementDofs; i++ ) {
    qA( i ) = 1e-4 * std::sin( 0.9 * i + 0.1 );
    qB( i ) = 1e-4 * std::sin( 1.3 * i + 0.7 );
  }

  auto run = [&]( E&                                                 el,
                  const Eigen::Matrix< double, nElementDofs, 1 >&    q,
                  std::array< double, nElementDofs >&                Pe,
                  std::array< double, nElementDofs * nElementDofs >& Ke ) {
    std::array< double, nElementDofs > Qa{}, dQa{};
    std::copy( q.data(), q.data() + nElementDofs, Qa.begin() );
    std::copy( q.data(), q.data() + nElementDofs, dQa.begin() );
    el.computeKernels( Qa.data(), dQa.data(), Pe.data(), Ke.data(), 0.0, 1.0 );
  };
  std::array< double, nElementDofs >                PeA{}, PeB{};
  std::array< double, nElementDofs * nElementDofs > KeA{}, KeB{};
  run( *elA, qA, PeA, KeA );
  run( *elB, qB, PeB, KeB );

  Eigen::MatrixXd K        = Eigen::MatrixXd::Zero( nGlobal, nGlobal );
  auto            assemble = [&]( const std::array< int, nElementDofs >&                   m,
                       const std::array< double, nElementDofs * nElementDofs >& Ke ) {
    for ( int i = 0; i < nElementDofs; i++ )
      for ( int j = 0; j < nElementDofs; j++ )
        K( m[i], m[j] ) += Ke[i + nElementDofs * j]; // element Ke is col-major
  };
  assemble( mA, KeA );
  assemble( mB, KeB );
  throwExceptionOnFailure( K.allFinite(), "assembled global K finite." );

  // Remove all 6 rigid-body modes via a 3-2-1 constraint on bottom displacement.
  // (disp dofs live at gdofBottomDisp(node,c) for bottom nodes.)
  std::vector< bool > fixed( nGlobal, false );
  fixed[gdofBottomDisp( 0, 0 )] = fixed[gdofBottomDisp( 0, 1 )] = fixed[gdofBottomDisp( 0, 2 )] = true;
  fixed[gdofBottomDisp( 2, 1 )] = fixed[gdofBottomDisp( 2, 2 )] = true;
  fixed[gdofBottomDisp( 3, 2 )]                                 = true;
  std::vector< int > freeDofs;
  for ( int i = 0; i < nGlobal; i++ )
    if ( !fixed[i] )
      freeDofs.push_back( i );

  Eigen::MatrixXd Kr( freeDofs.size(), freeDofs.size() );
  for ( size_t i = 0; i < freeDofs.size(); i++ )
    for ( size_t j = 0; j < freeDofs.size(); j++ )
      Kr( i, j ) = K( freeDofs[i], freeDofs[j] );

  Eigen::JacobiSVD< Eigen::MatrixXd > svd( Kr );
  const auto&                         s     = svd.singularValues();
  const double                        floor = 1e-10 * s( 0 );
  int                                 nz    = 0;
  for ( int i = 0; i < s.size(); i++ )
    if ( s( i ) < floor )
      nz++;
  std::cout << "two-element saddle system (" << Kr.rows() << "x" << Kr.cols() << "): smax=" << s( 0 )
            << " smin=" << s( s.size() - 1 ) << " near-zero=" << nz << "\n";
  throwExceptionOnFailure( s.allFinite(), "SVD finite." );
  throwExceptionOnFailure( nz == 0, "partial-mixed saddle system must be full rank after removing rigid-body modes." );
}

void Test_SingleElementEffectiveJumpStiffness()
{
  std::cout << "\n--- Test_SingleElementEffectiveJumpStiffness ---\n";
  // Solve the FULL coupled 48-dof element system with g and t as genuinely
  // FREE unknowns (only displacements prescribed), exactly as the global
  // solver would. Then verify the effective interface stiffness against the
  // analytic value: for a uniform normal jump delta with zero in-plane
  // strain, t_z = M * delta / h with the constrained modulus
  // M = E(1-nu)/((1+nu)(1-2nu)); for a tangential jump, t_x = G * delta / h.
  // This is the test that would have caught an interface that is orders of
  // magnitude too compliant in the assembled (mixed) system.
  const std::vector< double > props = { 8.0e4, 0.3, hMesh }; // E, nu, h ( = mesh separation )
  const double                Emod = 8.0e4, nu = 0.3;
  const double                M     = Emod * ( 1 - nu ) / ( ( 1 + nu ) * ( 1 - 2 * nu ) );
  const double                G     = Emod / ( 2 * ( 1 + nu ) );
  const double                delta = 1.0e-3;

  auto solveWithPrescribedJump = [&]( const Eigen::Vector3d& jump ) {
    // prescribed: bottom u = 0, top u = jump (all 4 top nodes)
    Eigen::Matrix< double, nElementDofs, 1 > q = Eigen::Matrix< double, nElementDofs, 1 >::Zero();
    for ( int A = 0; A < 4; A++ )
      q.segment< 3 >( E::offD + E::nSideDofU + 3 * A ) = jump;

    std::vector< int > freeDofs;
    for ( int i = E::offG; i < E::sizeLoadVector; i++ )
      freeDofs.push_back( i ); // g(12) + t(12) free
    const int nf = static_cast< int >( freeDofs.size() );

    const Eigen::Matrix< double, nElementDofs, 1 > zero = Eigen::Matrix< double, nElementDofs, 1 >::Zero();

    double resNorm = 1.0;
    for ( int it = 0; it < 12 && resNorm > 1e-11; it++ ) {
      const auto      ev = evaluate( "LINEARELASTIC", props, zero, q );
      Eigen::VectorXd Pf( nf );
      Eigen::MatrixXd Kf( nf, nf );
      for ( int i = 0; i < nf; i++ ) {
        Pf( i ) = ev.Pe( freeDofs[i] );
        for ( int j = 0; j < nf; j++ )
          Kf( i, j ) = ev.Ke( freeDofs[i], freeDofs[j] );
      }
      resNorm = Pf.norm();
      if ( resNorm <= 1e-11 )
        break;
      const Eigen::VectorXd dq = Kf.fullPivLu().solve( Pf ); // Pe = -R, Ke = dR/dq  =>  dq = K^-1 Pe
      for ( int i = 0; i < nf; i++ )
        q( freeDofs[i] ) += dq( i );
    }

    const auto evFinal = evaluate( "LINEARELASTIC", props, zero, q );
    struct Out {
      Eigen::Vector3d gNode0, tNode0;
      double          topReactionAlongJump;
      double          freeResidual;
    } out;
    out.gNode0       = q.segment< 3 >( E::offG );
    out.tNode0       = q.segment< 3 >( E::offT );
    out.freeResidual = 0.0;
    for ( int i = E::offG; i < E::sizeLoadVector; i++ )
      out.freeResidual = std::max( out.freeResidual, std::abs( evFinal.Pe( i ) ) );
    // reaction on the prescribed top dofs along the jump direction: R = -Pe
    const Eigen::Vector3d dir  = jump.normalized();
    double                reac = 0.0;
    for ( int A = 0; A < 4; A++ )
      reac += -evFinal.Pe.segment< 3 >( E::offD + E::nSideDofU + 3 * A ).dot( dir );
    out.topReactionAlongJump = reac;
    return out;
  };

  // ---- normal jump ----
  {
    const auto   r       = solveWithPrescribedJump( Eigen::Vector3d( 0, 0, delta ) );
    const double gExpect = delta / hMesh;     // 0.01
    const double tExpect = M * delta / hMesh; // ~1076.9
    const double gErr    = std::abs( r.gNode0( 2 ) - gExpect ) / gExpect;
    const double tErr    = std::abs( r.tNode0( 2 ) - tExpect ) / tExpect;
    const double reacErr = std::abs( r.topReactionAlongJump - tExpect ) / tExpect; // area = 1
    std::cout << "normal:    g_z=" << r.gNode0( 2 ) << " (expect " << gExpect << "),  t_z=" << r.tNode0( 2 )
              << " (expect " << tExpect << "),  top reaction=" << r.topReactionAlongJump
              << ",  residual=" << r.freeResidual << "\n";
    throwExceptionOnFailure( r.freeResidual < 1e-8, "free-dof residual must converge (normal jump)." );
    throwExceptionOnFailure( gErr < 1e-8, "g must equal [u]/h in the coupled solve (normal jump)." );
    throwExceptionOnFailure( tErr < 1e-6, "traction must equal M*delta/h -- effective normal stiffness wrong." );
    throwExceptionOnFailure( reacErr < 1e-6, "top reaction must equal t*Area -- equilibrium coupling wrong." );
  }

  // ---- tangential jump ----
  {
    const auto   r       = solveWithPrescribedJump( Eigen::Vector3d( delta, 0, 0 ) );
    const double gExpect = delta / hMesh;
    const double tExpect = G * delta / hMesh; // ~307.7
    const double gErr    = std::abs( r.gNode0( 0 ) - gExpect ) / gExpect;
    const double tErr    = std::abs( r.tNode0( 0 ) - tExpect ) / tExpect;
    const double reacErr = std::abs( r.topReactionAlongJump - tExpect ) / tExpect;
    std::cout << "tangential: g_x=" << r.gNode0( 0 ) << " (expect " << gExpect << "),  t_x=" << r.tNode0( 0 )
              << " (expect " << tExpect << "),  top reaction=" << r.topReactionAlongJump
              << ",  residual=" << r.freeResidual << "\n";
    throwExceptionOnFailure( r.freeResidual < 1e-8, "free-dof residual must converge (tangential jump)." );
    throwExceptionOnFailure( gErr < 1e-8, "g must equal [u]/h in the coupled solve (tangential jump)." );
    throwExceptionOnFailure( tErr < 1e-6, "traction must equal G*delta/h -- effective shear stiffness wrong." );
    throwExceptionOnFailure( reacErr < 1e-6, "top reaction must equal t*Area -- equilibrium coupling wrong." );
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    Test_DofLayout,
    Test_CompatibleAffinePatch,
    Test_TangentFD_Elastic,
    Test_TangentFD_Plastic,
    Test_EquivalenceCommonTraction,
    Test_TwoElementSharedNodeAndRank,
    Test_SingleElementEffectiveJumpStiffness,
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
