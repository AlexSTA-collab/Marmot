// Tests of DRUCKERPRAGER (libMarmotDRUCKERPRAGER.so), run at one material point:
//   1  phi = psi = 0 reproduces VONMISES (sigma_y = sqrt3 c0, H_J2 = 3 H): stress and tangent along a mixed path
//   2  consistency: every plastic step ends on the yield surface (cone) or at the apex
//   3  the consistent tangent against central finite differences, cone and apex, associated and non-associated
//   4  ellipticity along the interface path (shear with compression, alpha = 10 deg): the smallest det of the
//      acoustic tensor over all directions, and at the interface normal, relative to the elastic one
//   5  (with pt) the tension cap: cap-only and corner returns (p = pt, deviator elastic / on the cone), and a
//      random tensile path: consistency with both surfaces and the tangent against central finite differences
// usage: test_dp [E nu c0 phi psi H [pt]]   (defaults: a soil-like layer, E = 100 MPa; phi = 0 < psi: non-associated
// von Mises)
#include "Marmot/DruckerPrager.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <random>
#include <vector>

using namespace Marmot;
using Mat = MarmotMaterialHypoElastic;
using MarmotLibrary::MarmotMaterialHypoElasticFactory;

struct Point {
  std::unique_ptr< Mat > m;
  Mat::state3D           s;
  std::vector< double >  sv = std::vector< double >( 16, 0. );
  Point( const char* name, std::vector< double >& p )
    : m( MarmotMaterialHypoElasticFactory::createMaterial( name, p.data(), (int)p.size(), 1 ) )
  {
    s.stateVars = sv.data();
  }
  Matrix6d step( const Vector6d& dE )
  {
    Matrix6d C;
    m->computeStress( s, C, dE, { 1., 1. } );
    return C;
  }
  // stress after a trial increment from the current state, the state left untouched
  Vector6d trial( const Vector6d& dE, Matrix6d* C = nullptr ) const
  {
    Mat::state3D t     = s;
    auto         svCpy = sv;
    t.stateVars        = svCpy.data();
    Matrix6d Ct;
    m->computeStress( t, Ct, dE, { 1., 1. } );
    if ( C )
      *C = Ct;
    return t.stress;
  }
};

static int failures = 0;
static void check( bool ok, const char* what, double v )
{
  std::printf( "  %-62s %.3e  %s\n", what, v, ok ? "ok" : "FAIL" );
  failures += !ok;
}

static double fdTangentError( const Point& P, const Vector6d& dE )
{
  Matrix6d C;
  P.trial( dE, &C );
  Matrix6d     Cfd;
  const double h = 1e-7 * std::max( dE.norm(), 1e-6 );
  for ( int j = 0; j < 6; ++j ) {
    Vector6d e = Vector6d::Zero();
    e( j )     = h;
    Cfd.col( j ) = ( P.trial( dE + e ) - P.trial( dE - e ) ) / ( 2. * h );
  }
  return ( C - Cfd ).norm() / Cfd.norm();
}

// acoustic tensor Q_ik = D_ijkl n_j n_l of a Voigt tangent (order 11 22 33 12 13 23, engineering shear strain)
static Eigen::Matrix3d acoustic( const Matrix6d& D, const Eigen::Vector3d& n )
{
  const int       V[3][3] = { { 0, 3, 4 }, { 3, 1, 5 }, { 4, 5, 2 } };
  Eigen::Matrix3d Q       = Eigen::Matrix3d::Zero();
  for ( int i = 0; i < 3; ++i )
    for ( int k = 0; k < 3; ++k )
      for ( int j = 0; j < 3; ++j )
        for ( int l = 0; l < 3; ++l )
          Q( i, k ) += D( V[i][j], V[k][l] ) * n( j ) * n( l );
  return Q;
}

static double minDetRatio( const Matrix6d& D, const Matrix6d& Del, Eigen::Vector3d* worst = nullptr )
{
  double m = 1e300;
  for ( int a = 0; a <= 90; ++a )
    for ( int b = 0; b < 180; ++b ) {
      const double    th = a * M_PI / 180., ph = b * M_PI / 90.;
      Eigen::Vector3d n( std::sin( th ) * std::cos( ph ), std::sin( th ) * std::sin( ph ), std::cos( th ) );
      const double    r = acoustic( D, n ).determinant() / acoustic( Del, n ).determinant();
      if ( r < m ) {
        m = r;
        if ( worst )
          *worst = n;
      }
    }
  return m;
}

int main( int argc, char** argv )
{
  std::vector< double > dp = { 100., 0.3, 0.02, 30., 10., 1. };
  for ( int i = 1; i < argc && i <= 6; ++i )
    dp[i - 1] = std::atof( argv[i] );
  const bool   withCap = argc > 7;
  const double pt      = withCap ? std::atof( argv[7] ) : 0.;
  const double E = dp[0], nu = dp[1], c0 = dp[2], H = dp[5];
  const double sy = E * 1e-3; // yield strain of about 1e-3 for the J2 comparison
  const double sc = ( c0 / E ) / ( 0.02 / 100. ); // strain scale: the amplitudes below are for the default c0 / E

  std::printf( "1  phi = psi = 0 against VONMISES\n" );
  {
    std::vector< double > pdp = { E, nu, sy / std::sqrt( 3. ), 0., 0., 2. * E / 3. / 3. };
    std::vector< double > pvm = { E, nu, sy, 2. * E / 3., 0., 0. };
    Point                 a( "DRUCKERPRAGER", pdp ), b( "VONMISES", pvm );
    double                ds = 0., dc = 0.;
    for ( int k = 0; k < 40; ++k ) {
      Vector6d dE;
      dE << 1e-4, -2e-5, -3e-5, 2.5e-4 * std::cos( 0.2 * k ), 1e-4, -5e-5;
      const Matrix6d Ca = a.step( dE ), Cb = b.step( dE );
      ds = std::max( ds, ( a.s.stress - b.s.stress ).norm() / b.s.stress.norm() );
      dc = std::max( dc, ( Ca - Cb ).norm() / Cb.norm() );
    }
    check( ds < 1e-10, "largest relative stress difference over 40 steps", ds );
    check( dc < 1e-9, "largest relative tangent difference", dc );
  }

  const double et = Materials::DruckerPragerModel::eta( dp[3] ), x = Materials::DruckerPragerModel::xi( dp[3] );
  auto         yieldF = [&]( const Point& P ) {
    const Vector6d& S = P.s.stress;
    const double    p = S.head< 3 >().sum() / 3.;
    Vector6d        s = S;
    s.head< 3 >().array() -= p;
    const double J2 = 0.5 * ( s.head< 3 >().squaredNorm() + 2. * s.tail< 3 >().squaredNorm() );
    return std::sqrt( J2 ) + et * p - x * ( c0 + H * P.sv[0] );
  };

  std::printf( "2, 3  consistency and FD tangent, E = %g nu = %g c0 = %g phi = %g psi = %g H = %g\n", E, nu, c0, dp[3],
               dp[4], H );
  for ( double psi : { dp[3], dp[4] } ) {
    std::vector< double > p = { E, nu, c0, dp[3], psi, H };
    Point                 P( "DRUCKERPRAGER", p );
    double                fmax = 0., tmax = 0.;
    int                   nPlastic = 0;
    std::mt19937          rng( 7 );
    std::normal_distribution< double > N( 0., 1. );
    for ( int k = 0; k < 60; ++k ) {
      Vector6d dE;
      for ( int i = 0; i < 6; ++i )
        dE( i ) = sc * 2e-4 * N( rng );
      dE.head< 3 >().array() -= sc * 4e-5; // net compression
      tmax = std::max( tmax, fdTangentError( P, dE ) );
      const double kappaOld = P.sv[0];
      P.step( dE );
      if ( P.sv[0] > kappaOld ) { // a plastic step
        fmax = std::max( fmax, std::abs( yieldF( P ) ) / ( x * c0 ) );
        nPlastic++;
      }
    }
    char buf[128];
    std::snprintf( buf, sizeof buf, "psi = %4.1f: |f| / (xi c0) after %d plastic steps", psi, nPlastic );
    check( fmax < 1e-10 && nPlastic > 10, buf, fmax );
    std::snprintf( buf, sizeof buf, "psi = %4.1f: tangent vs central FD, worst over a random path", psi );
    check( tmax < 1e-5, buf, tmax );
    if ( psi > 0. && dp[3] > 0. ) { // apex: hydrostatic tension past the cone tip (no apex for phi = 0)
      Point Q( "DRUCKERPRAGER", p );
      Q.step( ( Vector6d() << 0, 0, 0, 0, 0, 0 ).finished() );
      Vector6d dE;
      dE << sc * 4e-4, sc * 4e-4, sc * 4e-4, sc * 1e-6, 0, 0;
      std::snprintf( buf, sizeof buf, "psi = %4.1f: apex return, tangent vs FD", psi );
      check( fdTangentError( Q, dE ) < 1e-5, buf, fdTangentError( Q, dE ) );
      Q.step( dE );
      check( Q.s.stress.tail< 3 >().norm() < 1e-12 * c0, "apex: deviatoric stress vanishes", Q.s.stress.tail< 3 >().norm() );
    }
    if ( psi == 0. && dp[3] > 0. ) { // tension cut-off: hydrostatic tension past the cone tip with non-dilatant flow
      Point Q( "DRUCKERPRAGER", p );
      Vector6d dE;
      dE << sc * 4e-4, sc * 4e-4, sc * 4e-4, sc * 1e-6, 0, 0;
      Matrix6d C;
      const Vector6d S      = Q.trial( dE, &C );
      const double   pApex  = x * c0 / et;
      const double   pError = std::abs( S.head< 3 >().sum() / 3. - pApex ) / pApex;
      std::snprintf( buf, sizeof buf, "psi = 0: cut-off puts p on the apex xi c0 / eta (rel. error)" );
      check( pError < 1e-12, buf, pError );
      check( S.tail< 3 >().norm() < 1e-12 * c0 && ( S.head< 3 >().array() - S( 0 ) ).abs().maxCoeff() < 1e-12 * c0,
             "psi = 0: cut-off stress is hydrostatic", S.tail< 3 >().norm() );
      Matrix6d     Cfd;
      const double hh = 1e-7 * dE.norm();
      for ( int j = 0; j < 6; ++j ) {
        Vector6d e = Vector6d::Zero();
        e( j )     = hh;
        Cfd.col( j ) = ( Q.trial( dE + e ) - Q.trial( dE - e ) ) / ( 2. * hh );
      }
      check( C.norm() == 0. && Cfd.norm() < 1e-6 * E, "psi = 0: cut-off tangent zero, FD agrees (|Cfd| / E)", Cfd.norm() / E );
    }
  }

  std::printf( "4  ellipticity along the interface path: gamma_xz = cos(10) d, eps_zz = -sin(10) d\n" );
  {
    std::vector< double > p   = { E, nu, c0, dp[3], dp[4], H };
    std::vector< double > pel = { E, nu, 1e30, dp[3], dp[4], H };
    Point                 P( "DRUCKERPRAGER", p ), El( "DRUCKERPRAGER", pel );
    const Matrix6d        Del = El.step( Vector6d::Zero() );
    const double          d   = 2e-5;
    double                worstAll = 1e300, worstNormal = 1e300;
    Eigen::Vector3d       nWorst;
    for ( int k = 0; k < 200; ++k ) {
      Vector6d dE = Vector6d::Zero();
      dE( 2 )     = -std::sin( M_PI / 18. ) * d;
      dE( 4 )     = std::cos( M_PI / 18. ) * d;
      const Matrix6d D = P.step( dE );
      if ( P.sv[0] == 0. )
        continue;
      Eigen::Vector3d n;
      const double    r = minDetRatio( D, Del, &n );
      if ( r < worstAll ) {
        worstAll = r;
        nWorst   = n;
      }
      worstNormal = std::min( worstNormal, acoustic( D, Eigen::Vector3d( 0, 0, 1 ) ).determinant() /
                                             acoustic( Del, Eigen::Vector3d( 0, 0, 1 ) ).determinant() );
    }
    std::printf( "  kappa at the end %.3e, sigma_xz %.4f, sigma_zz %.4f\n", P.sv[0], P.s.stress( 4 ), P.s.stress( 2 ) );
    check( worstAll > 0., "min over the path and all directions of det Q / det Q_el", worstAll );
    std::printf( "     reached at n = (%.3f, %.3f, %.3f)\n", nWorst( 0 ), nWorst( 1 ), nWorst( 2 ) );
    check( worstNormal > 0., "min over the path at the interface normal (<Q> of the reduction)", worstNormal );
  }
  if ( withCap ) {
    std::printf( "5  tension cap pt = %g\n", pt );
    std::vector< double > p = { E, nu, c0, dp[3], dp[4], H, pt };
    const double          K = E / ( 3. * ( 1. - 2. * nu ) ), G = E / ( 2. * ( 1. + nu ) );
    auto                  pOf  = []( const Vector6d& S ) { return S.head< 3 >().sum() / 3.; };
    auto                  sqJ2 = [&]( const Vector6d& S ) {
      Vector6d s = S;
      s.head< 3 >().array() -= pOf( S );
      return std::sqrt( 0.5 * s.head< 3 >().squaredNorm() + s.tail< 3 >().squaredNorm() );
    };
    {   // cap only: hydrostatic tension to 1.5 pt with a small shear
      Point    P( "DRUCKERPRAGER", p );
      Vector6d dE;
      dE << 1.5 * pt / ( 3. * K ), 1.5 * pt / ( 3. * K ), 1.5 * pt / ( 3. * K ), 0.05 * c0 / G, 0, 0;
      const double fd = fdTangentError( P, dE );
      P.step( dE );
      check( std::abs( pOf( P.s.stress ) - pt ) < 1e-12 * std::abs( pt ), "cap only: p = pt", std::abs( pOf( P.s.stress ) - pt ) );
      check( std::abs( P.s.stress( 3 ) - G * 0.05 * c0 / G ) < 1e-12 * c0 && P.sv[0] == 0., "cap only: deviator elastic, kappa unchanged",
             std::abs( P.s.stress( 3 ) - 0.05 * c0 ) );
      check( fd < 1e-5, "cap only: tangent vs FD", fd );
    }
    {   // tension to 1.5 pt with a shear well beyond the cone: a corner return (p = pt on the cone), or, if the dilatant
        // flow of the cone return already brings p below the cap, a cone return below the cap; either way the
        // state is admissible (f1 = 0, p <= pt) and kappa grows
      Point    P( "DRUCKERPRAGER", p );
      Vector6d dE;
      dE << 1.5 * pt / ( 3. * K ), 1.5 * pt / ( 3. * K ), 1.5 * pt / ( 3. * K ), 4. * c0 / G, 0, 0;
      const double fd = fdTangentError( P, dE );
      P.step( dE );
      const double f1     = sqJ2( P.s.stress ) + et * pOf( P.s.stress ) - x * ( c0 + H * P.sv[0] );
      const bool   corner = std::abs( pOf( P.s.stress ) - pt ) < 1e-12 * std::abs( pt );
      check( std::abs( f1 ) < 1e-10 * x * c0 && pOf( P.s.stress ) <= pt * ( 1. + 1e-12 ) && P.sv[0] > 0.,
             corner ? "tension + shear: corner (p = pt, on the cone, kappa grows)" : "tension + shear: cone below the cap (f1 = 0, p < pt)",
             std::abs( f1 ) );
      if ( dp[4] == 0. )
        check( corner, "psi = 0: no dilatancy to lower p, so the corner is reached", std::abs( pOf( P.s.stress ) - pt ) );
      check( fd < 1e-5, "tension + shear: tangent vs FD", fd );
    }
    {   // corner with dilatancy: the shear just beyond the cone at p = pt, so the cone return would end above the cap
      Point        P( "DRUCKERPRAGER", p );
      const double s0 = ( x * c0 - et * pt ) * 1.02; // sqrt(J2) slightly beyond the cone at p = pt
      Vector6d     dE;
      dE << 1.5 * pt / ( 3. * K ), 1.5 * pt / ( 3. * K ), 1.5 * pt / ( 3. * K ), s0 / G, 0, 0;
      const double fd = fdTangentError( P, dE );
      P.step( dE );
      const double f1 = sqJ2( P.s.stress ) + et * pOf( P.s.stress ) - x * ( c0 + H * P.sv[0] );
      check( std::abs( pOf( P.s.stress ) - pt ) < 1e-12 * std::abs( pt ) && std::abs( f1 ) < 1e-10 * x * c0 && P.sv[0] > 0.,
             "corner: p = pt, on the cone, kappa grows", std::abs( pOf( P.s.stress ) - pt ) + std::abs( f1 ) );
      check( fd < 1e-5, "corner: tangent vs FD", fd );
    }
    {   // random tensile path: both surfaces respected, tangent consistent
      Point                              P( "DRUCKERPRAGER", p );
      std::mt19937                       rng( 11 );
      std::normal_distribution< double > N( 0., 1. );
      double                             fmax = 0., pmax = -1e300, tmax = 0.;
      int                                nCap = 0;
      for ( int k = 0; k < 80; ++k ) {
        Vector6d dE;
        for ( int i = 0; i < 6; ++i )
          dE( i ) = sc * 2e-4 * N( rng );
        dE.head< 3 >().array() += 0.4 * pt / ( 3. * K ); // net tension
        tmax = std::max( tmax, fdTangentError( P, dE ) );
        P.step( dE );
        const double f1 = sqJ2( P.s.stress ) + et * pOf( P.s.stress ) - x * ( c0 + H * P.sv[0] );
        fmax            = std::max( fmax, f1 / ( x * c0 ) );
        pmax            = std::max( pmax, pOf( P.s.stress ) - pt );
        nCap += std::abs( pOf( P.s.stress ) - pt ) < 1e-12 * std::abs( pt );
      }
      check( fmax < 1e-10, "random path: max f1 / (xi c0) (cone never exceeded)", fmax );
      char buf[96];
      std::snprintf( buf, sizeof buf, "random path: max (p - pt), %d steps on the cap", nCap );
      check( pmax < 1e-12 * std::abs( pt ) && nCap > 5, buf, pmax );
      check( tmax < 1e-5, "random path: tangent vs FD, worst", tmax );
    }
  }
  std::printf( failures ? "%d FAILED\n" : "all passed\n", failures );
  return failures;
}
