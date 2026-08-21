/**
 * Tests for MarmotZStabPressureInterfaceMaterialHypoElastic -- the two-face
 * kernel with the gradient jump z external AND two independent pressures.
 *
 * The tests pin the structure the two-pressure form exists for:
 *   - the two blocks that must vanish identically (pbar cannot reach the
 *     traction jump; z cancels from the mean volumetric strain),
 *   - the two couplings the single-pressure form structurally cannot have,
 *     d r_z/d[p] = -(h/4) n  and  d r_[p]/dz = n,
 *   - the solvability of the (z.n, [p]) pair in the double limit K -> inf and
 *     n.Q_dev.n -> 0, where det -> h/4 rather than to zero,
 *   - consistency of every tangent block with finite differences.
 */
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotZStabPressureInterfaceMaterialHypoElastic.h"

#include <Eigen/Dense>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

namespace {

  using M          = MarmotZStabPressureInterfaceMaterialHypoElastic;
  constexpr int nX = M::nX, nZ = M::nZ;

  using VecX = Eigen::Matrix< double, nX, 1 >;
  using VecZ = Eigen::Matrix< double, nZ, 1 >;

  const double          hLayer = 0.1;
  const Eigen::Vector3d nrm( 0.0, 0.0, 1.0 );
  const Eigen::Vector3d sep( 0.0, 0.0, hLayer );

  std::vector< double > elasticProps( double nu = 0.3 )
  {
    return { 2.0e5, nu, hLayer };
  }
  std::vector< double > vmProps( double H )
  {
    return { 4.0e5, 0.3, hLayer, 5.0, H, 0.0, 0.0, 0.0 };
  }

  struct Out {
    VecX                                             pX;
    VecZ                                             pZ;
    double                                           rPm, rPj;
    Eigen::Matrix< double, nX, nX, Eigen::RowMajor > K_xx;
    Eigen::Matrix< double, nX, nZ, Eigen::RowMajor > K_xz;
    Eigen::Matrix< double, nX, 1 >                   K_xpm, K_xpj;
    Eigen::Matrix< double, nZ, nX, Eigen::RowMajor > K_zx;
    Eigen::Matrix< double, nZ, nZ, Eigen::RowMajor > K_zz;
    Eigen::Matrix< double, nZ, 1 >                   K_zpm, K_zpj;
    Eigen::Matrix< double, 1, nX >                   K_pmx, K_pjx;
    Eigen::Matrix< double, 1, nZ >                   K_pmz, K_pjz;
    double                                           K_pmpm, K_pjpj;
  };

  Out evaluate( M& mat, const std::vector< double >& svIn, const VecX& dX, const VecZ& dz, double dpm, double dpj )
  {
    std::vector< double >         sv = svIn;
    Eigen::Matrix< double, 6, 1 > dU;
    dU.setZero();
    dU.segment< 3 >( 0 ) = dX.segment< 3 >( 0 );
    Eigen::Matrix< double, 18, 1 > dA;
    dA.segment< 9 >( 0 ) = dX.segment< 9 >( 3 );
    dA.segment< 9 >( 9 ) = dX.segment< 9 >( 12 );
    Out              o;
    M::Response      r{ o.pX.data(), o.pX.data() + 3, o.pX.data() + 12, o.pZ.data(), &o.rPm, &o.rPj };
    M::Tangents      t{ o.K_xx.data(),
                   o.K_xz.data(),
                   o.K_xpm.data(),
                   o.K_xpj.data(),
                   o.K_zx.data(),
                   o.K_zz.data(),
                   o.K_zpm.data(),
                   o.K_zpj.data(),
                   o.K_pmx.data(),
                   o.K_pmz.data(),
                   &o.K_pmpm,
                   o.K_pjx.data(),
                   o.K_pjz.data(),
                   &o.K_pjpj };
    M::Deformation   d{ dU.data(), dA.data(), dz.data(), nrm.data(), sep.data(), dpm, dpj };
    M::TimeIncrement ti{ 0.0, 1.0 };
    mat.computeStress( sv.data(), r, t, d, ti );
    return o;
  }

  void testStructuralBlocks()
  {
    const auto p = elasticProps();
    M          mat( "LINEARELASTIC", p.data(), (int)p.size(), 0 );
    mat.setGradientJumpRegularization( 0.0 );
    std::vector< double > sv( mat.getNumberOfRequiredStateVars(), 0.0 );
    mat.initializeYourself( sv.data(), (int)sv.size() );

    VecX dX;
    for ( int i = 0; i < nX; ++i )
      dX( i ) = 1e-4 * std::sin( 0.4 + 0.8 * i );
    VecZ dz;
    dz << 1.1e-3, -7.0e-4, 5.0e-4;
    const Out o = evaluate( mat, sv, dX, dz, 3.0, -2.0 );

    std::cout << "  K_zpm (must be 0)      = " << o.K_zpm.transpose() << std::endl;
    std::cout << "  K_pmz (must be 0)      = " << o.K_pmz << std::endl;
    std::cout << "  K_zpj (must be -(h/4)n)= " << o.K_zpj.transpose() << "   expected "
              << ( -0.25 * hLayer * nrm ).transpose() << std::endl;
    std::cout << "  K_pjz (must be n^T)    = " << o.K_pjz << "   expected " << nrm.transpose() << std::endl;

    throwExceptionOnFailure( o.K_zpm.norm() < 1e-14, "pbar must not reach the traction jump." );
    throwExceptionOnFailure( o.K_pmz.norm() < 1e-14, "z must cancel from the MEAN volumetric strain." );
    throwExceptionOnFailure( ( o.K_zpj - Eigen::Vector3d( -0.25 * hLayer * nrm ) ).norm() < 1e-14,
                             "d r_z / d[p] must equal -(h/4) n." );
    throwExceptionOnFailure( ( o.K_pjz.transpose() - nrm ).norm() < 1e-14, "d r_[p] / dz must equal n." );
  }

  /** det of the (z.n, [p]) pair, in the double limit that breaks the other forms. */
  void testNormalPairStaysSolvable()
  {
    std::cout << "  (z.n,[p]) pair: det vs the h/4 floor" << std::endl;
    std::cout << "      nu        H         K          n.Qdev.n        det      det/(h/4)" << std::endl;
    double worst = 1e30;
    for ( double nu : { 0.3, 0.49, 0.4999 } )
      for ( double H : { 1000.0, 1.0, 0.0 } ) {
        auto p = vmProps( H );
        p[1]   = nu;
        M mat( "VONMISES", p.data(), (int)p.size(), 0 );
        mat.setGradientJumpRegularization( 0.0 );
        std::vector< double > sv( mat.getNumberOfRequiredStateVars(), 0.0 );
        mat.initializeYourself( sv.data(), (int)sv.size() );

        // drive both faces well past yield in shear
        VecX dX = VecX::Zero();
        dX( 0 ) = 3.0e-4;
        VecZ dz = VecZ::Zero();
        for ( int s = 0; s < 40; ++s ) {
          std::vector< double > next = sv;
          Out                   tmp;
          {
            std::vector< double >         scratch = sv;
            Eigen::Matrix< double, 6, 1 > dU;
            dU.setZero();
            dU.segment< 3 >( 0 ) = dX.segment< 3 >( 0 );
            Eigen::Matrix< double, 18, 1 > dA;
            dA.segment< 9 >( 0 ) = dX.segment< 9 >( 3 );
            dA.segment< 9 >( 9 ) = dX.segment< 9 >( 12 );
            M::Response    r{ tmp.pX.data(), tmp.pX.data() + 3, tmp.pX.data() + 12, tmp.pZ.data(), &tmp.rPm, &tmp.rPj };
            M::Tangents    t{ tmp.K_xx.data(),
                           tmp.K_xz.data(),
                           tmp.K_xpm.data(),
                           tmp.K_xpj.data(),
                           tmp.K_zx.data(),
                           tmp.K_zz.data(),
                           tmp.K_zpm.data(),
                           tmp.K_zpj.data(),
                           tmp.K_pmx.data(),
                           tmp.K_pmz.data(),
                           &tmp.K_pmpm,
                           tmp.K_pjx.data(),
                           tmp.K_pjz.data(),
                           &tmp.K_pjpj };
            M::Deformation d{ dU.data(), dA.data(), dz.data(), nrm.data(), sep.data(), 0.0, 0.0 };
            M::TimeIncrement ti{ double( s ), 1.0 };
            mat.computeStress( scratch.data(), r, t, d, ti );
            next = scratch;
          }
          sv = next;
        }
        const Out o = evaluate( mat, sv, dX, dz, 0.0, 0.0 );
        // 2x2 in (z.n, [p]) : rows (r_z . n, r_[p]) , cols (z.n, [p])
        const double a11    = nrm.transpose() * Eigen::Matrix3d( o.K_zz ) * nrm; // (h/4) n.Qdev.n
        const double a12    = nrm.dot( o.K_zpj );                                // -(h/4)
        const double a21    = o.K_pjz * nrm;                                     // 1
        const double a22    = o.K_pjpj;                                          // 1/K
        const double det    = a11 * a22 - a12 * a21;
        const double floorv = 0.25 * hLayer;
        worst               = std::min( worst, det / floorv );
        std::printf( "    %7.4f %8.1f %10.3e %12.4e %12.4e %10.4f\n",
                     nu,
                     H,
                     mat.getBulkModulus(),
                     a11 / floorv,
                     det,
                     det / floorv );
      }
    std::fflush( stdout );
    throwExceptionOnFailure( worst > 0.5,
                             "the (z.n,[p]) pair must stay solvable: det >= h/8 even as K -> inf and "
                             "the deviatoric acoustic response degenerates." );
  }

  void testTangentsAgainstFiniteDifferences( const std::string& name, const std::vector< double >& props )
  {
    M mat( name, props.data(), (int)props.size(), 0 );
    mat.setGradientJumpRegularization( 0.0 );
    std::vector< double > sv( mat.getNumberOfRequiredStateVars(), 0.0 );
    mat.initializeYourself( sv.data(), (int)sv.size() );

    VecX dX;
    for ( int i = 0; i < nX; ++i )
      dX( i ) = 1e-4 * std::sin( 1.0 + 0.7 * i );
    VecZ dz;
    for ( int i = 0; i < nZ; ++i )
      dz( i ) = 1e-4 * std::cos( 0.3 + 1.1 * i );
    const double pm = 2.0, pj = -1.5;
    const Out    b = evaluate( mat, sv, dX, dz, pm, pj );

    const double eps = 1e-9;
    double       eXX = 0, eXZ = 0, eZX = 0, eZZ = 0, ePX = 0, ePZ = 0, eXP = 0, eZP = 0;
    for ( int c = 0; c < nX; ++c ) {
      VecX d2 = dX;
      d2( c ) += eps;
      const Out o = evaluate( mat, sv, d2, dz, pm, pj );
      for ( int r = 0; r < nX; ++r )
        eXX = std::max( eXX, std::abs( ( o.pX( r ) - b.pX( r ) ) / eps - b.K_xx( r, c ) ) );
      for ( int r = 0; r < nZ; ++r )
        eZX = std::max( eZX, std::abs( ( o.pZ( r ) - b.pZ( r ) ) / eps - b.K_zx( r, c ) ) );
      ePX = std::max( ePX, std::abs( ( o.rPm - b.rPm ) / eps - b.K_pmx( c ) ) );
      ePX = std::max( ePX, std::abs( ( o.rPj - b.rPj ) / eps - b.K_pjx( c ) ) );
    }
    for ( int c = 0; c < nZ; ++c ) {
      VecZ d2 = dz;
      d2( c ) += eps;
      const Out o = evaluate( mat, sv, dX, d2, pm, pj );
      for ( int r = 0; r < nX; ++r )
        eXZ = std::max( eXZ, std::abs( ( o.pX( r ) - b.pX( r ) ) / eps - b.K_xz( r, c ) ) );
      for ( int r = 0; r < nZ; ++r )
        eZZ = std::max( eZZ, std::abs( ( o.pZ( r ) - b.pZ( r ) ) / eps - b.K_zz( r, c ) ) );
      ePZ = std::max( ePZ, std::abs( ( o.rPm - b.rPm ) / eps - b.K_pmz( c ) ) );
      ePZ = std::max( ePZ, std::abs( ( o.rPj - b.rPj ) / eps - b.K_pjz( c ) ) );
    }
    {
      const Out o = evaluate( mat, sv, dX, dz, pm + eps, pj );
      for ( int r = 0; r < nX; ++r )
        eXP = std::max( eXP, std::abs( ( o.pX( r ) - b.pX( r ) ) / eps - b.K_xpm( r ) ) );
      for ( int r = 0; r < nZ; ++r )
        eZP = std::max( eZP, std::abs( ( o.pZ( r ) - b.pZ( r ) ) / eps - b.K_zpm( r ) ) );
    }
    {
      const Out o = evaluate( mat, sv, dX, dz, pm, pj + eps );
      for ( int r = 0; r < nX; ++r )
        eXP = std::max( eXP, std::abs( ( o.pX( r ) - b.pX( r ) ) / eps - b.K_xpj( r ) ) );
      for ( int r = 0; r < nZ; ++r )
        eZP = std::max( eZP, std::abs( ( o.pZ( r ) - b.pZ( r ) ) / eps - b.K_zpj( r ) ) );
    }

    const double sc = b.K_xx.cwiseAbs().maxCoeff();
    std::cout << "  [" << name << "] FD rel errs: xx=" << eXX / sc << " xz=" << eXZ / sc << " zx=" << eZX / sc
              << " zz=" << eZZ / sc << " xp=" << eXP / sc << " zp=" << eZP / sc << " px=" << ePX << " pz=" << ePZ
              << std::endl;
    throwExceptionOnFailure( eXX < 1e-5 * sc && eXZ < 1e-5 * sc && eZX < 1e-5 * sc && eZZ < 1e-5 * sc,
                             "displacement/z tangent blocks inconsistent with finite differences." );
    throwExceptionOnFailure( eXP < 1e-5 * sc && eZP < 1e-5 * sc, "pressure coupling blocks inconsistent." );
    throwExceptionOnFailure( ePX < 1e-6 && ePZ < 1e-6, "volumetric residual gradients inconsistent." );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{
    [&]() { testStructuralBlocks(); },
    [&]() { testTangentsAgainstFiniteDifferences( "LINEARELASTIC", elasticProps() ); },
    [&]() { testTangentsAgainstFiniteDifferences( "VONMISES", vmProps( 0.1 ) ); },
    [&]() { testNormalPairStaysSolvable(); },
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
