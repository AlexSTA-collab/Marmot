#include "Marmot/MarmotPartialMixedInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

namespace {

  using Material = MarmotPartialMixedInterfaceMaterialHypoElastic;

  // x = [vec(Abar)(9); vec(DeltaA)(9); g(3)] = 21 ; y = [sAbar(9); sDeltaA(9); tMat(3)] = 21
  struct KernelEval {
    Eigen::Matrix< double, 21, 1 >  y;
    Eigen::Matrix< double, 21, 21 > H; // reduced tangent assembled from blocks
    Eigen::Vector3d                 tPlus;
    Eigen::Vector3d                 tMinus;
  };

  KernelEval runKernel( Material&                             material,
                        const std::vector< double >&          committedState,
                        const Eigen::Matrix< double, 21, 1 >& x,
                        double                                dT = 1.0 )
  {
    std::vector< double > scratch = committedState;

    Eigen::Matrix< double, 9, 1 > dAbar   = x.segment< 9 >( 0 );
    Eigen::Matrix< double, 9, 1 > dDeltaA = x.segment< 9 >( 9 );
    Eigen::Vector3d               dG      = x.segment< 3 >( 18 );
    const Eigen::Vector3d         normal( 0.0, 0.0, 1.0 );

    std::array< double, 9 >   sAbar{}, sDeltaA{};
    std::array< double, 3 >   tMat{};
    std::array< double, 324 > H_aa{}; // 18x18
    std::array< double, 54 >  H_ag{}; // 18x3
    std::array< double, 54 >  H_ga{}; // 3x18
    std::array< double, 9 >   H_gg{}; // 3x3

    Material::KernelInput in{ dAbar.data(), dDeltaA.data(), dG.data(), normal.data() };
    Material::KernelOutput
      out{ sAbar.data(), sDeltaA.data(), tMat.data(), H_aa.data(), H_ag.data(), H_ga.data(), H_gg.data() };
    Material::TimeIncrement ti{ 0.0, dT };

    material.computeMixedKernel( scratch.data(), in, out, ti );

    KernelEval e;
    e.y.segment< 9 >( 0 )  = Eigen::Map< Eigen::Matrix< double, 9, 1 > >( sAbar.data() );
    e.y.segment< 9 >( 9 )  = Eigen::Map< Eigen::Matrix< double, 9, 1 > >( sDeltaA.data() );
    e.y.segment< 3 >( 18 ) = Eigen::Map< Eigen::Matrix< double, 3, 1 > >( tMat.data() );

    e.H.setZero();
    e.H.block< 18, 18 >( 0, 0 ) = Eigen::Map< Eigen::Matrix< double, 18, 18, Eigen::RowMajor > >( H_aa.data() );
    e.H.block< 18, 3 >( 0, 18 ) = Eigen::Map< Eigen::Matrix< double, 18, 3, Eigen::RowMajor > >( H_ag.data() );
    e.H.block< 3, 18 >( 18, 0 ) = Eigen::Map< Eigen::Matrix< double, 3, 18, Eigen::RowMajor > >( H_ga.data() );
    e.H.block< 3, 3 >( 18, 18 ) = Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( H_gg.data() );

    // recover t+/t- from the committed sigma+/- for the equilibrium check
    Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > sp(
      material.stateLayout.getPtr( scratch.data(), "stressPlus" ) );
    Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > sm(
      material.stateLayout.getPtr( scratch.data(), "stressMinus" ) );
    e.tPlus  = sp * normal;
    e.tMinus = sm * normal;
    return e;
  }

  Eigen::Matrix< double, 21, 21 > fdTangent( Material&                             material,
                                             const std::vector< double >&          committedState,
                                             const Eigen::Matrix< double, 21, 1 >& x0,
                                             double                                rel = 1e-7 )
  {
    Eigen::Matrix< double, 21, 21 > J;
    for ( int k = 0; k < 21; k++ ) {
      const double                   eps = rel * std::max( 1.0, std::abs( x0( k ) ) );
      Eigen::Matrix< double, 21, 1 > xp = x0, xm = x0;
      xp( k ) += eps;
      xm( k ) -= eps;
      const auto yp = runKernel( material, committedState, xp ).y;
      const auto ym = runKernel( material, committedState, xm ).y;
      J.col( k )    = ( yp - ym ) / ( 2.0 * eps );
    }
    return J;
  }

} // namespace

void TestPartialMixedKernelElasticEquilibriumAndTangent()
{
  std::cout << "\n--- TestPartialMixedKernelElasticEquilibriumAndTangent ---\n";

  const std::vector< double > props = { 8.0e4, 0.3, 0.1 }; // E, nu, h  (LINEARELASTIC base)
  Material                    material( "LINEARELASTIC", props.data(), static_cast< int >( props.size() ), 1 );

  std::vector< double > committed( material.getNumberOfRequiredStateVars(), 0.0 );
  material.initializeYourself( committed.data(), static_cast< int >( committed.size() ) );

  Eigen::Matrix< double, 21, 1 > x;
  x.setZero();
  // generic small increment
  for ( int k = 0; k < 21; k++ )
    x( k ) = 1.0e-4 * std::sin( 0.7 * k + 0.4 );

  const auto base = runKernel( material, committed, x );

  const double equilibriumDefect = ( base.tPlus - base.tMinus ).norm() /
                                   std::max( 1.0, std::max( base.tPlus.norm(), base.tMinus.norm() ) );
  std::cout << "elastic |t+ - t-| (relative) = " << equilibriumDefect << "\n";
  throwExceptionOnFailure( equilibriumDefect < 1e-8, "local traction equilibrium t+=t- not satisfied (elastic)." );

  const auto   Jfd = fdTangent( material, committed, x );
  const double err = ( base.H - Jfd ).norm() / std::max( 1.0, Jfd.norm() );
  std::cout << "elastic reduced-tangent FD error = " << err << "\n";
  throwExceptionOnFailure( base.H.allFinite(), "elastic reduced tangent contains nan/inf." );
  throwExceptionOnFailure( err < 1e-5, "elastic reduced tangent does not match finite difference." );
}

void TestPartialMixedKernelPlasticTangent()
{
  std::cout << "\n--- TestPartialMixedKernelPlasticTangent ---\n";

  // E, nu, h, then VONMISES params (yield, ... hardening)
  const std::vector< double > props = { 210000., 0.3, 0.2, 120., 2100., 20., 20., 2400. };
  Material                    material( "VONMISES", props.data(), static_cast< int >( props.size() ), 1 );

  std::vector< double > committed( material.getNumberOfRequiredStateVars(), 0.0 );
  material.initializeYourself( committed.data(), static_cast< int >( committed.size() ) );

  // Drive into the plastic regime by committing a sequence of increments.
  auto commitStep = [&]( const Eigen::Matrix< double, 21, 1 >& x ) {
    Eigen::Matrix< double, 9, 1 > dAbar   = x.segment< 9 >( 0 );
    Eigen::Matrix< double, 9, 1 > dDeltaA = x.segment< 9 >( 9 );
    Eigen::Vector3d               dG      = x.segment< 3 >( 18 );
    const Eigen::Vector3d         normal( 0.0, 0.0, 1.0 );
    std::array< double, 9 >       sAbar{}, sDeltaA{};
    std::array< double, 3 >       tMat{};
    std::array< double, 324 >     H_aa{};
    std::array< double, 54 >      H_ag{}, H_ga{};
    std::array< double, 9 >       H_gg{};
    Material::KernelInput         in{ dAbar.data(), dDeltaA.data(), dG.data(), normal.data() };
    Material::KernelOutput
      out{ sAbar.data(), sDeltaA.data(), tMat.data(), H_aa.data(), H_ag.data(), H_ga.data(), H_gg.data() };
    Material::TimeIncrement ti{ 0.0, 1.0 };
    material.computeMixedKernel( committed.data(), in, out, ti );
  };

  Eigen::Matrix< double, 21, 1 > xHist;
  xHist.setZero();
  for ( int k = 0; k < 21; k++ )
    xHist( k ) = 3.0e-3 * std::sin( 0.9 * k + 0.2 );
  commitStep( xHist ); // commit a plastic step -> committed now holds a yielded state

  Eigen::Matrix< double, 21, 1 > x;
  x.setZero();
  for ( int k = 0; k < 21; k++ )
    x( k ) = 5.0e-4 * std::sin( 1.3 * k + 0.7 );

  const auto base = runKernel( material, committed, x );

  const double equilibriumDefect = ( base.tPlus - base.tMinus ).norm() /
                                   std::max( 1.0, std::max( base.tPlus.norm(), base.tMinus.norm() ) );
  std::cout << "plastic |t+ - t-| (relative) = " << equilibriumDefect << "\n";
  throwExceptionOnFailure( equilibriumDefect < 1e-8, "local traction equilibrium t+=t- not satisfied (plastic)." );

  const auto   Jfd = fdTangent( material, committed, x );
  const double err = ( base.H - Jfd ).norm() / std::max( 1.0, Jfd.norm() );
  std::cout << "plastic reduced-tangent FD error = " << err << "\n";
  throwExceptionOnFailure( base.H.allFinite(), "plastic reduced tangent contains nan/inf." );
  throwExceptionOnFailure( err < 2e-3, "plastic reduced tangent does not match finite difference." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    TestPartialMixedKernelElasticEquilibriumAndTangent,
    TestPartialMixedKernelPlasticTangent,
  };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
