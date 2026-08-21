/**
 * Tests for MarmotZInterfaceMaterialHypoElastic -- the two-sided interface
 * kernel that returns the gradient jump g as an EXTERNAL unknown instead of
 * condensing it.
 *
 * The tests pin the kernel to the paper's equations:
 *   - p_z equals (h/4)[t_i] of eq. (23), the term the condensed formulation
 *     is allowed to drop and this one is not;
 *   - K_zz equals (h/4)<Q> of eq. (13)/(16);
 *   - eliminating g by hand (solve p_z = 0, then Schur-complement) reproduces
 *     MarmotEquilibratedXInterfaceMaterialHypoElastic exactly, i.e. eq. (18)
 *     is recovered wherever it is well posed;
 *   - all four tangent blocks are consistent with finite differences.
 */
#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotZInterfaceMaterialHypoElastic.h"

#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

namespace {

  using Material = MarmotZInterfaceMaterialHypoElastic;

  constexpr int nX = Material::nX; // 21
  constexpr int nZ = Material::nZ; // 3

  using VectorX  = Eigen::Matrix< double, nX, 1 >;
  using VectorZ  = Eigen::Matrix< double, nZ, 1 >;
  using MatrixXX = Eigen::Matrix< double, nX, nX, Eigen::RowMajor >;
  using MatrixXZ = Eigen::Matrix< double, nX, nZ, Eigen::RowMajor >;
  using MatrixZX = Eigen::Matrix< double, nZ, nX, Eigen::RowMajor >;
  using MatrixZZ = Eigen::Matrix< double, nZ, nZ, Eigen::RowMajor >;

  const double          hLayer     = 0.1;
  const Eigen::Vector3d normal     = Eigen::Vector3d( 0.0, 0.0, 1.0 );
  const Eigen::Vector3d separation = Eigen::Vector3d( 0.0, 0.0, hLayer );

  std::vector< double > elasticProps()
  {
    return { 2.0e5, 0.3, hLayer };
  }
  // fy = 5, H = 0.1: the near-perfect-plastic corner of eq. (48), which is
  // where the condensed element loses the acoustic tensor.
  std::vector< double > vonMisesProps()
  {
    return { 4.0e5, 0.3, hLayer, 5.0, 0.1, 0.0, 0.0, 0.0 };
  }

  struct Result {
    VectorX  pX;
    VectorZ  pZ;
    MatrixXX K_xx;
    MatrixXZ K_xz;
    MatrixZX K_zx;
    MatrixZZ K_zz;
  };

  /** Evaluate on a scratch copy of stateIn, leaving stateIn untouched. */
  Result evaluate( Material&                    material,
                   const std::vector< double >& stateIn,
                   const VectorX&               dX,
                   const VectorZ&               dz,
                   std::vector< double >*       stateOut = nullptr )
  {
    std::vector< double > sv = stateIn;

    Eigen::Matrix< double, 6, 1 > dU;
    dU.setZero();
    dU.segment< 3 >( 0 ) = dX.segment< 3 >( 0 ); // u+ ; u- = 0 gives w = dX(0..2)
    Eigen::Matrix< double, 18, 1 > dSurf;
    dSurf.segment< 9 >( 0 ) = dX.segment< 9 >( 3 );
    dSurf.segment< 9 >( 9 ) = dX.segment< 9 >( 12 );

    Result                  r;
    Material::State         state{ r.pX.data(), r.pX.data() + 3, r.pX.data() + 12, r.pZ.data(), sv.data() };
    Material::Tangents      tangents{ r.K_xx.data(), r.K_xz.data(), r.K_zx.data(), r.K_zz.data() };
    Material::Deformation   deformation{ dU.data(), dSurf.data(), dz.data(), normal.data(), separation.data() };
    Material::TimeIncrement ti{ 0.0, 1.0 };

    material.computeStress( state, tangents, deformation, ti );

    if ( stateOut )
      *stateOut = sv;
    return r;
  }

  void testTangentAgainstFiniteDifferences( const std::string& name, const std::vector< double >& props )
  {
    Material material( name, props.data(), (int)props.size(), 0 );
    material.setCharacteristicElementLength( 1.0 );
    material.setGradientJumpRegularization( 0.0 );

    std::vector< double > sv( material.getNumberOfRequiredStateVars(), 0.0 );
    material.initializeYourself( sv.data(), (int)sv.size() );

    VectorX dX;
    VectorZ dz;
    for ( int i = 0; i < nX; ++i )
      dX( i ) = 1.0e-4 * std::sin( 1.0 + 0.7 * i );
    for ( int i = 0; i < nZ; ++i )
      dz( i ) = 1.0e-4 * std::cos( 0.3 + 1.1 * i );

    const Result base = evaluate( material, sv, dX, dz );

    const double eps      = 1.0e-9;
    double       maxErrXX = 0.0, maxErrXZ = 0.0, maxErrZX = 0.0, maxErrZZ = 0.0;

    for ( int c = 0; c < nX; ++c ) {
      VectorX dXp = dX;
      dXp( c ) += eps;
      const Result p = evaluate( material, sv, dXp, dz );
      for ( int r = 0; r < nX; ++r )
        maxErrXX = std::max( maxErrXX, std::abs( ( p.pX( r ) - base.pX( r ) ) / eps - base.K_xx( r, c ) ) );
      for ( int r = 0; r < nZ; ++r )
        maxErrZX = std::max( maxErrZX, std::abs( ( p.pZ( r ) - base.pZ( r ) ) / eps - base.K_zx( r, c ) ) );
    }

    for ( int c = 0; c < nZ; ++c ) {
      VectorZ dzp = dz;
      dzp( c ) += eps;
      const Result p = evaluate( material, sv, dX, dzp );
      for ( int r = 0; r < nX; ++r )
        maxErrXZ = std::max( maxErrXZ, std::abs( ( p.pX( r ) - base.pX( r ) ) / eps - base.K_xz( r, c ) ) );
      for ( int r = 0; r < nZ; ++r )
        maxErrZZ = std::max( maxErrZZ, std::abs( ( p.pZ( r ) - base.pZ( r ) ) / eps - base.K_zz( r, c ) ) );
    }

    const double scale = base.K_xx.cwiseAbs().maxCoeff();
    const double tol   = 1.0e-5 * scale;

    std::cout << "  [" << name << "] FD tangent errors (rel to " << scale << "): " << maxErrXX / scale << " "
              << maxErrXZ / scale << " " << maxErrZX / scale << " " << maxErrZZ / scale << std::endl;

    throwExceptionOnFailure( maxErrXX < tol, "K_xx inconsistent with finite differences." );
    throwExceptionOnFailure( maxErrXZ < tol, "K_xz inconsistent with finite differences." );
    throwExceptionOnFailure( maxErrZX < tol, "K_zx inconsistent with finite differences." );
    throwExceptionOnFailure( maxErrZZ < tol, "K_zz inconsistent with finite differences." );
  }

  /**
   * p_z must be exactly (h/4)[t_i] of eq. (23), with no regularization, and
   * K_zz must be exactly (h/4)<Q>. Both are checked against the definition
   * rebuilt from the material's own outputs.
   */
  void testResidualIsTheTractionImbalance()
  {
    const auto props = elasticProps();
    Material   material( "LINEARELASTIC", props.data(), (int)props.size(), 0 );
    material.setCharacteristicElementLength( 1.0 );
    material.setGradientJumpRegularization( 0.0 );

    std::vector< double > sv( material.getNumberOfRequiredStateVars(), 0.0 );
    material.initializeYourself( sv.data(), (int)sv.size() );

    VectorX dX = VectorX::Zero();
    VectorZ dz;
    dz << 1.0e-3, -2.0e-3, 5.0e-4;
    dX( 0 ) = 1.0e-4; // a normal jump as well, so g is not the only driver

    std::vector< double > svOut;
    const Result          r = evaluate( material, sv, dX, dz, &svOut );

    // The kernel commits both faces' Cauchy stress; rebuild [t] = (sigma+ - sigma-) n.
    const auto plusView  = material.getStateView( "stressPlus", svOut.data() );
    const auto minusView = material.getStateView( "stressMinus", svOut.data() );
    const auto toStress  = []( const double* v ) {
      Eigen::Matrix3d s;
      s << v[0], v[3], v[4], v[3], v[1], v[5], v[4], v[5], v[2];
      return s;
    };
    const Eigen::Vector3d tPlus  = toStress( plusView.stateLocation ) * normal;
    const Eigen::Vector3d tMinus = toStress( minusView.stateLocation ) * normal;

    const Eigen::Vector3d expected = 0.25 * hLayer * ( tPlus - tMinus );

    std::cout << "  p_z            = " << r.pZ.transpose() << std::endl;
    std::cout << "  (h/4)[t] eq 23 = " << expected.transpose() << std::endl;

    throwExceptionOnFailure( ( r.pZ - expected ).norm() <= 1.0e-10 * std::max( 1.0, expected.norm() ),
                             "p_z must equal the (h/4)[t_i] term of eq. (23)." );

    // K_zz = (h/4) <Q>. For an isotropic elastic layer <Q> = Q^e, so
    // K_zz = (h/4)( mu I + (lambda+mu) n (x) n ).
    const double    E = props[0], nu = props[1];
    const double    mu       = E / ( 2.0 * ( 1.0 + nu ) );
    const double    lambda   = E * nu / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
    Eigen::Matrix3d QElastic = mu * Eigen::Matrix3d::Identity() + ( lambda + mu ) * ( normal * normal.transpose() );
    const Eigen::Matrix3d expectedKzz = 0.25 * hLayer * QElastic;

    throwExceptionOnFailure( ( Eigen::Matrix3d( r.K_zz ) - expectedKzz ).norm() <= 1.0e-9 * expectedKzz.norm(),
                             "K_zz must equal (h/4)<Q> of eq. (13)/(16)." );
  }

  /**
   * The regularization is a spring on the z INCREMENT: it adds
   * (h/4) zeta Q^e dz to p_z and (h/4) zeta Q^e to K_zz, and NOTHING else.
   */
  void testRegularizationIsIncrementalAndLocalToZ()
  {
    const auto props = elasticProps();
    Material   plain( "LINEARELASTIC", props.data(), (int)props.size(), 0 );
    Material   regularized( "LINEARELASTIC", props.data(), (int)props.size(), 0 );
    plain.setGradientJumpRegularization( 0.0 );
    const double zeta = 1.0e-3;
    regularized.setGradientJumpRegularization( zeta );

    std::vector< double > sv( plain.getNumberOfRequiredStateVars(), 0.0 );
    plain.initializeYourself( sv.data(), (int)sv.size() );

    VectorX dX = VectorX::Zero();
    dX( 0 )    = 3.0e-4;
    VectorZ dz;
    dz << 1.0e-3, -2.0e-3, 5.0e-4;

    const Result a = evaluate( plain, sv, dX, dz );
    const Result b = evaluate( regularized, sv, dX, dz );

    const double          E = props[0], nu = props[1];
    const double          mu     = E / ( 2.0 * ( 1.0 + nu ) );
    const double          lambda = E * nu / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
    const Eigen::Matrix3d QReg   = 0.25 * hLayer * zeta *
                                 ( mu * Eigen::Matrix3d::Identity() +
                                   ( lambda + mu ) * ( normal * normal.transpose() ) );

    throwExceptionOnFailure( ( b.pX - a.pX ).norm() <= 1.0e-12 * std::max( 1.0, a.pX.norm() ),
                             "the regularization must not touch p_X." );
    throwExceptionOnFailure( ( b.K_xx - a.K_xx ).norm() <= 1.0e-12 * a.K_xx.norm() &&
                               ( b.K_xz - a.K_xz ).norm() <= 1.0e-12 * std::max( 1.0, a.K_xz.norm() ) &&
                               ( b.K_zx - a.K_zx ).norm() <= 1.0e-12 * std::max( 1.0, a.K_zx.norm() ),
                             "the regularization must not touch the coupling blocks." );
    throwExceptionOnFailure( ( Eigen::Vector3d( b.pZ - a.pZ ) - QReg * dz ).norm() <=
                               1.0e-10 * std::max( 1.0, ( QReg * dz ).norm() ),
                             "the regularization must add (h/4) zeta Q^e dz to p_z." );
    throwExceptionOnFailure( ( Eigen::Matrix3d( b.K_zz - a.K_zz ) - QReg ).norm() <= 1.0e-10 * QReg.norm(),
                             "the regularization must add (h/4) zeta Q^e to K_zz." );

    // Zero increment => zero regularization contribution: the spring never
    // pushes g back to zero, it only resists CHANGING it within the step.
    const Result c = evaluate( regularized, sv, dX, VectorZ::Zero() );
    const Result d = evaluate( plain, sv, dX, VectorZ::Zero() );
    throwExceptionOnFailure( ( c.pZ - d.pZ ).norm() <= 1.0e-14 * std::max( 1.0, d.pZ.norm() ),
                             "the regularization must vanish at zero z increment." );
  }

  /**
   * The decisive consistency check. With zeta = 0, solve p_z(z) = 0 by Newton
   * (which is eq. (45) driven by K_zz = (h/4)<Q>, i.e. eq. (13)), then form
   * the Schur complement of eq. (47) by hand. The result must reproduce
   * MarmotEquilibratedXInterfaceMaterialHypoElastic bit for bit -- i.e.
   * eq. (18) is recovered wherever it is well posed, and the difference
   * between the two elements is confined to HOW g is determined.
   */
  void testEquivalenceWithCondensedKernel( const std::string& name, const std::vector< double >& props )
  {
    Material zMaterial( name, props.data(), (int)props.size(), 0 );
    zMaterial.setCharacteristicElementLength( 1.0 );
    zMaterial.setGradientJumpRegularization( 0.0 );

    MarmotEquilibratedXInterfaceMaterialHypoElastic yMaterial( name, props.data(), (int)props.size(), 0 );
    yMaterial.setCharacteristicElementLength( 1.0 );

    std::vector< double > svZ( zMaterial.getNumberOfRequiredStateVars(), 0.0 );
    zMaterial.initializeYourself( svZ.data(), (int)svZ.size() );
    std::vector< double > svY( yMaterial.getNumberOfRequiredStateVars(), 0.0 );
    yMaterial.initializeYourself( svY.data(), (int)svY.size() );

    VectorX dX;
    for ( int i = 0; i < nX; ++i )
      dX( i ) = 2.0e-5 * std::sin( 0.4 + 0.9 * i );

    // --- local Newton on p_z = 0, exactly the statement the condensed kernel solves ---
    VectorZ z = VectorZ::Zero();
    Result  r;
    for ( int it = 0; it < 50; ++it ) {
      r = evaluate( zMaterial, svZ, dX, z );
      if ( r.pZ.norm() < 1.0e-14 )
        break;
      z -= Eigen::Matrix3d( r.K_zz ).fullPivLu().solve( Eigen::Vector3d( r.pZ ) );
    }
    r = evaluate( zMaterial, svZ, dX, z );
    throwExceptionOnFailure( r.pZ.norm() < 1.0e-10, "local Newton on p_z did not converge for this test state." );

    const MatrixXX schur = r.K_xx - r.K_xz * Eigen::Matrix3d( r.K_zz ).fullPivLu().solve( MatrixZX( r.K_zx ) );

    // --- the condensed kernel, same input ---
    Eigen::Vector3d                                force;
    Eigen::Matrix< double, 9, 1 >                  sPlus, sMinus;
    Eigen::Matrix< double, 3, 3, Eigen::RowMajor > Q_ww;
    Eigen::Matrix< double, 3, 9, Eigen::RowMajor > Q_wAp, Q_wAm;
    Eigen::Matrix< double, 9, 3, Eigen::RowMajor > Q_Apw, Q_Amw;
    Eigen::Matrix< double, 9, 9, Eigen::RowMajor > Q_ApAp, Q_ApAm, Q_AmAp, Q_AmAm;
    force.setZero();
    sPlus.setZero();
    sMinus.setZero();

    Eigen::Matrix< double, 6, 1 > dU;
    dU.setZero();
    dU.segment< 3 >( 0 ) = dX.segment< 3 >( 0 );
    Eigen::Matrix< double, 18, 1 > dSurf;
    dSurf.segment< 9 >( 0 ) = dX.segment< 9 >( 3 );
    dSurf.segment< 9 >( 9 ) = dX.segment< 9 >( 12 );

    MarmotEquilibratedXInterfaceMaterialHypoElastic::State         materialState{ force.data(),
                                                                          sPlus.data(),
                                                                          sMinus.data(),
                                                                          svY.data() };
    MarmotEquilibratedXInterfaceMaterialHypoElastic::Tangents      materialTangents{ Q_ww.data(),
                                                                                Q_wAp.data(),
                                                                                Q_wAm.data(),
                                                                                Q_Apw.data(),
                                                                                Q_ApAp.data(),
                                                                                Q_ApAm.data(),
                                                                                Q_Amw.data(),
                                                                                Q_AmAp.data(),
                                                                                Q_AmAm.data() };
    MarmotEquilibratedXInterfaceMaterialHypoElastic::Deformation   materialDeformation{ dU.data(),
                                                                                      dSurf.data(),
                                                                                      normal.data(),
                                                                                      separation.data() };
    MarmotEquilibratedXInterfaceMaterialHypoElastic::TimeIncrement materialTime{ 0.0, 1.0 };
    yMaterial.computeStress( materialState, materialTangents, materialDeformation, materialTime );

    MatrixXX condensed;
    condensed.setZero();
    condensed.block< 3, 3 >( 0, 0 )   = Q_ww;
    condensed.block< 3, 9 >( 0, 3 )   = Q_wAp;
    condensed.block< 3, 9 >( 0, 12 )  = Q_wAm;
    condensed.block< 9, 3 >( 3, 0 )   = Q_Apw;
    condensed.block< 9, 9 >( 3, 3 )   = Q_ApAp;
    condensed.block< 9, 9 >( 3, 12 )  = Q_ApAm;
    condensed.block< 9, 3 >( 12, 0 )  = Q_Amw;
    condensed.block< 9, 9 >( 12, 3 )  = Q_AmAp;
    condensed.block< 9, 9 >( 12, 12 ) = Q_AmAm;

    VectorX pCondensed;
    pCondensed.segment< 3 >( 0 )  = force;
    pCondensed.segment< 9 >( 3 )  = sPlus;
    pCondensed.segment< 9 >( 12 ) = sMinus;

    const double stressError  = ( r.pX - pCondensed ).norm() / std::max( 1.0, pCondensed.norm() );
    const double tangentError = ( schur - condensed ).norm() / condensed.norm();

    std::cout << "  [" << name << "] vs condensed kernel: |dp| = " << stressError << ", |dK| = " << tangentError
              << std::endl;

    throwExceptionOnFailure( stressError < 1.0e-9, "generalized stress must match the condensed kernel once p_z = 0." );
    throwExceptionOnFailure( tangentError < 1.0e-7,
                             "the Schur complement must reproduce eq. (18)'s condensed tangent." );
  }

  /**
   * Near-perfect plasticity: both faces yielding alike drives <Q> singular
   * (eq. 40/48). The condensed kernel must invert it; this kernel must not.
   * We check that K_zz becomes numerically singular in that state, and that
   * the regularization restores an invertible block without perturbing the
   * traction by more than O(zeta).
   */
  void testDegenerateAcousticTensorStaysUsable()
  {
    const auto props = vonMisesProps();

    Material material( "VONMISES", props.data(), (int)props.size(), 0 );
    material.setCharacteristicElementLength( 1.0 );
    material.setGradientJumpRegularization( 0.0 );

    std::vector< double > sv( material.getNumberOfRequiredStateVars(), 0.0 );
    material.initializeYourself( sv.data(), (int)sv.size() );

    // Drive both faces well past yield in identical simple shear on the
    // interface plane: dA+ = dA- = gamma e_x (x) e_x is not enough, so we
    // shear through the layer with a large tangential jump.
    VectorX               dX    = VectorX::Zero();
    VectorZ               z     = VectorZ::Zero();
    std::vector< double > state = sv;
    for ( int step = 0; step < 40; ++step ) {
      dX.setZero();
      dX( 0 ) = 2.0e-4; // tangential jump w_x, drives g_x through gbar = w/ell
      std::vector< double > next;
      evaluate( material, state, dX, z, &next );
      state = next;
    }

    const Result r = evaluate( material, state, dX, z );

    Eigen::JacobiSVD< Eigen::Matrix3d > svd( Eigen::Matrix3d( r.K_zz ) );
    const Eigen::Vector3d               s            = svd.singularValues();
    const double                        conditioning = s( 2 ) / s( 0 );
    std::cout << "  <Q> conditioning after yielding: sigma_min/sigma_max = " << conditioning << std::endl;

    Material regularized( "VONMISES", props.data(), (int)props.size(), 0 );
    regularized.setCharacteristicElementLength( 1.0 );
    regularized.setGradientJumpRegularization( 1.0e-6 );
    const Result rr = evaluate( regularized, state, dX, z );

    Eigen::JacobiSVD< Eigen::Matrix3d > svdReg( Eigen::Matrix3d( rr.K_zz ) );
    const Eigen::Vector3d               sReg = svdReg.singularValues();
    std::cout << "  with zeta = 1e-6:                sigma_min/sigma_max = " << sReg( 2 ) / sReg( 0 ) << std::endl;

    throwExceptionOnFailure( sReg( 2 ) > 0.0, "the regularized K_zz must be non-singular." );
    throwExceptionOnFailure( sReg( 2 ) / sReg( 0 ) >= conditioning,
                             "the regularization must not worsen the conditioning of K_zz." );
    // The regularization only ever acts through the z increment, which is zero here.
    throwExceptionOnFailure( ( rr.pZ - r.pZ ).norm() <= 1.0e-12 * std::max( 1.0, r.pZ.norm() ),
                             "the regularization must leave the residual untouched at dz = 0." );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{
    [&]() { testTangentAgainstFiniteDifferences( "LINEARELASTIC", elasticProps() ); },
    [&]() { testTangentAgainstFiniteDifferences( "VONMISES", vonMisesProps() ); },
    [&]() { testResidualIsTheTractionImbalance(); },
    [&]() { testRegularizationIsIncrementalAndLocalToZ(); },
    [&]() { testEquivalenceWithCondensedKernel( "LINEARELASTIC", elasticProps() ); },
    [&]() { testEquivalenceWithCondensedKernel( "VONMISES", vonMisesProps() ); },
    [&]() { testDegenerateAcousticTensorStaysUsable(); },
  };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
