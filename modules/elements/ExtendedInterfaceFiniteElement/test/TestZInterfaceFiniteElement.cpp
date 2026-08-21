/**
 * Tests for ZInterfaceFiniteElement (ZIQUAD4) -- the extended interface
 * element that carries the gradient jump g_k of eq. (12) as a nodal field
 * instead of condensing it at the quadrature point.
 *
 * The decisive test is testEquivalenceWithYIQUAD4: with the regularization
 * off, solving R_g = 0 for the nodal g and then condensing the 36x36 by hand
 * must reproduce YIQUAD4's residual and tangent exactly. With four Q1 nodes
 * and four stations the nodal-to-station map is invertible, so the weak
 * statement R_g = 0 is equivalent to pointwise traction continuity -- the new
 * element gives up nothing where the old one is well posed.
 */
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/YInterfaceFiniteElement.h"
#include "Marmot/ZInterfaceFiniteElement.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Elements;

namespace {

  using E                = ZInterfaceFiniteElement;
  using Y                = YInterfaceFiniteElement< 3, 8 >;
  constexpr int    nDofE = E::sizeLoadVector; // 36
  constexpr int    nDofY = 24;
  constexpr double hMesh = 0.1;               // top face z = hMesh, equal to the material h

  static std::array< double, 24 > flatCoordinates = {
    -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0,
    -0.5, -0.5, 0.1, 0.5, -0.5, 0.1, 0.5, 0.5, 0.1, -0.5, 0.5, 0.1,
  };
  static const std::array< std::array< double, 2 >, 4 > flatXY = {
    { { -0.5, -0.5 }, { 0.5, -0.5 }, { 0.5, 0.5 }, { -0.5, 0.5 } } };

  std::vector< double > elasticProps()
  {
    return { 2.0e5, 0.3, hMesh };
  }
  std::vector< double > vonMisesProps()
  {
    return { 4.0e5, 0.3, hMesh, 5.0, 0.1, 0.0, 0.0, 0.0 };
  }

  /** elementProperties = { out-of-plane thickness, zeta }. */
  /** A node pair separated OBLIQUELY: the top face is the bottom face rigidly
   * translated by d = (0.03, 0.02, 0.13), so the two faces stay congruent (as
   * the whole element family requires -- one gradN serves both) while
   * ell = 0.13 != h = 0.1 AND d_tau = (0.03, 0.02, 0) != 0. This exercises
   * both corrections of eq. (44) at once. */
  static std::array< double, 24 > obliqueCoordinates = {
    -0.5,  -0.5,  0.0,  0.5,  -0.5,  0.0,  0.5,  0.5,  0.0,  -0.5,  0.5,  0.0,
    -0.47, -0.48, 0.13, 0.53, -0.48, 0.13, 0.53, 0.52, 0.13, -0.47, 0.52, 0.13,
  };

  /** The same oblique connector, but on a WARPED (non-planar) bilinear quad:
   * the bottom face corners sit at four different heights, so the normal, the
   * tangent projector, ell and d_tau all vary from station to station. The top
   * face is the SAME warped surface rigidly translated by d, so the two faces
   * stay congruent, as the element family requires. */
  static std::array< double, 24 > warpedCoordinates = {
    -0.5,  -0.5,  0.00, 0.5,  -0.5,  0.06, 0.5,  0.5,  0.00, -0.5,  0.5,  0.06,
    -0.47, -0.48, 0.13, 0.53, -0.48, 0.19, 0.53, 0.52, 0.13, -0.47, 0.52, 0.19,
  };

  /** Affine field u = F x sampled on the given node positions, with a constant
   * nodal gradient jump g0. */
  Eigen::Matrix< double, nDofE, 1 > affineStateOn( const double*          coords,
                                                   const Eigen::Matrix3d& F,
                                                   const Eigen::Vector3d& g0 )
  {
    Eigen::Matrix< double, nDofE, 1 > q;
    q.setZero();
    for ( int A = 0; A < 4; A++ ) {
      const Eigen::Vector3d xB( coords[3 * A], coords[3 * A + 1], coords[3 * A + 2] );
      const Eigen::Vector3d xT( coords[12 + 3 * A], coords[12 + 3 * A + 1], coords[12 + 3 * A + 2] );
      q.segment< 3 >( E::offD + 3 * A )                = F * xB;
      q.segment< 3 >( E::offD + E::nSideDofU + 3 * A ) = F * xT;
      q.segment< 3 >( E::offG + 3 * A )                = g0;
    }
    return q;
  }

  Eigen::Matrix3d isotropicStress( const Eigen::Matrix3d& gradient, double E, double nu )
  {
    const double          mu     = E / ( 2.0 * ( 1.0 + nu ) );
    const double          lambda = E * nu / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
    const Eigen::Matrix3d eps    = 0.5 * ( gradient + gradient.transpose() );
    return lambda * eps.trace() * Eigen::Matrix3d::Identity() + 2.0 * mu * eps;
  }

  Eigen::Matrix3d toStressTensor( const double* v )
  {
    Eigen::Matrix3d s;
    s << v[0], v[3], v[4], v[3], v[1], v[5], v[4], v[5], v[2];
    return s;
  }

  std::unique_ptr< E > makeElement( const std::string&           mat,
                                    const std::vector< double >& props,
                                    double                       zeta,
                                    const double*                coords = flatCoordinates.data() )
  {
    auto el = std::make_unique< E >( 11,
                                     FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                     E::SectionType::Interface );
    el->assignNodeCoordinates( coords );
    static std::array< double, 2 > elp = { 1.0, 0.0 };
    elp[1]                             = zeta;
    ElementProperties ep( elp.data(), 2 );
    el->assignProperty( ep );
    el->assignMaterial( mat, props.data(), static_cast< int >( props.size() ) );
    return el;
  }

  std::unique_ptr< Y > makeYElement( const std::string&           mat,
                                     const std::vector< double >& props,
                                     const double*                coords = flatCoordinates.data() )
  {
    auto el = std::make_unique< Y >( 11,
                                     FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                     Y::SectionType::Interface );
    el->assignNodeCoordinates( coords );
    static std::array< double, 1 > elp = { 1.0 };
    ElementProperties              ep( elp.data(), 1 );
    el->assignProperty( ep );
    el->assignMaterial( mat, props.data(), static_cast< int >( props.size() ) );
    return el;
  }

  template < class El >
  void initState( El& el, std::vector< double >& sv )
  {
    sv.assign( el.getNumberOfRequiredStateVars(), 0.0 );
    el.assignStateVars( sv.data(), static_cast< int >( sv.size() ) );
    el.initializeYourself();
    el.setInitialConditions( MarmotElement::MarmotMaterialInitialization, nullptr );
  }

  struct Eval {
    Eigen::Matrix< double, nDofE, 1 >     Pe;
    Eigen::Matrix< double, nDofE, nDofE > Ke; // read column-major, as the element writes it
  };

  /** Evaluate on a scratch copy of stateIn; stateIn is never mutated. */
  Eval evaluate( E&                                       el,
                 const std::vector< double >&             stateIn,
                 const Eigen::Matrix< double, nDofE, 1 >& QTotal,
                 const Eigen::Matrix< double, nDofE, 1 >& dQ,
                 std::vector< double >*                   stateOut = nullptr )
  {
    std::vector< double > sv = stateIn;
    el.assignStateVars( sv.data(), static_cast< int >( sv.size() ) );

    std::array< double, nDofE >         Qa{}, dQa{}, Pe{};
    std::array< double, nDofE * nDofE > Ke{};
    std::copy( QTotal.data(), QTotal.data() + nDofE, Qa.begin() );
    std::copy( dQ.data(), dQ.data() + nDofE, dQa.begin() );
    el.computeKernels( Qa.data(), dQa.data(), Pe.data(), Ke.data(), 0.0, 1.0 );

    Eval e;
    e.Pe = Eigen::Map< Eigen::Matrix< double, nDofE, 1 > >( Pe.data() );
    e.Ke = Eigen::Map< Eigen::Matrix< double, nDofE, nDofE > >( Ke.data() ); // column-major
    if ( stateOut )
      *stateOut = sv;
    return e;
  }

  /** A single affine field u = c + F x applied to BOTH faces: the gradient
   * jump vanishes identically, so g = 0 is the exact answer. */
  Eigen::Matrix< double, nDofE, 1 > affineState( const Eigen::Vector3d& c, const Eigen::Matrix3d& F )
  {
    Eigen::Matrix< double, nDofE, 1 > q;
    q.setZero();
    for ( int A = 0; A < 4; A++ ) {
      const Eigen::Vector3d xB( flatXY[A][0], flatXY[A][1], 0.0 );
      const Eigen::Vector3d xT( flatXY[A][0], flatXY[A][1], hMesh );
      q.segment< 3 >( E::offD + 3 * A )                = c + F * xB; // u- bottom
      q.segment< 3 >( E::offD + E::nSideDofU + 3 * A ) = c + F * xT; // u+ top
      // g nodal stays zero
    }
    return q;
  }

  void testLayout()
  {
    throwExceptionOnFailure( E::sizeLoadVector == 36, "ZIQUAD4 must have 36 DOF." );
    throwExceptionOnFailure( E::offD == 0 && E::offG == 24, "offsets d/g = 0/24." );

    auto el = makeElement( "LINEARELASTIC", elasticProps(), 0.0 );
    throwExceptionOnFailure( el->getNDofPerElement() == 36, "getNDofPerElement 36." );

    const auto nf = el->getNodeFields();
    throwExceptionOnFailure( nf.size() == 8, "8 nodes." );
    for ( int i = 0; i < 4; i++ )
      throwExceptionOnFailure( nf[i].size() == 2 && nf[i][0] == "displacement" && nf[i][1] == "normalGradientJump",
                               "bottom nodes carry [displacement, normalGradientJump]." );
    for ( int i = 4; i < 8; i++ )
      throwExceptionOnFailure( nf[i].size() == 1 && nf[i][0] == "displacement", "top nodes carry [displacement]." );

    const auto         perm = el->getDofIndicesPermutationPattern();
    std::vector< int > s    = perm;
    std::sort( s.begin(), s.end() );
    throwExceptionOnFailure( perm.size() == 36, "perm size 36." );
    for ( int i = 0; i < 36; i++ )
      throwExceptionOnFailure( s[i] == i, "perm is a bijection onto [0,36)." );
  }

  /** A compatible affine field must produce zero traction imbalance, hence a
   * vanishing g-residual, for g = 0. This is eq. (23)'s dropped term being
   * genuinely zero where the theory says it is. */
  void testCompatibleAffineStateHasNoGradientJump()
  {
    auto                  el = makeElement( "LINEARELASTIC", elasticProps(), 0.0 );
    std::vector< double > sv;
    initState( *el, sv );

    Eigen::Matrix3d F;
    F << 1.0e-4, 2.0e-4, -3.0e-4, 5.0e-5, -1.5e-4, 2.0e-4, 3.0e-4, -2.0e-4, 4.0e-4;
    const auto q = affineState( Eigen::Vector3d( 1.0e-3, -2.0e-3, 5.0e-4 ), F );

    const Eval   e  = evaluate( *el, sv, q, q );
    const double Rg = e.Pe.segment< E::nDofG >( E::offG ).norm();
    const double Rd = e.Pe.segment< E::nDofD >( E::offD ).norm();

    std::cout << "  affine state: |R_g| = " << Rg << ", |R_d| = " << Rd << std::endl;
    throwExceptionOnFailure( Rg < 1.0e-10 * std::max( 1.0, Rd ),
                             "a compatible affine field must leave the gradient-jump residual at zero." );
  }

  /**
   * The residual, term by term against eq. (23).
   *
   * Rebuilt independently from the two faces' COMMITTED Cauchy stresses (read
   * back through getStateView) and the element's own raw geometric operators,
   * so it shares nothing with the Bx/Asrf machinery inside the material:
   *
   *   R_d = int_A ( t_i [uhat_i] + h <sigma>:Bhat^S_avg + (h/4)[sigma]:Bhat^S_jump ) dA
   *   R_g = int_A (h/4) [t_i] ghat_i dA                    (the term eq. 23 drops)
   *
   * with t_i = <sigma_ij> n_j of eq. (24). The mesh is collapsed (ell = h,
   * d_tau = 0), which is the configuration eq. (23) is written for.
   */
  void testResidualMatchesEquation23()
  {
    const auto            props = vonMisesProps();
    auto                  el    = makeElement( "VONMISES", props, 0.0 );
    std::vector< double > sv;
    initState( *el, sv );

    // a state with a genuine gradient jump AND a nonzero nodal g
    Eigen::Matrix< double, nDofE, 1 > q;
    q.setZero();
    Eigen::Matrix3d Fb, Ft;
    Fb << 2.0e-4, 0.0, 0.0, 0.0, -1.0e-4, 0.0, 3.0e-4, 1.0e-4, 0.0;
    Ft << -1.0e-4, 2.0e-4, 0.0, 1.5e-4, 3.0e-4, 0.0, -2.0e-4, 1.0e-4, 0.0;
    for ( int A = 0; A < 4; A++ ) {
      const Eigen::Vector3d x( flatXY[A][0], flatXY[A][1], 0.0 );
      q.segment< 3 >( E::offD + 3 * A )                = Fb * x;
      q.segment< 3 >( E::offD + E::nSideDofU + 3 * A ) = Ft * x + Eigen::Vector3d( 4.0e-4, -2.0e-4, 3.0e-4 );
      q.segment< 3 >( E::offG + 3 * A )                = Eigen::Vector3d( 1.0e-3, -5.0e-4, 2.0e-4 );
    }

    std::vector< double > svOut;
    const Eval            e = evaluate( *el, sv, q, q, &svOut );
    el->assignStateVars( svOut.data(), static_cast< int >( svOut.size() ) );

    const auto toStress = []( const double* v ) {
      Eigen::Matrix3d s;
      s << v[0], v[3], v[4], v[3], v[1], v[5], v[4], v[5], v[2];
      return s;
    };

    Eigen::Matrix< double, E::nDofD, 1 > Rd = Eigen::Matrix< double, E::nDofD, 1 >::Zero();
    Eigen::Matrix< double, E::nDofG, 1 > Rg = Eigen::Matrix< double, E::nDofG, 1 >::Zero();

    for ( int k = 0; k < el->getNumberOfQuadraturePoints(); ++k ) {
      const auto& qp = el->qps[k];

      const Eigen::Matrix3d sPlus  = toStress( el->getStateView( "stressPlus", k ).stateLocation );
      const Eigen::Matrix3d sMinus = toStress( el->getStateView( "stressMinus", k ).stateLocation );

      const Eigen::Matrix3d sAvg  = 0.5 * ( sPlus + sMinus ); // <sigma>
      const Eigen::Matrix3d sJump = sPlus - sMinus;           // [sigma]

      const Eigen::Vector3d t     = sAvg * qp.normal;         // eq. (24)
      const Eigen::Vector3d tJump = sJump * qp.normal;        // [t]

      // vec() in the element's row-major 3x3 order
      Eigen::Matrix< double, 9, 1 > vecAvg, vecJump;
      for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j ) {
          vecAvg( 3 * i + j )  = sAvg( i, j );
          vecJump( 3 * i + j ) = sJump( i, j );
        }

      // B for [uhat]^S_{,j} = A+ - A- : bottom nodes negative, top nodes positive
      Eigen::Matrix< double, 9, E::nDofD > BjumpSurface        = Eigen::Matrix< double, 9, E::nDofD >::Zero();
      BjumpSurface.block< 9, E::nSideDofU >( 0, 0 )            = -qp.BmatSide;
      BjumpSurface.block< 9, E::nSideDofU >( 0, E::nSideDofU ) = qp.BmatSide;

      Rd += ( qp.NmatJump.transpose() * t                               // t_i [uhat_i]
              + hMesh * ( qp.BmatAverage.transpose() * vecAvg )         // h sigma^S : uhat^S
              + 0.25 * hMesh * ( BjumpSurface.transpose() * vecJump ) ) // (h/4)[sigma]:[uhat]^S
            * qp.J0xW;

      Rg += qp.Ng.transpose() * ( 0.25 * hMesh * tJump ) * qp.J0xW; // (h/4)[t_i] ghat_i
    }

    // Pe holds MINUS the residual
    const Eigen::Matrix< double, E::nDofD, 1 > RdElement = -e.Pe.segment< E::nDofD >( E::offD );
    const Eigen::Matrix< double, E::nDofG, 1 > RgElement = -e.Pe.segment< E::nDofG >( E::offG );

    const double errD = ( RdElement - Rd ).norm() / std::max( 1.0, Rd.norm() );
    const double errG = ( RgElement - Rg ).norm() / std::max( 1.0, Rg.norm() );

    std::cout << "  R_d vs eq (23) terms : rel err = " << errD << "  (|R_d| = " << Rd.norm() << ")" << std::endl;
    std::cout << "  R_g vs (h/4)[t] ghat : rel err = " << errG << "  (|R_g| = " << Rg.norm() << ")" << std::endl;

    throwExceptionOnFailure( Rg.norm() > 1.0e-6 * Rd.norm(),
                             "this test state must have a genuinely nonzero traction imbalance." );
    throwExceptionOnFailure( errD < 1.0e-12, "R_d must be exactly the three terms of eq. (23)." );
    throwExceptionOnFailure( errG < 1.0e-12, "R_g must be exactly the (h/4)[t_i] term eq. (23) drops." );
  }

  /**
   * CLOSED FORM 1 -- the face reconstruction of eq. (44) is exact.
   *
   * For u = F x sampled on the OBLIQUE node pairs, the kinematics predict
   *
   *   Abar = A+ = A- = F T ,      w = F d ,
   *   gbar = ( w - Abar d_tau ) / ell = F ( d - d_tau ) / ell = F n ,
   *   G^pm = F T + ( F n +- g0/2 ) (x) n = F +- (g0/2) (x) n ,
   *
   * independently of ell and of d_tau -- both cancel. So each face's stress
   * must be the bulk response to sym( F +- (g0/2) (x) n ), a quantity computed
   * here from E and nu alone. A wrong SIGN on the d_tau correction, or ell
   * confused with h, breaks this; comparing against YIQUAD4 could not see it,
   * since both elements share the operator.
   */
  void testObliqueFaceReconstructionIsExact( const double* coords, const char* label )
  {
    const auto            props = elasticProps();
    auto                  el    = makeElement( "LINEARELASTIC", props, 0.0, coords );
    std::vector< double > sv;
    initState( *el, sv );

    Eigen::Matrix3d F;
    F << 3.0e-4, -1.0e-4, 2.0e-4, 1.5e-4, 2.5e-4, -3.0e-4, -2.0e-4, 4.0e-4, 1.0e-4;
    const Eigen::Vector3d g0( 1.2e-3, -8.0e-4, 6.0e-4 );

    const auto            q = affineStateOn( coords, F, g0 );
    std::vector< double > svOut;
    evaluate( *el, sv, q, q, &svOut );
    el->assignStateVars( svOut.data(), static_cast< int >( svOut.size() ) );

    double maxErr = 0.0, ellMin = 1e30, ellMax = -1e30, dTauMin = 1e30, dTauMax = -1e30;
    for ( int k = 0; k < el->getNumberOfQuadraturePoints(); ++k ) {
      // the normal, and with it ell and d_tau, varies per station on a warped face
      const Eigen::Vector3d n             = el->qps[k].normal;
      const Eigen::Matrix3d expectedPlus  = isotropicStress( Eigen::Matrix3d( F + 0.5 * g0 * n.transpose() ),
                                                            props[0],
                                                            props[1] );
      const Eigen::Matrix3d expectedMinus = isotropicStress( Eigen::Matrix3d( F - 0.5 * g0 * n.transpose() ),
                                                             props[0],
                                                             props[1] );

      const Eigen::Matrix3d sPlus  = toStressTensor( el->getStateView( "stressPlus", k ).stateLocation );
      const Eigen::Matrix3d sMinus = toStressTensor( el->getStateView( "stressMinus", k ).stateLocation );
      maxErr                       = std::max( maxErr, ( sPlus - expectedPlus ).norm() / expectedPlus.norm() );
      maxErr                       = std::max( maxErr, ( sMinus - expectedMinus ).norm() / expectedMinus.norm() );

      ellMin  = std::min( ellMin, el->qps[k].normalSeparation );
      ellMax  = std::max( ellMax, el->qps[k].normalSeparation );
      dTauMin = std::min( dTauMin, el->qps[k].tangentialSeparation.norm() );
      dTauMax = std::max( dTauMax, el->qps[k].tangentialSeparation.norm() );
    }

    std::cout << "  [" << label << "] ell in [" << ellMin << ", " << ellMax << "] (h = " << hMesh << "), |d_tau| in ["
              << dTauMin << ", " << dTauMax << "]" << std::endl;
    std::cout << "  [" << label << "] face stress vs closed form sym(F +- g0/2 (x) n): rel err = " << maxErr
              << std::endl;

    throwExceptionOnFailure( std::abs( ellMin - hMesh ) > 1.0e-3, "this mesh must have ell != h." );
    throwExceptionOnFailure( dTauMin > 1.0e-3, "this mesh must have a nonzero tangential separation." );
    throwExceptionOnFailure( maxErr < 1.0e-11, "eq. (44)'s face reconstruction must be exact for an affine field." );
  }

  /**
   * CLOSED FORM 2 -- the gradient-jump residual on the oblique mesh.
   *
   * For linear elasticity from a virgin state, [t]_i = C_ijkl n_j n_l g0_k
   * = Q^e_ik g0_k EXACTLY, so
   *
   *   R_g|_A = (h/4) ( int_A N_A dA ) Q^e g0 ,
   *
   * with no dependence on ell or d_tau at all. The nodal integrals are taken
   * from the element's own quadrature, everything else from E, nu and h.
   */
  void testObliqueGradientJumpResidualClosedForm()
  {
    const auto            props = elasticProps();
    auto                  el    = makeElement( "LINEARELASTIC", props, 0.0, obliqueCoordinates.data() );
    std::vector< double > sv;
    initState( *el, sv );

    Eigen::Matrix3d F;
    F << 3.0e-4, -1.0e-4, 2.0e-4, 1.5e-4, 2.5e-4, -3.0e-4, -2.0e-4, 4.0e-4, 1.0e-4;
    const Eigen::Vector3d g0( 1.2e-3, -8.0e-4, 6.0e-4 );

    const auto q = affineStateOn( obliqueCoordinates.data(), F, g0 );
    const Eval e = evaluate( *el, sv, q, q );

    const double          E = props[0], nu = props[1];
    const double          mu     = E / ( 2.0 * ( 1.0 + nu ) );
    const double          lambda = E * nu / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
    const Eigen::Vector3d n      = el->qps[0].normal;
    const Eigen::Matrix3d Qe     = mu * Eigen::Matrix3d::Identity() + ( lambda + mu ) * ( n * n.transpose() );

    Eigen::Matrix< double, E::nDofG, 1 > expected = Eigen::Matrix< double, E::nDofG, 1 >::Zero();
    for ( int k = 0; k < el->getNumberOfQuadraturePoints(); ++k ) {
      const auto& qp = el->qps[k];
      for ( int A = 0; A < 4; ++A )
        expected.segment< 3 >( 3 * A ) += 0.25 * hMesh * qp.N( A ) * ( Qe * g0 ) * qp.J0xW;
    }

    const Eigen::Matrix< double, E::nDofG, 1 > actual = -e.Pe.segment< E::nDofG >( E::offG );
    const double                               err    = ( actual - expected ).norm() / expected.norm();
    std::cout << "  R_g vs closed form (h/4)(int N) Q^e g0:        rel err = " << err << std::endl;

    throwExceptionOnFailure( err < 1.0e-12, "R_g must equal (h/4) int N^T Q^e g0 on the oblique mesh." );
  }

  /**
   * CLOSED FORM 3 -- the generalized eq. (23), derived from the energy alone.
   *
   * Varying phi = (h/2)( psi^+ + psi^- ) with G^pm = A^pm + (gbar +- z/2)(x)n
   * and gbar = ( w - Abar d_tau )/ell gives, using only the average identity
   * eq. (9) and t = <sigma> n:
   *
   *   delta phi = (h/ell) t_i dw_i
   *             + h <sigma>_ij dAbar_ij
   *             - (h/ell) t_i d^T_m dAbar_im            <- the d_tau correction
   *             + (h/4) [sigma]_ij dDeltaA_ij
   *             + (h/4) [t]_i dz_i .
   *
   * It reduces to eq. (23) for ell = h, d_tau = 0. The check below assembles
   * exactly these five terms from the element's raw geometric operators and
   * its committed face stresses -- nothing from Bx/Bz/Asrf, and nothing from
   * YIQUAD4.
   */
  void testGeneralizedResidualOnObliqueMesh()
  {
    const auto            props = vonMisesProps();
    auto                  el    = makeElement( "VONMISES", props, 0.0, obliqueCoordinates.data() );
    std::vector< double > sv;
    initState( *el, sv );

    Eigen::Matrix3d Fb, Ft;
    Fb << 2.0e-4, 0.0, 0.0, 0.0, -1.0e-4, 0.0, 3.0e-4, 1.0e-4, 0.0;
    Ft << -1.0e-4, 2.0e-4, 0.0, 1.5e-4, 3.0e-4, 0.0, -2.0e-4, 1.0e-4, 0.0;

    // deliberately NOT an affine field: the two faces get different gradients,
    // so DeltaA != 0 and every one of the five terms is populated
    Eigen::Matrix< double, nDofE, 1 > q;
    q.setZero();
    for ( int A = 0; A < 4; A++ ) {
      const Eigen::Vector3d xB( obliqueCoordinates[3 * A], obliqueCoordinates[3 * A + 1], 0.0 );
      const Eigen::Vector3d xT( obliqueCoordinates[12 + 3 * A], obliqueCoordinates[12 + 3 * A + 1], 0.0 );
      q.segment< 3 >( E::offD + 3 * A )                = Fb * xB;
      q.segment< 3 >( E::offD + E::nSideDofU + 3 * A ) = Ft * xT + Eigen::Vector3d( 4.0e-4, -2.0e-4, 3.0e-4 );
      q.segment< 3 >( E::offG + 3 * A )                = Eigen::Vector3d( 9.0e-4, -4.0e-4, 3.0e-4 );
    }

    std::vector< double > svOut;
    const Eval            e = evaluate( *el, sv, q, q, &svOut );
    el->assignStateVars( svOut.data(), static_cast< int >( svOut.size() ) );

    Eigen::Matrix< double, E::nDofD, 1 > Rd               = Eigen::Matrix< double, E::nDofD, 1 >::Zero();
    Eigen::Matrix< double, E::nDofG, 1 > Rg               = Eigen::Matrix< double, E::nDofG, 1 >::Zero();
    Eigen::Matrix< double, E::nDofD, 1 > dTauTermTotal    = Eigen::Matrix< double, E::nDofD, 1 >::Zero();
    double                               dTauContribution = 0.0;

    for ( int k = 0; k < el->getNumberOfQuadraturePoints(); ++k ) {
      const auto& qp = el->qps[k];

      const Eigen::Matrix3d sPlus  = toStressTensor( el->getStateView( "stressPlus", k ).stateLocation );
      const Eigen::Matrix3d sMinus = toStressTensor( el->getStateView( "stressMinus", k ).stateLocation );
      const Eigen::Matrix3d sAvg   = 0.5 * ( sPlus + sMinus );
      const Eigen::Matrix3d sJump  = sPlus - sMinus;

      const Eigen::Vector3d t     = sAvg * qp.normal;
      const Eigen::Vector3d tJump = sJump * qp.normal;

      const double          ell  = qp.normalSeparation;
      const Eigen::Vector3d dTau = qp.tangentialSeparation;

      // vec() in the element's row-major 3x3 order
      Eigen::Matrix< double, 9, 1 > vecAvg, vecJump, vecTdTau;
      for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j ) {
          vecAvg( 3 * i + j )   = sAvg( i, j );
          vecJump( 3 * i + j )  = sJump( i, j );
          vecTdTau( 3 * i + j ) = t( i ) * dTau( j ); // t (x) d^T
        }

      Eigen::Matrix< double, 9, E::nDofD > BjumpSurface        = Eigen::Matrix< double, 9, E::nDofD >::Zero();
      BjumpSurface.block< 9, E::nSideDofU >( 0, 0 )            = -qp.BmatSide;
      BjumpSurface.block< 9, E::nSideDofU >( 0, E::nSideDofU ) = qp.BmatSide;

      const auto dTauTerm = Eigen::Matrix< double, E::nDofD, 1 >( -( hMesh / ell ) *
                                                                  ( qp.BmatAverage.transpose() * vecTdTau ) * qp.J0xW );
      dTauTermTotal += dTauTerm;
      dTauContribution += dTauTerm.norm();

      Rd += ( hMesh / ell ) * ( qp.NmatJump.transpose() * t ) * qp.J0xW        // (h/ell) t_i dw_i
            + hMesh * ( qp.BmatAverage.transpose() * vecAvg ) * qp.J0xW        // h <sigma>:dAbar
            + dTauTerm                                                         // -(h/ell) t (x) d^T : dAbar
            + 0.25 * hMesh * ( BjumpSurface.transpose() * vecJump ) * qp.J0xW; // (h/4)[sigma]:dDeltaA

      Rg += qp.Ng.transpose() * ( 0.25 * hMesh * tJump ) * qp.J0xW;            // (h/4)[t]_i dz_i
    }

    const Eigen::Matrix< double, E::nDofD, 1 > RdElement = -e.Pe.segment< E::nDofD >( E::offD );
    const Eigen::Matrix< double, E::nDofG, 1 > RgElement = -e.Pe.segment< E::nDofG >( E::offG );

    const double errD = ( RdElement - Rd ).norm() / Rd.norm();
    const double errG = ( RgElement - Rg ).norm() / Rg.norm();

    // Negative controls: the same comparison with the d_tau correction DROPPED,
    // and with its SIGN FLIPPED. Both must fail loudly, otherwise the check
    // above has no power to detect an error in exactly the term it exists for.
    const double errNoDTau      = ( RdElement - ( Rd - dTauTermTotal ) ).norm() / Rd.norm();
    const double errFlippedDTau = ( RdElement - ( Rd - 2.0 * dTauTermTotal ) ).norm() / Rd.norm();

    std::cout << "  oblique R_d vs generalized eq (23): rel err = " << errD << std::endl;
    std::cout << "  oblique R_g vs (h/4)[t] ghat      : rel err = " << errG << std::endl;
    std::cout << "  the d_tau term contributes " << 100.0 * dTauContribution / Rd.norm() << "% of |R_d|" << std::endl;
    std::cout << "  negative control, d_tau dropped   : rel err = " << errNoDTau << std::endl;
    std::cout << "  negative control, d_tau sign flip : rel err = " << errFlippedDTau << std::endl;

    throwExceptionOnFailure( dTauContribution > 1.0e-3 * Rd.norm(),
                             "the d_tau term must be a non-negligible part of this test's residual." );
    throwExceptionOnFailure( errNoDTau > 1.0e-3 && errFlippedDTau > 1.0e-3,
                             "the check must be sensitive to the d_tau correction it is guarding." );
    throwExceptionOnFailure( errD < 1.0e-12, "R_d must equal the generalized eq. (23) on an oblique mesh." );
    throwExceptionOnFailure( errG < 1.0e-12, "R_g must equal (h/4)[t] on an oblique mesh." );
  }

  void testTangentAgainstFiniteDifferences( const std::string& mat, const std::vector< double >& props )
  {
    auto                  el = makeElement( mat, props, 0.0 );
    std::vector< double > sv;
    initState( *el, sv );

    Eigen::Matrix< double, nDofE, 1 > dQ;
    for ( int i = 0; i < nDofE; ++i )
      dQ( i ) = 1.0e-5 * std::sin( 0.6 + 0.83 * i );

    const Eval base = evaluate( *el, sv, dQ, dQ );

    const double eps    = 1.0e-9;
    double       maxErr = 0.0;
    for ( int c = 0; c < nDofE; ++c ) {
      Eigen::Matrix< double, nDofE, 1 > dQp = dQ;
      dQp( c ) += eps;
      const Eval p = evaluate( *el, sv, dQp, dQp );
      for ( int r = 0; r < nDofE; ++r ) {
        // Pe = -R, Ke = +dR/dq, so dPe/dq = -Ke
        maxErr = std::max( maxErr, std::abs( ( p.Pe( r ) - base.Pe( r ) ) / eps + base.Ke( r, c ) ) );
      }
    }

    const double scale = base.Ke.cwiseAbs().maxCoeff();
    std::cout << "  [" << mat << "] element FD tangent error (rel) = " << maxErr / scale << std::endl;
    throwExceptionOnFailure( base.Ke.allFinite(), "tangent finite." );
    throwExceptionOnFailure( maxErr < 1.0e-5 * scale, "element tangent inconsistent with finite differences." );
  }

  /**
   * With zeta = 0, solving R_g = 0 for the nodal g must reproduce YIQUAD4
   * exactly, both in the residual and -- after condensing g out of the 36x36
   * by hand -- in the 24x24 tangent.
   */
  void testEquivalenceWithYIQUAD4( const std::string&           mat,
                                   const std::vector< double >& props,
                                   const double*                coords = flatCoordinates.data(),
                                   const char*                  label  = "collapsed" )
  {
    auto                  el = makeElement( mat, props, 0.0, coords );
    std::vector< double > svZ;
    initState( *el, svZ );

    auto                  yel = makeYElement( mat, props, coords );
    std::vector< double > svY;
    initState( *yel, svY );

    // a displacement increment with a genuine gradient jump: the two faces
    // are given different in-plane gradients
    Eigen::Matrix< double, nDofE, 1 > q;
    q.setZero();
    Eigen::Matrix3d Fb, Ft;
    Fb << 1.0e-4, 0.0, 0.0, 0.0, -0.5e-4, 0.0, 2.0e-4, 1.0e-4, 0.0;
    Ft << -0.5e-4, 1.5e-4, 0.0, 1.0e-4, 2.0e-4, 0.0, -1.0e-4, 0.5e-4, 0.0;
    for ( int A = 0; A < 4; A++ ) {
      const Eigen::Vector3d x( flatXY[A][0], flatXY[A][1], 0.0 );
      q.segment< 3 >( E::offD + 3 * A )                = Fb * x;
      q.segment< 3 >( E::offD + E::nSideDofU + 3 * A ) = Ft * x + Eigen::Vector3d( 3.0e-4, -1.0e-4, 2.0e-4 );
    }

    // Newton on the g block only, holding the displacements fixed.
    Eval e;
    for ( int it = 0; it < 40; ++it ) {
      e                   = evaluate( *el, svZ, q, q );
      const auto residual = Eigen::Matrix< double, E::nDofG, 1 >( -e.Pe.segment< E::nDofG >( E::offG ) );
      if ( residual.norm() < 1.0e-14 )
        break;
      const auto Kgg = Eigen::Matrix< double, E::nDofG, E::nDofG >(
        e.Ke.block< E::nDofG, E::nDofG >( E::offG, E::offG ) );
      q.segment< E::nDofG >( E::offG ) -= Kgg.fullPivLu().solve( residual );
    }
    e                      = evaluate( *el, svZ, q, q );
    const double gResidual = e.Pe.segment< E::nDofG >( E::offG ).norm();
    throwExceptionOnFailure( gResidual < 1.0e-10, "the nodal-g Newton did not converge." );

    // YIQUAD4 on the same displacement increment
    std::array< double, nDofY >         Qy{}, dQy{}, Py{};
    std::array< double, nDofY * nDofY > Ky{};
    for ( int i = 0; i < nDofY; ++i ) {
      Qy[i]  = q( E::offD + i );
      dQy[i] = q( E::offD + i );
    }
    yel->computeKernels( Qy.data(), dQy.data(), Py.data(), Ky.data(), 0.0, 1.0 );
    const auto PeY = Eigen::Map< Eigen::Matrix< double, nDofY, 1 > >( Py.data() );
    const auto KeY = Eigen::Map< Eigen::Matrix< double, nDofY, nDofY > >( Ky.data() ); // column-major

    const auto PeZ = Eigen::Matrix< double, nDofY, 1 >( e.Pe.segment< nDofY >( E::offD ) );

    const auto Kdd = Eigen::Matrix< double, nDofY, nDofY >( e.Ke.block< nDofY, nDofY >( E::offD, E::offD ) );
    const auto Kdg = Eigen::Matrix< double, nDofY, E::nDofG >( e.Ke.block< nDofY, E::nDofG >( E::offD, E::offG ) );
    const auto Kgd = Eigen::Matrix< double, E::nDofG, nDofY >( e.Ke.block< E::nDofG, nDofY >( E::offG, E::offD ) );
    const auto Kgg = Eigen::Matrix< double, E::nDofG, E::nDofG >(
      e.Ke.block< E::nDofG, E::nDofG >( E::offG, E::offG ) );

    const Eigen::Matrix< double, nDofY, nDofY > condensed = Kdd - Kdg * Kgg.fullPivLu().solve( Kgd );

    const double residualError = ( PeZ - PeY ).norm() / std::max( 1.0, PeY.norm() );
    const double tangentError  = ( condensed - KeY ).norm() / KeY.norm();

    std::cout << "  [" << mat << ", " << label << "] vs YIQUAD4: |dPe| = " << residualError
              << ", |dKe| = " << tangentError << std::endl;

    throwExceptionOnFailure( residualError < 1.0e-9, "ZIQUAD4 residual must match YIQUAD4 once R_g = 0." );
    throwExceptionOnFailure( tangentError < 1.0e-7,
                             "the condensed ZIQUAD4 tangent must match YIQUAD4's condensed tangent." );
  }

  /**
   * The point of the element. Drive both faces past yield in the same shear
   * mode, so that the averaged acoustic tensor degenerates along m^T
   * (eq. 40/48). The g block must then be singular WITHOUT the regularization
   * -- which is exactly why YIQUAD4 needs the truncated inverse of eq. (51) --
   * and invertible WITH it.
   */
  /**
   * The term eq. (23) drops, measured where it actually matters.
   *
   * Under PERFECT plasticity (H = 0) the averaged acoustic tensor is fully
   * degenerate, cond(<Q>) ~ 1e-16, and the gauge of eqs. (49)-(51) is active.
   * The condensed kernel then commits a traction imbalance that lies ENTIRELY
   * along the undetermined direction m^T of eq. (49) -- no 3x3 solve can reach
   * a component of r that no g affects -- and its size grows QUADRATICALLY with
   * the mismatch alpha between the two faces, |[t].m^T|/|t| ~ 4.7e-3 alpha^2.
   *
   * That imbalance is the work conjugate to ghat. YIQUAD4 has no ghat degree of
   * freedom, so it does not "leave the remainder to the global Newton" as
   * section "condensation logic" states: it deletes it, and nothing minimises
   * it. ZIQUAD4 carries it as R_g and the global Newton drives it to zero.
   *
   * NOTE ON MEASUREMENT: judge degeneracy from the MATERIAL's
   * K_zz = (h/4)<Q>, never from the element block K_gg, which carries the Q1
   * mass matrix's own conditioning (~1/81) on top and reads a factor ~80 too
   * healthy against the gauge threshold tau = 1e-6.
   */
  std::vector< double > perfectlyPlasticProps()
  {
    return { 4.0e5, 0.3, hMesh, 5.0, 0.0, 0.0, 0.0, 0.0 };
  }

  /** Drive the condensed kernel to a perfectly-plastic state with a face
   * mismatch alpha, and return |[t].m^T| / |t| and |[t]| / |t|. */
  std::pair< double, double > measureCondensedImbalance( double alpha )
  {
    const auto                                      props = perfectlyPlasticProps();
    MarmotEquilibratedXInterfaceMaterialHypoElastic material( "VONMISES",
                                                              props.data(),
                                                              static_cast< int >( props.size() ),
                                                              0 );
    std::vector< double >                           sv( material.getNumberOfRequiredStateVars(), 0.0 );
    material.initializeYourself( sv.data(), static_cast< int >( sv.size() ) );

    const Eigen::Vector3d n( 0.0, 0.0, 1.0 );
    const Eigen::Vector3d sep( 0.0, 0.0, hMesh );

    Eigen::Matrix< double, 6, 1 >  dU = Eigen::Matrix< double, 6, 1 >::Zero();
    Eigen::Matrix< double, 18, 1 > dS = Eigen::Matrix< double, 18, 1 >::Zero();
    dU( 0 )                           = 3.0e-4; // shear jump through the layer
    Eigen::Matrix3d dFa, dMismatch;
    dFa << 0, 0, 0, 0, 0, 0, 1.0e-4, 0.5e-4, 0;
    dMismatch << 0, 1.0e-4, 0, 2.0e-4, 0, 0, -1.5e-4, 0.8e-4, 0;
    const Eigen::Matrix3d dFb = dFa + alpha * dMismatch;
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j ) {
        dS( 3 * i + j )     = dFa( i, j );
        dS( 9 + 3 * i + j ) = dFb( i, j );
      }

    for ( int step = 0; step < 60; ++step ) {
      Eigen::Vector3d                                                f  = Eigen::Vector3d::Zero();
      Eigen::Matrix< double, 9, 1 >                                  sp = Eigen::Matrix< double, 9, 1 >::Zero();
      Eigen::Matrix< double, 9, 1 >                                  sm = Eigen::Matrix< double, 9, 1 >::Zero();
      Eigen::Matrix< double, 3, 3, Eigen::RowMajor >                 Qww;
      Eigen::Matrix< double, 3, 9, Eigen::RowMajor >                 QwAp, QwAm;
      Eigen::Matrix< double, 9, 3, Eigen::RowMajor >                 QApw, QAmw;
      Eigen::Matrix< double, 9, 9, Eigen::RowMajor >                 QApAp, QApAm, QAmAp, QAmAm;
      MarmotEquilibratedXInterfaceMaterialHypoElastic::State         state{ f.data(), sp.data(), sm.data(), sv.data() };
      MarmotEquilibratedXInterfaceMaterialHypoElastic::Tangents      tangents{ Qww.data(),
                                                                          QwAp.data(),
                                                                          QwAm.data(),
                                                                          QApw.data(),
                                                                          QApAp.data(),
                                                                          QApAm.data(),
                                                                          QAmw.data(),
                                                                          QAmAp.data(),
                                                                          QAmAm.data() };
      MarmotEquilibratedXInterfaceMaterialHypoElastic::Deformation   def{ dU.data(), dS.data(), n.data(), sep.data() };
      MarmotEquilibratedXInterfaceMaterialHypoElastic::TimeIncrement ti{ double( step ), 1.0 };
      material.computeStress( state, tangents, def, ti );
    }

    const Eigen::Matrix3d sPlus  = toStressTensor( material.getStateView( "stressPlus", sv.data() ).stateLocation );
    const Eigen::Matrix3d sMinus = toStressTensor( material.getStateView( "stressMinus", sv.data() ).stateLocation );

    const auto unitDev = []( const Eigen::Matrix3d& s ) {
      const Eigen::Matrix3d d = s - ( s.trace() / 3.0 ) * Eigen::Matrix3d::Identity();
      return d.norm() > 0.0 ? Eigen::Matrix3d( d / d.norm() ) : Eigen::Matrix3d::Zero();
    };

    const Eigen::Vector3d tPlus  = sPlus * n;
    const Eigen::Vector3d tMinus = sMinus * n;
    const Eigen::Vector3d tJump  = tPlus - tMinus;
    const double          scale  = std::max( tPlus.norm(), tMinus.norm() );

    Eigen::Vector3d m = Eigen::Matrix3d( 0.5 * ( unitDev( sPlus ) + unitDev( sMinus ) ) ) * n;
    m -= m.dot( n ) * n; // eq. (49): the gauge direction is purely tangential
    const Eigen::Vector3d mT = m.norm() > 1.0e-14 ? Eigen::Vector3d( m.normalized() ) : Eigen::Vector3d::Zero();

    if ( scale <= 0.0 )
      return { 0.0, 0.0 };
    return { std::abs( tJump.dot( mT ) ) / scale, tJump.norm() / scale };
  }

  void testDroppedTermIsTheUndeterminedDirection()
  {
    std::cout << "  perfect plasticity (H = 0), condensed kernel:" << std::endl;
    std::cout << "      alpha      |[t]|/|t|     |[t].mT|/|t|   fraction along mT" << std::endl;

    double largest = 0.0;
    for ( double alpha : { 1.0e-3, 1.0e-2, 1.0e-1 } ) {
      const auto   r        = measureCondensedImbalance( alpha );
      const double fraction = r.second > 0.0 ? r.first / r.second : 0.0;
      largest               = std::max( largest, r.first );
      std::printf( "    %9.1e   %11.4e   %11.4e   %8.5f\n", alpha, r.second, r.first, fraction );

      // The imbalance the condensation cannot remove is the m^T component, and
      // ONLY the m^T component: the normal direction is always solved.
      throwExceptionOnFailure( fraction > 0.99,
                               "the committed traction imbalance must lie along the undetermined direction m^T." );
    }
    std::fflush( stdout );

    // ... and it is not roundoff: it grows quadratically with the face mismatch,
    // so two decades of alpha must move it by ~1e4.
    const double small = measureCondensedImbalance( 1.0e-3 ).first;
    const double big   = measureCondensedImbalance( 1.0e-1 ).first;
    std::cout << "  growth over two decades of mismatch: " << big / small << "x (quadratic => ~1e4)" << std::endl;

    throwExceptionOnFailure( largest > 1.0e-6,
                             "under perfect plasticity the condensed kernel must leave a traction imbalance "
                             "far above roundoff -- this is the (h/4)[t] ghat term eq. (23) drops." );
    throwExceptionOnFailure( big / small > 1.0e3,
                             "the dropped term must scale (quadratically) with the face mismatch, not sit at "
                             "roundoff." );
  }

  void testDegenerateStateStaysSolvable()
  {
    const auto props = vonMisesProps();

    auto                  el = makeElement( "VONMISES", props, 0.0 );
    std::vector< double > sv;
    initState( *el, sv );

    // load in tangential shear across the layer, in steps, updating the state
    Eigen::Matrix< double, nDofE, 1 > q;
    q.setZero();
    Eigen::Matrix< double, nDofE, 1 > dq;
    dq.setZero();
    for ( int A = 0; A < 4; A++ )
      dq( E::offD + E::nSideDofU + 3 * A + 0 ) = 4.0e-4; // uniform tangential jump w_x

    for ( int step = 0; step < 30; ++step ) {
      q += dq;
      std::vector< double > next;
      // one Newton pass on g per step keeps the layer near equilibrium
      for ( int it = 0; it < 6; ++it ) {
        const Eval s        = evaluate( *el, sv, q, dq );
        const auto residual = Eigen::Matrix< double, E::nDofG, 1 >( -s.Pe.segment< E::nDofG >( E::offG ) );
        const auto Kgg      = Eigen::Matrix< double, E::nDofG, E::nDofG >(
          s.Ke.block< E::nDofG, E::nDofG >( E::offG, E::offG ) );
        Eigen::FullPivLU< Eigen::Matrix< double, E::nDofG, E::nDofG > > lu( Kgg );
        if ( residual.norm() < 1.0e-12 || lu.rank() < E::nDofG )
          break;
        q.segment< E::nDofG >( E::offG ) -= lu.solve( residual );
      }
      evaluate( *el, sv, q, dq, &next );
      sv = next;
    }

    const Eval plain    = evaluate( *el, sv, q, dq );
    const auto KggPlain = Eigen::Matrix< double, E::nDofG, E::nDofG >(
      plain.Ke.block< E::nDofG, E::nDofG >( E::offG, E::offG ) );
    Eigen::JacobiSVD< Eigen::MatrixXd > svdPlain( ( Eigen::MatrixXd( KggPlain ) ) );
    const double condPlain = svdPlain.singularValues()( E::nDofG - 1 ) / svdPlain.singularValues()( 0 );

    auto                  regularized = makeElement( "VONMISES", props, 1.0e-6 );
    std::vector< double > svr         = sv;
    regularized->assignStateVars( svr.data(), static_cast< int >( svr.size() ) );
    regularized->initializeYourself();
    const Eval reg    = evaluate( *regularized, sv, q, dq );
    const auto KggReg = Eigen::Matrix< double, E::nDofG, E::nDofG >(
      reg.Ke.block< E::nDofG, E::nDofG >( E::offG, E::offG ) );
    Eigen::JacobiSVD< Eigen::MatrixXd > svdReg( ( Eigen::MatrixXd( KggReg ) ) );
    const double condReg = svdReg.singularValues()( E::nDofG - 1 ) / svdReg.singularValues()( 0 );

    std::cout << "  yielded layer: cond(K_gg) with zeta=0    : " << condPlain << std::endl;
    std::cout << "                 cond(K_gg) with zeta=1e-6 : " << condReg << std::endl;

    throwExceptionOnFailure( plain.Ke.allFinite() && reg.Ke.allFinite(),
                             "tangents must stay finite when <Q> degenerates." );
    throwExceptionOnFailure( condReg > condPlain, "the regularization must improve the conditioning of the g block." );
    throwExceptionOnFailure( condReg > 1.0e-9, "the regularized g block must be numerically invertible." );

    // How much of the cure comes from the coupling to the displacement DOFs
    // alone: count the near-null modes of the FULL 36x36 with zeta = 0.
    Eigen::JacobiSVD< Eigen::MatrixXd > svdPlainFull( ( Eigen::MatrixXd( plain.Ke ) ) );
    const auto                          svPlain36  = svdPlainFull.singularValues();
    int                                 nNearPlain = 0;
    for ( int i = 0; i < nDofE; ++i )
      if ( svPlain36( i ) < 1.0e-11 * svPlain36( 0 ) )
        ++nNearPlain;
    std::cout << "                 near-null modes of the 36x36, zeta=0: " << nNearPlain << std::endl;

    // and the whole 36x36 must have no null space beyond the 6 rigid-body modes
    Eigen::JacobiSVD< Eigen::MatrixXd > svdFull( ( Eigen::MatrixXd( reg.Ke ) ) );
    const auto                          sv36  = svdFull.singularValues();
    int                                 nNear = 0;
    for ( int i = 0; i < nDofE; ++i )
      if ( sv36( i ) < 1.0e-11 * sv36( 0 ) )
        ++nNear;
    std::cout << "                 near-null modes of the 36x36: " << nNear << " (6 rigid-body expected)" << std::endl;
    throwExceptionOnFailure( nNear <= 6, "no spurious null modes beyond rigid-body motion." );
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{
    [&]() { testLayout(); },
    [&]() { testCompatibleAffineStateHasNoGradientJump(); },
    [&]() { testResidualMatchesEquation23(); },
    [&]() { testTangentAgainstFiniteDifferences( "LINEARELASTIC", elasticProps() ); },
    [&]() { testTangentAgainstFiniteDifferences( "VONMISES", vonMisesProps() ); },
    [&]() { testEquivalenceWithYIQUAD4( "LINEARELASTIC", elasticProps() ); },
    [&]() { testEquivalenceWithYIQUAD4( "VONMISES", vonMisesProps() ); },
    // eq. (44)'s oblique-connector correction: ell != h and d_tau != 0
    [&]() { testEquivalenceWithYIQUAD4( "LINEARELASTIC", elasticProps(), obliqueCoordinates.data(), "oblique" ); },
    [&]() { testEquivalenceWithYIQUAD4( "VONMISES", vonMisesProps(), obliqueCoordinates.data(), "oblique" ); },
    [&]() { testObliqueFaceReconstructionIsExact( obliqueCoordinates.data(), "oblique" ); },
    [&]() { testObliqueFaceReconstructionIsExact( warpedCoordinates.data(), "warped+oblique" ); },
    [&]() { testObliqueGradientJumpResidualClosedForm(); },
    [&]() { testGeneralizedResidualOnObliqueMesh(); },
    [&]() { testDroppedTermIsTheUndeterminedDirection(); },
    [&]() { testDegenerateStateStaysSolvable(); },
  };

  executeTestsAndCollectExceptions( tests );
  return 0;
}
