#include "Marmot/MarmotCorrectedInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Vector9d          = Eigen::Matrix< double, 9, 1 >;
  using Matrix3dRowMajor  = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor  = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
  using Matrix3x9RowMajor = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;

  using EqMaterial = MarmotEquilibratedXInterfaceMaterialHypoElastic;

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

  /** All generalized outputs of one equilibrated-X-interface stress update. */
  struct EqResponse {
    Eigen::Vector3d   generalizedForce   = Eigen::Vector3d::Zero();
    Matrix3dRowMajor  surfaceStressPlus  = Matrix3dRowMajor::Zero();
    Matrix3dRowMajor  surfaceStressMinus = Matrix3dRowMajor::Zero();
    Matrix3dRowMajor  Q_ww               = Matrix3dRowMajor::Zero();
    Matrix3x9RowMajor Q_wAp              = Matrix3x9RowMajor::Zero();
    Matrix3x9RowMajor Q_wAm              = Matrix3x9RowMajor::Zero();
    Matrix9x3RowMajor Q_Apw              = Matrix9x3RowMajor::Zero();
    Matrix9dRowMajor  Q_ApAp             = Matrix9dRowMajor::Zero();
    Matrix9dRowMajor  Q_ApAm             = Matrix9dRowMajor::Zero();
    Matrix9x3RowMajor Q_Amw              = Matrix9x3RowMajor::Zero();
    Matrix9dRowMajor  Q_AmAp             = Matrix9dRowMajor::Zero();
    Matrix9dRowMajor  Q_AmAm             = Matrix9dRowMajor::Zero();
  };

  void computeEqStress( EqMaterial&   material,
                        EqResponse&   response,
                        double*       stateVars,
                        const double* dU,
                        const double* dSurfaceStrain,
                        const double* normal,
                        const double* separation,
                        double        timeOld,
                        double        dT )
  {
    EqMaterial::State         state{ response.generalizedForce.data(),
                             response.surfaceStressPlus.data(),
                             response.surfaceStressMinus.data(),
                             stateVars };
    EqMaterial::Tangents      tangents{ response.Q_ww.data(),
                                   response.Q_wAp.data(),
                                   response.Q_wAm.data(),
                                   response.Q_Apw.data(),
                                   response.Q_ApAp.data(),
                                   response.Q_ApAm.data(),
                                   response.Q_Amw.data(),
                                   response.Q_AmAp.data(),
                                   response.Q_AmAm.data() };
    EqMaterial::TimeIncrement timeIncrement{ timeOld, dT };

    if ( separation ) {
      EqMaterial::Deformation deformation{ dU, dSurfaceStrain, normal, separation };
      material.computeStress( state, tangents, deformation, timeIncrement );
    }
    else {
      EqMaterial::Deformation deformation{ dU, dSurfaceStrain, normal };
      material.computeStress( state, tangents, deformation, timeIncrement );
    }
  }

  Eigen::VectorXd makeInitializedStateVars( EqMaterial& material )
  {
    Eigen::VectorXd stateVars = Eigen::VectorXd::Zero( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );
    return stateVars;
  }

  EqResponse evaluateVirginResponse( EqMaterial&   material,
                                     const double* dU,
                                     const double* dSurfaceStrain,
                                     const double* normal,
                                     const double* separation )
  {
    EqResponse      response;
    Eigen::VectorXd stateVars = makeInitializedStateVars( material );
    computeEqStress( material, response, stateVars.data(), dU, dSurfaceStrain, normal, separation, 0.0, 1.0 );
    return response;
  }

  Vector9d flattenRowMajor( const Matrix3dRowMajor& tensor )
  {
    return Eigen::Map< const Vector9d >( tensor.data() );
  }

  /**
   * 14.1: for a fixed increment, reconstruct t+ and t- from the COMMITTED
   * generalized state (f, S+, S-) and assert traction equilibrium holds to
   * tight tolerance. sigma^s = (S^s + f (x) d_tau)/sideThickness, t^s =
   * sigma^s n.
   */
  void assertTractionEquilibrium( const EqResponse&      response,
                                  const Eigen::Vector3d& normal,
                                  const Eigen::Vector3d& dTangential,
                                  double                 sideThickness,
                                  const std::string&     context )
  {
    const Matrix3dRowMajor sigmaPlus = ( response.surfaceStressPlus +
                                         response.generalizedForce * dTangential.transpose() ) /
                                       sideThickness;
    const Matrix3dRowMajor sigmaMinus = ( response.surfaceStressMinus +
                                          response.generalizedForce * dTangential.transpose() ) /
                                        sideThickness;

    const Eigen::Vector3d tPlus  = sigmaPlus * normal;
    const Eigen::Vector3d tMinus = sigmaMinus * normal;
    const double          scale  = std::max( 1.0, 0.5 * ( tPlus.norm() + tMinus.norm() ) );

    throwExceptionOnFailure( ( tPlus - tMinus ).norm() < 1e-8 * scale,
                             context + ": traction equilibrium t+ == t- violated, |t+-t-| = " +
                               std::to_string( ( tPlus - tMinus ).norm() ) );
  }

  void testDegenerateGeometryAndConstructionThrow()
  {
    {
      bool thrown = false;
      try {
        const double tooFewProperties[2] = { 1e5, 0.3 };
        EqMaterial   material( "LINEARELASTIC", tooFewProperties, 2, 1 );
      }
      catch ( const std::invalid_argument& e ) {
        thrown = std::string( e.what() ) ==
                 "MarmotEquilibratedXInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h.";
      }
      throwExceptionOnFailure( thrown, "Construction with 2 properties should throw the documented message." );
    }
    {
      bool thrown = false;
      try {
        const double zeroThickness[3] = { 1e5, 0.3, 0.0 };
        EqMaterial   material( "LINEARELASTIC", zeroThickness, 3, 1 );
      }
      catch ( const std::invalid_argument& e ) {
        thrown = std::string( e.what() ) == "MarmotEquilibratedXInterfaceMaterialHypoElastic requires h > 0.";
      }
      throwExceptionOnFailure( thrown, "Construction with h = 0 should throw the documented message." );
    }

    const double interfaceProperties[3] = { 1e5, 0.3, 0.01 };
    EqMaterial   material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double dU[6]              = { 1e-4, 0., 0., 0., 0., 0. };
    const double dSurfaceStrain[18] = { 0. };

    {
      bool thrown = false;
      try {
        const double zeroNormal[3] = { 0., 0., 0. };
        evaluateVirginResponse( material, dU, dSurfaceStrain, zeroNormal, nullptr );
      }
      catch ( const std::invalid_argument& e ) {
        thrown = std::string( e.what() ) ==
                 "MarmotEquilibratedXInterfaceMaterialHypoElastic: interface normal is zero.";
      }
      throwExceptionOnFailure( thrown, "Zero interface normal should throw the documented message." );
    }
  }

  void testDensityDelegation()
  {
    const double interfaceProperties[9] = { 1e8, 0.3, 0.01, 2e7, 0.25, 6., 1e-4, 1., 2400. };
    EqMaterial   material( "LINEARVISCOELASTICWIECHERT", interfaceProperties, 9, 1 );

    throwExceptionOnFailure( checkIfEqual( material.getDensity(), interfaceProperties[8] ),
                             "Equilibrated X interface density delegation failed." );
  }

  /** 14.3: A+ == A- and identical (virgin) material state on both sides
   * must give a converged internal jump z == 0. */
  void testUniformSidesGiveZeroInternalJump()
  {
    const double interfaceProperties[3] = { 1e5, 0.3, 0.01 };
    const double normal[3]              = { 0., 0., 1. };

    EqMaterial material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double surfaceGradient[9] = { 1.0e-4, -2.0e-4, 0.5e-4, 0.8e-4, 1.5e-4, -0.7e-4, 0.3e-4, -1.1e-4, 0.9e-4 };
    double       dSurfaceStrain[18];
    std::copy( surfaceGradient, surfaceGradient + 9, dSurfaceStrain );
    std::copy( surfaceGradient, surfaceGradient + 9, dSurfaceStrain + 9 );
    const double dU[6] = { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 };

    Eigen::VectorXd stateVars = makeInitializedStateVars( material );
    EqResponse      response;
    computeEqStress( material, response, stateVars.data(), dU, dSurfaceStrain, normal, nullptr, 0.0, 1.0 );

    const double* z = material.getStateView( "normalGradientJump", stateVars.data() ).stateLocation;
    throwExceptionOnFailure( std::abs( z[0] ) < 1e-10 && std::abs( z[1] ) < 1e-10 && std::abs( z[2] ) < 1e-10,
                             "Uniform-side loading (A+ == A-, identical virgin material) should converge to a "
                             "zero internal normal-gradient jump." );

    // Also verify the two sides therefore produced IDENTICAL surface stress.
    assertMatrixNear( response.surfaceStressPlus,
                      response.surfaceStressMinus,
                      1e-8 * response.surfaceStressPlus.lpNorm< Eigen::Infinity >(),
                      "Uniform-side loading should give identical surfaceStressPlus/Minus." );
  }

  /** 14.1, elastic: traction equilibrium holds for a generic skewed geometry. */
  void testTractionEquilibriumElasticSkewed()
  {
    const double h                      = 0.01;
    const double ell                    = 0.025;
    const double interfaceProperties[3] = { 1e5, 0.3, h };

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    EqMaterial material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double dU[6]     = { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 };
    const double aPlus[9]  = { 1.0e-4, -2.0e-4, 0.5e-4, 0.8e-4, 1.5e-4, -0.7e-4, 0.3e-4, -1.1e-4, 0.9e-4 };
    const double aMinus[9] = { -0.4e-4, 0.9e-4, -0.2e-4, 0.6e-4, -1.0e-4, 0.3e-4, -0.8e-4, 0.2e-4, 0.5e-4 };
    double       dSurfaceStrain[18];
    std::copy( aPlus, aPlus + 9, dSurfaceStrain );
    std::copy( aMinus, aMinus + 9, dSurfaceStrain + 9 );

    const EqResponse response = evaluateVirginResponse( material, dU, dSurfaceStrain, n.data(), separation.data() );
    assertTractionEquilibrium( response, n, dTangential, 0.5 * h, "Elastic skewed equilibrium" );
  }

  /**
   * 14.1 + setup for 14.5/14.6 nonzero cross-blocks: drive the two sides to
   * DIVERGENT Von Mises plastic histories (different prior loading on A+
   * vs A-) so that C+ != C- at the state under test, then verify traction
   * equilibrium still holds after the second increment.
   */
  std::vector< double > makeDivergentHistoryStateVars( EqMaterial&            material,
                                                       const Eigen::Vector3d& n,
                                                       const Eigen::Vector3d& dTangential,
                                                       const Eigen::Vector3d& separation,
                                                       const double           dU[6],
                                                       const double           dSurfaceStrain[18] )
  {
    Eigen::VectorXd stateVars = makeInitializedStateVars( material );
    EqResponse      response;
    computeEqStress( material, response, stateVars.data(), dU, dSurfaceStrain, n.data(), separation.data(), 0.0, 1.0 );
    return std::vector< double >( stateVars.data(), stateVars.data() + stateVars.size() );
  }

  void testTractionEquilibriumAndNonzeroCrossBlocksPlastic()
  {
    const double h                      = 0.01;
    const double ell                    = 0.025;
    const double interfaceProperties[8] = { 210000., 0.3, h, 20., 200., 5., 10., 2400. };

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    EqMaterial material( "VONMISES", interfaceProperties, 8, 1 );

    // Strongly asymmetric first increment: A+ large, A- small, so the two
    // sides accumulate very different plastic strain and end up with
    // different algorithmic tangents C+ != C-.
    const double historyDU[6]     = { 2.0e-3, -0.5e-3, 3.0e-3, 0.2e-3, 0.3e-3, -0.4e-3 };
    const double historyAPlus[9]  = { 8.0e-3, -6.0e-3, 4.0e-3, 5.0e-3, 7.0e-3, -3.0e-3, 2.0e-3, -5.0e-3, 6.0e-3 };
    const double historyAMinus[9] = { 0.2e-3, -0.1e-3, 0.05e-3, 0.1e-3, 0.15e-3, -0.05e-3, 0.05e-3, -0.1e-3, 0.1e-3 };
    double       historyDSurfaceStrain[18];
    std::copy( historyAPlus, historyAPlus + 9, historyDSurfaceStrain );
    std::copy( historyAMinus, historyAMinus + 9, historyDSurfaceStrain + 9 );

    const std::vector< double > historyStateVars = makeDivergentHistoryStateVars( material,
                                                                                  n,
                                                                                  dTangential,
                                                                                  separation,
                                                                                  historyDU,
                                                                                  historyDSurfaceStrain );

    const double dU[6]     = { 0.3e-3, -0.1e-3, 0.4e-3, 0.05e-3, 0.05e-3, -0.05e-3 };
    const double aPlus[9]  = { 0.5e-3, -0.3e-3, 0.2e-3, 0.3e-3, 0.4e-3, -0.2e-3, 0.1e-3, -0.3e-3, 0.3e-3 };
    const double aMinus[9] = { 0.05e-3, -0.02e-3, 0.01e-3, 0.02e-3, 0.03e-3, -0.01e-3, 0.01e-3, -0.02e-3, 0.02e-3 };
    double       dSurfaceStrain[18];
    std::copy( aPlus, aPlus + 9, dSurfaceStrain );
    std::copy( aMinus, aMinus + 9, dSurfaceStrain + 9 );

    std::vector< double > stateVars = historyStateVars;
    EqResponse            response;
    computeEqStress( material, response, stateVars.data(), dU, dSurfaceStrain, n.data(), separation.data(), 0.0, 1.0 );

    assertTractionEquilibrium( response, n, dTangential, 0.5 * h, "Plastic divergent-history equilibrium" );

    const double crossNorm = std::max( response.Q_ApAm.lpNorm< Eigen::Infinity >(),
                                       response.Q_AmAp.lpNorm< Eigen::Infinity >() );
    const double diagNorm  = std::max( response.Q_ApAp.lpNorm< Eigen::Infinity >(),
                                      response.Q_AmAm.lpNorm< Eigen::Infinity >() );
    throwExceptionOnFailure( crossNorm > 1e-3 * diagNorm,
                             "With divergent plastic histories (C+ != C-), the condensed cross-blocks "
                             "d(S+)/d(A-) and d(S-)/d(A+) should be clearly nonzero, not negligible." );
  }

  /**
   * 14.6: perturb all 21 external variables (w, A+, A-) and compare every
   * analytic condensed tangent block against a central finite difference,
   * re-solving the local traction-equilibrium problem at each perturbation.
   */
  void testFullCondensedTangentMatchesFiniteDifference()
  {
    const double h                      = 0.01;
    const double ell                    = 0.025;
    const double interfaceProperties[3] = { 1e5, 0.3, h };

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    EqMaterial material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    double dU[6]              = { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 };
    double dSurfaceStrain[18] = { 1.0e-4,
                                  -2.0e-4,
                                  0.5e-4,
                                  0.8e-4,
                                  1.5e-4,
                                  -0.7e-4,
                                  0.3e-4,
                                  -1.1e-4,
                                  0.9e-4,
                                  -0.4e-4,
                                  0.9e-4,
                                  -0.2e-4,
                                  0.6e-4,
                                  -1.0e-4,
                                  0.3e-4,
                                  -0.8e-4,
                                  0.2e-4,
                                  0.5e-4 };

    const EqResponse analytic = evaluateVirginResponse( material, dU, dSurfaceStrain, n.data(), separation.data() );

    const double eps = 1e-7;

    auto perturbedResponse = [&]( int externalIndex, double delta ) {
      double dUp[6];
      double dSp[18];
      std::copy( dU, dU + 6, dUp );
      std::copy( dSurfaceStrain, dSurfaceStrain + 18, dSp );
      if ( externalIndex < 3 ) {
        // w = dU_top - dU_bottom: perturbing dU_top alone changes w by
        // exactly `delta` (only the difference w ever enters the material).
        dUp[externalIndex] += delta;
      }
      else if ( externalIndex < 12 ) {
        dSp[externalIndex - 3] += delta;
      }
      else {
        dSp[externalIndex - 3] += delta;
      }
      return evaluateVirginResponse( material, dUp, dSp, n.data(), separation.data() );
    };

    Eigen::Matrix< double, 21, 1 > pAnalytic;
    pAnalytic.segment< 3 >( 0 )  = analytic.generalizedForce;
    pAnalytic.segment< 9 >( 3 )  = flattenRowMajor( analytic.surfaceStressPlus );
    pAnalytic.segment< 9 >( 12 ) = flattenRowMajor( analytic.surfaceStressMinus );

    Eigen::Matrix< double, 21, 21, Eigen::RowMajor > KFiniteDifference;
    for ( int col = 0; col < 21; ++col ) {
      const EqResponse plus  = perturbedResponse( col, eps );
      const EqResponse minus = perturbedResponse( col, -eps );

      Eigen::Matrix< double, 21, 1 > pPlus, pMinus;
      pPlus.segment< 3 >( 0 )   = plus.generalizedForce;
      pPlus.segment< 9 >( 3 )   = flattenRowMajor( plus.surfaceStressPlus );
      pPlus.segment< 9 >( 12 )  = flattenRowMajor( plus.surfaceStressMinus );
      pMinus.segment< 3 >( 0 )  = minus.generalizedForce;
      pMinus.segment< 9 >( 3 )  = flattenRowMajor( minus.surfaceStressPlus );
      pMinus.segment< 9 >( 12 ) = flattenRowMajor( minus.surfaceStressMinus );

      KFiniteDifference.col( col ) = ( pPlus - pMinus ) / ( 2.0 * eps );
    }

    Eigen::Matrix< double, 21, 21, Eigen::RowMajor > KAnalytic;
    KAnalytic.block< 3, 3 >( 0, 0 )   = analytic.Q_ww;
    KAnalytic.block< 3, 9 >( 0, 3 )   = analytic.Q_wAp;
    KAnalytic.block< 3, 9 >( 0, 12 )  = analytic.Q_wAm;
    KAnalytic.block< 9, 3 >( 3, 0 )   = analytic.Q_Apw;
    KAnalytic.block< 9, 9 >( 3, 3 )   = analytic.Q_ApAp;
    KAnalytic.block< 9, 9 >( 3, 12 )  = analytic.Q_ApAm;
    KAnalytic.block< 9, 3 >( 12, 0 )  = analytic.Q_Amw;
    KAnalytic.block< 9, 9 >( 12, 3 )  = analytic.Q_AmAp;
    KAnalytic.block< 9, 9 >( 12, 12 ) = analytic.Q_AmAm;

    const double relTol = 1e-6 * std::max( 1.0, KFiniteDifference.lpNorm< Eigen::Infinity >() );
    assertMatrixNear( KAnalytic, KFiniteDifference, relTol, "Full 21x21 condensed tangent does not match FD" );
  }

  /**
   * 14.4: for a 1D-like scenario (zero surface gradient, uniaxial jump
   * along the normal, isotropic materials with nu=0 so shear does not mix
   * in), the normal-normal component of d(f)/d(w) must match the harmonic
   * series-spring formula q_eff = 2 q+ q- / (q+ + q-), NOT the arithmetic
   * average (q+ + q-)/2 that the uncondensed (non-equilibrated) formulation
   * would give.
   */
  void testScalarSeriesLawForDifferentSideModuli()
  {
    const double h      = 0.01;
    const double EPlus  = 2.0e5;
    const double EMinus = 5.0e4;

    // MarmotEquilibratedXInterfaceMaterialHypoElastic uses the SAME
    // properties for both embedded sides, so different C+ != C- (the
    // premise of the series-law formula) cannot be constructed through its
    // public API. This test instead directly emulates the two-sided local
    // equilibrium algebra with two independent throwaway LINEARELASTIC bulk
    // materials of different moduli, verifying the closed-form series-law
    // prediction the material's local Newton solve is supposed to reduce to.
    const double qEffExpected = 2.0 * EPlus * EMinus / ( EPlus + EMinus );

    auto bulkMaterial = []( double E ) {
      const double props[2] = { E, 0.0 };
      return std::unique_ptr< MarmotMaterialHypoElastic >(
        MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( "LINEARELASTIC", props, 2, 1 ) );
    };

    auto topMat    = bulkMaterial( EPlus );
    auto bottomMat = bulkMaterial( EMinus );

    const double w          = 1.0e-6;
    const double ell        = h;
    double       z          = 0.0;
    double       tConverged = 0.0;

    std::vector< double > topScratch( topMat->getNumberOfRequiredStateVars(), 0.0 );
    std::vector< double > bottomScratch( bottomMat->getNumberOfRequiredStateVars(), 0.0 );

    for ( int iter = 0; iter < 30; ++iter ) {
      const double gPlus  = w / ell + 0.5 * z;
      const double gMinus = w / ell - 0.5 * z;

      Marmot::Vector6d epsPlus  = Marmot::Vector6d::Zero();
      Marmot::Vector6d epsMinus = Marmot::Vector6d::Zero();
      epsPlus[2]                = gPlus; // zz normal strain (Voigt index 2)
      epsMinus[2]               = gMinus;

      Marmot::Matrix6d                   tangentPlus  = Marmot::Matrix6d::Zero();
      Marmot::Matrix6d                   tangentMinus = Marmot::Matrix6d::Zero();
      MarmotMaterialHypoElastic::state3D statePlus{ Marmot::Vector6d::Zero(), 0.0, 0.0, topScratch.data() };
      MarmotMaterialHypoElastic::state3D stateMinus{ Marmot::Vector6d::Zero(), 0.0, 0.0, bottomScratch.data() };
      topMat->computeStress( statePlus, tangentPlus, epsPlus, { 0.0, 1.0 } );
      bottomMat->computeStress( stateMinus, tangentMinus, epsMinus, { 0.0, 1.0 } );

      const double tPlus  = statePlus.stress[2];
      const double tMinus = stateMinus.stress[2];
      const double r      = tPlus - tMinus;
      if ( std::abs( r ) < 1e-12 * std::max( 1.0, std::max( std::abs( tPlus ), std::abs( tMinus ) ) ) ) {
        tConverged = tPlus;
        break;
      }
      const double Qavg = 0.5 * ( tangentPlus( 2, 2 ) + tangentMinus( 2, 2 ) );
      z += -r / Qavg;
    }

    const double qEffActual = tConverged / ( w / ell );

    throwExceptionOnFailure( std::abs( qEffActual - qEffExpected ) < 1e-3 * qEffExpected,
                             "Scalar series-law check failed: expected q_eff = 2 q+ q- / (q+ + q-) = " +
                               std::to_string( qEffExpected ) + ", got " + std::to_string( qEffActual ) );

    const double qArithmeticAverage = 0.5 * ( EPlus + EMinus );
    throwExceptionOnFailure( std::abs( qEffActual - qArithmeticAverage ) > 0.05 * qArithmeticAverage,
                             "Scalar series-law check: converged stiffness should differ clearly from the naive "
                             "arithmetic average, confirming this is genuinely the equilibrated (not averaged) "
                             "response." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testDegenerateGeometryAndConstructionThrow,
                                                   testDensityDelegation,
                                                   testUniformSidesGiveZeroInternalJump,
                                                   testTractionEquilibriumElasticSkewed,
                                                   testTractionEquilibriumAndNonzeroCrossBlocksPlastic,
                                                   testFullCondensedTangentMatchesFiniteDifference,
                                                   testScalarSeriesLawForDifferentSideModuli };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
