#include "Marmot/MarmotCorrectedInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/MarmotXInterfaceMaterialHypoElastic.h"

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

  using XMaterial = MarmotXInterfaceMaterialHypoElastic;

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

  /**
   * All generalized outputs of one X-interface stress update. The four
   * force/surfaceStress fields double as the persistent generalized state
   * between increments, exactly as in the element.
   */
  struct XResponse {
    Eigen::Vector3d   forcePlus          = Eigen::Vector3d::Zero();
    Eigen::Vector3d   forceMinus         = Eigen::Vector3d::Zero();
    Matrix3dRowMajor  surfaceStressPlus  = Matrix3dRowMajor::Zero();
    Matrix3dRowMajor  surfaceStressMinus = Matrix3dRowMajor::Zero();
    Matrix3dRowMajor  Qplus              = Matrix3dRowMajor::Zero();
    Matrix3dRowMajor  Qminus             = Matrix3dRowMajor::Zero();
    Matrix3x9RowMajor Hplus              = Matrix3x9RowMajor::Zero();
    Matrix3x9RowMajor Hminus             = Matrix3x9RowMajor::Zero();
    Matrix9x3RowMajor Kplus              = Matrix9x3RowMajor::Zero();
    Matrix9x3RowMajor Kminus             = Matrix9x3RowMajor::Zero();
    Matrix9dRowMajor  Zplus              = Matrix9dRowMajor::Zero();
    Matrix9dRowMajor  Zminus             = Matrix9dRowMajor::Zero();
  };

  void computeXStress( XMaterial&    material,
                       XResponse&    response,
                       double*       stateVars,
                       const double* dU,
                       const double* dSurfaceStrain,
                       const double* normal,
                       const double* separation,
                       double        timeOld,
                       double        dT )
  {
    XMaterial::State         state{ response.forcePlus.data(),
                            response.forceMinus.data(),
                            response.surfaceStressPlus.data(),
                            response.surfaceStressMinus.data(),
                            stateVars };
    XMaterial::Tangents      tangents{ response.Qplus.data(),
                                  response.Qminus.data(),
                                  response.Hplus.data(),
                                  response.Hminus.data(),
                                  response.Kplus.data(),
                                  response.Kminus.data(),
                                  response.Zplus.data(),
                                  response.Zminus.data() };
    XMaterial::TimeIncrement timeIncrement{ timeOld, dT };

    if ( separation ) {
      XMaterial::Deformation deformation{ dU, dSurfaceStrain, normal, separation };
      material.computeStress( state, tangents, deformation, timeIncrement );
    }
    else {
      XMaterial::Deformation deformation{ dU, dSurfaceStrain, normal };
      material.computeStress( state, tangents, deformation, timeIncrement );
    }
  }

  Eigen::VectorXd makeInitializedStateVars( XMaterial& material )
  {
    Eigen::VectorXd stateVars = Eigen::VectorXd::Zero( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );
    return stateVars;
  }

  /** Single virgin-state evaluation, used by the finite-difference checks. */
  XResponse evaluateVirginResponse( XMaterial&    material,
                                    const double* dU,
                                    const double* dSurfaceStrain,
                                    const double* normal,
                                    const double* separation )
  {
    XResponse       response;
    Eigen::VectorXd stateVars = makeInitializedStateVars( material );
    computeXStress( material, response, stateVars.data(), dU, dSurfaceStrain, normal, separation, 0.0, 1.0 );
    return response;
  }

  Vector9d flattenRowMajor( const Matrix3dRowMajor& tensor )
  {
    return Eigen::Map< const Vector9d >( tensor.data() );
  }

  /**
   * With A+ == A- (no through-thickness asymmetry) the X material must
   * reduce to the same generalized state as the Corrected material: each
   * side sees the same strain, integrates over half the thickness, and
   * together sum back to the single-averaged-gradient response.
   */
  void testUniformLoadingMatchesCorrectedMaterial()
  {
    const double interfaceProperties[3] = { 1e5, 0.3, 0.01 };
    const double normal[3]              = { 0., 0., 1. };

    XMaterial                                   xMaterial( "LINEARELASTIC", interfaceProperties, 3, 1 );
    MarmotCorrectedInterfaceMaterialHypoElastic correctedMaterial( "LINEARELASTIC", interfaceProperties, 3, 1 );

    Eigen::VectorXd xStateVars = makeInitializedStateVars( xMaterial );

    Eigen::VectorXd correctedStateVars = Eigen::VectorXd::Zero( correctedMaterial.getNumberOfRequiredStateVars() );
    correctedMaterial.initializeYourself( correctedStateVars.data(), static_cast< int >( correctedStateVars.size() ) );

    const double dU[6]              = { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 };
    const double surfaceGradient[9] = { 1.0e-4, -2.0e-4, 0.5e-4, 0.8e-4, 1.5e-4, -0.7e-4, 0.3e-4, -1.1e-4, 0.9e-4 };
    double       dSurfaceStrain[18];
    // A+ == A- == surfaceGradient: no through-thickness asymmetry.
    std::copy( surfaceGradient, surfaceGradient + 9, dSurfaceStrain );
    std::copy( surfaceGradient, surfaceGradient + 9, dSurfaceStrain + 9 );

    XResponse xResponse;
    computeXStress( xMaterial, xResponse, xStateVars.data(), dU, dSurfaceStrain, normal, nullptr, 0.0, 1.0 );

    Eigen::Vector3d   correctedForce         = Eigen::Vector3d::Zero();
    Matrix3dRowMajor  correctedSurfaceStress = Matrix3dRowMajor::Zero();
    Matrix3dRowMajor  Q                      = Matrix3dRowMajor::Zero();
    Matrix9dRowMajor  Z                      = Matrix9dRowMajor::Zero();
    Matrix3x9RowMajor H                      = Matrix3x9RowMajor::Zero();
    Matrix9x3RowMajor K                      = Matrix9x3RowMajor::Zero();

    MarmotCorrectedInterfaceMaterialHypoElastic::State    state{ correctedForce.data(),
                                                              correctedSurfaceStress.data(),
                                                              correctedStateVars.data() };
    MarmotCorrectedInterfaceMaterialHypoElastic::Tangents correctedTangents{ Q.data(), Z.data(), H.data(), K.data() };
    MarmotCorrectedInterfaceMaterialHypoElastic::Deformation   correctedDeformation{ dU, dSurfaceStrain, normal };
    MarmotCorrectedInterfaceMaterialHypoElastic::TimeIncrement correctedTimeIncrement{ 0.0, 1.0 };
    correctedMaterial.computeStress( state, correctedTangents, correctedDeformation, correctedTimeIncrement );

    const Eigen::Vector3d  combinedForce         = xResponse.forcePlus + xResponse.forceMinus;
    const Matrix3dRowMajor combinedSurfaceStress = xResponse.surfaceStressPlus + xResponse.surfaceStressMinus;

    assertMatrixNear( combinedForce,
                      correctedForce,
                      1e-12,
                      "Uniform-loading X material combined force differs from Corrected material" );
    assertMatrixNear( combinedSurfaceStress,
                      correctedSurfaceStress,
                      1e-12,
                      "Uniform-loading X material combined surface stress differs from Corrected material" );
  }

  /**
   * The decisive regression/validation test motivating the whole XIQUAD4
   * formulation: a cross-sectional-rotation mode with A_avg = n (x) t and
   * jump = -ell*t is an EXACT zero-strain (zero-energy) mode of the
   * Corrected material (single averaged surface gradient), because
   * G_avg = A_avg + (jump/ell) x n = n(x)t - t(x)n is antisymmetric and
   * therefore vanishes under symmetrization -- regardless of how A_avg is
   * split between the two faces.
   *
   * Choosing the split A+ = 2(n x t), A- = 0 keeps the SAME average
   * (n x t) but gives each side its own, individually nonzero, strain
   * (sym(G+) = -sym(G-) = (n(x)t + t(x)n)/2 -- a genuine bending/curvature
   * pattern). The X material evaluates the two sides' energies separately
   * before summing, so it does NOT see this cancellation.
   */
  void testRotationModeIsZeroForCorrectedButNotForXMaterial()
  {
    const double interfaceProperties[3] = { 1e5, 0.3, 0.01 };
    const double h                      = interfaceProperties[2];

    const Eigen::Vector3d n( 0., 0., 1. );
    const Eigen::Vector3d t( 1., 0., 0. );

    const Matrix3dRowMajor nOuterT = n * t.transpose();

    // Coincident faces: ell falls back to h, d_tau = 0.
    const double normal[3] = { n( 0 ), n( 1 ), n( 2 ) };
    const double dU[6]     = {
      0.,
      0.,
      0.,
      h * t( 0 ),
      h * t( 1 ),
      h * t( 2 ),
    }; // jump = dU_top - dU_bottom = -h*t

    double                         dSurfaceStrainCorrected[18];
    Eigen::Map< Matrix3dRowMajor > dSurfaceStrainCorrectedTop( dSurfaceStrainCorrected );
    Eigen::Map< Matrix3dRowMajor > dSurfaceStrainCorrectedBottom( dSurfaceStrainCorrected + 9 );
    dSurfaceStrainCorrectedTop    = nOuterT; // top surface gradient
    dSurfaceStrainCorrectedBottom = nOuterT; // bottom: A_avg = n(x)t

    double                         dSurfaceStrainX[18];
    Eigen::Map< Matrix3dRowMajor > dSurfaceStrainXPlus( dSurfaceStrainX );
    Eigen::Map< Matrix3dRowMajor > dSurfaceStrainXMinus( dSurfaceStrainX + 9 );
    dSurfaceStrainXPlus  = 2.0 * nOuterT;            // A+
    dSurfaceStrainXMinus = Matrix3dRowMajor::Zero(); // A-

    MarmotCorrectedInterfaceMaterialHypoElastic correctedMaterial( "LINEARELASTIC", interfaceProperties, 3, 1 );
    Eigen::VectorXd correctedStateVars = Eigen::VectorXd::Zero( correctedMaterial.getNumberOfRequiredStateVars() );
    correctedMaterial.initializeYourself( correctedStateVars.data(), static_cast< int >( correctedStateVars.size() ) );

    Eigen::Vector3d   correctedForce         = Eigen::Vector3d::Zero();
    Matrix3dRowMajor  correctedSurfaceStress = Matrix3dRowMajor::Zero();
    Matrix3dRowMajor  Q                      = Matrix3dRowMajor::Zero();
    Matrix9dRowMajor  Z                      = Matrix9dRowMajor::Zero();
    Matrix3x9RowMajor H                      = Matrix3x9RowMajor::Zero();
    Matrix9x3RowMajor K                      = Matrix9x3RowMajor::Zero();

    MarmotCorrectedInterfaceMaterialHypoElastic::State         state{ correctedForce.data(),
                                                              correctedSurfaceStress.data(),
                                                              correctedStateVars.data() };
    MarmotCorrectedInterfaceMaterialHypoElastic::Tangents      tangents{ Q.data(), Z.data(), H.data(), K.data() };
    MarmotCorrectedInterfaceMaterialHypoElastic::Deformation   deformation{ dU, dSurfaceStrainCorrected, normal };
    MarmotCorrectedInterfaceMaterialHypoElastic::TimeIncrement timeIncrement{ 0.0, 1.0 };
    correctedMaterial.computeStress( state, tangents, deformation, timeIncrement );

    throwExceptionOnFailure( correctedForce.template lpNorm< Eigen::Infinity >() < 1e-10,
                             "Regression check failed: the Corrected material should give (near-)zero force for the "
                             "rotation mode (A_avg = n(x)t, jump = -ell*t)." );
    throwExceptionOnFailure( correctedSurfaceStress.template lpNorm< Eigen::Infinity >() < 1e-10,
                             "Regression check failed: the Corrected material should give (near-)zero surface stress "
                             "for the rotation mode (A_avg = n(x)t, jump = -ell*t)." );

    XMaterial       xMaterial( "LINEARELASTIC", interfaceProperties, 3, 1 );
    XResponse       xResponse;
    Eigen::VectorXd xStateVars = makeInitializedStateVars( xMaterial );
    computeXStress( xMaterial, xResponse, xStateVars.data(), dU, dSurfaceStrainX, normal, nullptr, 0.0, 1.0 );

    // Genuine bending signature: the two sides carry equal and opposite
    // surface stress, and neither is (anywhere near) zero on its own.
    throwExceptionOnFailure( xResponse.surfaceStressPlus.template lpNorm< Eigen::Infinity >() > 1.0,
                             "X material should give a clearly nonzero surface stress on side + for the rotation "
                             "mode that the Corrected material misses entirely." );
    assertMatrixNear( xResponse.surfaceStressPlus,
                      ( -xResponse.surfaceStressMinus ).eval(),
                      1e-8 * xResponse.surfaceStressPlus.template lpNorm< Eigen::Infinity >(),
                      "X material rotation-mode response should be equal and opposite between the two sides "
                      "(a bending/curvature pattern), not zero on either side." );

    const Matrix3dRowMajor combinedSurfaceStress = xResponse.surfaceStressPlus + xResponse.surfaceStressMinus;
    throwExceptionOnFailure( combinedSurfaceStress.template lpNorm< Eigen::Infinity >() < 1e-8 * h,
                             "Sanity check failed: the two sides' surface stresses should still cancel in their sum "
                             "(same total membrane force as the Corrected material), only their INDIVIDUAL "
                             "contributions differ." );
  }

  /**
   * Finite-difference verification of all eight analytic tangent blocks,
   * for independently varying A+ and A-, and confirmation that the cross
   * blocks (d surfaceStressPlus / d A-, d surfaceStressMinus / d A+) are
   * exactly zero, as required by the energy split having no A+/A- coupling.
   */
  void testTangentsMatchFiniteDifferencesAndCrossBlocksAreZero()
  {
    const double h                      = 0.01;
    const double ell                    = 0.025;
    const double interfaceProperties[3] = { 1e5, 0.3, h };

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    XMaterial material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double dU[6]     = { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 };
    const double aPlus[9]  = { 1.0e-4, -2.0e-4, 0.5e-4, 0.8e-4, 1.5e-4, -0.7e-4, 0.3e-4, -1.1e-4, 0.9e-4 };
    const double aMinus[9] = { -0.4e-4, 0.9e-4, -0.2e-4, 0.6e-4, -1.0e-4, 0.3e-4, -0.8e-4, 0.2e-4, 0.5e-4 };
    double       dSurfaceStrain[18];
    std::copy( aPlus, aPlus + 9, dSurfaceStrain );
    std::copy( aMinus, aMinus + 9, dSurfaceStrain + 9 );

    const XResponse analytic = evaluateVirginResponse( material, dU, dSurfaceStrain, n.data(), separation.data() );

    const double eps = 1e-6;

    Matrix3dRowMajor  QplusFD = Matrix3dRowMajor::Zero(), QminusFD = Matrix3dRowMajor::Zero();
    Matrix9x3RowMajor KplusFD = Matrix9x3RowMajor::Zero(), KminusFD = Matrix9x3RowMajor::Zero();

    for ( int k = 0; k < 3; ++k ) {
      double dUPlus[6], dUMinus[6];
      std::copy( dU, dU + 6, dUPlus );
      std::copy( dU, dU + 6, dUMinus );
      dUPlus[k] += eps;
      dUMinus[k] -= eps;

      const XResponse plus  = evaluateVirginResponse( material, dUPlus, dSurfaceStrain, n.data(), separation.data() );
      const XResponse minus = evaluateVirginResponse( material, dUMinus, dSurfaceStrain, n.data(), separation.data() );

      QplusFD.col( k )  = ( plus.forcePlus - minus.forcePlus ) / ( 2.0 * eps );
      QminusFD.col( k ) = ( plus.forceMinus - minus.forceMinus ) / ( 2.0 * eps );
      KplusFD.col( k )  = ( flattenRowMajor( plus.surfaceStressPlus ) - flattenRowMajor( minus.surfaceStressPlus ) ) /
                         ( 2.0 * eps );
      KminusFD.col( k ) = ( flattenRowMajor( plus.surfaceStressMinus ) - flattenRowMajor( minus.surfaceStressMinus ) ) /
                          ( 2.0 * eps );
    }

    Matrix3x9RowMajor HplusFD = Matrix3x9RowMajor::Zero(), HminusFD = Matrix3x9RowMajor::Zero();
    Matrix9dRowMajor  ZplusFD = Matrix9dRowMajor::Zero(), ZminusFD = Matrix9dRowMajor::Zero();
    // Cross blocks: d(surfaceStressPlus)/d(A-) and d(surfaceStressMinus)/d(A+), expected to be zero.
    Matrix9dRowMajor crossPlusOnMinusFD = Matrix9dRowMajor::Zero();
    Matrix9dRowMajor crossMinusOnPlusFD = Matrix9dRowMajor::Zero();

    for ( int entry = 0; entry < 9; ++entry ) {
      // Perturb A+ only.
      double dSurfacePlus[18], dSurfaceMinus[18];
      std::copy( dSurfaceStrain, dSurfaceStrain + 18, dSurfacePlus );
      std::copy( dSurfaceStrain, dSurfaceStrain + 18, dSurfaceMinus );
      dSurfacePlus[entry] += eps;
      dSurfaceMinus[entry] -= eps;

      const XResponse plusA  = evaluateVirginResponse( material, dU, dSurfacePlus, n.data(), separation.data() );
      const XResponse minusA = evaluateVirginResponse( material, dU, dSurfaceMinus, n.data(), separation.data() );

      HplusFD.col( entry ) = ( plusA.forcePlus - minusA.forcePlus ) / ( 2.0 * eps );
      ZplusFD.col( entry ) = ( flattenRowMajor( plusA.surfaceStressPlus ) -
                               flattenRowMajor( minusA.surfaceStressPlus ) ) /
                             ( 2.0 * eps );
      crossMinusOnPlusFD.col( entry ) = ( flattenRowMajor( plusA.surfaceStressMinus ) -
                                          flattenRowMajor( minusA.surfaceStressMinus ) ) /
                                        ( 2.0 * eps );

      // Perturb A- only.
      double dSurfacePlus2[18], dSurfaceMinus2[18];
      std::copy( dSurfaceStrain, dSurfaceStrain + 18, dSurfacePlus2 );
      std::copy( dSurfaceStrain, dSurfaceStrain + 18, dSurfaceMinus2 );
      dSurfacePlus2[9 + entry] += eps;
      dSurfaceMinus2[9 + entry] -= eps;

      const XResponse plusB  = evaluateVirginResponse( material, dU, dSurfacePlus2, n.data(), separation.data() );
      const XResponse minusB = evaluateVirginResponse( material, dU, dSurfaceMinus2, n.data(), separation.data() );

      HminusFD.col( entry ) = ( plusB.forceMinus - minusB.forceMinus ) / ( 2.0 * eps );
      ZminusFD.col( entry ) = ( flattenRowMajor( plusB.surfaceStressMinus ) -
                                flattenRowMajor( minusB.surfaceStressMinus ) ) /
                              ( 2.0 * eps );
      crossPlusOnMinusFD.col( entry ) = ( flattenRowMajor( plusB.surfaceStressPlus ) -
                                          flattenRowMajor( minusB.surfaceStressPlus ) ) /
                                        ( 2.0 * eps );
    }

    const auto relativeTolerance = []( const auto& reference ) {
      return 1e-7 * std::max( 1.0, reference.template lpNorm< Eigen::Infinity >() );
    };

    assertMatrixNear( analytic.Qplus, QplusFD, relativeTolerance( QplusFD ), "Analytic Qplus does not match FD" );
    assertMatrixNear( analytic.Qminus, QminusFD, relativeTolerance( QminusFD ), "Analytic Qminus does not match FD" );
    assertMatrixNear( analytic.Kplus, KplusFD, relativeTolerance( KplusFD ), "Analytic Kplus does not match FD" );
    assertMatrixNear( analytic.Kminus, KminusFD, relativeTolerance( KminusFD ), "Analytic Kminus does not match FD" );
    assertMatrixNear( analytic.Hplus, HplusFD, relativeTolerance( HplusFD ), "Analytic Hplus does not match FD" );
    assertMatrixNear( analytic.Hminus, HminusFD, relativeTolerance( HminusFD ), "Analytic Hminus does not match FD" );
    assertMatrixNear( analytic.Zplus, ZplusFD, relativeTolerance( ZplusFD ), "Analytic Zplus does not match FD" );
    assertMatrixNear( analytic.Zminus, ZminusFD, relativeTolerance( ZminusFD ), "Analytic Zminus does not match FD" );

    const double crossTol = 1e-7 * std::max( 1.0, ZplusFD.template lpNorm< Eigen::Infinity >() );
    assertMatrixNear( crossMinusOnPlusFD,
                      Matrix9dRowMajor::Zero(),
                      crossTol,
                      "d(surfaceStressMinus)/d(A+) should be exactly zero (no cross-coupling)" );
    assertMatrixNear( crossPlusOnMinusFD,
                      Matrix9dRowMajor::Zero(),
                      crossTol,
                      "d(surfaceStressPlus)/d(A-) should be exactly zero (no cross-coupling)" );
  }

  template < typename Callable >
  void expectInvalidArgument( Callable&& callable, const std::string& expectedMessage, const std::string& context )
  {
    bool        thrown = false;
    std::string actualMessage;

    try {
      callable();
    }
    catch ( const std::invalid_argument& e ) {
      thrown        = true;
      actualMessage = e.what();
    }

    throwExceptionOnFailure( thrown, context + ": expected std::invalid_argument was not thrown." );
    throwExceptionOnFailure( actualMessage == expectedMessage,
                             context + ": unexpected message '" + actualMessage + "'." );
  }

  void testDegenerateGeometryAndConstructionThrow()
  {
    expectInvalidArgument(
      []() {
        const double tooFewProperties[2] = { 1e5, 0.3 };
        XMaterial    material( "LINEARELASTIC", tooFewProperties, 2, 1 );
      },
      "MarmotXInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h.",
      "Construction with 2 properties" );

    expectInvalidArgument(
      []() {
        const double zeroThickness[3] = { 1e5, 0.3, 0.0 };
        XMaterial    material( "LINEARELASTIC", zeroThickness, 3, 1 );
      },
      "MarmotXInterfaceMaterialHypoElastic requires h > 0.",
      "Construction with h = 0" );

    const double interfaceProperties[3] = { 1e5, 0.3, 0.01 };
    XMaterial    material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double dU[6]              = { 1e-4, 0., 0., 0., 0., 0. };
    const double dSurfaceStrain[18] = { 0. };

    expectInvalidArgument(
      [&]() {
        const double zeroNormal[3] = { 0., 0., 0. };
        evaluateVirginResponse( material, dU, dSurfaceStrain, zeroNormal, nullptr );
      },
      "MarmotXInterfaceMaterialHypoElastic: interface normal is zero.",
      "Stress update with zero interface normal" );

    const double normal[3] = { 0., 0., 1. };

    expectInvalidArgument(
      [&]() {
        const double tangentialOnlySeparation[3] = { 0.003, 0., 0. };
        evaluateVirginResponse( material, dU, dSurfaceStrain, normal, tangentialOnlySeparation );
      },
      "MarmotXInterfaceMaterialHypoElastic: the top-bottom connector must have a positive normal component.",
      "Stress update with purely tangential connector" );
  }

  void testDensityDelegation()
  {
    const double interfaceProperties[9] = { 1e8, 0.3, 0.01, 2e7, 0.25, 6., 1e-4, 1., 2400. };
    XMaterial    material( "LINEARVISCOELASTICWIECHERT", interfaceProperties, 9, 1 );

    throwExceptionOnFailure( checkIfEqual( material.getDensity(), interfaceProperties[8] ),
                             "X interface density delegation failed." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testUniformLoadingMatchesCorrectedMaterial,
                                                   testRotationModeIsZeroForCorrectedButNotForXMaterial,
                                                   testTangentsMatchFiniteDifferencesAndCrossBlocksAreZero,
                                                   testDegenerateGeometryAndConstructionThrow,
                                                   testDensityDelegation };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
