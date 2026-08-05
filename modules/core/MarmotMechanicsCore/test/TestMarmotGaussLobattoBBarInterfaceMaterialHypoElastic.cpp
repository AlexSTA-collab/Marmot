#include "Marmot/MarmotGaussLobattoBBarInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotGaussLobattoInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Vector9d          = Eigen::Matrix< double, 9, 1 >;
  using Matrix3dRowMajor  = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor  = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
  using Matrix3x9RowMajor = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;

  using GLMaterial     = MarmotGaussLobattoInterfaceMaterialHypoElastic;
  using GLBBarMaterial = MarmotGaussLobattoBBarInterfaceMaterialHypoElastic;

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

  struct GLResponse {
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

  void computeGLStress( GLMaterial&   material,
                        GLResponse&   response,
                        double*       stateVars,
                        const double* dU,
                        const double* dSurfaceStrain,
                        const double* normal,
                        const double* separation,
                        double        timeOld,
                        double        dT )
  {
    GLMaterial::State         state{ response.generalizedForce.data(),
                             response.surfaceStressPlus.data(),
                             response.surfaceStressMinus.data(),
                             stateVars };
    GLMaterial::Tangents      tangents{ response.Q_ww.data(),
                                   response.Q_wAp.data(),
                                   response.Q_wAm.data(),
                                   response.Q_Apw.data(),
                                   response.Q_ApAp.data(),
                                   response.Q_ApAm.data(),
                                   response.Q_Amw.data(),
                                   response.Q_AmAp.data(),
                                   response.Q_AmAm.data() };
    GLMaterial::TimeIncrement timeIncrement{ timeOld, dT };
    GLMaterial::Deformation   deformation{ dU, dSurfaceStrain, normal, separation };
    material.computeStress( state, tangents, deformation, timeIncrement );
  }

  void computeGLBBarStress( GLBBarMaterial& material,
                            GLResponse&     response,
                            double*         stateVars,
                            const double*   dU,
                            const double*   dSurfaceStrain,
                            const double*   normal,
                            const double*   separation,
                            const double*   traceCorrection,
                            double          timeOld,
                            double          dT )
  {
    GLBBarMaterial::State         state{ response.generalizedForce.data(),
                                 response.surfaceStressPlus.data(),
                                 response.surfaceStressMinus.data(),
                                 stateVars };
    GLBBarMaterial::Tangents      tangents{ response.Q_ww.data(),
                                       response.Q_wAp.data(),
                                       response.Q_wAm.data(),
                                       response.Q_Apw.data(),
                                       response.Q_ApAp.data(),
                                       response.Q_ApAm.data(),
                                       response.Q_Amw.data(),
                                       response.Q_AmAp.data(),
                                       response.Q_AmAm.data() };
    GLBBarMaterial::TimeIncrement timeIncrement{ timeOld, dT };
    GLBBarMaterial::Deformation   deformation{ dU, dSurfaceStrain, normal, separation, traceCorrection };
    material.computeStress( state, tangents, deformation, timeIncrement );
  }

  Vector9d flattenRowMajor( const Matrix3dRowMajor& tensor )
  {
    return Eigen::Map< const Vector9d >( tensor.data() );
  }

  // ------------------------------------------------------------------
  // Zero traceCorrection must reduce EXACTLY to the unmodified GLIQUAD4
  // material -- the core consistency check that the B-bar modification is
  // a pure additive correction with no other side effects.
  // ------------------------------------------------------------------
  void testZeroCorrectionMatchesOriginalElastic()
  {
    const double h                      = 0.01;
    const double ell                    = 0.025;
    const double interfaceProperties[3] = { 1e5, 0.3, h };

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    GLMaterial     glMaterial( "LINEARELASTIC", interfaceProperties, 3, 1 );
    GLBBarMaterial bbarMaterial( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double dU[6]     = { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 };
    const double aPlus[9]  = { 1.0e-4, -2.0e-4, 0.5e-4, 0.8e-4, 1.5e-4, -0.7e-4, 0.3e-4, -1.1e-4, 0.9e-4 };
    const double aMinus[9] = { -0.4e-4, 0.9e-4, -0.2e-4, 0.6e-4, -1.0e-4, 0.3e-4, -0.8e-4, 0.2e-4, 0.5e-4 };
    double       dSurfaceStrain[18];
    std::copy( aPlus, aPlus + 9, dSurfaceStrain );
    std::copy( aMinus, aMinus + 9, dSurfaceStrain + 9 );

    Eigen::VectorXd glStateVars = Eigen::VectorXd::Zero( glMaterial.getNumberOfRequiredStateVars() );
    glMaterial.initializeYourself( glStateVars.data(), static_cast< int >( glStateVars.size() ) );
    GLResponse glResponse;
    computeGLStress( glMaterial,
                     glResponse,
                     glStateVars.data(),
                     dU,
                     dSurfaceStrain,
                     n.data(),
                     separation.data(),
                     0.0,
                     1.0 );

    Eigen::VectorXd bbarStateVars = Eigen::VectorXd::Zero( bbarMaterial.getNumberOfRequiredStateVars() );
    bbarMaterial.initializeYourself( bbarStateVars.data(), static_cast< int >( bbarStateVars.size() ) );
    const double zeroCorrection[5] = { 0., 0., 0., 0., 0. };
    GLResponse   bbarResponse;
    computeGLBBarStress( bbarMaterial,
                         bbarResponse,
                         bbarStateVars.data(),
                         dU,
                         dSurfaceStrain,
                         n.data(),
                         separation.data(),
                         zeroCorrection,
                         0.0,
                         1.0 );

    const double tol = 1e-10 * std::max( 1.0, glResponse.generalizedForce.norm() );
    assertMatrixNear( bbarResponse.generalizedForce,
                      glResponse.generalizedForce,
                      tol,
                      "zero-correction: generalizedForce" );
    assertMatrixNear( bbarResponse.surfaceStressPlus,
                      glResponse.surfaceStressPlus,
                      tol,
                      "zero-correction: surfaceStressPlus" );
    assertMatrixNear( bbarResponse.surfaceStressMinus,
                      glResponse.surfaceStressMinus,
                      tol,
                      "zero-correction: surfaceStressMinus" );
    assertMatrixNear( bbarResponse.Q_ApAp, glResponse.Q_ApAp, tol, "zero-correction: Q_ApAp" );
    assertMatrixNear( bbarResponse.Q_ApAm, glResponse.Q_ApAm, tol, "zero-correction: Q_ApAm" );
    assertMatrixNear( bbarResponse.Q_AmAp, glResponse.Q_AmAp, tol, "zero-correction: Q_AmAp" );
    assertMatrixNear( bbarResponse.Q_AmAm, glResponse.Q_AmAm, tol, "zero-correction: Q_AmAm" );
  }

  void testZeroCorrectionMatchesOriginalPlasticWithHistory()
  {
    const double h                      = 0.01;
    const double ell                    = 0.025;
    const double interfaceProperties[8] = { 210000., 0.3, h, 20., 200., 5., 10., 2400. };

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    GLMaterial     glMaterial( "VONMISES", interfaceProperties, 8, 1 );
    GLBBarMaterial bbarMaterial( "VONMISES", interfaceProperties, 8, 1 );

    const double historyDU[6]     = { 2.0e-3, -0.5e-3, 3.0e-3, 0.2e-3, 0.3e-3, -0.4e-3 };
    const double historyAPlus[9]  = { 8.0e-3, -6.0e-3, 4.0e-3, 5.0e-3, 7.0e-3, -3.0e-3, 2.0e-3, -5.0e-3, 6.0e-3 };
    const double historyAMinus[9] = { 0.2e-3, -0.1e-3, 0.05e-3, 0.1e-3, 0.15e-3, -0.05e-3, 0.05e-3, -0.1e-3, 0.1e-3 };
    double       historyDS[18];
    std::copy( historyAPlus, historyAPlus + 9, historyDS );
    std::copy( historyAMinus, historyAMinus + 9, historyDS + 9 );

    const double zeroCorrection[5] = { 0., 0., 0., 0., 0. };

    Eigen::VectorXd glStateVars = Eigen::VectorXd::Zero( glMaterial.getNumberOfRequiredStateVars() );
    glMaterial.initializeYourself( glStateVars.data(), static_cast< int >( glStateVars.size() ) );
    GLResponse glHist;
    computeGLStress( glMaterial,
                     glHist,
                     glStateVars.data(),
                     historyDU,
                     historyDS,
                     n.data(),
                     separation.data(),
                     0.0,
                     1.0 );

    Eigen::VectorXd bbarStateVars = Eigen::VectorXd::Zero( bbarMaterial.getNumberOfRequiredStateVars() );
    bbarMaterial.initializeYourself( bbarStateVars.data(), static_cast< int >( bbarStateVars.size() ) );
    GLResponse bbarHist;
    computeGLBBarStress( bbarMaterial,
                         bbarHist,
                         bbarStateVars.data(),
                         historyDU,
                         historyDS,
                         n.data(),
                         separation.data(),
                         zeroCorrection,
                         0.0,
                         1.0 );

    const double dU[6]     = { 0.3e-3, -0.1e-3, 0.4e-3, 0.05e-3, 0.05e-3, -0.05e-3 };
    const double aPlus[9]  = { 0.5e-3, -0.3e-3, 0.2e-3, 0.3e-3, 0.4e-3, -0.2e-3, 0.1e-3, -0.3e-3, 0.3e-3 };
    const double aMinus[9] = { 0.05e-3, -0.02e-3, 0.01e-3, 0.02e-3, 0.03e-3, -0.01e-3, 0.01e-3, -0.02e-3, 0.02e-3 };
    double       dS[18];
    std::copy( aPlus, aPlus + 9, dS );
    std::copy( aMinus, aMinus + 9, dS + 9 );

    GLResponse glResponse;
    computeGLStress( glMaterial, glResponse, glStateVars.data(), dU, dS, n.data(), separation.data(), 0.0, 1.0 );

    GLResponse bbarResponse;
    computeGLBBarStress( bbarMaterial,
                         bbarResponse,
                         bbarStateVars.data(),
                         dU,
                         dS,
                         n.data(),
                         separation.data(),
                         zeroCorrection,
                         0.0,
                         1.0 );

    const double tol = 1e-9 * std::max( 1.0, glResponse.generalizedForce.norm() );
    assertMatrixNear( bbarResponse.generalizedForce,
                      glResponse.generalizedForce,
                      tol,
                      "zero-correction plastic: generalizedForce" );
    assertMatrixNear( bbarResponse.surfaceStressPlus,
                      glResponse.surfaceStressPlus,
                      tol,
                      "zero-correction plastic: surfaceStressPlus" );
    assertMatrixNear( bbarResponse.surfaceStressMinus,
                      glResponse.surfaceStressMinus,
                      tol,
                      "zero-correction plastic: surfaceStressMinus" );
    assertMatrixNear( bbarResponse.Q_ApAm, glResponse.Q_ApAm, tol, "zero-correction plastic: Q_ApAm (cross-block)" );
  }

  // ------------------------------------------------------------------
  // Nonzero traceCorrection: traction equilibrium across stations must
  // still hold (the local Newton / static condensation machinery is
  // untouched), and a nonzero correction must actually CHANGE the
  // response relative to zero correction (sanity that it does something).
  // ------------------------------------------------------------------
  void testNonzeroCorrectionChangesResponseAndKeepsEquilibrium()
  {
    const double h                      = 0.01;
    const double ell                    = 0.025;
    const double interfaceProperties[3] = { 1e5, 0.3, h };

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    GLBBarMaterial material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double dU[6]     = { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 };
    const double aPlus[9]  = { 1.0e-4, -2.0e-4, 0.5e-4, 0.8e-4, 1.5e-4, -0.7e-4, 0.3e-4, -1.1e-4, 0.9e-4 };
    const double aMinus[9] = { -0.4e-4, 0.9e-4, -0.2e-4, 0.6e-4, -1.0e-4, 0.3e-4, -0.8e-4, 0.2e-4, 0.5e-4 };
    double       dSurfaceStrain[18];
    std::copy( aPlus, aPlus + 9, dSurfaceStrain );
    std::copy( aMinus, aMinus + 9, dSurfaceStrain + 9 );

    const double zeroCorrection[5]    = { 0., 0., 0., 0., 0. };
    const double nonzeroCorrection[5] = { 2.0e-5, -1.0e-5, 0.5e-5, 1.5e-5, -0.8e-5 };

    Eigen::VectorXd stateVarsZero = Eigen::VectorXd::Zero( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVarsZero.data(), static_cast< int >( stateVarsZero.size() ) );
    GLResponse responseZero;
    computeGLBBarStress( material,
                         responseZero,
                         stateVarsZero.data(),
                         dU,
                         dSurfaceStrain,
                         n.data(),
                         separation.data(),
                         zeroCorrection,
                         0.0,
                         1.0 );

    Eigen::VectorXd stateVarsNonzero = Eigen::VectorXd::Zero( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVarsNonzero.data(), static_cast< int >( stateVarsNonzero.size() ) );
    GLResponse responseNonzero;
    computeGLBBarStress( material,
                         responseNonzero,
                         stateVarsNonzero.data(),
                         dU,
                         dSurfaceStrain,
                         n.data(),
                         separation.data(),
                         nonzeroCorrection,
                         0.0,
                         1.0 );

    const double diff = ( responseNonzero.generalizedForce - responseZero.generalizedForce ).norm();
    throwExceptionOnFailure( diff > 1e-6, "Nonzero traceCorrection should visibly change the generalized force." );

    // Traction equilibrium across all 5 stations must still hold exactly
    // (unaffected by the correction, since it only shifts the volumetric
    // part of A^(alpha), not the local Newton/condensation machinery).
    for ( int alpha = 0; alpha < 5; ++alpha ) {
      const double* stress9 = material
                                .getStateView( "stationStress" + std::to_string( alpha ), stateVarsNonzero.data() )
                                .stateLocation;
      const Eigen::Map< const Matrix3dRowMajor > sigma( stress9 );
      const Eigen::Vector3d                      t = sigma * n;
      const double* stress0 = material.getStateView( "stationStress0", stateVarsNonzero.data() ).stateLocation;
      const Eigen::Map< const Matrix3dRowMajor > sigma0( stress0 );
      const Eigen::Vector3d                      t0    = sigma0 * n;
      const double                               scale = std::max( 1.0, t0.norm() );
      throwExceptionOnFailure( ( t - t0 ).norm() < 1e-8 * scale,
                               "Nonzero-correction traction equilibrium violated at station " +
                                 std::to_string( alpha ) );
    }
  }

  // ------------------------------------------------------------------
  // Full condensed tangent vs finite difference, with a FIXED nonzero
  // traceCorrection held constant across all perturbations -- the
  // "frozen" B-bar tangent, matching the classical strain-projection
  // method (the correction is not itself differentiated with respect
  // to q).
  // ------------------------------------------------------------------
  void testFullCondensedTangentMatchesFiniteDifferenceFixedCorrection( const std::string&           materialName,
                                                                       const std::vector< double >& properties,
                                                                       bool                         withHistory )
  {
    const double ell = 0.025;

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    GLBBarMaterial material( materialName, properties.data(), static_cast< int >( properties.size() ), 1 );

    Eigen::VectorXd stateVars = Eigen::VectorXd::Zero( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );

    const double zeroCorrection[5] = { 0., 0., 0., 0., 0. };
    if ( withHistory ) {
      const double historyDU[6]     = { 2.0e-3, -0.5e-3, 3.0e-3, 0.2e-3, 0.3e-3, -0.4e-3 };
      const double historyAPlus[9]  = { 8.0e-3, -6.0e-3, 4.0e-3, 5.0e-3, 7.0e-3, -3.0e-3, 2.0e-3, -5.0e-3, 6.0e-3 };
      const double historyAMinus[9] = { 0.2e-3, -0.1e-3, 0.05e-3, 0.1e-3, 0.15e-3, -0.05e-3, 0.05e-3, -0.1e-3, 0.1e-3 };
      double       historyDS[18];
      std::copy( historyAPlus, historyAPlus + 9, historyDS );
      std::copy( historyAMinus, historyAMinus + 9, historyDS + 9 );
      GLResponse histResponse;
      computeGLBBarStress( material,
                           histResponse,
                           stateVars.data(),
                           historyDU,
                           historyDS,
                           n.data(),
                           separation.data(),
                           zeroCorrection,
                           0.0,
                           1.0 );
    }

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
    if ( withHistory ) {
      for ( int i = 0; i < 6; ++i )
        dU[i] *= 0.1;
      for ( int i = 0; i < 18; ++i )
        dSurfaceStrain[i] *= 0.1;
    }

    const double fixedCorrection[5] = { 1.5e-5, -0.7e-5, 0.3e-5, 0.9e-5, -1.2e-5 };

    Eigen::VectorXd baseStateVars = stateVars;
    GLResponse      analytic;
    computeGLBBarStress( material,
                         analytic,
                         baseStateVars.data(),
                         dU,
                         dSurfaceStrain,
                         n.data(),
                         separation.data(),
                         fixedCorrection,
                         0.0,
                         1.0 );

    const double eps = 1e-7;

    auto perturbedResponse = [&]( int externalIndex, double delta ) {
      double dUp[6];
      double dSp[18];
      std::copy( dU, dU + 6, dUp );
      std::copy( dSurfaceStrain, dSurfaceStrain + 18, dSp );
      if ( externalIndex < 3 ) {
        dUp[externalIndex] += delta;
      }
      else {
        dSp[externalIndex - 3] += delta;
      }
      Eigen::VectorXd perturbedStateVars = stateVars;
      GLResponse      response;
      computeGLBBarStress( material,
                           response,
                           perturbedStateVars.data(),
                           dUp,
                           dSp,
                           n.data(),
                           separation.data(),
                           fixedCorrection,
                           0.0,
                           1.0 );
      return response;
    };

    Eigen::Matrix< double, 21, 1 > pAnalytic;
    pAnalytic.segment< 3 >( 0 )  = analytic.generalizedForce;
    pAnalytic.segment< 9 >( 3 )  = flattenRowMajor( analytic.surfaceStressPlus );
    pAnalytic.segment< 9 >( 12 ) = flattenRowMajor( analytic.surfaceStressMinus );

    Eigen::Matrix< double, 21, 21, Eigen::RowMajor > KFiniteDifference;
    for ( int col = 0; col < 21; ++col ) {
      const GLResponse plus  = perturbedResponse( col, eps );
      const GLResponse minus = perturbedResponse( col, -eps );

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
    assertMatrixNear( KAnalytic,
                      KFiniteDifference,
                      relTol,
                      std::string( "Full 21x21 GL-BBar condensed tangent (fixed correction) vs FD (" ) + materialName +
                        ( withHistory ? ", with history)" : ", elastic)" ) );
  }

  void testFullCondensedTangentElasticSkewedFixedCorrection()
  {
    testFullCondensedTangentMatchesFiniteDifferenceFixedCorrection( "LINEARELASTIC", { 1e5, 0.3, 0.01 }, false );
  }

  void testFullCondensedTangentPlasticDivergentHistoryFixedCorrection()
  {
    testFullCondensedTangentMatchesFiniteDifferenceFixedCorrection( "VONMISES",
                                                                    { 210000., 0.3, 0.01, 20., 200., 5., 10., 2400. },
                                                                    true );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testZeroCorrectionMatchesOriginalElastic,
                                                   testZeroCorrectionMatchesOriginalPlasticWithHistory,
                                                   testNonzeroCorrectionChangesResponseAndKeepsEquilibrium,
                                                   testFullCondensedTangentElasticSkewedFixedCorrection,
                                                   testFullCondensedTangentPlasticDivergentHistoryFixedCorrection };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
