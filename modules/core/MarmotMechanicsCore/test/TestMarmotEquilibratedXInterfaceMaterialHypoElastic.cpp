#include "Marmot/LinearElastic.h"
#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using EqMaterial = MarmotEquilibratedXInterfaceMaterialHypoElastic;

  using Vector9d          = Eigen::Matrix< double, 9, 1 >;
  using Vector21d         = Eigen::Matrix< double, 21, 1 >;
  using Matrix3dRowMajor  = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor  = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;
  using Matrix3x9RowMajor = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;
  using Matrix21dRowMajor = Eigen::Matrix< double, 21, 21, Eigen::RowMajor >;

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

  /// All generalized outputs of one equilibrated interface stress update.
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

  /// Zero connector: the faces coincide.
  const double noSeparation[3] = { 0., 0., 0. };

  void computeEqStress( EqMaterial&                  material,
                        EqResponse&                  response,
                        double*                      stateVars,
                        const double*                dU,
                        const double*                dSurfaceStrain,
                        const double*                normal,
                        const double*                separation,
                        double                       timeOld,
                        double                       dT,
                        EqMaterial::IncrementMemory* memory = nullptr )
  {
    EqMaterial::State         state{ EqMaterial::Tensor3d( 0.0 ),
                             EqMaterial::Tensor33d( 0.0 ),
                             EqMaterial::Tensor33d( 0.0 ),
                             stateVars,
                             memory };
    EqMaterial::Tangents      tangents{ EqMaterial::Tensor33d( 0.0 ),
                                   EqMaterial::Tensor333d( 0.0 ),
                                   EqMaterial::Tensor333d( 0.0 ),
                                   EqMaterial::Tensor333d( 0.0 ),
                                   EqMaterial::Tensor3333d( 0.0 ),
                                   EqMaterial::Tensor3333d( 0.0 ),
                                   EqMaterial::Tensor333d( 0.0 ),
                                   EqMaterial::Tensor3333d( 0.0 ),
                                   EqMaterial::Tensor3333d( 0.0 ) };
    EqMaterial::Deformation   deformation{ EqMaterial::Tensor6d( dU ),
                                         EqMaterial::Tensor18d( dSurfaceStrain ),
                                         EqMaterial::Tensor3d( normal ),
                                         EqMaterial::Tensor3d( separation ) };
    EqMaterial::TimeIncrement timeIncrement{ timeOld, dT };

    material.computeStress( state, tangents, deformation, timeIncrement );

    response.generalizedForce   = Marmot::mapEigenToFastor( state.generalizedForce );
    response.surfaceStressPlus  = Marmot::mapEigenToFastor( state.surfaceStressPlus );
    response.surfaceStressMinus = Marmot::mapEigenToFastor( state.surfaceStressMinus );
    response.Q_ww               = Marmot::mapEigenToFastor( tangents.Q_ww );
    response.Q_wAp              = Marmot::mapEigenToFastor( tangents.Q_wAp );
    response.Q_wAm              = Marmot::mapEigenToFastor( tangents.Q_wAm );
    response.Q_Apw              = Eigen::Map< const Matrix9x3RowMajor >( tangents.Q_Apw.data() );
    response.Q_ApAp             = Marmot::mapEigenToFastor( tangents.Q_ApAp );
    response.Q_ApAm             = Marmot::mapEigenToFastor( tangents.Q_ApAm );
    response.Q_Amw              = Eigen::Map< const Matrix9x3RowMajor >( tangents.Q_Amw.data() );
    response.Q_AmAp             = Marmot::mapEigenToFastor( tangents.Q_AmAp );
    response.Q_AmAm             = Marmot::mapEigenToFastor( tangents.Q_AmAm );
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

  /// The generalized stress p = ( f, S+, S- ) of a response.
  Vector21d generalizedStress( const EqResponse& response )
  {
    Vector21d p;
    p.segment< 3 >( 0 )  = response.generalizedForce;
    p.segment< 9 >( 3 )  = flattenRowMajor( response.surfaceStressPlus );
    p.segment< 9 >( 12 ) = flattenRowMajor( response.surfaceStressMinus );
    return p;
  }

  /// The condensed tangent of a response, assembled over x = ( w, A+, A- ).
  Matrix21dRowMajor condensedTangent( const EqResponse& response )
  {
    Matrix21dRowMajor K;
    K.block< 3, 3 >( 0, 0 )   = response.Q_ww;
    K.block< 3, 9 >( 0, 3 )   = response.Q_wAp;
    K.block< 3, 9 >( 0, 12 )  = response.Q_wAm;
    K.block< 9, 3 >( 3, 0 )   = response.Q_Apw;
    K.block< 9, 9 >( 3, 3 )   = response.Q_ApAp;
    K.block< 9, 9 >( 3, 12 )  = response.Q_ApAm;
    K.block< 9, 3 >( 12, 0 )  = response.Q_Amw;
    K.block< 9, 9 >( 12, 3 )  = response.Q_AmAp;
    K.block< 9, 9 >( 12, 12 ) = response.Q_AmAm;
    return K;
  }

  /**
   * For a fixed increment, reconstruct t+ and t- from the COMMITTED generalized state (f, S+, S-) and assert
   * that the traction equilibrium holds to tight tolerance: sigma^s = (S^s + f (x) d_tau) / sideThickness,
   * t^s = sigma^s n.
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

  template < typename Callable >
  bool throwsInvalidArgument( Callable&& call )
  {
    try {
      call();
    }
    catch ( const std::invalid_argument& ) {
      return true;
    }
    return false;
  }

  void testConstructionAndGeometryThrow()
  {
    const double tooFewProperties[2] = { 1e5, 0.3 };
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { EqMaterial material( "LINEARELASTIC", tooFewProperties, 2, 1 ); } ),
                             "Construction with fewer than three properties must throw." );

    const double zeroThickness[3] = { 1e5, 0.3, 0.0 };
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { EqMaterial material( "LINEARELASTIC", zeroThickness, 3, 1 ); } ),
                             "Construction with h = 0 must throw." );

    const double properties[3] = { 1e5, 0.3, 0.01 };
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { EqMaterial material( "NOT_A_REGISTERED_MATERIAL", properties, 3, 1 ); } ),
                             "Construction with an unregistered base material must throw." );

    EqMaterial material( "LINEARELASTIC", properties, 3, 1 );

    const double dU[6]              = { 1e-4, 0., 0., 0., 0., 0. };
    const double dSurfaceStrain[18] = { 0. };

    const double zeroNormal[3] = { 0., 0., 0. };
    throwExceptionOnFailure( throwsInvalidArgument( [&] {
                               evaluateVirginResponse( material, dU, dSurfaceStrain, zeroNormal, noSeparation );
                             } ),
                             "A zero interface normal must throw." );

    // the connector points against the normal
    const double normal[3]    = { 0., 0., 1. };
    const double backwards[3] = { 0., 0., -0.02 };
    throwExceptionOnFailure( throwsInvalidArgument(
                               [&] { evaluateVirginResponse( material, dU, dSurfaceStrain, normal, backwards ); } ),
                             "A connector without a positive normal component must throw." );
  }

  void testDensityDelegation()
  {
    const double interfaceProperties[9] = { 1e8, 0.3, 0.01, 2e7, 0.25, 6., 1e-4, 1., 2400. };
    EqMaterial   material( "LINEARVISCOELASTICWIECHERT", interfaceProperties, 9, 1 );

    throwExceptionOnFailure( checkIfEqual( material.getDensity(), interfaceProperties[8] ),
                             "Equilibrated interface density delegation failed." );
  }

  void testCoincidentFacesAreAZeroConnector()
  {
    // A zero connector means coincident faces: the normal separation is the thickness h and the tangential
    // separation vanishes, the same as a connector of length h along the normal.
    const double h                      = 0.01;
    const double interfaceProperties[3] = { 1e5, 0.3, h };
    const double n[3]                   = { 0.6, 0.0, 0.8 };
    const double alongNormal[3]         = { h * 0.6, 0.0, h * 0.8 };

    EqMaterial material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double dU[6]              = { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 };
    const double dSurfaceStrain[18] = { 1.0e-4,
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

    const EqResponse coincident = evaluateVirginResponse( material, dU, dSurfaceStrain, n, noSeparation );
    const EqResponse connected  = evaluateVirginResponse( material, dU, dSurfaceStrain, n, alongNormal );

    assertMatrixNear( generalizedStress( coincident ),
                      generalizedStress( connected ),
                      1e-9 * generalizedStress( coincident ).lpNorm< Eigen::Infinity >(),
                      "Coincident faces must give the response of a connector of length h along the normal" );
    assertMatrixNear( condensedTangent( coincident ),
                      condensedTangent( connected ),
                      1e-9 * condensedTangent( coincident ).lpNorm< Eigen::Infinity >(),
                      "Coincident faces must give the tangent of a connector of length h along the normal" );
  }

  /// A+ == A- and identical (virgin) material state on both sides must give a converged internal jump z == 0.
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
    computeEqStress( material, response, stateVars.data(), dU, dSurfaceStrain, normal, noSeparation, 0.0, 1.0 );

    const double* z = material.getStateView( "normalGradientJump", stateVars.data() ).stateLocation;
    throwExceptionOnFailure( std::abs( z[0] ) < 1e-10 && std::abs( z[1] ) < 1e-10 && std::abs( z[2] ) < 1e-10,
                             "Uniform-side loading should converge to a zero internal normal-gradient jump." );

    assertMatrixNear( response.surfaceStressPlus,
                      response.surfaceStressMinus,
                      1e-8 * response.surfaceStressPlus.lpNorm< Eigen::Infinity >(),
                      "Uniform-side loading should give identical surfaceStressPlus/Minus" );
  }

  /// Geometry of the skewed interface used by several tests.
  struct SkewedInterface {
    double          h           = 0.01;
    double          ell         = 0.025;
    Eigen::Vector3d n           = Eigen::Vector3d( 0.6, 0.0, 0.8 );
    Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) + 0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    Eigen::Vector3d separation  = ell * n + dTangential;
  };

  void testTractionEquilibriumElasticSkewed()
  {
    const SkewedInterface geometry;
    const double          interfaceProperties[3] = { 1e5, 0.3, geometry.h };

    EqMaterial material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double dU[6]     = { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 };
    const double aPlus[9]  = { 1.0e-4, -2.0e-4, 0.5e-4, 0.8e-4, 1.5e-4, -0.7e-4, 0.3e-4, -1.1e-4, 0.9e-4 };
    const double aMinus[9] = { -0.4e-4, 0.9e-4, -0.2e-4, 0.6e-4, -1.0e-4, 0.3e-4, -0.8e-4, 0.2e-4, 0.5e-4 };
    double       dSurfaceStrain[18];
    std::copy( aPlus, aPlus + 9, dSurfaceStrain );
    std::copy( aMinus, aMinus + 9, dSurfaceStrain + 9 );

    const EqResponse response = evaluateVirginResponse( material,
                                                        dU,
                                                        dSurfaceStrain,
                                                        geometry.n.data(),
                                                        geometry.separation.data() );
    assertTractionEquilibrium( response,
                               geometry.n,
                               geometry.dTangential,
                               0.5 * geometry.h,
                               "Elastic skewed equilibrium" );
  }

  /**
   * Drive the two sides to divergent von Mises plastic histories (different prior loading on A+ and A-), so that
   * C+ != C- at the state under test, then verify that the traction equilibrium still holds after the second
   * increment and that the cross blocks of the condensed tangent are clearly nonzero.
   */
  void testTractionEquilibriumAndNonzeroCrossBlocksPlastic()
  {
    const SkewedInterface geometry;
    const double          interfaceProperties[8] = { 210000., 0.3, geometry.h, 20., 200., 5., 10., 2400. };

    EqMaterial material( "VONMISES", interfaceProperties, 8, 1 );

    // strongly asymmetric first increment: A+ large, A- small. It is plastic on both faces, and small enough
    // for the return mapping of the von Mises model to converge.
    const double historyScale     = 0.3;
    const double historyDU[6]     = { 2.0e-3, -0.5e-3, 3.0e-3, 0.2e-3, 0.3e-3, -0.4e-3 };
    const double historyAPlus[9]  = { 8.0e-3, -6.0e-3, 4.0e-3, 5.0e-3, 7.0e-3, -3.0e-3, 2.0e-3, -5.0e-3, 6.0e-3 };
    const double historyAMinus[9] = { 0.2e-3, -0.1e-3, 0.05e-3, 0.1e-3, 0.15e-3, -0.05e-3, 0.05e-3, -0.1e-3, 0.1e-3 };
    double       historyScaledDU[6];
    double       historyDSurfaceStrain[18];
    for ( int i = 0; i < 6; ++i )
      historyScaledDU[i] = historyScale * historyDU[i];
    for ( int i = 0; i < 9; ++i ) {
      historyDSurfaceStrain[i]     = historyScale * historyAPlus[i];
      historyDSurfaceStrain[9 + i] = historyScale * historyAMinus[i];
    }

    Eigen::VectorXd stateVars = makeInitializedStateVars( material );
    EqResponse      response;
    computeEqStress( material,
                     response,
                     stateVars.data(),
                     historyScaledDU,
                     historyDSurfaceStrain,
                     geometry.n.data(),
                     geometry.separation.data(),
                     0.0,
                     1.0 );

    const double dU[6]     = { 0.3e-3, -0.1e-3, 0.4e-3, 0.05e-3, 0.05e-3, -0.05e-3 };
    const double aPlus[9]  = { 0.5e-3, -0.3e-3, 0.2e-3, 0.3e-3, 0.4e-3, -0.2e-3, 0.1e-3, -0.3e-3, 0.3e-3 };
    const double aMinus[9] = { 0.05e-3, -0.02e-3, 0.01e-3, 0.02e-3, 0.03e-3, -0.01e-3, 0.01e-3, -0.02e-3, 0.02e-3 };
    double       dSurfaceStrain[18];
    std::copy( aPlus, aPlus + 9, dSurfaceStrain );
    std::copy( aMinus, aMinus + 9, dSurfaceStrain + 9 );

    computeEqStress( material,
                     response,
                     stateVars.data(),
                     dU,
                     dSurfaceStrain,
                     geometry.n.data(),
                     geometry.separation.data(),
                     1.0,
                     1.0 );

    assertTractionEquilibrium( response,
                               geometry.n,
                               geometry.dTangential,
                               0.5 * geometry.h,
                               "Plastic divergent-history equilibrium" );

    const double crossNorm = std::max( response.Q_ApAm.lpNorm< Eigen::Infinity >(),
                                       response.Q_AmAp.lpNorm< Eigen::Infinity >() );
    const double diagNorm  = std::max( response.Q_ApAp.lpNorm< Eigen::Infinity >(),
                                      response.Q_AmAm.lpNorm< Eigen::Infinity >() );
    throwExceptionOnFailure( crossNorm > 1e-3 * diagNorm,
                             "With divergent plastic histories (C+ != C-), the condensed cross blocks "
                             "d(S+)/d(A-) and d(S-)/d(A+) should be clearly nonzero." );
  }

  /**
   * Perturb all 21 external variables (w, A+, A-) and compare every analytic condensed tangent block against a
   * central finite difference, re-solving the local traction-equilibrium problem at each perturbation.
   */
  void testFullCondensedTangentMatchesFiniteDifference()
  {
    const SkewedInterface geometry;
    const double          interfaceProperties[3] = { 1e5, 0.3, geometry.h };

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

    const EqResponse analytic = evaluateVirginResponse( material,
                                                        dU,
                                                        dSurfaceStrain,
                                                        geometry.n.data(),
                                                        geometry.separation.data() );

    const double eps = 1e-7;

    auto perturbedResponse = [&]( int externalIndex, double delta ) {
      double dUp[6];
      double dSp[18];
      std::copy( dU, dU + 6, dUp );
      std::copy( dSurfaceStrain, dSurfaceStrain + 18, dSp );
      if ( externalIndex < 3 ) {
        // w = dU_top - dU_bottom: perturbing dU_top alone changes w by exactly delta
        dUp[externalIndex] += delta;
      }
      else {
        dSp[externalIndex - 3] += delta;
      }
      return evaluateVirginResponse( material, dUp, dSp, geometry.n.data(), geometry.separation.data() );
    };

    Matrix21dRowMajor KFiniteDifference;
    for ( int col = 0; col < 21; ++col ) {
      const EqResponse plus  = perturbedResponse( col, eps );
      const EqResponse minus = perturbedResponse( col, -eps );

      KFiniteDifference.col( col ) = ( generalizedStress( plus ) - generalizedStress( minus ) ) / ( 2.0 * eps );
    }

    const double relTol = 1e-6 * std::max( 1.0, KFiniteDifference.lpNorm< Eigen::Infinity >() );
    assertMatrixNear( condensedTangent( analytic ),
                      KFiniteDifference,
                      relTol,
                      "Full 21x21 condensed tangent does not match the finite difference" );
  }

  /// Parameters of the loading of the gauge tests: the top face slides along x, the interface normal is z.
  struct SlidingLoading {
    double normal[3]          = { 0., 0., 1. };
    double dSurfaceStrain[18] = { 0. };
    double separation[3]      = { 0., 0., 0.025 };

    /// sliding of the top face of an increment
    static void jump( double slide, double* dU )
    {
      std::fill( dU, dU + 6, 0.0 );
      dU[0] = slide;
    }
  };

  /// Near-perfectly plastic von Mises material: yield stress 20, hardening modulus 1e-3.
  const double perfectPlasticProperties[8] = { 210000., 0.3, 0.01, 20., 1e-3, 0., 0., 2400. };

  /// The same material with a considerable hardening.
  const double hardeningProperties[8] = { 210000., 0.3, 0.01, 20., 200., 5., 10., 2400. };

  /**
   * Pure sliding along the interface is the loading at which both faces flow in the same direction and the
   * averaged acoustic tensor loses its transverse stiffness. The gauge is then active; with a considerable
   * hardening it is not.
   */
  void testGaugeIsActiveForPlasticSlidingOnly()
  {
    SlidingLoading loading;
    double         dU[6];
    SlidingLoading::jump( 5e-6, dU );

    for ( const bool perfectPlasticity : { true, false } ) {
      EqMaterial material( "VONMISES", perfectPlasticity ? perfectPlasticProperties : hardeningProperties, 8, 1 );

      Eigen::VectorXd             stateVars = makeInitializedStateVars( material );
      EqMaterial::IncrementMemory memory;
      EqResponse                  response;
      computeEqStress( material,
                       response,
                       stateVars.data(),
                       dU,
                       loading.dSurfaceStrain,
                       loading.normal,
                       loading.separation,
                       0.0,
                       1.0,
                       &memory );

      throwExceptionOnFailure( memory.gaugeActive == perfectPlasticity,
                               perfectPlasticity
                                 ? "The gauge must be active for plastic sliding with a vanishing hardening."
                                 : "The gauge must be inactive for a hardening material." );
      throwExceptionOnFailure( generalizedStress( response ).allFinite() && condensedTangent( response ).allFinite(),
                               "The response with the gauge must be finite." );
    }
  }

  void testGaugeDecisionIsFrozenWithinAnIncrement()
  {
    SlidingLoading loading;
    double         loadingDU[6];
    double         unloadingDU[6];
    SlidingLoading::jump( 5e-6, loadingDU );
    SlidingLoading::jump( -5e-9, unloadingDU );

    EqMaterial material( "VONMISES", perfectPlasticProperties, 8, 1 );

    // the first increment [0, 1] is plastic and is committed
    Eigen::VectorXd committed = makeInitializedStateVars( material );
    {
      EqResponse response;
      computeEqStress( material,
                       response,
                       committed.data(),
                       loadingDU,
                       loading.dSurfaceStrain,
                       loading.normal,
                       loading.separation,
                       0.0,
                       1.0 );
    }

    // evaluations of the second increment [1, 2], each from the committed state, as the iterations of a global
    // Newton scheme are
    const auto evaluate = [&]( const double* dU, double timeOld, EqMaterial::IncrementMemory* memory ) {
      Eigen::VectorXd state = committed;
      EqResponse      response;
      computeEqStress( material,
                       response,
                       state.data(),
                       dU,
                       loading.dSurfaceStrain,
                       loading.normal,
                       loading.separation,
                       timeOld,
                       1.0,
                       memory );
      return response;
    };

    // an iteration that loads plastically, followed by one that unloads elastically: the unloading would not
    // activate the gauge by itself
    EqMaterial::IncrementMemory memory;
    evaluate( loadingDU, 1.0, &memory );
    throwExceptionOnFailure( memory.gaugeActive, "The plastic evaluation must activate the gauge." );

    const EqResponse frozen = evaluate( unloadingDU, 1.0, &memory );
    throwExceptionOnFailure( memory.gaugeActive, "The decision must stay active for the rest of the increment." );

    // Deciding anew, the unloading evaluation does not use the gauge. The frozen decision condenses through the
    // two-dimensional subspace of the gauge instead, which changes the tangent beyond rounding.
    EqMaterial::IncrementMemory freshMemory;
    const EqResponse            decidedAnew = evaluate( unloadingDU, 1.0, &freshMemory );
    throwExceptionOnFailure( !freshMemory.gaugeActive,
                             "The unloading evaluation must not activate the gauge by itself." );
    throwExceptionOnFailure( ( condensedTangent( frozen ) - condensedTangent( decidedAnew ) ).norm() >
                               1e-6 * condensedTangent( decidedAnew ).norm(),
                             "Freezing the decision must change the condensed tangent of the unloading evaluation." );

    // in the next increment the decision is taken anew, and the unloading evaluation does not activate the gauge
    evaluate( unloadingDU, 2.0, &memory );
    throwExceptionOnFailure( !memory.gaugeActive, "The decision must be renewed in the next increment." );
  }

  /// Linear elastic material that records the characteristic element length it is evaluated with.
  class LengthRecordingElastic : public Marmot::Materials::LinearElastic {
  public:
    using Marmot::Materials::LinearElastic::LinearElastic;

    /// the lengths of all evaluations
    static std::vector< double >& evaluatedLengths()
    {
      static std::vector< double > lengths;
      return lengths;
    }

    void computeStress( state3D&                state,
                        Marmot::Matrix6d&       dStressDDStrain,
                        const Marmot::Vector6d& dStrain,
                        const timeInfo&         timeInfo ) const override
    {
      evaluatedLengths().push_back( characteristicElementLength );
      Marmot::Materials::LinearElastic::computeStress( state, dStressDDStrain, dStrain, timeInfo );
    }
  };

  /// The base materials of both faces have to receive the characteristic length of the quadrature point, as the
  /// regularization of softening base materials depends on it.
  void testCharacteristicLengthReachesBothFaces()
  {
    static const bool
      registered = MarmotLibrary::MarmotMaterialHypoElasticFactory::registerMaterial< LengthRecordingElastic >(
        "LENGTHRECORDINGELASTIC" );
    throwExceptionOnFailure( registered, "The recording material could not be registered." );

    const double interfaceProperties[3] = { 1e5, 0.3, 0.01 };
    EqMaterial   material( "LENGTHRECORDINGELASTIC", interfaceProperties, 3, 1 );

    material.setCharacteristicElementLength( 0.37 );
    throwExceptionOnFailure( checkIfEqual( material.characteristicElementLength, 0.37 ),
                             "The material must keep the characteristic length." );

    const double normal[3]          = { 0., 0., 1. };
    const double dU[6]              = { 1e-4, 0., 0., 0., 0., 0. };
    const double dSurfaceStrain[18] = { 0. };

    LengthRecordingElastic::evaluatedLengths().clear();
    evaluateVirginResponse( material, dU, dSurfaceStrain, normal, noSeparation );

    // the local solve evaluates both faces, at least once each
    const auto& lengths = LengthRecordingElastic::evaluatedLengths();
    throwExceptionOnFailure( lengths.size() >= 2 && lengths.size() % 2 == 0,
                             "Both faces must be evaluated, and equally often." );
    throwExceptionOnFailure( std::all_of( lengths.begin(),
                                          lengths.end(),
                                          []( double length ) { return length == 0.37; } ),
                             "Every evaluation of the base materials must see the characteristic length." );
  }

  void testIncrementMemoryDecision()
  {
    EqMaterial::IncrementMemory memory;

    throwExceptionOnFailure( std::isnan( memory.timeOld ), "A new memory has taken no decision yet." );
    throwExceptionOnFailure( memory.decide( 0.0, true ), "The first decision of an increment is taken." );
    throwExceptionOnFailure( memory.decide( 0.0, false ), "The first decision holds for the whole increment." );
    throwExceptionOnFailure( !memory.decide( 1.0, false ), "A new increment takes a new decision." );
    throwExceptionOnFailure( !memory.decide( 1.0, true ), "The decision of the new increment holds." );
    throwExceptionOnFailure( memory.decide( 0.0, true ),
                             "Returning to an earlier time, as after a cutback of the step, takes a new decision." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testConstructionAndGeometryThrow,
                                                   testDensityDelegation,
                                                   testCoincidentFacesAreAZeroConnector,
                                                   testUniformSidesGiveZeroInternalJump,
                                                   testTractionEquilibriumElasticSkewed,
                                                   testTractionEquilibriumAndNonzeroCrossBlocksPlastic,
                                                   testFullCondensedTangentMatchesFiniteDifference,
                                                   testGaugeIsActiveForPlasticSlidingOnly,
                                                   testGaugeDecisionIsFrozenWithinAnIncrement,
                                                   testCharacteristicLengthReachesBothFaces,
                                                   testIncrementMemoryDecision };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
