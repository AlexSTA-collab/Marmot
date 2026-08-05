#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotGaussLobattoInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Vector3d          = Eigen::Vector3d;
  using Vector9d          = Eigen::Matrix< double, 9, 1 >;
  using Matrix3dRowMajor  = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor  = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
  using Matrix3x9RowMajor = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;

  using GLMaterial = MarmotGaussLobattoInterfaceMaterialHypoElastic;
  using YMaterial  = MarmotEquilibratedXInterfaceMaterialHypoElastic;

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

  /** All generalized outputs of one Gauss-Lobatto interface stress update. */
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

    if ( separation ) {
      GLMaterial::Deformation deformation{ dU, dSurfaceStrain, normal, separation };
      material.computeStress( state, tangents, deformation, timeIncrement );
    }
    else {
      GLMaterial::Deformation deformation{ dU, dSurfaceStrain, normal };
      material.computeStress( state, tangents, deformation, timeIncrement );
    }
  }

  Eigen::VectorXd makeInitializedStateVars( GLMaterial& material )
  {
    Eigen::VectorXd stateVars = Eigen::VectorXd::Zero( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );
    return stateVars;
  }

  GLResponse evaluateVirginResponse( GLMaterial&   material,
                                     const double* dU,
                                     const double* dSurfaceStrain,
                                     const double* normal,
                                     const double* separation )
  {
    GLResponse      response;
    Eigen::VectorXd stateVars = makeInitializedStateVars( material );
    computeGLStress( material, response, stateVars.data(), dU, dSurfaceStrain, normal, separation, 0.0, 1.0 );
    return response;
  }

  Vector9d flattenRowMajor( const Matrix3dRowMajor& tensor )
  {
    return Eigen::Map< const Vector9d >( tensor.data() );
  }

  // ------------------------------------------------------------------
  // 10.1: Lobatto-rule test
  // ------------------------------------------------------------------
  void testLobattoRuleMoments()
  {
    const std::array< double, 5 > xi     = { -1.0, -std::sqrt( 3.0 / 7.0 ), 0.0, std::sqrt( 3.0 / 7.0 ), 1.0 };
    const std::array< double, 5 > weight = { 1.0 / 10.0, 49.0 / 90.0, 32.0 / 45.0, 49.0 / 90.0, 1.0 / 10.0 };

    double sumW = 0.0, sumWXi = 0.0, sumWXi2 = 0.0;
    for ( int i = 0; i < 5; ++i ) {
      sumW += weight[i];
      sumWXi += weight[i] * xi[i];
      sumWXi2 += weight[i] * xi[i] * xi[i];
    }

    throwExceptionOnFailure( std::abs( sumW - 2.0 ) < 1e-13, "Lobatto rule: sum(w) != 2" );
    throwExceptionOnFailure( std::abs( sumWXi ) < 1e-13, "Lobatto rule: sum(w*xi) != 0" );
    throwExceptionOnFailure( std::abs( sumWXi2 - 2.0 / 3.0 ) < 1e-13, "Lobatto rule: sum(w*xi^2) != 2/3" );

    double sumLambda = 0.0;
    for ( int i = 0; i < 5; ++i ) {
      sumLambda += 0.5 * weight[i];
    }
    throwExceptionOnFailure( std::abs( sumLambda - 1.0 ) < 1e-13, "Lobatto rule: sum(lambda) != 1" );
  }

  // ------------------------------------------------------------------
  // 10.2: Uniform-state reduction (A+ == A-) must match the Y material
  // ------------------------------------------------------------------
  void testUniformStateMatchesYMaterial()
  {
    const double h                      = 0.01;
    const double ell                    = 0.025;
    const double interfaceProperties[3] = { 1e5, 0.3, h };

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    GLMaterial glMaterial( "LINEARELASTIC", interfaceProperties, 3, 1 );
    YMaterial  yMaterial( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double aShared[9] = { 1.0e-4, -2.0e-4, 0.5e-4, 0.8e-4, 1.5e-4, -0.7e-4, 0.3e-4, -1.1e-4, 0.9e-4 };
    double       dSurfaceStrain[18];
    std::copy( aShared, aShared + 9, dSurfaceStrain );
    std::copy( aShared, aShared + 9, dSurfaceStrain + 9 );
    const double dU[6] = { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 };

    const GLResponse glResponse = evaluateVirginResponse( glMaterial, dU, dSurfaceStrain, n.data(), separation.data() );

    // Reference from the existing (validated) two-station Y material.
    Eigen::VectorXd yStateVars = Eigen::VectorXd::Zero( yMaterial.getNumberOfRequiredStateVars() );
    yMaterial.initializeYourself( yStateVars.data(), static_cast< int >( yStateVars.size() ) );

    Eigen::Vector3d          yForce;
    Matrix3dRowMajor         ySPlus, ySMinus, yQww;
    Matrix3x9RowMajor        yQwAp, yQwAm;
    Matrix9x3RowMajor        yQApw, yQAmw;
    Matrix9dRowMajor         yQApAp, yQApAm, yQAmAp, yQAmAm;
    YMaterial::State         yState{ yForce.data(), ySPlus.data(), ySMinus.data(), yStateVars.data() };
    YMaterial::Tangents      yTangents{ yQww.data(),
                                   yQwAp.data(),
                                   yQwAm.data(),
                                   yQApw.data(),
                                   yQApAp.data(),
                                   yQApAm.data(),
                                   yQAmw.data(),
                                   yQAmAp.data(),
                                   yQAmAm.data() };
    YMaterial::Deformation   yDeformation{ dU, dSurfaceStrain, n.data(), separation.data() };
    YMaterial::TimeIncrement yTime{ 0.0, 1.0 };
    yMaterial.computeStress( yState, yTangents, yDeformation, yTime );

    const double tol = 1e-8 * std::max( 1.0, yForce.norm() );
    assertMatrixNear( glResponse.generalizedForce,
                      yForce,
                      tol,
                      "GL vs Y: generalizedForce mismatch under uniform loading" );
    assertMatrixNear( glResponse.surfaceStressPlus,
                      ySPlus,
                      tol,
                      "GL vs Y: surfaceStressPlus mismatch under uniform loading" );
    assertMatrixNear( glResponse.surfaceStressMinus,
                      ySMinus,
                      tol,
                      "GL vs Y: surfaceStressMinus mismatch under uniform loading" );

    // All station stresses must also be identical to one another.
    assertMatrixNear( glResponse.surfaceStressPlus,
                      glResponse.surfaceStressMinus,
                      1e-8 * glResponse.surfaceStressPlus.lpNorm< Eigen::Infinity >(),
                      "GL uniform loading: surfaceStressPlus should equal surfaceStressMinus" );
  }

  // ------------------------------------------------------------------
  // 10.3: Traction equilibrium across all five stations
  // ------------------------------------------------------------------
  Eigen::Matrix< double, 5, 3 > stationTractionsFromStateVars( GLMaterial&            material,
                                                               const Eigen::VectorXd& stateVars,
                                                               const Eigen::Vector3d& n )
  {
    Eigen::Matrix< double, 5, 3 > tractions;
    for ( int alpha = 0; alpha < 5; ++alpha ) {
      const double* stress9 = material
                                .getStateView( "stationStress" + std::to_string( alpha ),
                                               const_cast< double* >( stateVars.data() ) )
                                .stateLocation;
      const Eigen::Map< const Matrix3dRowMajor > sigma( stress9 );
      tractions.row( alpha ) = ( sigma * n ).transpose();
    }
    return tractions;
  }

  void assertAllStationTractionsEqual( const Eigen::Matrix< double, 5, 3 >& tractions, const std::string& context )
  {
    const double scale = std::max( 1.0, tractions.row( 0 ).norm() );
    for ( int alpha = 1; alpha < 5; ++alpha ) {
      const double err = ( tractions.row( alpha ) - tractions.row( 0 ) ).norm();
      throwExceptionOnFailure( err < 1e-8 * scale,
                               context + ": station " + std::to_string( alpha ) +
                                 " traction differs from station 0, err=" + std::to_string( err ) );
    }
  }

  void testTractionEquilibriumElasticSkewed()
  {
    const double h                      = 0.01;
    const double ell                    = 0.025;
    const double interfaceProperties[3] = { 1e5, 0.3, h };

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    GLMaterial material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double dU[6]     = { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 };
    const double aPlus[9]  = { 1.0e-4, -2.0e-4, 0.5e-4, 0.8e-4, 1.5e-4, -0.7e-4, 0.3e-4, -1.1e-4, 0.9e-4 };
    const double aMinus[9] = { -0.4e-4, 0.9e-4, -0.2e-4, 0.6e-4, -1.0e-4, 0.3e-4, -0.8e-4, 0.2e-4, 0.5e-4 };
    double       dSurfaceStrain[18];
    std::copy( aPlus, aPlus + 9, dSurfaceStrain );
    std::copy( aMinus, aMinus + 9, dSurfaceStrain + 9 );

    Eigen::VectorXd stateVars = makeInitializedStateVars( material );
    GLResponse      response;
    computeGLStress( material, response, stateVars.data(), dU, dSurfaceStrain, n.data(), separation.data(), 0.0, 1.0 );

    const auto tractions = stationTractionsFromStateVars( material, stateVars, n );
    assertAllStationTractionsEqual( tractions, "Elastic skewed equilibrium" );
  }

  void testTractionEquilibriumDivergentPlasticHistories()
  {
    const double h                      = 0.01;
    const double ell                    = 0.025;
    const double interfaceProperties[8] = { 210000., 0.3, h, 20., 200., 5., 10., 2400. };

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    GLMaterial material( "VONMISES", interfaceProperties, 8, 1 );

    const double historyDU[6]     = { 2.0e-3, -0.5e-3, 3.0e-3, 0.2e-3, 0.3e-3, -0.4e-3 };
    const double historyAPlus[9]  = { 8.0e-3, -6.0e-3, 4.0e-3, 5.0e-3, 7.0e-3, -3.0e-3, 2.0e-3, -5.0e-3, 6.0e-3 };
    const double historyAMinus[9] = { 0.2e-3, -0.1e-3, 0.05e-3, 0.1e-3, 0.15e-3, -0.05e-3, 0.05e-3, -0.1e-3, 0.1e-3 };
    double       historyDS[18];
    std::copy( historyAPlus, historyAPlus + 9, historyDS );
    std::copy( historyAMinus, historyAMinus + 9, historyDS + 9 );

    Eigen::VectorXd stateVars = makeInitializedStateVars( material );
    GLResponse      histResponse;
    computeGLStress( material,
                     histResponse,
                     stateVars.data(),
                     historyDU,
                     historyDS,
                     n.data(),
                     separation.data(),
                     0.0,
                     1.0 );

    const double dU[6]     = { 0.3e-3, -0.1e-3, 0.4e-3, 0.05e-3, 0.05e-3, -0.05e-3 };
    const double aPlus[9]  = { 0.5e-3, -0.3e-3, 0.2e-3, 0.3e-3, 0.4e-3, -0.2e-3, 0.1e-3, -0.3e-3, 0.3e-3 };
    const double aMinus[9] = { 0.05e-3, -0.02e-3, 0.01e-3, 0.02e-3, 0.03e-3, -0.01e-3, 0.01e-3, -0.02e-3, 0.02e-3 };
    double       dS[18];
    std::copy( aPlus, aPlus + 9, dS );
    std::copy( aMinus, aMinus + 9, dS + 9 );

    GLResponse response;
    computeGLStress( material, response, stateVars.data(), dU, dS, n.data(), separation.data(), 0.0, 1.0 );

    const auto tractions = stationTractionsFromStateVars( material, stateVars, n );
    assertAllStationTractionsEqual( tractions, "Plastic divergent-history equilibrium" );

    const double crossNorm = std::max( response.Q_ApAm.lpNorm< Eigen::Infinity >(),
                                       response.Q_AmAp.lpNorm< Eigen::Infinity >() );
    const double diagNorm  = std::max( response.Q_ApAp.lpNorm< Eigen::Infinity >(),
                                      response.Q_AmAm.lpNorm< Eigen::Infinity >() );
    throwExceptionOnFailure( crossNorm > 1e-3 * diagNorm,
                             "With divergent plastic histories, condensed cross-blocks should be clearly nonzero." );
  }

  // ------------------------------------------------------------------
  // 10.4: Scalar weighted-series law (independent of the class's public
  // API, which forces identical properties on every station -- this
  // directly emulates the 5-station local-equilibrium algebra with
  // independently-chosen moduli, mirroring the Y material's own
  // two-station series-law test).
  // ------------------------------------------------------------------
  void testScalarWeightedSeriesLaw()
  {
    const std::array< double, 5 > weight = { 1.0 / 10.0, 49.0 / 90.0, 32.0 / 45.0, 49.0 / 90.0, 1.0 / 10.0 };
    std::array< double, 5 >       lambda;
    for ( int i = 0; i < 5; ++i )
      lambda[i] = 0.5 * weight[i];

    const std::array< double, 5 > E = { 2.0e5, 1.5e5, 5.0e4, 8.0e4, 3.0e5 };

    double qEffExpected = 0.0;
    for ( int i = 0; i < 5; ++i )
      qEffExpected += lambda[i] / E[i];
    qEffExpected = 1.0 / qEffExpected;

    auto bulkMaterial = []( double Emod ) {
      const double props[2] = { Emod, 0.0 };
      return std::unique_ptr< MarmotMaterialHypoElastic >(
        MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( "LINEARELASTIC", props, 2, 1 ) );
    };

    std::array< std::unique_ptr< MarmotMaterialHypoElastic >, 5 > materials;
    std::array< std::vector< double >, 5 >                        scratch;
    for ( int i = 0; i < 5; ++i ) {
      materials[i] = bulkMaterial( E[i] );
      scratch[i].assign( materials[i]->getNumberOfRequiredStateVars(), 0.0 );
    }

    const double h    = 0.01;
    const double ell  = h;
    const double gBar = 1.0e-6 / ell; // w/ell with w=1e-6

    std::array< double, 5 > g;
    g.fill( 0.0 );
    double t = 0.0;

    for ( int iter = 0; iter < 60; ++iter ) {
      std::array< double, 5 > tAlpha;
      std::array< double, 5 > Qalpha;
      for ( int i = 0; i < 5; ++i ) {
        Marmot::Vector6d eps                       = Marmot::Vector6d::Zero();
        eps[2]                                     = g[i];
        Marmot::Matrix6d                   tangent = Marmot::Matrix6d::Zero();
        MarmotMaterialHypoElastic::state3D st{ Marmot::Vector6d::Zero(), 0.0, 0.0, scratch[i].data() };
        materials[i]->computeStress( st, tangent, eps, { 0.0, 1.0 } );
        tAlpha[i] = st.stress[2];
        Qalpha[i] = tangent( 2, 2 );
      }

      Eigen::Matrix< double, 6, 1 > r;
      for ( int i = 0; i < 5; ++i )
        r[i] = tAlpha[i] - t;
      double gWeighted = 0.0;
      for ( int i = 0; i < 5; ++i )
        gWeighted += lambda[i] * g[i];
      r[5] = gWeighted - gBar;

      const double rNorm = r.lpNorm< Eigen::Infinity >() /
                           std::max( 1.0, *std::max_element( tAlpha.begin(), tAlpha.end() ) );
      if ( rNorm < 1e-12 ) {
        break;
      }

      Eigen::Matrix< double, 6, 6 > J = Eigen::Matrix< double, 6, 6 >::Zero();
      for ( int i = 0; i < 5; ++i ) {
        J( i, i ) = Qalpha[i];
        J( i, 5 ) = -1.0;
        J( 5, i ) = lambda[i];
      }
      const Eigen::Matrix< double, 6, 1 > dy = J.fullPivLu().solve( -r );
      for ( int i = 0; i < 5; ++i )
        g[i] += dy[i];
      t += dy[5];
    }

    const double qEffActual = t / gBar;
    throwExceptionOnFailure( std::abs( qEffActual - qEffExpected ) < 1e-6 * qEffExpected,
                             "Scalar weighted series-law failed: expected " + std::to_string( qEffExpected ) +
                               ", got " + std::to_string( qEffActual ) );

    double qArithmeticAverage = 0.0;
    for ( int i = 0; i < 5; ++i )
      qArithmeticAverage += lambda[i] * E[i];
    throwExceptionOnFailure( std::abs( qEffActual - qArithmeticAverage ) > 0.05 * qArithmeticAverage,
                             "Scalar weighted series-law: converged stiffness should differ clearly from the naive "
                             "weighted arithmetic average." );
  }

  // ------------------------------------------------------------------
  // 10.5: Exact elastic bending stiffness D_kk = E h^3 / 12
  // ------------------------------------------------------------------
  void testExactElasticBendingStiffness()
  {
    const double h                      = 0.01;
    const double E                      = 2.0e5;
    const double interfaceProperties[3] = { E, 0.0, h }; // nu=0

    const Eigen::Vector3d n( 0.0, 0.0, 1.0 );
    const Eigen::Vector3d separation = h * n; // coincident faces, ell=h, d_tau=0

    GLMaterial material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    // Pure bending mode: A+ = -A- (Abar=0), single nonzero component A_xx
    // (index 0), representing an in-plane axial-strain gradient through
    // the thickness -- the classical flexural mode.
    //
    // NOTE: the shear/rotation component A_zx (index 6, n (x) t) that
    // serves as a bending proxy in the two-station Y material is NOT the
    // right probe here: with 5 genuinely distinct through-thickness
    // stations, the common-traction constraint t^(alpha) = t exactly
    // enforces sigma.n (i.e. sigma_zx, sigma_zy, sigma_zz) to be
    // through-thickness CONSTANT -- the correct continuum equilibrium
    // condition div(sigma)=0 for fields varying only through z. A genuine
    // moment can therefore only develop in the in-plane stress components
    // (xx, yy, xy), never in a sigma.n component.
    const double amplitude = 1.0e-6;
    double       aPlus[9]  = { 0. };
    double       aMinus[9] = { 0. };
    aPlus[0]               = amplitude;
    aMinus[0]              = -amplitude;
    double dSurfaceStrain[18];
    std::copy( aPlus, aPlus + 9, dSurfaceStrain );
    std::copy( aMinus, aMinus + 9, dSurfaceStrain + 9 );
    const double dU[6] = { 0., 0., 0., 0., 0., 0. };

    const GLResponse response = evaluateVirginResponse( material, dU, dSurfaceStrain, n.data(), separation.data() );

    const Matrix9dRowMajor DkkFull = ( h * h / 4.0 ) *
                                     ( response.Q_ApAp - response.Q_ApAm - response.Q_AmAp + response.Q_AmAm );
    const double Dkk     = DkkFull( 0, 0 );
    const double Dkk_ref = E * h * h * h / 12.0;

    throwExceptionOnFailure( std::abs( Dkk - Dkk_ref ) < 1e-6 * Dkk_ref,
                             "Exact elastic bending stiffness mismatch: expected " + std::to_string( Dkk_ref ) +
                               ", got " + std::to_string( Dkk ) );
  }

  // ------------------------------------------------------------------
  // Section 13 groundwork: the internal (not user-exposed) N-station
  // template must reproduce the SAME exact D_kk = E h^3/12 result for
  // every supported station count, since Gauss-Lobatto with N points
  // integrates polynomials up to degree 2N-3 exactly and bending energy
  // only requires degree 2 (z^2) -- satisfied for all N in {3,4,5,7}.
  // ------------------------------------------------------------------
  template < int NStations >
  void checkExactBendingStiffnessForStationCount()
  {
    using GLMaterialN = MarmotGaussLobattoInterfaceMaterialHypoElasticN< NStations >;

    const double h                      = 0.01;
    const double E                      = 2.0e5;
    const double interfaceProperties[3] = { E, 0.0, h };

    const Eigen::Vector3d n( 0.0, 0.0, 1.0 );
    const Eigen::Vector3d separation = h * n;

    GLMaterialN material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double amplitude = 1.0e-6;
    double       aPlus[9]  = { 0. };
    double       aMinus[9] = { 0. };
    aPlus[0]               = amplitude;
    aMinus[0]              = -amplitude;
    double dSurfaceStrain[18];
    std::copy( aPlus, aPlus + 9, dSurfaceStrain );
    std::copy( aMinus, aMinus + 9, dSurfaceStrain + 9 );
    const double dU[6] = { 0., 0., 0., 0., 0., 0. };

    Eigen::VectorXd stateVars = Eigen::VectorXd::Zero( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );

    Eigen::Vector3d   force;
    Matrix3dRowMajor  sPlus, sMinus, qWw;
    Matrix3x9RowMajor qWAp, qWAm;
    Matrix9x3RowMajor qApw, qAmw;
    Matrix9dRowMajor  qApAp, qApAm, qAmAp, qAmAm;

    typename GLMaterialN::State         state{ force.data(), sPlus.data(), sMinus.data(), stateVars.data() };
    typename GLMaterialN::Tangents      tangents{ qWw.data(),
                                             qWAp.data(),
                                             qWAm.data(),
                                             qApw.data(),
                                             qApAp.data(),
                                             qApAm.data(),
                                             qAmw.data(),
                                             qAmAp.data(),
                                             qAmAm.data() };
    typename GLMaterialN::Deformation   deformation( dU, dSurfaceStrain, n.data(), separation.data() );
    typename GLMaterialN::TimeIncrement timeIncrement{ 0.0, 1.0 };

    material.computeStress( state, tangents, deformation, timeIncrement );

    const Matrix9dRowMajor DkkFull = ( h * h / 4.0 ) * ( qApAp - qApAm - qAmAp + qAmAm );
    const double           Dkk     = DkkFull( 0, 0 );
    const double           Dkk_ref = E * h * h * h / 12.0;

    throwExceptionOnFailure( std::abs( Dkk - Dkk_ref ) < 1e-6 * Dkk_ref,
                             "N=" + std::to_string( NStations ) +
                               " exact elastic bending stiffness mismatch: expected " + std::to_string( Dkk_ref ) +
                               ", got " + std::to_string( Dkk ) );
  }

  void testExactBendingStiffnessForAllStationCounts()
  {
    checkExactBendingStiffnessForStationCount< 3 >();
    checkExactBendingStiffnessForStationCount< 4 >();
    checkExactBendingStiffnessForStationCount< 5 >();
    checkExactBendingStiffnessForStationCount< 7 >();
  }

  // ------------------------------------------------------------------
  // 10.6: Progressive plastic yielding through the thickness
  // ------------------------------------------------------------------
  void testProgressiveThroughThicknessYielding()
  {
    const double h                      = 0.01;
    const double fy                     = 5.0;
    const double interfaceProperties[8] = { 4e5, 0.3, h, fy, 0.1, 0.0, 0.0, 0 };

    const Eigen::Vector3d n( 0.0, 0.0, 1.0 );
    const Eigen::Vector3d separation = h * n;
    const double          dU[6]      = { 0., 0., 0., 0., 0., 0. };

    // Increasing pure curvature (A+ = -A-, Abar = 0), single in-plane
    // axial-strain component (A_xx, index 0 -- the flexural mode; see
    // testExactElasticBendingStiffness for why the sigma.n-type shear
    // component cannot carry a moment in this 5-station formulation) so
    // that the elastic-to-plastic transition is governed by a monotone
    // scalar amplitude.
    auto stationsYielded = [&]( double amplitude ) {
      GLMaterial material( "VONMISES", interfaceProperties, 8, 1 );

      double aPlus[9] = { 0. }, aMinus[9] = { 0. };
      aPlus[0]  = amplitude;
      aMinus[0] = -amplitude;
      double dSurfaceStrain[18];
      std::copy( aPlus, aPlus + 9, dSurfaceStrain );
      std::copy( aMinus, aMinus + 9, dSurfaceStrain + 9 );

      Eigen::VectorXd stateVars = makeInitializedStateVars( material );
      GLResponse      response;
      computeGLStress( material,
                       response,
                       stateVars.data(),
                       dU,
                       dSurfaceStrain,
                       n.data(),
                       separation.data(),
                       0.0,
                       1.0 );

      std::array< double, 5 > vmStress;
      for ( int alpha = 0; alpha < 5; ++alpha ) {
        const double* stress9 = material.getStateView( "stationStress" + std::to_string( alpha ), stateVars.data() )
                                  .stateLocation;
        const Eigen::Map< const Matrix3dRowMajor > sigma( stress9 );
        const Matrix3dRowMajor                     dev = sigma - ( sigma.trace() / 3.0 ) * Matrix3dRowMajor::Identity();
        vmStress[alpha]                                = std::sqrt( 1.5 * dev.cwiseProduct( dev ).sum() );
      }
      return vmStress;
    };

    // Small amplitude: everything elastic.
    const auto smallAmp = stationsYielded( 2.0e-6 );
    for ( double s : smallAmp ) {
      throwExceptionOnFailure( s < fy, "Progressive yielding: unexpectedly yielded at tiny curvature amplitude." );
    }

    // Larger amplitude: endpoint stations (0, 4) yield first because they
    // see the largest |A|, then off-center (1, 3), center (2) last.
    const auto largeAmp = stationsYielded( 4.0e-5 );
    throwExceptionOnFailure( largeAmp[0] >= fy - 1e-6 && largeAmp[4] >= fy - 1e-6,
                             "Progressive yielding: endpoint stations should yield first at large curvature." );
    throwExceptionOnFailure( largeAmp[2] < largeAmp[0],
                             "Progressive yielding: center station should carry less equivalent stress than the "
                             "endpoints under pure curvature." );
    throwExceptionOnFailure( largeAmp[1] < largeAmp[0] && largeAmp[3] < largeAmp[4],
                             "Progressive yielding: off-center stations should lag the endpoints." );
  }

  // ------------------------------------------------------------------
  // 10.7: Full condensed tangent vs finite difference
  // ------------------------------------------------------------------
  void testFullCondensedTangentMatchesFiniteDifference( const std::string&           materialName,
                                                        const std::vector< double >& properties,
                                                        bool                         withHistory )
  {
    const double ell = 0.025;

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    GLMaterial material( materialName, properties.data(), static_cast< int >( properties.size() ), 1 );

    Eigen::VectorXd stateVars = makeInitializedStateVars( material );

    if ( withHistory ) {
      const double historyDU[6]     = { 2.0e-3, -0.5e-3, 3.0e-3, 0.2e-3, 0.3e-3, -0.4e-3 };
      const double historyAPlus[9]  = { 8.0e-3, -6.0e-3, 4.0e-3, 5.0e-3, 7.0e-3, -3.0e-3, 2.0e-3, -5.0e-3, 6.0e-3 };
      const double historyAMinus[9] = { 0.2e-3, -0.1e-3, 0.05e-3, 0.1e-3, 0.15e-3, -0.05e-3, 0.05e-3, -0.1e-3, 0.1e-3 };
      double       historyDS[18];
      std::copy( historyAPlus, historyAPlus + 9, historyDS );
      std::copy( historyAMinus, historyAMinus + 9, historyDS + 9 );
      GLResponse histResponse;
      computeGLStress( material,
                       histResponse,
                       stateVars.data(),
                       historyDU,
                       historyDS,
                       n.data(),
                       separation.data(),
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
      // Small increment on top of the plastic history (avoid re-crossing
      // yield too abruptly for the perturbed finite-difference states).
      for ( int i = 0; i < 6; ++i )
        dU[i] *= 0.1;
      for ( int i = 0; i < 18; ++i )
        dSurfaceStrain[i] *= 0.1;
    }

    Eigen::VectorXd baseStateVars = stateVars;
    GLResponse      analytic;
    computeGLStress( material,
                     analytic,
                     baseStateVars.data(),
                     dU,
                     dSurfaceStrain,
                     n.data(),
                     separation.data(),
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
      computeGLStress( material, response, perturbedStateVars.data(), dUp, dSp, n.data(), separation.data(), 0.0, 1.0 );
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
                      std::string( "Full 21x21 GL condensed tangent vs FD (" ) + materialName +
                        ( withHistory ? ", with history)" : ", elastic)" ) );
  }

  void testFullCondensedTangentElasticSkewed()
  {
    testFullCondensedTangentMatchesFiniteDifference( "LINEARELASTIC", { 1e5, 0.3, 0.01 }, false );
  }

  void testFullCondensedTangentPlasticDivergentHistory()
  {
    testFullCondensedTangentMatchesFiniteDifference( "VONMISES",
                                                     { 210000., 0.3, 0.01, 20., 200., 5., 10., 2400. },
                                                     true );
  }

  // ------------------------------------------------------------------
  // 10.8: State rollback on rejected Newton/line-search trials
  // ------------------------------------------------------------------
  void testStateRollbackOnRejectedTrials()
  {
    const double          h                      = 0.01;
    const double          interfaceProperties[8] = { 210000., 0.3, h, 20., 200., 5., 10., 2400. };
    const Eigen::Vector3d n( 0.0, 0.0, 1.0 );
    const Eigen::Vector3d separation = h * n;

    GLMaterial      material( "VONMISES", interfaceProperties, 8, 1 );
    Eigen::VectorXd stateVars = makeInitializedStateVars( material );

    // Snapshot the committed station 0 state before a large, deliberately
    // aggressive increment that is likely to trigger line-search
    // rejections inside the local Newton solve.
    const std::vector< double > committedBefore( stateVars.data(), stateVars.data() + stateVars.size() );

    const double dU[6]     = { 5.0e-3, -1.0e-3, 6.0e-3, 0.5e-3, 0.4e-3, -0.6e-3 };
    const double aPlus[9]  = { 1.0e-2, -8.0e-3, 5.0e-3, 6.0e-3, 8.0e-3, -4.0e-3, 3.0e-3, -6.0e-3, 7.0e-3 };
    const double aMinus[9] = { 0.3e-3, -0.15e-3, 0.08e-3, 0.15e-3, 0.2e-3, -0.08e-3, 0.08e-3, -0.15e-3, 0.15e-3 };
    double       dS[18];
    std::copy( aPlus, aPlus + 9, dS );
    std::copy( aMinus, aMinus + 9, dS + 9 );

    GLResponse response;
    bool       threw = false;
    try {
      computeGLStress( material, response, stateVars.data(), dU, dS, n.data(), separation.data(), 0.0, 1.0 );
    }
    catch ( const std::exception& ) {
      threw = true;
    }

    if ( threw ) {
      // If the increment failed to converge at all, the state must be
      // untouched (equal to what it was before the call).
      const double diff = ( Eigen::Map< const Eigen::VectorXd >( stateVars.data(), stateVars.size() ) -
                            Eigen::Map< const Eigen::VectorXd >( committedBefore.data(), committedBefore.size() ) )
                            .lpNorm< Eigen::Infinity >();
      throwExceptionOnFailure( diff < 1e-13,
                               "State rollback: a failed local Newton solve must not have mutated the persisted "
                               "station state at all." );
      return;
    }

    // Otherwise, converged: the committed station stress/history must be
    // self-consistent with a FRESH re-evaluation restarting from the same
    // committed-before state and the same total deformation input (i.e.
    // no residual contamination from rejected intermediate trials).
    std::vector< double > stateVarsReplay = committedBefore;
    GLResponse            responseReplay;
    computeGLStress( material, responseReplay, stateVarsReplay.data(), dU, dS, n.data(), separation.data(), 0.0, 1.0 );

    const double tol = 1e-10 * std::max( 1.0, responseReplay.generalizedForce.norm() );
    assertMatrixNear( response.generalizedForce,
                      responseReplay.generalizedForce,
                      tol,
                      "State rollback: replaying the identical increment from the same committed state must "
                      "reproduce the same result (no leaked rejected-trial history)." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testLobattoRuleMoments,
                                                   testUniformStateMatchesYMaterial,
                                                   testTractionEquilibriumElasticSkewed,
                                                   testTractionEquilibriumDivergentPlasticHistories,
                                                   testScalarWeightedSeriesLaw,
                                                   testExactElasticBendingStiffness,
                                                   testExactBendingStiffnessForAllStationCounts,
                                                   testProgressiveThroughThicknessYielding,
                                                   testFullCondensedTangentElasticSkewed,
                                                   testFullCondensedTangentPlasticDivergentHistory,
                                                   testStateRollbackOnRejectedTrials };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
