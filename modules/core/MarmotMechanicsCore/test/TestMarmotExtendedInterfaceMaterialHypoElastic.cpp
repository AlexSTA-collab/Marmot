#include "Marmot/MarmotExtendedInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
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

  using ExtendedMaterial = MarmotExtendedInterfaceMaterialHypoElastic;

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
   * All generalized outputs of one extended-interface stress update.
   * force/surfaceStress double as the persistent generalized state between
   * increments, exactly as in the element.
   */
  struct ExtendedResponse {
    Eigen::Vector3d   force         = Eigen::Vector3d::Zero();
    Matrix3dRowMajor  surfaceStress = Matrix3dRowMajor::Zero();
    Matrix3dRowMajor  Q             = Matrix3dRowMajor::Zero();
    Matrix9dRowMajor  Z             = Matrix9dRowMajor::Zero();
    Matrix3x9RowMajor H             = Matrix3x9RowMajor::Zero();
    Matrix9x3RowMajor K             = Matrix9x3RowMajor::Zero();
  };

  std::unique_ptr< MarmotMaterialHypoElastic > createBulkMaterial( const std::string& materialName,
                                                                   const double*      properties,
                                                                   int                nProperties )
  {
    return std::unique_ptr< MarmotMaterialHypoElastic >(
      MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName, properties, nProperties, 1 ) );
  }

  /**
   * Run one extended-interface stress update. A null @p separation exercises
   * the backward-compatible 3-argument Deformation constructor (coincident
   * faces); otherwise the full separation-vector-aware path is used.
   */
  void computeExtendedStress( ExtendedMaterial& material,
                              ExtendedResponse& response,
                              double*           stateVars,
                              const double*     dU,
                              const double*     dSurfaceStrain,
                              const double*     normal,
                              const double*     separation,
                              double            timeOld,
                              double            dT )
  {
    ExtendedMaterial::State    state{ response.force.data(), response.surfaceStress.data(), stateVars };
    ExtendedMaterial::Tangents tangents{ response.Q.data(), response.Z.data(), response.H.data(), response.K.data() };
    ExtendedMaterial::TimeIncrement timeIncrement{ timeOld, dT };

    if ( separation ) {
      ExtendedMaterial::Deformation deformation{ dU, dSurfaceStrain, normal, separation };
      material.computeStress( state, tangents, deformation, timeIncrement );
    }
    else {
      ExtendedMaterial::Deformation deformation{ dU, dSurfaceStrain, normal };
      material.computeStress( state, tangents, deformation, timeIncrement );
    }
  }

  Eigen::VectorXd makeInitializedStateVars( ExtendedMaterial& material )
  {
    Eigen::VectorXd stateVars = Eigen::VectorXd::Zero( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );
    return stateVars;
  }

  /** Single virgin-state evaluation, used by the finite-difference checks. */
  ExtendedResponse evaluateVirginResponse( ExtendedMaterial& material,
                                           const double*     dU,
                                           const double*     dSurfaceStrain,
                                           const double*     normal,
                                           const double*     separation )
  {
    ExtendedResponse response;
    Eigen::VectorXd  stateVars = makeInitializedStateVars( material );
    computeExtendedStress( material, response, stateVars.data(), dU, dSurfaceStrain, normal, separation, 0.0, 1.0 );
    return response;
  }

  Vector9d flattenRowMajor( const Matrix3dRowMajor& tensor )
  {
    return Eigen::Map< const Vector9d >( tensor.data() );
  }

  /**
   * Zero (or omitted) separation vector must fall back to ell = h,
   * d_tau = 0: the generalized outputs then coincide with the plain
   * bulk response, force = sigma * n and surfaceStress = h * sigma.
   */
  void testZeroSeparationAgainstBulkMaterial( const std::string& materialName,
                                              const double*      interfaceProperties,
                                              int                nInterfaceProperties,
                                              const double*      bulkProperties,
                                              int                nBulkProperties )
  {
    const double h         = interfaceProperties[2];
    const double normal[3] = { 0., 0., 1. };

    auto interfaceMaterial = std::make_unique< ExtendedMaterial >( materialName,
                                                                   interfaceProperties,
                                                                   nInterfaceProperties,
                                                                   1 );
    auto bulkMaterial      = createBulkMaterial( materialName, bulkProperties, nBulkProperties );

    Eigen::VectorXd interfaceStateVars = makeInitializedStateVars( *interfaceMaterial );
    Eigen::VectorXd bulkStateVars( bulkMaterial->getNumberOfRequiredStateVars() );
    bulkMaterial->initializeYourself( bulkStateVars.data(), bulkStateVars.size() );

    ExtendedResponse response;
    Marmot::Vector6d bulkStress = Marmot::Vector6d::Zero();

    struct Increment {
      double dT;
      double jumpY;
      double surfaceShear;
    };
    const std::vector< Increment > increments = {
      { 0.01, 1e-4, 2e-4 },
      { 10.0, 0.0, 0.0 },
    };

    double timeOld = 0.0;
    for ( const auto& increment : increments ) {
      const double dU[6]              = { 0., increment.jumpY, 0., 0., 0., 0. };
      const double dSurfaceStrain[18] = { 0.,
                                          increment.surfaceShear,
                                          0.,
                                          increment.surfaceShear,
                                          0.,
                                          0.,
                                          0.,
                                          0.,
                                          0.,
                                          0.,
                                          increment.surfaceShear,
                                          0.,
                                          increment.surfaceShear,
                                          0.,
                                          0.,
                                          0.,
                                          0.,
                                          0. };

      computeExtendedStress( *interfaceMaterial,
                             response,
                             interfaceStateVars.data(),
                             dU,
                             dSurfaceStrain,
                             normal,
                             nullptr,
                             timeOld,
                             increment.dT );

      Marmot::Vector6d bulkStrainIncrement = Marmot::Vector6d::Zero();
      bulkStrainIncrement[3]               = 2. * increment.surfaceShear;
      bulkStrainIncrement[5]               = increment.jumpY / h;

      Marmot::Matrix6d                   bulkTangent = Marmot::Matrix6d::Zero();
      MarmotMaterialHypoElastic::state3D bulkState{ bulkStress, 0.0, 0.0, bulkStateVars.data() };
      bulkMaterial->computeStress( bulkState, bulkTangent, bulkStrainIncrement, { timeOld, increment.dT } );
      bulkStress = bulkState.stress;

      const Eigen::Matrix3d expectedStress = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( bulkStress );
      const Eigen::Vector3d expectedForce  = expectedStress * Eigen::Vector3d::UnitZ();

      assertMatrixNear( response.force,
                        expectedForce,
                        1e-10,
                        materialName + ": zero-separation interface force does not match bulk stress" );
      assertMatrixNear( response.surfaceStress,
                        ( h * expectedStress ).eval(),
                        1e-10,
                        materialName + ": zero-separation surface resultant does not match h * sigma" );
      throwExceptionOnFailure( checkIfEqual< double >( interfaceStateVars, bulkStateVars, 1e-10 ),
                               materialName + ": interface state variables do not match bulk state variables." );

      timeOld += increment.dT;
    }
  }

  void testZeroSeparationLinearElasticAgainstBulk()
  {
    const double interfaceProperties[3] = { 1e5, 0.3, 0.01 };
    const double bulkProperties[2]      = { 1e5, 0.3 };
    testZeroSeparationAgainstBulkMaterial( "LINEARELASTIC", interfaceProperties, 3, bulkProperties, 2 );
  }

  void testZeroSeparationVonMisesAgainstBulk()
  {
    const double interfaceProperties[8] = { 1e5, 0.3, 0.01, 100., 10., 0., 1., 2400. };
    const double bulkProperties[7]      = { 1e5, 0.3, 100., 10., 0., 1., 2400. };
    testZeroSeparationAgainstBulkMaterial( "VONMISES", interfaceProperties, 8, bulkProperties, 7 );
  }

  /**
   * The 3-argument Deformation constructor (coincident faces) and an explicit
   * zero separation vector must take the identical fallback branch.
   */
  void testExplicitZeroSeparationMatchesCoincidentFacePath()
  {
    const double     interfaceProperties[3] = { 1e5, 0.3, 0.01 };
    ExtendedMaterial material( "LINEARELASTIC", interfaceProperties, 3, 1 );

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
                                        -0.6e-4,
                                        1.2e-4,
                                        0.4e-4,
                                        -0.9e-4,
                                        0.7e-4,
                                        1.3e-4,
                                        -0.2e-4,
                                        0.5e-4,
                                        -1.4e-4 };
    const double normal[3]          = { 0., 0., 1. };
    const double zeroSeparation[3]  = { 0., 0., 0. };

    const ExtendedResponse coincident   = evaluateVirginResponse( material, dU, dSurfaceStrain, normal, nullptr );
    const ExtendedResponse explicitZero = evaluateVirginResponse( material,
                                                                  dU,
                                                                  dSurfaceStrain,
                                                                  normal,
                                                                  zeroSeparation );

    assertMatrixNear( explicitZero.force, coincident.force, 1e-14, "Explicit zero separation: force differs" );
    assertMatrixNear( explicitZero.surfaceStress,
                      coincident.surfaceStress,
                      1e-14,
                      "Explicit zero separation: surface resultant differs" );
    assertMatrixNear( explicitZero.Q, coincident.Q, 1e-14, "Explicit zero separation: Q differs" );
    assertMatrixNear( explicitZero.Z, coincident.Z, 1e-14, "Explicit zero separation: Z differs" );
    assertMatrixNear( explicitZero.H, coincident.H, 1e-14, "Explicit zero separation: H differs" );
    assertMatrixNear( explicitZero.K, coincident.K, 1e-14, "Explicit zero separation: K differs" );
  }

  /**
   * For coincident faces the extended material must reproduce the plain
   * interface material's generalized stress state for identical inputs.
   * (The tangent operators intentionally differ: the plain material uses a
   * condensed formulation, the extended one the full-gradient formulation.)
   */
  void testZeroSeparationMatchesPlainInterfaceMaterialStress()
  {
    const double interfaceProperties[3] = { 1e5, 0.3, 0.01 };

    ExtendedMaterial                   extendedMaterial( "LINEARELASTIC", interfaceProperties, 3, 1 );
    MarmotInterfaceMaterialHypoElastic plainMaterial( "LINEARELASTIC", interfaceProperties, 3, 1 );

    Eigen::VectorXd extendedStateVars = makeInitializedStateVars( extendedMaterial );
    Eigen::VectorXd plainStateVars    = Eigen::VectorXd::Zero( plainMaterial.getNumberOfRequiredStateVars() );
    plainMaterial.initializeYourself( plainStateVars.data(), plainStateVars.size() );

    ExtendedResponse extendedResponse;

    Eigen::Vector3d  plainForce         = Eigen::Vector3d::Zero();
    Matrix3dRowMajor plainSurfaceStress = Matrix3dRowMajor::Zero();

    const double normal[3] = { 0., 0., 1. };

    const double dUIncrements[2][6] = { { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 },
                                        { -0.5e-4, 0.8e-4, 1.0e-4, 0.2e-4, -0.3e-4, 0.4e-4 } };

    const double dSurfaceStrainIncrements[2][18] = { { 1.0e-4,
                                                       -2.0e-4,
                                                       0.5e-4,
                                                       0.8e-4,
                                                       1.5e-4,
                                                       -0.7e-4,
                                                       0.3e-4,
                                                       -1.1e-4,
                                                       0.9e-4,
                                                       -0.6e-4,
                                                       1.2e-4,
                                                       0.4e-4,
                                                       -0.9e-4,
                                                       0.7e-4,
                                                       1.3e-4,
                                                       -0.2e-4,
                                                       0.5e-4,
                                                       -1.4e-4 },
                                                     { -0.4e-4,
                                                       0.9e-4,
                                                       -0.2e-4,
                                                       0.6e-4,
                                                       -1.0e-4,
                                                       0.3e-4,
                                                       -0.8e-4,
                                                       0.2e-4,
                                                       0.5e-4,
                                                       0.7e-4,
                                                       -0.3e-4,
                                                       0.8e-4,
                                                       -0.5e-4,
                                                       0.4e-4,
                                                       -0.6e-4,
                                                       0.1e-4,
                                                       -0.9e-4,
                                                       0.2e-4 } };

    double timeOld = 0.0;
    for ( int increment = 0; increment < 2; ++increment ) {
      computeExtendedStress( extendedMaterial,
                             extendedResponse,
                             extendedStateVars.data(),
                             dUIncrements[increment],
                             dSurfaceStrainIncrements[increment],
                             normal,
                             nullptr,
                             timeOld,
                             1.0 );

      double Q[9]  = { 0. };
      double Z[81] = { 0. };
      double H[27] = { 0. };
      double Y[81] = { 0. };

      MarmotInterfaceMaterialHypoElastic::State         plainState{ plainForce.data(),
                                                            plainSurfaceStress.data(),
                                                            plainStateVars.data() };
      MarmotInterfaceMaterialHypoElastic::Tangents      plainTangents{ Q, Z, H, Y };
      MarmotInterfaceMaterialHypoElastic::Deformation   plainDeformation{ dUIncrements[increment],
                                                                        dSurfaceStrainIncrements[increment],
                                                                        normal };
      MarmotInterfaceMaterialHypoElastic::TimeIncrement plainTimeIncrement{ timeOld, 1.0 };
      plainMaterial.computeStress( plainState, plainTangents, plainDeformation, plainTimeIncrement );

      assertMatrixNear( extendedResponse.force,
                        plainForce,
                        1e-12,
                        "Zero-separation extended force differs from plain interface material" );
      assertMatrixNear( extendedResponse.surfaceStress,
                        plainSurfaceStress,
                        1e-12,
                        "Zero-separation extended surface resultant differs from plain interface material" );
      throwExceptionOnFailure( checkIfEqual< double >( extendedStateVars, plainStateVars, 1e-12 ),
                               "Zero-separation extended state variables differ from plain interface material." );

      timeOld += 1.0;
    }
  }

  /**
   * Nonzero separation vector with ell != h and a nonzero tangential
   * component d_tau: the reconstructed geometry
   *
   *   G = A + (1/ell) ( [u] - A d_tau ) \otimes n,
   *   force = (h/ell) sigma n,
   *   surfaceStress = h sigma - force \otimes d_tau,
   *
   * must feed through exactly, verified against a directly-driven bulk material.
   */
  void testNonzeroSeparationReconstructedGeometryMatchesReference()
  {
    const double h                      = 0.01;
    const double ell                    = 0.025;
    const double interfaceProperties[3] = { 1e5, 0.3, h };
    const double bulkProperties[2]      = { 1e5, 0.3 };

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    ExtendedMaterial interfaceMaterial( "LINEARELASTIC", interfaceProperties, 3, 1 );
    auto             bulkMaterial = createBulkMaterial( "LINEARELASTIC", bulkProperties, 2 );

    Eigen::VectorXd interfaceStateVars = makeInitializedStateVars( interfaceMaterial );
    Eigen::VectorXd bulkStateVars      = Eigen::VectorXd::Zero( bulkMaterial->getNumberOfRequiredStateVars() );
    bulkMaterial->initializeYourself( bulkStateVars.data(), bulkStateVars.size() );

    ExtendedResponse response;
    Marmot::Vector6d bulkStress = Marmot::Vector6d::Zero();

    const double dUIncrements[2][6] = { { 1.2e-4, -0.4e-4, 2.0e-4, 0.3e-4, 0.5e-4, -0.6e-4 },
                                        { -0.5e-4, 0.8e-4, 1.0e-4, 0.2e-4, -0.3e-4, 0.4e-4 } };

    const double dSurfaceStrainIncrements[2][18] = { { 1.0e-4,
                                                       -2.0e-4,
                                                       0.5e-4,
                                                       0.8e-4,
                                                       1.5e-4,
                                                       -0.7e-4,
                                                       0.3e-4,
                                                       -1.1e-4,
                                                       0.9e-4,
                                                       -0.6e-4,
                                                       1.2e-4,
                                                       0.4e-4,
                                                       -0.9e-4,
                                                       0.7e-4,
                                                       1.3e-4,
                                                       -0.2e-4,
                                                       0.5e-4,
                                                       -1.4e-4 },
                                                     { -0.4e-4,
                                                       0.9e-4,
                                                       -0.2e-4,
                                                       0.6e-4,
                                                       -1.0e-4,
                                                       0.3e-4,
                                                       -0.8e-4,
                                                       0.2e-4,
                                                       0.5e-4,
                                                       0.7e-4,
                                                       -0.3e-4,
                                                       0.8e-4,
                                                       -0.5e-4,
                                                       0.4e-4,
                                                       -0.6e-4,
                                                       0.1e-4,
                                                       -0.9e-4,
                                                       0.2e-4 } };

    double timeOld = 0.0;
    for ( int increment = 0; increment < 2; ++increment ) {
      computeExtendedStress( interfaceMaterial,
                             response,
                             interfaceStateVars.data(),
                             dUIncrements[increment],
                             dSurfaceStrainIncrements[increment],
                             n.data(),
                             separation.data(),
                             timeOld,
                             1.0 );

      // Reference: reconstruct the full displacement-gradient increment.
      const Eigen::Map< const Eigen::Matrix< double, 6, 1 > > dU( dUIncrements[increment] );
      const Eigen::Vector3d                                   jump = dU.segment< 3 >( 0 ) - dU.segment< 3 >( 3 );

      const Eigen::Map< const Matrix3dRowMajor > topGradient( dSurfaceStrainIncrements[increment] );
      const Eigen::Map< const Matrix3dRowMajor > bottomGradient( dSurfaceStrainIncrements[increment] + 9 );
      const Matrix3dRowMajor                     A = 0.5 * ( topGradient + bottomGradient );

      const Matrix3dRowMajor G                    = A + ( ( jump - A * dTangential ) / ell ) * n.transpose();
      const Eigen::Matrix3d  strainIncrement      = 0.5 * ( G + G.transpose() );
      const Marmot::Vector6d strainIncrementVoigt = Marmot::ContinuumMechanics::VoigtNotation::strainToVoigt(
        strainIncrement );

      Marmot::Matrix6d                   bulkTangent = Marmot::Matrix6d::Zero();
      MarmotMaterialHypoElastic::state3D bulkState{ bulkStress, 0.0, 0.0, bulkStateVars.data() };
      bulkMaterial->computeStress( bulkState, bulkTangent, strainIncrementVoigt, { timeOld, 1.0 } );
      bulkStress = bulkState.stress;

      const Eigen::Matrix3d sigma = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( bulkStress );

      const Eigen::Vector3d expectedForce         = ( h / ell ) * sigma * n;
      const Eigen::Matrix3d expectedSurfaceStress = h * sigma - expectedForce * dTangential.transpose();

      assertMatrixNear( response.force,
                        expectedForce,
                        1e-10,
                        "Nonzero-separation force does not match (h/ell) sigma n" );
      assertMatrixNear( response.surfaceStress,
                        expectedSurfaceStress,
                        1e-10,
                        "Nonzero-separation surface resultant does not match h sigma - force x d_tau" );

      timeOld += 1.0;
    }
  }

  /**
   * Finite-difference verification of all four analytic tangent blocks for a
   * nonzero separation vector (ell != h, d_tau != 0):
   *
   *   Q_ik       = d force_i / d [u]_k,
   *   H_i(kl)    = d force_i / d <u_{k,l}>_s,
   *   K_(ij)k    = d surfaceStress_ij / d [u]_k,
   *   Z_(ij)(kl) = d surfaceStress_ij / d <u_{k,l}>_s.
   */
  void testNonzeroSeparationTangentsMatchFiniteDifferences()
  {
    const double h                      = 0.01;
    const double ell                    = 0.025;
    const double interfaceProperties[3] = { 1e5, 0.3, h };

    const Eigen::Vector3d n( 0.6, 0.0, 0.8 );
    const Eigen::Vector3d dTangential = 0.004 * Eigen::Vector3d( 0.8, 0.0, -0.6 ) +
                                        0.003 * Eigen::Vector3d( 0.0, 1.0, 0.0 );
    const Eigen::Vector3d separation = ell * n + dTangential;

    ExtendedMaterial material( "LINEARELASTIC", interfaceProperties, 3, 1 );

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
                                        -0.6e-4,
                                        1.2e-4,
                                        0.4e-4,
                                        -0.9e-4,
                                        0.7e-4,
                                        1.3e-4,
                                        -0.2e-4,
                                        0.5e-4,
                                        -1.4e-4 };

    const ExtendedResponse analytic = evaluateVirginResponse( material,
                                                              dU,
                                                              dSurfaceStrain,
                                                              n.data(),
                                                              separation.data() );

    const double eps = 1e-6;

    Matrix3dRowMajor  QFiniteDifference = Matrix3dRowMajor::Zero();
    Matrix9x3RowMajor KFiniteDifference = Matrix9x3RowMajor::Zero();

    // Perturbing the top displacement entry k changes the jump [u]_k by +/- eps.
    for ( int k = 0; k < 3; ++k ) {
      double dUPlus[6];
      double dUMinus[6];
      std::copy( dU, dU + 6, dUPlus );
      std::copy( dU, dU + 6, dUMinus );
      dUPlus[k] += eps;
      dUMinus[k] -= eps;

      const ExtendedResponse plus = evaluateVirginResponse( material,
                                                            dUPlus,
                                                            dSurfaceStrain,
                                                            n.data(),
                                                            separation.data() );

      const ExtendedResponse minus = evaluateVirginResponse( material,
                                                             dUMinus,
                                                             dSurfaceStrain,
                                                             n.data(),
                                                             separation.data() );

      QFiniteDifference.col( k ) = ( plus.force - minus.force ) / ( 2.0 * eps );
      KFiniteDifference.col( k ) = ( flattenRowMajor( plus.surfaceStress ) - flattenRowMajor( minus.surfaceStress ) ) /
                                   ( 2.0 * eps );
    }

    Matrix3x9RowMajor HFiniteDifference = Matrix3x9RowMajor::Zero();
    Matrix9dRowMajor  ZFiniteDifference = Matrix9dRowMajor::Zero();

    // Perturbing the same entry on the top AND bottom surface gradients changes
    // the average surface gradient <u_{k,l}>_s by exactly +/- eps.
    for ( int entry = 0; entry < 9; ++entry ) {
      double dSurfacePlus[18];
      double dSurfaceMinus[18];
      std::copy( dSurfaceStrain, dSurfaceStrain + 18, dSurfacePlus );
      std::copy( dSurfaceStrain, dSurfaceStrain + 18, dSurfaceMinus );
      dSurfacePlus[entry] += eps;
      dSurfacePlus[9 + entry] += eps;
      dSurfaceMinus[entry] -= eps;
      dSurfaceMinus[9 + entry] -= eps;

      const ExtendedResponse plus = evaluateVirginResponse( material, dU, dSurfacePlus, n.data(), separation.data() );

      const ExtendedResponse minus = evaluateVirginResponse( material, dU, dSurfaceMinus, n.data(), separation.data() );

      HFiniteDifference.col( entry ) = ( plus.force - minus.force ) / ( 2.0 * eps );
      ZFiniteDifference.col(
        entry ) = ( flattenRowMajor( plus.surfaceStress ) - flattenRowMajor( minus.surfaceStress ) ) / ( 2.0 * eps );
    }

    const auto relativeTolerance = []( const auto& reference ) {
      return 1e-7 * std::max( 1.0, reference.template lpNorm< Eigen::Infinity >() );
    };

    assertMatrixNear( analytic.Q,
                      QFiniteDifference,
                      relativeTolerance( QFiniteDifference ),
                      "Analytic Q does not match d(force)/d(jump) finite difference" );
    assertMatrixNear( analytic.K,
                      KFiniteDifference,
                      relativeTolerance( KFiniteDifference ),
                      "Analytic K does not match d(surfaceStress)/d(jump) finite difference" );
    assertMatrixNear( analytic.H,
                      HFiniteDifference,
                      relativeTolerance( HFiniteDifference ),
                      "Analytic H does not match d(force)/d(surfaceGradient) finite difference" );
    assertMatrixNear( analytic.Z,
                      ZFiniteDifference,
                      relativeTolerance( ZFiniteDifference ),
                      "Analytic Z does not match d(surfaceStress)/d(surfaceGradient) finite difference" );
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
    // Construction guards.
    expectInvalidArgument(
      []() {
        const double     tooFewProperties[2] = { 1e5, 0.3 };
        ExtendedMaterial material( "LINEARELASTIC", tooFewProperties, 2, 1 );
      },
      "MarmotExtendedInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h.",
      "Construction with 2 properties" );

    expectInvalidArgument(
      []() {
        const double     zeroThickness[3] = { 1e5, 0.3, 0.0 };
        ExtendedMaterial material( "LINEARELASTIC", zeroThickness, 3, 1 );
      },
      "MarmotExtendedInterfaceMaterialHypoElastic requires h > 0.",
      "Construction with h = 0" );

    expectInvalidArgument(
      []() {
        const double     negativeThickness[3] = { 1e5, 0.3, -0.01 };
        ExtendedMaterial material( "LINEARELASTIC", negativeThickness, 3, 1 );
      },
      "MarmotExtendedInterfaceMaterialHypoElastic requires h > 0.",
      "Construction with h < 0" );

    // Stress-update geometry guards.
    const double     interfaceProperties[3] = { 1e5, 0.3, 0.01 };
    ExtendedMaterial material( "LINEARELASTIC", interfaceProperties, 3, 1 );

    const double dU[6]              = { 1e-4, 0., 0., 0., 0., 0. };
    const double dSurfaceStrain[18] = { 0. };

    expectInvalidArgument(
      [&]() {
        const double zeroNormal[3] = { 0., 0., 0. };
        evaluateVirginResponse( material, dU, dSurfaceStrain, zeroNormal, nullptr );
      },
      "MarmotExtendedInterfaceMaterialHypoElastic: interface normal is zero.",
      "Stress update with zero interface normal" );

    const double normal[3] = { 0., 0., 1. };

    expectInvalidArgument(
      [&]() {
        const double tangentialOnlySeparation[3] = { 0.003, 0., 0. };
        evaluateVirginResponse( material, dU, dSurfaceStrain, normal, tangentialOnlySeparation );
      },
      "MarmotExtendedInterfaceMaterialHypoElastic: the top-bottom connector must have a positive normal component.",
      "Stress update with purely tangential connector" );

    expectInvalidArgument(
      [&]() {
        const double invertedSeparation[3] = { 0.001, 0., -0.02 };
        evaluateVirginResponse( material, dU, dSurfaceStrain, normal, invertedSeparation );
      },
      "MarmotExtendedInterfaceMaterialHypoElastic: the top-bottom connector must have a positive normal component.",
      "Stress update with inverted connector" );
  }

  void testDensityDelegation()
  {
    const double     interfaceProperties[9] = { 1e8, 0.3, 0.01, 2e7, 0.25, 6., 1e-4, 1., 2400. };
    ExtendedMaterial material( "LINEARVISCOELASTICWIECHERT", interfaceProperties, 9, 1 );

    throwExceptionOnFailure( checkIfEqual( material.getDensity(), interfaceProperties[8] ),
                             "Extended interface density delegation failed." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testZeroSeparationLinearElasticAgainstBulk,
                                                   testZeroSeparationVonMisesAgainstBulk,
                                                   testExplicitZeroSeparationMatchesCoincidentFacePath,
                                                   testZeroSeparationMatchesPlainInterfaceMaterialStress,
                                                   testNonzeroSeparationReconstructedGeometryMatchesReference,
                                                   testNonzeroSeparationTangentsMatchFiniteDifferences,
                                                   testDegenerateGeometryAndConstructionThrow,
                                                   testDensityDelegation };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
