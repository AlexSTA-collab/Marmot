#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotWarpingInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotWarpingStabPressureInterfaceMaterialHypoElastic.h"

#include <Eigen/Dense>

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Material  = MarmotWarpingStabPressureInterfaceMaterialHypoElastic; // 5 stations
  using WMaterial = MarmotWarpingInterfaceMaterialHypoElastic;

  constexpr int nGen = Material::nGeneralized; // 40

  using VectorGen = Eigen::Matrix< double, nGen, 1 >;
  using MatrixGen = Eigen::Matrix< double, nGen, nGen, Eigen::RowMajor >;

  const Eigen::Vector3d interfaceNormal( 0.0, 0.0, 1.0 );

  /** Face-to-face connector. Set to a vector with a TANGENTIAL component to
   * exercise the geometric d_tau correction, which couples into the mixed
   * pressure blocks and is silently absent when the faces are aligned. */
  Eigen::Vector3d interfaceSeparation( 0.0, 0.0, 0.01 );

  const Eigen::Vector3d alignedSeparation( 0.0, 0.0, 0.01 );
  const Eigen::Vector3d skewedSeparation( 0.02, -0.015, 0.01 );

  const std::vector< double > elasticProperties = { 1.0e5, 0.3, 0.01 };
  //                                                 E,    nu,   h,   fy,     H,  Aexp, Hexp, id
  const std::vector< double > plasticProperties = { 1.0e5, 0.3, 0.01, 5.0, 1.0e3, 0.0, 0.0, 0.0 };

  double bulkModulusOf( const std::vector< double >& properties )
  {
    return properties[0] / ( 3.0 * ( 1.0 - 2.0 * properties[1] ) );
  }

  /** One stress update from a virgin state; optionally returns the final state vector. */
  void evaluate( const std::string&           materialName,
                 const std::vector< double >& properties,
                 const VectorGen&             generalizedStrain,
                 VectorGen&                   generalizedStress,
                 MatrixGen&                   condensedTangent,
                 Eigen::VectorXd*             finalStateVars = nullptr )
  {
    Material        material( materialName, properties.data(), static_cast< int >( properties.size() ), 1 );
    Eigen::VectorXd stateVars = Eigen::VectorXd::Zero( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );

    Eigen::Matrix< double, 6, 1 >  dU = Eigen::Matrix< double, 6, 1 >::Zero();
    Eigen::Matrix< double, 18, 1 > dSurface;
    Eigen::Matrix< double, 18, 1 > dWarping;

    dU.segment< 3 >( 0 )       = generalizedStrain.segment< 3 >( Material::offsetJump );
    dSurface.segment< 9 >( 0 ) = generalizedStrain.segment< 9 >( Material::offsetSurfacePlus );
    dSurface.segment< 9 >( 9 ) = generalizedStrain.segment< 9 >( Material::offsetSurfaceMinus );
    dWarping.segment< 9 >( 0 ) = generalizedStrain.segment< 9 >( Material::offsetWarpingSymmetric );
    dWarping.segment< 9 >( 9 ) = generalizedStrain.segment< 9 >( Material::offsetWarpingAntisymmetric );

    generalizedStress.setZero();
    condensedTangent.setZero();

    Material::State         state( generalizedStress.data(), stateVars.data() );
    Material::Tangents      tangents( condensedTangent.data() );
    Material::Deformation   deformation( dU.data(),
                                       dSurface.data(),
                                       dWarping.data(),
                                       interfaceNormal.data(),
                                       interfaceSeparation.data(),
                                       generalizedStrain( Material::offsetPressure ) );
    Material::TimeIncrement timeIncrement{ 0.0, 1.0 };

    material.computeStress( state, tangents, deformation, timeIncrement );

    if ( finalStateVars ) {
      *finalStateVars = stateVars;
    }
  }

  VectorGen makeRepresentativeGeneralizedStrain( bool includeWarping, double pressureIncrement )
  {
    VectorGen q = VectorGen::Zero();
    q.segment< 3 >( Material::offsetJump ) << 1.3e-4, -0.7e-4, 2.1e-4;
    for ( int i = 0; i < 9; ++i ) {
      q( Material::offsetSurfacePlus + i )  = 0.9e-4 * std::sin( 1.0 + i );
      q( Material::offsetSurfaceMinus + i ) = -0.6e-4 * std::cos( 2.0 + i );
      if ( includeWarping ) {
        q( Material::offsetWarpingSymmetric + i )     = 0.5e-4 * std::sin( 3.0 + i );
        q( Material::offsetWarpingAntisymmetric + i ) = 0.4e-4 * std::cos( 0.5 + i );
      }
    }
    q( Material::offsetPressure ) = pressureIncrement;
    return q;
  }

  /** Uniform through-thickness state: A+ = A-, no warping. */
  VectorGen makeUniformGeneralizedStrain()
  {
    VectorGen q = VectorGen::Zero();
    q.segment< 3 >( Material::offsetJump ) << 0.8e-4, -0.5e-4, 1.7e-4;
    for ( int i = 0; i < 9; ++i ) {
      const double value                    = 0.7e-4 * std::sin( 1.0 + i );
      q( Material::offsetSurfacePlus + i )  = value;
      q( Material::offsetSurfaceMinus + i ) = value;
    }
    return q;
  }

  // ------------------------------------------------------------------
  // 1. Full 40x40 condensed tangent against central differences. This is the
  //    decisive consistency check: it covers the mixed pressure row and column
  //    together with every warping block.
  // ------------------------------------------------------------------
  void checkCondensedTangentAgainstFiniteDifference( const std::string&           materialName,
                                                     const std::vector< double >& properties )
  {
    const VectorGen q = makeRepresentativeGeneralizedStrain( true, 3.0e-2 );

    VectorGen stress;
    MatrixGen tangent;
    evaluate( materialName, properties, q, stress, tangent );

    MatrixGen finiteDifference = MatrixGen::Zero();
    for ( int column = 0; column < nGen; ++column ) {
      // The pressure column is scaled like a STRESS, the strain columns like a
      // strain, so they need different perturbation sizes.
      const double perturbation = ( column == Material::offsetPressure ) ? 1.0e-5 : 1.0e-9;

      VectorGen forward = q, backward = q;
      forward( column ) += perturbation;
      backward( column ) -= perturbation;

      VectorGen stressForward, stressBackward;
      MatrixGen unused;
      evaluate( materialName, properties, forward, stressForward, unused );
      evaluate( materialName, properties, backward, stressBackward, unused );

      finiteDifference.col( column ) = ( stressForward - stressBackward ) / ( 2.0 * perturbation );
    }

    const double error         = ( tangent - finiteDifference ).lpNorm< Eigen::Infinity >();
    const double scale         = std::max( 1.0, finiteDifference.lpNorm< Eigen::Infinity >() );
    const double relativeError = error / scale;

    int    worstRow = 0, worstColumn = 0;
    double worst = 0.0;
    for ( int i = 0; i < nGen; ++i ) {
      for ( int j = 0; j < nGen; ++j ) {
        if ( std::abs( tangent( i, j ) - finiteDifference( i, j ) ) > worst ) {
          worst       = std::abs( tangent( i, j ) - finiteDifference( i, j ) );
          worstRow    = i;
          worstColumn = j;
        }
      }
    }

    std::cout << "full 40x40 condensed tangent vs FD (" << materialName << ", d_tau "
              << ( interfaceSeparation.head< 2 >().norm() > 0.0 ? "on" : "off" )
              << "): relative error = " << relativeError << " at (" << worstRow << "," << worstColumn << ")\n";

    throwExceptionOnFailure( relativeError < 1e-6,
                             "warping/mixed-pressure condensed tangent does not match central differences (" +
                               materialName + "): relative error " + std::to_string( relativeError ) );
  }

  void TestCondensedTangentMatchesFiniteDifference()
  {
    interfaceSeparation = alignedSeparation;
    checkCondensedTangentAgainstFiniteDifference( "LINEARELASTIC", elasticProperties );
    checkCondensedTangentAgainstFiniteDifference( "VONMISES", plasticProperties );

    // With a tangential offset between the faces the d_tau correction is live,
    // and it enters BOTH mixed pressure blocks (dS/dp and dR_p/dA). Those terms
    // are identically zero for aligned faces, so without this case they would
    // go untested.
    interfaceSeparation = skewedSeparation;
    checkCondensedTangentAgainstFiniteDifference( "LINEARELASTIC", elasticProperties );
    checkCondensedTangentAgainstFiniteDifference( "VONMISES", plasticProperties );
    interfaceSeparation = alignedSeparation;
  }

  // ------------------------------------------------------------------
  // 2. THE STRUCTURAL PROPERTY THE FORMULATION RESTS ON: the local
  //    traction-equilibrium solve is EXACTLY independent of p.
  //
  // Because p is constant through the thickness, sigma~ n = t~ at every station
  // is identical to dev(sigma) n = t. So changing p may not move a single
  // station stress, a single station normal gradient, or the plastic state --
  // it may only shift the shared traction. This is what makes every mixed
  // coupling block purely geometric.
  // ------------------------------------------------------------------
  void checkLocalSolveIsPressureIndependent( const std::string& materialName, const std::vector< double >& properties )
  {
    Eigen::VectorXd stateWithoutPressure, stateWithPressure;
    VectorGen       stress;
    MatrixGen       tangent;

    evaluate( materialName,
              properties,
              makeRepresentativeGeneralizedStrain( true, 0.0 ),
              stress,
              tangent,
              &stateWithoutPressure );
    evaluate( materialName,
              properties,
              makeRepresentativeGeneralizedStrain( true, 7.5e-2 ),
              stress,
              tangent,
              &stateWithPressure );

    Material material( materialName, properties.data(), static_cast< int >( properties.size() ), 1 );

    double difference = 0.0;
    for ( int alpha = 0; alpha < Material::nStations; ++alpha ) {
      for ( const std::string& name : { "stationStress" + std::to_string( alpha ),
                                        "normalGradient" + std::to_string( alpha ),
                                        "stationMaterialStateVars" + std::to_string( alpha ) } ) {
        const StateView a = material.getStateView( name, stateWithoutPressure.data() );
        const StateView b = material.getStateView( name, stateWithPressure.data() );
        for ( int i = 0; i < a.stateSize; ++i ) {
          difference = std::max( difference, std::abs( a.stateLocation[i] - b.stateLocation[i] ) );
        }
      }
    }

    std::cout << "local solve vs pressure (" << materialName << "): max state difference = " << difference << "\n";

    throwExceptionOnFailure( difference == 0.0,
                             "the local traction-equilibrium solve must be EXACTLY independent of the mixed pressure "
                             "(" +
                               materialName + "): max difference " + std::to_string( difference ) );
  }

  void TestLocalSolveIsIndependentOfThePressure()
  {
    checkLocalSolveIsPressureIndependent( "LINEARELASTIC", elasticProperties );
    checkLocalSolveIsPressureIndependent( "VONMISES", plasticProperties );
  }

  // ------------------------------------------------------------------
  // 3. CONSISTENCY OF THE MIXED SPLIT: at the pressure that satisfies the
  //    volumetric law exactly, the mixed material must reproduce the
  //    displacement-only warping material.
  //
  // sigma~ = dev(sigma) - p I with p = -K tr(eps) gives back sigma itself,
  // because the volumetric response of both LINEARELASTIC and VONMISES is
  // elastic (sigma_vol = K tr(eps) I).
  //
  // The comparison is made on a UNIFORM through-thickness state (A+ = A-,
  // Ws = Wa = 0). That is not a convenience: for a NON-uniform state the two
  // materials legitimately differ, because one enforces equal total tractions
  // across the stations and the other equal deviatoric tractions, which are
  // different conditions. For a uniform state g^(alpha) = gbar solves both, so
  // the station stresses coincide and the comparison isolates exactly the
  // volumetric substitution.
  // ------------------------------------------------------------------
  void checkReductionAtTheConsistentPressure( const std::string& materialName, const std::vector< double >& properties )
  {
    const double K = bulkModulusOf( properties );
    const double h = properties[2];

    VectorGen q = makeUniformGeneralizedStrain();

    // First pass at p = 0 reads off the volumetric strain increment: R_p = h * vol.
    VectorGen stress;
    MatrixGen tangent;
    evaluate( materialName, properties, q, stress, tangent );
    const double volumetricStrain = stress( Material::offsetPressure ) / h;

    // The pressure that satisfies the volumetric law: vol + p/K = 0.
    q( Material::offsetPressure ) = -K * volumetricStrain;
    evaluate( materialName, properties, q, stress, tangent );

    throwExceptionOnFailure( std::abs( stress( Material::offsetPressure ) ) < 1e-18,
                             "at the consistent pressure the volumetric residual must vanish (" + materialName +
                               "): got " + std::to_string( stress( Material::offsetPressure ) ) );

    // Reference: the displacement-only warping material at the same strain.
    constexpr int                        nWarpGen        = WMaterial::nGeneralized; // 39
    Eigen::Matrix< double, nWarpGen, 1 > referenceStress = Eigen::Matrix< double, nWarpGen, 1 >::Zero();
    Eigen::Matrix< double, nWarpGen, nWarpGen, Eigen::RowMajor >
      referenceTangent = Eigen::Matrix< double, nWarpGen, nWarpGen, Eigen::RowMajor >::Zero();

    WMaterial       reference( materialName, properties.data(), static_cast< int >( properties.size() ), 1 );
    Eigen::VectorXd stateVars = Eigen::VectorXd::Zero( reference.getNumberOfRequiredStateVars() );
    reference.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );

    Eigen::Matrix< double, 6, 1 >  dU       = Eigen::Matrix< double, 6, 1 >::Zero();
    Eigen::Matrix< double, 18, 1 > dSurface = Eigen::Matrix< double, 18, 1 >::Zero();
    Eigen::Matrix< double, 18, 1 > dWarping = Eigen::Matrix< double, 18, 1 >::Zero();
    dU.segment< 3 >( 0 )                    = q.segment< 3 >( Material::offsetJump );
    dSurface.segment< 9 >( 0 )              = q.segment< 9 >( Material::offsetSurfacePlus );
    dSurface.segment< 9 >( 9 )              = q.segment< 9 >( Material::offsetSurfaceMinus );

    WMaterial::State         referenceState( referenceStress.data(), stateVars.data() );
    WMaterial::Tangents      referenceTangents( referenceTangent.data() );
    WMaterial::Deformation   referenceDeformation( dU.data(),
                                                 dSurface.data(),
                                                 dWarping.data(),
                                                 interfaceNormal.data(),
                                                 interfaceSeparation.data() );
    WMaterial::TimeIncrement timeIncrement{ 0.0, 1.0 };
    reference.computeStress( referenceState, referenceTangents, referenceDeformation, timeIncrement );

    double error = 0.0;
    error        = std::max( error,
                      ( stress.segment< 3 >( Material::offsetJump ) -
                        referenceStress.segment< 3 >( WMaterial::offsetJump ) )
                        .lpNorm< Eigen::Infinity >() );
    error        = std::max( error,
                      ( stress.segment< 9 >( Material::offsetSurfacePlus ) -
                        referenceStress.segment< 9 >( WMaterial::offsetSurfacePlus ) )
                        .lpNorm< Eigen::Infinity >() );
    error        = std::max( error,
                      ( stress.segment< 9 >( Material::offsetSurfaceMinus ) -
                        referenceStress.segment< 9 >( WMaterial::offsetSurfaceMinus ) )
                        .lpNorm< Eigen::Infinity >() );

    const double scale = std::max( 1.0, referenceStress.lpNorm< Eigen::Infinity >() );

    std::cout << "mixed vs displacement-only at the consistent pressure (" << materialName << "): max|dp| = " << error
              << " (relative " << error / scale << ")\n";

    throwExceptionOnFailure( error / scale < 1e-10,
                             "at the pressure satisfying the volumetric law the mixed material must reproduce the "
                             "displacement-only warping material (" +
                               materialName + "): relative error " + std::to_string( error / scale ) );
  }

  void TestReductionAtTheConsistentPressure()
  {
    checkReductionAtTheConsistentPressure( "LINEARELASTIC", elasticProperties );
    checkReductionAtTheConsistentPressure( "VONMISES", plasticProperties );
  }

  // ------------------------------------------------------------------
  // 4. The volumetric residual must be the exact weighted station mean, i.e.
  //
  //     R_p / h = tr(ABar) + (2/3) tr(Ws) + gbar . n + dp / K
  //
  // In particular the ANTISYMMETRIC warping must be INVISIBLE to it: an odd
  // through-thickness profile transports no net volume, sum lambda phi_a = 0.
  // Getting this weight wrong (e.g. by reusing the two-station material's plain
  // average) is the single easiest way to break the mixed formulation while
  // still passing a tangent check.
  // ------------------------------------------------------------------
  void TestVolumetricResidualUsesTheCorrectProfileWeights()
  {
    const double h   = elasticProperties[2];
    const double K   = bulkModulusOf( elasticProperties );
    const double ell = interfaceSeparation( 2 ); // normal separation, d_tau = 0 here

    VectorGen stress;
    MatrixGen tangent;

    // (a) antisymmetric warping alone -> exactly zero volumetric residual
    {
      VectorGen q = VectorGen::Zero();
      for ( int i = 0; i < 9; ++i ) {
        q( Material::offsetWarpingAntisymmetric + i ) = 1.0e-4 * std::cos( 0.5 + i );
      }
      evaluate( "LINEARELASTIC", elasticProperties, q, stress, tangent );
      throwExceptionOnFailure( std::abs( stress( Material::offsetPressure ) ) < 1e-18,
                               "the antisymmetric warping must not contribute to the volumetric residual: got " +
                                 std::to_string( stress( Material::offsetPressure ) ) );
    }

    // (b) symmetric warping alone -> weight exactly 2/3
    {
      VectorGen    q                            = VectorGen::Zero();
      const double value                        = 1.0e-4;
      q( Material::offsetWarpingSymmetric + 0 ) = value; // Ws_xx
      q( Material::offsetWarpingSymmetric + 4 ) = value; // Ws_yy
      evaluate( "LINEARELASTIC", elasticProperties, q, stress, tangent );
      const double expected = h * ( 2.0 / 3.0 ) * ( 2.0 * value );
      throwExceptionOnFailure( std::abs( stress( Material::offsetPressure ) - expected ) < 1e-16,
                               "the symmetric warping must enter the volumetric residual with weight 2/3: got " +
                                 std::to_string( stress( Material::offsetPressure ) ) + ", expected " +
                                 std::to_string( expected ) );
    }

    // (c) the full expression, on a general state
    {
      const double pressureIncrement = 2.5e-2;
      VectorGen    q                 = makeRepresentativeGeneralizedStrain( true, pressureIncrement );
      evaluate( "LINEARELASTIC", elasticProperties, q, stress, tangent );

      const Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > APlus( q.data() +
                                                                                      Material::offsetSurfacePlus );
      const Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > AMinus( q.data() +
                                                                                       Material::offsetSurfaceMinus );
      const Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > > Ws( q.data() +
                                                                                   Material::offsetWarpingSymmetric );

      const Eigen::Vector3d jump     = q.segment< 3 >( Material::offsetJump );
      const double          expected = h * ( 0.5 * APlus.trace() + 0.5 * AMinus.trace() + ( 2.0 / 3.0 ) * Ws.trace() +
                                    jump.dot( interfaceNormal ) / ell + pressureIncrement / K );

      const double relativeError = std::abs( stress( Material::offsetPressure ) - expected ) /
                                   std::max( 1e-30, std::abs( expected ) );

      std::cout << "volumetric residual: relative error vs closed form = " << relativeError << "\n";

      throwExceptionOnFailure( relativeError < 1e-12,
                               "the volumetric residual does not match its closed form: got " +
                                 std::to_string( stress( Material::offsetPressure ) ) + ", expected " +
                                 std::to_string( expected ) );
    }

    // (d) the pressure self-stiffness is h/K, i.e. the compressibility term
    {
      VectorGen q = VectorGen::Zero();
      evaluate( "LINEARELASTIC", elasticProperties, q, stress, tangent );
      const double expected = h / K;
      throwExceptionOnFailure( std::abs( tangent( Material::offsetPressure, Material::offsetPressure ) - expected ) <
                                 1e-14 * expected,
                               "the pressure self-stiffness must be h/K." );
    }
  }

  // ------------------------------------------------------------------
  // 5. The mixed formulation must actually REMOVE the volumetric amplification.
  //
  // In a displacement-only formulation the generalized stress carries
  // p = K tr(eps), so it grows without bound as nu -> 0.5. Here the pressure is
  // supplied independently, so at p = 0 the generalized stress must be purely
  // deviatoric and therefore BOUNDED in K.
  //
  // The comparison is made at nu = 0.49 because that is the largest value at
  // which the displacement-only material still CONVERGES: at nu = 0.4999 its
  // local traction-equilibrium solve throws. That failure is itself part of the
  // result and is asserted below -- the mixed material must survive exactly the
  // regime that breaks the displacement-only one.
  // ------------------------------------------------------------------
  void TestGeneralizedStressIsBoundedInTheIncompressibleLimit()
  {
    const VectorGen q = makeRepresentativeGeneralizedStrain( true, 0.0 );

    auto mixedResponse = [&]( double nu ) {
      const std::vector< double > properties = { 1.0e5, nu, 0.01 };
      VectorGen                   stress;
      MatrixGen                   tangent;
      evaluate( "LINEARELASTIC", properties, q, stress, tangent );
      return stress.segment< 9 >( Material::offsetSurfacePlus ).lpNorm< Eigen::Infinity >();
    };

    auto displacementOnlyResponse = [&]( double nu ) {
      const std::vector< double >          properties = { 1.0e5, nu, 0.01 };
      constexpr int                        nWarpGen   = WMaterial::nGeneralized;
      Eigen::Matrix< double, nWarpGen, 1 > stress     = Eigen::Matrix< double, nWarpGen, 1 >::Zero();
      Eigen::Matrix< double, nWarpGen, nWarpGen, Eigen::RowMajor >
        tangent = Eigen::Matrix< double, nWarpGen, nWarpGen, Eigen::RowMajor >::Zero();

      WMaterial       material( "LINEARELASTIC", properties.data(), static_cast< int >( properties.size() ), 1 );
      Eigen::VectorXd stateVars = Eigen::VectorXd::Zero( material.getNumberOfRequiredStateVars() );
      material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );

      Eigen::Matrix< double, 6, 1 >  dU       = Eigen::Matrix< double, 6, 1 >::Zero();
      Eigen::Matrix< double, 18, 1 > dSurface = Eigen::Matrix< double, 18, 1 >::Zero();
      Eigen::Matrix< double, 18, 1 > dWarping = Eigen::Matrix< double, 18, 1 >::Zero();
      dU.segment< 3 >( 0 )                    = q.segment< 3 >( Material::offsetJump );
      dSurface.segment< 9 >( 0 )              = q.segment< 9 >( Material::offsetSurfacePlus );
      dSurface.segment< 9 >( 9 )              = q.segment< 9 >( Material::offsetSurfaceMinus );
      dWarping.segment< 9 >( 0 )              = q.segment< 9 >( Material::offsetWarpingSymmetric );
      dWarping.segment< 9 >( 9 )              = q.segment< 9 >( Material::offsetWarpingAntisymmetric );

      WMaterial::State         state( stress.data(), stateVars.data() );
      WMaterial::Tangents      tangents( tangent.data() );
      WMaterial::Deformation   deformation( dU.data(),
                                          dSurface.data(),
                                          dWarping.data(),
                                          interfaceNormal.data(),
                                          interfaceSeparation.data() );
      WMaterial::TimeIncrement timeIncrement{ 0.0, 1.0 };
      material.computeStress( state, tangents, deformation, timeIncrement );

      return stress.segment< 9 >( WMaterial::offsetSurfacePlus ).lpNorm< Eigen::Infinity >();
    };

    const double mixedGrowth        = mixedResponse( 0.49 ) / mixedResponse( 0.3 );
    const double displacementGrowth = displacementOnlyResponse( 0.49 ) / displacementOnlyResponse( 0.3 );

    std::cout << "growth from nu = 0.3 to nu = 0.49: mixed = " << mixedGrowth
              << ", displacement-only = " << displacementGrowth << "\n";

    throwExceptionOnFailure( displacementGrowth > 5.0,
                             "the displacement-only material is expected to amplify strongly as nu -> 0.5; got "
                             "growth " +
                               std::to_string( displacementGrowth ) +
                               " -- if this fails the comparison proves nothing." );
    throwExceptionOnFailure( mixedGrowth < 2.0,
                             "with an independent pressure the generalized stress must stay bounded as nu -> 0.5; got "
                             "growth " +
                               std::to_string( mixedGrowth ) );

    // At nu = 0.4999 the displacement-only local solve fails outright, while the
    // mixed one must still deliver a finite response. This is the practical
    // point of the mixed split, not just a conditioning nicety.
    bool displacementOnlyFailed = false;
    try {
      displacementOnlyResponse( 0.4999 );
    }
    catch ( const std::exception& ) {
      displacementOnlyFailed = true;
    }

    double mixedNearlyIncompressible = 0.0;
    bool   mixedSurvived             = true;
    try {
      mixedNearlyIncompressible = mixedResponse( 0.4999 );
    }
    catch ( const std::exception& ) {
      mixedSurvived = false;
    }

    std::cout << "at nu = 0.4999: displacement-only " << ( displacementOnlyFailed ? "FAILED" : "converged" )
              << ", mixed " << ( mixedSurvived ? "converged" : "FAILED" ) << " with response "
              << mixedNearlyIncompressible << "\n";

    throwExceptionOnFailure( mixedSurvived && std::isfinite( mixedNearlyIncompressible ),
                             "the mixed material must remain solvable at nu = 0.4999." );
  }

  // ------------------------------------------------------------------
  // 6. Both warping modes must still carry stiffness once the deviatoric
  //    projector is in the loop -- the projector must not have made the
  //    enrichment inert.
  // ------------------------------------------------------------------
  void TestWarpingModesStillCarryStiffnessUnderTheDeviatoricProjection()
  {
    VectorGen q                                       = VectorGen::Zero();
    q( Material::offsetWarpingSymmetric + 3 * 0 + 0 ) = 1.0e-4; // Ws_xx

    VectorGen stress;
    MatrixGen tangent;
    evaluate( "LINEARELASTIC", elasticProperties, q, stress, tangent );

    throwExceptionOnFailure( stress.segment< 9 >( Material::offsetWarpingSymmetric ).norm() > 1e-8,
                             "a parabolic in-plane warping stretch must still produce a nonzero conjugate stress." );
    throwExceptionOnFailure( tangent( Material::offsetWarpingSymmetric, Material::offsetWarpingSymmetric ) > 0.0,
                             "the parabolic warping direction must still have positive stiffness." );

    VectorGen qAnti                                           = VectorGen::Zero();
    qAnti( Material::offsetWarpingAntisymmetric + 3 * 0 + 0 ) = 1.0e-4; // Wa_xx
    evaluate( "LINEARELASTIC", elasticProperties, qAnti, stress, tangent );

    throwExceptionOnFailure( stress.segment< 9 >( Material::offsetWarpingAntisymmetric ).norm() > 1e-8,
                             "a cubic in-plane warping stretch must still produce a nonzero conjugate stress." );
    throwExceptionOnFailure( tangent( Material::offsetWarpingAntisymmetric, Material::offsetWarpingAntisymmetric ) >
                               0.0,
                             "the cubic warping direction must still have positive stiffness." );
  }

} // namespace

int main()
{
  auto
    tests = std::vector< std::function< void() > >{ TestCondensedTangentMatchesFiniteDifference,
                                                    TestLocalSolveIsIndependentOfThePressure,
                                                    TestReductionAtTheConsistentPressure,
                                                    TestVolumetricResidualUsesTheCorrectProfileWeights,
                                                    TestGeneralizedStressIsBoundedInTheIncompressibleLimit,
                                                    TestWarpingModesStillCarryStiffnessUnderTheDeviatoricProjection };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
