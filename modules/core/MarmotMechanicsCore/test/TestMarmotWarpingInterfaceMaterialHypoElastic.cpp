#include "Marmot/MarmotGaussLobattoInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotTesting.h"
#include "Marmot/MarmotWarpingInterfaceMaterialHypoElastic.h"

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using WMaterial  = MarmotWarpingInterfaceMaterialHypoElastic; // 5 stations
  using GLMaterial = MarmotGaussLobattoInterfaceMaterialHypoElastic;

  constexpr int nGen = WMaterial::nGeneralized; // 39

  using VectorGen        = Eigen::Matrix< double, nGen, 1 >;
  using MatrixGen        = Eigen::Matrix< double, nGen, nGen, Eigen::RowMajor >;
  using Matrix3dRowMajor = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;

  const Eigen::Vector3d interfaceNormal( 0.0, 0.0, 1.0 );
  const Eigen::Vector3d interfaceSeparation( 0.0, 0.0, 0.01 );

  const std::vector< double > elasticProperties = { 1.0e5, 0.3, 0.01 };
  //                                      E,     nu,  h,    fy,   H,     Aexp, Hexp, id
  const std::vector< double > plasticProperties = { 1.0e5, 0.3, 0.01, 50.0, 1.0e3, 0.0, 0.0, 0.0 };

  /** One warping-material stress update from a virgin state. */
  void evaluateWarpingMaterial( const std::string&           materialName,
                                const std::vector< double >& properties,
                                const VectorGen&             generalizedStrain,
                                VectorGen&                   generalizedStress,
                                MatrixGen&                   condensedTangent )
  {
    WMaterial       material( materialName, properties.data(), static_cast< int >( properties.size() ), 1 );
    Eigen::VectorXd stateVars = Eigen::VectorXd::Zero( material.getNumberOfRequiredStateVars() );
    material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );

    Eigen::Matrix< double, 6, 1 >  dU = Eigen::Matrix< double, 6, 1 >::Zero();
    Eigen::Matrix< double, 18, 1 > dSurface;
    Eigen::Matrix< double, 18, 1 > dWarping;

    // The material forms w = u+ - u-; place the whole jump on the + side.
    dU.segment< 3 >( 0 )       = generalizedStrain.segment< 3 >( WMaterial::offsetJump );
    dSurface.segment< 9 >( 0 ) = generalizedStrain.segment< 9 >( WMaterial::offsetSurfacePlus );
    dSurface.segment< 9 >( 9 ) = generalizedStrain.segment< 9 >( WMaterial::offsetSurfaceMinus );
    dWarping.segment< 9 >( 0 ) = generalizedStrain.segment< 9 >( WMaterial::offsetWarpingSymmetric );
    dWarping.segment< 9 >( 9 ) = generalizedStrain.segment< 9 >( WMaterial::offsetWarpingAntisymmetric );

    generalizedStress.setZero();
    condensedTangent.setZero();

    WMaterial::State         state( generalizedStress.data(), stateVars.data() );
    WMaterial::Tangents      tangents( condensedTangent.data() );
    WMaterial::Deformation   deformation( dU.data(),
                                        dSurface.data(),
                                        dWarping.data(),
                                        interfaceNormal.data(),
                                        interfaceSeparation.data() );
    WMaterial::TimeIncrement timeIncrement{ 0.0, 1.0 };

    material.computeStress( state, tangents, deformation, timeIncrement );
  }

  VectorGen makeRepresentativeGeneralizedStrain( bool includeWarping )
  {
    VectorGen q = VectorGen::Zero();
    q.segment< 3 >( WMaterial::offsetJump ) << 1.3e-4, -0.7e-4, 2.1e-4;
    for ( int i = 0; i < 9; ++i ) {
      q( WMaterial::offsetSurfacePlus + i )  = 0.9e-4 * std::sin( 1.0 + i );
      q( WMaterial::offsetSurfaceMinus + i ) = -0.6e-4 * std::cos( 2.0 + i );
      if ( includeWarping ) {
        q( WMaterial::offsetWarpingSymmetric + i )     = 0.5e-4 * std::sin( 3.0 + i );
        q( WMaterial::offsetWarpingAntisymmetric + i ) = 0.4e-4 * std::cos( 0.5 + i );
      }
    }
    return q;
  }

  // ------------------------------------------------------------------
  // 1. The through-thickness quadrature identities the formulation rests on.
  //
  // The warping profiles must have ZERO weighted-mean derivative, otherwise
  // their normal-gradient content would not be spanned by the free station
  // gradients g^(alpha) and the enrichment would be redundant (see the
  // material header). They must also be mutually orthogonal, so that the
  // symmetric and antisymmetric modes do not couple.
  // ------------------------------------------------------------------
  template < int NStations >
  void checkQuadratureIdentitiesForStationCount()
  {
    using Rule        = Marmot::Detail::GaussLobattoInterface::LobattoRule< NStations >;
    const auto xi     = Rule::xi();
    const auto weight = Rule::weight();

    double sumWeights = 0.0, sumDSym = 0.0, sumDAnti = 0.0, sumSym = 0.0, sumAnti = 0.0, sumCross = 0.0;
    for ( int alpha = 0; alpha < NStations; ++alpha ) {
      const double lambda = 0.5 * weight[alpha];
      const double x      = xi[alpha];
      sumWeights += lambda;
      sumSym += lambda * WMaterial::warpingProfileSymmetric( x );
      sumAnti += lambda * WMaterial::warpingProfileAntisymmetric( x );
      sumDSym += lambda * ( -2.0 * x );
      sumDAnti += lambda * ( 1.0 - 3.0 * x * x );
      sumCross += lambda * WMaterial::warpingProfileSymmetric( x ) * WMaterial::warpingProfileAntisymmetric( x );
    }

    const std::string tag = "Lobatto rule n=" + std::to_string( NStations ) + ": ";
    throwExceptionOnFailure( std::abs( sumWeights - 1.0 ) < 1e-13, tag + "normalized weights must sum to 1." );
    throwExceptionOnFailure( std::abs( sumDSym ) < 1e-13,
                             tag + "sum(lambda phi_s') must vanish -- otherwise the symmetric warping duplicates the "
                                   "station normal gradients." );
    throwExceptionOnFailure( std::abs( sumDAnti ) < 1e-13,
                             tag + "sum(lambda phi_a') must vanish -- otherwise the antisymmetric warping duplicates "
                                   "the station normal gradients." );
    throwExceptionOnFailure( std::abs( sumAnti ) < 1e-13, tag + "sum(lambda phi_a) must vanish (odd profile)." );
    throwExceptionOnFailure( std::abs( sumSym - 2.0 / 3.0 ) < 1e-13, tag + "sum(lambda phi_s) must equal 2/3." );
    throwExceptionOnFailure( std::abs( sumCross ) < 1e-13,
                             tag + "the two warping profiles must be orthogonal through the thickness." );
  }

  void TestWarpingProfileQuadratureIdentities()
  {
    checkQuadratureIdentitiesForStationCount< 4 >();
    checkQuadratureIdentitiesForStationCount< 5 >();
    checkQuadratureIdentitiesForStationCount< 7 >();

    // The profiles must vanish at both faces, or the warping would change the
    // face displacements and stop being condensable.
    throwExceptionOnFailure( std::abs( WMaterial::warpingProfileSymmetric( 1.0 ) ) < 1e-15 &&
                               std::abs( WMaterial::warpingProfileSymmetric( -1.0 ) ) < 1e-15,
                             "the symmetric warping profile must vanish at both faces." );
    throwExceptionOnFailure( std::abs( WMaterial::warpingProfileAntisymmetric( 1.0 ) ) < 1e-15 &&
                               std::abs( WMaterial::warpingProfileAntisymmetric( -1.0 ) ) < 1e-15,
                             "the antisymmetric warping profile must vanish at both faces." );
  }

  // ------------------------------------------------------------------
  // 2. With Ws = Wa = 0 the material must reproduce the Gauss-Lobatto
  //    material EXACTLY -- same generalized stresses, same 21x21 tangent.
  // ------------------------------------------------------------------
  void checkReductionToGaussLobatto( const std::string& materialName, const std::vector< double >& properties )
  {
    const VectorGen q = makeRepresentativeGeneralizedStrain( false );

    VectorGen warpingStress;
    MatrixGen warpingTangent;
    evaluateWarpingMaterial( materialName, properties, q, warpingStress, warpingTangent );

    GLMaterial      reference( materialName, properties.data(), static_cast< int >( properties.size() ), 1 );
    Eigen::VectorXd stateVars = Eigen::VectorXd::Zero( reference.getNumberOfRequiredStateVars() );
    reference.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );

    Eigen::Vector3d                                force              = Eigen::Vector3d::Zero();
    Matrix3dRowMajor                               surfaceStressPlus  = Matrix3dRowMajor::Zero();
    Matrix3dRowMajor                               surfaceStressMinus = Matrix3dRowMajor::Zero();
    Matrix3dRowMajor                               Q_ww               = Matrix3dRowMajor::Zero();
    Eigen::Matrix< double, 3, 9, Eigen::RowMajor > Q_wAp  = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >::Zero();
    Eigen::Matrix< double, 3, 9, Eigen::RowMajor > Q_wAm  = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >::Zero();
    Eigen::Matrix< double, 9, 3, Eigen::RowMajor > Q_Apw  = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >::Zero();
    Eigen::Matrix< double, 9, 3, Eigen::RowMajor > Q_Amw  = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >::Zero();
    Matrix9dRowMajor                               Q_ApAp = Matrix9dRowMajor::Zero(), Q_ApAm = Matrix9dRowMajor::Zero();
    Matrix9dRowMajor                               Q_AmAp = Matrix9dRowMajor::Zero(), Q_AmAm = Matrix9dRowMajor::Zero();

    Eigen::Matrix< double, 6, 1 >  dU       = Eigen::Matrix< double, 6, 1 >::Zero();
    Eigen::Matrix< double, 18, 1 > dSurface = Eigen::Matrix< double, 18, 1 >::Zero();
    dU.segment< 3 >( 0 )                    = q.segment< 3 >( WMaterial::offsetJump );
    dSurface.segment< 9 >( 0 )              = q.segment< 9 >( WMaterial::offsetSurfacePlus );
    dSurface.segment< 9 >( 9 )              = q.segment< 9 >( WMaterial::offsetSurfaceMinus );

    GLMaterial::State    state{ force.data(), surfaceStressPlus.data(), surfaceStressMinus.data(), stateVars.data() };
    GLMaterial::Tangents tangents{ Q_ww.data(),
                                   Q_wAp.data(),
                                   Q_wAm.data(),
                                   Q_Apw.data(),
                                   Q_ApAp.data(),
                                   Q_ApAm.data(),
                                   Q_Amw.data(),
                                   Q_AmAp.data(),
                                   Q_AmAm.data() };
    GLMaterial::Deformation   deformation{ dU.data(),
                                         dSurface.data(),
                                         interfaceNormal.data(),
                                         interfaceSeparation.data() };
    GLMaterial::TimeIncrement timeIncrement{ 0.0, 1.0 };
    reference.computeStress( state, tangents, deformation, timeIncrement );

    double stressError = ( warpingStress.segment< 3 >( WMaterial::offsetJump ) - force ).lpNorm< Eigen::Infinity >();
    stressError        = std::max( stressError,
                            ( warpingStress.segment< 9 >( WMaterial::offsetSurfacePlus ) -
                              Eigen::Map< Eigen::Matrix< double, 9, 1 > >( surfaceStressPlus.data() ) )
                              .lpNorm< Eigen::Infinity >() );
    stressError        = std::max( stressError,
                            ( warpingStress.segment< 9 >( WMaterial::offsetSurfaceMinus ) -
                              Eigen::Map< Eigen::Matrix< double, 9, 1 > >( surfaceStressMinus.data() ) )
                              .lpNorm< Eigen::Infinity >() );

    double tangentError = 0.0;
    auto   compare      = [&]( const auto& actual, const auto& expected ) {
      tangentError = std::max( tangentError, ( actual - expected ).template lpNorm< Eigen::Infinity >() );
    };
    compare( warpingTangent.block< 3, 3 >( 0, 0 ), Q_ww );
    compare( warpingTangent.block< 3, 9 >( 0, 3 ), Q_wAp );
    compare( warpingTangent.block< 3, 9 >( 0, 12 ), Q_wAm );
    compare( warpingTangent.block< 9, 3 >( 3, 0 ), Q_Apw );
    compare( warpingTangent.block< 9, 9 >( 3, 3 ), Q_ApAp );
    compare( warpingTangent.block< 9, 9 >( 3, 12 ), Q_ApAm );
    compare( warpingTangent.block< 9, 3 >( 12, 0 ), Q_Amw );
    compare( warpingTangent.block< 9, 9 >( 12, 3 ), Q_AmAp );
    compare( warpingTangent.block< 9, 9 >( 12, 12 ), Q_AmAm );

    std::cout << "reduction to Gauss-Lobatto (" << materialName << "): max|dp| = " << stressError
              << ", max|dK| = " << tangentError << "\n";

    throwExceptionOnFailure( stressError < 1e-14,
                             "with Ws = Wa = 0 the warping material must reproduce the Gauss-Lobatto generalized "
                             "stresses exactly (" +
                               materialName + "): max error " + std::to_string( stressError ) );
    throwExceptionOnFailure( tangentError < 1e-9,
                             "with Ws = Wa = 0 the warping material must reproduce the Gauss-Lobatto tangent "
                             "exactly (" +
                               materialName + "): max error " + std::to_string( tangentError ) );
  }

  void TestReductionToGaussLobattoWithoutWarping()
  {
    checkReductionToGaussLobatto( "LINEARELASTIC", elasticProperties );
    checkReductionToGaussLobatto( "VONMISES", plasticProperties );
  }

  // ------------------------------------------------------------------
  // 3. Full 39x39 condensed tangent against central differences.
  // ------------------------------------------------------------------
  void checkCondensedTangentAgainstFiniteDifference( const std::string&           materialName,
                                                     const std::vector< double >& properties )
  {
    const VectorGen q = makeRepresentativeGeneralizedStrain( true );

    VectorGen stress;
    MatrixGen tangent;
    evaluateWarpingMaterial( materialName, properties, q, stress, tangent );

    MatrixGen    finiteDifference = MatrixGen::Zero();
    const double perturbation     = 1.0e-9;
    for ( int column = 0; column < nGen; ++column ) {
      VectorGen forward = q, backward = q;
      forward( column ) += perturbation;
      backward( column ) -= perturbation;

      VectorGen stressForward, stressBackward;
      MatrixGen unused;
      evaluateWarpingMaterial( materialName, properties, forward, stressForward, unused );
      evaluateWarpingMaterial( materialName, properties, backward, stressBackward, unused );

      finiteDifference.col( column ) = ( stressForward - stressBackward ) / ( 2.0 * perturbation );
    }

    const double error         = ( tangent - finiteDifference ).lpNorm< Eigen::Infinity >();
    const double scale         = std::max( 1.0, finiteDifference.lpNorm< Eigen::Infinity >() );
    const double relativeError = error / scale;

    std::cout << "full 39x39 condensed tangent vs FD (" << materialName << "): relative error = " << relativeError
              << "\n";

    throwExceptionOnFailure( relativeError < 1e-6,
                             "warping material condensed tangent does not match central differences (" + materialName +
                               "): relative error " + std::to_string( relativeError ) );
  }

  void TestCondensedTangentMatchesFiniteDifference()
  {
    checkCondensedTangentAgainstFiniteDifference( "LINEARELASTIC", elasticProperties );
    checkCondensedTangentAgainstFiniteDifference( "VONMISES", plasticProperties );
  }

  // ------------------------------------------------------------------
  // 4. The warping modes must actually carry energy -- and exactly the
  //    directions the design predicts must NOT.
  //
  // With n = e_z the surface gradient has a structurally zero normal column,
  // so only the six tangential components of each slot are meaningful. Of
  // those, the material is provably blind to:
  //   * the in-plane skew part (only sym(.) enters)   -- for A+/A- alike,
  //   * (grad_s w_a)_{z j}, exactly absorbed by the station gradients because
  //     sum_alpha lambda_alpha phi_a(xi_alpha) = 0.
  // Expected ranks of the 6x6 tangential sub-blocks: A+ 5, Ws 5, Wa 3.
  // ------------------------------------------------------------------
  void TestWarpingBlocksCarryEnergyWithTheExpectedRank()
  {
    const VectorGen q = VectorGen::Zero();

    VectorGen stress;
    MatrixGen tangent;
    evaluateWarpingMaterial( "LINEARELASTIC", elasticProperties, q, stress, tangent );

    std::array< int, 6 > tangential{};
    int                  k = 0;
    for ( int i = 0; i < 3; ++i ) {
      for ( int j = 0; j < 2; ++j ) {
        tangential[k++] = 3 * i + j;
      }
    }

    auto numericalRank = [&]( int offset ) {
      Eigen::Matrix< double, 6, 6 > block;
      for ( int a = 0; a < 6; ++a ) {
        for ( int b = 0; b < 6; ++b ) {
          block( a, b ) = tangent( offset + tangential[a], offset + tangential[b] );
        }
      }
      Eigen::JacobiSVD< Eigen::Matrix< double, 6, 6 > > svd( block );
      const double                                      largest = svd.singularValues()( 0 );
      int                                               rank    = 0;
      for ( int i = 0; i < 6; ++i ) {
        if ( svd.singularValues()( i ) > 1.0e-10 * largest ) {
          ++rank;
        }
      }
      return std::make_pair( rank, largest );
    };

    const auto [rankSurfacePlus, scaleSurfacePlus]     = numericalRank( WMaterial::offsetSurfacePlus );
    const auto [rankSymmetric, scaleSymmetric]         = numericalRank( WMaterial::offsetWarpingSymmetric );
    const auto [rankAntisymmetric, scaleAntisymmetric] = numericalRank( WMaterial::offsetWarpingAntisymmetric );

    std::cout << "tangential sub-block ranks: A+ = " << rankSurfacePlus << ", Ws = " << rankSymmetric
              << ", Wa = " << rankAntisymmetric << "\n";

    throwExceptionOnFailure( scaleSymmetric > 0.0 && scaleAntisymmetric > 0.0,
                             "both warping stiffness blocks must be nonzero -- the modes must carry energy." );
    throwExceptionOnFailure( rankSurfacePlus == 5,
                             "unexpected rank of the A+ tangential sub-block: " + std::to_string( rankSurfacePlus ) );
    throwExceptionOnFailure( rankSymmetric == 5,
                             "the symmetric warping block must have the same rank as A+ (only the in-plane skew "
                             "direction is inert), got " +
                               std::to_string( rankSymmetric ) );
    throwExceptionOnFailure( rankAntisymmetric == 3,
                             "the antisymmetric warping block must have rank 3 (in-plane skew plus the two "
                             "normal-row directions absorbed by the station gradients), got " +
                               std::to_string( rankAntisymmetric ) );

    // The two warping modes are orthogonal through the thickness, so their
    // cross-coupling block must vanish identically.
    const double crossCoupling = tangent
                                   .block< 9, 9 >( WMaterial::offsetWarpingSymmetric,
                                                   WMaterial::offsetWarpingAntisymmetric )
                                   .lpNorm< Eigen::Infinity >();
    throwExceptionOnFailure( crossCoupling < 1e-12,
                             "the symmetric and antisymmetric warping modes must not couple: max|K_WsWa| = " +
                               std::to_string( crossCoupling ) );
  }

  // ------------------------------------------------------------------
  // 5. THE DECISIVE TEST: a through-thickness-nonlinear in-plane strain is
  //    invisible to the Gauss-Lobatto material and visible here.
  //
  // Apply a pure symmetric warping in-plane stretch Ws_xx. The Gauss-Lobatto
  // material has no generalized strain that can represent it at all -- its
  // in-plane strain is linear in zeta by construction. The warping material
  // must respond with a nonzero conjugate stress and a positive second
  // variation along that direction.
  // ------------------------------------------------------------------
  void TestParabolicInPlaneWarpingCarriesStiffness()
  {
    VectorGen q                                        = VectorGen::Zero();
    q( WMaterial::offsetWarpingSymmetric + 3 * 0 + 0 ) = 1.0e-4; // Ws_xx

    VectorGen stress;
    MatrixGen tangent;
    evaluateWarpingMaterial( "LINEARELASTIC", elasticProperties, q, stress, tangent );

    const double conjugate = stress.segment< 9 >( WMaterial::offsetWarpingSymmetric ).norm();
    const double stiffness = tangent( WMaterial::offsetWarpingSymmetric, WMaterial::offsetWarpingSymmetric );

    std::cout << "pure Ws_xx: |S_Ws| = " << conjugate << ", diagonal stiffness = " << stiffness << "\n";

    throwExceptionOnFailure( conjugate > 1e-6,
                             "a parabolic in-plane warping stretch must produce a nonzero conjugate stress." );
    throwExceptionOnFailure( stiffness > 0.0, "the parabolic warping direction must have positive stiffness." );

    // The antisymmetric slot must stay silent: the profiles are orthogonal.
    throwExceptionOnFailure( stress.segment< 9 >( WMaterial::offsetWarpingAntisymmetric ).norm() < 1e-14,
                             "a purely symmetric warping strain must not excite the antisymmetric conjugate stress." );

    // And the antisymmetric in-plane counterpart must carry stiffness too.
    VectorGen qAnti                                            = VectorGen::Zero();
    qAnti( WMaterial::offsetWarpingAntisymmetric + 3 * 0 + 0 ) = 1.0e-4; // Wa_xx
    VectorGen stressAnti;
    MatrixGen tangentAnti;
    evaluateWarpingMaterial( "LINEARELASTIC", elasticProperties, qAnti, stressAnti, tangentAnti );

    throwExceptionOnFailure( stressAnti.segment< 9 >( WMaterial::offsetWarpingAntisymmetric ).norm() > 1e-8,
                             "a cubic in-plane warping stretch must produce a nonzero conjugate stress." );
    throwExceptionOnFailure( tangentAnti( WMaterial::offsetWarpingAntisymmetric,
                                          WMaterial::offsetWarpingAntisymmetric ) > 0.0,
                             "the cubic warping direction must have positive stiffness." );
  }

  // ------------------------------------------------------------------
  // 6. The transverse-shear directions of the ANTISYMMETRIC field are exactly
  //    inert -- this is what proves the enrichment is not double-counting the
  //    station normal gradients. Any nonzero response here would mean the
  //    formulation had become redundant (and the element block singular for a
  //    reason other than the documented one).
  // ------------------------------------------------------------------
  void TestAntisymmetricTransverseShearIsAbsorbedByStationGradients()
  {
    for ( int tangentialDirection = 0; tangentialDirection < 2; ++tangentialDirection ) {
      VectorGen q = VectorGen::Zero();
      // (grad_s w_a)_{z j}, j tangential
      q( WMaterial::offsetWarpingAntisymmetric + 3 * 2 + tangentialDirection ) = 1.0e-4;

      VectorGen stress;
      MatrixGen tangent;
      evaluateWarpingMaterial( "LINEARELASTIC", elasticProperties, q, stress, tangent );

      const double response = stress.lpNorm< Eigen::Infinity >();
      throwExceptionOnFailure( response < 1e-12,
                               "the antisymmetric warping's transverse-shear direction must be exactly absorbed by "
                               "the station normal gradients (response " +
                                 std::to_string( response ) + " for direction " +
                                 std::to_string( tangentialDirection ) + ")." );
    }
  }

} // namespace

int main()
{
  auto tests = std::vector< std::function< void() > >{ TestWarpingProfileQuadratureIdentities,
                                                       TestReductionToGaussLobattoWithoutWarping,
                                                       TestCondensedTangentMatchesFiniteDifference,
                                                       TestWarpingBlocksCarryEnergyWithTheExpectedRank,
                                                       TestParabolicInPlaneWarpingCarriesStiffness,
                                                       TestAntisymmetricTransverseShearIsAbsorbedByStationGradients };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
