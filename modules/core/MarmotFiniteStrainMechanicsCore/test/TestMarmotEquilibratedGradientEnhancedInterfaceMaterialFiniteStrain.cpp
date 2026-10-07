#include "Marmot/MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain.h"
#include "Marmot/MarmotEquilibratedInterfaceMaterialFiniteStrain.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace Marmot::Testing;

namespace {

  using Material  = MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain;
  using Vector21d = Eigen::Matrix< double, 21, 1 >;
  using Matrix21d = Eigen::Matrix< double, 21, 21 >;

  constexpr double pi = 3.14159265358979323846;

  // bulk modulus, shear modulus, thickness, damage threshold, softening parameter, nonlocal radius, density
  const std::vector< double > card = { 3500.0, 1500.0, 0.0985, 1.0e-3, 5.0e-3, 0.1, 1.0 };

  const Eigen::Vector3d normal( -std::sin( 10.0 * pi / 180.0 ), 0.0, std::cos( 10.0 * pi / 180.0 ) );
  const Eigen::Vector3d connector( 0.0, 0.0, 0.1 );

  struct Response {
    Vector21d p;
    Matrix21d tangent;
    Vector21d dStress_dNonlocal;
    Vector21d dLocalField_dStrain;
    double    dLocalField_dNonlocal;
    double    localField;
    double    nonlocalRadius;
  };

  /// A generalized strain with two different faces, so that the normal-gradient jump is nonzero.
  Vector21d generalizedStrain()
  {
    Vector21d x = Vector21d::Zero();
    x( 0 )      = 0.004;
    x( 1 )      = 0.002;
    x( 2 )      = 0.003;
    x( 3 )      = 0.021;
    x( 7 )      = 0.014;
    x( 4 )      = 0.009;
    x( 12 )     = -0.017;
    x( 16 )     = -0.011;
    x( 13 )     = 0.006;
    return x;
  }

  /// One stress update from a freshly initialized state, in which the top face may already be damaged up to the
  /// history value topKappa of its bulk material.
  Response evaluate( const Vector21d& x, double nonlocalField, double topKappa = 0.0 )
  {
    Material material( "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE",
                       card.data(),
                       static_cast< int >( card.size() ),
                       1 );

    std::vector< double > stateVars( material.getNumberOfRequiredStateVars(), 0.0 );
    material.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );
    material.getStateView( "topMaterialStateVars", stateVars.data() ).stateLocation[0] = topKappa;

    std::array< double, 6 >  U{};
    std::array< double, 18 > surfaceGradient{};
    for ( int i = 0; i < 3; ++i )
      U[i] = x( i );
    for ( int i = 0; i < 9; ++i ) {
      surfaceGradient[i]     = x( 3 + i );
      surfaceGradient[9 + i] = x( 12 + i );
    }

    Material::State       state{ { Material::Tensor3d( 0.0 ),
                                   Material::Tensor33d( 0.0 ),
                                   Material::Tensor33d( 0.0 ),
                                   stateVars.data() },
                           0.0,
                           0.0 };
    Material::Tangents    tangents{ { Material::Tensor33d( 0.0 ),
                                      Material::Tensor333d( 0.0 ),
                                      Material::Tensor333d( 0.0 ),
                                      Material::Tensor333d( 0.0 ),
                                      Material::Tensor3333d( 0.0 ),
                                      Material::Tensor3333d( 0.0 ),
                                      Material::Tensor333d( 0.0 ),
                                      Material::Tensor3333d( 0.0 ),
                                      Material::Tensor3333d( 0.0 ) },
                                 Material::Tensor21d( 0.0 ),
                                 Material::Tensor21d( 0.0 ),
                                 0.0 };
    Material::Deformation deformation{ { Material::Tensor6d( U.data() ),
                                         Material::Tensor18d( surfaceGradient.data() ),
                                         Material::Tensor3d( normal.data() ),
                                         Material::Tensor3d( connector.data() ) },
                                       nonlocalField };

    material.computeStress( state, tangents, deformation, { 0.0, 1.0 } );

    Response response;
    response.p.segment< 3 >( 0 )  = Eigen::Map< const Eigen::Vector3d >( state.generalizedForce.data() );
    response.p.segment< 9 >( 3 )  = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( state.surfaceStressPlus.data() );
    response.p.segment< 9 >( 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >(
      state.surfaceStressMinus.data() );

    using RowMajor = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >;
    response.tangent.setZero();
    response.tangent.block< 3, 3 >( 0, 0 )   = Eigen::Map< const RowMajor >( tangents.Q_ww.data(), 3, 3 );
    response.tangent.block< 3, 9 >( 0, 3 )   = Eigen::Map< const RowMajor >( tangents.Q_wAp.data(), 3, 9 );
    response.tangent.block< 3, 9 >( 0, 12 )  = Eigen::Map< const RowMajor >( tangents.Q_wAm.data(), 3, 9 );
    response.tangent.block< 9, 3 >( 3, 0 )   = Eigen::Map< const RowMajor >( tangents.Q_Apw.data(), 9, 3 );
    response.tangent.block< 9, 9 >( 3, 3 )   = Eigen::Map< const RowMajor >( tangents.Q_ApAp.data(), 9, 9 );
    response.tangent.block< 9, 9 >( 3, 12 )  = Eigen::Map< const RowMajor >( tangents.Q_ApAm.data(), 9, 9 );
    response.tangent.block< 9, 3 >( 12, 0 )  = Eigen::Map< const RowMajor >( tangents.Q_Amw.data(), 9, 3 );
    response.tangent.block< 9, 9 >( 12, 3 )  = Eigen::Map< const RowMajor >( tangents.Q_AmAp.data(), 9, 9 );
    response.tangent.block< 9, 9 >( 12, 12 ) = Eigen::Map< const RowMajor >( tangents.Q_AmAm.data(), 9, 9 );

    response.dStress_dNonlocal     = Eigen::Map< const Vector21d >( tangents.dGeneralizedStress_dNonlocal.data() );
    response.dLocalField_dStrain   = Eigen::Map< const Vector21d >( tangents.dLocalField_dGeneralizedStrain.data() );
    response.dLocalField_dNonlocal = tangents.dLocalField_dNonlocal;
    response.localField            = state.localField;
    response.nonlocalRadius        = state.nonlocalRadius;
    return response;
  }

  /// Below the damage threshold the stress does not depend on the nonlocal field, and the response is the one of the
  /// neo-Hookean material.
  void testMatchesTheLocalMaterialBelowTheThreshold()
  {
    const Vector21d x = generalizedStrain();

    const auto enhanced = evaluate( x, 0.0 );

    // the local material with the same bulk and shear modulus and thickness
    const std::vector< double >                     localCard = { 3500.0, 1500.0, 0.0985 };
    MarmotEquilibratedInterfaceMaterialFiniteStrain local( "COMPRESSIBLENEOHOOKE", localCard.data(), 3, 1 );

    std::vector< double > stateVars( local.getNumberOfRequiredStateVars(), 0.0 );
    local.initializeYourself( stateVars.data(), static_cast< int >( stateVars.size() ) );

    std::array< double, 6 >  U{};
    std::array< double, 18 > surfaceGradient{};
    for ( int i = 0; i < 3; ++i )
      U[i] = x( i );
    for ( int i = 0; i < 9; ++i ) {
      surfaceGradient[i]     = x( 3 + i );
      surfaceGradient[9 + i] = x( 12 + i );
    }

    using L = MarmotEquilibratedInterfaceMaterialFiniteStrain;
    L::State       state{ L::Tensor3d( 0.0 ), L::Tensor33d( 0.0 ), L::Tensor33d( 0.0 ), stateVars.data() };
    L::Tangents    tangents{ L::Tensor33d( 0.0 ),
                          L::Tensor333d( 0.0 ),
                          L::Tensor333d( 0.0 ),
                          L::Tensor333d( 0.0 ),
                          L::Tensor3333d( 0.0 ),
                          L::Tensor3333d( 0.0 ),
                          L::Tensor333d( 0.0 ),
                          L::Tensor3333d( 0.0 ),
                          L::Tensor3333d( 0.0 ) };
    L::Deformation deformation{ L::Tensor6d( U.data() ),
                                L::Tensor18d( surfaceGradient.data() ),
                                L::Tensor3d( normal.data() ),
                                L::Tensor3d( connector.data() ) };
    local.computeStress( state, tangents, deformation, { 0.0, 1.0 } );

    Vector21d pLocal;
    pLocal.segment< 3 >( 0 )  = Eigen::Map< const Eigen::Vector3d >( state.generalizedForce.data() );
    pLocal.segment< 9 >( 3 )  = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( state.surfaceStressPlus.data() );
    pLocal.segment< 9 >( 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( state.surfaceStressMinus.data() );

    throwExceptionOnFailure( pLocal.norm() > 0.0, "The state should load the layer." );
    throwExceptionOnFailure( ( pLocal - enhanced.p ).norm() < 1e-12 * pLocal.norm(),
                             "Below the damage threshold the stress is the one of the local material." );
    throwExceptionOnFailure( enhanced.dStress_dNonlocal.norm() == 0.0,
                             "Below the damage threshold the stress does not depend on the nonlocal field." );
    throwExceptionOnFailure( enhanced.localField > 0.0 && enhanced.nonlocalRadius == 0.1,
                             "The local field and the nonlocal radius are reported." );
  }

  /// The condensed derivatives with respect to the generalized strain and the nonlocal field are the derivatives of
  /// the generalized stress and the local field, with the damage active.
  void checkTangentsAreConsistent( double topKappa, const std::string& context )
  {
    const Vector21d x             = generalizedStrain();
    const double    nonlocalField = 3.0e-3;
    const auto      base          = evaluate( x, nonlocalField, topKappa );

    throwExceptionOnFailure( base.dStress_dNonlocal.norm() > 0.0,
                             context + ": the damage should couple the stress to the field." );

    const double eps = 1e-7;

    Matrix21d finiteDifference;
    Vector21d localFieldGradient;
    for ( int j = 0; j < 21; ++j ) {
      Vector21d plus = x, minus = x;
      plus( j ) += eps;
      minus( j ) -= eps;

      const auto rp = evaluate( plus, nonlocalField, topKappa );
      const auto rm = evaluate( minus, nonlocalField, topKappa );

      finiteDifference.col( j ) = ( rp.p - rm.p ) / ( 2.0 * eps );
      localFieldGradient( j )   = ( rp.localField - rm.localField ) / ( 2.0 * eps );
    }

    const double tangentError = ( finiteDifference - base.tangent ).norm() / std::max( 1.0, finiteDifference.norm() );
    throwExceptionOnFailure( tangentError < 1e-6,
                             context + ": the condensed tangent is not the derivative of the generalized stress, " +
                               std::to_string( tangentError ) );

    const double gradientError = ( localFieldGradient - base.dLocalField_dStrain ).norm() /
                                 std::max( 1.0, localFieldGradient.norm() );
    throwExceptionOnFailure( gradientError < 1e-6,
                             context + ": the derivative of the local field is not consistent, " +
                               std::to_string( gradientError ) );

    const double epsN = 1e-8;
    const auto   np   = evaluate( x, nonlocalField + epsN, topKappa );
    const auto   nm   = evaluate( x, nonlocalField - epsN, topKappa );

    const Vector21d stressDerivative = ( np.p - nm.p ) / ( 2.0 * epsN );
    const double    stressError      = ( stressDerivative - base.dStress_dNonlocal ).norm() /
                               std::max( 1.0, stressDerivative.norm() );
    throwExceptionOnFailure( stressError < 1e-6,
                             context +
                               ": the derivative of the stress with respect to the nonlocal field is not consistent, " +
                               std::to_string( stressError ) );

    const double localDerivative = ( np.localField - nm.localField ) / ( 2.0 * epsN );
    throwExceptionOnFailure( std::abs( localDerivative - base.dLocalField_dNonlocal ) <
                               1e-6 * std::max( 1.0, std::abs( localDerivative ) ),
                             "The derivative of the local field with respect to the nonlocal field is not "
                             "consistent." );
  }

  /**
   * With the same history on both faces, the traction that the nonlocal field changes is the same on both of them, and
   * the path of the nonlocal field through z vanishes. With a different history it does not, and the condensation of
   * the derivatives with respect to the nonlocal field matters.
   */
  void testTangentsAreConsistent()
  {
    checkTangentsAreConsistent( 0.0, "same history" );

    // the top face is damaged already, and unloads while the bottom face loads
    checkTangentsAreConsistent( 5.0e-3, "different histories" );
  }

  void testInvalidCardsAreRejected()
  {
    const auto rejects = []( const std::string& name, const std::vector< double >& values ) {
      try {
        Material material( name, values.data(), static_cast< int >( values.size() ), 1 );
      }
      catch ( const std::invalid_argument& ) {
        return true;
      }
      return false;
    };

    throwExceptionOnFailure( rejects( "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE", { 3500.0, 1500.0 } ),
                             "A card without the thickness is rejected." );
    throwExceptionOnFailure( rejects( "GRADIENTENHANCEDCOMPRESSIBLENEOHOOKEDAMAGE",
                                      { 3500.0, 1500.0, -1.0, 1e-3, 5e-3, 0.1, 1.0 } ),
                             "A negative thickness is rejected." );
    throwExceptionOnFailure( rejects( "NOSUCHMATERIAL", { 3500.0, 1500.0, 0.1 } ), "An unknown material is rejected." );
  }

} // namespace

int main()
{
  std::vector< std::function< void() > > tests = { testMatchesTheLocalMaterialBelowTheThreshold,
                                                   testTangentsAreConsistent,
                                                   testInvalidCardsAreRejected };
  executeTestsAndCollectExceptions( tests );
  return 0;
}
