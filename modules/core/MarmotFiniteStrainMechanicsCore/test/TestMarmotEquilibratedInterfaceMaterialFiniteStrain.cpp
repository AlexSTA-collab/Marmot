#include "Marmot/MarmotEquilibratedInterfaceMaterialFiniteStrain.h"
#include "Marmot/MarmotTesting.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;

namespace {

  constexpr double pi = 3.14159265358979323846;

  using Vector3d         = Eigen::Matrix< double, 3, 1 >;
  using Vector21d        = Eigen::Matrix< double, 21, 1 >;
  using Matrix3dRowMajor = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix21d        = Eigen::Matrix< double, 21, 21 >;

  constexpr double h  = 0.1;
  constexpr double E  = 210000.0;
  constexpr double nu = 0.3;

  const double K = E / ( 3.0 * ( 1.0 - 2.0 * nu ) );
  const double G = E / ( 2.0 * ( 1.0 + nu ) );

  /** The geometry of the study meshes: an interface inclined by alpha, its faces separated by
   * the nominal thickness offset VERTICALLY, so the connector has a tangential part. */
  struct Geometry {
    Vector3d normal;
    Vector3d connector;
    double   hNormal;
  };

  Geometry inclinedGeometry( double alpha )
  {
    Geometry geometry;
    geometry.normal    = Vector3d( -std::sin( alpha ), 0.0, std::cos( alpha ) );
    geometry.connector = Vector3d( 0.0, 0.0, h );
    geometry.hNormal   = geometry.connector.dot( geometry.normal );
    return geometry;
  }

  struct Response {
    Vector21d             p;
    Matrix21d             tangent;
    std::vector< double > stateVars;
    Matrix3dRowMajor      PPlus;
    Matrix3dRowMajor      PMinus;
    Vector3d              z;
  };

  /** One stress update from a freshly initialised state, so that repeated calls -- as a finite
   * difference makes -- all start from the same committed state and re-solve z from scratch. */
  Response evaluate( const std::string&    materialName,
                     std::vector< double > materialProperties,
                     const Geometry&       geometry,
                     const Vector21d&      x,
                     double                dT = 1.0 )
  {
    MarmotEquilibratedInterfaceMaterialFiniteStrain material( materialName,
                                                              materialProperties.data(),
                                                              static_cast< int >( materialProperties.size() ),
                                                              1 );

    Response response;
    response.stateVars.assign( material.getNumberOfRequiredStateVars(), 0.0 );
    material.initializeYourself( response.stateVars.data(), static_cast< int >( response.stateVars.size() ) );

    // the material takes the two faces' TOTAL displacements, so w = x(0..2) is realised as
    // u+ = w and u- = 0
    std::array< double, 6 >  U{};
    std::array< double, 18 > surfaceGradient{};

    for ( int i = 0; i < 3; ++i )
      U[i] = x( i );

    for ( int i = 0; i < 9; ++i ) {
      surfaceGradient[i]     = x( 3 + i );
      surfaceGradient[9 + i] = x( 12 + i );
    }

    using Material = MarmotEquilibratedInterfaceMaterialFiniteStrain;

    Material::State         state{ Material::Tensor3d( 0.0 ),
                           Material::Tensor33d( 0.0 ),
                           Material::Tensor33d( 0.0 ),
                           response.stateVars.data() };
    Material::Tangents      tangents{ Material::Tensor33d( 0.0 ),
                                 Material::Tensor333d( 0.0 ),
                                 Material::Tensor333d( 0.0 ),
                                 Material::Tensor333d( 0.0 ),
                                 Material::Tensor3333d( 0.0 ),
                                 Material::Tensor3333d( 0.0 ),
                                 Material::Tensor333d( 0.0 ),
                                 Material::Tensor3333d( 0.0 ),
                                 Material::Tensor3333d( 0.0 ) };
    Material::Deformation   deformation{ Material::Tensor6d( U.data() ),
                                       Material::Tensor18d( surfaceGradient.data() ),
                                       Material::Tensor3d( geometry.normal.data() ),
                                       Material::Tensor3d( geometry.connector.data() ) };
    Material::TimeIncrement timeIncrement{ 0.0, dT };

    material.computeStress( state, tangents, deformation, timeIncrement );

    response.p.segment< 3 >( 0 )  = Eigen::Map< const Vector3d >( state.generalizedForce.data() );
    response.p.segment< 9 >( 3 )  = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( state.surfaceStressPlus.data() );
    response.p.segment< 9 >( 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >(
      state.surfaceStressMinus.data() );

    response.tangent.setZero();
    response.tangent.block< 3, 3 >( 0, 0 ) = Eigen::Map< const Matrix3dRowMajor >( tangents.Q_ww.data() );
    response.tangent.block< 3, 9 >( 0, 3 ) = Eigen::Map< const Eigen::Matrix< double, 3, 9, Eigen::RowMajor > >(
      tangents.Q_wAp.data() );
    response.tangent.block< 3, 9 >( 0, 12 ) = Eigen::Map< const Eigen::Matrix< double, 3, 9, Eigen::RowMajor > >(
      tangents.Q_wAm.data() );
    response.tangent.block< 9, 3 >( 3, 0 ) = Eigen::Map< const Eigen::Matrix< double, 9, 3, Eigen::RowMajor > >(
      tangents.Q_Apw.data() );
    response.tangent.block< 9, 9 >( 3, 3 ) = Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >(
      tangents.Q_ApAp.data() );
    response.tangent.block< 9, 9 >( 3, 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >(
      tangents.Q_ApAm.data() );
    response.tangent.block< 9, 3 >( 12, 0 ) = Eigen::Map< const Eigen::Matrix< double, 9, 3, Eigen::RowMajor > >(
      tangents.Q_Amw.data() );
    response.tangent.block< 9, 9 >( 12, 3 ) = Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >(
      tangents.Q_AmAp.data() );
    response.tangent.block< 9, 9 >( 12, 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >(
      tangents.Q_AmAm.data() );

    const double* zStore = response.stateVars.data();
    response.z           = Eigen::Map< const Vector3d >( zStore );
    response.PPlus       = Eigen::Map< const Matrix3dRowMajor >( zStore + 3 );
    response.PMinus      = Eigen::Map< const Matrix3dRowMajor >( zStore + 12 );

    return response;
  }

  /** The generalized strain of a finite rigid rotation of the whole layer. */
  Vector21d rigidRotationStrain( const Geometry& geometry, const Vector3d& axis, double theta )
  {
    const Matrix3dRowMajor R = Eigen::AngleAxisd( theta, axis.normalized() ).toRotationMatrix();
    const Matrix3dRowMajor T = Matrix3dRowMajor::Identity() - geometry.normal * geometry.normal.transpose();

    const Matrix3dRowMajor A = ( R - Matrix3dRowMajor::Identity() ) * T;
    const Vector3d         w = ( R - Matrix3dRowMajor::Identity() ) * geometry.connector;

    Vector21d x;
    x.segment< 3 >( 0 ) = w;

    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j ) {
        x( 3 + 3 * i + j )  = A( i, j );
        x( 12 + 3 * i + j ) = A( i, j );
      }

    return x;
  }

} // namespace

void TestEquilibratedFiniteStrainInterfaceIsObjective()
{
  std::cout << "\n--- TestEquilibratedFiniteStrainInterfaceIsObjective ---\n";

  for ( const double alphaDegrees : { 0.0, 10.0 } ) {
    const double   alpha    = alphaDegrees * pi / 180.0;
    const Geometry geometry = inclinedGeometry( alpha );

    const std::vector< double > materialProperties = { K, G, geometry.hNormal };

    std::cout << "  alpha = " << alphaDegrees << " deg, tangential offset / h_normal = " << std::tan( alpha ) << "\n";

    for ( const double thetaDegrees : { 1.0, 15.0, 40.0 } ) {
      const auto x        = rigidRotationStrain( geometry, Vector3d( 0.0, 1.0, 0.0 ), thetaDegrees * pi / 180.0 );
      const auto response = evaluate( "COMPRESSIBLENEOHOOKE", materialProperties, geometry, x );

      const double stress = response.p.lpNorm< Eigen::Infinity >() / ( E * h );
      const double zNorm  = response.z.norm();

      std::cout << "      rotation " << thetaDegrees << " deg -> |p| / (E h) = " << stress << ", |z| = " << zNorm
                << "\n";

      throwExceptionOnFailure( stress < 1e-10,
                               "The equilibrated finite-strain interface carries stress through a rigid rotation." );
      throwExceptionOnFailure( zNorm < 1e-10,
                               "The equilibrated finite-strain interface opens a normal-gradient jump under a rigid "
                               "rotation." );
    }
  }
}

void TestEquilibratedFiniteStrainInterfaceTangentIsConsistent()
{
  std::cout << "\n--- TestEquilibratedFiniteStrainInterfaceTangentIsConsistent ---\n";

  const Geometry              geometry           = inclinedGeometry( 10.0 * pi / 180.0 );
  const std::vector< double > materialProperties = { K, G, geometry.hNormal };

  // a state well away from the reference configuration: a 20 degree rotation of the layer, a
  // sliding and opening jump, and two DIFFERENT face gradients, so that z is genuinely nonzero
  Vector21d x = rigidRotationStrain( geometry, Vector3d( 0.0, 1.0, 0.0 ), 20.0 * pi / 180.0 );

  x( 0 ) += 0.004;
  x( 1 ) += 0.002;
  x( 2 ) += 0.003;

  x( 3 ) += 0.021;  // A+_xx
  x( 7 ) += 0.014;  // A+_yy
  x( 4 ) += 0.009;  // A+_xy
  x( 12 ) -= 0.017; // A-_xx
  x( 16 ) -= 0.011; // A-_yy
  x( 13 ) += 0.006; // A-_xy

  const auto base = evaluate( "COMPRESSIBLENEOHOOKE", materialProperties, geometry, x );

  std::cout << "  |z| at the solved state = " << base.z.norm() << "\n";

  throwExceptionOnFailure( base.z.norm() > 1e-6, "The chosen state should open a normal-gradient jump." );

  Matrix21d finiteDifference;
  finiteDifference.setZero();

  const double eps = 1e-7;

  for ( int j = 0; j < 21; ++j ) {
    Vector21d plus  = x;
    Vector21d minus = x;

    plus( j ) += eps;
    minus( j ) -= eps;

    finiteDifference.col( j ) = ( evaluate( "COMPRESSIBLENEOHOOKE", materialProperties, geometry, plus ).p -
                                  evaluate( "COMPRESSIBLENEOHOOKE", materialProperties, geometry, minus ).p ) /
                                ( 2.0 * eps );
  }

  const double error = ( base.tangent - finiteDifference ).lpNorm< Eigen::Infinity >();
  const double scale = std::max( 1.0, finiteDifference.lpNorm< Eigen::Infinity >() );

  std::cout << "  relative tangent error against central differences = " << error / scale << "\n";

  throwExceptionOnFailure( error / scale < 1e-6,
                           "The condensed finite-strain tangent is not the derivative of the generalized stress." );
}

void TestEquilibratedFiniteStrainInterfaceEquilibratesTheNominalTraction()
{
  std::cout << "\n--- TestEquilibratedFiniteStrainInterfaceEquilibratesTheNominalTraction ---\n";

  const Geometry              geometry           = inclinedGeometry( 10.0 * pi / 180.0 );
  const std::vector< double > materialProperties = { K, G, geometry.hNormal };

  Vector21d x;
  x.setZero();

  x( 0 )  = 0.002;
  x( 2 )  = 0.004;
  x( 3 )  = 0.030; // the two faces are stretched differently, which is what drives z
  x( 12 ) = -0.020;

  const auto response = evaluate( "COMPRESSIBLENEOHOOKE", materialProperties, geometry, x );

  const Vector3d tPlus  = response.PPlus * geometry.normal;
  const Vector3d tMinus = response.PMinus * geometry.normal;

  const double imbalance = ( tPlus - tMinus ).norm() / std::max( 1.0, tPlus.norm() );

  std::cout << "  |z| = " << response.z.norm() << ", |P+ N - P- N| / |P+ N| = " << imbalance << "\n";

  throwExceptionOnFailure( response.z.norm() > 1e-6, "The chosen state should open a normal-gradient jump." );
  throwExceptionOnFailure( imbalance < 1e-9, "The local solve does not equilibrate the nominal traction." );
}

void TestEquilibratedFiniteStrainInterfaceCarriesARateMaterial()
{
  std::cout << "\n--- TestEquilibratedFiniteStrainInterfaceCarriesARateMaterial ---\n";

  const Geometry geometry = inclinedGeometry( 10.0 * pi / 180.0 );

  // baseModel = NeoHooke, onlyShearCreep = off, K, G, one Maxwell branch -- the property
  // list the material's own test uses, with this study's moduli
  const std::vector< double > materialProperties = { 0.0, 0.0, geometry.hNormal, K, G, 1.0, 0.3, 10.0 };

  Vector21d x;
  x.setZero();
  x( 0 )  = 0.002;
  x( 2 )  = 0.004;
  x( 3 )  = 0.010;
  x( 12 ) = -0.008;

  Vector3d previous = Vector3d::Zero();

  for ( const double dT : { 0.1, 1.0, 10.0, 100.0 } ) {
    const auto response = evaluate( "COMPRESSIBLEFINITESTRAINLINEARVISCOELASTICITY",
                                    materialProperties,
                                    geometry,
                                    x,
                                    dT );

    const Vector3d force = response.p.segment< 3 >( 0 );

    std::cout << "  dT = " << dT << " -> |f| = " << force.norm() << ", |z| = " << response.z.norm() << "\n";

    throwExceptionOnFailure( force.allFinite(), "The viscoelastic response is not finite." );

    // relaxation: a longer step over the same total deformation must carry less force
    if ( previous.norm() > 0.0 )
      throwExceptionOnFailure( force.norm() < previous.norm() * ( 1.0 + 1e-12 ),
                               "The viscoelastic interface does not relax as the step grows." );

    previous = force;
  }
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    TestEquilibratedFiniteStrainInterfaceIsObjective,
    TestEquilibratedFiniteStrainInterfaceTangentIsConsistent,
    TestEquilibratedFiniteStrainInterfaceEquilibratesTheNominalTraction,
    TestEquilibratedFiniteStrainInterfaceCarriesARateMaterial,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
