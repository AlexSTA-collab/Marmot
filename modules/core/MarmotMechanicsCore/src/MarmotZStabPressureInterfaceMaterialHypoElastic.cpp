/* ---------------------------------------------------------------------
 *  marmot - MAteRialMOdellingToolbox
 *  Alexandros Stathas alexandros.stathas@boku.ac.at
 *  LGPL 2.1 or later; see LICENSE.md at the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#include "Marmot/MarmotZStabPressureInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot;

namespace {

  using Vector3d           = Eigen::Matrix< double, 3, 1 >;
  using Vector9d           = Eigen::Matrix< double, 9, 1 >;
  using Vector21d          = Eigen::Matrix< double, 21, 1 >;
  using Matrix3dRowMajor   = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor   = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor  = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;
  using Matrix3x9RowMajor  = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;
  using Matrix9x21RowMajor = Eigen::Matrix< double, 9, 21, Eigen::RowMajor >;
  using Matrix21x3RowMajor = Eigen::Matrix< double, 21, 3, Eigen::RowMajor >;
  using Matrix3x21RowMajor = Eigen::Matrix< double, 3, 21, Eigen::RowMajor >;
  using Matrix21dRowMajor  = Eigen::Matrix< double, 21, 21, Eigen::RowMajor >;

  constexpr int flatIndex( int i, int j )
  {
    return 3 * i + j;
  }

  Matrix9dRowMajor fullGradientTangent( const Matrix6d& tangentVoigt )
  {
    Matrix9dRowMajor T = Matrix9dRowMajor::Zero();
    for ( int k = 0; k < 3; ++k )
      for ( int l = 0; l < 3; ++l ) {
        Matrix3dRowMajor dG = Matrix3dRowMajor::Zero();
        dG( k, l )          = 1.0;
        const Vector6d dEps = ContinuumMechanics::VoigtNotation::strainToVoigt(
          Matrix3dRowMajor( 0.5 * ( dG + dG.transpose() ) ) );
        const Eigen::Matrix3d dS = ContinuumMechanics::VoigtNotation::voigtToStress( Vector6d( tangentVoigt * dEps ) );
        for ( int i = 0; i < 3; ++i )
          for ( int j = 0; j < 3; ++j )
            T( flatIndex( i, j ), flatIndex( k, l ) ) = dS( i, j );
      }
    return T;
  }

  /** 9x9 deviatoric projector on vec(sigma): dev = sigma - (tr sigma / 3) I. */
  Matrix9dRowMajor deviatoricProjector()
  {
    Matrix9dRowMajor P = Matrix9dRowMajor::Identity();
    for ( int i = 0; i < 3; ++i )
      for ( int k = 0; k < 3; ++k )
        P( flatIndex( i, i ), flatIndex( k, k ) ) -= 1.0 / 3.0;
    return P;
  }

  struct InterfaceGeometry {
    double   normalSeparation;
    Vector3d tangentialSeparation;
  };

  InterfaceGeometry evaluateInterfaceGeometry( const Vector3d& n, const Vector3d& d, double h )
  {
    constexpr double tol = 1.0e-12;
    if ( h <= 0.0 )
      throw std::invalid_argument( "MarmotZStabPressureInterfaceMaterialHypoElastic: h must be positive." );
    if ( d.norm() <= tol )
      return { h, Vector3d::Zero() };
    const double ell = d.dot( n );
    if ( ell <= tol )
      throw std::invalid_argument(
        "MarmotZStabPressureInterfaceMaterialHypoElastic: connector must have a positive normal component." );
    return { ell, d - ell * n };
  }

  Matrix9x3RowMajor jumpToGradient( const Vector3d& n, double ell )
  {
    Matrix9x3RowMajor m = Matrix9x3RowMajor::Zero();
    for ( int k = 0; k < 3; ++k )
      for ( int l = 0; l < 3; ++l )
        m( flatIndex( k, l ), k ) = n( l ) / ell;
    return m;
  }

  Matrix9dRowMajor surfaceToGradient( const Vector3d& n, const Vector3d& dTau, double ell )
  {
    Matrix9dRowMajor m = Matrix9dRowMajor::Zero();
    for ( int k = 0; k < 3; ++k )
      for ( int l = 0; l < 3; ++l )
        for ( int b = 0; b < 3; ++b )
          m( flatIndex( k, l ), flatIndex( k, b ) ) = ( l == b ? 1.0 : 0.0 ) - dTau( b ) * n( l ) / ell;
    return m;
  }

  struct SideTrial {
    Matrix3dRowMajor      stress;
    Matrix9dRowMajor      CFull;
    std::vector< double > trialStateVars;
  };

  SideTrial evaluateSideTrial( MarmotMaterialHypoElastic&                 material,
                               const double*                              oldSv,
                               int                                        nSv,
                               const Matrix3dRowMajor&                    stressCurrent,
                               const Vector6d&                            dEpsVoigt,
                               const MarmotMaterialHypoElastic::timeInfo& ti )
  {
    SideTrial t;
    t.trialStateVars.assign( oldSv, oldSv + nSv );
    const Vector6d sV = ContinuumMechanics::VoigtNotation::stressToVoigt(
      Eigen::Matrix3d( 0.5 * ( stressCurrent + stressCurrent.transpose() ) ) );
    Matrix6d                           CV = Matrix6d::Zero();
    MarmotMaterialHypoElastic::state3D s{ sV, 0.0, 0.0, t.trialStateVars.data() };
    material.computeStress( s, CV, dEpsVoigt, ti );
    t.stress = ContinuumMechanics::VoigtNotation::voigtToStress( s.stress );
    t.CFull  = fullGradientTangent( CV );
    return t;
  }

  double defaultRegularization()
  {
    static const double v = []() {
      const char* e = std::getenv( "MARMOT_ZIFACE_REG" );
      return e ? std::atof( e ) : 1.0e-6;
    }();
    return v;
  }

} // namespace

MarmotZStabPressureInterfaceMaterialHypoElastic::MarmotZStabPressureInterfaceMaterialHypoElastic(
  const std::string& materialName,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : materialProperties( matProperties_ ),
    nMaterialProperties( nMaterialProperties_ ),
    regularization( defaultRegularization() ),
    materialNumber( materialNumber_ )
{
  if ( nMaterialProperties < 3 )
    throw std::invalid_argument( "MarmotZStabPressureInterfaceMaterialHypoElastic requires at least E, nu and h." );

  h = materialProperties[2];
  if ( h <= 0.0 )
    throw std::invalid_argument( "MarmotZStabPressureInterfaceMaterialHypoElastic requires h > 0." );

  const double E  = materialProperties[0];
  const double nu = std::min( std::max( materialProperties[1], -0.999999 ), 0.499999 );
  shearModulus    = E / ( 2.0 * ( 1.0 + nu ) );
  lameParameter   = E * nu / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );
  bulkModulus     = E / ( 3.0 * ( 1.0 - 2.0 * nu ) );

  baseMaterialProperties.reserve( nMaterialProperties - 1 );
  baseMaterialProperties.push_back( materialProperties[0] );
  baseMaterialProperties.push_back( materialProperties[1] );
  baseMaterialProperties.insert( baseMaterialProperties.end(),
                                 materialProperties + 3,
                                 materialProperties + nMaterialProperties );

  topMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     baseMaterialProperties.data(),
                                                                     (int)baseMaterialProperties.size(),
                                                                     materialNumber ) );
  bottomMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     baseMaterialProperties.data(),
                                                                     (int)baseMaterialProperties.size(),
                                                                     materialNumber ) );

  stateLayout.add( "topMaterialStateVars", topMaterial->getNumberOfRequiredStateVars() );
  stateLayout.add( "bottomMaterialStateVars", bottomMaterial->getNumberOfRequiredStateVars() );
  stateLayout.add( "stressPlus", 6 );
  stateLayout.add( "stressMinus", 6 );
  stateLayout.finalize();
}

void MarmotZStabPressureInterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  if ( topMaterial )
    topMaterial->setCharacteristicElementLength( length );
  if ( bottomMaterial )
    bottomMaterial->setCharacteristicElementLength( length );
}

void MarmotZStabPressureInterfaceMaterialHypoElastic::setGradientJumpRegularization( double zeta )
{
  if ( zeta < 0.0 )
    throw std::invalid_argument( "MarmotZStabPressureInterfaceMaterialHypoElastic: zeta must be >= 0." );
  regularization = zeta;
}

void MarmotZStabPressureInterfaceMaterialHypoElastic::computeStress( double*              stateVars,
                                                                     Response&            response,
                                                                     Tangents&            tangents,
                                                                     const Deformation&   deformation,
                                                                     const TimeIncrement& timeIncrement )
{
  if ( !topMaterial || !bottomMaterial )
    throw std::logic_error( "MarmotZStabPressureInterfaceMaterialHypoElastic has no base material." );

  Vector3d n = Eigen::Map< const Vector3d >( deformation.normal );
  if ( n.norm() <= 1.0e-12 )
    throw std::invalid_argument( "MarmotZStabPressureInterfaceMaterialHypoElastic: normal is zero." );
  n /= n.norm();

  const Vector3d          d   = Eigen::Map< const Vector3d >( deformation.separationVector );
  const InterfaceGeometry geo = evaluateInterfaceGeometry( n, d, h );
  const double            ell = geo.normalSeparation;
  const double            st  = 0.5 * h; // side thickness

  const Eigen::Map< const Eigen::Matrix< double, 6, 1 > >  dU( deformation.dU );
  const Eigen::Map< const Eigen::Matrix< double, 18, 1 > > dA( deformation.dSurfaceStrain );
  const Vector3d dz = Eigen::Map< const Vector3d >( deformation.dNormalGradientJump );

  const double dpBar  = deformation.dPressureMean;
  const double dpJump = deformation.dPressureJump;

  // --- geometric operators ---
  const Matrix9x3RowMajor Bw   = jumpToGradient( n, ell );
  const Matrix9x3RowMajor Bn   = jumpToGradient( n, 1.0 );
  const Matrix3x9RowMajor Rt   = Bn.transpose();
  const Matrix9dRowMajor  Asrf = surfaceToGradient( n, geo.tangentialSeparation, ell );
  const Matrix9dRowMajor  I9   = Matrix9dRowMajor::Identity();
  const Matrix9dRowMajor  hIp  = 0.5 * ( I9 + Asrf );
  const Matrix9dRowMajor  hIm  = 0.5 * ( Asrf - I9 );
  const Matrix9dRowMajor  Pdev = deviatoricProjector();
  const Vector9d          eI   = Eigen::Map< const Vector9d >( Matrix3dRowMajor::Identity().eval().data() );

  Matrix9x21RowMajor BxP = Matrix9x21RowMajor::Zero(), BxM = Matrix9x21RowMajor::Zero();
  BxP.block< 9, 3 >( 0, offW )  = Bw;
  BxP.block< 9, 9 >( 0, offAp ) = hIp;
  BxP.block< 9, 9 >( 0, offAm ) = hIm;
  BxM.block< 9, 3 >( 0, offW )  = Bw;
  BxM.block< 9, 9 >( 0, offAp ) = hIm;
  BxM.block< 9, 9 >( 0, offAm ) = hIp;
  const Matrix9x3RowMajor BzP = 0.5 * Bn, BzM = -0.5 * Bn;

  Vector21d x             = Vector21d::Zero();
  x.segment< 3 >( offW )  = dU.segment< 3 >( 0 ) - dU.segment< 3 >( 3 );
  x.segment< 9 >( offAp ) = dA.segment< 9 >( 0 );
  x.segment< 9 >( offAm ) = dA.segment< 9 >( 9 );

  const Vector9d         dGP = BxP * x + BzP * dz;
  const Vector9d         dGM = BxM * x + BzM * dz;
  const Matrix3dRowMajor GP  = Eigen::Map< const Matrix3dRowMajor >( dGP.data() );
  const Matrix3dRowMajor GM  = Eigen::Map< const Matrix3dRowMajor >( dGM.data() );

  // --- committed state, one constitutive call per face ---
  double* topSv   = stateLayout.getPtr( stateVars, "topMaterialStateVars" );
  double* botSv   = stateLayout.getPtr( stateVars, "bottomMaterialStateVars" );
  double* sPStore = stateLayout.getPtr( stateVars, "stressPlus" );
  double* sMStore = stateLayout.getPtr( stateVars, "stressMinus" );

  const Matrix3dRowMajor sCurP = ContinuumMechanics::VoigtNotation::voigtToStress(
    Vector6d( Eigen::Map< const Vector6d >( sPStore ) ) );
  const Matrix3dRowMajor sCurM = ContinuumMechanics::VoigtNotation::voigtToStress(
    Vector6d( Eigen::Map< const Vector6d >( sMStore ) ) );

  const MarmotMaterialHypoElastic::timeInfo ti{ timeIncrement.timeOld, timeIncrement.dT };

  const SideTrial tp = evaluateSideTrial( *topMaterial,
                                          topSv,
                                          topMaterial->getNumberOfRequiredStateVars(),
                                          sCurP,
                                          ContinuumMechanics::VoigtNotation::strainToVoigt(
                                            Matrix3dRowMajor( 0.5 * ( GP + GP.transpose() ) ) ),
                                          ti );
  const SideTrial tm = evaluateSideTrial( *bottomMaterial,
                                          botSv,
                                          bottomMaterial->getNumberOfRequiredStateVars(),
                                          sCurM,
                                          ContinuumMechanics::VoigtNotation::strainToVoigt(
                                            Matrix3dRowMajor( 0.5 * ( GM + GM.transpose() ) ) ),
                                          ti );

  std::copy( tp.trialStateVars.begin(), tp.trialStateVars.end(), topSv );
  std::copy( tm.trialStateVars.begin(), tm.trialStateVars.end(), botSv );
  Eigen::Map< Vector6d > sPMap( sPStore ), sMMap( sMStore );
  sPMap = ContinuumMechanics::VoigtNotation::stressToVoigt(
    Eigen::Matrix3d( 0.5 * ( tp.stress + tp.stress.transpose() ) ) );
  sMMap = ContinuumMechanics::VoigtNotation::stressToVoigt(
    Eigen::Matrix3d( 0.5 * ( tm.stress + tm.stress.transpose() ) ) );

  // --- sigma~^pm = Pdev sigma^pm - p^pm I,  p^pm = pbar +- [p]/2 ---
  const double   pPlus = dpBar + 0.5 * dpJump, pMinus = dpBar - 0.5 * dpJump;
  const Vector9d stP = Pdev * Eigen::Map< const Vector9d >( tp.stress.data() ) - pPlus * eI;
  const Vector9d stM = Pdev * Eigen::Map< const Vector9d >( tm.stress.data() ) - pMinus * eI;

  const Matrix9dRowMajor CTp = Pdev * tp.CFull, CTm = Pdev * tm.CFull;

  // --- regularization on the z increment (see the Z kernel) ---
  const Matrix3dRowMajor Qe = shearModulus * Matrix3dRowMajor::Identity() +
                              ( lameParameter + shearModulus ) * ( n * n.transpose() );
  const Matrix3dRowMajor QReg = ( 0.25 * h * regularization ) * Qe;

  // --- generalized stresses ---
  const Vector21d pX = st * ( BxP.transpose() * stP + BxM.transpose() * stM );
  const Vector3d  pZ = st * ( BzP.transpose() * stP + BzM.transpose() * stM ) + QReg * dz;

  Eigen::Map< Vector3d >( response.generalizedForce )   = pX.segment< 3 >( offW );
  Eigen::Map< Vector9d >( response.surfaceStressPlus )  = pX.segment< 9 >( offAp );
  Eigen::Map< Vector9d >( response.surfaceStressMinus ) = pX.segment< 9 >( offAm );
  Eigen::Map< Vector3d >( response.tractionImbalance )  = pZ;

  response.volumetricResidualMean[0] = 0.5 * ( GP.trace() + GM.trace() ) + dpBar / bulkModulus;
  response.volumetricResidualJump[0] = ( GP.trace() - GM.trace() ) + dpJump / bulkModulus;

  // --- tangents ---
  Eigen::Map< Matrix21dRowMajor >( tangents.K_xx ) = st * ( BxP.transpose() * CTp * BxP + BxM.transpose() * CTm * BxM );
  Eigen::Map< Matrix21x3RowMajor >( tangents.K_xz ) = st *
                                                      ( BxP.transpose() * CTp * BzP + BxM.transpose() * CTm * BzM );
  Eigen::Map< Vector21d >( tangents.K_xpm ) = -st * ( ( BxP + BxM ).transpose() * eI );
  Eigen::Map< Vector21d >( tangents.K_xpj ) = -0.5 * st * ( ( BxP - BxM ).transpose() * eI );

  Eigen::Map< Matrix3x21RowMajor >( tangents.K_zx ) = st *
                                                      ( BzP.transpose() * CTp * BxP + BzM.transpose() * CTm * BxM );
  Eigen::Map< Matrix3dRowMajor >( tangents.K_zz ) = st * ( BzP.transpose() * CTp * BzP + BzM.transpose() * CTm * BzM ) +
                                                    QReg;
  // BzP + BzM == 0, so pbar cannot reach the traction jump: this block is exactly zero.
  Eigen::Map< Vector3d >( tangents.K_zpm ) = Vector3d::Zero();
  // -(h/4) n : the coupling the single-pressure form structurally cannot have
  Eigen::Map< Vector3d >( tangents.K_zpj ) = -0.5 * st * ( ( BzP - BzM ).transpose() * eI );

  Eigen::Map< Eigen::Matrix< double, 1, 21 > >( tangents.K_pmx ) = 0.5 * ( eI.transpose() * ( BxP + BxM ) );
  // eI^T (BzP + BzM) == 0: z cancels from the MEAN volumetric strain.
  Eigen::Map< Eigen::Matrix< double, 1, 3 > >( tangents.K_pmz ) = eI.transpose() * ( BzP + BzM );
  tangents.K_pmpm[0]                                            = 1.0 / bulkModulus;

  Eigen::Map< Eigen::Matrix< double, 1, 21 > >( tangents.K_pjx ) = eI.transpose() * ( BxP - BxM );
  // eI^T Bn == n^T : z.n IS the jump in volumetric strain
  Eigen::Map< Eigen::Matrix< double, 1, 3 > >( tangents.K_pjz ) = eI.transpose() * ( BzP - BzM );
  tangents.K_pjpj[0]                                            = 1.0 / bulkModulus;

  (void)Rt;
}

void MarmotZStabPressureInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int nStateVars )
{
  for ( int i = 0; i < nStateVars; ++i )
    stateVars[i] = 0.0;
  if ( !topMaterial || !bottomMaterial )
    return;
  topMaterial->initializeYourself( stateLayout.getPtr( stateVars, "topMaterialStateVars" ),
                                   topMaterial->getNumberOfRequiredStateVars() );
  bottomMaterial->initializeYourself( stateLayout.getPtr( stateVars, "bottomMaterialStateVars" ),
                                      bottomMaterial->getNumberOfRequiredStateVars() );
}

double MarmotZStabPressureInterfaceMaterialHypoElastic::getDensity()
{
  return topMaterial ? topMaterial->getDensity( nullptr ) : -1;
}
