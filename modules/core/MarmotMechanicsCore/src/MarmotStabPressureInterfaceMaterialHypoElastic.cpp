/* ---------------------------------------------------------------------
 *  Marmot - MarmotStabPressureInterfaceMaterialHypoElastic
 *  Alexandros Stathas alexandros.stathas@boku.ac.at
 * --------------------------------------------------------------------- */

#include "Marmot/MarmotStabPressureInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
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
  constexpr int    maxLocalNewtonIterations = 40;
  constexpr double localNewtonTolerance     = 1.0e-10;

  Matrix9dRowMajor fullGradientTangent( const Matrix6d& tangentVoigt )
  {
    Matrix9dRowMajor T = Matrix9dRowMajor::Zero();
    for ( int k = 0; k < 3; ++k )
      for ( int l = 0; l < 3; ++l ) {
        Matrix3dRowMajor dG        = Matrix3dRowMajor::Zero();
        dG( k, l )                 = 1.0;
        const Matrix3dRowMajor de  = 0.5 * ( dG + dG.transpose() );
        const Vector6d         dev = ContinuumMechanics::VoigtNotation::strainToVoigt( de );
        const Vector6d         dsv = tangentVoigt * dev;
        const Eigen::Matrix3d  ds  = ContinuumMechanics::VoigtNotation::voigtToStress( dsv );
        for ( int i = 0; i < 3; ++i )
          for ( int j = 0; j < 3; ++j )
            T( flatIndex( i, j ), flatIndex( k, l ) ) = ds( i, j );
      }
    return T;
  }

  /** 9x9 deviatoric projector acting on vec(sigma): dev = sigma - (tr sigma/3) I */
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

  InterfaceGeometry evaluateInterfaceGeometry( const Vector3d& n, const Vector3d& sep, double hh )
  {
    constexpr double tol = 1.0e-12;
    if ( hh <= 0.0 )
      throw std::invalid_argument( "StabPressureInterface: h must be positive." );
    if ( sep.norm() <= tol )
      return { hh, Vector3d::Zero() };
    const double ns = sep.dot( n );
    if ( ns <= tol )
      throw std::invalid_argument( "StabPressureInterface: non-positive normal separation." );
    return { ns, sep - ns * n };
  }

  Matrix9x3RowMajor jumpToGradient( const Vector3d& n, double ns )
  {
    Matrix9x3RowMajor m = Matrix9x3RowMajor::Zero();
    for ( int k = 0; k < 3; ++k )
      for ( int l = 0; l < 3; ++l )
        m( flatIndex( k, l ), k ) = n( l ) / ns;
    return m;
  }

  Matrix9dRowMajor surfaceToGradient( const Vector3d& n, const Vector3d& dtau, double ns )
  {
    Matrix9dRowMajor m = Matrix9dRowMajor::Zero();
    for ( int k = 0; k < 3; ++k )
      for ( int l = 0; l < 3; ++l )
        for ( int b = 0; b < 3; ++b )
          m( flatIndex( k, l ), flatIndex( k, b ) ) = ( l == b ? 1.0 : 0.0 ) - dtau( b ) * n( l ) / ns;
    return m;
  }

  Matrix3x9RowMajor stressToForce( const Vector3d& n )
  {
    Matrix3x9RowMajor m = Matrix3x9RowMajor::Zero();
    for ( int i = 0; i < 3; ++i )
      for ( int b = 0; b < 3; ++b )
        m( i, flatIndex( i, b ) ) = n( b );
    return m;
  }

  struct SideTrial {
    Matrix3dRowMajor      stress;
    Matrix9dRowMajor      CFull;
    std::vector< double > trialStateVars;
  };

  SideTrial evaluateSideTrial( MarmotMaterialHypoElastic&                 mat,
                               const double*                              oldSv,
                               int                                        nSv,
                               const Matrix3dRowMajor&                    stressCur,
                               const Vector6d&                            dEpsVoigt,
                               const MarmotMaterialHypoElastic::timeInfo& ti )
  {
    SideTrial t;
    t.trialStateVars.assign( oldSv, oldSv + nSv );
    const Eigen::Matrix3d              sSym( 0.5 * ( stressCur + stressCur.transpose() ) );
    Matrix6d                           C = Matrix6d::Zero();
    MarmotMaterialHypoElastic::state3D st( ContinuumMechanics::VoigtNotation::stressToVoigt( sSym ),
                                           0.0,
                                           0.0,
                                           t.trialStateVars.data() );
    mat.computeStress( st, C, dEpsVoigt, ti );
    t.stress = ContinuumMechanics::VoigtNotation::voigtToStress( st.stress );
    t.CFull  = fullGradientTangent( C );
    return t;
  }

} // namespace

MarmotStabPressureInterfaceMaterialHypoElastic::MarmotStabPressureInterfaceMaterialHypoElastic(
  const std::string& materialName,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : materialProperties( matProperties_ ), nMaterialProperties( nMaterialProperties_ ), materialNumber( materialNumber_ )
{
  if ( nMaterialProperties < 3 )
    throw std::invalid_argument( "StabPressureInterface requires at least E, nu, h." );
  const double E  = materialProperties[0];
  const double nu = materialProperties[1];
  h               = materialProperties[2];
  if ( h <= 0.0 )
    throw std::invalid_argument( "StabPressureInterface requires h > 0." );
  if ( std::abs( 1.0 - 2.0 * nu ) < 1.0e-12 )
    throw std::invalid_argument( "StabPressureInterface: nu = 0.5 (K infinite) is not supported." );
  bulkModulus  = E / ( 3.0 * ( 1.0 - 2.0 * nu ) );
  shearModulus = E / ( 2.0 * ( 1.0 + nu ) );

  baseMaterialProperties.reserve( nMaterialProperties - 1 );
  baseMaterialProperties.push_back( E );
  baseMaterialProperties.push_back( nu );
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
  stateLayout.add( "normalGradientJump", 3 );
  stateLayout.add( "stressPlus", 9 );
  stateLayout.add( "stressMinus", 9 );
  stateLayout.add( "pressure", 1 );
  stateLayout.finalize();
}

void MarmotStabPressureInterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  if ( topMaterial )
    topMaterial->setCharacteristicElementLength( length );
  if ( bottomMaterial )
    bottomMaterial->setCharacteristicElementLength( length );
}

void MarmotStabPressureInterfaceMaterialHypoElastic::computeStress( double*              stateVars,
                                                                    Response&            response,
                                                                    Tangents&            tangents,
                                                                    const Deformation&   deformation,
                                                                    const TimeIncrement& timeIncrement )
{
  if ( !topMaterial || !bottomMaterial )
    throw std::logic_error( "StabPressureInterface: no base material." );

  Vector3d n = Eigen::Map< const Vector3d >( deformation.normal );
  if ( n.norm() <= 1.0e-12 )
    throw std::invalid_argument( "StabPressureInterface: zero normal." );
  n.normalize();
  const Vector3d sep = Eigen::Map< const Vector3d >( deformation.separationVector );

  const auto   geo  = evaluateInterfaceGeometry( n, sep, h );
  const double ell  = geo.normalSeparation;
  const auto&  dtau = geo.tangentialSeparation;
  const double st   = 0.5 * h; // side thickness

  const Eigen::Map< const Eigen::Matrix< double, 6, 1 > >  dU( deformation.dU );
  const Eigen::Map< const Eigen::Matrix< double, 18, 1 > > dS( deformation.dSurfaceStrain );
  const Vector3d                                           w = dU.segment< 3 >( 0 ) - dU.segment< 3 >( 3 );
  const Eigen::Map< const Matrix3dRowMajor >               APlus( dS.data() );
  const Eigen::Map< const Matrix3dRowMajor >               AMinus( dS.data() + 9 );
  const Matrix3dRowMajor                                   ABar = 0.5 * ( APlus + AMinus );
  const Vector3d                                           gBar = ( w - ABar * dtau ) / ell;

  double*   topSv    = stateLayout.getPtr( stateVars, "topMaterialStateVars" );
  double*   botSv    = stateLayout.getPtr( stateVars, "bottomMaterialStateVars" );
  double*   gammaPtr = stateLayout.getPtr( stateVars, "normalGradientJump" );
  double*   spPtr    = stateLayout.getPtr( stateVars, "stressPlus" );
  double*   smPtr    = stateLayout.getPtr( stateVars, "stressMinus" );
  double*   pPtr     = stateLayout.getPtr( stateVars, "pressure" );
  const int nTopSv   = topMaterial->getNumberOfRequiredStateVars();
  const int nBotSv   = bottomMaterial->getNumberOfRequiredStateVars();

  // total pressure at the end of the step
  const double pTotal = pPtr[0] + deformation.dPressure;

  const Matrix3dRowMajor sCurPlus  = Eigen::Map< const Matrix3dRowMajor >( spPtr );
  const Matrix3dRowMajor sCurMinus = Eigen::Map< const Matrix3dRowMajor >( smPtr );

  const MarmotMaterialHypoElastic::timeInfo ti{ timeIncrement.timeOld, timeIncrement.dT };

  const Matrix3x9RowMajor Rt   = stressToForce( n );
  const Matrix9x3RowMajor Bn   = jumpToGradient( n, 1.0 );
  const Matrix9dRowMajor  Pdev = deviatoricProjector();
  const Vector9d          eI   = Eigen::Map< const Vector9d >( Matrix3dRowMajor::Identity().eval().data() );

  struct Local {
    SideTrial p, m;
    Vector3d  r;
    double    rNorm;
  };

  // sigma~ = Pdev sigma - p I ; the -p I cancels in (sigma~+ - sigma~-)
  auto evaluate = [&]( const Vector3d& gam ) -> Local {
    const Vector3d         gP = gBar + 0.5 * gam, gM = gBar - 0.5 * gam;
    const Matrix3dRowMajor GP  = APlus + gP * n.transpose();
    const Matrix3dRowMajor GM  = AMinus + gM * n.transpose();
    const Vector6d         deP = ContinuumMechanics::VoigtNotation::strainToVoigt( 0.5 * ( GP + GP.transpose() ) );
    const Vector6d         deM = ContinuumMechanics::VoigtNotation::strainToVoigt( 0.5 * ( GM + GM.transpose() ) );
    Local                  L;
    L.p                 = evaluateSideTrial( *topMaterial, topSv, nTopSv, sCurPlus, deP, ti );
    L.m                 = evaluateSideTrial( *bottomMaterial, botSv, nBotSv, sCurMinus, deM, ti );
    const Vector9d devP = Pdev * Eigen::Map< const Vector9d >( L.p.stress.data() );
    const Vector9d devM = Pdev * Eigen::Map< const Vector9d >( L.m.stress.data() );
    L.r                 = Rt * ( devP - devM );
    L.rNorm             = L.r.norm() / std::max( 1.0, std::max( ( Rt * devP ).norm(), ( Rt * devM ).norm() ) );
    return L;
  };

  Vector3d gamma = Eigen::Map< const Vector3d >( gammaPtr );
  Local    cur   = evaluate( gamma );
  bool     conv  = cur.rNorm <= localNewtonTolerance;

  for ( int it = 0; it < maxLocalNewtonIterations && !conv; ++it ) {
    const Matrix3dRowMajor Qp   = Rt * ( Pdev * cur.p.CFull ) * Bn;
    const Matrix3dRowMajor Qm   = Rt * ( Pdev * cur.m.CFull ) * Bn;
    const Matrix3dRowMajor Qavg = 0.5 * ( Qp + Qm );
    const Vector3d         dg   = Qavg.fullPivLu().solve( -cur.r );
    double                 al   = 1.0;
    Local                  cand;
    bool                   acc = false;
    for ( int ls = 0; ls < 12; ++ls ) {
      cand = evaluate( gamma + al * dg );
      if ( cand.rNorm <= localNewtonTolerance || cand.rNorm < cur.rNorm ) {
        acc = true;
        break;
      }
      al *= 0.5;
    }
    gamma += al * dg;
    cur = acc ? cand : evaluate( gamma );
    if ( cur.rNorm <= localNewtonTolerance )
      conv = true;
  }
  if ( !conv )
    throw Marmot::StressUpdateFailed( "MarmotStabPressureInterfaceMaterialHypoElastic: local traction-equilibrium "
                                      "solve did not converge." );

  // commit
  Eigen::Map< Vector3d >         gMap( gammaPtr );
  Eigen::Map< Matrix3dRowMajor > spMap( spPtr ), smMap( smPtr );
  gMap    = gamma;
  spMap   = cur.p.stress;
  smMap   = cur.m.stress;
  pPtr[0] = pTotal;
  std::copy( cur.p.trialStateVars.begin(), cur.p.trialStateVars.end(), topSv );
  std::copy( cur.m.trialStateVars.begin(), cur.m.trialStateVars.end(), botSv );

  // ---- geometric maps (independent of state) ----
  const Matrix9x3RowMajor Bw  = jumpToGradient( n, ell );
  const Matrix9dRowMajor  Asr = surfaceToGradient( n, dtau, ell );
  const Matrix9dRowMajor  I9  = Matrix9dRowMajor::Identity();
  const Matrix9dRowMajor  hp  = 0.5 * ( I9 + Asr );
  const Matrix9dRowMajor  hm  = 0.5 * ( Asr - I9 );

  Matrix9x21RowMajor BxP = Matrix9x21RowMajor::Zero(), BxM = Matrix9x21RowMajor::Zero();
  BxP.block< 9, 3 >( 0, 0 )   = Bw;
  BxP.block< 9, 9 >( 0, 3 )   = hp;
  BxP.block< 9, 9 >( 0, 12 )  = hm;
  BxM.block< 9, 3 >( 0, 0 )   = Bw;
  BxM.block< 9, 9 >( 0, 3 )   = hm;
  BxM.block< 9, 9 >( 0, 12 )  = hp;
  const Matrix9x3RowMajor BzP = 0.5 * Bn, BzM = -0.5 * Bn;

  // ---- projected material tangents ----
  const Matrix9dRowMajor CTp = Pdev * cur.p.CFull;
  const Matrix9dRowMajor CTm = Pdev * cur.m.CFull;

  const Matrix3dRowMajor   Qavg = 0.5 * ( Rt * ( CTp + CTm ) * Bn ); // = rgam
  const Matrix3x21RowMajor rx   = Rt * ( CTp * BxP - CTm * BxM );

  const Matrix21dRowMajor  px   = st * ( BxP.transpose() * CTp * BxP + BxM.transpose() * CTm * BxM );
  const Matrix21x3RowMajor pgam = st * ( BxP.transpose() * CTp * BzP + BxM.transpose() * CTm * BzM );
  const Matrix21dRowMajor  Kc   = px - pgam * Qavg.fullPivLu().solve( rx );

  // ---- generalized stresses from sigma~ ----
  const Vector9d  stP = Pdev * Eigen::Map< const Vector9d >( cur.p.stress.data() ) - pTotal * eI;
  const Vector9d  stM = Pdev * Eigen::Map< const Vector9d >( cur.m.stress.data() ) - pTotal * eI;
  const Vector21d y   = st * ( BxP.transpose() * stP + BxM.transpose() * stM );

  Eigen::Map< Vector3d >( response.generalizedForce )          = y.segment< 3 >( 0 );
  Eigen::Map< Matrix3dRowMajor >( response.surfaceStressPlus ) = Eigen::Map< const Matrix3dRowMajor >(
    y.segment< 9 >( 3 ).eval().data() );
  Eigen::Map< Matrix3dRowMajor >( response.surfaceStressMinus ) = Eigen::Map< const Matrix3dRowMajor >(
    y.segment< 9 >( 12 ).eval().data() );

  // ---- volumetric residual: <tr eps> + p/K   (independent of gamma) ----
  const Vector3d         gP = gBar + 0.5 * gamma, gM = gBar - 0.5 * gamma;
  const Matrix3dRowMajor GP = APlus + gP * n.transpose();
  const Matrix3dRowMajor GM = AMinus + gM * n.transpose();
  // INCREMENTAL volumetric law: d<tr eps> + dp/K = 0.  Both the volumetric
  // strain and p accumulate from zero, so enforcing the increment each step is
  // equivalent to the total relation p = -K tr(eps) -- and it is consistent with
  // the hypoelastic (rate) form used for the deviatoric response. Mixing the
  // strain INCREMENT with the TOTAL pressure here would under-predict p by the
  // number of load steps.
  response.volumetricResidual[0] = 0.5 * ( GP.trace() + GM.trace() ) + deformation.dPressure / bulkModulus;

  // ---- tangents ----
  Eigen::Map< Matrix3dRowMajor >( tangents.Q_ww )   = Kc.block< 3, 3 >( 0, 0 );
  Eigen::Map< Matrix3x9RowMajor >( tangents.Q_wAp ) = Kc.block< 3, 9 >( 0, 3 );
  Eigen::Map< Matrix3x9RowMajor >( tangents.Q_wAm ) = Kc.block< 3, 9 >( 0, 12 );
  Eigen::Map< Matrix9x3RowMajor >( tangents.Q_Apw ) = Kc.block< 9, 3 >( 3, 0 );
  Eigen::Map< Matrix9dRowMajor >( tangents.Q_ApAp ) = Kc.block< 9, 9 >( 3, 3 );
  Eigen::Map< Matrix9dRowMajor >( tangents.Q_ApAm ) = Kc.block< 9, 9 >( 3, 12 );
  Eigen::Map< Matrix9x3RowMajor >( tangents.Q_Amw ) = Kc.block< 9, 3 >( 12, 0 );
  Eigen::Map< Matrix9dRowMajor >( tangents.Q_AmAp ) = Kc.block< 9, 9 >( 12, 3 );
  Eigen::Map< Matrix9dRowMajor >( tangents.Q_AmAm ) = Kc.block< 9, 9 >( 12, 12 );

  // dy/dp = -st ( BxP^T e + BxM^T e )     (purely geometric; dgamma/dp = 0)
  const Vector21d dydp                     = -st * ( BxP.transpose() * eI + BxM.transpose() * eI );
  Eigen::Map< Vector3d >( tangents.Q_wp )  = dydp.segment< 3 >( 0 );
  Eigen::Map< Vector9d >( tangents.Q_App ) = dydp.segment< 9 >( 3 );
  Eigen::Map< Vector9d >( tangents.Q_Amp ) = dydp.segment< 9 >( 12 );

  // dr_p/dx = 1/2 e^T ( BxP + BxM )       (purely geometric; dr_p/dgamma = 0)
  const Eigen::Matrix< double, 1, 21 > drdx                     = 0.5 * ( eI.transpose() * ( BxP + BxM ) );
  Eigen::Map< Eigen::Matrix< double, 1, 3 > >( tangents.Q_pw )  = drdx.block< 1, 3 >( 0, 0 );
  Eigen::Map< Eigen::Matrix< double, 1, 9 > >( tangents.Q_pAp ) = drdx.block< 1, 9 >( 0, 3 );
  Eigen::Map< Eigen::Matrix< double, 1, 9 > >( tangents.Q_pAm ) = drdx.block< 1, 9 >( 0, 12 );
  tangents.Q_pp[0]                                              = 1.0 / bulkModulus;
}

void MarmotStabPressureInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int nStateVars )
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

double MarmotStabPressureInterfaceMaterialHypoElastic::getDensity()
{
  return topMaterial ? topMaterial->getDensity( nullptr ) : -1;
}
