#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::FastorStandardTensors;

namespace {

  using Vector3d           = Eigen::Matrix< double, 3, 1 >;
  using Vector9d           = Eigen::Matrix< double, 9, 1 >;
  using Vector21d          = Eigen::Matrix< double, 21, 1 >;
  using Matrix3dRowMajor   = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor   = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor  = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;
  using Matrix3x9RowMajor  = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;
  using Matrix9x21RowMajor = Eigen::Matrix< double, 9, 21, Eigen::RowMajor >;
  using Matrix21x9RowMajor = Eigen::Matrix< double, 21, 9, Eigen::RowMajor >;
  using Matrix21x3RowMajor = Eigen::Matrix< double, 21, 3, Eigen::RowMajor >;
  using Matrix3x21RowMajor = Eigen::Matrix< double, 3, 21, Eigen::RowMajor >;
  using Matrix21dRowMajor  = Eigen::Matrix< double, 21, 21, Eigen::RowMajor >;

  constexpr int flatIndex( int i, int j )
  {
    return 3 * i + j;
  }

  constexpr int    maxLocalNewtonIterations = 40;
  constexpr double localNewtonTolerance     = 1.0e-10;

  /** Identical to the corresponding helper in MarmotXInterfaceMaterialHypoElastic.cpp. */
  Matrix9dRowMajor fullGradientTangent( const Matrix6d& tangentVoigt )
  {
    Matrix9dRowMajor tangentFull = Matrix9dRowMajor::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        Matrix3dRowMajor dGradient = Matrix3dRowMajor::Zero();
        dGradient( k, l )          = 1.0;

        const Matrix3dRowMajor dStrain      = 0.5 * ( dGradient + dGradient.transpose() );
        const Vector6d         dStrainVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt( dStrain );
        const Vector6d         dStressVoigt = tangentVoigt * dStrainVoigt;
        const Eigen::Matrix3d  dStress      = ContinuumMechanics::VoigtNotation::voigtToStress( dStressVoigt );

        const int column = flatIndex( k, l );
        for ( int i = 0; i < 3; ++i ) {
          for ( int j = 0; j < 3; ++j ) {
            tangentFull( flatIndex( i, j ), column ) = dStress( i, j );
          }
        }
      }
    }

    return tangentFull;
  }

  struct InterfaceGeometry {
    double   normalSeparation;
    Vector3d tangentialSeparation;
  };

  InterfaceGeometry evaluateInterfaceGeometry( const Vector3d& normal,
                                               const Vector3d& separationVector,
                                               double          constitutiveThickness )
  {
    constexpr double tolerance = 1.0e-12;

    if ( constitutiveThickness <= 0.0 ) {
      throw std::invalid_argument(
        "MarmotEquilibratedXInterfaceMaterialHypoElastic: interface thickness h must be positive." );
    }

    if ( separationVector.norm() <= tolerance ) {
      return { constitutiveThickness, Vector3d::Zero() };
    }

    const double normalSeparation = separationVector.dot( normal );
    if ( normalSeparation <= tolerance ) {
      throw std::invalid_argument( "MarmotEquilibratedXInterfaceMaterialHypoElastic: the top-bottom connector must "
                                   "have a positive normal component." );
    }

    return { normalSeparation, separationVector - normalSeparation * normal };
  }

  /** Maps a 3-vector g to vec(g \otimes n) / normalSeparation. With
   * normalSeparation=1 this is B_n; with normalSeparation=ell this is B_w. */
  Matrix9x3RowMajor jumpToGradient( const Vector3d& normal, double normalSeparation )
  {
    Matrix9x3RowMajor map = Matrix9x3RowMajor::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        map( flatIndex( k, l ), k ) = normal( l ) / normalSeparation;
      }
    }

    return map;
  }

  /** Maps vec(A) to vec(A) - (1/normalSeparation) vec((A d_tau) \otimes n), i.e. I_A + B_tau. */
  Matrix9dRowMajor surfaceToGradient( const Vector3d& normal,
                                      const Vector3d& tangentialSeparation,
                                      double          normalSeparation )
  {
    Matrix9dRowMajor map = Matrix9dRowMajor::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        const int row = flatIndex( k, l );

        for ( int a = 0; a < 3; ++a ) {
          for ( int b = 0; b < 3; ++b ) {
            if ( k != a )
              continue;

            const double identityPart     = l == b ? 1.0 : 0.0;
            map( row, flatIndex( a, b ) ) = identityPart - tangentialSeparation( b ) * normal( l ) / normalSeparation;
          }
        }
      }
    }

    return map;
  }

  /** Maps vec(sigma) to (sideThickness/normalSeparation) * sigma n. With
   * (sideThickness=1, normalSeparation=1) this is R_t. */
  Matrix3x9RowMajor stressToForce( const Vector3d& normal, double sideThickness, double normalSeparation )
  {
    Matrix3x9RowMajor map   = Matrix3x9RowMajor::Zero();
    const double      scale = sideThickness / normalSeparation;

    for ( int i = 0; i < 3; ++i ) {
      for ( int b = 0; b < 3; ++b ) {
        map( i, flatIndex( i, b ) ) = scale * normal( b );
      }
    }

    return map;
  }

  /** One side's trial constitutive response, evaluated from a scratch copy of
   * the OLD (committed) state -- never mutates the persistent state vars. */
  struct SideTrial {
    Matrix3dRowMajor      stress;
    Matrix9dRowMajor      CFull;
    std::vector< double > trialStateVars;
  };

  SideTrial evaluateSideTrial( MarmotMaterialHypoElastic&                 material,
                               const double*                              oldStateVars,
                               int                                        nStateVars,
                               const Matrix3dRowMajor&                    stressCurrent,
                               const Vector6d&                            strainIncrementVoigt,
                               const MarmotMaterialHypoElastic::timeInfo& timeInfo )
  {
    SideTrial trial;
    trial.trialStateVars.assign( oldStateVars, oldStateVars + nStateVars );

    const Eigen::Matrix3d stressCurrentSym( 0.5 * ( stressCurrent + stressCurrent.transpose() ) );
    const Vector6d        stressVoigt = ContinuumMechanics::VoigtNotation::stressToVoigt( stressCurrentSym );

    Matrix6d                           tangentVoigt = Matrix6d::Zero();
    MarmotMaterialHypoElastic::state3D baseState{ stressVoigt, 0.0, 0.0, trial.trialStateVars.data() };

    material.computeStress( baseState, tangentVoigt, strainIncrementVoigt, timeInfo );

    const Eigen::Matrix3d stressUpdatedEigen = ContinuumMechanics::VoigtNotation::voigtToStress( baseState.stress );
    trial.stress                             = stressUpdatedEigen;
    trial.CFull                              = fullGradientTangent( tangentVoigt );

    return trial;
  }

} // namespace

MarmotEquilibratedXInterfaceMaterialHypoElastic::MarmotEquilibratedXInterfaceMaterialHypoElastic(
  const std::string& materialName,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : materialProperties( matProperties_ ), nMaterialProperties( nMaterialProperties_ ), materialNumber( materialNumber_ )
{
  if ( nMaterialProperties < 3 ) {
    throw std::invalid_argument(
      "MarmotEquilibratedXInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h." );
  }

  h = materialProperties[2];
  if ( h <= 0.0 ) {
    throw std::invalid_argument( "MarmotEquilibratedXInterfaceMaterialHypoElastic requires h > 0." );
  }

  baseMaterialProperties.reserve( nMaterialProperties - 1 );
  baseMaterialProperties.push_back( materialProperties[0] );
  baseMaterialProperties.push_back( materialProperties[1] );
  baseMaterialProperties.insert( baseMaterialProperties.end(),
                                 materialProperties + 3,
                                 materialProperties + nMaterialProperties );

  topMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     baseMaterialProperties.data(),
                                                                     static_cast< int >(
                                                                       baseMaterialProperties.size() ),
                                                                     materialNumber ) );
  bottomMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     baseMaterialProperties.data(),
                                                                     static_cast< int >(
                                                                       baseMaterialProperties.size() ),
                                                                     materialNumber ) );

  stateLayout.add( "topMaterialStateVars", topMaterial->getNumberOfRequiredStateVars() );
  stateLayout.add( "bottomMaterialStateVars", bottomMaterial->getNumberOfRequiredStateVars() );
  stateLayout.add( "normalGradientJump", 3 );
  stateLayout.finalize();
}

void MarmotEquilibratedXInterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  if ( topMaterial ) {
    topMaterial->setCharacteristicElementLength( length );
  }
  if ( bottomMaterial ) {
    bottomMaterial->setCharacteristicElementLength( length );
  }
}

void MarmotEquilibratedXInterfaceMaterialHypoElastic::computeStress( State&               state,
                                                                     Tangents&            tangents,
                                                                     const Deformation&   deformation,
                                                                     const TimeIncrement& timeIncrement )
{
  if ( !topMaterial || !bottomMaterial ) {
    throw std::logic_error( "MarmotEquilibratedXInterfaceMaterialHypoElastic has no base material." );
  }

  const Eigen::Map< const Vector3d > normalMap( deformation.normal.data() );
  const Eigen::Map< const Vector3d > separationMap( deformation.separationVector.data() );

  Vector3d     normal     = normalMap;
  const double normalNorm = normal.norm();
  if ( normalNorm <= 1.0e-12 ) {
    throw std::invalid_argument( "MarmotEquilibratedXInterfaceMaterialHypoElastic: interface normal is zero." );
  }
  normal /= normalNorm;

  const InterfaceGeometry geometry      = evaluateInterfaceGeometry( normal, separationMap, h );
  const double            ell           = geometry.normalSeparation;
  const Vector3d&         dTangential   = geometry.tangentialSeparation;
  const double            sideThickness = 0.5 * h;

  const Eigen::Map< const Eigen::Matrix< double, 6, 1 > >  dU( deformation.dU.data() );
  const Eigen::Map< const Eigen::Matrix< double, 18, 1 > > dSurfaceGradient( deformation.dSurfaceStrain.data() );

  const Vector3d w = dU.template segment< 3 >( 0 ) - dU.template segment< 3 >( 3 );

  const Eigen::Map< const Matrix3dRowMajor > APlus( dSurfaceGradient.data() );
  const Eigen::Map< const Matrix3dRowMajor > AMinus( dSurfaceGradient.data() + 9 );
  const Matrix3dRowMajor                     ABar = 0.5 * ( APlus + AMinus );

  const Vector3d gBar = ( w - ABar * dTangential ) / ell;

  // Recover both sides' OLD (committed, at start of increment) Cauchy stress
  // from the persisted equilibrated generalized state. A single stored
  // generalizedForce is used for both sides' reconstruction -- consistent
  // with the equilibrium constraint (t+ == t- == f * ell / h at every
  // committed state).
  const Eigen::Map< const Vector3d >         forceOld( state.generalizedForce.data() );
  const Eigen::Map< const Matrix3dRowMajor > surfaceStressPlusOld( state.surfaceStressPlus.data() );
  const Eigen::Map< const Matrix3dRowMajor > surfaceStressMinusOld( state.surfaceStressMinus.data() );

  const Matrix3dRowMajor stressCurrentPlus = ( surfaceStressPlusOld + forceOld * dTangential.transpose() ) /
                                             sideThickness;
  const Matrix3dRowMajor stressCurrentMinus = ( surfaceStressMinusOld + forceOld * dTangential.transpose() ) /
                                              sideThickness;

  double*   topStateVars       = stateLayout.getPtr( state.stateVars, "topMaterialStateVars" );
  double*   bottomStateVars    = stateLayout.getPtr( state.stateVars, "bottomMaterialStateVars" );
  double*   normalGradientJump = stateLayout.getPtr( state.stateVars, "normalGradientJump" );
  const int nTopStateVars      = topMaterial->getNumberOfRequiredStateVars();
  const int nBottomStateVars   = bottomMaterial->getNumberOfRequiredStateVars();

  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld, timeIncrement.dT };

  // Local Newton solve for the normal-gradient jump z, enforcing
  // sigma+(z) n == sigma-(z) n. Every trial evaluation restarts from the OLD
  // committed state (via evaluateSideTrial's scratch copy) -- state is never
  // accumulated across local iterations.
  const Matrix3x9RowMajor Rt = stressToForce( normal, 1.0, 1.0 );
  const Matrix9x3RowMajor Bn = jumpToGradient( normal, 1.0 );

  struct LocalEvaluation {
    SideTrial plus;
    SideTrial minus;
    Vector3d  r;
    double    rNorm;
  };

  auto evaluateLocal = [&]( const Vector3d& zTrial ) -> LocalEvaluation {
    const Vector3d gPlus  = gBar + 0.5 * zTrial;
    const Vector3d gMinus = gBar - 0.5 * zTrial;

    const Matrix3dRowMajor GPlus  = APlus + gPlus * normal.transpose();
    const Matrix3dRowMajor GMinus = AMinus + gMinus * normal.transpose();

    const Matrix3dRowMajor strainIncrementPlus  = 0.5 * ( GPlus + GPlus.transpose() );
    const Matrix3dRowMajor strainIncrementMinus = 0.5 * ( GMinus + GMinus.transpose() );

    const Vector6d strainIncrementPlusVoigt  = ContinuumMechanics::VoigtNotation::strainToVoigt( strainIncrementPlus );
    const Vector6d strainIncrementMinusVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt( strainIncrementMinus );

    LocalEvaluation evaluation;
    evaluation.plus  = evaluateSideTrial( *topMaterial,
                                         topStateVars,
                                         nTopStateVars,
                                         stressCurrentPlus,
                                         strainIncrementPlusVoigt,
                                         timeInfo );
    evaluation.minus = evaluateSideTrial( *bottomMaterial,
                                          bottomStateVars,
                                          nBottomStateVars,
                                          stressCurrentMinus,
                                          strainIncrementMinusVoigt,
                                          timeInfo );

    const Vector3d tPlus  = evaluation.plus.stress * normal;
    const Vector3d tMinus = evaluation.minus.stress * normal;
    evaluation.r          = tPlus - tMinus;
    evaluation.rNorm      = evaluation.r.norm() / std::max( 1.0, std::max( tPlus.norm(), tMinus.norm() ) );
    return evaluation;
  };

  Vector3d z = Eigen::Map< const Vector3d >( normalGradientJump );

  LocalEvaluation current   = evaluateLocal( z );
  bool            converged = current.rNorm <= localNewtonTolerance;

  for ( int iteration = 0; iteration < maxLocalNewtonIterations && !converged; ++iteration ) {
    const Matrix3dRowMajor Qplus  = Rt * current.plus.CFull * Bn;
    const Matrix3dRowMajor Qminus = Rt * current.minus.CFull * Bn;
    const Matrix3dRowMajor Qavg   = 0.5 * ( Qplus + Qminus );

    const Vector3d dz = Qavg.fullPivLu().solve( -current.r );

    // Backtracking line search: a full Newton step can overshoot when a
    // trial state crosses a yield-surface kink, so only accept a step that
    // reduces the residual (or is already converged), halving otherwise.
    double          alpha = 1.0;
    LocalEvaluation candidate;
    bool            accepted = false;
    for ( int lineSearchIter = 0; lineSearchIter < 12; ++lineSearchIter ) {
      candidate = evaluateLocal( z + alpha * dz );
      if ( candidate.rNorm <= localNewtonTolerance || candidate.rNorm < current.rNorm ) {
        accepted = true;
        break;
      }
      alpha *= 0.5;
    }

    z += alpha * dz;
    current = accepted ? candidate : evaluateLocal( z );

    if ( current.rNorm <= localNewtonTolerance ) {
      converged = true;
    }
  }

  const SideTrial& plusTrial  = current.plus;
  const SideTrial& minusTrial = current.minus;

  if ( !converged ) {
    throw Marmot::StressUpdateFailed(
      "MarmotEquilibratedXInterfaceMaterialHypoElastic: local traction-equilibrium solve did not converge." );
  }

  // Commit: persist the converged internal variable (warm start for the
  // next increment) and the two sides' updated material state.
  Eigen::Map< Vector3d > normalGradientJumpMap( normalGradientJump );
  normalGradientJumpMap = z;
  std::copy( plusTrial.trialStateVars.begin(), plusTrial.trialStateVars.end(), topStateVars );
  std::copy( minusTrial.trialStateVars.begin(), minusTrial.trialStateVars.end(), bottomStateVars );

  // Purely geometric operators (independent of z and of the material
  // state), built from the same primitives as
  // MarmotXInterfaceMaterialHypoElastic/MarmotCorrectedInterfaceMaterialHypoElastic.
  // (Rt, Bn already declared above, ahead of the local Newton loop.)
  const Matrix9x3RowMajor Bw                = jumpToGradient( normal, ell );
  const Matrix9dRowMajor  Asrf              = surfaceToGradient( normal, dTangential, ell ); // = I_A + B_tau
  const Matrix9dRowMajor  I9                = Matrix9dRowMajor::Identity();
  const Matrix9dRowMajor  half_I_plus_Asrf  = 0.5 * ( I9 + Asrf );
  const Matrix9dRowMajor  half_Asrf_minus_I = 0.5 * ( Asrf - I9 );

  Matrix9x21RowMajor BxPlus      = Matrix9x21RowMajor::Zero();
  Matrix9x21RowMajor BxMinus     = Matrix9x21RowMajor::Zero();
  BxPlus.block< 9, 3 >( 0, 0 )   = Bw;
  BxPlus.block< 9, 9 >( 0, 3 )   = half_I_plus_Asrf;
  BxPlus.block< 9, 9 >( 0, 12 )  = half_Asrf_minus_I;
  BxMinus.block< 9, 3 >( 0, 0 )  = Bw;
  BxMinus.block< 9, 9 >( 0, 3 )  = half_Asrf_minus_I;
  BxMinus.block< 9, 9 >( 0, 12 ) = half_I_plus_Asrf;

  const Matrix9x3RowMajor BzPlus  = 0.5 * Bn;
  const Matrix9x3RowMajor BzMinus = -0.5 * Bn;

  const Matrix3dRowMajor Qplus  = Rt * plusTrial.CFull * Bn;
  const Matrix3dRowMajor Qminus = Rt * minusTrial.CFull * Bn;
  const Matrix3dRowMajor Qavg   = 0.5 * ( Qplus + Qminus );

  const Matrix3x21RowMajor rx = Rt * ( plusTrial.CFull * BxPlus - minusTrial.CFull * BxMinus );

  const Matrix21dRowMajor  px = sideThickness * ( BxPlus.transpose() * plusTrial.CFull * BxPlus +
                                                 BxMinus.transpose() * minusTrial.CFull * BxMinus );
  const Matrix21x3RowMajor pz = sideThickness * ( BxPlus.transpose() * plusTrial.CFull * BzPlus +
                                                  BxMinus.transpose() * minusTrial.CFull * BzMinus );

  const Matrix3x21RowMajor QavgInvRx = Qavg.fullPivLu().solve( rx );
  const Matrix21dRowMajor  Kcond     = px - pz * QavgInvRx;

  const Eigen::Map< const Vector9d > vecSigmaPlus( plusTrial.stress.data() );
  const Eigen::Map< const Vector9d > vecSigmaMinus( minusTrial.stress.data() );
  const Vector21d p = sideThickness * ( BxPlus.transpose() * vecSigmaPlus + BxMinus.transpose() * vecSigmaMinus );

  Eigen::Map< Vector3d >( state.generalizedForce.data() )          = p.segment< 3 >( 0 );
  Eigen::Map< Matrix3dRowMajor >( state.surfaceStressPlus.data() ) = Eigen::Map< const Matrix3dRowMajor >(
    p.segment< 9 >( 3 ).eval().data() );
  Eigen::Map< Matrix3dRowMajor >( state.surfaceStressMinus.data() ) = Eigen::Map< const Matrix3dRowMajor >(
    p.segment< 9 >( 12 ).eval().data() );

  Eigen::Map< Matrix3dRowMajor >( tangents.Q_ww.data() )   = Kcond.block< 3, 3 >( 0, 0 );
  Eigen::Map< Matrix3x9RowMajor >( tangents.Q_wAp.data() ) = Kcond.block< 3, 9 >( 0, 3 );
  Eigen::Map< Matrix3x9RowMajor >( tangents.Q_wAm.data() ) = Kcond.block< 3, 9 >( 0, 12 );

  Eigen::Map< Matrix9x3RowMajor >( tangents.Q_Apw.data() ) = Kcond.block< 9, 3 >( 3, 0 );
  Eigen::Map< Matrix9dRowMajor >( tangents.Q_ApAp.data() ) = Kcond.block< 9, 9 >( 3, 3 );
  Eigen::Map< Matrix9dRowMajor >( tangents.Q_ApAm.data() ) = Kcond.block< 9, 9 >( 3, 12 );

  Eigen::Map< Matrix9x3RowMajor >( tangents.Q_Amw.data() ) = Kcond.block< 9, 3 >( 12, 0 );
  Eigen::Map< Matrix9dRowMajor >( tangents.Q_AmAp.data() ) = Kcond.block< 9, 9 >( 12, 3 );
  Eigen::Map< Matrix9dRowMajor >( tangents.Q_AmAm.data() ) = Kcond.block< 9, 9 >( 12, 12 );
}

void MarmotEquilibratedXInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int nStateVars )
{
  if ( !topMaterial || !bottomMaterial ) {
    for ( int i = 0; i < nStateVars; ++i ) {
      stateVars[i] = 0.0;
    }
    return;
  }

  topMaterial->initializeYourself( stateLayout.getPtr( stateVars, "topMaterialStateVars" ),
                                   topMaterial->getNumberOfRequiredStateVars() );
  bottomMaterial->initializeYourself( stateLayout.getPtr( stateVars, "bottomMaterialStateVars" ),
                                      bottomMaterial->getNumberOfRequiredStateVars() );

  double* normalGradientJump = stateLayout.getPtr( stateVars, "normalGradientJump" );
  normalGradientJump[0]      = 0.0;
  normalGradientJump[1]      = 0.0;
  normalGradientJump[2]      = 0.0;
}

double MarmotEquilibratedXInterfaceMaterialHypoElastic::getDensity()
{
  if ( !topMaterial ) {
    return -1;
  }

  return topMaterial->getDensity( nullptr );
}
