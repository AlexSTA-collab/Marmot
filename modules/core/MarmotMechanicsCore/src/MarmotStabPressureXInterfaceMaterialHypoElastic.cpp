#include "Marmot/MarmotStabPressureXInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtils.h"
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
  using Matrix3x21RowMajor = Eigen::Matrix< double, 3, 21, Eigen::RowMajor >;
  using Matrix21x3RowMajor = Eigen::Matrix< double, 21, 3, Eigen::RowMajor >;
  using Matrix21dRowMajor  = Eigen::Matrix< double, 21, 21, Eigen::RowMajor >;

  /// tolerance of the geometry of the connector and of the interface normal
  constexpr double geometryTolerance = 1.0e-12;

  /**
   * @brief Maximum number of iterations of the local Newton solve for the normal-gradient jump.
   */
  constexpr int maxLocalIterations = 40;

  /**
   * @brief Tolerance of the relative traction mismatch of the local Newton solve.
   */
  constexpr double localTolerance = 1.0e-10;

  /**
   * @brief Flat row-major index of the entry (i, j) of a 3x3 tensor.
   *
   * @param[in] i row
   * @param[in] j column
   * @return the flat index
   */
  constexpr int flatIndex( int i, int j )
  {
    return 3 * i + j;
  }

  /**
   * @brief Geometry of the connector between the two faces of the interface layer.
   */
  struct InterfaceGeometry {
    double   normalSeparation;     ///< component of the connector along the normal
    Vector3d tangentialSeparation; ///< component of the connector in the tangent plane
  };

  /**
   * @brief Split the connector between the two faces into its normal and its tangential part.
   *
   * @param[in] normal unit normal of the interface
   * @param[in] separationVector connector from the bottom to the top face, zero for coincident faces
   * @param[in] constitutiveThickness thickness h, used as the normal separation of coincident faces
   * @return the normal separation and the tangential separation
   *
   * @throws std::invalid_argument if the connector has no positive normal component
   */
  InterfaceGeometry evaluateInterfaceGeometry( const Vector3d& normal,
                                               const Vector3d& separationVector,
                                               double          constitutiveThickness )
  {
    if ( separationVector.norm() <= geometryTolerance ) {
      return { constitutiveThickness, Vector3d::Zero() };
    }

    const double normalSeparation = separationVector.dot( normal );
    if ( normalSeparation <= geometryTolerance ) {
      throw std::invalid_argument( "MarmotStabPressureXInterfaceMaterialHypoElastic: the top-bottom connector must "
                                   "have a positive normal component." );
    }

    return { normalSeparation, separationVector - normalSeparation * normal };
  }

  /**
   * @brief Map a 3-vector g to vec(g (x) n) / normalSeparation.
   *
   * With a normal separation of 1 this is B_n, with the normal separation of the connector it is B_w.
   *
   * @param[in] normal unit normal of the interface
   * @param[in] normalSeparation normal separation of the two faces
   * @return the 9x3 map
   */
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

  /**
   * @brief Map vec(A) to vec(A) - vec((A d_tau) (x) n) / normalSeparation, i.e. I_A + B_tau.
   *
   * @param[in] normal unit normal of the interface
   * @param[in] tangentialSeparation tangential part of the connector
   * @param[in] normalSeparation normal separation of the two faces
   * @return the 9x9 map
   */
  Matrix9dRowMajor surfaceToGradient( const Vector3d& normal,
                                      const Vector3d& tangentialSeparation,
                                      double          normalSeparation )
  {
    Matrix9dRowMajor map = Matrix9dRowMajor::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        const int row = flatIndex( k, l );

        for ( int b = 0; b < 3; ++b ) {
          const double identityPart     = l == b ? 1.0 : 0.0;
          map( row, flatIndex( k, b ) ) = identityPart - tangentialSeparation( b ) * normal( l ) / normalSeparation;
        }
      }
    }

    return map;
  }

  /**
   * @brief Map vec(sigma) to (sideThickness / normalSeparation) sigma n.
   *
   * With a side thickness and a normal separation of 1 this is R_t.
   *
   * @param[in] normal unit normal of the interface
   * @param[in] sideThickness thickness of one face
   * @param[in] normalSeparation normal separation of the two faces
   * @return the 3x9 map
   */
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

  /**
   * @brief Trial constitutive response of one face, evaluated from a scratch copy of the old state.
   *
   * The persistent state variables are never modified by a trial.
   */
  struct SideTrial {
    Matrix3dRowMajor      stress;         ///< updated Cauchy stress
    Matrix9dRowMajor      CFull;          ///< tangent d(sigma_ij)/d(G_kl) with respect to the gradient
    std::vector< double > trialStateVars; ///< updated state variables of the base material
  };

  /**
   * @brief Evaluate the trial response of one face.
   *
   * @param[in,out] material base material of the face
   * @param[in] oldStateVars committed state variables of the face at the beginning of the increment
   * @param[in] nStateVars number of state variables of the face
   * @param[in] stressCurrent committed Cauchy stress of the face
   * @param[in] strainIncrementVoigt strain increment of the face in Voigt notation
   * @param[in] timeInfo time and time increment
   * @return the trial response
   */
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

    // the stress depends on the symmetric part of the gradient only, so d(sigma_ij)/d(G_kl) is the stiffness tensor
    const Tensor3333d stiffness = ContinuumMechanics::VoigtNotation::voigtToStiffnessFastor( tangentVoigt );
    trial.CFull                 = mapEigenToFastor( stiffness );

    return trial;
  }

  /**
   * @brief Projector on the deviator, acting on vec(sigma).
   *
   * @return the 9x9 matrix of sigma -> sigma - tr(sigma) / 3 I
   */
  Matrix9dRowMajor deviatoricProjector()
  {
    Matrix9dRowMajor projector = Matrix9dRowMajor::Identity();

    for ( int i = 0; i < 3; ++i ) {
      for ( int k = 0; k < 3; ++k ) {
        projector( flatIndex( i, i ), flatIndex( k, k ) ) -= 1.0 / 3.0;
      }
    }

    return projector;
  }

} // namespace

MarmotStabPressureXInterfaceMaterialHypoElastic::MarmotStabPressureXInterfaceMaterialHypoElastic(
  const std::string& materialName,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : materialProperties( matProperties_ ), nMaterialProperties( nMaterialProperties_ ), materialNumber( materialNumber_ )
{
  h = checkedMaterialProperty( matProperties_, nMaterialProperties_, 2 );

  if ( h <= 0.0 ) {
    throw std::invalid_argument( "MarmotStabPressureXInterfaceMaterialHypoElastic requires h > 0." );
  }

  const double E  = materialProperties[0];
  const double nu = materialProperties[1];

  if ( std::abs( 1.0 - 2.0 * nu ) < 1.0e-12 ) {
    throw std::invalid_argument( "MarmotStabPressureXInterfaceMaterialHypoElastic: nu = 0.5 (K infinite) is not "
                                 "supported." );
  }

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
  // Each face carries its own committed Cauchy stress, see MarmotEquilibratedXInterfaceMaterialHypoElastic.
  stateLayout.add( "stressPlus", 6 );
  stateLayout.add( "stressMinus", 6 );
  stateLayout.add( "pressure", 1 );
  stateLayout.finalize();
}

void MarmotStabPressureXInterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  topMaterial->setCharacteristicElementLength( length );
  bottomMaterial->setCharacteristicElementLength( length );
}

void MarmotStabPressureXInterfaceMaterialHypoElastic::computeStress( State&               state,
                                                                     Tangents&            tangents,
                                                                     const Deformation&   deformation,
                                                                     const TimeIncrement& timeIncrement )
{
  const Eigen::Map< const Vector3d > normalMap( deformation.normal.data() );
  const Eigen::Map< const Vector3d > separationMap( deformation.separationVector.data() );

  Vector3d     normal     = normalMap;
  const double normalNorm = normal.norm();
  if ( normalNorm <= geometryTolerance ) {
    throw std::invalid_argument( "MarmotStabPressureXInterfaceMaterialHypoElastic: interface normal is zero." );
  }
  normal /= normalNorm;

  const InterfaceGeometry geometry      = evaluateInterfaceGeometry( normal, separationMap, h );
  const double            ell           = geometry.normalSeparation;
  const Vector3d&         dTangential   = geometry.tangentialSeparation;
  const double            sideThickness = 0.5 * h;

  const Eigen::Map< const Vector3d > dUPlus( deformation.dU.data() );
  const Eigen::Map< const Vector3d > dUMinus( deformation.dU.data() + 3 );

  const Vector3d w = dUPlus - dUMinus;

  const Eigen::Map< const Matrix3dRowMajor > APlus( deformation.dSurfaceStrain.data() );
  const Eigen::Map< const Matrix3dRowMajor > AMinus( deformation.dSurfaceStrain.data() + 9 );
  const Matrix3dRowMajor                     ABar = 0.5 * ( APlus + AMinus );

  const Vector3d gBar = ( w - ABar * dTangential ) / ell;

  double*   topStateVars       = stateLayout.getPtr( state.stateVars, "topMaterialStateVars" );
  double*   bottomStateVars    = stateLayout.getPtr( state.stateVars, "bottomMaterialStateVars" );
  double*   normalGradientJump = stateLayout.getPtr( state.stateVars, "normalGradientJump" );
  double*   stressPlusStore    = stateLayout.getPtr( state.stateVars, "stressPlus" );
  double*   stressMinusStore   = stateLayout.getPtr( state.stateVars, "stressMinus" );
  double*   pressureStore      = stateLayout.getPtr( state.stateVars, "pressure" );
  const int nTopStateVars      = topMaterial->getNumberOfRequiredStateVars();
  const int nBottomStateVars   = bottomMaterial->getNumberOfRequiredStateVars();

  // total pressure at the end of the increment
  const double pressureTotal = pressureStore[0] + deformation.dPressure;

  const Matrix3dRowMajor stressCurrentPlus = Matrix3dRowMajor(
    ContinuumMechanics::VoigtNotation::voigtToStress( Vector6d( Eigen::Map< const Vector6d >( stressPlusStore ) ) ) );
  const Matrix3dRowMajor stressCurrentMinus = Matrix3dRowMajor(
    ContinuumMechanics::VoigtNotation::voigtToStress( Vector6d( Eigen::Map< const Vector6d >( stressMinusStore ) ) ) );

  // the base materials expect the current time, i.e. the time at the end of the increment
  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld + timeIncrement.dT, timeIncrement.dT };

  const Matrix3x9RowMajor Rt   = stressToForce( normal, 1.0, 1.0 );
  const Matrix9x3RowMajor Bn   = jumpToGradient( normal, 1.0 );
  const Matrix9dRowMajor  Pdev = deviatoricProjector();
  const Vector9d          eI   = Eigen::Map< const Vector9d >( Matrix3dRowMajor::Identity().eval().data() );

  struct LocalEvaluation {
    SideTrial plus;
    SideTrial minus;
    Vector3d  r;
    double    rNorm;
  };

  // The equilibrated stress is dev(sigma) - p I. The term -p I cancels in the difference of the two faces, so the
  // traction mismatch is that of the deviators.
  auto evaluateLocal = [&]( const Vector3d& zTrial ) -> LocalEvaluation {
    const Vector3d gPlus  = gBar + 0.5 * zTrial;
    const Vector3d gMinus = gBar - 0.5 * zTrial;

    const Matrix3dRowMajor GPlus  = APlus + gPlus * normal.transpose();
    const Matrix3dRowMajor GMinus = AMinus + gMinus * normal.transpose();

    const Vector6d strainIncrementPlusVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt(
      Eigen::Matrix3d( 0.5 * ( GPlus + GPlus.transpose() ) ) );
    const Vector6d strainIncrementMinusVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt(
      Eigen::Matrix3d( 0.5 * ( GMinus + GMinus.transpose() ) ) );

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

    const Vector9d deviatorPlus  = Pdev * Eigen::Map< const Vector9d >( evaluation.plus.stress.data() );
    const Vector9d deviatorMinus = Pdev * Eigen::Map< const Vector9d >( evaluation.minus.stress.data() );

    evaluation.r     = Rt * ( deviatorPlus - deviatorMinus );
    evaluation.rNorm = evaluation.r.norm() /
                       std::max( 1.0, std::max( ( Rt * deviatorPlus ).norm(), ( Rt * deviatorMinus ).norm() ) );
    return evaluation;
  };

  // z enters the strain of both faces as an INCREMENT, like gBar, so the solve starts from zero. The committed z of the
  // previous increment is not a starting guess for this one: at a load reversal it has the wrong sign, and at a zero
  // increment it is applied as a strain increment to faces that should not move, which can make the return mapping of
  // the base material fail before the first iteration. It is still stored for output.
  Vector3d        z         = Vector3d::Zero();
  LocalEvaluation current   = evaluateLocal( z );
  bool            converged = current.rNorm <= localTolerance;

  for ( int iteration = 0; iteration < maxLocalIterations && !converged; ++iteration ) {
    const Matrix3dRowMajor Qplus  = Rt * ( Pdev * current.plus.CFull ) * Bn;
    const Matrix3dRowMajor Qminus = Rt * ( Pdev * current.minus.CFull ) * Bn;
    const Matrix3dRowMajor Qavg   = 0.5 * ( Qplus + Qminus );
    const Vector3d         dz     = Qavg.fullPivLu().solve( -current.r );

    // Backtracking: only a step that reduces the residual is accepted.
    double          alpha = 1.0;
    LocalEvaluation candidate;
    bool            accepted = false;
    for ( int lineSearchIter = 0; lineSearchIter < 12; ++lineSearchIter ) {
      candidate = evaluateLocal( z + alpha * dz );
      if ( candidate.rNorm <= localTolerance || candidate.rNorm < current.rNorm ) {
        accepted = true;
        break;
      }
      alpha *= 0.5;
    }

    z += alpha * dz;
    current   = accepted ? candidate : evaluateLocal( z );
    converged = current.rNorm <= localTolerance;
  }

  if ( !converged ) {
    throw Marmot::StressUpdateFailed( "MarmotStabPressureXInterfaceMaterialHypoElastic: local traction-equilibrium "
                                      "solve did not converge." );
  }

  const SideTrial& plusTrial  = current.plus;
  const SideTrial& minusTrial = current.minus;

  // Commit the converged internal variable, the state of both faces, their stresses and the pressure.
  Eigen::Map< Vector3d > normalGradientJumpMap( normalGradientJump );
  normalGradientJumpMap = z;
  std::copy( plusTrial.trialStateVars.begin(), plusTrial.trialStateVars.end(), topStateVars );
  std::copy( minusTrial.trialStateVars.begin(), minusTrial.trialStateVars.end(), bottomStateVars );
  {
    const Eigen::Matrix3d  plusSym( 0.5 * ( plusTrial.stress + plusTrial.stress.transpose() ) );
    const Eigen::Matrix3d  minusSym( 0.5 * ( minusTrial.stress + minusTrial.stress.transpose() ) );
    Eigen::Map< Vector6d > plusStoreMap( stressPlusStore );
    Eigen::Map< Vector6d > minusStoreMap( stressMinusStore );
    plusStoreMap  = ContinuumMechanics::VoigtNotation::stressToVoigt( plusSym );
    minusStoreMap = ContinuumMechanics::VoigtNotation::stressToVoigt( minusSym );
  }
  pressureStore[0] = pressureTotal;

  // Purely geometric operators, independent of z and of the material state.
  const Matrix9x3RowMajor Bw                = jumpToGradient( normal, ell );
  const Matrix9dRowMajor  Asrf              = surfaceToGradient( normal, dTangential, ell );
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

  // tangents of the deviatoric stress
  const Matrix9dRowMajor CTplus  = Pdev * plusTrial.CFull;
  const Matrix9dRowMajor CTminus = Pdev * minusTrial.CFull;

  const Matrix3dRowMajor   Qavg = 0.5 * ( Rt * ( CTplus + CTminus ) * Bn );
  const Matrix3x21RowMajor rx   = Rt * ( CTplus * BxPlus - CTminus * BxMinus );

  const Matrix21dRowMajor px = sideThickness *
                               ( BxPlus.transpose() * CTplus * BxPlus + BxMinus.transpose() * CTminus * BxMinus );
  const Matrix21x3RowMajor pz = sideThickness *
                                ( BxPlus.transpose() * CTplus * BzPlus + BxMinus.transpose() * CTminus * BzMinus );

  const Matrix21dRowMajor Kcond = px - pz * Qavg.fullPivLu().solve( rx );

  // generalized stress from the equilibrated stress dev(sigma) - p I
  const Vector9d stressPlusTilde  = Pdev * Eigen::Map< const Vector9d >( plusTrial.stress.data() ) - pressureTotal * eI;
  const Vector9d stressMinusTilde = Pdev * Eigen::Map< const Vector9d >( minusTrial.stress.data() ) -
                                    pressureTotal * eI;
  const Vector21d p = sideThickness * ( BxPlus.transpose() * stressPlusTilde + BxMinus.transpose() * stressMinusTilde );

  state.generalizedForce   = Tensor3d( p.data() );
  state.surfaceStressPlus  = Tensor33d( p.data() + 3 );
  state.surfaceStressMinus = Tensor33d( p.data() + 12 );

  // Volumetric residual, independent of z. The law is incremental: mixing the strain increment with the total
  // pressure would underpredict p by the number of increments.
  const Matrix3dRowMajor GPlus  = APlus + ( gBar + 0.5 * z ) * normal.transpose();
  const Matrix3dRowMajor GMinus = AMinus + ( gBar - 0.5 * z ) * normal.transpose();
  state.volumetricResidual      = 0.5 * ( GPlus.trace() + GMinus.trace() ) + deformation.dPressure / bulkModulus;

  mapEigenToFastor( tangents.Q_ww )  = Kcond.block< 3, 3 >( 0, 0 );
  mapEigenToFastor( tangents.Q_wAp ) = Kcond.block< 3, 9 >( 0, 3 );
  mapEigenToFastor( tangents.Q_wAm ) = Kcond.block< 3, 9 >( 0, 12 );

  Eigen::Map< Matrix9x3RowMajor >( tangents.Q_Apw.data() ) = Kcond.block< 9, 3 >( 3, 0 );
  mapEigenToFastor( tangents.Q_ApAp )                      = Kcond.block< 9, 9 >( 3, 3 );
  mapEigenToFastor( tangents.Q_ApAm )                      = Kcond.block< 9, 9 >( 3, 12 );

  Eigen::Map< Matrix9x3RowMajor >( tangents.Q_Amw.data() ) = Kcond.block< 9, 3 >( 12, 0 );
  mapEigenToFastor( tangents.Q_AmAp )                      = Kcond.block< 9, 9 >( 12, 3 );
  mapEigenToFastor( tangents.Q_AmAm )                      = Kcond.block< 9, 9 >( 12, 12 );

  // d(y)/d(p) = -sideThickness ( BxPlus^T e + BxMinus^T e ), purely geometric because dz/dp = 0
  const Vector21d dydp = -sideThickness * ( BxPlus.transpose() * eI + BxMinus.transpose() * eI );
  Eigen::Map< Vector3d >( tangents.Q_wp.data() )  = dydp.segment< 3 >( 0 );
  Eigen::Map< Vector9d >( tangents.Q_App.data() ) = dydp.segment< 9 >( 3 );
  Eigen::Map< Vector9d >( tangents.Q_Amp.data() ) = dydp.segment< 9 >( 12 );

  // d(r_p)/d(x) = 1/2 e^T ( BxPlus + BxMinus ), purely geometric because d(r_p)/dz = 0
  const Eigen::Matrix< double, 1, 21 > drdx                           = 0.5 * ( eI.transpose() * ( BxPlus + BxMinus ) );
  Eigen::Map< Eigen::Matrix< double, 1, 3 > >( tangents.Q_pw.data() ) = drdx.block< 1, 3 >( 0, 0 );
  Eigen::Map< Eigen::Matrix< double, 1, 9 > >( tangents.Q_pAp.data() ) = drdx.block< 1, 9 >( 0, 3 );
  Eigen::Map< Eigen::Matrix< double, 1, 9 > >( tangents.Q_pAm.data() ) = drdx.block< 1, 9 >( 0, 12 );
  tangents.Q_pp                                                        = 1.0 / bulkModulus;
}

void MarmotStabPressureXInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int nStateVars )
{
  for ( int i = 0; i < nStateVars; ++i ) {
    stateVars[i] = 0.0;
  }

  topMaterial->initializeYourself( stateLayout.getPtr( stateVars, "topMaterialStateVars" ),
                                   topMaterial->getNumberOfRequiredStateVars() );
  bottomMaterial->initializeYourself( stateLayout.getPtr( stateVars, "bottomMaterialStateVars" ),
                                      bottomMaterial->getNumberOfRequiredStateVars() );
}

double MarmotStabPressureXInterfaceMaterialHypoElastic::getDensity()
{
  return topMaterial->getDensity( nullptr );
}
