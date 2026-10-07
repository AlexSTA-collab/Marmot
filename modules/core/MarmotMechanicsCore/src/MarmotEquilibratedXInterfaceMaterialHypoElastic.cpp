#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"

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
   * @brief Ratio of the smallest to the largest singular value of the averaged acoustic tensor, below which the
   * layer counts as degenerate.
   *
   * Elastic and viscoelastic tangents sit at a ratio of 0.2857, near-perfect plasticity at 1e-7. Six decades
   * separate them, so the threshold is a separation and not a tuned constant.
   */
  constexpr double gaugeDegeneracyThreshold = 1.0e-6;

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
      throw std::invalid_argument( "MarmotEquilibratedXInterfaceMaterialHypoElastic: the top-bottom connector must "
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
    Matrix3dRowMajor      stress;               ///< updated Cauchy stress
    Matrix9dRowMajor      CFull;                ///< tangent d(sigma_ij)/d(G_kl) with respect to the gradient
    std::vector< double > trialStateVars;       ///< updated state variables of the base material
    bool                  stateChanged = false; ///< the material updated its state variables
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

    // An elastic step leaves the internal variables untouched. This is material-agnostic, but NOT a plasticity
    // test on its own: a viscoelastic material updates its Kelvin strains every increment. The degeneracy
    // gate of the gauge supplies that distinction.
    trial.stateChanged = !std::equal( trial.trialStateVars.begin(), trial.trialStateVars.end(), oldStateVars );

    return trial;
  }

  /**
   * @brief Unit deviator of a stress.
   *
   * @param[in] sigma stress
   * @return the deviator divided by its norm, or zero for a vanishing deviator
   */
  Matrix3dRowMajor unitDeviator( const Matrix3dRowMajor& sigma )
  {
    Matrix3dRowMajor d = sigma;
    const double     m = sigma.trace() / 3.0;
    d( 0, 0 ) -= m;
    d( 1, 1 ) -= m;
    d( 2, 2 ) -= m;
    const double norm = d.norm();
    return Matrix3dRowMajor( norm > 0.0 ? Matrix3dRowMajor( d / norm ) : Matrix3dRowMajor( Matrix3dRowMajor::Zero() ) );
  }

} // namespace

MarmotEquilibratedXInterfaceMaterialHypoElastic::MarmotEquilibratedXInterfaceMaterialHypoElastic(
  const std::string& materialName,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : materialProperties( matProperties_ ),
    nMaterialProperties( nMaterialProperties_ ),
    h( checkedMaterialProperty( matProperties_, nMaterialProperties_, 2 ) ),
    materialNumber( materialNumber_ )
{
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
  // Each face carries its own committed Cauchy stress. Recovering them instead from the single generalized force
  // is exact only where t+ == t- holds identically; any residual splits the two faces and compounds over increments.
  stateLayout.add( "stressPlus", 6 );
  stateLayout.add( "stressMinus", 6 );
  stateLayout.finalize();
}

void MarmotEquilibratedXInterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  topMaterial->setCharacteristicElementLength( length );
  bottomMaterial->setCharacteristicElementLength( length );
}

void MarmotEquilibratedXInterfaceMaterialHypoElastic::computeStress( State&               state,
                                                                     Tangents&            tangents,
                                                                     const Deformation&   deformation,
                                                                     const TimeIncrement& timeIncrement )
{
  const Eigen::Map< const Vector3d > normalMap( deformation.normal.data() );
  const Eigen::Map< const Vector3d > separationMap( deformation.separationVector.data() );

  Vector3d     normal     = normalMap;
  const double normalNorm = normal.norm();
  if ( normalNorm <= geometryTolerance ) {
    throw std::invalid_argument( "MarmotEquilibratedXInterfaceMaterialHypoElastic: interface normal is zero." );
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

  // Recover both sides' OLD (committed, at the beginning of the increment) Cauchy stress from the persisted state.
  double* topStateVars       = stateLayout.getPtr( state.stateVars, "topMaterialStateVars" );
  double* bottomStateVars    = stateLayout.getPtr( state.stateVars, "bottomMaterialStateVars" );
  double* normalGradientJump = stateLayout.getPtr( state.stateVars, "normalGradientJump" );
  double* stressPlusStore    = stateLayout.getPtr( state.stateVars, "stressPlus" );
  double* stressMinusStore   = stateLayout.getPtr( state.stateVars, "stressMinus" );

  const Vector6d        stressPlusVoigtOld     = Eigen::Map< const Vector6d >( stressPlusStore );
  const Vector6d        stressMinusVoigtOld    = Eigen::Map< const Vector6d >( stressMinusStore );
  const Eigen::Matrix3d stressCurrentPlusEigen = ContinuumMechanics::VoigtNotation::voigtToStress( stressPlusVoigtOld );
  const Eigen::Matrix3d stressCurrentMinusEigen = ContinuumMechanics::VoigtNotation::voigtToStress(
    stressMinusVoigtOld );
  const Matrix3dRowMajor stressCurrentPlus  = stressCurrentPlusEigen;
  const Matrix3dRowMajor stressCurrentMinus = stressCurrentMinusEigen;
  const int              nTopStateVars      = topMaterial->getNumberOfRequiredStateVars();
  const int              nBottomStateVars   = bottomMaterial->getNumberOfRequiredStateVars();

  // the base materials expect the current time, i.e. the time at the end of the increment
  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld + timeIncrement.dT, timeIncrement.dT };

  // Local Newton solve for the normal-gradient jump z, enforcing sigma+(z) n == sigma-(z) n. Every trial
  // evaluation restarts from the OLD committed state (via the scratch copy of evaluateSideTrial), so that the
  // state is never accumulated across local iterations.
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

  // The local solve plays the same role for the layer that a return map plays for a material point: it makes the
  // condensed state consistent. It is kept.
  //
  // What is NOT kept is the demand that it reach traction continuity exactly. The reduction imposes the
  // LINEARISED condition [dt_i] = 0, one solve with Qbar, which is invertible. Near perfect plasticity the
  // NONLINEAR condition t+(z) == t-(z) has no solution at all: both faces sit on the yield surface, their
  // deviatoric tractions stop responding to z, and an initial deviatoric difference cannot be removed by any z.
  // The iteration therefore floors at a small residual instead of converging, and the best iterate is accepted;
  // the global Newton iteration carries the remaining imbalance.
  Vector3d z = Eigen::Map< const Vector3d >( normalGradientJump );

  LocalEvaluation current = evaluateLocal( z );

  // Gauge on the undetermined component of g.
  //
  // Once both faces flow, the averaged acoustic tensor Qbar loses its transverse stiffness, and the component of g
  // along
  //
  //     m_i = 1/2 ( N+_ij + N-_ij ) n_j ,      N = dev(sigma) / ||dev(sigma)||
  //
  // is undetermined. The local solve is then bistable: it converges equally well to two states that differ by a
  // jump along m, and alternates between them from one global iteration to the next. The traction is the same on
  // both branches; the coupling blocks of the condensation, assembled from the two faces separately, are not, so
  // the global residual becomes discontinuous and the iteration stalls.
  //
  // Three details matter:
  //   - the TANGENTIAL part of m is constrained, so the constraint is purely deviatoric and never touches the
  //     volumetric response;
  //   - the direction comes from the material, not from a singular vector, whose sign is arbitrary and would
  //     make the gauge itself two-valued;
  //   - the constraint SET is fixed once per increment, see IncrementMemory. Decided anew in every iteration, a
  //     handful of points of several thousand change state and hold the global Newton iteration in a limit cycle.
  bool gaugeActive = current.plus.stateChanged && current.minus.stateChanged;
  if ( gaugeActive ) {
    const Matrix3dRowMajor              QEntry = 0.5 * ( Rt * current.plus.CFull * Bn + Rt * current.minus.CFull * Bn );
    const Eigen::Matrix3d               QEntryDense = QEntry;
    Eigen::JacobiSVD< Eigen::Matrix3d > svd( QEntryDense );
    const Vector3d                      sv = svd.singularValues();
    gaugeActive                            = sv( 0 ) > 0.0 && sv( 2 ) / sv( 0 ) < gaugeDegeneracyThreshold;
  }
  if ( state.incrementMemory ) {
    gaugeActive = state.incrementMemory->decide( timeIncrement.timeOld, gaugeActive );
  }

  Eigen::Matrix< double, 3, 2 > gaugeBasis = Eigen::Matrix< double, 3, 2 >::Zero();
  if ( gaugeActive ) {
    const Matrix3dRowMajor flowPlus  = unitDeviator( current.plus.stress );
    const Matrix3dRowMajor flowMinus = unitDeviator( current.minus.stress );
    Vector3d               direction = Matrix3dRowMajor( 0.5 * ( flowPlus + flowMinus ) ) * normal;
    if ( direction.norm() <= 1.0e-10 )
      direction = flowPlus * normal;               // exactly antipodal faces
    direction -= direction.dot( normal ) * normal; // tangential: the constraint is deviatoric
    if ( direction.norm() > geometryTolerance ) {
      const Vector3d gaugeDirection = direction.normalized();
      const Vector3d trial          = ( std::abs( gaugeDirection[0] ) < 0.9 ) ? Vector3d::UnitX() : Vector3d::UnitY();
      gaugeBasis.col( 0 )           = ( trial - trial.dot( gaugeDirection ) * gaugeDirection ).normalized();
      gaugeBasis.col( 1 )           = gaugeDirection.cross( gaugeBasis.col( 0 ) ).normalized();
      z -= gaugeDirection.dot( z ) * gaugeDirection; // start on the gauge surface
      current = evaluateLocal( z );
    }
    else {
      gaugeActive = false;
    }
  }

  // Norm of the part of the residual that can be removed by the iteration.
  auto solvableNorm = [&]( const LocalEvaluation& e ) {
    return gaugeActive
             ? ( gaugeBasis.transpose() * e.r ).norm() /
                 std::max( 1.0, std::max( ( e.plus.stress * normal ).norm(), ( e.minus.stress * normal ).norm() ) )
             : e.rNorm;
  };

  Vector3d        bestZ    = z;
  LocalEvaluation bestEval = current;

  for ( int iteration = 0; iteration < maxLocalIterations && solvableNorm( current ) > localTolerance; ++iteration ) {
    const Matrix3dRowMajor Qplus  = Rt * current.plus.CFull * Bn;
    const Matrix3dRowMajor Qminus = Rt * current.minus.CFull * Bn;
    const Matrix3dRowMajor Qavg   = 0.5 * ( Qplus + Qminus );

    Vector3d dz;
    if ( gaugeActive ) {
      // two directions the layer still resists; the near-null one is never inverted
      const Eigen::Matrix2d restrictedQ = gaugeBasis.transpose() * Qavg * gaugeBasis;
      dz = gaugeBasis * restrictedQ.fullPivLu().solve( Eigen::Vector2d( -( gaugeBasis.transpose() * current.r ) ) );
    }
    else {
      dz = Qavg.fullPivLu().solve( -current.r );
    }

    // Backtracking: a full step can overshoot when a trial state crosses a yield-surface kink, so only a step that
    // reduces the residual is accepted.
    double          alpha = 1.0;
    LocalEvaluation candidate;
    bool            accepted = false;
    for ( int lineSearchIter = 0; lineSearchIter < 12; ++lineSearchIter ) {
      candidate = evaluateLocal( z + alpha * dz );
      if ( solvableNorm( candidate ) <= localTolerance || solvableNorm( candidate ) < solvableNorm( current ) ) {
        accepted = true;
        break;
      }
      alpha *= 0.5;
    }
    if ( !accepted )
      break;

    z += alpha * dz;
    current = candidate;
    if ( solvableNorm( current ) < solvableNorm( bestEval ) ) {
      bestEval = current;
      bestZ    = z;
    }
  }

  z       = bestZ;
  current = bestEval;

  const SideTrial& plusTrial  = current.plus;
  const SideTrial& minusTrial = current.minus;

  // Commit: persist the converged internal variable (warm start for the next increment) and the updated material
  // state of both sides.
  Eigen::Map< Vector3d > normalGradientJumpMap( normalGradientJump );
  normalGradientJumpMap = z;
  std::copy( plusTrial.trialStateVars.begin(), plusTrial.trialStateVars.end(), topStateVars );
  std::copy( minusTrial.trialStateVars.begin(), minusTrial.trialStateVars.end(), bottomStateVars );

  // Commit each face's own Cauchy stress, so the next increment restarts from it directly.
  {
    const Eigen::Matrix3d  plusSym( 0.5 * ( plusTrial.stress + plusTrial.stress.transpose() ) );
    const Eigen::Matrix3d  minusSym( 0.5 * ( minusTrial.stress + minusTrial.stress.transpose() ) );
    Eigen::Map< Vector6d > plusStoreMap( stressPlusStore );
    Eigen::Map< Vector6d > minusStoreMap( stressMinusStore );
    plusStoreMap  = ContinuumMechanics::VoigtNotation::stressToVoigt( plusSym );
    minusStoreMap = ContinuumMechanics::VoigtNotation::stressToVoigt( minusSym );
  }

  // Purely geometric operators, independent of z and of the material state.
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

  // The tangent must condense through the subspace the solve used, or Kcond stops being the derivative of the
  // residual the element returns.
  Matrix3x21RowMajor QavgInvRx;
  if ( gaugeActive ) {
    const Eigen::Matrix2d restrictedQ = gaugeBasis.transpose() * Qavg * gaugeBasis;
    QavgInvRx = gaugeBasis * restrictedQ.fullPivLu().solve( ( gaugeBasis.transpose() * rx ).eval() );
  }
  else {
    QavgInvRx = Qavg.fullPivLu().solve( rx );
  }

  const Matrix21dRowMajor Kcond = px - pz * QavgInvRx;

  const Eigen::Map< const Vector9d > vecSigmaPlus( plusTrial.stress.data() );
  const Eigen::Map< const Vector9d > vecSigmaMinus( minusTrial.stress.data() );
  const Vector21d p = sideThickness * ( BxPlus.transpose() * vecSigmaPlus + BxMinus.transpose() * vecSigmaMinus );

  state.generalizedForce   = Tensor3d( p.data() );
  state.surfaceStressPlus  = Tensor33d( p.data() + 3 );
  state.surfaceStressMinus = Tensor33d( p.data() + 12 );

  mapEigenToFastor( tangents.Q_ww )  = Kcond.block< 3, 3 >( 0, 0 );
  mapEigenToFastor( tangents.Q_wAp ) = Kcond.block< 3, 9 >( 0, 3 );
  mapEigenToFastor( tangents.Q_wAm ) = Kcond.block< 3, 9 >( 0, 12 );

  // the blocks of the stress of a face with respect to the jump are 9x3 matrices, stored as third order tensors
  Eigen::Map< Matrix9x3RowMajor >( tangents.Q_Apw.data() ) = Kcond.block< 9, 3 >( 3, 0 );
  mapEigenToFastor( tangents.Q_ApAp )                      = Kcond.block< 9, 9 >( 3, 3 );
  mapEigenToFastor( tangents.Q_ApAm )                      = Kcond.block< 9, 9 >( 3, 12 );

  Eigen::Map< Matrix9x3RowMajor >( tangents.Q_Amw.data() ) = Kcond.block< 9, 3 >( 12, 0 );
  mapEigenToFastor( tangents.Q_AmAp )                      = Kcond.block< 9, 9 >( 12, 3 );
  mapEigenToFastor( tangents.Q_AmAm )                      = Kcond.block< 9, 9 >( 12, 12 );
}

void MarmotEquilibratedXInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int nStateVars )
{
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
  return topMaterial->getDensity( nullptr );
}
