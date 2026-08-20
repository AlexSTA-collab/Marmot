#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <mutex>
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
    bool                  stateChanged = false; // the material updated its state
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

    // An elastic step leaves the internal variables untouched. Material-agnostic, but
    // NOT a plasticity test on its own: a viscoelastic material updates its Kelvin
    // strains every increment. The degeneracy gate below supplies that distinction.
    trial.stateChanged = !std::equal( trial.trialStateVars.begin(), trial.trialStateVars.end(), oldStateVars );

    return trial;
  }

  // ---------------------------------------------------------------------------
  // Gauge on the undetermined component of g.
  //
  // Once both faces flow, the averaged acoustic tensor Qbar loses its transverse
  // stiffness -- lambda_min = G H / (3G + H) -- and the component of g along
  //
  //     m_i = 1/2 ( N+_ij + N-_ij ) n_j ,      N = dev(sigma)/||dev(sigma)||
  //
  // is undetermined. The local solve is then bistable: it converges equally well to two
  // states differing by a jump along m, and alternates between them from one global
  // iteration to the next. The traction is the same on both branches; the coupling blocks
  // of the condensation, assembled from the two faces separately, are not, so the global
  // residual becomes discontinuous and the iteration stalls.
  //
  // Three details matter, each established by measurement:
  //   - the TANGENTIAL part of m is constrained, so the constraint is purely deviatoric
  //     and never touches the volumetric response;
  //   - the direction comes from the material, not from a singular vector, whose sign is
  //     arbitrary and would make the gauge itself two-valued;
  //   - the constraint SET is fixed once per increment. Re-decided every iteration, a
  //     handful of points of several thousand change state and hold the global Newton in
  //     a limit cycle. This is a property of any binary criterion, not of the measure.
  // ---------------------------------------------------------------------------
  bool envFlagOn( const char* name, bool fallback )
  {
    const char* e = std::getenv( name );
    if ( !e )
      return fallback;
    return std::string( e ) != "0";
  }

  bool gaugeEnabled()
  {
    static const bool v = envFlagOn( "MARMOT_IFACE_GAUGE", true );
    return v;
  }
  bool gaugeFreeze()
  {
    static const bool v = envFlagOn( "MARMOT_IFACE_GAUGE_FREEZE", true );
    return v;
  }

  // Elastic and viscoelastic tangents sit at sigma_min/sigma_max = 0.2857; near-perfect
  // plasticity at 1e-7. Six decades separate them, so this is a separation rather than a
  // tuned constant.
  double gaugeDegeneracyThreshold()
  {
    static const double v = []() {
      const char* e = std::getenv( "MARMOT_IFACE_GAUGE_TAU" );
      return e ? std::atof( e ) : 1.0e-6;
    }();
    return v;
  }

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

  // first decision taken in an increment wins, for the whole increment and across cutbacks
  bool frozenGaugeDecision( const double* stateVars, double timeOld, bool proposed )
  {
    const auto                                                   key = reinterpret_cast< std::uintptr_t >( stateVars );
    static std::mutex                                            mutex;
    static std::map< std::uintptr_t, std::pair< double, bool > > decisions;
    std::lock_guard< std::mutex >                                lock( mutex );
    auto                                                         it = decisions.find( key );
    if ( it != decisions.end() && it->second.first == timeOld )
      return it->second.second;
    decisions[key] = { timeOld, proposed };
    return proposed;
  }

  // Permanent monitors, behind MARMOT_IFACE_MONITOR.
  //   -(g.n)/tr D  predicts 1.000000 -- the layer's near-incompressibility expressed
  //   inside the condensation, and the best-conditioned check on the whole reduction.
  void monitor( double gDotN, double trD, bool gaugeActive )
  {
    if ( !std::getenv( "MARMOT_IFACE_MONITOR" ) )
      return;
    static std::mutex             mutex;
    static long                   calls = 0, active = 0;
    std::lock_guard< std::mutex > lock( mutex );
    ++calls;
    if ( gaugeActive )
      ++active;
    if ( calls % 100000 == 0 ) {
      std::printf( "  [iface] calls=%ld gaugeActive=%ld (%.1f%%)  -(g.n)/trD=%.6f\n",
                   calls,
                   active,
                   100.0 * double( active ) / double( calls ),
                   std::abs( trD ) > 0.0 ? -gDotN / trD : 0.0 );
      std::fflush( stdout );
    }
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
  // Each face carries its own committed Cauchy stress. Recovering them instead from the
  // single generalized force is exact only where t+ == t- holds identically; any residual
  // splits the two faces and compounds over increments.
  stateLayout.add( "stressPlus", 6 );
  stateLayout.add( "stressMinus", 6 );
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

  // The local solve plays the same role for the layer that a return map plays for a
  // material point: it makes the condensed state consistent. It is kept.
  //
  // What is NOT kept is the demand that it reach traction continuity exactly. The
  // reduction imposes the LINEARISED condition [dt_i] = 0, one solve with Qbar, which is
  // invertible. Near perfect plasticity the NONLINEAR condition t+(z) == t-(z) has no
  // solution at all: both faces sit on the yield surface, their deviatoric tractions stop
  // responding to z, and an initial deviatoric difference cannot be removed by any z. The
  // iteration therefore floors at a small residual instead of converging, and the best
  // iterate is accepted; the global Newton carries the remaining imbalance.
  constexpr int    maxLocalIterations = 40;
  constexpr double localTolerance     = 1.0e-10;

  Vector3d z = Eigen::Map< const Vector3d >( normalGradientJump );

  LocalEvaluation current = evaluateLocal( z );

  // --- gauge: is the layer degenerate here, and in which direction? ---
  bool gaugeActive = gaugeEnabled() && current.plus.stateChanged && current.minus.stateChanged;
  if ( gaugeActive ) {
    const Matrix3dRowMajor              QEntry = 0.5 * ( Rt * current.plus.CFull * Bn + Rt * current.minus.CFull * Bn );
    const Eigen::Matrix3d               QEntryDense = QEntry;
    Eigen::JacobiSVD< Eigen::Matrix3d > svd( QEntryDense );
    const Vector3d                      sv = svd.singularValues();
    gaugeActive                            = sv( 0 ) > 0.0 && sv( 2 ) / sv( 0 ) < gaugeDegeneracyThreshold();
  }
  if ( gaugeFreeze() )
    gaugeActive = frozenGaugeDecision( state.stateVars, timeIncrement.timeOld, gaugeActive );

  Vector3d                      gaugeDirection = Vector3d::Zero();
  Eigen::Matrix< double, 3, 2 > gaugeBasis     = Eigen::Matrix< double, 3, 2 >::Zero();
  if ( gaugeActive ) {
    const Matrix3dRowMajor flowPlus  = unitDeviator( current.plus.stress );
    const Matrix3dRowMajor flowMinus = unitDeviator( current.minus.stress );
    Vector3d               direction = Matrix3dRowMajor( 0.5 * ( flowPlus + flowMinus ) ) * normal;
    if ( direction.norm() <= 1.0e-10 )
      direction = flowPlus * normal;               // exactly antipodal faces
    direction -= direction.dot( normal ) * normal; // tangential: constraint is deviatoric
    if ( direction.norm() > 1.0e-12 ) {
      gaugeDirection       = direction.normalized();
      const Vector3d trial = ( std::abs( gaugeDirection[0] ) < 0.9 ) ? Vector3d::UnitX() : Vector3d::UnitY();
      gaugeBasis.col( 0 )  = ( trial - trial.dot( gaugeDirection ) * gaugeDirection ).normalized();
      gaugeBasis.col( 1 )  = gaugeDirection.cross( gaugeBasis.col( 0 ) ).normalized();
      z -= gaugeDirection.dot( z ) * gaugeDirection; // start on the gauge surface
      current = evaluateLocal( z );
    }
    else {
      gaugeActive = false;
    }
  }

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

    // Backtracking: a full step can overshoot when a trial state crosses a yield-surface
    // kink, so only a step that reduces the residual is accepted.
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

  // Commit: persist the converged internal variable (warm start for the
  // next increment) and the two sides' updated material state.
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

  // The tangent must condense through the subspace the solve used, or Kcond stops being
  // the derivative of the residual the element returns.
  Matrix3x21RowMajor QavgInvRx;
  if ( gaugeActive ) {
    const Eigen::Matrix2d restrictedQ = gaugeBasis.transpose() * Qavg * gaugeBasis;
    QavgInvRx = gaugeBasis * restrictedQ.fullPivLu().solve( ( gaugeBasis.transpose() * rx ).eval() );
  }
  else {
    QavgInvRx = Qavg.fullPivLu().solve( rx );
  }

  {
    const Matrix3dRowMajor D = APlus - AMinus;
    monitor( z.dot( normal ), D.trace(), gaugeActive );
  }
  const Matrix21dRowMajor Kcond = px - pz * QavgInvRx;

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
