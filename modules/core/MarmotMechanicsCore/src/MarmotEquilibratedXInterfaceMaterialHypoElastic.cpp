#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <mutex>
#include <set>
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
    bool                  stateChanged = false; // the material updated its state: inelastic step
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

    // Material-agnostic test for an inelastic step: an elastic step leaves the
    // internal variables untouched.
    trial.stateChanged = !std::equal( trial.trialStateVars.begin(), trial.trialStateVars.end(), oldStateVars );

    return trial;
  }

  // Relative truncation tolerance for the guarded inverse of the averaged
  // acoustic tensor. MARMOT_IFACE_QGUARD_RELTOL overrides it; a value of 0
  // retains every direction and reduces the guard to a plain pseudo-inverse,
  // which is the A/B control.
  double qGuardRelativeTolerance()
  {
    static const double value = []() {
      const char*  s = std::getenv( "MARMOT_IFACE_QGUARD_RELTOL" );
      const double v = s ? std::atof( s ) : 0.0; // truncation off unless asked for
      return v > 0.0 ? v : 0.0;
    }();
    return value;
  }

  // Which residual the local iteration converges on. "full" restores the original
  // criterion (the whole traction imbalance) while leaving the truncation active, so
  // the truncation and the projected criterion can be attributed separately.
  bool qGuardFullResidualConvergence()
  {
    static const bool value = []() {
      const char* s = std::getenv( "MARMOT_IFACE_QGUARD_CONVERGENCE" );
      return s && std::string( s ) == "full";
    }();
    return value;
  }

  // The gauge is on unless MARMOT_IFACE_GAUGE=0.
  bool gaugeEnabled()
  {
    static const bool value = []() {
      const char* s = std::getenv( "MARMOT_IFACE_GAUGE" );
      return !( s && std::string( s ) == "0" );
    }();
    return value;
  }

  // Forces the gauge onto the surface g.m = offset instead of g.m = 0. Used to test
  // whether the discarded component affects the reconstructed face states at all.
  double gaugeOffset()
  {
    static const double value = []() {
      const char* s = std::getenv( "MARMOT_IFACE_GAUGE_OFFSET" );
      return s ? std::atof( s ) : 0.0;
    }();
    return value;
  }

  // MARMOT_IFACE_GAUGE_TANGENTIAL=0 constrains m itself rather than its tangential
  // part, which is the earlier form of the gauge; used to attribute behaviour to the
  // tangential projection.
  bool gaugeTangential()
  {
    static const bool value = []() {
      const char* s = std::getenv( "MARMOT_IFACE_GAUGE_TANGENTIAL" );
      return !( s && std::string( s ) == "0" );
    }();
    return value;
  }

  // Degeneracy gate. A changing internal state is not evidence of degeneracy: a
  // viscoelastic material updates its Kelvin strains every increment while Qbar stays
  // positive definite at sigma_min/sigma_max = 0.2857. The gauge must engage on the
  // degeneracy itself. Elastic and viscoelastic sit at 0.2857, near-perfect plasticity at
  // 1e-7, so the separation is six decades and the gate is not a tuned quantity.
  double gaugeDegeneracyThreshold()
  {
    static const double value = []() {
      const char* s = std::getenv( "MARMOT_IFACE_GAUGE_TAU" );
      return s ? std::atof( s ) : 1.0e-6;
    }();
    return value;
  }

  bool gaugeFreeze()
  {
    static const bool value = []() {
      const char* s = std::getenv( "MARMOT_IFACE_GAUGE_FREEZE" );
      return s && std::string( s ) == "1";
    }();
    return value;
  }

  // first decision taken at a given increment wins, for the whole increment
  bool frozenGaugeDecision( unsigned long quadraturePointId, double timeOld, bool proposed )
  {
    static std::mutex                                           mutex;
    static std::map< unsigned long, std::pair< double, bool > > decisions;
    std::lock_guard< std::mutex >                               lock( mutex );
    auto                                                        it = decisions.find( quadraturePointId );
    if ( it != decisions.end() && it->second.first == timeOld )
      return it->second.second;
    decisions[quadraturePointId] = { timeOld, proposed };
    return proposed;
  }

  double gaugeSelfCheckDelta()
  {
    static const double value = []() {
      const char* s = std::getenv( "MARMOT_IFACE_GAUGE_SELFCHECK" );
      return s ? std::atof( s ) : 0.0;
    }();
    return value;
  }

  unsigned long quadraturePointIdOf( const double* stateVars )
  {
    return static_cast< unsigned long >( ( reinterpret_cast< std::uintptr_t >( stateVars ) >> 4 ) & 0xffffffffu );
  }

  // Unit deviatoric direction of a stress state: the plastic flow direction for J2.
  Matrix3dRowMajor unitDeviator( const Matrix3dRowMajor& sigma )
  {
    Matrix3dRowMajor deviator = sigma;
    const double     mean     = sigma.trace() / 3.0;
    deviator( 0, 0 ) -= mean;
    deviator( 1, 1 ) -= mean;
    deviator( 2, 2 ) -= mean;
    const double norm = deviator.norm();
    return Matrix3dRowMajor( norm > 0.0 ? Matrix3dRowMajor( deviator / norm )
                                        : Matrix3dRowMajor( Matrix3dRowMajor::Zero() ) );
  }

  bool criterionDetS()
  {
    static const bool value = []() {
      const char* e = std::getenv( "MARMOT_IFACE_CRITERION" );
      return e && std::string( e ) == "dets";
    }();
    return value;
  }

  // theta and thetabar read out of the algorithmic tangent:
  //   C:X = 2 mu theta X   for X deviatoric, symmetric, orthogonal to M
  //   M:C:M = 2 mu ( theta - thetabar )
  void probeAlgorithmicTheta( const Matrix3dRowMajor& sigma,
                              const Matrix9dRowMajor& C,
                              double                  mu,
                              double&                 theta,
                              double&                 thetaBar )
  {
    const Matrix3dRowMajor M = unitDeviator( sigma );
    Matrix3dRowMajor       X;
    X << 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0;
    X -= ( X.array() * M.array() ).sum() * M;
    const double                       xn = X.norm();
    const Eigen::Map< const Vector9d > vX( X.data() );
    const Eigen::Map< const Vector9d > vM( M.data() );
    theta    = xn > 1.0e-12 ? double( ( vX.transpose() * C * vX )( 0, 0 ) ) / ( 2.0 * mu * xn * xn ) : 1.0;
    thetaBar = theta - double( ( vM.transpose() * C * vM )( 0, 0 ) ) / ( 2.0 * mu );
  }

  // Sweep accounting: how often the gauge engages, whether the near-antipodal
  // fallback is ever exercised, and whether the tangential constraint really is
  // orthogonal to the normal.
  void reportGaugeSetup( double flowSumNorm, bool fallbackUsed, double tangentialDotNormal )
  {
    static std::mutex mutex;
    static long       activations = 0, fallbacks = 0;
    static double     minFlowSum = 1.0e300, maxTangentialDotNormal = 0.0;

    std::lock_guard< std::mutex > lock( mutex );
    ++activations;
    if ( fallbackUsed )
      ++fallbacks;
    minFlowSum             = std::min( minFlowSum, flowSumNorm );
    maxTangentialDotNormal = std::max( maxTangentialDotNormal, std::abs( tangentialDotNormal ) );

    if ( activations % 100000 == 0 ) {
      std::printf( "  [sweep] gaugeActivations=%ld  fallbacks=%ld  min||N+ + N-||=%.6e  max|m^T.n|=%.3e\n",
                   activations,
                   fallbacks,
                   minFlowSum,
                   maxTangentialDotNormal );
      std::fflush( stdout );
    }
  }

  // total calls, so activation can be reported as a fraction
  void reportCall()
  {
    static std::atomic< long > calls{ 0 };
    const long                 n = ++calls;
    if ( n % 500000 == 0 ) {
      std::printf( "  [sweep] materialCalls=%ld\n", n );
      std::fflush( stdout );
    }
  }

  // Continuous replacement for the binary activation. w = 1 well inside the degenerate
  // regime, w = 0 well outside, with a C1 ramp between, so no quadrature point can change
  // the local problem discontinuously. A single point toggling a binary criterion was
  // enough to hold the global Newton in a limit cycle.
  double gaugeWeight( double ratio )
  {
    static const double lo = []() {
      const char* s = std::getenv( "MARMOT_IFACE_GAUGE_W1" );
      return s ? std::atof( s ) : 1.0e-7;
    }();
    static const double hi = []() {
      const char* s = std::getenv( "MARMOT_IFACE_GAUGE_W0" );
      return s ? std::atof( s ) : 1.0e-5;
    }();
    if ( !( ratio > 0.0 ) || ratio <= lo )
      return 1.0;
    if ( ratio >= hi )
      return 0.0;
    const double x = ( std::log10( ratio ) - std::log10( lo ) ) / ( std::log10( hi ) - std::log10( lo ) );
    return 1.0 - x * x * ( 3.0 - 2.0 * x ); // smoothstep: C1 at both ends
  }

  // A global iteration is one sweep over every interface quadrature point. The sweep
  // boundary is detected by a point being visited twice; at that moment the counts for
  // the finished iteration are emitted and reset. This counts ALL points, not a sample.
  void reportIteration( unsigned long quadraturePointId, double weight, double timeOld )
  {
    static std::mutex                mutex;
    static std::set< unsigned long > seen;
    static double                    weightSum = 0.0;
    static long                      total = 0, iteration = 0;

    std::lock_guard< std::mutex > lock( mutex );
    if ( seen.count( quadraturePointId ) ) {
      ++iteration;
      std::printf( "[iter] %ld %.6e sumW=%.6f of %ld\n", iteration, timeOld, weightSum, total );
      std::fflush( stdout );
      seen.clear();
      weightSum = 0.0;
      total     = 0;
    }
    seen.insert( quadraturePointId );
    ++total;
    weightSum += weight;
  }

  // Guarded inverse of the averaged acoustic tensor Qbar.
  //
  // The model establishes that Qbar is invertible whenever the condensation is
  // non-trivial, and that the only configuration degenerating it is the
  // homogeneous one, f+ parallel f- and g+ parallel g-, in which the internal
  // variable it would solve for is itself zero. That statement is qualitative:
  // it bounds no singular value away from zero. Once both faces yield at small
  // hardening the smallest singular value of Qbar falls to ~H/2 against an
  // elastic scale K + 4G/3, and the plain inverse then amplifies the part of
  // the traction imbalance that no z can remove.
  //
  // The degeneracy is directional, so the treatment is too: solve in the
  // retained subspace, assign zero along the degenerate directions -- the value
  // the theory gives them -- and leave the irreducible residual to the global
  // Newton.
  struct GuardedAcousticInverse {
    Eigen::JacobiSVD< Eigen::Matrix3d > svd;
    Vector3d                            invSingular       = Vector3d::Zero();
    Matrix3dRowMajor                    retainedProjector = Matrix3dRowMajor::Zero();
    int                                 nRetained         = 0;
    Vector3d                            nullDirection     = Vector3d::Zero();
    double                              sMax              = 0.0;
    double                              sMinRetained      = 0.0;

    GuardedAcousticInverse( const Matrix3dRowMajor& Q, double relativeTolerance )
      : svd( Eigen::Matrix3d( Q ), Eigen::ComputeFullU | Eigen::ComputeFullV )
    {
      const Vector3d singular = svd.singularValues();
      sMax                    = singular( 0 );
      nullDirection           = svd.matrixV().col( 2 ); // right-singular vector of the SMALLEST value
      const double threshold  = relativeTolerance * sMax;
      for ( int i = 0; i < 3; ++i ) {
        if ( singular( i ) > threshold && singular( i ) > 0.0 ) {
          invSingular( i ) = 1.0 / singular( i );
          retainedProjector += svd.matrixU().col( i ) * svd.matrixU().col( i ).transpose();
          sMinRetained = singular( i );
          ++nRetained;
        }
      }
    }

    template < typename Derived >
    Eigen::Matrix< double, 3, Derived::ColsAtCompileTime > solve( const Eigen::MatrixBase< Derived >& rhs ) const
    {
      return svd.matrixV() * ( invSingular.asDiagonal() * ( svd.matrixU().transpose() * rhs ) );
    }

    // Part of a traction imbalance that some z can still remove.
    double reducibleNorm( const Vector3d& r, double scale ) const
    {
      if ( qGuardFullResidualConvergence() )
        return r.norm() / scale;
      return ( retainedProjector * r ).norm() / scale;
    }

    // with a gauge in force the solvable part is the projection on its basis
    static double gaugedNorm( const Eigen::Matrix< double, 3, 2 >& basis, const Vector3d& r, double scale )
    {
      return ( basis.transpose() * r ).norm() / scale;
    }
  };

  // Diagnostic: how often the guard truncates, and how far Qbar degenerates.
  void reportGuard( const GuardedAcousticInverse& guard )
  {
    static std::atomic< long > calls{ 0 };
    static std::atomic< long > truncatedCalls{ 0 };

    const long n = ++calls;
    if ( guard.nRetained < 3 )
      ++truncatedCalls;

    if ( n % 100000 == 0 ) {
      const double ratio = guard.sMax > 0.0 ? guard.sMinRetained / guard.sMax : 0.0;
      std::printf( "  [qguard] call=%ld truncated=%ld  relTol=%.1e  sMinRetained/sMax=%.3e  nRetained=%d\n",
                   n,
                   truncatedCalls.load(),
                   qGuardRelativeTolerance(),
                   ratio,
                   guard.nRetained );
      std::fflush( stdout );
    }
  }

  // Physics probe. Answers two questions with numbers instead of inference:
  //
  //   1. do the two faces actually carry the same elastoplastic tangent?
  //      Reported as kappa+ , kappa- and the tangent mismatch ||C+ - C-||/||C+||.
  //   2. which half of the traction continuity condition can the local solve
  //      actually enforce? The imbalance r = t+ - t- is split into the component
  //      along the normal, which the bulk response can still move, and the shear
  //      component, which the yield surface pins once both faces flow. Both are
  //      reported on entry to the local solve and on exit from it.
  //
  // Sampled one call in a thousand, printed every hundred samples.
  void reportPhysics( double        kappaPlus,
                      double        kappaMinus,
                      double        tangentMismatch,
                      double        rNormalIn,
                      double        rShearIn,
                      double        rNormalOut,
                      double        rShearOut,
                      double        cosTheta,
                      double        condRatio,
                      double        deviatorNormPlus,
                      double        deviatorNormMinus,
                      unsigned long quadraturePointId,
                      double        timeOld,
                      const double* zVector,
                      const double* nullDirection,
                      double        condensedTangentNorm,
                      bool          gaugeReport,
                      double        gaugeAlong,
                      double        gaugePerp,
                      double        stressNormPlus,
                      double        stressNormMinus,
                      double        tangentNormPlus,
                      double        tangentNormMinus,
                      const double* tractionPlus,
                      const double* tractionMinus,
                      int           localIterations,
                      int           lineSearchTrials,
                      double        dEpsTracePlus,
                      double        dEpsNormPlus,
                      double        pressurePlus,
                      double        pressureMinus )
  {
    static std::atomic< long > calls{ 0 };
    const long                 n = ++calls;

    // Full history for a sparse, FIXED subset of quadrature points, so that an
    // opposed flow direction can be told apart from one that flickers: physical
    // opposition persists at the same point across increments, noise does not.
    if ( gaugeReport && quadraturePointId % 256u == 0u ) {
      static std::mutex             gaugeMutex;
      std::lock_guard< std::mutex > gaugeLock( gaugeMutex );
      std::printf( "[gauge] %lu %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e "
                   "%.17e %.17e %.17e %.17e %.17e\n",
                   quadraturePointId,
                   timeOld,
                   gaugeAlong,
                   gaugePerp,
                   kappaPlus,
                   kappaMinus,
                   stressNormPlus,
                   stressNormMinus,
                   tangentNormPlus,
                   tangentNormMinus,
                   tractionPlus[0],
                   tractionPlus[1],
                   tractionPlus[2],
                   tractionMinus[0],
                   tractionMinus[1] );
      std::fflush( stdout );
    }

    if ( quadraturePointId % 16u == 0u ) {
      static std::mutex             grindMutex;
      std::lock_guard< std::mutex > grindLock( grindMutex );
      std::printf( "[grind] %lu %.6e %d %d %d %.9e %.9e %.9e %.9e %.9e %.9e\n",
                   quadraturePointId,
                   timeOld,
                   gaugeReport ? 1 : 0,
                   localIterations,
                   lineSearchTrials,
                   dEpsTracePlus,
                   dEpsNormPlus,
                   pressurePlus,
                   deviatorNormPlus,
                   pressureMinus,
                   deviatorNormMinus );
      std::fflush( stdout );
    }

    if ( quadraturePointId % 256u == 0u ) {
      static std::mutex             traceMutex;
      std::lock_guard< std::mutex > traceLock( traceMutex );
      std::printf( "[trace] %lu %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
                   quadraturePointId,
                   timeOld,
                   cosTheta,
                   deviatorNormPlus,
                   deviatorNormMinus,
                   kappaPlus,
                   kappaMinus,
                   gaugeReport ? 1.0 : 0.0,
                   zVector[1],
                   zVector[2],
                   nullDirection[0],
                   nullDirection[1],
                   nullDirection[2],
                   condensedTangentNorm );
      std::fflush( stdout );
    }

    if ( n % 200 != 0 )
      return;

    static std::mutex mutex;
    static long       samples    = 0;
    static double     sKappaPlus = 0.0, sKappaMinus = 0.0, sMismatch = 0.0, maxMismatch = 0.0;
    static double     sRnIn = 0.0, sRsIn = 0.0, sRnOut = 0.0, sRsOut = 0.0;
    static long       bothPlastic = 0;

    std::lock_guard< std::mutex > lock( mutex );
    ++samples;

    // one row per sample: orientation of plastic flow at the two faces against the
    // difference in their hardening state, so the two can be correlated afterwards.
    const double rIn  = std::sqrt( rNormalIn * rNormalIn + rShearIn * rShearIn );
    const double rOut = std::sqrt( rNormalOut * rNormalOut + rShearOut * rShearOut );
    std::printf( "[row] %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %lu %.6e\n",
                 kappaPlus,
                 kappaMinus,
                 kappaPlus - kappaMinus,
                 cosTheta,
                 tangentMismatch,
                 condRatio,
                 rIn,
                 rOut,
                 deviatorNormPlus,
                 deviatorNormMinus,
                 quadraturePointId,
                 timeOld );

    sKappaPlus += kappaPlus;
    sKappaMinus += kappaMinus;
    sMismatch += tangentMismatch;
    maxMismatch = std::max( maxMismatch, tangentMismatch );
    sRnIn += rNormalIn;
    sRsIn += rShearIn;
    sRnOut += rNormalOut;
    sRsOut += rShearOut;
    if ( kappaPlus > 0.0 && kappaMinus > 0.0 )
      ++bothPlastic;

    if ( samples % 100 == 0 ) {
      const double inv = 1.0 / static_cast< double >( samples );
      std::printf( "  [phys] n=%ld  kappa+=%.3e kappa-=%.3e  |dC|/|C|: mean=%.3e max=%.3e  "
                   "bothPlastic=%.0f%%  r_normal %.3e -> %.3e   r_shear %.3e -> %.3e\n",
                   n,
                   sKappaPlus * inv,
                   sKappaMinus * inv,
                   sMismatch * inv,
                   maxMismatch,
                   100.0 * static_cast< double >( bothPlastic ) * inv,
                   sRnIn * inv,
                   sRnOut * inv,
                   sRsIn * inv,
                   sRsOut * inv );
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
    double    scale;
    Vector6d  dEpsPlus  = Vector6d::Zero(); // input handed to the + face material
    Vector6d  dEpsMinus = Vector6d::Zero();
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
    evaluation.dEpsPlus  = strainIncrementPlusVoigt;
    evaluation.dEpsMinus = strainIncrementMinusVoigt;
    evaluation.plus      = evaluateSideTrial( *topMaterial,
                                         topStateVars,
                                         nTopStateVars,
                                         stressCurrentPlus,
                                         strainIncrementPlusVoigt,
                                         timeInfo );
    evaluation.minus     = evaluateSideTrial( *bottomMaterial,
                                          bottomStateVars,
                                          nBottomStateVars,
                                          stressCurrentMinus,
                                          strainIncrementMinusVoigt,
                                          timeInfo );

    const Vector3d tPlus  = evaluation.plus.stress * normal;
    const Vector3d tMinus = evaluation.minus.stress * normal;
    evaluation.r          = tPlus - tMinus;
    evaluation.scale      = std::max( 1.0, std::max( tPlus.norm(), tMinus.norm() ) );
    evaluation.rNorm      = evaluation.r.norm() / evaluation.scale;
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
  constexpr int maxLocalIterations = 40;
  // MARMOT_IFACE_LOCAL_TOL overrides it, so the Step 0 residual can be tested against the
  // local solver's own convergence rather than assumed to be a modelling error.
  static const double localTolerance = []() {
    const char* e = std::getenv( "MARMOT_IFACE_LOCAL_TOL" );
    return e ? std::atof( e ) : 1.0e-10;
  }();

  Vector3d z = Eigen::Map< const Vector3d >( normalGradientJump );

  LocalEvaluation current = evaluateLocal( z );

  const LocalEvaluation entryEval = current;
  reportCall();

  // ---------------------------------------------------------------------------
  // Gauge. Once both faces flow, the layer stops resisting a motion of g along
  //
  //     a_i = 1/2 ( N+_ij + N-_ij ) n_j ,   N = dev(sigma)/||dev(sigma)||,
  //
  // and the local solve becomes bistable: it converges equally well to two states
  // that differ by a jump along that direction, alternating between them from one
  // global iteration to the next. The traction is the same on both branches -- it
  // is the coupling blocks of the condensation, which are built from the two faces
  // SEPARATELY, that then differ, and the global residual becomes discontinuous.
  //
  // The direction is taken from the material rather than from an SVD of Qbar: the
  // singular vector is worst determined exactly here, and its sign is arbitrary,
  // which would make the gauge itself bistable. Only the projector I - a (x) a is
  // used, so the sign of a is irrelevant.
  // ---------------------------------------------------------------------------
  double gaugeWeightCurrent         = 0.0;
  double gaugeAlongBeforeProjection = 0.0;
  double gaugePerpBeforeProjection  = 0.0;
  bool   gaugeActive                = gaugeEnabled() && entryEval.plus.stateChanged && entryEval.minus.stateChanged;

  // ... and only where the averaged acoustic tensor has actually gone degenerate.
  //
  // Patch B: the criterion is det S, built from the material state
  //
  //   S = diag( 1/(mu thbar+) , 1/(mu thbar-) ) - U^T Qe^-1 U ,   U = [ m+ , m- ]
  //
  // rather than sigma_min/sigma_max of an SVD of Qbar. The spectral ratio wanders in the
  // solution, which let one quadrature point of 3600 toggle between a 2-DOF and a 3-DOF
  // local problem and hold the global Newton in a limit cycle. det S is smooth in the
  // state. MARMOT_IFACE_CRITERION=dets selects it.
  if ( gaugeActive ) {
    if ( criterionDetS() ) {
      const double E2  = baseMaterialProperties[0];
      const double nu2 = baseMaterialProperties[1];
      const double mu2 = E2 / ( 2.0 * ( 1.0 + nu2 ) );
      const double K2  = E2 / ( 3.0 * ( 1.0 - 2.0 * nu2 ) );

      double thP2, thBP2, thM2, thBM2;
      probeAlgorithmicTheta( entryEval.plus.stress, entryEval.plus.CFull, mu2, thP2, thBP2 );
      probeAlgorithmicTheta( entryEval.minus.stress, entryEval.minus.CFull, mu2, thM2, thBM2 );
      const double thAvg2 = 0.5 * ( thP2 + thM2 );

      const Matrix3dRowMajor nn2 = normal * normal.transpose();
      const Matrix3dRowMajor Qei = ( Matrix3dRowMajor::Identity() - nn2 ) / std::max( mu2 * thAvg2, 1.0e-300 ) +
                                   nn2 / ( K2 + 4.0 * mu2 * thAvg2 / 3.0 );
      Eigen::Matrix< double, 3, 2 > Uc;
      Uc.col( 0 )              = unitDeviator( entryEval.plus.stress ) * normal;
      Uc.col( 1 )              = unitDeviator( entryEval.minus.stress ) * normal;
      const Eigen::Matrix2d Gc = Uc.transpose() * Qei * Uc;
      Eigen::Matrix2d       Sc = -Gc;
      Sc( 0, 0 ) += 1.0 / std::max( mu2 * thBP2, 1.0e-300 );
      Sc( 1, 1 ) += 1.0 / std::max( mu2 * thBM2, 1.0e-300 );
      // Dimensionless: |det S| / |S11 S22|. No external scale to get wrong, and it is
      // the natural measure of how close a 2x2 is to singular.
      const double detSc  = Sc.determinant();
      const double diagSc = Sc( 0, 0 ) * Sc( 1, 1 );
      const double ratioS = std::abs( detSc ) / std::max( std::abs( diagSc ), 1.0e-300 );
      gaugeActive         = ratioS < gaugeDegeneracyThreshold();

      if ( std::getenv( "MARMOT_IFACE_CRITERION_LOG" ) && ( quadraturePointIdOf( state.stateVars ) % 256u ) == 0u ) {
        const Matrix3dRowMajor QEntry2 = 0.5 * ( Rt * entryEval.plus.CFull * Bn + Rt * entryEval.minus.CFull * Bn );
        const GuardedAcousticInverse  sp2( QEntry2, 0.0 );
        const double                  spectral = sp2.sMax > 0.0 ? sp2.sMinRetained / sp2.sMax : 0.0;
        static std::mutex             mc;
        std::lock_guard< std::mutex > lc( mc );
        std::printf( "[crit] thbP=%.6e thbM=%.6e detS=%.6e S11=%.6e S22=%.6e "
                     "scaleOld=%.6e ratioDimless=%.6e ratioOld=%.6e spectral=%.6e\n",
                     thBP2,
                     thBM2,
                     detSc,
                     Sc( 0, 0 ),
                     Sc( 1, 1 ),
                     1.0 / std::max( mu2 * mu2 * thBP2 * thBM2, 1.0e-300 ),
                     ratioS,
                     std::abs( detSc ) * mu2 * mu2 * thBP2 * thBM2,
                     spectral );
        std::fflush( stdout );
      }
    }
    else {
      const Matrix3dRowMajor       QEntry = 0.5 * ( Rt * entryEval.plus.CFull * Bn + Rt * entryEval.minus.CFull * Bn );
      const GuardedAcousticInverse entrySpectrum( QEntry, 0.0 );
      const double ratio = entrySpectrum.sMax > 0.0 ? entrySpectrum.sMinRetained / entrySpectrum.sMax : 0.0;
      gaugeActive        = ratio < gaugeDegeneracyThreshold();
    }
  }

  // MARMOT_IFACE_GAUGE_FREEZE=1 evaluates the criterion once per increment and holds it
  // through the global iterations. If the grind survives a frozen set, the toggling point
  // is a symptom rather than the cause.
  if ( gaugeFreeze() )
    gaugeActive = frozenGaugeDecision( quadraturePointIdOf( state.stateVars ), timeIncrement.timeOld, gaugeActive );
  Vector3d                      gaugeDirection = Vector3d::Zero();
  Eigen::Matrix< double, 3, 2 > gaugeBasis     = Eigen::Matrix< double, 3, 2 >::Zero();

  if ( gaugeActive ) {
    const Matrix3dRowMajor flowPlus  = unitDeviator( entryEval.plus.stress );
    const Matrix3dRowMajor flowMinus = unitDeviator( entryEval.minus.stress );

    const double flowSumNorm = Matrix3dRowMajor( flowPlus + flowMinus ).norm();

    Vector3d direction = Matrix3dRowMajor( 0.5 * ( flowPlus + flowMinus ) ) * normal;
    // exactly antiparallel faces: the average vanishes and either face spans the
    // same line, so fall back to one of them
    bool fallbackUsed = false;
    if ( direction.norm() <= 1.0e-10 ) {
      direction    = flowPlus * normal;
      fallbackUsed = true;
    }

    // Constrain only the TANGENTIAL part. With m^T . n = 0 the perturbation
    // delta_eps_kk = delta (m^T . n) vanishes identically, so the gauge never touches
    // the volumetric response -- including at the minority of points where m itself
    // carries a normal component (p95 1.8e-4 -> 7.7e-2, max 1.9e-1).
    if ( gaugeTangential() )
      direction -= direction.dot( normal ) * normal;

    reportGaugeSetup( flowSumNorm,
                      fallbackUsed,
                      direction.norm() > 1.0e-12 ? direction.normalized().dot( normal ) : 0.0 );

    if ( direction.norm() > 1.0e-12 ) {
      gaugeDirection       = direction.normalized();
      const Vector3d trial = ( std::abs( gaugeDirection[0] ) < 0.9 ) ? Vector3d::UnitX() : Vector3d::UnitY();
      const Vector3d first = ( trial - trial.dot( gaugeDirection ) * gaugeDirection ).normalized();
      gaugeBasis.col( 0 )  = first;
      gaugeBasis.col( 1 )  = gaugeDirection.cross( first ).normalized();
    }
    else {
      gaugeActive = false;
    }
  }

  if ( gaugeActive ) {
    // the discarded component, and the retained magnitude, as they stood before the
    // projection: evidence that what is thrown away was genuinely free
    gaugeAlongBeforeProjection = gaugeDirection.dot( z );
    gaugePerpBeforeProjection  = ( z - gaugeAlongBeforeProjection * gaugeDirection ).norm();

    z -= ( gaugeDirection.dot( z ) - gaugeOffset() ) * gaugeDirection;
    current = evaluateLocal( z );

    // Instantaneous test that the discarded component is genuinely free: evaluate the
    // same state displaced along m and compare stress, tangent and hardening. The
    // solution is not affected -- these are throwaway evaluations.
    const double delta = gaugeSelfCheckDelta();
    if ( delta > 0.0 && ( quadraturePointIdOf( state.stateVars ) % 256u ) == 0u ) {
      const LocalEvaluation up       = evaluateLocal( z + delta * gaugeDirection );
      const LocalEvaluation down     = evaluateLocal( z - delta * gaugeDirection );
      auto                  relative = []( double a, double b ) {
        const double scale = std::max( std::abs( a ), std::abs( b ) );
        return scale > 0.0 ? std::abs( a - b ) / scale : 0.0;
      };
      static std::mutex             checkMutex;
      std::lock_guard< std::mutex > checkLock( checkMutex );
      auto                          pressureOf = []( const Matrix3dRowMajor& sigma ) { return sigma.trace() / 3.0; };
      auto                          deviatorNormOf = []( const Matrix3dRowMajor& sigma ) {
        Matrix3dRowMajor d = sigma;
        const double     m = sigma.trace() / 3.0;
        d( 0, 0 ) -= m;
        d( 1, 1 ) -= m;
        d( 2, 2 ) -= m;
        return d.norm();
      };
      // check 2: how much of the flow direction, of the material gauge direction, and of
      // the SVD null direction points along the interface normal
      const Matrix3dRowMajor       NPlus    = unitDeviator( current.plus.stress );
      const Matrix3dRowMajor       NMinus   = unitDeviator( current.minus.stress );
      const double                 n33Plus  = normal.dot( NPlus * normal );
      const double                 n33Minus = normal.dot( NMinus * normal );
      const double                 mDotN    = gaugeDirection.dot( normal );
      const GuardedAcousticInverse svdGuard( Matrix3dRowMajor(
                                               0.5 * ( Rt * current.plus.CFull * Bn + Rt * current.minus.CFull * Bn ) ),
                                             0.0 );
      const double                 aDotN = svdGuard.nullDirection.dot( normal );
      std::printf( "[check2] N33p=%.6e N33m=%.6e  m.n=%.6e  a.n=%.6e  sMinRatio=%.6e  "
                   "p=%.6e  dp_pred=%.6e  delta=%.3e\n",
                   n33Plus,
                   n33Minus,
                   mDotN,
                   aDotN,
                   svdGuard.sMax > 0.0 ? svdGuard.sMinRetained / svdGuard.sMax : 0.0,
                   current.plus.stress.trace() / 3.0,
                   -( materialProperties[0] / ( 3.0 * ( 1.0 - 2.0 * materialProperties[1] ) ) ) * 0.5 * delta * mDotN,
                   delta );
      const Vector3d tCur = current.plus.stress * normal;
      const Vector3d tUp  = up.plus.stress * normal;
      const Vector3d tDn  = down.plus.stress * normal;
      std::printf( "[decomp] dev %.3e %.3e   p %.3e %.3e   |t| %.3e %.3e   dt/|t| %.3e %.3e\n",
                   relative( deviatorNormOf( up.plus.stress ), deviatorNormOf( current.plus.stress ) ),
                   relative( deviatorNormOf( down.plus.stress ), deviatorNormOf( current.plus.stress ) ),
                   relative( pressureOf( up.plus.stress ), pressureOf( current.plus.stress ) ),
                   relative( pressureOf( down.plus.stress ), pressureOf( current.plus.stress ) ),
                   relative( tUp.norm(), tCur.norm() ),
                   relative( tDn.norm(), tCur.norm() ),
                   ( tUp - tCur ).norm() / std::max( tCur.norm(), 1e-30 ),
                   ( tDn - tCur ).norm() / std::max( tCur.norm(), 1e-30 ) );
      std::printf( "[selfcheck] %.6e %.6e   dSigmaP=%.3e %.3e  dSigmaM=%.3e %.3e  "
                   "dCP=%.3e %.3e  dCM=%.3e %.3e  dKapP=%.6e %.6e  dKapM=%.6e %.6e  delta=%.3e\n",
                   timeIncrement.timeOld,
                   ( current.plus.stress * normal ).norm(),
                   relative( up.plus.stress.norm(), current.plus.stress.norm() ),
                   relative( down.plus.stress.norm(), current.plus.stress.norm() ),
                   relative( up.minus.stress.norm(), current.minus.stress.norm() ),
                   relative( down.minus.stress.norm(), current.minus.stress.norm() ),
                   relative( up.plus.CFull.norm(), current.plus.CFull.norm() ),
                   relative( down.plus.CFull.norm(), current.plus.CFull.norm() ),
                   relative( up.minus.CFull.norm(), current.minus.CFull.norm() ),
                   relative( down.minus.CFull.norm(), current.minus.CFull.norm() ),
                   up.plus.trialStateVars.empty() ? 0.0 : up.plus.trialStateVars[0] - current.plus.trialStateVars[0],
                   down.plus.trialStateVars.empty() ? 0.0
                                                    : down.plus.trialStateVars[0] - current.plus.trialStateVars[0],
                   up.minus.trialStateVars.empty() ? 0.0 : up.minus.trialStateVars[0] - current.minus.trialStateVars[0],
                   down.minus.trialStateVars.empty() ? 0.0
                                                     : down.minus.trialStateVars[0] - current.minus.trialStateVars[0],
                   delta );
      std::fflush( stdout );
    }
  }

  // Convergence, the line search and the best-iterate bookkeeping are all measured
  // on the REDUCIBLE part of the imbalance, the component the retained directions
  // of Qbar can still act on. The rest cannot be removed by any z and is carried by
  // the global Newton, which is what the condensation asks of it.
  Vector3d        bestZ         = z;
  LocalEvaluation bestEval      = current;
  double          bestReducible = std::numeric_limits< double >::max();

  auto guardAt = [&]( const LocalEvaluation& evaluation ) {
    const Matrix3dRowMajor Qplus  = Rt * evaluation.plus.CFull * Bn;
    const Matrix3dRowMajor Qminus = Rt * evaluation.minus.CFull * Bn;
    const Matrix3dRowMajor Qavg   = 0.5 * ( Qplus + Qminus );
    return GuardedAcousticInverse( Qavg, qGuardRelativeTolerance() );
  };

  auto recordIfBetter =
    [&]( const GuardedAcousticInverse& guard, const LocalEvaluation& evaluation, const Vector3d& zTrial ) {
      const double reducible = gaugeActive
                                 ? GuardedAcousticInverse::gaugedNorm( gaugeBasis, evaluation.r, evaluation.scale )
                                 : guard.reducibleNorm( evaluation.r, evaluation.scale );
      if ( reducible < bestReducible ) {
        bestReducible = reducible;
        bestEval      = evaluation;
        bestZ         = zTrial;
      }
      return reducible;
    };

  int localIterationsUsed  = 0;
  int lineSearchTrialsUsed = 0;
  for ( int iteration = 0; iteration < maxLocalIterations; ++iteration ) {
    ++localIterationsUsed;
    const GuardedAcousticInverse guard = guardAt( current );
    reportGuard( guard );

    // Fully degenerate: every direction has stopped responding, which is the
    // configuration in which the theory gives z = 0 outright.
    if ( guard.nRetained == 0 )
      break;

    const double currentReducible = recordIfBetter( guard, current, z );
    if ( currentReducible <= localTolerance )
      break;

    Vector3d dz;
    if ( gaugeActive ) {
      const Matrix3dRowMajor QavgLocal   = 0.5 * ( Rt * current.plus.CFull * Bn + Rt * current.minus.CFull * Bn );
      const Eigen::Matrix2d  restrictedQ = gaugeBasis.transpose() * QavgLocal * gaugeBasis;
      const Eigen::Vector2d  restrictedR = -( gaugeBasis.transpose() * current.r );
      dz                                 = gaugeBasis * restrictedQ.fullPivLu().solve( restrictedR );
    }
    else {
      dz = guard.solve( -current.r );
    }

    // Backtracking: a full step can overshoot when a trial state crosses a yield-surface
    // kink, so only a step that reduces the reducible residual is accepted.
    double          alpha = 1.0;
    LocalEvaluation candidate;
    bool            accepted = false;
    for ( int lineSearchIter = 0; lineSearchIter < 12; ++lineSearchIter ) {
      candidate                     = evaluateLocal( z + alpha * dz );
      const double candidateReduced = gaugeActive
                                        ? GuardedAcousticInverse::gaugedNorm( gaugeBasis, candidate.r, candidate.scale )
                                        : guard.reducibleNorm( candidate.r, candidate.scale );
      if ( candidateReduced <= localTolerance || candidateReduced < currentReducible ) {
        accepted = true;
        break;
      }
      alpha *= 0.5;
      ++lineSearchTrialsUsed;
    }
    if ( !accepted )
      break;

    z += alpha * dz;
    current = candidate;
  }

  // The iterate the loop ends on has not been ranked yet.
  recordIfBetter( guardAt( current ), current, z );

  z       = bestZ;
  current = bestEval;

  double        probeKappaPlus = 0.0, probeKappaMinus = 0.0, probeMismatch = 0.0;
  double        probeRnIn = 0.0, probeRsIn = 0.0, probeRnOut = 0.0, probeRsOut = 0.0;
  double        probeCosTheta = 0.0, probeCondRatio = 0.0, probeDevPlus = 0.0, probeDevMinus = 0.0;
  unsigned long probeQpId        = 0;
  Vector3d      probeNull        = Vector3d::Zero();
  Vector3d      probeZ           = Vector3d::Zero();
  bool          probeGaugeActive = false;
  double        probeGaugeAlong = 0.0, probeGaugePerp = 0.0;
  double        probeStressPlus = 0.0, probeStressMinus = 0.0;
  double        probeTangentPlus = 0.0, probeTangentMinus = 0.0;
  Vector3d      probeTractionPlus = Vector3d::Zero(), probeTractionMinus = Vector3d::Zero();
  double        probeDEpsTrace = 0.0, probeDEpsNorm = 0.0;
  double        probePressurePlus = 0.0, probePressureMinus = 0.0;

  // ---------------------------------------------------------------------------
  // Step 0 of the exact-condensation spec. Diagnostic only.
  //
  // Rebuild Qbar^ep as the rank-two downdate of Q^e and Phi_i from
  //
  //   Phi_i = lam (tr D) n_i + mu d_i - alpha+ A+_i + alpha- A-_i ,   d_i = n_k D_ki
  //
  // and test (13),  Qbar^ep g + Phi = 0,  at the converged g of the unmodified solve.
  // It is an identity, so it must sit at solver tolerance everywhere. The mu d_i term is
  // the one every earlier attempt dropped; its magnitude is logged separately.
  //
  // <Du_{k,l}> is ambiguous between the average SURFACE gradient and the full average
  // gradient including gBar (x) n, so both readings are evaluated and reported.
  if ( std::getenv( "MARMOT_IFACE_STEP0" ) ) {
    const double E   = baseMaterialProperties[0];
    const double nu  = baseMaterialProperties[1];
    const double Hh  = baseMaterialProperties[3];
    const double mu  = E / ( 2.0 * ( 1.0 + nu ) );
    const double lam = E * nu / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu ) );

    // associated J2: a = b = 2 mu M , A = B = 2 mu M n , Dpm = H + 2 mu
    const Matrix3dRowMajor Mp = unitDeviator( current.plus.stress );
    const Matrix3dRowMajor Mm = unitDeviator( current.minus.stress );
    const Vector3d         Ap = 2.0 * mu * ( Mp * normal );
    const Vector3d         Am = 2.0 * mu * ( Mm * normal );
    const double           Dp = Hh + 2.0 * mu;
    const double           Dm = Hh + 2.0 * mu;

    const Matrix3dRowMajor nn    = normal * normal.transpose();
    const Matrix3dRowMajor Qe    = mu * Matrix3dRowMajor::Identity() + ( lam + mu ) * nn;
    const Matrix3dRowMajor Qeinv = ( Matrix3dRowMajor::Identity() - nn ) / mu + nn / ( lam + 2.0 * mu );

    Eigen::Matrix< double, 3, 2 > Umat, Vmat;
    Umat.col( 0 ) = Ap;
    Umat.col( 1 ) = Am;
    Vmat.col( 0 ) = Ap;
    Vmat.col( 1 ) = Am; // associated: B = A

    const Matrix3dRowMajor QbarSpec = Qe - 0.5 * ( Ap * Ap.transpose() / Dp + Am * Am.transpose() / Dm );

    const Eigen::Matrix2d Gmat = Vmat.transpose() * Qeinv * Umat;
    Eigen::Matrix2d       Smat;
    Smat << 2.0 * Dp - Gmat( 0, 0 ), -Gmat( 0, 1 ), -Gmat( 1, 0 ), 2.0 * Dm - Gmat( 1, 1 );
    const double detS = Smat.determinant();

    // ---- 2a / 2b: is the mismatch geometric, or is CFull the ALGORITHMIC tangent? ----
    // theta and thetabar are read straight out of CFull:
    //   C:X = 2 mu theta X          for X deviatoric, symmetric, orthogonal to M
    //   M:C:M = 2 mu ( theta - thetabar )
    auto probeTheta = [&]( const Matrix9dRowMajor& C, const Matrix3dRowMajor& M, double& theta, double& thetaBar ) {
      Matrix3dRowMajor X;
      X << 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0; // symmetric, traceless
      X -= ( X.array() * M.array() ).sum() * M;          // remove the M component
      const double                       xn = X.norm();
      const Eigen::Map< const Vector9d > vX( X.data() );
      const Eigen::Map< const Vector9d > vM( M.data() );
      theta    = xn > 1.0e-12 ? double( ( vX.transpose() * C * vX )( 0, 0 ) ) / ( 2.0 * mu * xn * xn ) : 1.0;
      thetaBar = theta - double( ( vM.transpose() * C * vM )( 0, 0 ) ) / ( 2.0 * mu );
    };
    double thetaP, thetaBarP, thetaM, thetaBarM;
    probeTheta( current.plus.CFull, Mp, thetaP, thetaBarP );
    probeTheta( current.minus.CFull, Mm, thetaM, thetaBarM );

    // ---- section 3: theta two ways, and the size of the strain increment ----
    // theta measured from the stresses themselves: s^tr = dev(sigma_n) + 2 mu dev(deps)
    auto devOf = []( const Matrix3dRowMajor& t ) {
      Matrix3dRowMajor d = t;
      const double     m = t.trace() / 3.0;
      d( 0, 0 ) -= m;
      d( 1, 1 ) -= m;
      d( 2, 2 ) -= m;
      return d;
    };
    const Vector3d         gP3       = gBar + 0.5 * z;
    const Vector3d         gM3       = gBar - 0.5 * z;
    const Matrix3dRowMajor GP3       = APlus + gP3 * normal.transpose();
    const Matrix3dRowMajor GM3       = AMinus + gM3 * normal.transpose();
    const Matrix3dRowMajor epsP      = 0.5 * ( GP3 + GP3.transpose() );
    const Matrix3dRowMajor epsM      = 0.5 * ( GM3 + GM3.transpose() );
    const Matrix3dRowMajor strP      = devOf( stressCurrentPlus ) + 2.0 * mu * devOf( epsP );
    const Matrix3dRowMajor strM      = devOf( stressCurrentMinus ) + 2.0 * mu * devOf( epsM );
    const double           thetaDirP = strP.norm() > 1.0e-300 ? devOf( current.plus.stress ).norm() / strP.norm() : 1.0;
    const double thetaDirM   = strM.norm() > 1.0e-300 ? devOf( current.minus.stress ).norm() / strM.norm() : 1.0;
    const double yieldStrain = baseMaterialProperties[2] / E;
    const double depsRatio   = epsP.norm() / std::max( yieldStrain, 1.0e-300 );

    const double   Kbulk = E / ( 3.0 * ( 1.0 - 2.0 * nu ) );
    const Vector3d mpv = Mp * normal, mmv = Mm * normal;
    // Corrected: Qbar is the AVERAGE of the two per-face acoustic tensors
    //   Q^alg+- = (K + mu th/3) n n + mu th I - 2 mu thbar m (x) m ,
    // so averaging halves the downdate. All the strength lives in mu thetabar; there is
    // no separate A/B/D triple for associated J2, and routing through one double-counted.
    const double           thAvg   = 0.5 * ( thetaP + thetaM );
    const Matrix3dRowMajor QbarAlg = ( Kbulk + mu * thAvg / 3.0 ) * nn + mu * thAvg * Matrix3dRowMajor::Identity() -
                                     mu * thetaBarP * ( mpv * mpv.transpose() ) -
                                     mu * thetaBarM * ( mmv * mmv.transpose() );

    // corrected criterion: S = Ccal^-1 - G , G_ab = m^a . Qe^-1 . m^b
    const Matrix3dRowMajor QeinvAlg = ( Matrix3dRowMajor::Identity() - nn ) / ( mu * thAvg ) +
                                      nn / ( Kbulk + 4.0 * mu * thAvg / 3.0 );
    Eigen::Matrix< double, 3, 2 > Umb;
    Umb.col( 0 )              = mpv;
    Umb.col( 1 )              = mmv;
    const Eigen::Matrix2d Gmb = Umb.transpose() * QeinvAlg * Umb;
    Eigen::Matrix2d       Smb = -Gmb;
    Smb( 0, 0 ) += 1.0 / std::max( mu * thetaBarP, 1.0e-300 );
    Smb( 1, 1 ) += 1.0 / std::max( mu * thetaBarM, 1.0e-300 );
    const double detSmb = Smb.determinant();

    const Matrix3dRowMajor QbarCode = 0.5 * ( Rt * current.plus.CFull * Bn + Rt * current.minus.CFull * Bn );

    // Phi, built directly from the face contractions.
    //
    //   F^{e,alg}(th, X)_i = ( K - 2/3 mu th ) tr(X) n_i + mu th ( X_ij n_j + X_ji n_j )
    //   F^{p,alg}+-(X)_i   = 2 mu thbar+- m+-_i ( M+- : X )
    //   F+- = F^{e,alg}(th+-) - F^{p,alg}+-
    //   Phi = [F]( <Du> ) + <F>( D )
    //
    // Two coefficients that earlier versions got wrong: the elastic prefactor is
    // K - 2/3 mu <theta>, NOT lam = K - 2/3 mu, and the d_i term carries mu <theta>,
    // not mu. Both matter because theta ~ 0.018.
    const Matrix3dRowMajor Dj    = APlus - AMinus;                   // D_kl
    const Vector3d         dvec  = Dj.transpose() * normal;          // d_i = n_k D_ki
    const double           trD   = Dj.trace();
    const Matrix3dRowMajor DuAvg = ABar + gBar * normal.transpose(); // FULL average gradient

    // MARMOT_IFACE_PHI_LAMBDA=elastic restores K - 2/3 mu (the control); default is the
    // theta-corrected K - 2/3 mu <theta>.
    static const bool lambdaElastic = []() {
      const char* e = std::getenv( "MARMOT_IFACE_PHI_LAMBDA" );
      return e && std::string( e ) == "elastic";
    }();
    auto Fe = [&]( double th, const Matrix3dRowMajor& X ) {
      const double coeff = lambdaElastic ? ( Kbulk - 2.0 / 3.0 * mu ) : ( Kbulk - 2.0 / 3.0 * mu * th );
      return Vector3d( coeff * X.trace() * normal + mu * th * ( X * normal + X.transpose() * normal ) );
    };
    auto Fp = [&]( double thBar, const Matrix3dRowMajor& M, const Vector3d& mv, const Matrix3dRowMajor& X ) {
      return Vector3d( 2.0 * mu * thBar * ( M.array() * X.array() ).sum() * mv );
    };

    const Vector3d FeJumpDu = Fe( thetaP, DuAvg ) - Fe( thetaM, DuAvg );
    const Vector3d FpJumpDu = Fp( thetaBarP, Mp, mpv, DuAvg ) - Fp( thetaBarM, Mm, mmv, DuAvg );
    const Vector3d FeAvgD   = Fe( thAvg, Dj ); // linear in theta
    const Vector3d FpAvgD   = 0.5 * ( Fp( thetaBarP, Mp, mpv, Dj ) + Fp( thetaBarM, Mm, mmv, Dj ) );

    const Vector3d PhiFull = FeJumpDu - FpJumpDu + FeAvgD - FpAvgD;
    const Vector3d PhiSurf = PhiFull; // the surface-only reading is retired

    // The residual must be formed with the tensor the element actually assembles.
    // QbarSpec is the CONTINUUM tensor, disproved earlier (60-75 % out in the transverse
    // block); using it here made every reported Step 0 residual meaningless.
    const Vector3d resSurf = QbarSpec * z + PhiFull; // kept only as the old, wrong measure
    const Vector3d resFull = QbarCode * z + PhiFull; // the identity, correctly formed
    const double   scal    = std::max( { ( QbarCode * z ).norm(), PhiFull.norm(), 1.0e-300 } );

    // ---- Woodbury solve, run alongside the existing local Newton for comparison ----
    //   phi  = Qe^-1 Phi ,  beta = U^T phi ,  g = -phi - Qe^-1 U S^-1 beta
    // Diagnostic only: g from this closed form is compared against z from the local
    // Newton. The old path stays until they agree.
    const Vector3d        phiW    = QeinvAlg * PhiFull;
    const Eigen::Vector2d betaW   = Umb.transpose() * phiW;
    const Vector3d        gWood   = -phiW - QeinvAlg * Umb * Smb.fullPivLu().solve( betaW );
    const double          woodRel = ( gWood - z ).norm() / std::max( z.norm(), 1.0e-300 );
    const double woodDot = ( gWood.norm() > 0 && z.norm() > 0 ) ? gWood.dot( z ) / ( gWood.norm() * z.norm() ) : 0.0;

    // ---- 1a: how much of the traction imbalance is IRREDUCIBLE ----
    // bhat is the left null direction of Qbar: no g_k removes a residual component along
    // it. Fredholm says the local system is solvable iff bhat.Phi = 0, which holds for
    // identical faces in layer-plane shear and fails for dissimilar ones. What is handed
    // to the global Newton has never been measured.
    Eigen::JacobiSVD< Eigen::Matrix3d > svdQ( Eigen::Matrix3d( QbarCode ), Eigen::ComputeFullU | Eigen::ComputeFullV );
    const Vector3d                      bHat    = svdQ.matrixU().col( 2 );
    const double                        rAlongB = std::abs( bHat.dot( current.r ) );
    const double                        rPerpB  = ( current.r - bHat.dot( current.r ) * bHat ).norm();
    const double phiAlongB = std::abs( bHat.dot( PhiFull ) ) / std::max( PhiFull.norm(), 1.0e-300 );
    const double tScale    = std::max(
      { ( current.plus.stress * normal ).norm(), ( current.minus.stress * normal ).norm(), 1.0e-300 } );

    // Section 2 bracket: factor K out of the n-component before subtracting, so the
    // near-cancellation is resolved at full precision rather than costing seven digits.
    const double gDotN   = z.dot( normal );
    const double bracket = ( 1.0 + mu * thAvg / ( 3.0 * Kbulk ) ) * gDotN +
                           ( 1.0 - 2.0 * mu * thAvg / ( 3.0 * Kbulk ) ) * trD;
    const double bracketRel = std::abs( bracket ) / std::max( { std::abs( gDotN ), std::abs( trD ), 1.0e-300 } );

    const double qMismatch = ( QbarSpec - QbarCode ).norm() / std::max( QbarCode.norm(), 1.0e-300 );

    // ---- Run 2: which direction does the residual live in? ----
    const Vector3d resid = QbarCode * z + PhiFull;

    // ---- sec 1: is <theta> actually alive at the point Phi is assembled? ----
    const double coeffTheta   = Kbulk - 2.0 / 3.0 * mu * thAvg;
    const double coeffElastic = Kbulk - 2.0 / 3.0 * mu;
    const double coeffRatio   = coeffTheta / coeffElastic;
    const bool   thetaSane    = ( thAvg > 0.0 && thAvg <= 1.0 );

    // ---- sec 2: the three contractions of D. D_il n_l must vanish. ----
    const Vector3d Dn    = Dj * normal;             // D_il n_l  -> must be ZERO
    const Vector3d nD    = Dj.transpose() * normal; // n_k D_ki = d_i
    const double   DnRel = Dn.norm() / std::max( Dj.norm(), 1.0e-300 );

    // ---- sec 3: residual per TERM, not as a fraction of ||Phi|| ----
    const double termN = std::abs( coeffTheta * trD );
    const double termD = ( mu * thAvg * dvec ).norm();
    const double termP = FpAvgD.norm() + FpJumpDu.norm();
    const double relN  = std::abs( resid.dot( normal ) ) / std::max( termN, 1.0e-300 );
    const double relD  = ( dvec.norm() > 1e-300 ? std::abs( resid.dot( dvec ) ) / dvec.norm() : 0.0 ) /
                        std::max( termD, 1.0e-300 );
    const double relP = ( mpv.norm() > 1e-300 ? std::abs( resid.dot( mpv ) ) / mpv.norm() : 0.0 ) /
                        std::max( termP, 1.0e-300 );
    auto share = [&]( const Vector3d& dir ) {
      const double dn = dir.norm();
      return dn > 1.0e-300 ? std::abs( resid.dot( dir ) ) / dn / std::max( resid.norm(), 1.0e-300 ) : 0.0;
    };
    const double shN    = share( normal );
    const double shMp   = share( mpv );
    const double shMm   = share( mmv );
    const double shD    = share( dvec );
    const double trDrel = std::abs( ( Kbulk - 2.0 / 3.0 * mu * thAvg ) * trD ) / std::max( PhiFull.norm(), 1.0e-300 );
    // ---- 2a(3): thetabar read from CFull vs reconstructed from theta and H ----
    const double thetaBarReconP = 1.0 / ( 1.0 + Hh / ( 3.0 * mu ) ) - ( 1.0 - thetaP );
    const double thetaBarReconM = 1.0 / ( 1.0 + Hh / ( 3.0 * mu ) ) - ( 1.0 - thetaM );
    const double tbRelP         = std::abs( thetaBarReconP - thetaBarP ) / std::max( std::abs( thetaBarP ), 1.0e-300 );
    // the quantity that actually sets the transverse eigenvalue
    const double diffMeasP  = thetaP - thetaBarP;
    const double diffReconP = thetaP - thetaBarReconP; // = H/(3 mu + H) exactly

    // ---- 2b, redone in a proper ORTHONORMAL tangential basis. Projecting a 3x3 with
    // (I - n n^T) leaves a structural zero along n, which the previous extraction picked
    // up instead of the smallest transverse eigenvalue. ----
    Vector3d t1 = ( std::abs( normal[0] ) < 0.9 ) ? Vector3d::UnitX() : Vector3d::UnitY();
    t1 -= t1.dot( normal ) * normal;
    t1.normalize();
    const Vector3d                t2 = normal.cross( t1 );
    Eigen::Matrix< double, 3, 2 > Tb;
    Tb.col( 0 ) = t1;
    Tb.col( 1 ) = t2;

    const Eigen::Matrix2d QttCode     = Tb.transpose() * QbarCode * Tb;
    const Eigen::Matrix2d QttAlg      = Tb.transpose() * QbarAlg * Tb;
    const double          lamPerpMeas = QttCode.eigenvalues().real().minCoeff();
    const double          lamPerpAlg  = QttAlg.eigenvalues().real().minCoeff();
    const double          lamPerpPred = mu * Hh / ( 3.0 * mu + Hh );
    const double          dQtt        = ( QttAlg - QttCode ).norm() / std::max( QttCode.norm(), 1.0e-300 );

    const Eigen::Vector2d QntCode = Tb.transpose() * ( QbarCode * normal );
    const Eigen::Vector2d QntAlg  = Tb.transpose() * ( QbarAlg * normal );
    const double          dQnt    = QntCode.norm();
    const double          Qnn     = std::abs( normal.dot( QbarCode * normal ) );
    const double          dQnn    = std::abs( normal.dot( ( QbarAlg - QbarCode ) * normal ) );
    const double          Qtt     = QntAlg.norm();
    // hardening-driven vs misalignment-driven degeneracy
    const Vector3d mhatP     = mpv.norm() > 1e-300 ? Vector3d( mpv / mpv.norm() ) : Vector3d::Zero();
    const Vector3d mhatM     = mmv.norm() > 1e-300 ? Vector3d( mmv / mmv.norm() ) : Vector3d::Zero();
    const double   cosPhiVec = mhatP.dot( mhatM );
    const double   lamH      = mu * Hh / ( 3.0 * mu + Hh );
    const double   lamPhi    = 0.5 * mu * ( 1.0 - std::abs( cosPhiVec ) );
    const double   sigMaxAlg = Kbulk + 4.0 * mu * thAvg / 3.0;

    const double n33P = normal.dot( Mp * normal );
    const double n33M = normal.dot( Mm * normal );

    const double algMismatch = ( QbarAlg - QbarCode ).norm() / std::max( QbarCode.norm(), 1.0e-300 );
    const double normRatio   = QbarCode.norm() / std::max( QbarSpec.norm(), 1.0e-300 );

    if ( current.plus.stateChanged && current.minus.stateChanged &&
         ( quadraturePointIdOf( state.stateVars ) % 64u ) == 0u ) {
      static std::mutex             m0;
      std::lock_guard< std::mutex > l0( m0 );
      std::printf( "[w0] resSurf=%.6e resFull=%.6e scale=%.6e qMis=%.6e detS=%.6e sc4D=%.6e "
                   "T_lam=%.6e T_mud=%.6e T_aP=%.6e T_aM=%.6e Qg=%.6e "
                   "thP=%.6e thBP=%.6e thM=%.6e thBM=%.6e algMis=%.6e nRat=%.6e "
                   "thDirP=%.6e thDirM=%.6e depsRat=%.6e tbRelP=%.6e "
                   "dMeas=%.9e dRecon=%.9e dQnn=%.6e dQtt=%.6e dQnt=%.6e Qnn=%.6e Qtt=%.6e "
                   "lamMeas=%.6e lamPred=%.6e lamAlg=%.6e n33P=%.6e n33M=%.6e QntAlg=%.6e "
                   "cosPhiV=%.9f lamH=%.6e lamPhi=%.6e sigMaxA=%.6e detSmb=%.6e shN=%.6e shMp=%.6e shMm=%.6e shD=%.6e "
                   "trDrel=%.6e trD=%.6e "
                   "rAlongB=%.6e rPerpB=%.6e phiAlongB=%.6e tScale=%.6e rNormAbs=%.6e woodRel=%.6e woodDot=%.6f "
                   "woodN=%.6e zN=%.6e gDotN=%.9e trD2=%.9e brRel=%.6e thSane=%d cRatio=%.6e DnRel=%.6e Dn=%.6e "
                   "nD=%.6e relN=%.6e relD=%.6e relP=%.6e\n",
                   resSurf.norm(),
                   resFull.norm(),
                   scal,
                   qMismatch,
                   detS,
                   4.0 * Dp * Dm,
                   FeAvgD.norm(),
                   FeJumpDu.norm(),
                   FpAvgD.norm(),
                   FpJumpDu.norm(),
                   ( QbarCode * z ).norm(),
                   thetaP,
                   thetaBarP,
                   thetaM,
                   thetaBarM,
                   algMismatch,
                   normRatio,
                   thetaDirP,
                   thetaDirM,
                   depsRatio,
                   tbRelP,
                   diffMeasP,
                   diffReconP,
                   dQnn,
                   dQtt,
                   dQnt,
                   Qnn,
                   Qtt,
                   lamPerpMeas,
                   lamPerpPred,
                   lamPerpAlg,
                   n33P,
                   n33M,
                   Qtt,
                   cosPhiVec,
                   lamH,
                   lamPhi,
                   sigMaxAlg,
                   detSmb,
                   shN,
                   shMp,
                   shMm,
                   shD,
                   trDrel,
                   trD,
                   rAlongB,
                   rPerpB,
                   phiAlongB,
                   tScale,
                   current.r.norm(),
                   woodRel,
                   woodDot,
                   gWood.norm(),
                   z.norm(),
                   gDotN,
                   trD,
                   bracketRel,
                   thetaSane ? 1 : 0,
                   coeffRatio,
                   DnRel,
                   Dn.norm(),
                   nD.norm(),
                   relN,
                   relD,
                   relP );
      std::fflush( stdout );
    }
  }

  {
    // normal / shear split of the traction imbalance, entering and leaving the solve
    auto split = [&]( const LocalEvaluation& e, double& normalPart, double& shearPart ) {
      const double alongNormal = e.r.dot( normal );
      normalPart               = std::abs( alongNormal ) / e.scale;
      shearPart                = ( e.r - alongNormal * normal ).norm() / e.scale;
    };
    double rnIn, rsIn, rnOut, rsOut;
    split( entryEval, rnIn, rsIn );
    split( current, rnOut, rsOut );

    const double kappaPlus  = nTopStateVars > 0 ? current.plus.trialStateVars[0] : 0.0;
    const double kappaMinus = nBottomStateVars > 0 ? current.minus.trialStateVars[0] : 0.0;

    const double normPlus = current.plus.CFull.norm();
    const double mismatch = normPlus > 0.0 ? ( current.plus.CFull - current.minus.CFull ).norm() / normPlus : 0.0;

    // Orientation of plastic flow at each face: the unit deviatoric stress direction,
    // which is the direction of the plastic flow for J2. cos(theta) = N+ : N- is 1 when
    // the two faces flow the same way.
    auto flowDirection = []( const Matrix3dRowMajor& sigma ) {
      Matrix3dRowMajor deviator = sigma;
      const double     mean     = sigma.trace() / 3.0;
      deviator( 0, 0 ) -= mean;
      deviator( 1, 1 ) -= mean;
      deviator( 2, 2 ) -= mean;
      const double norm = deviator.norm();
      return Matrix3dRowMajor( norm > 0.0 ? Matrix3dRowMajor( deviator / norm )
                                          : Matrix3dRowMajor( Matrix3dRowMajor::Zero() ) );
    };
    const Matrix3dRowMajor flowPlus  = flowDirection( current.plus.stress );
    const Matrix3dRowMajor flowMinus = flowDirection( current.minus.stress );
    const double           cosTheta  = ( flowPlus.array() * flowMinus.array() ).sum();

    const GuardedAcousticInverse finalGuard = guardAt( current );
    const double                 condRatio  = finalGuard.sMax > 0.0 ? finalGuard.sMinRetained / finalGuard.sMax : 0.0;

    auto deviatorNorm = []( const Matrix3dRowMajor& sigma ) {
      Matrix3dRowMajor deviator = sigma;
      const double     mean     = sigma.trace() / 3.0;
      deviator( 0, 0 ) -= mean;
      deviator( 1, 1 ) -= mean;
      deviator( 2, 2 ) -= mean;
      return deviator.norm();
    };

    // identity of this quadrature point: its state block is persistent across increments
    const unsigned long quadraturePointId = static_cast< unsigned long >(
      ( reinterpret_cast< std::uintptr_t >( state.stateVars ) >> 4 ) & 0xffffffffu );

    probeKappaPlus     = kappaPlus;
    probeKappaMinus    = kappaMinus;
    probeMismatch      = mismatch;
    probeRnIn          = rnIn;
    probeRsIn          = rsIn;
    probeRnOut         = rnOut;
    probeRsOut         = rsOut;
    probeCosTheta      = cosTheta;
    probeCondRatio     = condRatio;
    probeDevPlus       = deviatorNorm( current.plus.stress );
    probeDevMinus      = deviatorNorm( current.minus.stress );
    probeQpId          = quadraturePointId;
    probeNull          = finalGuard.nullDirection;
    probeZ             = z;
    probeGaugeActive   = gaugeActive;
    probeGaugeAlong    = gaugeAlongBeforeProjection;
    probeGaugePerp     = gaugeActive ? gaugePerpBeforeProjection : z.norm();
    probeStressPlus    = current.plus.stress.norm();
    probeStressMinus   = current.minus.stress.norm();
    probeTangentPlus   = current.plus.CFull.norm();
    probeTangentMinus  = current.minus.CFull.norm();
    probeTractionPlus  = current.plus.stress * normal;
    probeTractionMinus = current.minus.stress * normal;
    probeDEpsTrace     = current.dEpsPlus[0] + current.dEpsPlus[1] + current.dEpsPlus[2];
    probeDEpsNorm      = current.dEpsPlus.norm();
    probePressurePlus  = current.plus.stress.trace() / 3.0;
    probePressureMinus = current.minus.stress.trace() / 3.0;

    reportIteration( quadraturePointId, gaugeActive ? 1.0 : 0.0, timeIncrement.timeOld );
  }

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

  // The same guard, and the same truncation, as the local solve: the direction the
  // iteration refuses to solve in is the direction the tangent refuses to condense
  // through. Using a different rule here would leave Kcond no longer the derivative
  // of the residual the element returns.
  Matrix3x21RowMajor QavgInvRx;
  if ( gaugeActive ) {
    const Eigen::Matrix2d restrictedQ = gaugeBasis.transpose() * Qavg * gaugeBasis;
    QavgInvRx = gaugeBasis * restrictedQ.fullPivLu().solve( ( gaugeBasis.transpose() * rx ).eval() );
  }
  else {
    const GuardedAcousticInverse condensationGuard( Qavg, qGuardRelativeTolerance() );
    QavgInvRx = condensationGuard.solve( rx );
  }
  const Matrix21dRowMajor Kcond = px - pz * QavgInvRx;

  reportPhysics( probeKappaPlus,
                 probeKappaMinus,
                 probeMismatch,
                 probeRnIn,
                 probeRsIn,
                 probeRnOut,
                 probeRsOut,
                 probeCosTheta,
                 probeCondRatio,
                 probeDevPlus,
                 probeDevMinus,
                 probeQpId,
                 timeIncrement.timeOld,
                 probeZ.data(),
                 probeNull.data(),
                 Kcond.norm(),
                 probeGaugeActive,
                 probeGaugeAlong,
                 probeGaugePerp,
                 probeStressPlus,
                 probeStressMinus,
                 probeTangentPlus,
                 probeTangentMinus,
                 probeTractionPlus.data(),
                 probeTractionMinus.data(),
                 localIterationsUsed,
                 lineSearchTrialsUsed,
                 probeDEpsTrace,
                 probeDEpsNorm,
                 probePressurePlus,
                 probePressureMinus );

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
