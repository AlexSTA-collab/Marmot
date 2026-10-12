#include "Marmot/MarmotGradedXInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotThicknessMicroproblem.h"
#include "Marmot/MarmotTwoSidedInterfaceGeometry.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::Materials::TwoSidedInterface;

namespace {

  using Matrix3x6  = Eigen::Matrix< double, 3, 6 >;
  using Matrix3x9  = Eigen::Matrix< double, 3, 9 >;
  using Matrix3x21 = Eigen::Matrix< double, 3, 21 >;
  using Matrix6x9  = Eigen::Matrix< double, 6, 9 >;
  using Matrix6x21 = Eigen::Matrix< double, 6, 21 >;
  using Matrix9x21 = Eigen::Matrix< double, 9, 21 >;
  using Matrix21d  = Eigen::Matrix< double, 21, 21 >;
  using Vector21d  = Eigen::Matrix< double, 21, 1 >;

  /// Maximum number of passes of the iteration for the derivative of the correction of the closure.
  constexpr int maxClosurePasses = 8;

  /// Largest change of the derivative of the correction, relative to its norm, which ends the closure iteration.
  constexpr double closureTolerance = 1e-10;

  /**
   * @brief Terms of the closure of the tangential gradient of the correction at all integration points.
   */
  struct ClosureTerms {
    std::vector< Vector6d >        strain; ///< strain of the tangential gradient of the correction, in Voigt notation
    std::vector< Eigen::Vector3d > force;  ///< change of the traction from the bottom to the point
  };

} // namespace

void MarmotGradedXInterfaceMaterialHypoElastic::computeStress( State&               state,
                                                               Tangents&            tangents,
                                                               const Deformation&   deformation,
                                                               const TimeIncrement& timeIncrement )
{
  Eigen::Vector3d normal = Eigen::Map< const Eigen::Vector3d >( deformation.normal.data() );
  if ( normal.norm() <= 1.0e-12 )
    throw std::invalid_argument( "MarmotGradedXInterfaceMaterialHypoElastic: interface normal is zero." );
  normal.normalize();

  const InterfaceGeometry geometry = evaluateInterfaceGeometry( normal,
                                                                Eigen::Map< const Eigen::Vector3d >(
                                                                  deformation.separationVector.data() ),
                                                                h );
  const double            ell      = geometry.normalSeparation;

  // generalized strain increment x = ( w, A+, A- )
  Vector21d x;
  {
    const Eigen::Map< const Eigen::Matrix< double, 6, 1 > > dU( deformation.dU.data() );
    x.segment< 3 >( 0 )  = dU.segment< 3 >( 0 ) - dU.segment< 3 >( 3 );
    x.segment< 9 >( 3 )  = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( deformation.dSurfaceStrain.data() );
    x.segment< 9 >( 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( deformation.dSurfaceStrain.data() + 9 );
  }
  const Eigen::Matrix3d dA = Eigen::Map< const Matrix3dRowMajor >( x.data() + 3 ) -
                             Eigen::Map< const Matrix3dRowMajor >( x.data() + 12 );

  const auto S  = normalStrainOperator( normal );
  const auto St = S.transpose().eval();
  const auto E  = gradientToVoigtStrain();
  const auto Ts = voigtStressToVector();

  const Matrix9x9 I9    = Matrix9x9::Identity();
  const Matrix9x9 Hsrf  = surfaceToGradient( normal, geometry.tangentialSeparation, ell ) - I9;
  const Matrix9x3 Bw    = jumpToGradient( normal, ell );
  const int       nq    = numberOfThicknessPoints();
  const double    thick = h;

  // gradient of an integration point as a function of the generalized strain, and its strain
  std::vector< Matrix9x21 > Bx( nq );
  std::vector< Matrix6x21 > EBx( nq );
  std::vector< Vector6d >   strain( nq );
  std::vector< double >     weights( nq );
  for ( int q = 0; q < nq; q++ ) {
    const double y = points[q].position;
    Bx[q].setZero();
    Bx[q].block< 9, 3 >( 0, 0 )  = Bw;
    Bx[q].block< 9, 9 >( 0, 3 )  = y * I9 + 0.5 * Hsrf;
    Bx[q].block< 9, 9 >( 0, 12 ) = ( 1. - y ) * I9 + 0.5 * Hsrf;
    EBx[q]                       = E * Bx[q];
    strain[q]                    = EBx[q] * x;
    weights[q]                   = points[q].weight;
  }

  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld + timeIncrement.dT, timeIncrement.dT };

  // Terms of the closure for the tangential gradient of the correction. They are linear in dA, the difference of the
  // surface gradients, for the derivative of the correction w with respect to the strain and the tangent of the
  // points.
  const auto closure = [&]( const std::vector< Matrix3x6 >& dw,
                            const std::vector< Matrix6d >&  D,
                            const Eigen::Matrix3d&          surfaceGradientDifference ) {
    ClosureTerms terms;
    terms.strain.assign( nq, Vector6d::Zero() );
    terms.force.assign( nq, Eigen::Vector3d::Zero() );

    // tangential derivative of the strain along e_r, and the response of the corrections and the stresses to it
    std::vector< Eigen::Matrix3d > derivativeOfCorrection( nq, Eigen::Matrix3d::Zero() ); // column r: d_r w_q
    std::vector< Eigen::Vector3d > source( nq, Eigen::Vector3d::Zero() );
    for ( int r = 0; r < 3; r++ ) {
      const Eigen::Vector3d column         = surfaceGradientDifference.col( r );
      const Eigen::Matrix3d strainGradient = 0.5 * ( column * normal.transpose() + normal * column.transpose() ) / ell;
      const Vector6d        s              = ContinuumMechanics::VoigtNotation::strainToVoigt( strainGradient );
      for ( int q = 0; q < nq; q++ ) {
        const Eigen::Vector3d dwr          = dw[q] * s;
        derivativeOfCorrection[q].col( r ) = dwr;
        const Eigen::Matrix3d dSigma       = ContinuumMechanics::VoigtNotation::voigtToStress(
          Vector6d( D[q] * ( s + S * dwr ) ) );
        source[q] += thick * dSigma.col( r );
      }
    }

    // cumulative integrals from the bottom: the force of the traction and the tangential gradient of u1
    for ( int q = 0; q < nq; q++ ) {
      Eigen::Matrix3d gradient = Eigen::Matrix3d::Zero();
      for ( int p = 0; p < nq; p++ ) {
        terms.force[q] += cumulativeWeights( q, p ) * source[p];
        gradient += cumulativeWeights( q, p ) * thick * derivativeOfCorrection[p];
      }
      const Eigen::Matrix3d symmetric = 0.5 * ( gradient + gradient.transpose() );
      terms.strain[q]                 = ContinuumMechanics::VoigtNotation::strainToVoigt( symmetric );
    }
    return terms;
  };

  // evaluate the integration point q for the correction w, with the closure terms of the pass
  const auto makeEvaluate = [&]( const ClosureTerms& terms ) {
    return [&, terms]( int q, const Eigen::Vector3d& w ) {
      auto bulk = evaluateBulk( q, state.stateVars, strain[q] + terms.strain[q] + S * w, timeInfo );

      PointEvaluation evaluation;
      evaluation.w     = w;
      evaluation.sigma = bulk.sigma;
      evaluation.D     = bulk.D;
      evaluation.state = std::move( bulk.state );
      evaluation.a     = St * evaluation.sigma + terms.force[q];

      evaluation.Qinv = ( St * evaluation.D * S ).inverse();
      if ( !evaluation.Qinv.allFinite() )
        throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__
                                               << ": the acoustic tensor of an integration point through the "
                                                  "thickness is singular" );
      return evaluation;
    };
  };

  const auto condenseOverStrain = [&]( const std::vector< PointEvaluation >& evaluations ) {
    std::vector< Matrix3x6 > B( nq );
    for ( int q = 0; q < nq; q++ )
      B[q] = St * evaluations[q].D;
    return Marmot::Materials::ThicknessMicroproblem::condensedCorrectionTangents< 6 >( weights, evaluations, B );
  };

  double strainScale = 0.;
  for ( int q = 0; q < nq; q++ )
    strainScale = std::max( strainScale, strain[q].cwiseAbs().maxCoeff() );
  const double tolerance = 1e-11 * strainScale + 1e-15;

  const std::vector< Eigen::Vector3d > initialCorrections( nq, Eigen::Vector3d::Zero() );

  // the closure is iterated, because the derivative of the correction depends on the converged state
  const bool closureActive = tangentialGradient && dA.cwiseAbs().maxCoeff() > 0.;

  ClosureTerms                   terms;
  std::vector< Matrix3x6 >       dwStrain( nq );
  std::vector< Matrix6d >        tangentPoint( nq );
  std::vector< PointEvaluation > evaluations;

  if ( closureActive ) {
    // the derivative of the correction of the state at the beginning of the increment starts the iteration
    ClosureTerms none;
    none.strain.assign( nq, Vector6d::Zero() );
    none.force.assign( nq, Eigen::Vector3d::Zero() );

    const auto                     evaluateNone = makeEvaluate( none );
    std::vector< PointEvaluation > start;
    start.reserve( nq );
    for ( int q = 0; q < nq; q++ )
      start.push_back( evaluateNone( q, Eigen::Vector3d::Zero() ) );

    dwStrain = condenseOverStrain( start );
    for ( int q = 0; q < nq; q++ )
      tangentPoint[q] = start[q].D;
    terms = closure( dwStrain, tangentPoint, dA );
  }
  else {
    terms.strain.assign( nq, Vector6d::Zero() );
    terms.force.assign( nq, Eigen::Vector3d::Zero() );
  }

  for ( int pass = 0; pass < maxClosurePasses; pass++ ) {
    evaluations = Marmot::Materials::ThicknessMicroproblem::solve< PointEvaluation >( weights,
                                                                                      makeEvaluate( terms ),
                                                                                      initialCorrections,
                                                                                      tolerance );
    dwStrain    = condenseOverStrain( evaluations );
    for ( int q = 0; q < nq; q++ )
      tangentPoint[q] = evaluations[q].D;

    if ( !closureActive )
      break;

    const ClosureTerms updated = closure( dwStrain, tangentPoint, dA );

    double change = 0.;
    double scale  = 0.;
    for ( int q = 0; q < nq; q++ ) {
      change = std::max( change, ( updated.strain[q] - terms.strain[q] ).cwiseAbs().maxCoeff() );
      change = std::max( change, ( updated.force[q] - terms.force[q] ).cwiseAbs().maxCoeff() );
      scale  = std::max( scale, updated.strain[q].cwiseAbs().maxCoeff() );
      scale  = std::max( scale, updated.force[q].cwiseAbs().maxCoeff() );
    }
    terms = updated;
    if ( change <= closureTolerance * scale + 1e-300 )
      break;
  }

  // commit the state of the converged evaluations
  for ( int q = 0; q < nq; q++ )
    std::copy( evaluations[q].state.begin(), evaluations[q].state.end(), state.stateVars + points[q].offset );

  // linear dependence of the closure terms on the surface gradients A+ and A-, with the derivative of the correction
  // of the converged state frozen
  std::vector< Matrix6x21 > closureStrainTangent( nq, Matrix6x21::Zero() );
  std::vector< Matrix3x21 > closureForceTangent( nq, Matrix3x21::Zero() );
  for ( int k = 0; k < ( tangentialGradient ? 9 : 0 ); k++ ) {
    Eigen::Matrix3d unit        = Eigen::Matrix3d::Zero();
    unit( k / 3, k % 3 )        = 1.0;
    const ClosureTerms response = closure( dwStrain, tangentPoint, unit );
    for ( int q = 0; q < nq; q++ ) {
      closureStrainTangent[q].col( 3 + k )  = response.strain[q];
      closureStrainTangent[q].col( 12 + k ) = -response.strain[q];
      closureForceTangent[q].col( 3 + k )   = response.force[q];
      closureForceTangent[q].col( 12 + k )  = -response.force[q];
    }
  }

  // condensation of the corrections over the generalized strain
  std::vector< Matrix3x21 > B( nq );
  for ( int q = 0; q < nq; q++ )
    B[q] = St * evaluations[q].D * ( EBx[q] + closureStrainTangent[q] ) + closureForceTangent[q];
  const auto dw = Marmot::Materials::ThicknessMicroproblem::condensedCorrectionTangents< 21 >( weights,
                                                                                               evaluations,
                                                                                               B );

  Vector21d p = Vector21d::Zero();
  Matrix21d K = Matrix21d::Zero();
  for ( int q = 0; q < nq; q++ ) {
    const double                         scale  = thick * weights[q];
    const Eigen::Matrix< double, 9, 21 > dSigma = Ts * evaluations[q].D *
                                                  ( EBx[q] + closureStrainTangent[q] + S * dw[q] );
    p += scale * Bx[q].transpose() * ( Ts * evaluations[q].sigma );
    K += scale * Bx[q].transpose() * dSigma;
  }

  state.generalizedForce   = Tensor3d( p.data() );
  state.surfaceStressPlus  = Tensor33d( p.data() + 3 );
  state.surfaceStressMinus = Tensor33d( p.data() + 12 );

  mapEigenToFastor( tangents.Q_ww )  = K.block< 3, 3 >( 0, 0 );
  mapEigenToFastor( tangents.Q_wAp ) = K.block< 3, 9 >( 0, 3 );
  mapEigenToFastor( tangents.Q_wAm ) = K.block< 3, 9 >( 0, 12 );

  // the blocks of the stress of a face with respect to the jump are 9x3 matrices, stored as third order tensors
  Eigen::Map< Eigen::Matrix< double, 9, 3, Eigen::RowMajor > >( tangents.Q_Apw.data() ) = K.block< 9, 3 >( 3, 0 );
  mapEigenToFastor( tangents.Q_ApAp )                                                   = K.block< 9, 9 >( 3, 3 );
  mapEigenToFastor( tangents.Q_ApAm )                                                   = K.block< 9, 9 >( 3, 12 );

  Eigen::Map< Eigen::Matrix< double, 9, 3, Eigen::RowMajor > >( tangents.Q_Amw.data() ) = K.block< 9, 3 >( 12, 0 );
  mapEigenToFastor( tangents.Q_AmAp )                                                   = K.block< 9, 9 >( 12, 3 );
  mapEigenToFastor( tangents.Q_AmAm )                                                   = K.block< 9, 9 >( 12, 12 );
}
