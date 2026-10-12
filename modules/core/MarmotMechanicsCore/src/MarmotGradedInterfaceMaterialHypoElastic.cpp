#include "Marmot/MarmotGradedInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotThicknessMicroproblem.h"
#include "Marmot/MarmotThicknessProfile.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtils.h"
#include "Marmot/MarmotVoigt.h"

#include "Fastor/Fastor.h"
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>

using namespace Marmot;
using namespace Marmot::FastorIndices;
using namespace Marmot::FastorStandardTensors;

void MarmotGradedInterfaceMaterialHypoElastic::computeStress( State&               state,
                                                              Tangents&            tangents,
                                                              const Deformation&   deformation,
                                                              const TimeIncrement& timeIncrement )
{
  using namespace Marmot::Materials::InterfaceMaterialHelperFunctions;

  const auto& normal         = deformation.normal;
  auto        dU             = deformation.dU;
  auto        dSurfaceStrain = deformation.dSurfaceStrain;

  // strain of the macroscopic kinematics, as in MarmotInterfaceMaterialHypoElastic
  const Tensor3d jumpU = dU( Fastor::seq( 0, 3 ) ) - dU( Fastor::seq( 3, Fastor::last ) );

  const Tensor9d  averageSurfaceGradientFlat = 0.5 * ( dSurfaceStrain( Fastor::seq( 0, 9 ) ) +
                                                      dSurfaceStrain( Fastor::seq( 9, Fastor::last ) ) );
  const Tensor33d averageSurfaceGradient( Fastor::reshape< 3, 3 >( averageSurfaceGradientFlat ) );
  const Tensor33d displacementGradient = ( 1. / h ) * Fastor::einsum< i, j, to_ij >( jumpU, normal ) +
                                         averageSurfaceGradient;
  const Tensor33d strainIncrement = 0.5 * ( displacementGradient + Fastor::transpose( displacementGradient ) );

  const Vector6d dE0 = ContinuumMechanics::VoigtNotation::strainToVoigt( mapEigenToFastor( strainIncrement ) );

  const Eigen::Vector3d n( normal( 0 ), normal( 1 ), normal( 2 ) );
  const auto            S  = normalStrainOperator( n );
  const auto            St = S.transpose().eval();

  const int nq = static_cast< int >( points.size() );

  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld + timeIncrement.dT, timeIncrement.dT };

  // Evaluate the point q for the trial value w of the correction, starting from the state of the previous increment.
  const auto evaluate = [&]( int q, const Eigen::Vector3d& w ) {
    auto bulk = evaluateBulk( q, state.stateVars, dE0 + S * w, timeInfo );

    PointEvaluation evaluation;
    evaluation.w     = w;
    evaluation.sigma = bulk.sigma;
    evaluation.D     = bulk.D;
    evaluation.state = std::move( bulk.state );
    evaluation.a     = St * evaluation.sigma;

    evaluation.Qinv = ( St * evaluation.D * S ).inverse();
    if ( !evaluation.Qinv.allFinite() )
      throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__
                                             << ": the acoustic tensor of an integration point through the "
                                                "thickness is singular" );
    return evaluation;
  };

  std::vector< double > weights( nq );
  for ( int q = 0; q < nq; q++ )
    weights[q] = points[q].weight;

  const double tolerance = 1e-11 * dE0.norm() + 1e-15;

  const std::vector< Eigen::Vector3d > initialCorrections( nq, Eigen::Vector3d::Zero() );
  const auto evaluations = Marmot::Materials::ThicknessMicroproblem::solve< PointEvaluation >( weights,
                                                                                               evaluate,
                                                                                               initialCorrections,
                                                                                               tolerance );

  // commit the state of the converged evaluation
  for ( int q = 0; q < nq; q++ )
    std::copy( evaluations[q].state.begin(), evaluations[q].state.end(), state.stateVars + points[q].offset );

  // condensed tangent of the layer and through-thickness mean of the stress
  std::vector< Eigen::Matrix< double, 3, 6 > > B( nq );
  Vector6d                                     meanStress = Vector6d::Zero();
  for ( int q = 0; q < nq; q++ ) {
    B[q] = St * evaluations[q].D;
    meanStress += weights[q] * evaluations[q].sigma;
  }
  const auto dw = Marmot::Materials::ThicknessMicroproblem::condensedCorrectionTangents< 6 >( weights, evaluations, B );

  Matrix6d tangent = Matrix6d::Zero();
  for ( int q = 0; q < nq; q++ )
    tangent += weights[q] * evaluations[q].D * ( Matrix6d::Identity() + S * dw[q] );

  auto [Z, Q, H, Y] = calculateInterfaceMaterialParameters( normal, tangent );
  tangents.Q_ij     = ( 1. / h ) * Q;
  tangents.Z_ijkl   = h * Z;
  tangents.H_ijk    = H;
  tangents.Y_ijkl   = h * Y;

  mapEigenToFastor( state.surfaceStress ) = h * ContinuumMechanics::VoigtNotation::voigtToStress( meanStress );
  state.force                             = ( 1. / h ) * Fastor::einsum< ij, j, to_i >( state.surfaceStress, normal );
}
