#include "Marmot/MarmotGradedInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
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

namespace {

  constexpr int maxMicroIterations = 60; ///< iterations of the microproblem
  constexpr int maxBacktracks      = 12; ///< halvings of a step after a failed update of a bulk material

  /**
   * @brief Nodes and weights of the Gauss-Legendre rule on [-1, 1], found by Newton's method on the Legendre
   * polynomial.
   *
   * @param[in] n number of points
   * @return the nodes and the weights, which sum up to 2
   */
  std::pair< std::vector< double >, std::vector< double > > gaussLegendre( int n )
  {
    std::vector< double > nodes( n ), weights( n );
    for ( int i = 0; i < n; i++ ) {
      double x  = std::cos( M_PI * ( i + 0.75 ) / ( n + 0.5 ) );
      double dP = 1.0;
      for ( int iteration = 0; iteration < 100; iteration++ ) {
        double p0 = 1.0, p1 = x;
        for ( int k = 2; k <= n; k++ ) {
          const double p2 = ( ( 2. * k - 1. ) * x * p1 - ( k - 1. ) * p0 ) / k;
          p0              = p1;
          p1              = p2;
        }
        const double pn = ( n == 0 ) ? 1.0 : p1;
        const double pm = ( n == 1 ) ? 1.0 : p0;
        dP              = n * ( x * pn - pm ) / ( x * x - 1. );
        const double dx = pn / dP;
        x -= dx;
        if ( std::abs( dx ) < 1e-15 )
          break;
      }
      nodes[i]   = x;
      weights[i] = 2. / ( ( 1. - x * x ) * dP * dP );
    }
    return { nodes, weights };
  }

  /**
   * @brief Operator S that maps a vector w to the strain sym( w (x) n ) in Voigt notation.
   *
   * With the engineering shear strains of the Voigt notation, the transpose of S maps a stress in Voigt notation
   * to the traction sigma n.
   *
   * @param[in] n the normal
   * @return the 6x3 operator
   */
  Eigen::Matrix< double, 6, 3 > normalStrainOperator( const Eigen::Vector3d& n )
  {
    Eigen::Matrix< double, 6, 3 > S;
    for ( int j = 0; j < 3; j++ ) {
      Eigen::Vector3d e            = Eigen::Vector3d::Zero();
      e( j )                       = 1.0;
      const Eigen::Matrix3d strain = 0.5 * ( e * n.transpose() + n * e.transpose() );
      S.col( j )                   = ContinuumMechanics::VoigtNotation::strainToVoigt( strain );
    }
    return S;
  }

  /// Response of an integration point through the thickness to the strain of a trial value of w.
  struct PointEvaluation {
    Eigen::Vector3d       w = Eigen::Vector3d::Zero(); ///< first-order strain correction w = du1/dy
    Vector6d              sigma;                       ///< stress
    Matrix6d              D;                           ///< tangent
    Eigen::Matrix3d       Qinv;                        ///< inverse of the acoustic tensor
    Eigen::Vector3d       a;                           ///< traction sigma n
    std::vector< double > state;                       ///< stress and state variables after the update
  };

} // namespace

MarmotGradedInterfaceMaterialHypoElastic::MarmotGradedInterfaceMaterialHypoElastic( const std::string& materialName,
                                                                                    const double*      matProperties,
                                                                                    int nMaterialProperties,
                                                                                    int materialNumber_ )
  : materialNumber( materialNumber_ ), h( checkedMaterialProperty( matProperties, nMaterialProperties, 0 ) )
{
  if ( h <= 0. )
    throw std::invalid_argument( "MarmotGradedInterfaceMaterialHypoElastic: the thickness must be positive" );

  const int mode   = static_cast< int >( checkedMaterialProperty( matProperties, nMaterialProperties, 1 ) );
  const int nGauss = static_cast< int >( checkedMaterialProperty( matProperties, nMaterialProperties, 2 ) );
  const int K      = static_cast< int >( checkedMaterialProperty( matProperties, nMaterialProperties, 3 ) );
  const int m      = static_cast< int >( checkedMaterialProperty( matProperties, nMaterialProperties, 4 ) );

  if ( mode != 0 && mode != 1 )
    throw std::invalid_argument( "MarmotGradedInterfaceMaterialHypoElastic: the mode must be 0 (layers) or 1 "
                                 "(interpolated stations)" );
  if ( nGauss < 1 || m < 1 || K < ( mode == 0 ? 1 : 2 ) )
    throw std::invalid_argument( "MarmotGradedInterfaceMaterialHypoElastic: nGauss and m must be positive, and K at "
                                 "least one layer or two stations" );
  if ( nMaterialProperties != 5 + K * ( 1 + m ) )
    throw std::invalid_argument( "MarmotGradedInterfaceMaterialHypoElastic: " + std::to_string( 5 + K * ( 1 + m ) ) +
                                 " material properties are required for K = " + std::to_string( K ) +
                                 " and m = " + std::to_string( m ) + ", but " + std::to_string( nMaterialProperties ) +
                                 " are given" );

  const auto position = [&]( int k ) { return matProperties[5 + k * ( 1 + m )]; };
  const auto property = [&]( int k, int i ) { return matProperties[5 + k * ( 1 + m ) + 1 + i]; };

  // segments [a, b] of the thickness, with the properties at both ends
  struct Segment {
    double a, b;
    int    kA, kB; // stations at both ends, kA == kB for a layer
  };
  std::vector< Segment > segments;

  if ( mode == 0 ) {
    double sum = 0.;
    for ( int k = 0; k < K; k++ ) {
      if ( position( k ) <= 0. )
        throw std::invalid_argument( "MarmotGradedInterfaceMaterialHypoElastic: the layer thicknesses must be "
                                     "positive" );
      sum += position( k );
    }
    double a = 0.;
    for ( int k = 0; k < K; k++ ) {
      const double b = ( k == K - 1 ) ? 1.0 : a + position( k ) / sum;
      segments.push_back( { a, b, k, k } );
      a = b;
    }
  }
  else {
    if ( position( 0 ) != 0. || position( K - 1 ) != 1. )
      throw std::invalid_argument( "MarmotGradedInterfaceMaterialHypoElastic: the first station must be at 0 and "
                                   "the last at 1" );
    for ( int k = 0; k < K - 1; k++ ) {
      if ( position( k + 1 ) <= position( k ) )
        throw std::invalid_argument( "MarmotGradedInterfaceMaterialHypoElastic: the stations must be increasing" );
      segments.push_back( { position( k ), position( k + 1 ), k, k + 1 } );
    }
  }

  const auto [nodes, weights] = gaussLegendre( nGauss );

  std::size_t offset = 0;
  for ( const auto& segment : segments ) {
    for ( int g = 0; g < nGauss; g++ ) {
      ThicknessPoint point;
      point.weight   = 0.5 * weights[g] * ( segment.b - segment.a );
      point.position = segment.a + 0.5 * ( nodes[g] + 1. ) * ( segment.b - segment.a );

      point.properties.resize( m );
      const double s = ( segment.kA == segment.kB ) ? 0. : ( point.position - segment.a ) / ( segment.b - segment.a );
      for ( int i = 0; i < m; i++ )
        point.properties[i] = ( 1. - s ) * property( segment.kA, i ) + s * property( segment.kB, i );

      point.material = std::unique_ptr< MarmotMaterialHypoElastic >(
        MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                         point.properties.data(),
                                                                         m,
                                                                         materialNumber ) );
      if ( !point.material )
        throw std::invalid_argument( "Unknown base material for MarmotGradedInterfaceMaterialHypoElastic: " +
                                     materialName );

      point.nBaseStateVars = point.material->getNumberOfRequiredStateVars();
      point.offset         = offset;
      offset += 6 + point.nBaseStateVars;
      points.push_back( std::move( point ) );
    }
  }

  for ( std::size_t q = 0; q < points.size(); q++ ) {
    stateLayout.add( "stress" + std::to_string( q ), 6 );
    stateLayout.add( "baseMaterialStateVars" + std::to_string( q ), points[q].nBaseStateVars );
  }
  stateLayout.finalize();
}

void MarmotGradedInterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  for ( auto& point : points )
    point.material->setCharacteristicElementLength( length );
}

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

  const double tolerance = 1e-11 * dE0.norm() + 1e-15;

  // Evaluate the point q for the trial value w of the correction, starting from the state of the previous increment.
  const auto evaluate = [&]( int q, const Eigen::Vector3d& w ) {
    const auto& point = points[q];

    PointEvaluation evaluation;
    evaluation.w = w;
    evaluation.state.assign( state.stateVars + point.offset,
                             state.stateVars + point.offset + 6 + point.nBaseStateVars );

    MarmotMaterialHypoElastic::state3D baseState;
    baseState.stress    = Eigen::Map< const Vector6d >( evaluation.state.data() );
    baseState.stateVars = evaluation.state.data() + 6;

    evaluation.D.setZero();
    point.material->computeStress( baseState, evaluation.D, dE0 + S * w, timeInfo );

    Eigen::Map< Vector6d > storedStress( evaluation.state.data() );
    storedStress     = baseState.stress;
    evaluation.sigma = baseState.stress;
    evaluation.a     = St * evaluation.sigma;

    evaluation.Qinv = ( St * evaluation.D * S ).inverse();
    if ( !evaluation.Qinv.allFinite() )
      throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__
                                             << ": the acoustic tensor of an integration point through the "
                                                "thickness is singular" );
    return evaluation;
  };

  const auto evaluateAll =
    [&]( const std::vector< PointEvaluation >& base, double alpha, const std::vector< Eigen::Vector3d >& steps ) {
      std::vector< PointEvaluation > result( nq );
      for ( int q = 0; q < nq; q++ )
        result[q] = evaluate( q, base[q].w + alpha * steps[q] );
      return result;
    };

  std::vector< PointEvaluation > evaluations( nq );
  for ( int q = 0; q < nq; q++ )
    evaluations[q] = evaluate( q, Eigen::Vector3d::Zero() );

  // The unknowns are the corrections w of all points and the common traction t. The residual consists of the
  // misfit a_q - t of the tractions and of the mean of w, which has to vanish. Its Jacobian has the blocks Q_q
  // and -1, and the weights, so that the Newton step condenses to the formulas of the class description. The step
  // is damped by a backtracking line search on the squared norm of the residual, where the tractions are scaled
  // to strains with the stiffest acoustic tensor.
  double stiffness = 0.;
  for ( const auto& evaluation : evaluations )
    stiffness = std::max( stiffness, evaluation.Qinv.inverse().cwiseAbs().maxCoeff() );

  Eigen::Matrix3d M   = Eigen::Matrix3d::Zero();
  Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
  for ( int q = 0; q < nq; q++ ) {
    M += points[q].weight * evaluations[q].Qinv;
    rhs += points[q].weight * evaluations[q].Qinv * evaluations[q].a;
  }
  Eigen::Vector3d t = M.partialPivLu().solve( rhs );

  const auto merit = [&]( const std::vector< PointEvaluation >& trial, const Eigen::Vector3d& tTrial ) {
    double          sum  = 0.;
    Eigen::Vector3d mean = Eigen::Vector3d::Zero();
    for ( int q = 0; q < nq; q++ ) {
      sum += points[q].weight * ( ( trial[q].a - tTrial ) / stiffness ).squaredNorm();
      mean += points[q].weight * trial[q].w;
    }
    return 0.5 * ( sum + mean.squaredNorm() );
  };

  bool converged = false;
  for ( int iteration = 0; iteration < maxMicroIterations && !converged; iteration++ ) {

    Eigen::Vector3d mean = Eigen::Vector3d::Zero();
    M.setZero();
    rhs.setZero();
    for ( int q = 0; q < nq; q++ ) {
      mean += points[q].weight * evaluations[q].w;
      M += points[q].weight * evaluations[q].Qinv;
      rhs += points[q].weight * evaluations[q].Qinv * ( evaluations[q].a - t );
    }
    const Eigen::Vector3d dt = M.partialPivLu().solve( rhs - mean );

    std::vector< Eigen::Vector3d > steps( nq );
    double                         maxStep = 0.;
    for ( int q = 0; q < nq; q++ ) {
      steps[q] = evaluations[q].Qinv * ( dt - ( evaluations[q].a - t ) );
      maxStep  = std::max( maxStep, steps[q].cwiseAbs().maxCoeff() );
    }

    if ( maxStep <= tolerance ) {
      converged = true;
      break;
    }

    const double meritBefore = merit( evaluations, t );

    double alpha = 1.0;
    for ( int backtrack = 0;; backtrack++ ) {
      try {
        auto trial = evaluateAll( evaluations, alpha, steps );
        if ( merit( trial, t + alpha * dt ) <= ( 1. - 2e-4 * alpha ) * meritBefore || backtrack >= maxBacktracks ) {
          evaluations = std::move( trial );
          t += alpha * dt;
          break;
        }
      }
      catch ( const std::exception& ) {
        // a bulk material could not update the strain of this trial: take a shorter step
        if ( backtrack >= maxBacktracks )
          throw;
      }
      alpha *= 0.5;
    }
  }

  if ( !converged )
    throw StressUpdateFailed( MakeString()
                              << __PRETTY_FUNCTION__ << ": the microproblem through the thickness did not converge" );

  // commit the state of the converged evaluation
  for ( int q = 0; q < nq; q++ )
    std::copy( evaluations[q].state.begin(), evaluations[q].state.end(), state.stateVars + points[q].offset );

  // condensed tangent of the layer and through-thickness mean of the stress
  M.setZero();
  Eigen::Matrix< double, 3, 6 > G          = Eigen::Matrix< double, 3, 6 >::Zero();
  Vector6d                      meanStress = Vector6d::Zero();
  for ( int q = 0; q < nq; q++ ) {
    const double                        weight = points[q].weight;
    const Eigen::Matrix< double, 3, 6 > B      = St * evaluations[q].D;
    M += weight * evaluations[q].Qinv;
    G += weight * evaluations[q].Qinv * B;
    meanStress += weight * evaluations[q].sigma;
  }
  const Eigen::Matrix< double, 3, 6 > dt = M.partialPivLu().solve( G );

  Matrix6d tangent = Matrix6d::Zero();
  for ( int q = 0; q < nq; q++ ) {
    const Eigen::Matrix< double, 3, 6 > B  = St * evaluations[q].D;
    const Eigen::Matrix< double, 3, 6 > dw = evaluations[q].Qinv * ( dt - B );
    tangent += points[q].weight * evaluations[q].D * ( Matrix6d::Identity() + S * dw );
  }

  auto [Z, Q, H, Y] = calculateInterfaceMaterialParameters( normal, tangent );
  tangents.Q_ij     = ( 1. / h ) * Q;
  tangents.Z_ijkl   = h * Z;
  tangents.H_ijk    = H;
  tangents.Y_ijkl   = h * Y;

  mapEigenToFastor( state.surfaceStress ) = h * ContinuumMechanics::VoigtNotation::voigtToStress( meanStress );
  state.force                             = ( 1. / h ) * Fastor::einsum< ij, j, to_i >( state.surfaceStress, normal );
}

void MarmotGradedInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int )
{
  for ( auto& point : points ) {
    std::fill( stateVars + point.offset, stateVars + point.offset + 6, 0.0 );
    point.material->initializeYourself( stateVars + point.offset + 6, static_cast< int >( point.nBaseStateVars ) );
  }
}

double MarmotGradedInterfaceMaterialHypoElastic::getDensity()
{
  double density = 0.;
  for ( auto& point : points )
    density += point.weight * point.material->getDensity( nullptr );
  return density;
}
