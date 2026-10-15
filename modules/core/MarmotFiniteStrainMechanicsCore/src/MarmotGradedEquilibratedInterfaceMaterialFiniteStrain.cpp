#include "Marmot/MarmotGradedEquilibratedInterfaceMaterialFiniteStrain.h"

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotTwoSidedGradedLayer.h"
#include "Marmot/MarmotTwoSidedInterfaceGeometry.h"

#include <Eigen/Dense>

#include <algorithm>
#include <stdexcept>
#include <vector>

using namespace Marmot;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::Materials::TwoSidedInterface;
using namespace Marmot::Materials::TwoSidedGradedLayer;

void MarmotGradedEquilibratedInterfaceMaterialFiniteStrain::computeStress( State&               state,
                                                                           Tangents&            tangents,
                                                                           const Deformation&   deformation,
                                                                           const TimeIncrement& timeIncrement )
{
  Eigen::Vector3d normal = Eigen::Map< const Eigen::Vector3d >( deformation.normal.data() );
  if ( normal.norm() <= 1.0e-12 )
    throw std::invalid_argument( "MarmotGradedEquilibratedInterfaceMaterialFiniteStrain: interface normal is zero." );
  normal.normalize();

  const InterfaceGeometry geometry = evaluateInterfaceGeometry( normal,
                                                                Eigen::Map< const Eigen::Vector3d >(
                                                                  deformation.connector.data() ),
                                                                h );

  const int nq = numberOfThicknessPoints();

  Layer< 9 > layer;
  layer.normal               = normal;
  layer.normalSeparation     = geometry.normalSeparation;
  layer.tangentialSeparation = geometry.tangentialSeparation;
  layer.thickness            = h;
  layer.cumulativeWeights    = &cumulativeWeights;
  layer.tangentialGradient   = tangentialGradient;
  layer.gradientToStrain     = Eigen::Matrix< double, 9, 9 >::Identity();
  layer.stressToVector       = Eigen::Matrix< double, 9, 9 >::Identity();
  layer.strainOffset         = Eigen::Matrix< double, 9, 1 >::Zero();
  for ( int i = 0; i < 3; i++ )
    layer.strainOffset( 3 * i + i ) = 1.0;
  for ( int q = 0; q < nq; q++ ) {
    layer.weights.push_back( points[q].weight );
    layer.positions.push_back( points[q].position );
  }

  // generalized strain x = ( w, A+, A- ), total
  Vector21d x;
  {
    const Eigen::Map< const Eigen::Matrix< double, 6, 1 > > U( deformation.U.data() );
    x.segment< 3 >( 0 )  = U.segment< 3 >( 0 ) - U.segment< 3 >( 3 );
    x.segment< 9 >( 3 )  = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( deformation.surfaceGradient.data() );
    x.segment< 9 >( 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( deformation.surfaceGradient.data() + 9 );
  }

  const MarmotMaterialFiniteStrain::TimeIncrement materialTimeIncrement{ timeIncrement.timeOld, timeIncrement.dT };

  const auto evaluateBulkPoint = [&]( int q, const Eigen::Matrix< double, 9, 1 >& F ) {
    auto bulk = evaluateBulk( q, state.stateVars, Tensor33d( F.data() ), materialTimeIncrement );

    PointResponse< 9 > response;
    response.stress = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( bulk.P.data() );
    response.D      = Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >( bulk.A.data() );
    response.state  = std::move( bulk.state );
    return response;
  };

  // the corrections of the last converged state start the iteration
  std::vector< Eigen::Vector3d > initialCorrections( nq );
  for ( int q = 0; q < nq; q++ )
    initialCorrections[q] = Eigen::Map< const Eigen::Vector3d >( state.stateVars + points[q].offset );

  const double tolerance = 1e-11 * x.segment< 18 >( 3 ).cwiseAbs().maxCoeff() + 1e-14;

  auto solution = solve< 9 >( layer, x, initialCorrections, tolerance, evaluateBulkPoint );

  // commit the corrections and the state of the converged evaluations
  for ( int q = 0; q < nq; q++ ) {
    auto& pointState = solution.points[q].state;
    std::copy( solution.points[q].w.data(), solution.points[q].w.data() + 3, pointState.begin() );
    std::copy( pointState.begin(), pointState.end(), state.stateVars + points[q].offset );
  }

  const Vector21d& p = solution.p;
  resultant_         = solution.resultant;
  const Matrix21d& K = solution.K;

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
