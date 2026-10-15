#include "Marmot/MarmotGradedXInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotTwoSidedGradedLayer.h"
#include "Marmot/MarmotTwoSidedInterfaceGeometry.h"
#include "Marmot/MarmotTypedefs.h"

#include <Eigen/Dense>

#include <algorithm>
#include <stdexcept>
#include <vector>

using namespace Marmot;
using namespace Marmot::Materials::TwoSidedInterface;
using namespace Marmot::Materials::TwoSidedGradedLayer;

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

  const int nq = numberOfThicknessPoints();

  Layer< 6 > layer;
  layer.normal               = normal;
  layer.normalSeparation     = geometry.normalSeparation;
  layer.tangentialSeparation = geometry.tangentialSeparation;
  layer.thickness            = h;
  layer.cumulativeWeights    = &cumulativeWeights;
  layer.tangentialGradient   = tangentialGradient;
  layer.gradientToStrain     = gradientToVoigtStrain();
  layer.stressToVector       = voigtStressToVector();
  layer.strainOffset         = Vector6d::Zero();
  for ( int q = 0; q < nq; q++ ) {
    layer.weights.push_back( points[q].weight );
    layer.positions.push_back( points[q].position );
  }

  // generalized strain increment x = ( w, A+, A- )
  Vector21d x;
  {
    const Eigen::Map< const Eigen::Matrix< double, 6, 1 > > dU( deformation.dU.data() );
    x.segment< 3 >( 0 )  = dU.segment< 3 >( 0 ) - dU.segment< 3 >( 3 );
    x.segment< 9 >( 3 )  = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( deformation.dSurfaceStrain.data() );
    x.segment< 9 >( 12 ) = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( deformation.dSurfaceStrain.data() + 9 );
  }

  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld + timeIncrement.dT, timeIncrement.dT };

  const auto evaluateBulkPoint = [&]( int q, const Vector6d& strain ) {
    auto               bulk = evaluateBulk( q, state.stateVars, strain, timeInfo );
    PointResponse< 6 > response;
    response.stress = bulk.sigma;
    response.D      = bulk.D;
    response.state  = std::move( bulk.state );
    return response;
  };

  // the corrections of the increment start from zero
  const std::vector< Eigen::Vector3d > initialCorrections( nq, Eigen::Vector3d::Zero() );

  // the corrections are of the order of the strain of the layer
  const double tolerance = 1e-11 * x.cwiseAbs().maxCoeff() + 1e-15;

  auto solution = solve< 6 >( layer, x, initialCorrections, tolerance, evaluateBulkPoint );

  // commit the state of the converged evaluations
  for ( int q = 0; q < nq; q++ )
    std::copy( solution.points[q].state.begin(), solution.points[q].state.end(), state.stateVars + points[q].offset );

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
