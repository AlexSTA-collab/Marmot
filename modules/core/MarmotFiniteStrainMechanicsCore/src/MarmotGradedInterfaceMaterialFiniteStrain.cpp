#include "Marmot/MarmotGradedInterfaceMaterialFiniteStrain.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceFiniteStrainHelperFunctions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"
#include "Marmot/MarmotThicknessMicroproblem.h"
#include "Marmot/MarmotThicknessProfile.h"
#include "Marmot/MarmotUtils.h"

#include "Fastor/Fastor.h"
#include <Eigen/Dense>

#include <algorithm>
#include <stdexcept>
#include <string>

using namespace Marmot;
using namespace Marmot::FastorIndices;
using namespace Marmot::FastorStandardTensors;

namespace {

  /// Response of an integration point through the thickness to the deformation gradient of a trial value of w.
  struct PointEvaluation {
    Eigen::Vector3d       w = Eigen::Vector3d::Zero(); ///< correction w of the normal derivative
    Tensor33d             P;                           ///< first Piola-Kirchhoff stress
    Tensor3333d           A;                           ///< derivative of the stress with respect to F
    Eigen::Matrix3d       Qinv;                        ///< inverse of the acoustic tensor A[N, N]
    Eigen::Vector3d       a;                           ///< traction P N
    std::vector< double > state;                       ///< correction and state variables after the update
  };

} // namespace

MarmotGradedInterfaceMaterialFiniteStrain::MarmotGradedInterfaceMaterialFiniteStrain( const std::string& materialName,
                                                                                      const double*      matProperties,
                                                                                      int nMaterialProperties,
                                                                                      int materialNumber_ )
  : materialNumber( materialNumber_ ), h( checkedMaterialProperty( matProperties, nMaterialProperties, 0 ) )
{
  auto profile = Marmot::Materials::ThicknessProfile::fromMaterialProperties( matProperties, nMaterialProperties );

  std::size_t offset = 0;
  for ( auto& profilePoint : profile.points ) {
    ThicknessPoint point;
    point.weight     = profilePoint.weight;
    point.position   = profilePoint.position;
    point.properties = std::move( profilePoint.properties );

    point.material = std::unique_ptr< MarmotMaterialFiniteStrain >(
      MarmotLibrary::MarmotMaterialFiniteStrainFactory::createMaterial( materialName,
                                                                        point.properties.data(),
                                                                        static_cast< int >( point.properties.size() ),
                                                                        materialNumber ) );
    if ( !point.material )
      throw std::invalid_argument( "Unknown base material for MarmotGradedInterfaceMaterialFiniteStrain: " +
                                   materialName );

    point.nBaseStateVars = point.material->getNumberOfRequiredStateVars();
    point.offset         = offset;
    offset += 3 + point.nBaseStateVars;
    points.push_back( std::move( point ) );
  }

  for ( std::size_t q = 0; q < points.size(); q++ ) {
    stateLayout.add( "correction" + std::to_string( q ), 3 );
    stateLayout.add( "baseMaterialStateVars" + std::to_string( q ), points[q].nBaseStateVars );
  }
  stateLayout.finalize();
}

void MarmotGradedInterfaceMaterialFiniteStrain::setCharacteristicElementLength( double length )
{
  for ( auto& point : points )
    point.material->setCharacteristicElementLength( length );
}

void MarmotGradedInterfaceMaterialFiniteStrain::computeStress( State&               state,
                                                               Tangents&            tangents,
                                                               const Deformation&   deformation,
                                                               const TimeIncrement& timeIncrement )
{
  using namespace Marmot::Materials::InterfaceFiniteStrainHelperFunctions;

  const Tensor3d& normal = deformation.normal;

  const NormalFibre fibre = normalFibre( deformation.U, deformation.surfaceGradient, normal, deformation.connector );

  const int nq = static_cast< int >( points.size() );

  const MarmotMaterialFiniteStrain::TimeIncrement materialTimeIncrement{ timeIncrement.timeOld, timeIncrement.dT };

  // Evaluate the point q for the trial value w of the correction, starting from the state of the previous increment.
  const auto evaluate = [&]( int q, const Eigen::Vector3d& w ) {
    const auto& point = points[q];

    PointEvaluation evaluation;
    evaluation.w = w;
    evaluation.state.assign( state.stateVars + point.offset,
                             state.stateVars + point.offset + 3 + point.nBaseStateVars );

    Tensor3d wTensor;
    for ( int i = 0; i < 3; i++ )
      wTensor( i ) = w( i );

    const Tensor33d F = fibre.F + einsum< i, j, to_ij >( wTensor, normal );

    MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ),
                                                                    0.0,
                                                                    0.0,
                                                                    evaluation.state.data() + 3 );
    MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 >    materialTangents;
    MarmotMaterialFiniteStrain::Deformation< 3 >          materialDeformation{ F };

    point.material->computeStress( response, materialTangents, materialDeformation, materialTimeIncrement );

    firstPiolaKirchhoff( response.tau, materialTangents.dTau_dF, F, evaluation.P, evaluation.A );

    const Tensor3d  traction = einsum< ij, j, to_i >( evaluation.P, normal );
    const Tensor33d ANN      = einsum< ijkl, j, l, to_ik >( evaluation.A, normal, normal );

    evaluation.a    = Eigen::Vector3d( traction( 0 ), traction( 1 ), traction( 2 ) );
    evaluation.Qinv = Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( ANN.data() ).inverse();

    if ( !evaluation.Qinv.allFinite() )
      throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__
                                             << ": the acoustic tensor of an integration point through the "
                                                "thickness is singular" );

    std::copy( w.data(), w.data() + 3, evaluation.state.begin() );
    return evaluation;
  };

  std::vector< double >          weights( nq );
  std::vector< Eigen::Vector3d > initialCorrections( nq );
  for ( int q = 0; q < nq; q++ ) {
    weights[q]            = points[q].weight;
    initialCorrections[q] = Eigen::Map< const Eigen::Vector3d >( state.stateVars + points[q].offset );
  }

  // the corrections are of the order of the strain of the layer
  Tensor33d strain = fibre.F;
  for ( int i = 0; i < 3; i++ )
    strain( i, i ) -= 1.0;
  const double strainNorm = std::sqrt( Fastor::inner( strain, strain ) );
  const double tolerance  = 1e-11 * strainNorm + 1e-14;

  const auto evaluations = Marmot::Materials::ThicknessMicroproblem::solve< PointEvaluation >( weights,
                                                                                               evaluate,
                                                                                               initialCorrections,
                                                                                               tolerance );

  // commit the state of the converged evaluation
  for ( int q = 0; q < nq; q++ )
    std::copy( evaluations[q].state.begin(), evaluations[q].state.end(), state.stateVars + points[q].offset );

  // derivative of the traction of every point with respect to F^0, B(i, (k, l)) = A_ijkl N_j
  std::vector< Eigen::Matrix< double, 3, 9 > > B( nq );
  for ( int q = 0; q < nq; q++ ) {
    const Tensor333d Bq = einsum< ijkl, j, to_ikl >( evaluations[q].A, normal );
    B[q]                = Eigen::Map< const Eigen::Matrix< double, 3, 9, Eigen::RowMajor > >( Bq.data() );
  }

  const auto dw = Marmot::Materials::ThicknessMicroproblem::condensedCorrectionTangents< 9 >( weights, evaluations, B );

  // mean stress and condensed mean tangent: A_eff_ijmn = sum( weight * ( A_ijmn + A_ijkl N_l dw_k,mn ) )
  Tensor33d   PMean( 0.0 );
  Tensor3333d AMean( 0.0 );
  for ( int q = 0; q < nq; q++ ) {
    PMean += weights[q] * evaluations[q].P;
    AMean += weights[q] * evaluations[q].A;
    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ )
        for ( int m = 0; m < 3; m++ )
          for ( int n = 0; n < 3; n++ ) {
            double sum = 0.;
            for ( int k = 0; k < 3; k++ )
              for ( int l = 0; l < 3; l++ )
                sum += evaluations[q].A( i, j, k, l ) * normal( l ) * dw[q]( k, 3 * m + n );
            AMean( i, j, m, n ) += weights[q] * sum;
          }
  }

  const ReducedResponse reduced = reduceToInterface( PMean, AMean, normal, fibre, h );

  state.force         = reduced.force;
  state.surfaceStress = reduced.surfaceStress;

  tangents.Q_ik   = reduced.Q_ik;
  tangents.H_ikR  = reduced.H_ikR;
  tangents.H_iRk  = reduced.H_iRk;
  tangents.A_iRkP = reduced.A_iRkP;
}

void MarmotGradedInterfaceMaterialFiniteStrain::initializeYourself( double* stateVars, int )
{
  for ( auto& point : points ) {
    std::fill( stateVars + point.offset, stateVars + point.offset + 3, 0.0 );
    point.material->initializeYourself( stateVars + point.offset + 3, static_cast< int >( point.nBaseStateVars ) );
  }
}

double MarmotGradedInterfaceMaterialFiniteStrain::getDensity() const
{
  double density = 0.;
  for ( const auto& point : points )
    density += point.weight * point.material->getDensity( nullptr );
  return density;
}
