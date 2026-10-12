#include "Marmot/MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceFiniteStrainHelperFunctions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialGradientEnhancedFiniteStrainFactory.h"
#include "Marmot/MarmotThicknessMicroproblem.h"
#include "Marmot/MarmotThicknessProfile.h"
#include "Marmot/MarmotUtils.h"

#include "Fastor/Fastor.h"

#include <algorithm>
#include <stdexcept>
#include <string>

using namespace Marmot;
using namespace Marmot::FastorIndices;
using namespace Marmot::FastorStandardTensors;

namespace {

  /// Response of a station to the deformation gradient of a trial value of w.
  struct StationEvaluation {
    Eigen::Vector3d       w = Eigen::Vector3d::Zero(); ///< correction w of the normal derivative
    Tensor33d             P;                           ///< first Piola-Kirchhoff stress
    Tensor3333d           A;                           ///< derivative of the stress with respect to F
    Tensor33d             dP_dN;                       ///< derivative of the stress with respect to the nonlocal field
    Tensor33d             dL_dF;                       ///< derivative of the local field with respect to F
    double                dL_dN;                       ///< derivative of the local field with respect to the field
    double                L;                           ///< local driving field
    double                radius;                      ///< nonlocal radius
    Eigen::Matrix3d       Qinv;                        ///< inverse of the acoustic tensor A[N, N]
    Eigen::Vector3d       a;                           ///< traction P N
    std::vector< double > state;                       ///< correction and state variables after the update
  };

} // namespace

MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain::MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain(
  const std::string& materialName,
  const double*      matProperties,
  int                nMaterialProperties,
  int                materialNumber_ )
  : materialNumber( materialNumber_ ), h( checkedMaterialProperty( matProperties, nMaterialProperties, 0 ) )
{
  if ( static_cast< int >( checkedMaterialProperty( matProperties, nMaterialProperties, 1 ) ) != 2 )
    throw std::invalid_argument( "MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain requires the mode 2 of "
                                 "the thickness profile: the nonlocal fields live at the stations" );

  auto profile = Marmot::Materials::ThicknessProfile::fromMaterialProperties( matProperties, nMaterialProperties );

  std::size_t offset = 0;
  for ( auto& profilePoint : profile.points ) {
    Station station;
    station.weight     = profilePoint.weight;
    station.position   = profilePoint.position;
    station.properties = std::move( profilePoint.properties );

    station.material = std::unique_ptr< MarmotMaterialGradientEnhancedFiniteStrain >(
      MarmotLibrary::MarmotMaterialGradientEnhancedFiniteStrainFactory::createMaterial( materialName,
                                                                                        station.properties.data(),
                                                                                        static_cast< int >(
                                                                                          station.properties.size() ),
                                                                                        materialNumber ) );
    if ( !station.material )
      throw std::invalid_argument(
        "Unknown base material for MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain: " + materialName );

    station.nBaseStateVars = station.material->getNumberOfRequiredStateVars();
    station.offset         = offset;
    offset += 3 + station.nBaseStateVars;
    points.push_back( std::move( station ) );
  }

  for ( std::size_t q = 0; q < points.size(); q++ ) {
    stateLayout.add( "correction" + std::to_string( q ), 3 );
    stateLayout.add( "baseMaterialStateVars" + std::to_string( q ), points[q].nBaseStateVars );
  }
  stateLayout.finalize();
}

void MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain::setCharacteristicElementLength( double length )
{
  for ( auto& station : points )
    station.material->setCharacteristicElementLength( length );
}

void MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain::computeStress( State&               state,
                                                                               Tangents&            tangents,
                                                                               const Deformation&   deformation,
                                                                               const TimeIncrement& timeIncrement )
{
  using namespace Marmot::Materials::InterfaceFiniteStrainHelperFunctions;

  const int nq = static_cast< int >( points.size() );

  if ( deformation.nonlocalFields.size() != nq )
    throw std::invalid_argument( "MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain: " + std::to_string( nq ) +
                                 " nonlocal fields are required, but " +
                                 std::to_string( deformation.nonlocalFields.size() ) + " are given" );

  const Tensor3d& normal = deformation.normal;

  const NormalFibre fibre = normalFibre( deformation.U, deformation.surfaceGradient, normal, deformation.connector );

  const MarmotMaterialGradientEnhancedFiniteStrain::TimeIncrement materialTimeIncrement{ timeIncrement.timeOld,
                                                                                         timeIncrement.dT };

  // Evaluate the station q for the trial value w of the correction, starting from the state of the previous
  // increment.
  const auto evaluate = [&]( int q, const Eigen::Vector3d& w ) {
    const auto& station = points[q];

    StationEvaluation evaluation;
    evaluation.w = w;
    evaluation.state.assign( state.stateVars + station.offset,
                             state.stateVars + station.offset + 3 + station.nBaseStateVars );

    Tensor3d wTensor;
    for ( int i = 0; i < 3; i++ )
      wTensor( i ) = w( i );

    const Tensor33d F = fibre.F + einsum< i, j, to_ij >( wTensor, normal );

    MarmotMaterialGradientEnhancedFiniteStrain::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ),
                                                                                    0.0,
                                                                                    0.0,
                                                                                    0.0,
                                                                                    0.0,
                                                                                    evaluation.state.data() + 3 );
    MarmotMaterialGradientEnhancedFiniteStrain::AlgorithmicModuli< 3 >    materialTangents;
    MarmotMaterialGradientEnhancedFiniteStrain::Deformation< 3 >          materialDeformation{ F,
                                                                                      deformation.nonlocalFields( q ) };

    station.material->computeStress( response, materialTangents, materialDeformation, materialTimeIncrement );

    firstPiolaKirchhoff( response.tau, materialTangents.dTau_dF, F, evaluation.P, evaluation.A );

    const Tensor33d FInvT = transpose( inverse( F ) );

    evaluation.dP_dN  = einsum< ij, jk, to_ik >( materialTangents.dTau_dN, FInvT );
    evaluation.dL_dF  = materialTangents.dL_dF;
    evaluation.dL_dN  = materialTangents.dL_dN;
    evaluation.L      = response.L;
    evaluation.radius = response.nonLocalRadius;

    const Tensor3d  traction = einsum< ij, j, to_i >( evaluation.P, normal );
    const Tensor33d ANN      = einsum< ijkl, j, l, to_ik >( evaluation.A, normal, normal );

    evaluation.a    = Eigen::Vector3d( traction( 0 ), traction( 1 ), traction( 2 ) );
    evaluation.Qinv = Eigen::Map< const Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( ANN.data() ).inverse();

    if ( !evaluation.Qinv.allFinite() )
      throw StressUpdateFailed( MakeString() << __PRETTY_FUNCTION__
                                             << ": the acoustic tensor of a station through the thickness is "
                                                "singular" );

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
  const double tolerance = 1e-11 * std::sqrt( Fastor::inner( strain, strain ) ) + 1e-14;

  const auto evaluations = Marmot::Materials::ThicknessMicroproblem::solve< StationEvaluation >( weights,
                                                                                                 evaluate,
                                                                                                 initialCorrections,
                                                                                                 tolerance );

  // commit the state of the converged evaluation
  for ( int q = 0; q < nq; q++ )
    std::copy( evaluations[q].state.begin(), evaluations[q].state.end(), state.stateVars + points[q].offset );

  // Derivative of the traction of every station with respect to the strain of the layer, which consists of the
  // deformation gradient F^0 (9 columns, B(i, (k, l)) = A_ijkl N_j) and the nonlocal fields (nq columns, only the
  // field of the station enters its traction).
  const int nColumns = 9 + nq;

  std::vector< Eigen::Matrix< double, 3, Eigen::Dynamic > >
    B( nq, Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, nColumns ) );
  for ( int q = 0; q < nq; q++ ) {
    const Tensor333d Bq = einsum< ijkl, j, to_ikl >( evaluations[q].A, normal );
    B[q].leftCols( 9 )  = Eigen::Map< const Eigen::Matrix< double, 3, 9, Eigen::RowMajor > >( Bq.data() );

    const Tensor3d c = einsum< ij, j, to_i >( evaluations[q].dP_dN, normal );
    for ( int i = 0; i < 3; i++ )
      B[q]( i, 9 + q ) = c( i );
  }

  const auto dw = Marmot::Materials::ThicknessMicroproblem::condensedCorrectionTangents< Eigen::Dynamic >( weights,
                                                                                                           evaluations,
                                                                                                           B );

  // mean stress and its condensed derivatives with respect to F^0 and to the nonlocal fields
  Tensor33d                PMean( 0.0 );
  Tensor3333d              AMean( 0.0 );
  std::vector< Tensor33d > dPMean_dN( nq, Tensor33d( 0.0 ) );

  for ( int q = 0; q < nq; q++ ) {
    PMean += weights[q] * evaluations[q].P;
    AMean += weights[q] * evaluations[q].A;

    for ( int i = 0; i < 3; i++ )
      for ( int j = 0; j < 3; j++ ) {
        // the stress of the station through its own derivative with respect to the strain of the layer
        for ( int m = 0; m < 3; m++ )
          for ( int n = 0; n < 3; n++ ) {
            double sum = 0.;
            for ( int k = 0; k < 3; k++ )
              for ( int l = 0; l < 3; l++ )
                sum += evaluations[q].A( i, j, k, l ) * normal( l ) * dw[q]( k, 3 * m + n );
            AMean( i, j, m, n ) += weights[q] * sum;
          }

        for ( int f = 0; f < nq; f++ ) {
          double sum = 0.;
          for ( int k = 0; k < 3; k++ )
            for ( int l = 0; l < 3; l++ )
              sum += evaluations[q].A( i, j, k, l ) * normal( l ) * dw[q]( k, 9 + f );
          dPMean_dN[f]( i, j ) += weights[q] * sum;
        }
      }

    dPMean_dN[q] += weights[q] * evaluations[q].dP_dN;
  }

  const ReducedResponse reduced = reduceToInterface( PMean, AMean, normal, fibre, h );

  state.force         = reduced.force;
  state.surfaceStress = reduced.surfaceStress;

  tangents.Q_ik   = reduced.Q_ik;
  tangents.H_ikR  = reduced.H_ikR;
  tangents.H_iRk  = reduced.H_iRk;
  tangents.A_iRkP = reduced.A_iRkP;

  // coupling to the nonlocal fields
  tangents.dForce_dNonlocal         = Eigen::MatrixXd::Zero( 3, nq );
  tangents.dSurfaceStress_dNonlocal = Eigen::MatrixXd::Zero( 9, nq );
  tangents.dLocal_dJump             = Eigen::MatrixXd::Zero( nq, 3 );
  tangents.dLocal_dAverage          = Eigen::MatrixXd::Zero( nq, 9 );
  tangents.dLocal_dNonlocal         = Eigen::MatrixXd::Zero( nq, nq );

  state.localFields   = Eigen::VectorXd::Zero( nq );
  state.nonlocalRadii = Eigen::VectorXd::Zero( nq );

  for ( int f = 0; f < nq; f++ ) {
    Tensor3d  dForce( 0.0 );
    Tensor33d dSurface( 0.0 );
    reduceStress( dPMean_dN[f], normal, fibre, h, dForce, dSurface );

    for ( int i = 0; i < 3; i++ ) {
      tangents.dForce_dNonlocal( i, f ) = dForce( i );
      for ( int R = 0; R < 3; R++ )
        tangents.dSurfaceStress_dNonlocal( 3 * i + R, f ) = dSurface( i, R );
    }
  }

  for ( int q = 0; q < nq; q++ ) {
    state.localFields( q )   = evaluations[q].L;
    state.nonlocalRadii( q ) = evaluations[q].radius;

    // derivative of the local field of the station with respect to F^0 and to the fields
    Tensor33d G( evaluations[q].dL_dF );
    for ( int m = 0; m < 3; m++ )
      for ( int n = 0; n < 3; n++ ) {
        double sum = 0.;
        for ( int k = 0; k < 3; k++ )
          for ( int l = 0; l < 3; l++ )
            sum += evaluations[q].dL_dF( k, l ) * normal( l ) * dw[q]( k, 3 * m + n );
        G( m, n ) += sum;
      }

    Tensor3d  dJump( 0.0 );
    Tensor33d dAverage( 0.0 );
    reduceDerivative( G, normal, fibre, dJump, dAverage );

    for ( int i = 0; i < 3; i++ ) {
      tangents.dLocal_dJump( q, i ) = dJump( i );
      for ( int R = 0; R < 3; R++ )
        tangents.dLocal_dAverage( q, 3 * i + R ) = dAverage( i, R );
    }

    for ( int f = 0; f < nq; f++ ) {
      double sum = ( f == q ) ? evaluations[q].dL_dN : 0.;
      for ( int k = 0; k < 3; k++ )
        for ( int l = 0; l < 3; l++ )
          sum += evaluations[q].dL_dF( k, l ) * normal( l ) * dw[q]( k, 9 + f );
      tangents.dLocal_dNonlocal( q, f ) = sum;
    }
  }
}

void MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain::initializeYourself( double* stateVars, int )
{
  for ( auto& station : points ) {
    std::fill( stateVars + station.offset, stateVars + station.offset + 3, 0.0 );
    station.material->initializeYourself( stateVars + station.offset + 3,
                                          static_cast< int >( station.nBaseStateVars ) );
  }
}

double MarmotGradedGradientEnhancedInterfaceMaterialFiniteStrain::getDensity( const double* stateVars ) const
{
  double density = 0.;
  for ( const auto& station : points )
    density += station.weight * station.material->getDensity( stateVars + station.offset + 3 );
  return density;
}
