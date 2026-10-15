#include "Marmot/MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain.h"

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceFiniteStrainHelperFunctions.h"
#include "Marmot/MarmotTwoSidedGradedLayer.h"
#include "Marmot/MarmotTwoSidedInterfaceGeometry.h"
#include "Marmot/MarmotUtils.h"

#include "Fastor/Fastor.h"
#include <Eigen/Dense>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot;
using namespace Marmot::FastorIndices;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::Materials::TwoSidedInterface;
using namespace Marmot::Materials::TwoSidedGradedLayer;

MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain::
  MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain( const std::string& materialName,
                                                                         const double*      matProperties,
                                                                         int                nMaterialProperties,
                                                                         int                materialNumber_ )
  : MarmotGradedFiniteStrainLayer( materialName, matProperties, nMaterialProperties, materialNumber_ )
{
  if ( static_cast< int >( checkedMaterialProperty( matProperties, nMaterialProperties, 1 ) ) != 2 )
    throw std::invalid_argument( "MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain requires the "
                                 "mode 2 of the thickness profile: the nonlocal fields live at the stations" );
}

void MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain::computeStress(
  State&               state,
  Tangents&            tangents,
  const Deformation&   deformation,
  const TimeIncrement& timeIncrement )
{
  const int nq = numberOfThicknessPoints();

  if ( deformation.nonlocalFields.size() != nq )
    throw std::invalid_argument(
      "MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain: " + std::to_string( nq ) +
      " nonlocal fields are required, but " + std::to_string( deformation.nonlocalFields.size() ) + " are given" );

  Eigen::Vector3d normal = Eigen::Map< const Eigen::Vector3d >( deformation.normal.data() );
  if ( normal.norm() <= 1.0e-12 )
    throw std::invalid_argument(
      "MarmotGradedEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain: interface normal is zero." );
  normal.normalize();

  const InterfaceGeometry geometry = evaluateInterfaceGeometry( normal,
                                                                Eigen::Map< const Eigen::Vector3d >(
                                                                  deformation.connector.data() ),
                                                                h );

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

  const MarmotMaterialGradientEnhancedFiniteStrain::TimeIncrement materialTimeIncrement{ timeIncrement.timeOld,
                                                                                         timeIncrement.dT };

  // shift of the nonlocal fields of the stations, for the derivative of the tangent of a point
  std::vector< double > nonlocalShift( nq, 0.0 );

  // keep the state variables of the points of the committed state for every trial
  const auto evaluateBulkPoint = [&]( int q, const Eigen::Matrix< double, 9, 1 >& Fvector ) {
    const auto& point = points[q];

    PointResponse< 9 > response;
    response.state.assign( state.stateVars + point.offset, state.stateVars + point.offset + 3 + point.nBaseStateVars );

    const Tensor33d F( Fvector.data() );

    MarmotMaterialGradientEnhancedFiniteStrain::ConstitutiveResponse< 3 > bulk( Tensor33d( 0.0 ),
                                                                                0.0,
                                                                                0.0,
                                                                                0.0,
                                                                                0.0,
                                                                                response.state.data() + 3 );
    MarmotMaterialGradientEnhancedFiniteStrain::AlgorithmicModuli< 3 >    moduli;
    MarmotMaterialGradientEnhancedFiniteStrain::Deformation< 3 >          materialDeformation{ F,
                                                                                      deformation.nonlocalFields( q ) +
                                                                                        nonlocalShift[q] };

    point.material->computeStress( bulk, moduli, materialDeformation, materialTimeIncrement );

    Tensor33d   P;
    Tensor3333d A;
    Marmot::Materials::InterfaceFiniteStrainHelperFunctions::firstPiolaKirchhoff( bulk.tau, moduli.dTau_dF, F, P, A );

    const Tensor33d FInvT = transpose( inverse( F ) );
    const Tensor33d dP_dN = einsum< ij, jk, to_ik >( moduli.dTau_dN, FInvT );

    response.stress = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( P.data() );
    response.D      = Eigen::Map< const Eigen::Matrix< double, 9, 9, Eigen::RowMajor > >( A.data() );

    response.dStress_dNonlocal = Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( dP_dN.data() );
    response.dLocal_dStrain    = Eigen::Map< const Eigen::Matrix< double, 1, 9 > >( moduli.dL_dF.data() );
    response.dLocal_dNonlocal  = moduli.dL_dN;
    response.local             = bulk.L;
    response.radius            = bulk.nonLocalRadius;
    return response;
  };

  // the corrections of the last converged state start the iteration
  std::vector< Eigen::Vector3d > initialCorrections( nq );
  for ( int q = 0; q < nq; q++ )
    initialCorrections[q] = Eigen::Map< const Eigen::Vector3d >( state.stateVars + points[q].offset );

  const double tolerance = 1e-11 * x.segment< 18 >( 3 ).cwiseAbs().maxCoeff() + 1e-14;

  LayerSolver< 9 > solver( layer, x );
  solver.solve( initialCorrections, tolerance, evaluateBulkPoint );

  const auto& solved = solver.points();

  // derivative of the tangent of a point with respect to the nonlocal field of its station, from the committed
  // state: the response of the microproblem depends on it through the tangent
  constexpr double                             fieldStep = 1e-7;
  std::vector< Eigen::Matrix< double, 9, 9 > > dTangent_dNonlocal( nq, Eigen::Matrix< double, 9, 9 >::Zero() );
  if ( tangentialGradient ) {
    for ( int q = 0; q < nq; q++ ) {
      const Eigen::Matrix< double, 9, 1 > Fvector = layer.strainOffset + solved[q].gradient;
      Eigen::Matrix< double, 9, 9 >       sides[2];
      for ( int side = 0; side < 2; side++ ) {
        nonlocalShift[q] = ( side == 0 ? -fieldStep : fieldStep );
        sides[side]      = evaluateBulkPoint( q, Fvector ).D;
      }
      nonlocalShift[q]      = 0.0;
      dTangent_dNonlocal[q] = ( sides[1] - sides[0] ) / ( 2. * fieldStep );
    }
  }

  // commit the corrections and the state of the converged evaluations
  for ( int q = 0; q < nq; q++ ) {
    std::vector< double > pointState = solved[q].state;
    std::copy( solved[q].w.data(), solved[q].w.data() + 3, pointState.begin() );
    std::copy( pointState.begin(), pointState.end(), state.stateVars + points[q].offset );
  }

  // the traction of a station depends on the nonlocal field of the station only
  Eigen::Matrix< double, 3, 9 > Rt = Eigen::Matrix< double, 3, 9 >::Zero();
  for ( int i = 0; i < 3; i++ )
    for ( int b = 0; b < 3; b++ )
      Rt( i, 3 * i + b ) = normal( b );

  std::vector< Eigen::MatrixXd > extraTraction( nq, Eigen::MatrixXd::Zero( 3, nq ) );
  for ( int q = 0; q < nq; q++ )
    extraTraction[q].col( q ) = Rt * solved[q].dStress_dNonlocal;

  const auto dw = solver.condensedCorrections( extraTraction );

  const Vector21d p = solver.stress( dw );
  const Matrix21d K = solver.tangent( dw );

  // derivative of the generalized stress with respect to the nonlocal fields
  const Eigen::Matrix< double, 9, 3 >& N9 = solver.correctionGradientOperator();

  tangents.dStress_dNonlocal = Eigen::MatrixXd::Zero( 21, nq );
  tangents.dLocal_dStrain    = Eigen::MatrixXd::Zero( nq, 21 );
  tangents.dLocal_dNonlocal  = Eigen::MatrixXd::Zero( nq, nq );

  // the test operator depends on the nonlocal field through the tangent of the point, if there is a closure
  for ( int f = 0; f < nq; f++ )
    tangents.dStress_dNonlocal.col(
      f ) += solver.stressDerivativeThroughMicroResponse( f, dTangent_dNonlocal[f], fieldStep );

  state.localFields   = Eigen::VectorXd::Zero( nq );
  state.nonlocalRadii = Eigen::VectorXd::Zero( nq );

  for ( int q = 0; q < nq; q++ ) {
    const double scale = h * layer.weights[q];

    for ( int f = 0; f < nq; f++ ) {
      Eigen::Matrix< double, 9, 1 > dStress = solved[q].D * N9 * dw[q].col( 21 + f );
      if ( f == q )
        dStress += solved[q].dStress_dNonlocal;
      tangents.dStress_dNonlocal.col( f ) += scale * solver.testGradientOperator( q, dw ).transpose() * dStress;

      tangents.dLocal_dNonlocal( q, f ) = ( f == q ? solved[q].dLocal_dNonlocal : 0. ) +
                                          ( solved[q].dLocal_dStrain * N9 * dw[q].col( 21 + f ) )( 0, 0 );
    }

    tangents.dLocal_dStrain.row( q ) = solved[q].dLocal_dStrain *
                                       ( solver.strainOperator( q ) + N9 * dw[q].leftCols( 21 ) );

    state.localFields( q )   = solved[q].local;
    state.nonlocalRadii( q ) = solved[q].radius;
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
