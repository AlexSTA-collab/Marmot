/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Alexandros Stathas alexandros.stathas@boku.ac.at
 * Matthias Neuner matthias.neuner@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#include "Marmot/MarmotEquilibratedInterfaceMaterialFiniteStrainBase.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotUtils.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace Marmot;
using namespace Marmot::FastorStandardTensors;
using namespace Marmot::FastorIndices;

namespace {

  using Vector3d           = Eigen::Matrix< double, 3, 1 >;
  using Vector9d           = Eigen::Matrix< double, 9, 1 >;
  using Vector21d          = Eigen::Matrix< double, 21, 1 >;
  using Matrix3dRowMajor   = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor   = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor  = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;
  using Matrix3x9RowMajor  = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;
  using Matrix9x21RowMajor = Eigen::Matrix< double, 9, 21, Eigen::RowMajor >;
  using Matrix21x3RowMajor = Eigen::Matrix< double, 21, 3, Eigen::RowMajor >;
  using Matrix3x21RowMajor = Eigen::Matrix< double, 3, 21, Eigen::RowMajor >;
  using Matrix21dRowMajor  = Eigen::Matrix< double, 21, 21, Eigen::RowMajor >;

  constexpr int flatIndex( int i, int j )
  {
    return 3 * i + j;
  }

  /**
   * @brief Operator that maps a vector g to vec( g (x) n ) / normalSeparation.
   *
   * @param[in] normal normal of the interface
   * @param[in] normalSeparation separation of the faces along the normal
   * @return the 9x3 operator
   */
  Matrix9x3RowMajor jumpToGradient( const Vector3d& normal, double normalSeparation )
  {
    Matrix9x3RowMajor map = Matrix9x3RowMajor::Zero();

    for ( int k = 0; k < 3; ++k )
      for ( int l = 0; l < 3; ++l )
        map( flatIndex( k, l ), k ) = normal( l ) / normalSeparation;

    return map;
  }

  /**
   * @brief Operator that maps vec( A ) to vec( A ) - vec( ( A d_tau ) (x) n ) / normalSeparation.
   *
   * @param[in] normal normal of the interface
   * @param[in] tangentialSeparation tangential part of the connector between the faces
   * @param[in] normalSeparation separation of the faces along the normal
   * @return the 9x9 operator
   */
  Matrix9dRowMajor surfaceToGradient( const Vector3d& normal,
                                      const Vector3d& tangentialSeparation,
                                      double          normalSeparation )
  {
    Matrix9dRowMajor map = Matrix9dRowMajor::Zero();

    for ( int k = 0; k < 3; ++k )
      for ( int l = 0; l < 3; ++l )
        for ( int b = 0; b < 3; ++b )
          map( flatIndex( k, l ), flatIndex( k, b ) ) = ( l == b ? 1.0 : 0.0 ) -
                                                        tangentialSeparation( b ) * normal( l ) / normalSeparation;

    return map;
  }

  /**
   * @brief Operator that maps vec( P ) to the traction P n.
   *
   * @param[in] normal normal of the interface
   * @return the 3x9 operator
   */
  Matrix3x9RowMajor stressToTraction( const Vector3d& normal )
  {
    Matrix3x9RowMajor map = Matrix3x9RowMajor::Zero();

    for ( int i = 0; i < 3; ++i )
      for ( int b = 0; b < 3; ++b )
        map( i, flatIndex( i, b ) ) = normal( b );

    return map;
  }

  /// Split of the connector between the faces into its normal and tangential part.
  struct InterfaceGeometry {
    double   normalSeparation;     ///< component of the connector along the normal
    Vector3d tangentialSeparation; ///< component of the connector in the tangential plane
  };

  /**
   * @brief Split the connector between the faces into its normal and tangential part.
   *
   * @param[in] normal unit normal of the interface
   * @param[in] connector connector between the paired faces
   * @return the normal separation and the tangential separation
   */
  InterfaceGeometry evaluateInterfaceGeometry( const Vector3d& normal, const Vector3d& connector )
  {
    InterfaceGeometry geometry;
    geometry.normalSeparation     = connector.dot( normal );
    geometry.tangentialSeparation = connector - geometry.normalSeparation * normal;
    return geometry;
  }

} // namespace

MarmotEquilibratedInterfaceMaterialFiniteStrainBase::MarmotEquilibratedInterfaceMaterialFiniteStrainBase(
  const double* matProperties_,
  int           nMaterialProperties_,
  int           materialNumber_,
  const char*   className )
  : materialProperties( matProperties_ ),
    nMaterialProperties( nMaterialProperties_ ),
    h( checkedMaterialProperty( matProperties_, nMaterialProperties_, 2 ) ),
    materialNumber( materialNumber_ )
{
  if ( h <= 0.0 )
    throw std::invalid_argument( std::string( className ) + " requires a positive interface thickness h." );

  // The card is [ p0, p1, h, ... ], the layout of MarmotEquilibratedXInterfaceMaterialHypoElastic, so that a study
  // deck moves to the finite-strain element by changing the element type alone. p0 and p1 are the first two
  // properties of the bulk law, h is the thickness of the layer normal to the interface and is removed from the list
  // handed to the bulk law.
  baseMaterialProperties.reserve( nMaterialProperties - 1 );
  baseMaterialProperties.push_back( materialProperties[0] );
  baseMaterialProperties.push_back( materialProperties[1] );
  baseMaterialProperties.insert( baseMaterialProperties.end(),
                                 materialProperties + 3,
                                 materialProperties + nMaterialProperties );
}

void MarmotEquilibratedInterfaceMaterialFiniteStrainBase::finalizeStateLayout( int nTopStateVars, int nBottomStateVars )
{
  stateLayout.add( "normalGradientJump", 3 );
  stateLayout.add( "firstPiolaPlus", 9 );
  stateLayout.add( "firstPiolaMinus", 9 );
  stateLayout.add( "topMaterialStateVars", nTopStateVars );
  stateLayout.add( "bottomMaterialStateVars", nBottomStateVars );
  stateLayout.finalize();
}

MarmotEquilibratedInterfaceMaterialFiniteStrainBase::SideTrial MarmotEquilibratedInterfaceMaterialFiniteStrainBase::
  firstPiolaFromKirchhoff( const Tensor33d& F, const Tensor33d& tau, const Tensor3333d& dTau_dF )
{
  using namespace Fastor;

  const Tensor33d   FInvT = transpose( inverse( F ) );
  const Tensor33d   P     = einsum< ij, jk, to_ik >( tau, FInvT );
  const Tensor3333d A     = einsum< imkl, mj, to_ijkl >( dTau_dF, FInvT ) -
                        einsum< im, kj, ml, to_ijkl >( tau, FInvT, FInvT );

  SideTrial trial;
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j ) {
      trial.P( i, j ) = P( i, j );
      for ( int k = 0; k < 3; ++k )
        for ( int l = 0; l < 3; ++l )
          trial.A( flatIndex( i, j ), flatIndex( k, l ) ) = A( i, j, k, l );
    }

  return trial;
}

MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Solution MarmotEquilibratedInterfaceMaterialFiniteStrainBase::
  solve( double*                                                            stateVars,
         const Deformation&                                                 deformation,
         const std::function< SideTrial( Side, const Matrix3dRowMajor& ) >& evaluateSide,
         bool                                                               withNonlocal ) const
{
  Vector3d     normal     = Eigen::Map< const Vector3d >( deformation.normal.data() );
  const double normalNorm = normal.norm();
  if ( normalNorm <= 1.0e-12 ) {
    throw std::invalid_argument( "MarmotEquilibratedInterfaceMaterialFiniteStrain: interface normal is zero." );
  }
  normal /= normalNorm;

  const Vector3d          connector = Eigen::Map< const Vector3d >( deformation.connector.data() );
  const InterfaceGeometry geometry  = evaluateInterfaceGeometry( normal, connector );

  const double    ell         = geometry.normalSeparation;
  const Vector3d& dTangential = geometry.tangentialSeparation;

  if ( ell <= 0.0 ) {
    throw std::invalid_argument(
      "MarmotEquilibratedInterfaceMaterialFiniteStrain: the paired faces must be separated along the interface "
      "normal, and the normal must point from the first side to the second. The connector between them measures "
      "the normal fibre the reduction differentiates along." );
  }

  const double sideThickness = 0.5 * h;

  const Eigen::Map< const Eigen::Matrix< double, 6, 1 > > U( deformation.U.data() );

  const Vector3d w = U.segment< 3 >( 0 ) - U.segment< 3 >( 3 );

  const Eigen::Map< const Matrix3dRowMajor > APlus( deformation.surfaceGradient.data() );
  const Eigen::Map< const Matrix3dRowMajor > AMinus( deformation.surfaceGradient.data() + 9 );
  const Matrix3dRowMajor                     ABar = 0.5 * ( APlus + AMinus );

  const Vector3d gBar = ( w - ABar * dTangential ) / ell;

  double* normalGradientJump = stateLayout.getPtr( stateVars, "normalGradientJump" );
  double* firstPiolaPlus     = stateLayout.getPtr( stateVars, "firstPiolaPlus" );
  double* firstPiolaMinus    = stateLayout.getPtr( stateVars, "firstPiolaMinus" );
  double* topStateVars       = stateLayout.getPtr( stateVars, "topMaterialStateVars" );
  double* bottomStateVars    = stateLayout.getPtr( stateVars, "bottomMaterialStateVars" );

  const Matrix3x9RowMajor Rt = stressToTraction( normal );
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

    LocalEvaluation evaluation;
    evaluation.plus  = evaluateSide( Side::Top, GPlus );
    evaluation.minus = evaluateSide( Side::Bottom, GMinus );

    const Vector3d tPlus  = evaluation.plus.P * normal;
    const Vector3d tMinus = evaluation.minus.P * normal;

    evaluation.r     = tPlus - tMinus;
    evaluation.rNorm = evaluation.r.norm() / std::max( 1.0, std::max( tPlus.norm(), tMinus.norm() ) );
    return evaluation;
  };

  // Local Newton for z, enforcing P+ N == P- N. Every trial restarts from the committed state.
  constexpr int    maxLocalIterations = 40;
  constexpr double localTolerance     = 1.0e-10;

  Vector3d z = Eigen::Map< const Vector3d >( normalGradientJump );

  LocalEvaluation current = evaluateLocal( z );

  Vector3d        bestZ    = z;
  LocalEvaluation bestEval = current;

  for ( int iteration = 0; iteration < maxLocalIterations && current.rNorm > localTolerance; ++iteration ) {
    const Matrix3dRowMajor Qplus  = Rt * current.plus.A * Bn;
    const Matrix3dRowMajor Qminus = Rt * current.minus.A * Bn;
    const Matrix3dRowMajor Qavg   = 0.5 * ( Qplus + Qminus );

    const Vector3d dz = Qavg.fullPivLu().solve( -current.r );

    // Backtracking: a full step can overshoot when a trial state crosses a yield-surface kink,
    // so only a step that reduces the residual is accepted. Without the guard the base material
    // can be asked for a state it cannot integrate and throws, which surfaces as a global
    // cutback -- the same failure the small-strain element documents.
    double          alpha = 1.0;
    LocalEvaluation candidate;
    bool            accepted = false;

    for ( int lineSearchIteration = 0; lineSearchIteration < 12; ++lineSearchIteration ) {
      try {
        candidate = evaluateLocal( z + alpha * dz );
      }
      catch ( const Marmot::StressUpdateFailed& ) {
        alpha *= 0.5;
        continue;
      }

      if ( candidate.rNorm <= localTolerance || candidate.rNorm < current.rNorm ) {
        accepted = true;
        break;
      }
      alpha *= 0.5;
    }

    if ( !accepted )
      break;

    z += alpha * dz;
    current = candidate;

    if ( current.rNorm < bestEval.rNorm ) {
      bestEval = current;
      bestZ    = z;
    }
  }

  z       = bestZ;
  current = bestEval;

  const SideTrial& plusTrial  = current.plus;
  const SideTrial& minusTrial = current.minus;

  // Commit: the converged internal variable warm-starts the next increment, and each side's
  // material state and first Piola stress are persisted.
  Eigen::Map< Vector3d >         normalGradientJumpMap( normalGradientJump );
  Eigen::Map< Matrix3dRowMajor > firstPiolaPlusMap( firstPiolaPlus );
  Eigen::Map< Matrix3dRowMajor > firstPiolaMinusMap( firstPiolaMinus );

  normalGradientJumpMap = z;
  firstPiolaPlusMap     = plusTrial.P;
  firstPiolaMinusMap    = minusTrial.P;
  std::copy( plusTrial.trialStateVars.begin(), plusTrial.trialStateVars.end(), topStateVars );
  std::copy( minusTrial.trialStateVars.begin(), minusTrial.trialStateVars.end(), bottomStateVars );

  // Purely geometric operators: dF/dx and dF/dz, constant because F is affine in both.
  const Matrix9x3RowMajor Bw                = jumpToGradient( normal, ell );
  const Matrix9dRowMajor  Asrf              = surfaceToGradient( normal, dTangential, ell );
  const Matrix9dRowMajor  I9                = Matrix9dRowMajor::Identity();
  const Matrix9dRowMajor  half_I_plus_Asrf  = 0.5 * ( I9 + Asrf );
  const Matrix9dRowMajor  half_Asrf_minus_I = 0.5 * ( Asrf - I9 );

  Matrix9x21RowMajor BxPlus  = Matrix9x21RowMajor::Zero();
  Matrix9x21RowMajor BxMinus = Matrix9x21RowMajor::Zero();

  BxPlus.block< 9, 3 >( 0, 0 )   = Bw;
  BxPlus.block< 9, 9 >( 0, 3 )   = half_I_plus_Asrf;
  BxPlus.block< 9, 9 >( 0, 12 )  = half_Asrf_minus_I;
  BxMinus.block< 9, 3 >( 0, 0 )  = Bw;
  BxMinus.block< 9, 9 >( 0, 3 )  = half_Asrf_minus_I;
  BxMinus.block< 9, 9 >( 0, 12 ) = half_I_plus_Asrf;

  const Matrix9x3RowMajor BzPlus  = 0.5 * Bn;
  const Matrix9x3RowMajor BzMinus = -0.5 * Bn;

  const Matrix3dRowMajor Qplus  = Rt * plusTrial.A * Bn;
  const Matrix3dRowMajor Qminus = Rt * minusTrial.A * Bn;
  const Matrix3dRowMajor Qavg   = 0.5 * ( Qplus + Qminus );

  const Matrix3x21RowMajor rx = Rt * ( plusTrial.A * BxPlus - minusTrial.A * BxMinus );

  const Matrix21dRowMajor  px = sideThickness * ( BxPlus.transpose() * plusTrial.A * BxPlus +
                                                 BxMinus.transpose() * minusTrial.A * BxMinus );
  const Matrix21x3RowMajor pz = sideThickness * ( BxPlus.transpose() * plusTrial.A * BzPlus +
                                                  BxMinus.transpose() * minusTrial.A * BzMinus );

  const auto QavgDecomposition = Qavg.fullPivLu();

  const Matrix21dRowMajor Kcond = px - pz * QavgDecomposition.solve( rx );

  const Eigen::Map< const Vector9d > vecPPlus( plusTrial.P.data() );
  const Eigen::Map< const Vector9d > vecPMinus( minusTrial.P.data() );

  const Vector21d p = sideThickness * ( BxPlus.transpose() * vecPPlus + BxMinus.transpose() * vecPMinus );

  Solution solution;
  solution.z       = z;
  solution.p       = p;
  solution.tangent = Kcond;

  if ( withNonlocal ) {
    // The nonlocal field enters both faces, so z depends on it too, and all quantities below share the same <Q>^-1.
    Vector9d vecdPPlus_dNonlocal, vecdPMinus_dNonlocal;
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j ) {
        vecdPPlus_dNonlocal( flatIndex( i, j ) )  = plusTrial.dP_dNonlocal( i, j );
        vecdPMinus_dNonlocal( flatIndex( i, j ) ) = minusTrial.dP_dNonlocal( i, j );
      }

    const Vector3d ra = Rt * ( vecdPPlus_dNonlocal - vecdPMinus_dNonlocal );

    const Vector21d pa = sideThickness *
                         ( BxPlus.transpose() * vecdPPlus_dNonlocal + BxMinus.transpose() * vecdPMinus_dNonlocal );

    // Lbar = ( L+ + L- ) / 2, the two faces carrying half the layer each
    const Eigen::Matrix< double, 1, 21 > Lx = 0.5 * ( plusTrial.dL_dF.transpose() * BxPlus +
                                                      minusTrial.dL_dF.transpose() * BxMinus );
    const Eigen::Matrix< double, 1, 3 >  Lz = 0.5 * ( plusTrial.dL_dF.transpose() * BzPlus +
                                                     minusTrial.dL_dF.transpose() * BzMinus );

    const Matrix3x21RowMajor QavgInverseRx = QavgDecomposition.solve( rx );
    const Vector3d           QavgInverseRa = QavgDecomposition.solve( ra );

    solution.dStress_dNonlocal     = pa - pz * QavgInverseRa;
    solution.dLocalField_dStrain   = Lx - Lz * QavgInverseRx;
    solution.dLocalField_dNonlocal = 0.5 * ( plusTrial.dL_dN + minusTrial.dL_dN ) - ( Lz * QavgInverseRa )( 0, 0 );
    solution.localField            = 0.5 * ( plusTrial.L + minusTrial.L );
    solution.nonlocalRadius        = 0.5 * ( plusTrial.nonlocalRadius + minusTrial.nonlocalRadius );
  }

  solution.plus  = plusTrial;
  solution.minus = minusTrial;
  return solution;
}

void MarmotEquilibratedInterfaceMaterialFiniteStrainBase::writeResponse( const Solution& solution,
                                                                         State&          state,
                                                                         Tangents&       tangents )
{
  using Matrix3x9RowMajor = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;

  const Vector21d& p     = solution.p;
  const Matrix21d& Kcond = solution.tangent;

  Eigen::Map< Vector3d >( state.generalizedForce.data() )   = p.segment< 3 >( 0 );
  Eigen::Map< Vector9d >( state.surfaceStressPlus.data() )  = p.segment< 9 >( 3 );
  Eigen::Map< Vector9d >( state.surfaceStressMinus.data() ) = p.segment< 9 >( 12 );

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
