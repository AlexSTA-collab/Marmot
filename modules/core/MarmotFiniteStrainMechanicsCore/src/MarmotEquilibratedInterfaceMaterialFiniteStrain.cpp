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

#include "Marmot/MarmotEquilibratedInterfaceMaterialFiniteStrain.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialFiniteStrainFromGradientEnhanced.h"
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

  /**
   * @brief Trial response of one side at a total displacement gradient, evaluated from a scratch copy of the
   * committed state.
   *
   * The persistent state variables are never touched: a local Newton iteration re-evaluates the same increment many
   * times, and a material that mutated its state in place would integrate it once per trial.
   */
  struct SideTrial {
    Matrix3dRowMajor      P;              ///< first Piola-Kirchhoff stress
    Matrix9dRowMajor      A;              ///< derivative of the first Piola-Kirchhoff stress with respect to F
    std::vector< double > trialStateVars; ///< state variables after the trial
  };

  /**
   * @brief Evaluate the bulk material of one side at a total displacement gradient.
   *
   * @param[in] material bulk material of the side
   * @param[in] committedStateVars committed state variables of the side
   * @param[in] nStateVars number of state variables of the side
   * @param[in] G total displacement gradient
   * @param[in] timeIncrement time at the beginning of the increment and time increment
   * @return the first Piola-Kirchhoff stress, its tangent and the updated state variables
   */
  SideTrial evaluateSide( MarmotMaterialFiniteStrain&                      material,
                          const double*                                    committedStateVars,
                          int                                              nStateVars,
                          const Matrix3dRowMajor&                          G,
                          const MarmotMaterialFiniteStrain::TimeIncrement& timeIncrement )
  {
    using namespace Fastor;

    SideTrial trial;
    trial.trialStateVars.assign( committedStateVars, committedStateVars + nStateVars );

    Tensor33d F( 0.0 );
    F.eye();
    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j )
        F( i, j ) += G( i, j );

    MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ),
                                                                    0.0,
                                                                    0.0,
                                                                    trial.trialStateVars.data() );
    MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 >    moduli;
    MarmotMaterialFiniteStrain::Deformation< 3 >          deformation{ F };

    material.computeStress( response, moduli, deformation, timeIncrement );

    const Tensor33d   FInvT = transpose( inverse( F ) );
    const Tensor33d   P     = einsum< ij, jk, to_ik >( response.tau, FInvT );
    const Tensor3333d A     = einsum< imkl, mj, to_ijkl >( moduli.dTau_dF, FInvT ) -
                          einsum< im, kj, ml, to_ijkl >( response.tau, FInvT, FInvT );

    for ( int i = 0; i < 3; ++i )
      for ( int j = 0; j < 3; ++j ) {
        trial.P( i, j ) = P( i, j );
        for ( int k = 0; k < 3; ++k )
          for ( int l = 0; l < 3; ++l )
            trial.A( flatIndex( i, j ), flatIndex( k, l ) ) = A( i, j, k, l );
      }

    return trial;
  }

} // namespace

MarmotEquilibratedInterfaceMaterialFiniteStrain::MarmotEquilibratedInterfaceMaterialFiniteStrain(
  const std::string& materialName,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : materialProperties( matProperties_ ),
    nMaterialProperties( nMaterialProperties_ ),
    h( checkedMaterialProperty( matProperties_, nMaterialProperties_, 2 ) ),
    materialNumber( materialNumber_ )
{
  if ( h <= 0.0 )
    throw std::invalid_argument( "MarmotEquilibratedInterfaceMaterialFiniteStrain requires a positive interface "
                                 "thickness h." );

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

  // A gradient-enhanced material is accepted too, with its nonlocal field held fixed. The adapter refuses a material
  // whose stress depends on the field, so the substitution is never silent.
  const auto createBulkMaterial = [&]() {
    return MarmotLibrary::createFiniteStrainMaterialAllowingGradientEnhanced( materialName,
                                                                              baseMaterialProperties.data(),
                                                                              static_cast< int >(
                                                                                baseMaterialProperties.size() ),
                                                                              materialNumber );
  };

  topMaterial    = createBulkMaterial();
  bottomMaterial = createBulkMaterial();

  stateLayout.add( "normalGradientJump", 3 );
  stateLayout.add( "firstPiolaPlus", 9 );
  stateLayout.add( "firstPiolaMinus", 9 );
  stateLayout.add( "topMaterialStateVars", topMaterial->getNumberOfRequiredStateVars() );
  stateLayout.add( "bottomMaterialStateVars", bottomMaterial->getNumberOfRequiredStateVars() );
  stateLayout.finalize();
}

void MarmotEquilibratedInterfaceMaterialFiniteStrain::setCharacteristicElementLength( double length )
{
  topMaterial->setCharacteristicElementLength( length );
  bottomMaterial->setCharacteristicElementLength( length );
}

void MarmotEquilibratedInterfaceMaterialFiniteStrain::computeStress( State&               state,
                                                                     Tangents&            tangents,
                                                                     const Deformation&   deformation,
                                                                     const TimeIncrement& timeIncrement )
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

  double* normalGradientJump = stateLayout.getPtr( state.stateVars, "normalGradientJump" );
  double* firstPiolaPlus     = stateLayout.getPtr( state.stateVars, "firstPiolaPlus" );
  double* firstPiolaMinus    = stateLayout.getPtr( state.stateVars, "firstPiolaMinus" );
  double* topStateVars       = stateLayout.getPtr( state.stateVars, "topMaterialStateVars" );
  double* bottomStateVars    = stateLayout.getPtr( state.stateVars, "bottomMaterialStateVars" );

  const int nTopStateVars    = topMaterial->getNumberOfRequiredStateVars();
  const int nBottomStateVars = bottomMaterial->getNumberOfRequiredStateVars();

  const MarmotMaterialFiniteStrain::TimeIncrement materialTime{ timeIncrement.timeOld, timeIncrement.dT };

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
    evaluation.plus  = evaluateSide( *topMaterial, topStateVars, nTopStateVars, GPlus, materialTime );
    evaluation.minus = evaluateSide( *bottomMaterial, bottomStateVars, nBottomStateVars, GMinus, materialTime );

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

  const Matrix21dRowMajor Kcond = px - pz * Qavg.fullPivLu().solve( rx );

  const Eigen::Map< const Vector9d > vecPPlus( plusTrial.P.data() );
  const Eigen::Map< const Vector9d > vecPMinus( minusTrial.P.data() );

  const Vector21d p = sideThickness * ( BxPlus.transpose() * vecPPlus + BxMinus.transpose() * vecPMinus );

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

void MarmotEquilibratedInterfaceMaterialFiniteStrain::initializeYourself( double* stateVars, int nStateVars )
{
  for ( int i = 0; i < nStateVars; ++i )
    stateVars[i] = 0.0;

  topMaterial->initializeYourself( stateLayout.getPtr( stateVars, "topMaterialStateVars" ),
                                   topMaterial->getNumberOfRequiredStateVars() );
  bottomMaterial->initializeYourself( stateLayout.getPtr( stateVars, "bottomMaterialStateVars" ),
                                      bottomMaterial->getNumberOfRequiredStateVars() );
}

double MarmotEquilibratedInterfaceMaterialFiniteStrain::getDensity() const
{
  return topMaterial->getDensity( nullptr );
}
