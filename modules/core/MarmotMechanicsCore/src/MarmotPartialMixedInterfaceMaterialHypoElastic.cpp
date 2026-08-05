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

#include "Marmot/MarmotPartialMixedInterfaceMaterialHypoElastic.h"

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Marmot;

namespace {

  using Vector3d          = Eigen::Matrix< double, 3, 1 >;
  using Vector9d          = Eigen::Matrix< double, 9, 1 >;
  using Matrix3dRowMajor  = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor  = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;
  using Matrix3x9RowMajor = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;

  constexpr int flatIndex( int i, int j )
  {
    return 3 * i + j;
  }

  constexpr int    maxLocalNewtonIterations = 40;
  constexpr double localNewtonTolerance     = 1.0e-10;

  Matrix9dRowMajor fullGradientTangent( const Matrix6d& tangentVoigt )
  {
    Matrix9dRowMajor tangentFull = Matrix9dRowMajor::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        Matrix3dRowMajor dGradient = Matrix3dRowMajor::Zero();
        dGradient( k, l )          = 1.0;

        const Matrix3dRowMajor dStrain      = 0.5 * ( dGradient + dGradient.transpose() );
        const Vector6d         dStrainVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt( dStrain );
        const Vector6d         dStressVoigt = tangentVoigt * dStrainVoigt;
        const Eigen::Matrix3d  dStress      = ContinuumMechanics::VoigtNotation::voigtToStress( dStressVoigt );

        const int column = flatIndex( k, l );
        for ( int i = 0; i < 3; ++i )
          for ( int j = 0; j < 3; ++j )
            tangentFull( flatIndex( i, j ), column ) = dStress( i, j );
      }
    }

    return tangentFull;
  }

  /** Maps a 3-vector v to vec(v (x) n). */
  Matrix9x3RowMajor jumpToGradient( const Vector3d& normal )
  {
    Matrix9x3RowMajor map = Matrix9x3RowMajor::Zero();
    for ( int k = 0; k < 3; ++k )
      for ( int l = 0; l < 3; ++l )
        map( flatIndex( k, l ), k ) = normal( l );
    return map;
  }

  /** Maps vec(sigma) to sigma n. */
  Matrix3x9RowMajor stressToForce( const Vector3d& normal )
  {
    Matrix3x9RowMajor map = Matrix3x9RowMajor::Zero();
    for ( int i = 0; i < 3; ++i )
      for ( int b = 0; b < 3; ++b )
        map( i, flatIndex( i, b ) ) = normal( b );
    return map;
  }

  struct SideTrial {
    Matrix3dRowMajor      stress;
    Matrix9dRowMajor      CFull;
    std::vector< double > trialStateVars;
  };

  SideTrial evaluateSideTrial( MarmotMaterialHypoElastic&                 material,
                               const double*                              oldStateVars,
                               int                                        nStateVars,
                               const Matrix3dRowMajor&                    stressCurrent,
                               const Vector6d&                            strainIncrementVoigt,
                               const MarmotMaterialHypoElastic::timeInfo& timeInfo )
  {
    SideTrial trial;
    trial.trialStateVars.assign( oldStateVars, oldStateVars + nStateVars );

    const Eigen::Matrix3d stressCurrentSym( 0.5 * ( stressCurrent + stressCurrent.transpose() ) );
    const Vector6d        stressVoigt = ContinuumMechanics::VoigtNotation::stressToVoigt( stressCurrentSym );

    Matrix6d                           tangentVoigt = Matrix6d::Zero();
    MarmotMaterialHypoElastic::state3D baseState{ stressVoigt, 0.0, 0.0, trial.trialStateVars.data() };

    material.computeStress( baseState, tangentVoigt, strainIncrementVoigt, timeInfo );

    trial.stress = ContinuumMechanics::VoigtNotation::voigtToStress( baseState.stress );
    trial.CFull  = fullGradientTangent( tangentVoigt );
    return trial;
  }

} // namespace

MarmotPartialMixedInterfaceMaterialHypoElastic::MarmotPartialMixedInterfaceMaterialHypoElastic(
  const std::string& materialName,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : materialProperties( matProperties_ ), nMaterialProperties( nMaterialProperties_ ), materialNumber( materialNumber_ )
{
  if ( nMaterialProperties < 3 )
    throw std::invalid_argument(
      "MarmotPartialMixedInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h." );

  h = materialProperties[2];
  if ( h <= 0.0 )
    throw std::invalid_argument( "MarmotPartialMixedInterfaceMaterialHypoElastic requires h > 0." );

  baseMaterialProperties.reserve( nMaterialProperties - 1 );
  baseMaterialProperties.push_back( materialProperties[0] );
  baseMaterialProperties.push_back( materialProperties[1] );
  baseMaterialProperties.insert( baseMaterialProperties.end(),
                                 materialProperties + 3,
                                 materialProperties + nMaterialProperties );

  topMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     baseMaterialProperties.data(),
                                                                     static_cast< int >(
                                                                       baseMaterialProperties.size() ),
                                                                     materialNumber ) );
  bottomMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
    MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                     baseMaterialProperties.data(),
                                                                     static_cast< int >(
                                                                       baseMaterialProperties.size() ),
                                                                     materialNumber ) );

  stateLayout.add( "topMaterialStateVars", topMaterial->getNumberOfRequiredStateVars() );
  stateLayout.add( "bottomMaterialStateVars", bottomMaterial->getNumberOfRequiredStateVars() );
  stateLayout.add( "normalGradientJump", 3 );
  stateLayout.add( "stressPlus", 9 );
  stateLayout.add( "stressMinus", 9 );
  stateLayout.finalize();
}

void MarmotPartialMixedInterfaceMaterialHypoElastic::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  if ( topMaterial )
    topMaterial->setCharacteristicElementLength( length );
  if ( bottomMaterial )
    bottomMaterial->setCharacteristicElementLength( length );
}

void MarmotPartialMixedInterfaceMaterialHypoElastic::computeMixedKernel( double*              stateVars,
                                                                         const KernelInput&   input,
                                                                         KernelOutput&        output,
                                                                         const TimeIncrement& timeIncrement )
{
  if ( !topMaterial || !bottomMaterial )
    throw std::logic_error( "MarmotPartialMixedInterfaceMaterialHypoElastic has no base material." );

  Vector3d     normal     = Eigen::Map< const Vector3d >( input.normal );
  const double normalNorm = normal.norm();
  if ( normalNorm <= 1.0e-12 )
    throw std::invalid_argument( "MarmotPartialMixedInterfaceMaterialHypoElastic: interface normal is zero." );
  normal /= normalNorm;

  const double sideThickness = 0.5 * h;
  (void)sideThickness;

  const Eigen::Map< const Matrix3dRowMajor > dAbar( input.dAbar );
  const Eigen::Map< const Matrix3dRowMajor > dDeltaA( input.dDeltaA );
  const Vector3d                             dG = Eigen::Map< const Vector3d >( input.dG );

  const Matrix3dRowMajor dAPlus  = dAbar + 0.5 * dDeltaA;
  const Matrix3dRowMajor dAMinus = dAbar - 0.5 * dDeltaA;

  double*   topStateVars       = stateLayout.getPtr( stateVars, "topMaterialStateVars" );
  double*   bottomStateVars    = stateLayout.getPtr( stateVars, "bottomMaterialStateVars" );
  double*   normalGradientJump = stateLayout.getPtr( stateVars, "normalGradientJump" );
  double*   stressPlusPtr      = stateLayout.getPtr( stateVars, "stressPlus" );
  double*   stressMinusPtr     = stateLayout.getPtr( stateVars, "stressMinus" );
  const int nTopStateVars      = topMaterial->getNumberOfRequiredStateVars();
  const int nBottomStateVars   = bottomMaterial->getNumberOfRequiredStateVars();

  const Matrix3dRowMajor stressCurrentPlus  = Eigen::Map< const Matrix3dRowMajor >( stressPlusPtr );
  const Matrix3dRowMajor stressCurrentMinus = Eigen::Map< const Matrix3dRowMajor >( stressMinusPtr );

  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld, timeIncrement.dT };

  const Matrix3x9RowMajor Rt = stressToForce( normal );  // vec(sigma) -> sigma n
  const Matrix9x3RowMajor Bn = jumpToGradient( normal ); // v -> vec(v (x) n)

  struct LocalEvaluation {
    SideTrial plus;
    SideTrial minus;
    Vector3d  r;
    double    rNorm;
  };

  // gamma is the increment of the normal-gradient jump for THIS step (warm
  // started from the committed value). Reconstruction uses the increment
  // gradients and the increment normal gradients g^pm = dG +/- gamma/2.
  auto evaluateLocal = [&]( const Vector3d& gammaTrial ) -> LocalEvaluation {
    const Vector3d gPlus  = dG + 0.5 * gammaTrial;
    const Vector3d gMinus = dG - 0.5 * gammaTrial;

    const Matrix3dRowMajor GPlus  = dAPlus + gPlus * normal.transpose();
    const Matrix3dRowMajor GMinus = dAMinus + gMinus * normal.transpose();

    const Matrix3dRowMajor strainIncrementPlus  = 0.5 * ( GPlus + GPlus.transpose() );
    const Matrix3dRowMajor strainIncrementMinus = 0.5 * ( GMinus + GMinus.transpose() );

    const Vector6d strainIncrementPlusVoigt  = ContinuumMechanics::VoigtNotation::strainToVoigt( strainIncrementPlus );
    const Vector6d strainIncrementMinusVoigt = ContinuumMechanics::VoigtNotation::strainToVoigt( strainIncrementMinus );

    LocalEvaluation evaluation;
    evaluation.plus  = evaluateSideTrial( *topMaterial,
                                         topStateVars,
                                         nTopStateVars,
                                         stressCurrentPlus,
                                         strainIncrementPlusVoigt,
                                         timeInfo );
    evaluation.minus = evaluateSideTrial( *bottomMaterial,
                                          bottomStateVars,
                                          nBottomStateVars,
                                          stressCurrentMinus,
                                          strainIncrementMinusVoigt,
                                          timeInfo );

    const Vector3d tPlus  = evaluation.plus.stress * normal;
    const Vector3d tMinus = evaluation.minus.stress * normal;
    evaluation.r          = tPlus - tMinus;
    evaluation.rNorm      = evaluation.r.norm() / std::max( 1.0, std::max( tPlus.norm(), tMinus.norm() ) );
    return evaluation;
  };

  Vector3d gamma = Eigen::Map< const Vector3d >( normalGradientJump );

  LocalEvaluation current   = evaluateLocal( gamma );
  bool            converged = current.rNorm <= localNewtonTolerance;

  for ( int iteration = 0; iteration < maxLocalNewtonIterations && !converged; ++iteration ) {
    const Matrix3dRowMajor Qplus  = Rt * current.plus.CFull * Bn;
    const Matrix3dRowMajor Qminus = Rt * current.minus.CFull * Bn;
    const Matrix3dRowMajor Qavg   = 0.5 * ( Qplus + Qminus );

    const Vector3d dgamma = Qavg.fullPivLu().solve( -current.r );

    double          alpha = 1.0;
    LocalEvaluation candidate;
    bool            accepted = false;
    for ( int lineSearchIter = 0; lineSearchIter < 12; ++lineSearchIter ) {
      candidate = evaluateLocal( gamma + alpha * dgamma );
      if ( candidate.rNorm <= localNewtonTolerance || candidate.rNorm < current.rNorm ) {
        accepted = true;
        break;
      }
      alpha *= 0.5;
    }

    gamma += alpha * dgamma;
    current = accepted ? candidate : evaluateLocal( gamma );

    if ( current.rNorm <= localNewtonTolerance )
      converged = true;
  }

  if ( !converged )
    throw Marmot::StressUpdateFailed(
      "MarmotPartialMixedInterfaceMaterialHypoElastic: local traction-equilibrium solve did not converge." );

  const SideTrial& plusTrial  = current.plus;
  const SideTrial& minusTrial = current.minus;

  // Commit warm-start and both sides' state + stress.
  Eigen::Map< Vector3d >         gammaMap( normalGradientJump );
  Eigen::Map< Matrix3dRowMajor > stressPlusMap( stressPlusPtr );
  Eigen::Map< Matrix3dRowMajor > stressMinusMap( stressMinusPtr );
  gammaMap       = gamma;
  stressPlusMap  = plusTrial.stress;
  stressMinusMap = minusTrial.stress;
  std::copy( plusTrial.trialStateVars.begin(), plusTrial.trialStateVars.end(), topStateVars );
  std::copy( minusTrial.trialStateVars.begin(), minusTrial.trialStateVars.end(), bottomStateVars );

  // -----------------------------------------------------------------
  // Generalized stress outputs.
  //   s_Abar   = <sigma>       s_DeltaA = 1/4 [sigma]     t^mat = <sigma> n
  // -----------------------------------------------------------------
  const Eigen::Map< const Vector9d > vecSigmaPlus( plusTrial.stress.data() );
  const Eigen::Map< const Vector9d > vecSigmaMinus( minusTrial.stress.data() );

  const Vector9d sAbar   = 0.5 * ( vecSigmaPlus + vecSigmaMinus );
  const Vector9d sDeltaA = 0.25 * ( vecSigmaPlus - vecSigmaMinus );
  const Vector3d tMat    = Rt * ( 0.5 * ( vecSigmaPlus + vecSigmaMinus ) );

  Eigen::Map< Vector9d >( output.sAbar )   = sAbar;
  Eigen::Map< Vector9d >( output.sDeltaA ) = sDeltaA;
  Eigen::Map< Vector3d >( output.tMat )    = tMat;

  // -----------------------------------------------------------------
  // Consistent reduced tangent H^red = y,x - y,gamma <Q>^{-1} F_gamma,x
  // over x = (Abar[9], DeltaA[9], g[3]), y = (s_Abar[9], s_DeltaA[9], t[3]).
  // -----------------------------------------------------------------
  const Matrix9dRowMajor& Cp    = plusTrial.CFull;
  const Matrix9dRowMajor& Cm    = minusTrial.CFull;
  const Matrix9dRowMajor  Csum  = Cp + Cm;
  const Matrix9dRowMajor  Cdiff = Cp - Cm;

  const Matrix9x3RowMajor CsumBn  = Csum * Bn;  // 9x3
  const Matrix9x3RowMajor CdiffBn = Cdiff * Bn; // 9x3

  Eigen::Matrix< double, 21, 21 > yx  = Eigen::Matrix< double, 21, 21 >::Zero();
  Eigen::Matrix< double, 21, 3 >  yg  = Eigen::Matrix< double, 21, 3 >::Zero();
  Eigen::Matrix< double, 3, 21 >  Fgx = Eigen::Matrix< double, 3, 21 >::Zero();

  // y,x rows: s_Abar (0..8), s_DeltaA (9..17), t (18..20); cols: Abar (0..8), DeltaA (9..17), g (18..20)
  yx.block< 9, 9 >( 0, 0 )  = 0.5 * Csum;
  yx.block< 9, 9 >( 0, 9 )  = 0.25 * Cdiff;
  yx.block< 9, 3 >( 0, 18 ) = 0.5 * CsumBn;

  yx.block< 9, 9 >( 9, 0 )  = 0.25 * Cdiff;
  yx.block< 9, 9 >( 9, 9 )  = 0.125 * Csum;
  yx.block< 9, 3 >( 9, 18 ) = 0.25 * CdiffBn;

  yx.block< 3, 9 >( 18, 0 )  = 0.5 * ( Rt * Csum );
  yx.block< 3, 9 >( 18, 9 )  = 0.25 * ( Rt * Cdiff );
  yx.block< 3, 3 >( 18, 18 ) = 0.5 * ( Rt * CsumBn );

  // y,gamma
  yg.block< 9, 3 >( 0, 0 )  = 0.25 * CdiffBn;
  yg.block< 9, 3 >( 9, 0 )  = 0.125 * CsumBn;
  yg.block< 3, 3 >( 18, 0 ) = 0.25 * ( Rt * CdiffBn );

  // F_gamma,x
  Fgx.block< 3, 9 >( 0, 0 )  = Rt * Cdiff;
  Fgx.block< 3, 9 >( 0, 9 )  = 0.5 * ( Rt * Csum );
  Fgx.block< 3, 3 >( 0, 18 ) = Rt * CdiffBn;

  // F_gamma,gamma = <Q>
  const Matrix3dRowMajor Qavg = 0.5 * ( Rt * CsumBn );

  const Eigen::Matrix< double, 3, 21 >  QinvFgx = Qavg.fullPivLu().solve( Fgx );
  const Eigen::Matrix< double, 21, 21 > Hred    = yx - yg * QinvFgx;

  // Partition: a = (Abar, DeltaA) -> indices 0..17, g -> 18..20.
  Eigen::Map< Eigen::Matrix< double, 18, 18, Eigen::RowMajor > >( output.H_aa ) = Hred.block< 18, 18 >( 0, 0 );
  Eigen::Map< Eigen::Matrix< double, 18, 3, Eigen::RowMajor > >( output.H_ag )  = Hred.block< 18, 3 >( 0, 18 );
  Eigen::Map< Eigen::Matrix< double, 3, 18, Eigen::RowMajor > >( output.H_ga )  = Hred.block< 3, 18 >( 18, 0 );
  Eigen::Map< Eigen::Matrix< double, 3, 3, Eigen::RowMajor > >( output.H_gg )   = Hred.block< 3, 3 >( 18, 18 );
}

void MarmotPartialMixedInterfaceMaterialHypoElastic::initializeYourself( double* stateVars, int nStateVars )
{
  for ( int i = 0; i < nStateVars; ++i )
    stateVars[i] = 0.0;

  if ( !topMaterial || !bottomMaterial )
    return;

  topMaterial->initializeYourself( stateLayout.getPtr( stateVars, "topMaterialStateVars" ),
                                   topMaterial->getNumberOfRequiredStateVars() );
  bottomMaterial->initializeYourself( stateLayout.getPtr( stateVars, "bottomMaterialStateVars" ),
                                      bottomMaterial->getNumberOfRequiredStateVars() );
}

double MarmotPartialMixedInterfaceMaterialHypoElastic::getDensity()
{
  if ( !topMaterial )
    return -1;
  return topMaterial->getDensity( nullptr );
}
