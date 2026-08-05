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

#pragma once

// Reuses the LobattoRule<N> table and the Marmot::Detail::GaussLobattoInterface
// helper functions (fullGradientTangent, evaluateInterfaceGeometry,
// jumpToGradient, stressToForce, tangentialContraction, evaluateStationTrial,
// naming helpers) from the validated GLIQUAD4 material -- that file is
// included, not modified.
#include "Marmot/MarmotGaussLobattoInterfaceMaterialHypoElastic.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * B-bar variant of MarmotGaussLobattoInterfaceMaterialHypoElasticN.
 *
 * Identical formulation in every respect (kinematics, per-station Lobatto
 * integration, local traction-equilibrium Newton solve with the
 * warm-start/fresh-restart fallback, static condensation) EXCEPT for one
 * change: at every Lobatto station alpha, the tangential surface gradient
 * A^(alpha) is augmented by a caller-supplied per-station volumetric
 * correction,
 *
 *   A^(alpha) <- A^(alpha) + (traceCorrection[alpha] / 3) * I,
 *
 * before it enters the (otherwise unmodified) strain reconstruction
 * G^(alpha) = A^(alpha) + g^(alpha) (x) n. The correction is computed by
 * the element (GaussLobattoBBarInterfaceFiniteElement) as
 *
 *   traceCorrection[alpha] = averageTrace[alpha] - trace(A^(alpha)_raw),
 *
 * i.e. it replaces the pointwise volumetric part of the surface-gradient
 * strain by its element-average value (over the surface Gauss points,
 * SEPARATELY per Lobatto station) -- the classical B-bar projection,
 * applied independently at each through-thickness station.
 *
 * The correction is treated as a fixed external input for the purposes of
 * this call's own local Newton and static condensation -- i.e. the
 * analytic tangent here is the standard "frozen" B-bar tangent (the
 * correction is not differentiated with respect to q), matching the
 * classical B-bar strain-displacement projection used in solid mechanics:
 * the same B-bar operator is used for both the strain and the stiffness
 * assembly, and the cross-Gauss-point dependence of the average is not
 * carried through the linearization.
 */
template < int NStations >
class MarmotGaussLobattoBBarInterfaceMaterialHypoElasticN {

public:
  static constexpr int nStations = NStations;

protected:
  const double*                                                         materialProperties;
  const int                                                             nMaterialProperties;
  double                                                                h = 0.0;
  std::vector< double >                                                 baseMaterialProperties;
  std::array< std::unique_ptr< MarmotMaterialHypoElastic >, NStations > stationMaterials;

public:
  using TensorMap3d    = Marmot::FastorStandardTensors::TensorMap3d;
  using TensorMap33d   = Marmot::FastorStandardTensors::TensorMap33d;
  using TensorMap333d  = Marmot::FastorStandardTensors::TensorMap333d;
  using TensorMap3333d = Marmot::FastorStandardTensors::TensorMap3333d;
  using TensorMap6d    = Marmot::FastorStandardTensors::TensorMap6d;
  using TensorMap18d   = Marmot::FastorStandardTensors::TensorMap18d;

  const int materialNumber;

  MarmotGaussLobattoBBarInterfaceMaterialHypoElasticN( const std::string& materialName,
                                                       const double*      matProperties_,
                                                       int                nMaterialProperties_,
                                                       int                materialNumber_ );

  virtual ~MarmotGaussLobattoBBarInterfaceMaterialHypoElasticN() = default;

  MarmotStateLayoutDynamic stateLayout;

  double characteristicElementLength;

  void setCharacteristicElementLength( double length );

  struct State {
    TensorMap3d  generalizedForce;
    TensorMap33d surfaceStressPlus;
    TensorMap33d surfaceStressMinus;
    double*      stateVars;
  };

  struct Tangents {
    TensorMap33d   Q_ww;
    TensorMap333d  Q_wAp;
    TensorMap333d  Q_wAm;
    TensorMap333d  Q_Apw;
    TensorMap3333d Q_ApAp;
    TensorMap3333d Q_ApAm;
    TensorMap333d  Q_Amw;
    TensorMap3333d Q_AmAp;
    TensorMap3333d Q_AmAm;
  };

  struct Deformation {
    TensorMap6d  dU;
    TensorMap18d dSurfaceStrain;
    TensorMap3d  normal;
    TensorMap3d  separationVector;

    /** Per-station volumetric (trace) correction, NStations scalars:
     * traceCorrection[alpha] = averageTrace[alpha] - trace(A^(alpha)_raw).
     * Owned by the caller; must outlive the computeStress call. */
    const double* traceCorrection;

    Deformation( const double* dU_,
                 const double* dSurfaceStrain_,
                 const double* normal_,
                 const double* separationVector_,
                 const double* traceCorrection_ )
      : dU( const_cast< double* >( dU_ ) ),
        dSurfaceStrain( const_cast< double* >( dSurfaceStrain_ ) ),
        normal( const_cast< double* >( normal_ ) ),
        separationVector( const_cast< double* >( separationVector_ ) ),
        traceCorrection( traceCorrection_ )
    {
    }

    /** Backward-compatible constructor for coincident interface faces. */
    Deformation( const double* dU_,
                 const double* dSurfaceStrain_,
                 const double* normal_,
                 const double* traceCorrection_ )
      : Deformation( dU_, dSurfaceStrain_, normal_, zeroSeparationVector(), traceCorrection_ )
    {
    }

  private:
    static double* zeroSeparationVector()
    {
      static double zero[3] = { 0.0, 0.0, 0.0 };
      return zero;
    }
  };

  struct TimeIncrement {
    double timeOld;
    double dT;
  };

  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );

  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    return stateLayout.getStateView( stateVars, stateName );
  }

  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  virtual void initializeYourself( double* stateVars, int nStateVars );

  virtual double getDensity();
};

/** Preferred production rule: fixed five-point Gauss-Lobatto, B-bar variant. */
using MarmotGaussLobattoBBarInterfaceMaterialHypoElastic = MarmotGaussLobattoBBarInterfaceMaterialHypoElasticN< 5 >;

// ---------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------

template < int NStations >
MarmotGaussLobattoBBarInterfaceMaterialHypoElasticN< NStations >::MarmotGaussLobattoBBarInterfaceMaterialHypoElasticN(
  const std::string& materialName,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : materialProperties( matProperties_ ), nMaterialProperties( nMaterialProperties_ ), materialNumber( materialNumber_ )
{
  using namespace Marmot::Detail::GaussLobattoInterface;

  if ( nMaterialProperties < 3 ) {
    throw std::invalid_argument(
      "MarmotGaussLobattoBBarInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h." );
  }

  h = materialProperties[2];
  if ( h <= 0.0 ) {
    throw std::invalid_argument( "MarmotGaussLobattoBBarInterfaceMaterialHypoElastic requires h > 0." );
  }

  baseMaterialProperties.reserve( nMaterialProperties - 1 );
  baseMaterialProperties.push_back( materialProperties[0] );
  baseMaterialProperties.push_back( materialProperties[1] );
  baseMaterialProperties.insert( baseMaterialProperties.end(),
                                 materialProperties + 3,
                                 materialProperties + nMaterialProperties );

  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    stationMaterials[alpha] = std::unique_ptr< MarmotMaterialHypoElastic >(
      MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                       baseMaterialProperties.data(),
                                                                       static_cast< int >(
                                                                         baseMaterialProperties.size() ),
                                                                       materialNumber ) );
  }

  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    stateLayout.add( stationStateVarsName( alpha ), stationMaterials[alpha]->getNumberOfRequiredStateVars() );
    stateLayout.add( normalGradientName( alpha ), 3 );
    stateLayout.add( stationStressName( alpha ), 9 );
  }
  stateLayout.add( "commonTraction", 3 );
  stateLayout.finalize();
}

template < int NStations >
void MarmotGaussLobattoBBarInterfaceMaterialHypoElasticN< NStations >::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  for ( auto& material : stationMaterials ) {
    if ( material ) {
      material->setCharacteristicElementLength( length );
    }
  }
}

template < int NStations >
void MarmotGaussLobattoBBarInterfaceMaterialHypoElasticN< NStations >::computeStress(
  State&               state,
  Tangents&            tangents,
  const Deformation&   deformation,
  const TimeIncrement& timeIncrement )
{
  using namespace Marmot::Detail::GaussLobattoInterface;

  constexpr int nLocal = 3 * NStations + 3;

  using Vector9d            = Eigen::Matrix< double, 9, 1 >;
  using VectorLocal         = Eigen::Matrix< double, nLocal, 1 >;
  using Matrix3dRowMajor    = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor    = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor   = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;
  using Matrix3x9RowMajor   = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;
  using MatrixLocalRowMajor = Eigen::Matrix< double, nLocal, nLocal, Eigen::RowMajor >;
  using MatrixLocalx21      = Eigen::Matrix< double, nLocal, 21, Eigen::RowMajor >;
  using Matrix21xLocal      = Eigen::Matrix< double, 21, nLocal, Eigen::RowMajor >;
  using Matrix21dRowMajor   = Eigen::Matrix< double, 21, 21, Eigen::RowMajor >;

  for ( auto& material : stationMaterials ) {
    if ( !material ) {
      throw std::logic_error( "MarmotGaussLobattoBBarInterfaceMaterialHypoElastic has no base material." );
    }
  }

  const Eigen::Map< const Vector3d > normalMap( deformation.normal.data() );
  const Eigen::Map< const Vector3d > separationMap( deformation.separationVector.data() );

  Vector3d     normal     = normalMap;
  const double normalNorm = normal.norm();
  if ( normalNorm <= 1.0e-12 ) {
    throw std::invalid_argument( "MarmotGaussLobattoBBarInterfaceMaterialHypoElastic: interface normal is zero." );
  }
  normal /= normalNorm;

  const InterfaceGeometry geometry    = evaluateInterfaceGeometry( normal, separationMap, h );
  const double            ell         = geometry.normalSeparation;
  const Vector3d&         dTangential = geometry.tangentialSeparation;

  const Eigen::Map< const Eigen::Matrix< double, 6, 1 > >  dU( deformation.dU.data() );
  const Eigen::Map< const Eigen::Matrix< double, 18, 1 > > dSurfaceGradient( deformation.dSurfaceStrain.data() );

  const Vector3d w = dU.template segment< 3 >( 0 ) - dU.template segment< 3 >( 3 );

  const Eigen::Map< const Matrix3dRowMajor > APlus( dSurfaceGradient.data() );
  const Eigen::Map< const Matrix3dRowMajor > AMinus( dSurfaceGradient.data() + 9 );
  const Matrix3dRowMajor                     ABar = 0.5 * ( APlus + AMinus );

  const Vector3d gBar = ( w - ABar * dTangential ) / ell;

  const auto lobattoXi     = LobattoRule< NStations >::xi();
  const auto lobattoWeight = LobattoRule< NStations >::weight();

  // The ONE functional difference from MarmotGaussLobattoInterfaceMaterialHypoElasticN:
  // each station's tangential surface gradient gets an added B-bar
  // volumetric correction, supplied by the element as
  // deformation.traceCorrection[alpha]. Treated as a fixed external input
  // for this call (not differentiated with respect to q) -- the standard
  // "frozen" B-bar tangent.
  std::array< Matrix3dRowMajor, NStations > AStation;
  std::array< double, NStations >           lambda, Nplus, Nminus;
  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    lambda[alpha]   = 0.5 * lobattoWeight[alpha];
    Nplus[alpha]    = 0.5 * ( 1.0 + lobattoXi[alpha] );
    Nminus[alpha]   = 0.5 * ( 1.0 - lobattoXi[alpha] );
    AStation[alpha] = Nminus[alpha] * AMinus + Nplus[alpha] * APlus +
                      ( deformation.traceCorrection[alpha] / 3.0 ) * Matrix3dRowMajor::Identity();
  }

  std::array< double*, NStations > stationStateVars;
  std::array< int, NStations >     nStationStateVars;
  std::array< double*, NStations > normalGradientPtr;
  std::array< double*, NStations > stationStressPtr;

  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    stationStateVars[alpha]  = stateLayout.getPtr( state.stateVars, stationStateVarsName( alpha ) );
    nStationStateVars[alpha] = stationMaterials[alpha]->getNumberOfRequiredStateVars();
    normalGradientPtr[alpha] = stateLayout.getPtr( state.stateVars, normalGradientName( alpha ) );
    stationStressPtr[alpha]  = stateLayout.getPtr( state.stateVars, stationStressName( alpha ) );
  }
  double* commonTractionPtr = stateLayout.getPtr( state.stateVars, "commonTraction" );

  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld, timeIncrement.dT };

  const Matrix3x9RowMajor Rt = stressToForce( normal, 1.0, 1.0 );
  const Matrix9x3RowMajor Bn = jumpToGradient( normal, 1.0 );

  const double Qc = std::max( 1.0, std::abs( materialProperties[0] ) );

  struct LocalEvaluation {
    std::array< StationTrial, NStations > station;
    std::array< Vector3d, NStations >     tAlpha;
    VectorLocal                           R;
    double                                rNorm;
  };

  auto evaluateLocal = [&]( const VectorLocal& yTrial ) -> LocalEvaluation {
    LocalEvaluation evaluation;
    const Vector3d  t = yTrial.template segment< 3 >( 3 * NStations );

    Vector3d gWeightedSum  = Vector3d::Zero();
    double   tractionScale = 1.0;

    for ( int alpha = 0; alpha < NStations; ++alpha ) {
      const Vector3d         g                    = yTrial.template segment< 3 >( 3 * alpha );
      const Matrix3dRowMajor G                    = AStation[alpha] + g * normal.transpose();
      const Matrix3dRowMajor strainIncrement      = 0.5 * ( G + G.transpose() );
      const Marmot::Vector6d strainIncrementVoigt = Marmot::ContinuumMechanics::VoigtNotation::strainToVoigt(
        strainIncrement );
      const Matrix3dRowMajor stressCurrent = Eigen::Map< const Matrix3dRowMajor >( stationStressPtr[alpha] );

      evaluation.station[alpha] = evaluateStationTrial( *stationMaterials[alpha],
                                                        stationStateVars[alpha],
                                                        nStationStateVars[alpha],
                                                        stressCurrent,
                                                        strainIncrementVoigt,
                                                        timeInfo );

      evaluation.tAlpha[alpha]                        = evaluation.station[alpha].stress * normal;
      evaluation.R.template segment< 3 >( 3 * alpha ) = evaluation.tAlpha[alpha] - t;
      gWeightedSum += lambda[alpha] * g;
      tractionScale = std::max( tractionScale, evaluation.tAlpha[alpha].norm() );
    }

    evaluation.R.template segment< 3 >( 3 * NStations ) = gWeightedSum - gBar;

    const double tractionResidualNorm = evaluation.R.template segment< 3 * NStations >( 0 )
                                          .template lpNorm< Eigen::Infinity >();
    const double compatResidualNorm = evaluation.R.template segment< 3 >( 3 * NStations ).norm();
    evaluation.rNorm                = std::max( tractionResidualNorm, Qc * compatResidualNorm ) / tractionScale;

    return evaluation;
  };

  auto attemptLocalNewton = [&]( const VectorLocal& yInitial ) {
    struct Attempt {
      VectorLocal     yCurrent;
      LocalEvaluation current;
      bool            converged;
    };

    VectorLocal     yCurrent  = yInitial;
    LocalEvaluation current   = evaluateLocal( yCurrent );
    bool            converged = current.rNorm <= localNewtonTolerance;

    for ( int iteration = 0; iteration < maxLocalNewtonIterations && !converged; ++iteration ) {
      MatrixLocalRowMajor J = MatrixLocalRowMajor::Zero();
      for ( int alpha = 0; alpha < NStations; ++alpha ) {
        const Matrix3dRowMajor Qalpha                        = Rt * current.station[alpha].CFull * Bn;
        J.template block< 3, 3 >( 3 * alpha, 3 * alpha )     = Qalpha;
        J.template block< 3, 3 >( 3 * alpha, 3 * NStations ) = -Matrix3dRowMajor::Identity();
        J.template block< 3, 3 >( 3 * NStations, 3 * alpha ) = lambda[alpha] * Matrix3dRowMajor::Identity();
      }

      const VectorLocal dy = J.fullPivLu().solve( -current.R );

      double          alphaStep = 1.0;
      LocalEvaluation candidate;
      bool            accepted = false;
      bool            stuck    = true;
      for ( int lineSearchIter = 0; lineSearchIter < 24; ++lineSearchIter ) {
        candidate = evaluateLocal( yCurrent + alphaStep * dy );
        if ( candidate.rNorm <= localNewtonTolerance || candidate.rNorm < current.rNorm ) {
          accepted = true;
          stuck    = false;
          break;
        }
        alphaStep *= 0.5;
      }

      if ( stuck ) {
        break;
      }

      yCurrent += alphaStep * dy;
      current = accepted ? candidate : evaluateLocal( yCurrent );

      if ( current.rNorm <= localNewtonTolerance ) {
        converged = true;
      }
    }

    return Attempt{ yCurrent, current, converged };
  };

  VectorLocal yWarmStart;
  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    yWarmStart.template segment< 3 >( 3 * alpha ) = Eigen::Map< const Vector3d >( normalGradientPtr[alpha] );
  }
  yWarmStart.template segment< 3 >( 3 * NStations ) = Eigen::Map< const Vector3d >( commonTractionPtr );

  auto attempt = attemptLocalNewton( yWarmStart );

  if ( !attempt.converged ) {
    VectorLocal yFresh;
    for ( int alpha = 0; alpha < NStations; ++alpha ) {
      yFresh.template segment< 3 >( 3 * alpha ) = gBar;
    }
    yFresh.template segment< 3 >( 3 * NStations ) = Vector3d::Zero();

    const auto freshAttempt = attemptLocalNewton( yFresh );
    if ( freshAttempt.converged || freshAttempt.current.rNorm < attempt.current.rNorm ) {
      attempt = freshAttempt;
    }
  }

  VectorLocal&     yCurrent  = attempt.yCurrent;
  LocalEvaluation& current   = attempt.current;
  const bool       converged = attempt.converged;

  if ( !converged ) {
    throw Marmot::StressUpdateFailed(
      "MarmotGaussLobattoBBarInterfaceMaterialHypoElastic: local traction-equilibrium solve did not converge." );
  }

  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    Eigen::Map< Vector3d > normalGradientMap( normalGradientPtr[alpha] );
    normalGradientMap = yCurrent.template segment< 3 >( 3 * alpha );
    std::copy( current.station[alpha].trialStateVars.begin(),
               current.station[alpha].trialStateVars.end(),
               stationStateVars[alpha] );
    Eigen::Map< Matrix3dRowMajor > stationStressMap( stationStressPtr[alpha] );
    stationStressMap = current.station[alpha].stress;
  }
  Eigen::Map< Vector3d > commonTractionMap( commonTractionPtr );
  commonTractionMap = yCurrent.template segment< 3 >( 3 * NStations );

  const Matrix3x9RowMajor Mdtau = tangentialContraction( dTangential );

  MatrixLocalRowMajor Ry = MatrixLocalRowMajor::Zero();
  MatrixLocalx21      Rq = MatrixLocalx21::Zero();

  std::array< Matrix9dRowMajor, NStations >  CFullConverged;
  std::array< Matrix3x9RowMajor, NStations > RtCFull;
  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    CFullConverged[alpha] = current.station[alpha].CFull;
    RtCFull[alpha]        = Rt * CFullConverged[alpha];

    const Matrix3dRowMajor Qalpha                         = RtCFull[alpha] * Bn;
    Ry.template block< 3, 3 >( 3 * alpha, 3 * alpha )     = Qalpha;
    Ry.template block< 3, 3 >( 3 * alpha, 3 * NStations ) = -Matrix3dRowMajor::Identity();
    Ry.template block< 3, 3 >( 3 * NStations, 3 * alpha ) = lambda[alpha] * Matrix3dRowMajor::Identity();

    Rq.template block< 3, 9 >( 3 * alpha, 3 )  = Nplus[alpha] * RtCFull[alpha];
    Rq.template block< 3, 9 >( 3 * alpha, 12 ) = Nminus[alpha] * RtCFull[alpha];
  }
  Rq.template block< 3, 3 >( 3 * NStations, 0 )  = -Matrix3dRowMajor::Identity() / ell;
  Rq.template block< 3, 9 >( 3 * NStations, 3 )  = ( 0.5 / ell ) * Mdtau;
  Rq.template block< 3, 9 >( 3 * NStations, 12 ) = ( 0.5 / ell ) * Mdtau;

  Matrix21dRowMajor pq = Matrix21dRowMajor::Zero();
  Matrix21xLocal    py = Matrix21xLocal::Zero();

  Matrix9dRowMajor pApAp = Matrix9dRowMajor::Zero();
  Matrix9dRowMajor pApAm = Matrix9dRowMajor::Zero();
  Matrix9dRowMajor pAmAm = Matrix9dRowMajor::Zero();
  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    pApAp += ( h * lambda[alpha] * Nplus[alpha] * Nplus[alpha] ) * CFullConverged[alpha];
    pApAm += ( h * lambda[alpha] * Nplus[alpha] * Nminus[alpha] ) * CFullConverged[alpha];
    pAmAm += ( h * lambda[alpha] * Nminus[alpha] * Nminus[alpha] ) * CFullConverged[alpha];
  }
  pq.template block< 9, 9 >( 3, 3 )   = pApAp;
  pq.template block< 9, 9 >( 3, 12 )  = pApAm;
  pq.template block< 9, 9 >( 12, 3 )  = pApAm;
  pq.template block< 9, 9 >( 12, 12 ) = pAmAm;

  py.template block< 3, 3 >( 0, 3 * NStations ) = ( h / ell ) * Matrix3dRowMajor::Identity();
  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    const Matrix9x3RowMajor CB                 = CFullConverged[alpha] * Bn;
    py.template block< 9, 3 >( 3, 3 * alpha )  = ( h * lambda[alpha] * Nplus[alpha] ) * CB;
    py.template block< 9, 3 >( 12, 3 * alpha ) = ( h * lambda[alpha] * Nminus[alpha] ) * CB;
  }
  py.template block< 9, 3 >( 3, 3 * NStations )  = -( 0.5 * h / ell ) * Mdtau.transpose();
  py.template block< 9, 3 >( 12, 3 * NStations ) = -( 0.5 * h / ell ) * Mdtau.transpose();

  const MatrixLocalx21    RyInvRq = Ry.fullPivLu().solve( Rq );
  const Matrix21dRowMajor Kcond   = pq - py * RyInvRq;

  const Vector3d tShared = yCurrent.template segment< 3 >( 3 * NStations );

  Vector9d sPlusSum  = Vector9d::Zero();
  Vector9d sMinusSum = Vector9d::Zero();
  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    const Eigen::Map< const Vector9d > stressVec( current.station[alpha].stress.data() );
    sPlusSum += ( h * lambda[alpha] * Nplus[alpha] ) * stressVec;
    sMinusSum += ( h * lambda[alpha] * Nminus[alpha] ) * stressVec;
  }
  const Vector9d mdtauTtShared = Mdtau.transpose() * tShared;
  sPlusSum -= ( 0.5 * h / ell ) * mdtauTtShared;
  sMinusSum -= ( 0.5 * h / ell ) * mdtauTtShared;

  Eigen::Map< Vector3d >( state.generalizedForce.data() )   = ( h / ell ) * tShared;
  Eigen::Map< Vector9d >( state.surfaceStressPlus.data() )  = sPlusSum;
  Eigen::Map< Vector9d >( state.surfaceStressMinus.data() ) = sMinusSum;

  Eigen::Map< Matrix3dRowMajor >( tangents.Q_ww.data() )   = Kcond.template block< 3, 3 >( 0, 0 );
  Eigen::Map< Matrix3x9RowMajor >( tangents.Q_wAp.data() ) = Kcond.template block< 3, 9 >( 0, 3 );
  Eigen::Map< Matrix3x9RowMajor >( tangents.Q_wAm.data() ) = Kcond.template block< 3, 9 >( 0, 12 );

  Eigen::Map< Matrix9x3RowMajor >( tangents.Q_Apw.data() ) = Kcond.template block< 9, 3 >( 3, 0 );
  Eigen::Map< Matrix9dRowMajor >( tangents.Q_ApAp.data() ) = Kcond.template block< 9, 9 >( 3, 3 );
  Eigen::Map< Matrix9dRowMajor >( tangents.Q_ApAm.data() ) = Kcond.template block< 9, 9 >( 3, 12 );

  Eigen::Map< Matrix9x3RowMajor >( tangents.Q_Amw.data() ) = Kcond.template block< 9, 3 >( 12, 0 );
  Eigen::Map< Matrix9dRowMajor >( tangents.Q_AmAp.data() ) = Kcond.template block< 9, 9 >( 12, 3 );
  Eigen::Map< Matrix9dRowMajor >( tangents.Q_AmAm.data() ) = Kcond.template block< 9, 9 >( 12, 12 );
}

template < int NStations >
void MarmotGaussLobattoBBarInterfaceMaterialHypoElasticN< NStations >::initializeYourself( double* stateVars,
                                                                                           int     nStateVars )
{
  using namespace Marmot::Detail::GaussLobattoInterface;
  (void)nStateVars;

  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    if ( !stationMaterials[alpha] ) {
      continue;
    }
    stationMaterials[alpha]->initializeYourself( stateLayout.getPtr( stateVars, stationStateVarsName( alpha ) ),
                                                 stationMaterials[alpha]->getNumberOfRequiredStateVars() );

    double* g = stateLayout.getPtr( stateVars, normalGradientName( alpha ) );
    g[0] = g[1] = g[2] = 0.0;

    double* s = stateLayout.getPtr( stateVars, stationStressName( alpha ) );
    std::fill( s, s + 9, 0.0 );
  }

  double* t = stateLayout.getPtr( stateVars, "commonTraction" );
  t[0] = t[1] = t[2] = 0.0;
}

template < int NStations >
double MarmotGaussLobattoBBarInterfaceMaterialHypoElasticN< NStations >::getDensity()
{
  if ( !stationMaterials[0] ) {
    return -1;
  }

  return stationMaterials[0]->getDensity( nullptr );
}
