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

#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotGaussLobattoInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"
#include "Marmot/MarmotStateHelpers.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @file MarmotWarpingInterfaceMaterialHypoElastic.h
 * @brief Multi-station interface material with SYMMETRIC (parabolic) and
 * ANTISYMMETRIC (cubic) through-thickness WARPING micro-fields.
 *
 * WHAT THIS ADDS OVER THE GAUSS-LOBATTO MATERIAL
 * ----------------------------------------------
 * MarmotGaussLobattoInterfaceMaterialHypoElasticN resolves the through-
 * thickness direction with NStations constitutive points, but its TANGENTIAL
 * surface gradient is interpolated LINEARLY between the two faces,
 *
 *     A^(alpha) = Nminus_alpha A- + Nplus_alpha A+ ,
 *
 * so the in-plane strain is forced to vary linearly (constant + bending)
 * across the layer. Its per-station normal gradients g^(alpha) are, by
 * contrast, already free (subject to one weighted-mean constraint), so the
 * TRANSVERSE-SHEAR profile is unrestricted there.
 *
 * This material lifts the remaining restriction by adding two warping
 * micro-fields to the through-thickness displacement expansion. Writing
 * zeta in [-1,1] for the normalized thickness coordinate,
 *
 *     u(x_s, zeta) = ubar(x_s) + (zeta/2) [u](x_s)
 *                    + phi_s(zeta) w_s(x_s)  +  phi_a(zeta) w_a(x_s),
 *
 *     phi_s(zeta) = 1 - zeta^2      (EVEN  -> symmetric,     parabolic)
 *     phi_a(zeta) = zeta - zeta^3   (ODD   -> antisymmetric, cubic)
 *
 * Both profiles vanish at zeta = +/-1, so the warping fields leave the FACE
 * displacements -- and therefore the mesh jump [u] and the coupling to the
 * surrounding bulk elements -- completely untouched. That is what makes them
 * strictly internal, hence condensable at element level.
 *
 * WHY ONLY THE TANGENTIAL GRADIENT IS ENRICHED (the crux -- do not "fix" this)
 * ---------------------------------------------------------------------------
 * The full gradient of a warping term phi(zeta) w(x_s) is
 *
 *     grad( phi w ) = phi(zeta) grad_s w  +  (2/h) phi'(zeta) w (x) n .
 *
 * The second (normal-gradient) part would be REDUNDANT here: the Gauss-
 * Lobatto station gradients g^(alpha) are already free apart from the single
 * weighted-compatibility constraint sum_alpha lambda_alpha g^(alpha) = gbar,
 * and BOTH warping profiles have exactly zero weighted mean derivative,
 *
 *     sum_alpha lambda_alpha phi_s'(xi_alpha) = sum_alpha lambda_alpha (-2 xi_alpha)     = 0,
 *     sum_alpha lambda_alpha phi_a'(xi_alpha) = sum_alpha lambda_alpha (1 - 3 xi_alpha^2) = 0,
 *
 * (the first because the Lobatto abscissae/weights are symmetric, the second
 * because sum lambda xi^2 = 1/3 exactly for every rule used here). So the
 * normal-gradient content of both warping modes lies entirely inside the
 * constraint-preserving subspace that g^(alpha) already spans, and the local
 * traction-equilibrium solve picks it up on its own. Feeding it in a second
 * time through w would make the combined local/internal system exactly
 * singular. Nothing is lost by omitting it; something is broken by adding it.
 *
 * Consequently the material's NEW generalized strains are the two warping
 * SURFACE GRADIENTS
 *
 *     Ws = grad_s w_s ,      Wa = grad_s w_a       (3x3 each),
 *
 * entering the station tangential gradient as
 *
 *     A^(alpha) = Nminus_alpha A- + Nplus_alpha A+ + phi_s(xi_alpha) Ws
 *                                                  + phi_a(xi_alpha) Wa .
 *
 * The element supplies Ws and Wa; a constant warping amplitude produces
 * Ws = Wa = 0 and is therefore (correctly) invisible here.
 *
 * STATION COUNT
 * -------------
 * NStations >= 4 is REQUIRED. The three-point Lobatto rule samples
 * xi = {-1, 0, +1}, where phi_a = zeta - zeta^3 vanishes at every single
 * station -- the antisymmetric warping mode would be identically invisible
 * and its stiffness block exactly singular. Five stations is the production
 * default, matching the Gauss-Lobatto material.
 *
 * REDUCTION PROPERTY
 * ------------------
 * With Ws = Wa = 0 this material reproduces
 * MarmotGaussLobattoInterfaceMaterialHypoElasticN<NStations> exactly (same
 * local problem, same generalized stresses, same 21x21 sub-block of the
 * condensed tangent). The geometric d_tau correction deliberately keeps
 * using the FACE-interpolated mean ABar = (A+ + A-)/2, exactly as the
 * Gauss-Lobatto material does, so that the reduction is bit-exact.
 *
 * GENERALIZED PAIR
 * ----------------
 *     q = ( w, A+, A-, Ws, Wa )       3 + 9 + 9 + 9 + 9 = 39
 *     p = ( f, S+, S-, S_Ws, S_Wa )
 *
 * with, for each tensor slot X in {A+, A-, Ws, Wa} and its through-thickness
 * profile coefficient c^X_alpha in {Nplus, Nminus, phi_s, phi_a},
 *
 *     S_X = h sum_alpha lambda_alpha c^X_alpha sigma^(alpha)
 *           - [X in {A+,A-}] (h/(2 ell)) Mdtau^T t .
 *
 * The tangent is returned as ONE dense condensed 39x39 block rather than as
 * 25 named sub-blocks; the element slices it.
 */
template < int NStations = 5 >
class MarmotWarpingInterfaceMaterialHypoElasticN {

  static_assert( NStations >= 4,
                 "MarmotWarpingInterfaceMaterialHypoElastic requires at least 4 through-thickness stations: the "
                 "3-point Gauss-Lobatto rule samples xi = {-1,0,1}, where the antisymmetric warping profile "
                 "phi_a = xi - xi^3 vanishes identically, leaving its stiffness block singular." );

public:
  static constexpr int nStations = NStations;

  /** Number of 3x3 tensor slots in the generalized strain: A+, A-, Ws, Wa. */
  static constexpr int nTensorSlots = 4;

  static constexpr int nGeneralized = 3 + 9 * nTensorSlots; // 39

  static constexpr int offsetJump                 = 0;
  static constexpr int offsetSurfacePlus          = 3;
  static constexpr int offsetSurfaceMinus         = 12;
  static constexpr int offsetWarpingSymmetric     = 21;
  static constexpr int offsetWarpingAntisymmetric = 30;

  /** Offsets of the four tensor slots, in the slot order used everywhere below. */
  static constexpr std::array< int, nTensorSlots > slotOffset = { offsetSurfacePlus,
                                                                  offsetSurfaceMinus,
                                                                  offsetWarpingSymmetric,
                                                                  offsetWarpingAntisymmetric };

  using CondensedTangent  = Eigen::Matrix< double, nGeneralized, nGeneralized, Eigen::RowMajor >;
  using GeneralizedVector = Eigen::Matrix< double, nGeneralized, 1 >;

protected:
  const double*                                                         materialProperties;
  const int                                                             nMaterialProperties;
  double                                                                h = 0.0;
  std::vector< double >                                                 baseMaterialProperties;
  std::array< std::unique_ptr< MarmotMaterialHypoElastic >, NStations > stationMaterials;

public:
  const int materialNumber;

  MarmotStateLayoutDynamic stateLayout;

  double characteristicElementLength = 0.0;

  /** Generalized stress p, laid out as ( f(3), S+(9), S-(9), S_Ws(9), S_Wa(9) ). */
  struct State {
    Eigen::Map< GeneralizedVector > generalizedStress;
    double*                         stateVars;

    State( double* generalizedStress_, double* stateVars_ )
      : generalizedStress( generalizedStress_ ), stateVars( stateVars_ )
    {
    }
  };

  /** Condensed d(p)/d(q), row-major, over the 39-component generalized pair. */
  struct Tangents {
    Eigen::Map< CondensedTangent > condensed;

    explicit Tangents( double* condensed_ ) : condensed( condensed_ ) {}
  };

  struct Deformation {
    Eigen::Map< const Eigen::Matrix< double, 6, 1 > >  dU;             /**< (u+, u-) increments. */
    Eigen::Map< const Eigen::Matrix< double, 18, 1 > > dSurfaceStrain; /**< (A+, A-) increments. */
    Eigen::Map< const Eigen::Matrix< double, 18, 1 > > dWarpingStrain; /**< (Ws, Wa) increments. */
    Eigen::Map< const Eigen::Matrix< double, 3, 1 > >  normal;
    Eigen::Map< const Eigen::Matrix< double, 3, 1 > >  separationVector;

    Deformation( const double* dU_,
                 const double* dSurfaceStrain_,
                 const double* dWarpingStrain_,
                 const double* normal_,
                 const double* separationVector_ )
      : dU( dU_ ),
        dSurfaceStrain( dSurfaceStrain_ ),
        dWarpingStrain( dWarpingStrain_ ),
        normal( normal_ ),
        separationVector( separationVector_ )
    {
    }

    /** Coincident-face convenience overload: d = 0, ell = h, d_tau = 0. */
    Deformation( const double* dU_,
                 const double* dSurfaceStrain_,
                 const double* dWarpingStrain_,
                 const double* normal_ )
      : Deformation( dU_, dSurfaceStrain_, dWarpingStrain_, normal_, zeroVector() )
    {
    }

  private:
    static const double* zeroVector()
    {
      static const double zero[3] = { 0.0, 0.0, 0.0 };
      return zero;
    }
  };

  struct TimeIncrement {
    double timeOld;
    double dT;
  };

  MarmotWarpingInterfaceMaterialHypoElasticN( const std::string& materialName,
                                              const double*      matProperties_,
                                              int                nMaterialProperties_,
                                              int                materialNumber_ );

  virtual ~MarmotWarpingInterfaceMaterialHypoElasticN() = default;

  void setCharacteristicElementLength( double length );

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

  /** Through-thickness profile of the symmetric (parabolic) warping mode. */
  static constexpr double warpingProfileSymmetric( double xi ) { return 1.0 - xi * xi; }

  /** Through-thickness profile of the antisymmetric (cubic) warping mode. */
  static constexpr double warpingProfileAntisymmetric( double xi ) { return xi - xi * xi * xi; }
};

/** Production rule: five Gauss-Lobatto stations, matching the GL material. */
using MarmotWarpingInterfaceMaterialHypoElastic = MarmotWarpingInterfaceMaterialHypoElasticN< 5 >;

// ---------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------

template < int NStations >
MarmotWarpingInterfaceMaterialHypoElasticN< NStations >::MarmotWarpingInterfaceMaterialHypoElasticN(
  const std::string& materialName,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : materialProperties( matProperties_ ), nMaterialProperties( nMaterialProperties_ ), materialNumber( materialNumber_ )
{
  using namespace Marmot::Detail::GaussLobattoInterface;

  if ( nMaterialProperties < 3 ) {
    throw std::invalid_argument(
      "MarmotWarpingInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h." );
  }

  h = materialProperties[2];
  if ( h <= 0.0 ) {
    throw std::invalid_argument( "MarmotWarpingInterfaceMaterialHypoElastic requires h > 0." );
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
void MarmotWarpingInterfaceMaterialHypoElasticN< NStations >::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  for ( auto& material : stationMaterials ) {
    if ( material ) {
      material->setCharacteristicElementLength( length );
    }
  }
}

template < int NStations >
void MarmotWarpingInterfaceMaterialHypoElasticN< NStations >::computeStress( State&               state,
                                                                             Tangents&            tangents,
                                                                             const Deformation&   deformation,
                                                                             const TimeIncrement& timeIncrement )
{
  using namespace Marmot::Detail::GaussLobattoInterface;

  constexpr int nLocal = 3 * NStations + 3;
  constexpr int nGen   = nGeneralized;

  using Vector3d            = Eigen::Matrix< double, 3, 1 >;
  using Vector9d            = Eigen::Matrix< double, 9, 1 >;
  using VectorLocal         = Eigen::Matrix< double, nLocal, 1 >;
  using Matrix3dRowMajor    = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
  using Matrix9dRowMajor    = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;
  using Matrix9x3RowMajor   = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >;
  using Matrix3x9RowMajor   = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >;
  using MatrixLocalRowMajor = Eigen::Matrix< double, nLocal, nLocal, Eigen::RowMajor >;
  using MatrixLocalxGen     = Eigen::Matrix< double, nLocal, nGen, Eigen::RowMajor >;
  using MatrixGenxLocal     = Eigen::Matrix< double, nGen, nLocal, Eigen::RowMajor >;

  for ( auto& material : stationMaterials ) {
    if ( !material ) {
      throw std::logic_error( "MarmotWarpingInterfaceMaterialHypoElastic has no base material." );
    }
  }

  Vector3d     normal     = deformation.normal;
  const double normalNorm = normal.norm();
  if ( normalNorm <= 1.0e-12 ) {
    throw std::invalid_argument( "MarmotWarpingInterfaceMaterialHypoElastic: interface normal is zero." );
  }
  normal /= normalNorm;

  const InterfaceGeometry geometry = evaluateInterfaceGeometry( normal, Vector3d( deformation.separationVector ), h );
  const double            ell      = geometry.normalSeparation;
  const Vector3d&         dTangential = geometry.tangentialSeparation;

  const Vector3d w = deformation.dU.template segment< 3 >( 0 ) - deformation.dU.template segment< 3 >( 3 );

  const Eigen::Map< const Matrix3dRowMajor > APlus( deformation.dSurfaceStrain.data() );
  const Eigen::Map< const Matrix3dRowMajor > AMinus( deformation.dSurfaceStrain.data() + 9 );
  const Eigen::Map< const Matrix3dRowMajor > WarpSymmetric( deformation.dWarpingStrain.data() );
  const Eigen::Map< const Matrix3dRowMajor > WarpAntisymmetric( deformation.dWarpingStrain.data() + 9 );

  // The d_tau geometric correction intentionally keeps the FACE-interpolated
  // mean, exactly as the Gauss-Lobatto material defines it, so that Ws=Wa=0
  // reproduces that material bit-for-bit.
  const Matrix3dRowMajor ABar = 0.5 * ( APlus + AMinus );
  const Vector3d         gBar = ( w - ABar * dTangential ) / ell;

  const auto lobattoXi     = LobattoRule< NStations >::xi();
  const auto lobattoWeight = LobattoRule< NStations >::weight();

  // Through-thickness profile coefficient of every tensor slot at every
  // station, in slot order { A+, A-, Ws, Wa }.
  std::array< std::array< double, nTensorSlots >, NStations > slotCoefficient;
  std::array< Matrix3dRowMajor, NStations >                   AStation;
  std::array< double, NStations >                             lambda;

  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    const double xi = lobattoXi[alpha];

    lambda[alpha] = 0.5 * lobattoWeight[alpha];

    slotCoefficient[alpha][0] = 0.5 * ( 1.0 + xi ); // Nplus  -> A+
    slotCoefficient[alpha][1] = 0.5 * ( 1.0 - xi ); // Nminus -> A-
    slotCoefficient[alpha][2] = warpingProfileSymmetric( xi );
    slotCoefficient[alpha][3] = warpingProfileAntisymmetric( xi );

    AStation[alpha] = slotCoefficient[alpha][0] * APlus + slotCoefficient[alpha][1] * AMinus +
                      slotCoefficient[alpha][2] * WarpSymmetric + slotCoefficient[alpha][3] * WarpAntisymmetric;
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

  // ---- local problem y = ( g^(1..N), t ): unchanged from the Gauss-Lobatto
  // material. The warping only shifts A^(alpha), it adds no local unknown.
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
      bool            stuck = true;
      for ( int lineSearchIter = 0; lineSearchIter < 24; ++lineSearchIter ) {
        candidate = evaluateLocal( yCurrent + alphaStep * dy );
        if ( candidate.rNorm <= localNewtonTolerance || candidate.rNorm < current.rNorm ) {
          stuck = false;
          break;
        }
        alphaStep *= 0.5;
      }

      if ( stuck ) {
        break;
      }

      yCurrent += alphaStep * dy;
      current = candidate;

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

  VectorLocal&     yCurrent = attempt.yCurrent;
  LocalEvaluation& current  = attempt.current;

  if ( !attempt.converged ) {
    throw Marmot::StressUpdateFailed(
      "MarmotWarpingInterfaceMaterialHypoElastic: local traction-equilibrium solve did not converge." );
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

  // ---- condensation of the local unknowns onto the 39-component pair ----
  const Matrix3x9RowMajor Mdtau = tangentialContraction( dTangential );

  MatrixLocalRowMajor Ry = MatrixLocalRowMajor::Zero();
  MatrixLocalxGen     Rq = MatrixLocalxGen::Zero();

  std::array< Matrix9dRowMajor, NStations >  CFullConverged;
  std::array< Matrix3x9RowMajor, NStations > RtCFull;

  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    CFullConverged[alpha] = current.station[alpha].CFull;
    RtCFull[alpha]        = Rt * CFullConverged[alpha];

    const Matrix3dRowMajor Qalpha                         = RtCFull[alpha] * Bn;
    Ry.template block< 3, 3 >( 3 * alpha, 3 * alpha )     = Qalpha;
    Ry.template block< 3, 3 >( 3 * alpha, 3 * NStations ) = -Matrix3dRowMajor::Identity();
    Ry.template block< 3, 3 >( 3 * NStations, 3 * alpha ) = lambda[alpha] * Matrix3dRowMajor::Identity();

    for ( int slot = 0; slot < nTensorSlots; ++slot ) {
      Rq.template block< 3, 9 >( 3 * alpha, slotOffset[slot] ) = slotCoefficient[alpha][slot] * RtCFull[alpha];
    }
  }

  // The weighted-compatibility row sees only w and the FACE gradients: gBar is
  // built from the mesh jump and ABar, and the warping does not enter it.
  Rq.template block< 3, 3 >( 3 * NStations, offsetJump )         = -Matrix3dRowMajor::Identity() / ell;
  Rq.template block< 3, 9 >( 3 * NStations, offsetSurfacePlus )  = ( 0.5 / ell ) * Mdtau;
  Rq.template block< 3, 9 >( 3 * NStations, offsetSurfaceMinus ) = ( 0.5 / ell ) * Mdtau;

  CondensedTangent pq = CondensedTangent::Zero();
  MatrixGenxLocal  py = MatrixGenxLocal::Zero();

  for ( int rowSlot = 0; rowSlot < nTensorSlots; ++rowSlot ) {
    for ( int colSlot = 0; colSlot < nTensorSlots; ++colSlot ) {
      Matrix9dRowMajor block = Matrix9dRowMajor::Zero();
      for ( int alpha = 0; alpha < NStations; ++alpha ) {
        block += ( h * lambda[alpha] * slotCoefficient[alpha][rowSlot] * slotCoefficient[alpha][colSlot] ) *
                 CFullConverged[alpha];
      }
      pq.template block< 9, 9 >( slotOffset[rowSlot], slotOffset[colSlot] ) = block;
    }
  }

  py.template block< 3, 3 >( offsetJump, 3 * NStations ) = ( h / ell ) * Matrix3dRowMajor::Identity();
  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    const Matrix9x3RowMajor CB = CFullConverged[alpha] * Bn;
    for ( int slot = 0; slot < nTensorSlots; ++slot ) {
      py.template block< 9, 3 >( slotOffset[slot], 3 * alpha ) = ( h * lambda[alpha] * slotCoefficient[alpha][slot] ) *
                                                                 CB;
    }
  }
  // Only the face slots inherit the d_tau term, because only they appear in gBar.
  py.template block< 9, 3 >( offsetSurfacePlus, 3 * NStations )  = -( 0.5 * h / ell ) * Mdtau.transpose();
  py.template block< 9, 3 >( offsetSurfaceMinus, 3 * NStations ) = -( 0.5 * h / ell ) * Mdtau.transpose();

  const MatrixLocalxGen RyInvRq = Ry.fullPivLu().solve( Rq );

  tangents.condensed = pq - py * RyInvRq;

  // ---- generalized stresses ----
  const Vector3d tShared = yCurrent.template segment< 3 >( 3 * NStations );

  state.generalizedStress.setZero();
  state.generalizedStress.template segment< 3 >( offsetJump ) = ( h / ell ) * tShared;

  for ( int slot = 0; slot < nTensorSlots; ++slot ) {
    Vector9d sum = Vector9d::Zero();
    for ( int alpha = 0; alpha < NStations; ++alpha ) {
      const Eigen::Map< const Vector9d > stressVec( current.station[alpha].stress.data() );
      sum += ( h * lambda[alpha] * slotCoefficient[alpha][slot] ) * stressVec;
    }
    if ( slotOffset[slot] == offsetSurfacePlus || slotOffset[slot] == offsetSurfaceMinus ) {
      sum -= ( 0.5 * h / ell ) * ( Mdtau.transpose() * tShared );
    }
    state.generalizedStress.template segment< 9 >( slotOffset[slot] ) = sum;
  }
}

template < int NStations >
void MarmotWarpingInterfaceMaterialHypoElasticN< NStations >::initializeYourself( double* stateVars, int nStateVars )
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
double MarmotWarpingInterfaceMaterialHypoElasticN< NStations >::getDensity()
{
  if ( !stationMaterials[0] ) {
    return -1;
  }

  return stationMaterials[0]->getDensity( nullptr );
}
