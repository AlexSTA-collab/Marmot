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
 * @file MarmotWarpingStabPressureInterfaceMaterialHypoElastic.h
 * @brief Multi-station interface material carrying BOTH through-thickness
 * WARPING micro-fields AND an independent (mixed) hydrostatic pressure.
 *
 * WHY BOTH, AND WHY THEY DO NOT OVERLAP
 * -------------------------------------
 * The two enrichments were measured to fix two different errors, and neither
 * fixes the other's:
 *
 *   * WARPING (MarmotWarpingInterfaceMaterialHypoElastic) removes a KINEMATIC
 *     error. The Gauss-Lobatto material forces the in-plane strain to vary
 *     LINEARLY across the layer; the interior of a real Cauchy layer does not.
 *     Adding the parabolic/cubic micro-fields cut the displacement-jump error
 *     38% and the overall L2 error 44% at h = 0.01 -- but it left the pressure
 *     checkerboard essentially untouched (oscillation metric 1.400e-1 vs
 *     1.516e-1 for GLIQUAD4, i.e. 7.6%).
 *
 *   * MIXED PRESSURE (MarmotStabPressureInterfaceMaterialHypoElastic) removes a
 *     CONSTITUTIVE/INF-SUP error. Near the incompressible limit reached by
 *     near-perfectly-plastic von Mises (K/mu_plastic ~ 1e7) a displacement-only
 *     formulation evaluates p = K tr(eps) with tr(eps) ~ 1e-4, amplifying
 *     discretisation noise into an element-to-element checkerboard. Making p an
 *     independent nodal field collapsed the same metric to 1.227e-4, i.e. by
 *     three orders of magnitude.
 *
 * They are complementary because they act on ORTHOGONAL parts of the strain:
 * warping enriches the through-thickness profile of the DEVIATORIC response,
 * the mixed field replaces the VOLUMETRIC one. This material composes them.
 *
 * FORMULATION
 * -----------
 * Stress used for equilibrium, at every through-thickness station:
 *
 *     sigma~^(alpha) = dev( sigma^(alpha) )  -  p I ,
 *
 * with p the independent pressure -- ONE value per interface point, constant
 * across the thickness (it is a midsurface nodal field). The through-thickness
 * equilibrium condition therefore reads sigma~^(alpha) n = t~ for every alpha,
 * which -- because p I n = p n is the SAME at every station -- is identical to
 *
 *     dev( sigma^(alpha) ) n  =  t ,        t~ = t - p n .
 *
 * So the local traction-equilibrium solve is EXACTLY INDEPENDENT OF p; p only
 * shifts the shared traction. This is the multi-station generalisation of the
 * d gamma / dp = 0 property of the two-station mixed material, and it is what
 * keeps the mixed coupling blocks purely geometric.
 *
 * The volumetric law is returned as a residual for the element to enforce
 * weakly (the Brezzi-Pitkaranta pressure-gradient stabilisation lives in the
 * element, since it needs the mesh size):
 *
 *     R_p = h [ sum_alpha lambda_alpha tr( d eps^(alpha) )  +  dp / K ] .
 *
 * INCREMENTAL, exactly as in the two-station material: both the volumetric
 * strain and p accumulate from zero, so enforcing the increment each step is
 * equivalent to the total relation p = -K tr(eps), and it is consistent with
 * the hypoelastic (rate) form of the deviatoric response. Mixing the strain
 * INCREMENT with the TOTAL pressure would under-predict p by the number of
 * load steps.
 *
 * Note that the weighted station sum collapses to a purely geometric
 * expression. With sum_alpha lambda_alpha of the four profile coefficients
 * equal to ( 1/2, 1/2, 2/3, 0 ) for ( A+, A-, Ws, Wa ), and with the station
 * gradients constrained by sum_alpha lambda_alpha g^(alpha) = gbar,
 *
 *     sum_alpha lambda_alpha tr( d eps^(alpha) )
 *          = tr(ABar) + (2/3) tr(Ws) + gbar . n .
 *
 * The ANTISYMMETRIC warping is invisible to the volumetric residual (its
 * weighted mean vanishes) -- as it must be, since an odd profile transports no
 * net volume. The SYMMETRIC one is not, and its 2/3 weight is the reason this
 * material cannot be assembled by simply calling the two existing ones.
 *
 * RESTRICTION
 * -----------
 * Overriding the constitutive pressure is admissible only because the von Mises
 * yield surface is PRESSURE-INDEPENDENT: the deviatoric return mapping is
 * unaffected by p. For a pressure-DEPENDENT law (Drucker-Prager, cap models)
 * p would have to be fed INTO the return mapping and this construction is
 * invalid.
 *
 * GENERALIZED PAIR
 * ----------------
 *     q = ( w, A+, A-, Ws, Wa, p )      3 + 9 + 9 + 9 + 9 + 1 = 40
 *     p = ( f, S+, S-, S_Ws, S_Wa, R_p )
 *
 * with t~ = t - p n and, for each tensor slot X with profile coefficient
 * c^X_alpha in { Nplus, Nminus, phi_s, phi_a },
 *
 *     f   = (h/ell) t~
 *     S_X = h sum_alpha lambda_alpha c^X_alpha sigma~^(alpha)
 *           - [X in {A+,A-}] (h/(2 ell)) Mdtau^T t~ .
 *
 * Unlike the two-station material, R_p is returned ALREADY MULTIPLIED BY h, so
 * that the whole 40x40 tangent is one consistent dense block the element can
 * slice without rescaling any row. The element adds only the stabilisation.
 */
template < int NStations = 5 >
class MarmotWarpingStabPressureInterfaceMaterialHypoElasticN {

  static_assert( NStations >= 4,
                 "MarmotWarpingStabPressureInterfaceMaterialHypoElastic requires at least 4 through-thickness "
                 "stations: the 3-point Gauss-Lobatto rule samples xi = {-1,0,1}, where the antisymmetric warping "
                 "profile phi_a = xi - xi^3 vanishes identically, leaving its stiffness block singular." );

public:
  static constexpr int nStations = NStations;

  /** Number of 3x3 tensor slots in the generalized strain: A+, A-, Ws, Wa. */
  static constexpr int nTensorSlots = 4;

  /** ( w, A+, A-, Ws, Wa, p ). */
  static constexpr int nGeneralized = 3 + 9 * nTensorSlots + 1; // 40

  static constexpr int offsetJump                 = 0;
  static constexpr int offsetSurfacePlus          = 3;
  static constexpr int offsetSurfaceMinus         = 12;
  static constexpr int offsetWarpingSymmetric     = 21;
  static constexpr int offsetWarpingAntisymmetric = 30;
  static constexpr int offsetPressure             = 39;

  static constexpr std::array< int, nTensorSlots > slotOffset = { offsetSurfacePlus,
                                                                  offsetSurfaceMinus,
                                                                  offsetWarpingSymmetric,
                                                                  offsetWarpingAntisymmetric };

  using CondensedTangent  = Eigen::Matrix< double, nGeneralized, nGeneralized, Eigen::RowMajor >;
  using GeneralizedVector = Eigen::Matrix< double, nGeneralized, 1 >;

protected:
  const double*                                                         materialProperties;
  const int                                                             nMaterialProperties;
  double                                                                h            = 0.0;
  double                                                                bulkModulus  = 0.0;
  double                                                                shearModulus = 0.0;
  std::vector< double >                                                 baseMaterialProperties;
  std::array< std::unique_ptr< MarmotMaterialHypoElastic >, NStations > stationMaterials;

public:
  const int materialNumber;

  MarmotStateLayoutDynamic stateLayout;

  double characteristicElementLength = 0.0;

  double getInterfaceThickness() const { return h; }
  double getBulkModulus() const { return bulkModulus; }
  double getShearModulus() const { return shearModulus; }

  /** Generalized stress p, laid out as ( f(3), S+(9), S-(9), S_Ws(9), S_Wa(9), R_p(1) ). */
  struct State {
    Eigen::Map< GeneralizedVector > generalizedStress;
    double*                         stateVars;

    State( double* generalizedStress_, double* stateVars_ )
      : generalizedStress( generalizedStress_ ), stateVars( stateVars_ )
    {
    }
  };

  /** Condensed d(p)/d(q), row-major, over the 40-component generalized pair. */
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
    double                                             dPressure; /**< increment of the independent pressure. */

    Deformation( const double* dU_,
                 const double* dSurfaceStrain_,
                 const double* dWarpingStrain_,
                 const double* normal_,
                 const double* separationVector_,
                 double        dPressure_ )
      : dU( dU_ ),
        dSurfaceStrain( dSurfaceStrain_ ),
        dWarpingStrain( dWarpingStrain_ ),
        normal( normal_ ),
        separationVector( separationVector_ ),
        dPressure( dPressure_ )
    {
    }

    /** Coincident-face convenience overload: d = 0, ell = h, d_tau = 0. */
    Deformation( const double* dU_,
                 const double* dSurfaceStrain_,
                 const double* dWarpingStrain_,
                 const double* normal_,
                 double        dPressure_ )
      : Deformation( dU_, dSurfaceStrain_, dWarpingStrain_, normal_, zeroVector(), dPressure_ )
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

  MarmotWarpingStabPressureInterfaceMaterialHypoElasticN( const std::string& materialName,
                                                          const double*      matProperties_,
                                                          int                nMaterialProperties_,
                                                          int                materialNumber_ );

  virtual ~MarmotWarpingStabPressureInterfaceMaterialHypoElasticN() = default;

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

/** Production rule: five Gauss-Lobatto stations, matching the warping material. */
using MarmotWarpingStabPressureInterfaceMaterialHypoElastic = MarmotWarpingStabPressureInterfaceMaterialHypoElasticN<
  5 >;

namespace Marmot::Detail::WarpingStabPressureInterface {

  /** 9x9 deviatoric projector acting on vec(sigma): dev = sigma - (tr sigma / 3) I. */
  inline Eigen::Matrix< double, 9, 9, Eigen::RowMajor > deviatoricProjector()
  {
    using Marmot::Detail::GaussLobattoInterface::flatIndex;

    Eigen::Matrix< double, 9, 9, Eigen::RowMajor >
      projector = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >::Identity();
    for ( int i = 0; i < 3; ++i ) {
      for ( int k = 0; k < 3; ++k ) {
        projector( flatIndex( i, i ), flatIndex( k, k ) ) -= 1.0 / 3.0;
      }
    }
    return projector;
  }

} // namespace Marmot::Detail::WarpingStabPressureInterface

// ---------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------

template < int NStations >
MarmotWarpingStabPressureInterfaceMaterialHypoElasticN<
  NStations >::MarmotWarpingStabPressureInterfaceMaterialHypoElasticN( const std::string& materialName,
                                                                       const double*      matProperties_,
                                                                       int                nMaterialProperties_,
                                                                       int                materialNumber_ )
  : materialProperties( matProperties_ ), nMaterialProperties( nMaterialProperties_ ), materialNumber( materialNumber_ )
{
  using namespace Marmot::Detail::GaussLobattoInterface;

  if ( nMaterialProperties < 3 ) {
    throw std::invalid_argument(
      "MarmotWarpingStabPressureInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h." );
  }

  const double E  = materialProperties[0];
  const double nu = materialProperties[1];

  h = materialProperties[2];
  if ( h <= 0.0 ) {
    throw std::invalid_argument( "MarmotWarpingStabPressureInterfaceMaterialHypoElastic requires h > 0." );
  }
  if ( std::abs( 1.0 - 2.0 * nu ) < 1.0e-12 ) {
    throw std::invalid_argument(
      "MarmotWarpingStabPressureInterfaceMaterialHypoElastic: nu = 0.5 (infinite K) is not supported." );
  }

  bulkModulus  = E / ( 3.0 * ( 1.0 - 2.0 * nu ) );
  shearModulus = E / ( 2.0 * ( 1.0 + nu ) );

  baseMaterialProperties.reserve( nMaterialProperties - 1 );
  baseMaterialProperties.push_back( E );
  baseMaterialProperties.push_back( nu );
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
  stateLayout.add( "pressure", 1 );
  stateLayout.finalize();
}

template < int NStations >
void MarmotWarpingStabPressureInterfaceMaterialHypoElasticN< NStations >::setCharacteristicElementLength(
  double length )
{
  characteristicElementLength = length;
  for ( auto& material : stationMaterials ) {
    if ( material ) {
      material->setCharacteristicElementLength( length );
    }
  }
}

template < int NStations >
void MarmotWarpingStabPressureInterfaceMaterialHypoElasticN< NStations >::computeStress(
  State&               state,
  Tangents&            tangents,
  const Deformation&   deformation,
  const TimeIncrement& timeIncrement )
{
  using namespace Marmot::Detail::GaussLobattoInterface;
  using Marmot::Detail::WarpingStabPressureInterface::deviatoricProjector;

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
      throw std::logic_error( "MarmotWarpingStabPressureInterfaceMaterialHypoElastic has no base material." );
    }
  }

  Vector3d     normal     = deformation.normal;
  const double normalNorm = normal.norm();
  if ( normalNorm <= 1.0e-12 ) {
    throw std::invalid_argument( "MarmotWarpingStabPressureInterfaceMaterialHypoElastic: interface normal is zero." );
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

  const Matrix3dRowMajor ABar = 0.5 * ( APlus + AMinus );
  const Vector3d         gBar = ( w - ABar * dTangential ) / ell;

  const auto lobattoXi     = LobattoRule< NStations >::xi();
  const auto lobattoWeight = LobattoRule< NStations >::weight();

  std::array< std::array< double, nTensorSlots >, NStations > slotCoefficient;
  std::array< Matrix3dRowMajor, NStations >                   AStation;
  std::array< double, NStations >                             lambda;

  /** sum_alpha lambda_alpha c^X_alpha -- ( 1/2, 1/2, 2/3, 0 ) for the exact rules. */
  std::array< double, nTensorSlots > meanSlotCoefficient{};

  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    const double xi = lobattoXi[alpha];

    lambda[alpha] = 0.5 * lobattoWeight[alpha];

    slotCoefficient[alpha][0] = 0.5 * ( 1.0 + xi ); // Nplus  -> A+
    slotCoefficient[alpha][1] = 0.5 * ( 1.0 - xi ); // Nminus -> A-
    slotCoefficient[alpha][2] = warpingProfileSymmetric( xi );
    slotCoefficient[alpha][3] = warpingProfileAntisymmetric( xi );

    AStation[alpha] = slotCoefficient[alpha][0] * APlus + slotCoefficient[alpha][1] * AMinus +
                      slotCoefficient[alpha][2] * WarpSymmetric + slotCoefficient[alpha][3] * WarpAntisymmetric;

    for ( int slot = 0; slot < nTensorSlots; ++slot ) {
      meanSlotCoefficient[slot] += lambda[alpha] * slotCoefficient[alpha][slot];
    }
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
  double* pressurePtr       = stateLayout.getPtr( state.stateVars, "pressure" );

  const double pressureTotal = pressurePtr[0] + deformation.dPressure;

  const MarmotMaterialHypoElastic::timeInfo timeInfo{ timeIncrement.timeOld, timeIncrement.dT };

  const Matrix3x9RowMajor Rt          = stressToForce( normal, 1.0, 1.0 );
  const Matrix9x3RowMajor Bn          = jumpToGradient( normal, 1.0 );
  const Matrix9dRowMajor  Pdev        = deviatoricProjector();
  const Vector9d          identityVec = Eigen::Map< const Vector9d >( Matrix3dRowMajor::Identity().eval().data() );

  const double Qc = std::max( 1.0, std::abs( materialProperties[0] ) );

  // ---- local problem y = ( g^(1..N), t ) ----
  // Equilibrium is imposed on the DEVIATORIC traction: sigma~^(alpha) n = t~ with
  // sigma~ = Pdev sigma - p I and p constant through the thickness is identical
  // to Pdev sigma^(alpha) n = t, with t~ = t - p n. The local solve is therefore
  // exactly p-independent -- the multi-station form of the two-station
  // material's d gamma / dp = 0.
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

      const Vector9d deviatoric = Pdev * Eigen::Map< const Vector9d >( evaluation.station[alpha].stress.data() );

      evaluation.tAlpha[alpha]                        = Rt * deviatoric;
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
        const Matrix3dRowMajor Qalpha                        = Rt * ( Pdev * current.station[alpha].CFull ) * Bn;
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
    throw Marmot::StressUpdateFailed( "MarmotWarpingStabPressureInterfaceMaterialHypoElastic: local "
                                      "traction-equilibrium solve did not converge." );
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
  pressurePtr[0]    = pressureTotal;

  // ---- condensation of the local unknowns onto the 40-component pair ----
  const Matrix3x9RowMajor Mdtau = tangentialContraction( dTangential );

  MatrixLocalRowMajor Ry = MatrixLocalRowMajor::Zero();
  MatrixLocalxGen     Rq = MatrixLocalxGen::Zero();

  std::array< Matrix9dRowMajor, NStations >  CDeviatoric;
  std::array< Matrix3x9RowMajor, NStations > RtCDeviatoric;

  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    CDeviatoric[alpha]   = Pdev * current.station[alpha].CFull;
    RtCDeviatoric[alpha] = Rt * CDeviatoric[alpha];

    const Matrix3dRowMajor Qalpha                         = RtCDeviatoric[alpha] * Bn;
    Ry.template block< 3, 3 >( 3 * alpha, 3 * alpha )     = Qalpha;
    Ry.template block< 3, 3 >( 3 * alpha, 3 * NStations ) = -Matrix3dRowMajor::Identity();
    Ry.template block< 3, 3 >( 3 * NStations, 3 * alpha ) = lambda[alpha] * Matrix3dRowMajor::Identity();

    for ( int slot = 0; slot < nTensorSlots; ++slot ) {
      Rq.template block< 3, 9 >( 3 * alpha, slotOffset[slot] ) = slotCoefficient[alpha][slot] * RtCDeviatoric[alpha];
    }
  }

  // The pressure column of Rq stays exactly zero: the local problem does not see p.
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
                 CDeviatoric[alpha];
      }
      pq.template block< 9, 9 >( slotOffset[rowSlot], slotOffset[colSlot] ) = block;
    }
  }

  py.template block< 3, 3 >( offsetJump, 3 * NStations ) = ( h / ell ) * Matrix3dRowMajor::Identity();
  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    const Matrix9x3RowMajor CB = CDeviatoric[alpha] * Bn;
    for ( int slot = 0; slot < nTensorSlots; ++slot ) {
      py.template block< 9, 3 >( slotOffset[slot], 3 * alpha ) = ( h * lambda[alpha] * slotCoefficient[alpha][slot] ) *
                                                                 CB;
    }
  }
  py.template block< 9, 3 >( offsetSurfacePlus, 3 * NStations )  = -( 0.5 * h / ell ) * Mdtau.transpose();
  py.template block< 9, 3 >( offsetSurfaceMinus, 3 * NStations ) = -( 0.5 * h / ell ) * Mdtau.transpose();

  const MatrixLocalxGen RyInvRq = Ry.fullPivLu().solve( Rq );

  tangents.condensed = pq - py * RyInvRq;

  // ---- mixed pressure blocks: purely geometric, since dy/dp = 0 ----
  // f   = (h/ell) ( t - p n )
  // S_X = h sum lambda c^X ( Pdev sigma - p I ) - [face] (h/(2 ell)) Mdtau^T ( t - p n )
  tangents.condensed.template block< 3, 1 >( offsetJump, offsetPressure ) = -( h / ell ) * normal;
  for ( int slot = 0; slot < nTensorSlots; ++slot ) {
    Vector9d column = -( h * meanSlotCoefficient[slot] ) * identityVec;
    if ( slotOffset[slot] == offsetSurfacePlus || slotOffset[slot] == offsetSurfaceMinus ) {
      column += ( 0.5 * h / ell ) * ( Mdtau.transpose() * normal );
    }
    tangents.condensed.template block< 9, 1 >( slotOffset[slot], offsetPressure ) = column;
  }

  // R_p = h [ sum lambda tr(d eps^(alpha)) + dp/K ]
  //     = h [ sum_slot meanCoefficient_slot tr(X_slot) + gbar . n + dp/K ]
  // The station gradients enter only through sum lambda g = gbar, which is a
  // constraint of the local problem, so this row is geometric: dR_p/dy = 0.
  tangents.condensed.template block< 1, 3 >( offsetPressure, offsetJump ) = ( h / ell ) * normal.transpose();
  for ( int slot = 0; slot < nTensorSlots; ++slot ) {
    Eigen::Matrix< double, 1, 9 > row = ( h * meanSlotCoefficient[slot] ) * identityVec.transpose();
    if ( slotOffset[slot] == offsetSurfacePlus || slotOffset[slot] == offsetSurfaceMinus ) {
      row -= ( 0.5 * h / ell ) * ( normal.transpose() * Mdtau );
    }
    tangents.condensed.template block< 1, 9 >( offsetPressure, slotOffset[slot] ) = row;
  }
  tangents.condensed( offsetPressure, offsetPressure ) = h / bulkModulus;

  // ---- generalized stresses ----
  const Vector3d tShared   = yCurrent.template segment< 3 >( 3 * NStations );
  const Vector3d tModified = tShared - pressureTotal * normal;

  state.generalizedStress.setZero();
  state.generalizedStress.template segment< 3 >( offsetJump ) = ( h / ell ) * tModified;

  for ( int slot = 0; slot < nTensorSlots; ++slot ) {
    Vector9d sum = Vector9d::Zero();
    for ( int alpha = 0; alpha < NStations; ++alpha ) {
      const Vector9d deviatoricStress = Pdev * Eigen::Map< const Vector9d >( current.station[alpha].stress.data() );
      sum += ( h * lambda[alpha] * slotCoefficient[alpha][slot] ) * deviatoricStress;
    }
    sum -= ( h * meanSlotCoefficient[slot] * pressureTotal ) * identityVec;
    if ( slotOffset[slot] == offsetSurfacePlus || slotOffset[slot] == offsetSurfaceMinus ) {
      sum -= ( 0.5 * h / ell ) * ( Mdtau.transpose() * tModified );
    }
    state.generalizedStress.template segment< 9 >( slotOffset[slot] ) = sum;
  }

  double volumetricStrainIncrement = gBar.dot( normal );
  volumetricStrainIncrement += meanSlotCoefficient[0] * APlus.trace();
  volumetricStrainIncrement += meanSlotCoefficient[1] * AMinus.trace();
  volumetricStrainIncrement += meanSlotCoefficient[2] * WarpSymmetric.trace();
  volumetricStrainIncrement += meanSlotCoefficient[3] * WarpAntisymmetric.trace();

  state.generalizedStress( offsetPressure ) = h * ( volumetricStrainIncrement + deformation.dPressure / bulkModulus );
}

template < int NStations >
void MarmotWarpingStabPressureInterfaceMaterialHypoElasticN< NStations >::initializeYourself( double* stateVars,
                                                                                              int     nStateVars )
{
  using namespace Marmot::Detail::GaussLobattoInterface;

  for ( int i = 0; i < nStateVars; ++i ) {
    stateVars[i] = 0.0;
  }

  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    if ( !stationMaterials[alpha] ) {
      continue;
    }
    stationMaterials[alpha]->initializeYourself( stateLayout.getPtr( stateVars, stationStateVarsName( alpha ) ),
                                                 stationMaterials[alpha]->getNumberOfRequiredStateVars() );
  }
}

template < int NStations >
double MarmotWarpingStabPressureInterfaceMaterialHypoElasticN< NStations >::getDensity()
{
  if ( !stationMaterials[0] ) {
    return -1;
  }

  return stationMaterials[0]->getDensity( nullptr );
}
