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
#include "Marmot/MarmotFastorTensorBasics.h"
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

namespace Marmot::Detail::GaussLobattoInterface {

  /**
   * Fixed Gauss-Lobatto rules on xi in [-1,1] for the station counts
   * supported by the through-thickness convergence study (Section 13):
   * n=3 (insufficient for plastic bending -- its only interior point sits
   * at xi=0, zero lever arm), n=4 (minimum useful rule), n=5 (preferred
   * production rule), n=7 (convergence check). NOT exposed as a user
   * input: the station count is a compile-time template parameter, chosen
   * internally.
   */
  template < int NStations >
  struct LobattoRule;

  template <>
  struct LobattoRule< 3 > {
    static std::array< double, 3 > xi() { return { -1.0, 0.0, 1.0 }; }
    static std::array< double, 3 > weight() { return { 1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0 }; }
  };

  template <>
  struct LobattoRule< 4 > {
    static std::array< double, 4 > xi()
    {
      const double a = 1.0 / std::sqrt( 5.0 );
      return { -1.0, -a, a, 1.0 };
    }
    static std::array< double, 4 > weight() { return { 1.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0 }; }
  };

  template <>
  struct LobattoRule< 5 > {
    static std::array< double, 5 > xi()
    {
      const double a = std::sqrt( 3.0 / 7.0 );
      return { -1.0, -a, 0.0, a, 1.0 };
    }
    static std::array< double, 5 > weight()
    {
      return { 1.0 / 10.0, 49.0 / 90.0, 32.0 / 45.0, 49.0 / 90.0, 1.0 / 10.0 };
    }
  };

  template <>
  struct LobattoRule< 7 > {
    static std::array< double, 7 > xi()
    {
      return { -1.0, -0.8302238962785670, -0.4688487934707142, 0.0, 0.4688487934707142, 0.8302238962785670, 1.0 };
    }
    static std::array< double, 7 > weight()
    {
      return { 0.0476190476190476,
               0.2768260473615659,
               0.4317453812098627,
               0.4876190476190476,
               0.4317453812098627,
               0.2768260473615659,
               0.0476190476190476 };
    }
  };

} // namespace Marmot::Detail::GaussLobattoInterface

/**
 * Multi-station (Gauss-Lobatto) through-thickness interface material.
 *
 * Generalizes MarmotEquilibratedXInterfaceMaterialHypoElastic's two-sided
 * (+/-) construction to a fixed NStations-point Gauss-Lobatto rule through
 * the constitutive thickness h. The external generalized kinematics/statics
 * and the full nine-block tangent structure are IDENTICAL to the
 * equilibrated two-sided material -- q = (w, A+, A-), p = (f, S+, S-) --
 * so this material is a drop-in replacement wherever
 * MarmotEquilibratedXInterfaceMaterialHypoElastic is used (in particular by
 * GaussLobattoInterfaceFiniteElement, cloned from YInterfaceFiniteElement
 * with only the material type changed).
 *
 * At station alpha (xi_alpha in [-1,1], normalized weight lambda_alpha,
 * sum lambda_alpha = 1) the tangential surface gradient is LINEARLY
 * interpolated (no new kinematics, no bubble):
 *
 *   A^(alpha) = Nminus_alpha A- + Nplus_alpha A+,
 *   Nminus_alpha = (1 - xi_alpha)/2,  Nplus_alpha = (1 + xi_alpha)/2.
 *
 * Each station carries its own independent embedded bulk-material instance
 * and its own persistent history, and its own internal normal-gradient
 * unknown g^(alpha), subject to the single weighted-compatibility
 * constraint
 *
 *   sum_alpha lambda_alpha g^(alpha) = gbar(q,d,n),
 *
 * where gbar is the SAME corrected average normal gradient used by the
 * equilibrated two-sided material (including the tangential-separation
 * d_tau correction). A shared (Lagrange-multiplier) traction t enforces
 * equal station tractions t^(alpha) = t for all alpha. The local unknown
 * is y = (g^(1),...,g^(NStations), t), 3*NStations+3 scalars; after local
 * Newton convergence the generalized force p and the full condensed 21x21
 * tangent are obtained by static condensation of y.
 *
 * NStations is a template parameter (compile-time), NOT a runtime material
 * property -- per Section 13, the station count is not user-exposed. The
 * default alias `MarmotGaussLobattoInterfaceMaterialHypoElastic` below is
 * the fixed five-point production rule.
 */
template < int NStations >
class MarmotGaussLobattoInterfaceMaterialHypoElasticN {

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

  MarmotGaussLobattoInterfaceMaterialHypoElasticN( const std::string& materialName,
                                                   const double*      matProperties_,
                                                   int                nMaterialProperties_,
                                                   int                materialNumber_ );

  virtual ~MarmotGaussLobattoInterfaceMaterialHypoElasticN() = default;

  MarmotStateLayoutDynamic stateLayout;

  double characteristicElementLength;

  void setCharacteristicElementLength( double length );

  struct State {
    /** Equilibrated generalized force conjugate to the raw mesh jump w=[u]. */
    TensorMap3d generalizedForce;

    /** Generalized surface resultant conjugate to A+ (NOT the physical
     * stress at the top station -- see the "stationStress0"/"stationStressN-1"
     * state views for that). */
    TensorMap33d surfaceStressPlus;

    /** Generalized surface resultant conjugate to A- (NOT the physical
     * stress at the bottom station). */
    TensorMap33d surfaceStressMinus;

    double* stateVars;
  };

  /**
   * Full 3x3 block tangent over the external generalized strain
   * x = (w, A+, A-) and generalized stress p = (f, S+, S-), after static
   * condensation of the internal per-station normal gradients and the
   * shared traction. No major symmetry is assumed.
   */
  struct Tangents {
    TensorMap33d   Q_ww;   /**< d(f)/d(w).    3x3 */
    TensorMap333d  Q_wAp;  /**< d(f)/d(A+).   3x9 */
    TensorMap333d  Q_wAm;  /**< d(f)/d(A-).   3x9 */
    TensorMap333d  Q_Apw;  /**< d(S+)/d(w).   9x3 */
    TensorMap3333d Q_ApAp; /**< d(S+)/d(A+).  9x9 */
    TensorMap3333d Q_ApAm; /**< d(S+)/d(A-).  9x9 */
    TensorMap333d  Q_Amw;  /**< d(S-)/d(w).   9x3 */
    TensorMap3333d Q_AmAp; /**< d(S-)/d(A+).  9x9 */
    TensorMap3333d Q_AmAm; /**< d(S-)/d(A-).  9x9 */
  };

  struct Deformation {
    TensorMap6d  dU;
    TensorMap18d dSurfaceStrain;
    TensorMap3d  normal;
    TensorMap3d  separationVector;

    Deformation( const double* dU_,
                 const double* dSurfaceStrain_,
                 const double* normal_,
                 const double* separationVector_ )
      : dU( const_cast< double* >( dU_ ) ),
        dSurfaceStrain( const_cast< double* >( dSurfaceStrain_ ) ),
        normal( const_cast< double* >( normal_ ) ),
        separationVector( const_cast< double* >( separationVector_ ) )
    {
    }

    /**
     * Backward-compatible constructor for coincident interface faces.
     * The material then uses d = 0, ell = h, and d_tau = 0.
     */
    Deformation( const double* dU_, const double* dSurfaceStrain_, const double* normal_ )
      : Deformation( dU_, dSurfaceStrain_, normal_, zeroSeparationVector() )
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

/** Preferred production rule (Section 3): fixed five-point Gauss-Lobatto. */
using MarmotGaussLobattoInterfaceMaterialHypoElastic = MarmotGaussLobattoInterfaceMaterialHypoElasticN< 5 >;

// ---------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------

namespace Marmot::Detail::GaussLobattoInterface {

  using Vector3d = Eigen::Matrix< double, 3, 1 >;

  constexpr int flatIndex( int i, int j )
  {
    return 3 * i + j;
  }

  constexpr double localNewtonTolerance     = 1.0e-10;
  constexpr int    maxLocalNewtonIterations = 60;

  /** Identical construction to the corresponding helper in
   * MarmotEquilibratedXInterfaceMaterialHypoElastic.cpp. */
  inline Eigen::Matrix< double, 9, 9, Eigen::RowMajor > fullGradientTangent( const Marmot::Matrix6d& tangentVoigt )
  {
    using Matrix3dRowMajor = Eigen::Matrix< double, 3, 3, Eigen::RowMajor >;
    using Matrix9dRowMajor = Eigen::Matrix< double, 9, 9, Eigen::RowMajor >;

    Matrix9dRowMajor tangentFull = Matrix9dRowMajor::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        Matrix3dRowMajor dGradient = Matrix3dRowMajor::Zero();
        dGradient( k, l )          = 1.0;

        const Matrix3dRowMajor dStrain      = 0.5 * ( dGradient + dGradient.transpose() );
        const Marmot::Vector6d dStrainVoigt = Marmot::ContinuumMechanics::VoigtNotation::strainToVoigt( dStrain );
        const Marmot::Vector6d dStressVoigt = tangentVoigt * dStrainVoigt;
        const Eigen::Matrix3d  dStress      = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress( dStressVoigt );

        const int column = flatIndex( k, l );
        for ( int i = 0; i < 3; ++i ) {
          for ( int j = 0; j < 3; ++j ) {
            tangentFull( flatIndex( i, j ), column ) = dStress( i, j );
          }
        }
      }
    }

    return tangentFull;
  }

  struct InterfaceGeometry {
    double   normalSeparation;
    Vector3d tangentialSeparation;
  };

  inline InterfaceGeometry evaluateInterfaceGeometry( const Vector3d& normal,
                                                      const Vector3d& separationVector,
                                                      double          constitutiveThickness )
  {
    constexpr double tolerance = 1.0e-12;

    if ( constitutiveThickness <= 0.0 ) {
      throw std::invalid_argument(
        "MarmotGaussLobattoInterfaceMaterialHypoElastic: interface thickness h must be positive." );
    }

    if ( separationVector.norm() <= tolerance ) {
      return { constitutiveThickness, Vector3d::Zero() };
    }

    const double normalSeparation = separationVector.dot( normal );
    if ( normalSeparation <= tolerance ) {
      throw std::invalid_argument( "MarmotGaussLobattoInterfaceMaterialHypoElastic: the top-bottom connector must "
                                   "have a positive normal component." );
    }

    return { normalSeparation, separationVector - normalSeparation * normal };
  }

  /** Maps a 3-vector g to vec(g \otimes n) / normalSeparation. */
  inline Eigen::Matrix< double, 9, 3, Eigen::RowMajor > jumpToGradient( const Vector3d& normal,
                                                                        double          normalSeparation )
  {
    Eigen::Matrix< double, 9, 3, Eigen::RowMajor > map = Eigen::Matrix< double, 9, 3, Eigen::RowMajor >::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int l = 0; l < 3; ++l ) {
        map( flatIndex( k, l ), k ) = normal( l ) / normalSeparation;
      }
    }

    return map;
  }

  /** Maps vec(sigma) to (sideThickness/normalSeparation) * sigma n. */
  inline Eigen::Matrix< double, 3, 9, Eigen::RowMajor > stressToForce( const Vector3d& normal,
                                                                       double          sideThickness,
                                                                       double          normalSeparation )
  {
    Eigen::Matrix< double, 3, 9, Eigen::RowMajor > map   = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >::Zero();
    const double                                   scale = sideThickness / normalSeparation;

    for ( int i = 0; i < 3; ++i ) {
      for ( int b = 0; b < 3; ++b ) {
        map( i, flatIndex( i, b ) ) = scale * normal( b );
      }
    }

    return map;
  }

  /** Maps a 9-vector dA direction to d(A . tangentialSeparation), a 3x9
   * matrix: M(k,(a,b)) = delta(k,a) * tangentialSeparation(b). */
  inline Eigen::Matrix< double, 3, 9, Eigen::RowMajor > tangentialContraction( const Vector3d& tangentialSeparation )
  {
    Eigen::Matrix< double, 3, 9, Eigen::RowMajor > map = Eigen::Matrix< double, 3, 9, Eigen::RowMajor >::Zero();

    for ( int k = 0; k < 3; ++k ) {
      for ( int b = 0; b < 3; ++b ) {
        map( k, flatIndex( k, b ) ) = tangentialSeparation( b );
      }
    }

    return map;
  }

  /** One station's trial constitutive response, evaluated from a scratch
   * copy of the OLD (committed) state -- never mutates the persistent
   * state vars. */
  struct StationTrial {
    Eigen::Matrix< double, 3, 3, Eigen::RowMajor > stress;
    Eigen::Matrix< double, 9, 9, Eigen::RowMajor > CFull;
    std::vector< double >                          trialStateVars;
  };

  inline StationTrial evaluateStationTrial( MarmotMaterialHypoElastic&                            material,
                                            const double*                                         oldStateVars,
                                            int                                                   nStateVars,
                                            const Eigen::Matrix< double, 3, 3, Eigen::RowMajor >& stressCurrent,
                                            const Marmot::Vector6d&                               strainIncrementVoigt,
                                            const MarmotMaterialHypoElastic::timeInfo&            timeInfo )
  {
    StationTrial trial;
    trial.trialStateVars.assign( oldStateVars, oldStateVars + nStateVars );

    const Eigen::Matrix3d  stressCurrentSym( 0.5 * ( stressCurrent + stressCurrent.transpose() ) );
    const Marmot::Vector6d stressVoigt = Marmot::ContinuumMechanics::VoigtNotation::stressToVoigt( stressCurrentSym );

    Marmot::Matrix6d                   tangentVoigt = Marmot::Matrix6d::Zero();
    MarmotMaterialHypoElastic::state3D baseState{ stressVoigt, 0.0, 0.0, trial.trialStateVars.data() };

    material.computeStress( baseState, tangentVoigt, strainIncrementVoigt, timeInfo );

    const Eigen::Matrix3d stressUpdatedEigen = Marmot::ContinuumMechanics::VoigtNotation::voigtToStress(
      baseState.stress );
    trial.stress = stressUpdatedEigen;
    trial.CFull  = fullGradientTangent( tangentVoigt );

    return trial;
  }

  inline std::string stationStateVarsName( int alpha )
  {
    return "stationMaterialStateVars" + std::to_string( alpha );
  }

  inline std::string normalGradientName( int alpha )
  {
    return "normalGradient" + std::to_string( alpha );
  }

  inline std::string stationStressName( int alpha )
  {
    return "stationStress" + std::to_string( alpha );
  }

} // namespace Marmot::Detail::GaussLobattoInterface

template < int NStations >
MarmotGaussLobattoInterfaceMaterialHypoElasticN< NStations >::MarmotGaussLobattoInterfaceMaterialHypoElasticN(
  const std::string& materialName,
  const double*      matProperties_,
  int                nMaterialProperties_,
  int                materialNumber_ )
  : materialProperties( matProperties_ ), nMaterialProperties( nMaterialProperties_ ), materialNumber( materialNumber_ )
{
  using namespace Marmot::Detail::GaussLobattoInterface;

  if ( nMaterialProperties < 3 ) {
    throw std::invalid_argument(
      "MarmotGaussLobattoInterfaceMaterialHypoElastic requires at least E, nu, and interface thickness h." );
  }

  h = materialProperties[2];
  if ( h <= 0.0 ) {
    throw std::invalid_argument( "MarmotGaussLobattoInterfaceMaterialHypoElastic requires h > 0." );
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
void MarmotGaussLobattoInterfaceMaterialHypoElasticN< NStations >::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
  for ( auto& material : stationMaterials ) {
    if ( material ) {
      material->setCharacteristicElementLength( length );
    }
  }
}

template < int NStations >
void MarmotGaussLobattoInterfaceMaterialHypoElasticN< NStations >::computeStress( State&               state,
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
      throw std::logic_error( "MarmotGaussLobattoInterfaceMaterialHypoElastic has no base material." );
    }
  }

  const Eigen::Map< const Vector3d > normalMap( deformation.normal.data() );
  const Eigen::Map< const Vector3d > separationMap( deformation.separationVector.data() );

  Vector3d     normal     = normalMap;
  const double normalNorm = normal.norm();
  if ( normalNorm <= 1.0e-12 ) {
    throw std::invalid_argument( "MarmotGaussLobattoInterfaceMaterialHypoElastic: interface normal is zero." );
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

  std::array< Matrix3dRowMajor, NStations > AStation;
  std::array< double, NStations >           lambda, Nplus, Nminus;
  for ( int alpha = 0; alpha < NStations; ++alpha ) {
    lambda[alpha]   = 0.5 * lobattoWeight[alpha];
    Nplus[alpha]    = 0.5 * ( 1.0 + lobattoXi[alpha] );
    Nminus[alpha]   = 0.5 * ( 1.0 - lobattoXi[alpha] );
    AStation[alpha] = Nminus[alpha] * AMinus + Nplus[alpha] * APlus;
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

  // Attempt the local Newton solve from a given initial guess for y. Tried
  // first from the warm-started (persisted) guess; if that attempt gets
  // stuck (e.g. a bad warm start pointing away from the solution near a
  // yield-surface kink -- the raw Newton direction can have no improving
  // step at ANY line-search fraction, so the outer loop would otherwise
  // burn all iterations making no progress), it is retried from a neutral
  // "fresh" guess. This is a standard multi-start robustness technique.
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
        // No improving step along the Newton direction at any tested
        // fraction: further iterating from here is pointless. Report
        // non-convergence immediately so the caller can retry from a
        // different starting point instead of burning the remaining
        // iteration budget in place.
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
    // Fresh, neutral restart: every station gets the shared compatible
    // normal gradient gBar (trivially satisfying the weighted-compatibility
    // constraint) and a zero traction guess.
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
      "MarmotGaussLobattoInterfaceMaterialHypoElastic: local traction-equilibrium solve did not converge." );
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
void MarmotGaussLobattoInterfaceMaterialHypoElasticN< NStations >::initializeYourself( double* stateVars,
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
double MarmotGaussLobattoInterfaceMaterialHypoElasticN< NStations >::getDensity()
{
  if ( !stationMaterials[0] ) {
    return -1;
  }

  return stationMaterials[0]->getDensity( nullptr );
}
