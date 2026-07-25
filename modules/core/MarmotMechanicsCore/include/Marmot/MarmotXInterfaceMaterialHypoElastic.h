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

#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"

#include <memory>
#include <string>
#include <vector>

/**
 * Hypoelastic interface material resolving the top and bottom surface
 * strains of the interface layer independently (a "+ / -" split), instead
 * of collapsing them into one averaged surface gradient.
 *
 * The layer's virtual work is evaluated as the thin-interface weak form
 *
 *   h \int_I < C_{ijkl} u_{k,l} \hat u_{i,j} > dA,
 *   <.> := 1/2 [ (.)^+ + (.)^- ],
 *
 * i.e. the average of two SEPARATELY contracted one-sided energies, each
 * built from its own face's surface gradient A^+ / A^- plus the shared
 * normal-jump correction:
 *
 *   q^s   = [u] - A^s d_tau,
 *   G^s   = A^s + (1/ell) q^s \otimes n,
 *   eps^s = sym(G^s),        s in {+, -}.
 *
 * Each eps^s drives its own embedded bulk-material instance (own
 * plastic/history state), each integrated over half the constitutive
 * thickness h/2. Unlike the single-averaged-gradient formulation, a
 * cross-sectional-rotation mode (A^+ = -A^-) does not produce zero strain
 * on both sides simultaneously, so the layer is not blind to it.
 *
 * The interface kinematics use the actual connector between paired mesh
 * points,
 *
 *   d = x_top - x_bottom = ell n + d_tau,
 *
 * exactly as in MarmotCorrectedInterfaceMaterialHypoElastic. For coincident
 * interface faces, ell falls back to the constitutive interface thickness h
 * and d_tau is set to zero.
 */
class MarmotXInterfaceMaterialHypoElastic {

protected:
  const double*                                materialProperties;
  const int                                    nMaterialProperties;
  double                                       h = 0.0;
  std::vector< double >                        baseMaterialProperties;
  std::unique_ptr< MarmotMaterialHypoElastic > topMaterial;
  std::unique_ptr< MarmotMaterialHypoElastic > bottomMaterial;

public:
  using TensorMap3d    = Marmot::FastorStandardTensors::TensorMap3d;
  using TensorMap33d   = Marmot::FastorStandardTensors::TensorMap33d;
  using TensorMap333d  = Marmot::FastorStandardTensors::TensorMap333d;
  using TensorMap3333d = Marmot::FastorStandardTensors::TensorMap3333d;
  using TensorMap6d    = Marmot::FastorStandardTensors::TensorMap6d;
  using TensorMap18d   = Marmot::FastorStandardTensors::TensorMap18d;

  const int materialNumber;

  MarmotXInterfaceMaterialHypoElastic( const std::string& materialName,
                                       const double*      matProperties_,
                                       int                nMaterialProperties_,
                                       int                materialNumber_ );

  virtual ~MarmotXInterfaceMaterialHypoElastic() = default;

  MarmotStateLayoutDynamic stateLayout;

  double characteristicElementLength;

  void setCharacteristicElementLength( double length );

  struct State {
    /** Generalized force conjugate to [u], contribution from side +. */
    TensorMap3d forcePlus;

    /** Generalized force conjugate to [u], contribution from side -. */
    TensorMap3d forceMinus;

    /** Generalized surface resultant conjugate to A^+. */
    TensorMap33d surfaceStressPlus;

    /** Generalized surface resultant conjugate to A^-. */
    TensorMap33d surfaceStressMinus;

    double* stateVars;
  };

  struct Tangents {
    /** d(forcePlus_i) / d([u_k]). */
    TensorMap33d Q_plus;

    /** d(forceMinus_i) / d([u_k]). */
    TensorMap33d Q_minus;

    /** d(forcePlus_i) / d(A^+_{kl}). */
    TensorMap333d H_plus;

    /** d(forceMinus_i) / d(A^-_{kl}). */
    TensorMap333d H_minus;

    /** d(surfaceStressPlus_ij) / d([u_k]). */
    TensorMap333d K_plus;

    /** d(surfaceStressMinus_ij) / d([u_k]). */
    TensorMap333d K_minus;

    /** d(surfaceStressPlus_ij) / d(A^+_{kl}). */
    TensorMap3333d Z_plus;

    /** d(surfaceStressMinus_ij) / d(A^-_{kl}). */
    TensorMap3333d Z_minus;
  };

  struct Deformation {
    TensorMap6d  dU;
    TensorMap18d dSurfaceStrain;
    TensorMap3d  normal;
    TensorMap3d  separationVector;

    // Fastor's const TensorMap cannot be used with slicing and norm operations.
    // These views are therefore mutable types but are exposed through const Deformation&.
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
