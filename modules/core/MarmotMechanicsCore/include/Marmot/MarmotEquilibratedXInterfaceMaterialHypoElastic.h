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
 * Equilibrated two-sided hypoelastic interface material.
 *
 * Like MarmotXInterfaceMaterialHypoElastic, the top and bottom surface
 * gradients A+/A- are kept separate and drive independent embedded
 * bulk-material instances. Unlike that material, this one does NOT impose
 * the same normal-gradient contribution g = [u]/ell on both sides. Instead
 * it introduces the normal-gradient jump
 *
 *   z = g+ - g-
 *
 * as an internal variable, condensed at every stress update by a local
 * Newton iteration enforcing physical traction equilibrium
 *
 *   sigma+ n = sigma- n.
 *
 * The average surface gradient
 *
 *   Abar = (A+ + A-)/2
 *
 * together with the connector geometry provides the shared normal-gradient
 * average
 *
 *   gbar = (1/ell) ( w - Abar d_tau ),    w = [u],
 *
 * and
 *
 *   g+ = gbar + z/2,   g- = gbar - z/2.
 *
 * The local Jacobian of the equilibrium residual r(z) = t+(z) - t-(z) is the
 * AVERAGE one-sided acoustic tensor
 *
 *   <Q> = 1/2 ( n C+ n + n C- n ),
 *
 * which remains invertible even where Q+ or Q- individually degenerate
 * (e.g. at yield), so only <Q> is ever inverted, never Q+ or Q- separately.
 *
 * After local convergence, the generalized stress and the FULLY CONDENSED
 * tangent (including the generally nonzero A+/A- cross-coupling blocks
 * introduced by the equilibrium constraint) are returned to the element.
 */
class MarmotEquilibratedXInterfaceMaterialHypoElastic {

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

  MarmotEquilibratedXInterfaceMaterialHypoElastic( const std::string& materialName,
                                                   const double*      matProperties_,
                                                   int                nMaterialProperties_,
                                                   int                materialNumber_ );

  virtual ~MarmotEquilibratedXInterfaceMaterialHypoElastic() = default;

  MarmotStateLayoutDynamic stateLayout;

  double characteristicElementLength;

  void setCharacteristicElementLength( double length );

  struct State {
    /** Equilibrated generalized force conjugate to the raw mesh jump w=[u]. */
    TensorMap3d generalizedForce;

    /** Generalized surface resultant conjugate to A+. */
    TensorMap33d surfaceStressPlus;

    /** Generalized surface resultant conjugate to A-. */
    TensorMap33d surfaceStressMinus;

    double* stateVars;
  };

  /**
   * Full 3x3 block tangent over the external generalized strain
   * x = (w, A+, A-) and generalized stress p = (f, S+, S-), after static
   * condensation of the internal normal-gradient jump z. No major symmetry
   * is assumed and, in general, the cross blocks are nonzero:
   * d(S+)/d(A-) != 0, d(S-)/d(A+) != 0.
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
