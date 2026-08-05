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

#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotStateHelpers.h"

#include <memory>
#include <string>
#include <vector>

/**
 * Partial-mixed interface constitutive kernel in the (Abar, DeltaA, g)
 * parameterization, for the partial Hu-Washizu interface element with a
 * nodal common-traction field.
 *
 * Unlike MarmotEquilibratedXInterfaceMaterialHypoElastic (which derives the
 * average normal gradient gbar from the displacement jump w = [u]/h), this
 * kernel takes the average normal gradient g as an INDEPENDENT input and
 * keeps the tangential average/jump gradients (Abar, DeltaA) as its other
 * inputs. The local normal-gradient jump gamma is solved from traction
 * equilibrium exactly as before.
 *
 * Per-side gradient reconstruction (spec section 1, d_tau = 0):
 *   G^pm = Abar +/- (1/2) DeltaA + ( g +/- (1/2) gamma ) (x) n .
 * Interface energy per unit midsurface area is h * phi, with
 *   phi = 1/2 ( psi^+(G^+) + psi^-(G^-) ).
 *
 * Local solve (spec section 3): gamma from
 *   F_gamma = ( sigma^+ - sigma^- ) n = 0,
 * Newton matrix <Q> = 1/2 ( Q^+ + Q^- ), Q^pm_ik = C^pm_ijkl n_j n_l.
 *
 * Reduced outputs (spec sections 2, 4), with <.> = 1/2( (.)^+ + (.)^- ),
 * [.] = (.)^+ - (.)^-:
 *   s_Abar   = <sigma>          (9, tangential 2nd index used by the element)
 *   s_DeltaA = 1/4 [sigma]      (9)
 *   t^mat    = <sigma> n        (3, common constitutive traction)
 * The reduced constitutive tangent over x = (Abar, DeltaA, g) is condensed
 * consistently (spec section 9):
 *   H^red = y_,x - y_,gamma (F_gamma,gamma)^{-1} F_gamma,x ,   F_gamma,gamma = <Q>,
 * and returned as the blocks H_aa (18x18), H_ag (18x3), H_ga (3x18),
 * H_gg (3x3), with a = (Abar, DeltaA).
 *
 * State: two embedded bulk-material state blocks (independent top/bottom
 * plastic histories), the committed normal-gradient jump gamma (warm start),
 * and the committed one-sided Cauchy stresses sigma^+/sigma^- (needed as the
 * hypoelastic increment base).
 */
class MarmotPartialMixedInterfaceMaterialHypoElastic {

protected:
  const double*                                materialProperties;
  const int                                    nMaterialProperties;
  double                                       h = 0.0;
  std::vector< double >                        baseMaterialProperties;
  std::unique_ptr< MarmotMaterialHypoElastic > topMaterial;
  std::unique_ptr< MarmotMaterialHypoElastic > bottomMaterial;

public:
  const int materialNumber;

  MarmotPartialMixedInterfaceMaterialHypoElastic( const std::string& materialName,
                                                  const double*      matProperties_,
                                                  int                nMaterialProperties_,
                                                  int                materialNumber_ );

  virtual ~MarmotPartialMixedInterfaceMaterialHypoElastic() = default;

  MarmotStateLayoutDynamic stateLayout;

  double characteristicElementLength = 0.0;

  void setCharacteristicElementLength( double length );

  double getInterfaceThickness() const { return h; }

  /** Increment inputs (over the current step). Abar/DeltaA are full 3x3
   * tangential gradients stored row-major (vec), g is the average normal
   * gradient, normal is the (unit) interface normal. */
  struct KernelInput {
    const double* dAbar;   // 9
    const double* dDeltaA; // 9
    const double* dG;      // 3
    const double* normal;  // 3
  };

  /** Current (total) generalized stresses and the reduced constitutive
   * tangent blocks. All row-major. a = (Abar, DeltaA) has size 18. */
  struct KernelOutput {
    double* sAbar;   // 9
    double* sDeltaA; // 9
    double* tMat;    // 3
    double* H_aa;    // 18 x 18
    double* H_ag;    // 18 x 3
    double* H_ga;    // 3  x 18
    double* H_gg;    // 3  x 3
  };

  struct TimeIncrement {
    double timeOld;
    double dT;
  };

  /** Solve the local gamma-equilibrium and return reduced generalized
   * stresses and the consistent reduced tangent. Mutates the two bulk
   * state blocks, gamma, and sigma^+/sigma^- in stateVars (call on a scratch
   * copy for finite-difference perturbations). */
  virtual void computeMixedKernel( double*              stateVars,
                                   const KernelInput&   input,
                                   KernelOutput&        output,
                                   const TimeIncrement& timeIncrement );

  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    // Expose the embedded bulk materials' internal variables (e.g. the Von
    // Mises hardening variable kappa) for postprocessing: "kappaTop" /
    // "kappaBottom" resolve to "kappa" inside the respective side's state
    // block. Any other embedded-state name works the same way with the
    // "Top"/"Bottom" suffix (e.g. "kappaTop" -> topMaterial's "kappa").
    const auto delegate =
      [&]( const std::string& base, MarmotMaterialHypoElastic* mat, const char* block ) -> StateView {
      return mat->getStateView( base, stateLayout.getPtr( stateVars, block ) );
    };
    const std::size_t n = stateName.size();
    if ( n > 3 && stateName.compare( n - 3, 3, "Top" ) == 0 && topMaterial )
      return delegate( stateName.substr( 0, n - 3 ), topMaterial.get(), "topMaterialStateVars" );
    if ( n > 6 && stateName.compare( n - 6, 6, "Bottom" ) == 0 && bottomMaterial )
      return delegate( stateName.substr( 0, n - 6 ), bottomMaterial.get(), "bottomMaterialStateVars" );

    return stateLayout.getStateView( stateVars, stateName );
  }

  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  virtual void initializeYourself( double* stateVars, int nStateVars );

  virtual double getDensity();
};
