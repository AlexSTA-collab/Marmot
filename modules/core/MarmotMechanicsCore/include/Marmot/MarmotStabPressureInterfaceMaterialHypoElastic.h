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
 * Two-station equilibrated interface material with a MIXED (independent)
 * hydrostatic pressure.
 *
 * WHY
 * ---
 * With near-perfectly-plastic von Mises (H/E -> 0) the interface layer reaches
 * the incompressible limit: K/mu_plastic ~ 1e7. In a displacement-only
 * formulation the pressure is then evaluated as
 *      p = K * tr(eps)
 * where tr(eps) is driven to ~0 (measured: tr(eps) ~ 1e-4 against a deviatoric
 * strain of ~0.55, four orders of magnitude smaller). Multiplying a near-zero
 * quantity by K ~ 3e5 amplifies discretisation noise into a large spurious
 * element-to-element pressure checkerboard (measured p2p ~ 38 with the von
 * Mises stress uniform to 5 significant figures). Confirmed by a three-point
 * Poisson-ratio test: pressure oscillation ~ K^0.70 over a 10x range in K while
 * the deviatoric solution stays bit-identical.
 *
 * This material removes the amplification by making the pressure an INDEPENDENT
 * unknown p (supplied by the element as a nodal field) instead of a byproduct of
 * the volumetric strain. The stress used for equilibrium is
 *
 *      sigma~ = dev( sigma(eps) )  -  p * I
 *
 * and the volumetric constitutive law is returned as a residual to be enforced
 * weakly by the element (together with a pressure-gradient stabilisation, which
 * lives in the element, not here):
 *
 *      r_p = <tr(eps)>  +  p / K            ( <.> = mean of the two stations )
 *
 * RESTRICTION
 * -----------
 * Overriding the constitutive pressure is only admissible because the von Mises
 * yield surface is PRESSURE-INDEPENDENT: the deviatoric return-mapping is
 * unaffected by p, so taking the material's deviator and replacing its
 * volumetric part is consistent. For a pressure-DEPENDENT law (Drucker-Prager,
 * cap models, ...) p would have to be fed INTO the return-mapping and this
 * construction would be invalid.
 *
 * STRUCTURE OF THE TANGENT (derived, not approximated)
 * ---------------------------------------------------
 * With sigma~ = Pdev sigma - p I, which is linear in sigma and in p:
 *   d sigma~ / d vec(G) = Pdev * CFull   (=: CT),      d sigma~ / dp = -vec(I)
 * Two terms vanish identically, which keeps the mixed coupling purely geometric:
 *   - F_gamma = Rt ( sigma~+ - sigma~- ) is INDEPENDENT of p  (the -p I cancels)
 *       => d gamma / dp = 0
 *   - the mean volumetric strain is INDEPENDENT of gamma  (Bz+ + Bz- = 0)
 *       => d r_p / d gamma = 0
 * Hence, with x = (w, A+, A-) (21) and y = (f, S+, S-) (21):
 *   dy/dx = px - pgam * rgam^-1 * rx        (same condensation as the Y material)
 *   dy/dp = -sideThickness ( Bx+^T e + Bx-^T e )        (purely geometric)
 *   dr_p/dx = 1/2 e^T ( Bx+ + Bx- )                      (purely geometric)
 *   dr_p/dp = 1/K
 */
class MarmotStabPressureInterfaceMaterialHypoElastic {

protected:
  const double*                                materialProperties;
  const int                                    nMaterialProperties;
  double                                       h            = 0.0;
  double                                       bulkModulus  = 0.0;
  double                                       shearModulus = 0.0;
  std::vector< double >                        baseMaterialProperties;
  std::unique_ptr< MarmotMaterialHypoElastic > topMaterial;
  std::unique_ptr< MarmotMaterialHypoElastic > bottomMaterial;

public:
  const int materialNumber;

  MarmotStabPressureInterfaceMaterialHypoElastic( const std::string& materialName,
                                                  const double*      matProperties_,
                                                  int                nMaterialProperties_,
                                                  int                materialNumber_ );

  virtual ~MarmotStabPressureInterfaceMaterialHypoElastic() = default;

  MarmotStateLayoutDynamic stateLayout;

  double characteristicElementLength = 0.0;

  void setCharacteristicElementLength( double length );

  double getInterfaceThickness() const { return h; }
  double getBulkModulus() const { return bulkModulus; }
  double getShearModulus() const { return shearModulus; }

  /** Increments over the step. dSurfaceStrain holds [A+ (9); A- (9)]. */
  struct Deformation {
    const double* dU;               // 6  (top, bottom) -> w = dU[0:3]-dU[3:6]
    const double* dSurfaceStrain;   // 18
    const double* normal;           // 3
    const double* separationVector; // 3
    double        dPressure;        // increment of the independent pressure p
  };

  /** Totals returned to the element. All row-major. */
  struct Response {
    double* generalizedForce;   // 3   f      (from sigma~)
    double* surfaceStressPlus;  // 9   S+     (from sigma~)
    double* surfaceStressMinus; // 9   S-     (from sigma~)
    double* volumetricResidual; // 1   r_p = <tr eps> + p/K
  };

  /** Tangents. x = (w, A+, A-). */
  struct Tangents {
    double* Q_ww;   // 3x3
    double* Q_wAp;  // 3x9
    double* Q_wAm;  // 3x9
    double* Q_Apw;  // 9x3
    double* Q_ApAp; // 9x9
    double* Q_ApAm; // 9x9
    double* Q_Amw;  // 9x3
    double* Q_AmAp; // 9x9
    double* Q_AmAm; // 9x9
    double* Q_wp;   // 3x1   d f  / dp
    double* Q_App;  // 9x1   d S+ / dp
    double* Q_Amp;  // 9x1   d S- / dp
    double* Q_pw;   // 1x3   d r_p / dw
    double* Q_pAp;  // 1x9   d r_p / dA+
    double* Q_pAm;  // 1x9   d r_p / dA-
    double* Q_pp;   // 1x1   d r_p / dp = 1/K
  };

  struct TimeIncrement {
    double timeOld;
    double dT;
  };

  virtual void computeStress( double*              stateVars,
                              Response&            response,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );

  StateView getStateView( const std::string& stateName, double* stateVars ) const
  {
    const std::size_t n = stateName.size();
    if ( n > 3 && stateName.compare( n - 3, 3, "Top" ) == 0 && topMaterial )
      return topMaterial->getStateView( stateName.substr( 0, n - 3 ),
                                        stateLayout.getPtr( stateVars, "topMaterialStateVars" ) );
    if ( n > 6 && stateName.compare( n - 6, 6, "Bottom" ) == 0 && bottomMaterial )
      return bottomMaterial->getStateView( stateName.substr( 0, n - 6 ),
                                           stateLayout.getPtr( stateVars, "bottomMaterialStateVars" ) );
    return stateLayout.getStateView( stateVars, stateName );
  }

  virtual int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

  virtual void initializeYourself( double* stateVars, int nStateVars );

  virtual double getDensity();
};
