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
 * Two-face interface kernel with the gradient jump z as an EXTERNAL unknown and
 * TWO independent pressures.
 *
 * WHY TWO
 * -------
 * MarmotStabPressureInterfaceMaterialHypoElastic carries one mixed pressure and
 * documents two exact cancellations. Both survive when z is promoted to an
 * external field, and together they leave a hole:
 *
 *   d r_p / dz = 0   because tr( sym(z (x) n) ) = z.n enters the two faces as
 *                    +z.n/2 and -z.n/2 and cancels from the MEAN <tr eps>;
 *   d r_z / dp = 0   because a single p shared by both faces cancels from the
 *                    traction jump [ t~ ].
 *
 * So a single mean pressure does not see z at all. The volumetric mode the
 * gradient jump carries lives in the JUMP of the volumetric strain,
 *
 *   [ tr eps ] = z.n + [ tr sym(A) ] ,
 *
 * and nothing constrains it. Under isochoric plastic flow BOTH faces must be
 * volume preserving, not just their average, so p = K tr eps is amplified on
 * each face separately exactly as it was before the pressure was made mixed.
 * That is the traction pollution measured for ZIQUAD4 (4.9x the condensed
 * element's traction-jump error in the hard-yielding family, unchanged by the
 * gradient-jump regularization zeta).
 *
 * This kernel therefore carries the pressure as a pair
 *
 *   p^pm = pbar +- [p]/2 ,
 *
 * one per face, parameterized by its mean and its jump, mirroring the
 * (Abar, DeltaA) and (gbar, z) split the rest of the formulation already uses.
 * The stress driving equilibrium is, per face,
 *
 *   sigma~^pm = dev( sigma^pm ) - p^pm I ,
 *
 * and the two volumetric laws are returned as residuals for the element to
 * enforce weakly,
 *
 *   r_pbar = < tr eps > + pbar / K ,
 *   r_pjump = [ tr eps ] + [p] / K .
 *
 * The pairing is now complete and the two new couplings are exactly the ones
 * the single-pressure form was missing:
 *
 *   d r_z    / d[p]  = -(h/4) n ,
 *   d r_pjump / dz   = n .
 *
 * while d r_z/d pbar and d r_pbar/dz remain identically zero, so pbar keeps its
 * original role untouched.
 *
 * RESTRICTION
 * -----------
 * As in the single-pressure kernel, overriding the constitutive pressure is
 * admissible only because the von Mises yield surface is PRESSURE-INDEPENDENT:
 * the deviatoric return mapping does not see p, so replacing the volumetric
 * part is consistent. For a pressure-dependent law (Drucker-Prager, cap models)
 * p would have to enter the return mapping and this construction is invalid.
 *
 * INF-SUP
 * -------
 * (z.n, [p]) is an equal-order Q1/Q1 pair on the same midsurface nodes, exactly
 * like (u, pbar). Both need the element's pressure-gradient stabilisation; the
 * element applies it to BOTH fields. Nothing here stabilises anything -- this
 * kernel is constitutive only.
 *
 * Material properties: E, nu, h, then the base material's own.
 */
class MarmotZStabPressureInterfaceMaterialHypoElastic {

protected:
  const double*                                materialProperties;
  const int                                    nMaterialProperties;
  double                                       h             = 0.0;
  double                                       bulkModulus   = 0.0;
  double                                       shearModulus  = 0.0;
  double                                       lameParameter = 0.0;
  std::vector< double >                        baseMaterialProperties;
  std::unique_ptr< MarmotMaterialHypoElastic > topMaterial;
  std::unique_ptr< MarmotMaterialHypoElastic > bottomMaterial;
  double                                       regularization;

public:
  static constexpr int nX = 21; /**< x = (w, A+, A-) */
  static constexpr int nZ = 3;  /**< z = [[u_{k,r}]] n_r */

  static constexpr int offW  = 0;
  static constexpr int offAp = 3;
  static constexpr int offAm = 12;

  const int materialNumber;

  MarmotZStabPressureInterfaceMaterialHypoElastic( const std::string& materialName,
                                                   const double*      matProperties_,
                                                   int                nMaterialProperties_,
                                                   int                materialNumber_ );

  virtual ~MarmotZStabPressureInterfaceMaterialHypoElastic() = default;

  MarmotStateLayoutDynamic stateLayout;

  double characteristicElementLength = 0.0;

  void setCharacteristicElementLength( double length );

  double getInterfaceThickness() const { return h; }
  double getBulkModulus() const { return bulkModulus; }
  double getShearModulus() const { return shearModulus; }

  double getGradientJumpRegularization() const { return regularization; }
  void   setGradientJumpRegularization( double zeta );

  /** Increments over the step. */
  struct Deformation {
    const double* dU;                  /**< 6:  [ du+ ; du- ] */
    const double* dSurfaceStrain;      /**< 18: [ dA+ ; dA- ] */
    const double* dNormalGradientJump; /**< 3:  dz */
    const double* normal;              /**< 3 */
    const double* separationVector;    /**< 3 */
    double        dPressureMean;       /**< d pbar */
    double        dPressureJump;       /**< d [p]  */
  };

  /** Totals returned to the element. All row-major. */
  struct Response {
    double* generalizedForce;       /**< 3  f,  from sigma~ */
    double* surfaceStressPlus;      /**< 9  S+, from sigma~ */
    double* surfaceStressMinus;     /**< 9  S-, from sigma~ */
    double* tractionImbalance;      /**< 3  r_z = (h/4)([dev sigma] n - [p] n) + reg */
    double* volumetricResidualMean; /**< 1  r_pbar  = <tr eps> + pbar/K */
    double* volumetricResidualJump; /**< 1  r_pjump = [tr eps] + [p]/K  */
  };

  /**
   * Hessian over ( x(21), z(3), pbar(1), [p](1) ). Nothing is condensed.
   * The two blocks that vanish identically are still written (as zeros) so the
   * element never has to know which they are.
   */
  struct Tangents {
    double* K_xx;   /**< 21 x 21 */
    double* K_xz;   /**< 21 x  3 */
    double* K_xpm;  /**< 21 x  1 */
    double* K_xpj;  /**< 21 x  1 */
    double* K_zx;   /**<  3 x 21 */
    double* K_zz;   /**<  3 x  3 */
    double* K_zpm;  /**<  3 x  1  == 0 */
    double* K_zpj;  /**<  3 x  1  == -(h/4) n */
    double* K_pmx;  /**<  1 x 21 */
    double* K_pmz;  /**<  1 x  3  == 0 */
    double* K_pmpm; /**<  1 x  1  == 1/K */
    double* K_pjx;  /**<  1 x 21 */
    double* K_pjz;  /**<  1 x  3  == n */
    double* K_pjpj; /**<  1 x  1  == 1/K */
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
