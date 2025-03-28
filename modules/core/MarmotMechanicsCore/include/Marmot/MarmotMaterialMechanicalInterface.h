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
 * Matthias Neuner matthias.neuner@uibk.ac.at
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
#include "Marmot/MarmotMaterial.h"
#include "Fastor/Fastor.h"

/**
 *  Abstract basic class for Mechanical materials with Interface.
 *  'Mechanical' is meant in the 'most general sense', i.e., any material which describes a mechanical (cauchy)
 *  stress - deformation relationship (e.g, hyperelastic, hypoelastic, elasto-plastic, visco-elastic materials)
 *
 *  force = g (force, surface_stress, dx, dxdX, t, .. ),
 *  surface_stress = h(force, surface_stress, dx, dxdX, t, ..)
 *
 *  formulated incrementally as:
 *  force_np = f (force_n, surface_stress_n, dx_n, dx_np, dxdX_n, dxdX_np, Δt, t_n, .. )
 *
 *  Algorithmic tangent with format: 
 *                 [[[dsurface_stress dsurface_strain]_{ijkl},[dforce dsurface_strain]_{ijk}, [0]],
 *                  [[0]                                     ,[dforce d dispalcement]_{ij}  , [0]],
 *                  [[0]                                     ,[0]                           , [dsurface_stress dsurface_strain]_{ijkl} ]] 
 * it can be interpreted as a Row major 2D Tensor, a matrix of shape 21x21.
 */
class MarmotMaterialMechanicalInterface : public MarmotMaterial {

public:
  using MarmotMaterial::MarmotMaterial;

  using Tensor1D = Fastor::Tensor<double,3>;
  using Tensor2D = Fastor::Tensor<double,3,3>;

  virtual void computeStress(
                             Tensor1D&  force,
                             Tensor2D&  surface_stress,
                             Fastor::Tensor<double, 21,21>& dStress_dStrain,
                             const Fastor::Tensor<double, 6,1>& dU,
                             const Fastor::Tensor<double, 18,1>& dSurface_strain,
                             const Tensor1D& normal,
                             const double* timeOld,
                             const double  dT,
                             double&       pNewDT ) = 0;

};
