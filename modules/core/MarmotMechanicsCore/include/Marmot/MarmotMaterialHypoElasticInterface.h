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
#include <ios>
#include <iostream>
#include "Marmot/MarmotMaterialMechanicalInterface.h"
#include "Fastor/Fastor.h"

/**
 *
 * Derived abstract base class for elastic materials expressed purely in rate form.
 *
 * In general, the nominal stress rate tensor \f$ \sigRate \f$ can be written as a function of the nominal stress tensor
 * \f$ \sig \f$, the stretching rate tensor \f$ \epsRate \f$ and the time \f$ t \f$.
 *
 * \f[  \displaystyle \sigRate = f( \sig, \epsRate, t, ...) \f]
 *
 * In course of numerical time integration, this relation will be formulated incrementally as
 *
 * \f[  \displaystyle \Delta \sig = f ( \sig_n, \Delta\eps, \Delta t, t_n, ...) \f]
 *
 * with
 *
 * \f[  \displaystyle \Delta\eps =  \epsRate\, \Delta t \f]
 *
 * and the algorithmic tangent
 *
 * \f[ \displaystyle \frac{d \sig }{d \eps } =  \frac{d \Delta \sig }{d \Delta \eps } \f]
 *
 * This formulation is compatible with an Abaqus interface.
 */
class MarmotMaterialHypoElasticInterface : public MarmotMaterialMechanicalInterface {

public:
  using MarmotMaterialMechanicalInterface::MarmotMaterialMechanicalInterface;
  using Tensor1D = Fastor::Tensor<double,3>;
  using Tensor2D = Fastor::Tensor<double,3,3>;

  /// Characteristic element length
  double characteristicElementLength;
  /**
   * Set the characteristic element length at the considered quadrature point.
   * It is needed for the regularization of materials with softening behavior based on the mesh-adjusted softening
   * modulus.
   *
   * @param[in] length characteristic length; will be assigned to @ref characteristicElementLength
   */
  void setCharacteristicElementLength( double length );

  /**
   * For a given linearized strain increment \f$\Delta\boldsymbol{\varepsilon}\f$ at the old and the current time,
   * compute the Cauchy stress and the algorithmic tangent
   * \f$\frac{\partial\boldsymbol{\sigma}^{(n+1)}}{\partial\boldsymbol{\varepsilon}^{(n+1)}}\f$.
   *
   * @param[in,out]	force           conjugate force due to displacement jump
   * @param[in,out]	surface_stress  conjugate surface_stress due to conjugate average surface strain
   * @param[in,out]	dStressDDstrain	Algorithmic tangent representing the derivatives of:
   *                                    1)  The Cauchy surface stress tensor with respect to:
   *                                        a) the linearized average surface strain (the homogenization procedure -Gu et al 2011- 
   *                                           returns two material tensors Z_{ijkl}, Yn_H_inv_nF_{ijkl}), 
   *                                        b) the linearized displacement jump (H_inv_nF_{ijk})
   *                                    2)  The force vector with respect to:
   *                                        a) the linearized displacement jump (H_inv_{ij})
   * @param[in]	dU linearized displacement increment on "top (0)" and "bottom (1)" sides of the interface
   * @param[in] dSurface_strain linearized surface strain increment on top and bottom sides of the interface
   * @param[in] normal The normal to the interface positive to the direction of the "top (0)" side
   * @param[in]	timeOld	Old (pseudo-)time
   * @param[in]	dt	(Pseudo-)time increment from the old (pseudo-)time to the current (pseudo-)time
   * @param[in,out]	pNewDT	Suggestion for a new time increment
   */
  virtual void computeStress( double*  force,
                        double*  surface_stress,
                        double* dStress_dStrain,
                        const double* dU,
                        const double* dSurface_strain,
                        const double* normal,
                        const double* timeOld,
                        const double  dT,
                        double&       pNewDT)=0;



  //virtual void computeStress(
  //                           Tensor1D&  force,
  //                           Tensor2D&  surface_stress,
  //                           Fastor::Tensor<double, 21,21>& dStress_dStrain,
  //                           const Fastor::Tensor<double, 6,1>& dU,
  //                           const Fastor::Tensor<double, 18,1>& dSurface_strain,
  //                           const Tensor1D& normal,
  //                           const double* timeOld,
  //                           const double  dT,
  //                           double&       pNewDT );
};
