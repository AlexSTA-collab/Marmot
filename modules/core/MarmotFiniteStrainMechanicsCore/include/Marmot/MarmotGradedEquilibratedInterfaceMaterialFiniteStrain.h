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

#include "Marmot/MarmotEquilibratedInterfaceMaterialFiniteStrainBase.h"
#include "Marmot/MarmotGradedInterfaceMaterialFiniteStrainBase.h"

#include <string>

/**
 * @brief Two-sided finite-strain interface material with a stiffness that varies through the thickness of the layer.
 *
 * This material is the finite-strain counterpart of MarmotGradedXInterfaceMaterialHypoElastic and the generalization
 * of the equilibrated two-sided material (MarmotEquilibratedInterfaceMaterialFiniteStrain), which has the two faces
 * as the only integration points through the thickness, to an arbitrary number of integration points with different
 * bulk materials. It has the contract of the equilibrated material, so that it is used by the same element
 * formulation: it receives the total displacements and the total surface gradients @f$A^+@f$ and @f$A^-@f$ of the
 * two faces, and returns the generalized force conjugate to the jump and the two generalized surface stresses
 * conjugate to @f$A^+@f$ and @f$A^-@f$, based on the first Piola-Kirchhoff stress.
 *
 * With the layer coordinate @f$y \in [0,1]@f$ from the bottom to the top, the deformation gradient at the
 * integration point @f$q@f$ is
 * @f[
 *   F_q = I + \frac{1}{\ell} [\![u]\!] \otimes N + (1 - y_q) A^- + y_q A^+ + B_\tau \bar{A} + g_q \otimes N , \qquad
 *   \sum_q \omega_q g_q = 0 ,
 * @f]
 * with the normal separation @f$\ell@f$ and the tangential separation of the faces from the connector, and the
 * correction @f$g_q = \partial_Y u^1@f$ of the microproblem, which requires that the traction @f$P_q N@f$ is the same
 * at all points. The generalized stress is @f$p = h \sum_q \omega_q B_q^T P_q@f$ with the derivative @f$B_q@f$ of the
 * gradient with respect to @f$([\![u]\!], A^+, A^-)@f$, and the tangent is the one condensed over the corrections.
 * The tangential gradient of the correction is kept in the way described in
 * MarmotGradedXInterfaceMaterialHypoElastic, with the derivative of the correction of the first Piola-Kirchhoff
 * tangent and the divergence of the first Piola-Kirchhoff stress in the reference configuration.
 *
 * The material properties and the state variables are those of MarmotGradedInterfaceMaterialFiniteStrain, with the
 * optional last property of Marmot::Materials::ThicknessProfile, which switches the closure of the tangential
 * gradient off. With two integration points at the faces and without the closure, the response is the one of the
 * equilibrated material.
 */
class MarmotGradedEquilibratedInterfaceMaterialFiniteStrain : public MarmotGradedInterfaceMaterialFiniteStrainBase {

public:
  using Tensor3d    = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tensor3d;    ///< first order tensor in 3D
  using Tensor33d   = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tensor33d;   ///< second order tensor in 3D
  using Tensor333d  = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tensor333d;  ///< third order tensor in 3D
  using Tensor3333d = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tensor3333d; ///< fourth order tensor in 3D
  using Tensor6d    = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tensor6d;    ///< vector of six components
  using Tensor18d   = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tensor18d;   ///< vector of 18 components

  using State         = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::State;         ///< response, as there
  using Tangents      = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Tangents;      ///< block tangent
  using Deformation   = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::Deformation;   ///< kinematics
  using TimeIncrement = MarmotEquilibratedInterfaceMaterialFiniteStrainBase::TimeIncrement; ///< time data

  /**
   * @brief Construct the material from a registered finite-strain bulk material.
   *
   * @param[in] materialName name of the bulk material, registered in the MarmotMaterialFiniteStrainFactory
   * @param[in] matProperties material properties, see Marmot::Materials::ThicknessProfile
   * @param[in] nMaterialProperties number of material properties
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if the properties are inconsistent or the bulk material is not registered
   */
  MarmotGradedEquilibratedInterfaceMaterialFiniteStrain( const std::string& materialName,
                                                         const double*      matProperties,
                                                         int                nMaterialProperties,
                                                         int                materialNumber_ )
    : MarmotGradedInterfaceMaterialFiniteStrainBase( materialName, matProperties, nMaterialProperties, materialNumber_ )
  {
  }

  /**
   * @brief Solve the microproblem of the layer and compute the generalized stress and the condensed tangent.
   *
   * @param[in,out] state generalized stress and state variables
   * @param[out] tangents condensed block tangent over the jump and the two surface gradients
   * @param[in] deformation total displacements and surface gradients of both faces, the normal and the connector
   * @param[in] timeIncrement time at the beginning of the increment and time increment
   *
   * @throws Marmot::StressUpdateFailed if the iteration for the microproblem does not converge, or if the acoustic
   * tensor of an integration point is singular
   * @throws std::invalid_argument if the normal is zero, or the connector has no positive normal component
   */
  virtual void computeStress( State&               state,
                              Tangents&            tangents,
                              const Deformation&   deformation,
                              const TimeIncrement& timeIncrement );
};
