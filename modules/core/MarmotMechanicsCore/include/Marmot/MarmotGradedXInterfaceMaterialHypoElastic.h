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

#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotGradedInterfaceMaterialHypoElasticBase.h"

#include <string>

/**
 * @brief Two-sided hypoelastic interface material with a stiffness that varies through the thickness of the layer.
 *
 * This material is the generalization of the equilibrated two-sided interface material
 * (MarmotEquilibratedXInterfaceMaterialHypoElastic), which has the two faces as the only integration points through
 * the thickness, to an arbitrary number of integration points, whose bulk materials differ. It has the contract of
 * the equilibrated material, so that it is used by the same element formulation: the surface gradients @f$A^+@f$ and
 * @f$A^-@f$ of the top and of the bottom face are given separately, and the response is the generalized force
 * conjugate to the jump @f$w = [\![u]\!]@f$ and the two generalized surface stresses conjugate to @f$A^+@f$ and
 * @f$A^-@f$.
 *
 * With the layer coordinate @f$y \in [0,1]@f$ from the bottom to the top, the strain gradient at the integration
 * point @f$q@f$ is
 * @f[
 *   G_q = \frac{1}{\ell} w \otimes n + (1 - y_q) A^- + y_q A^+ + B_\tau \bar{A} + g_q \otimes n , \qquad
 *   \sum_q \omega_q g_q = 0 ,
 * @f]
 * where @f$\ell@f$ is the normal separation and @f$B_\tau \bar{A}@f$ the correction of the tangential separation of
 * the faces, exactly as for the equilibrated material, and @f$g_q = \partial_y u^1@f$ is the correction of the
 * microproblem. The surface gradient is interpolated linearly between the faces, which keeps the terms of the
 * bending type that the two surface gradients carry. The microproblem requires that the traction
 * @f$\sigma_q n@f$ is the same at all points (see Marmot::Materials::ThicknessMicroproblem), and is solved at every
 * call. The generalized stress is
 * @f[
 *   p = h \sum_q \omega_q \hat{B}_q^T \sigma_q , \qquad
 *   \hat{B}_q = \frac{\partial G_q}{\partial (w, A^+, A^-)} ,
 * @f]
 * where the derivative of the gradient of the point contains the response of the microproblem and the tangential
 * gradient of the correction below. The test function is therefore built from the same operator as the strain, the
 * generalized stress is the derivative of the energy of the layer, and the tangent, condensed over the corrections,
 * is symmetric for an elastic layer. With two integration points at the faces (mode 2 of the profile) and without
 * the closure below, the response is identical to the one of the equilibrated material.
 *
 * @par Tangential gradient of the correction
 * The variational form has the terms @f$h \, \nabla_T u^1@f$, which are kept. The correction @f$u^1@f$ varies along
 * the surface because the strain of the layer does, through the jump of the displacement: its tangential gradient
 * is @f$\nabla_T [\![u]\!] = A^+ - A^-@f$. The material therefore uses the closure
 * @f[
 *   \partial_r \varepsilon_q = \frac{1}{\ell} \mathrm{sym} \left( ( A^+ - A^- ) e_r \otimes n \right) , \qquad
 *   \partial_r w_q = \frac{\partial w_q}{\partial \varepsilon} \, \partial_r \varepsilon_q ,
 * @f]
 * with the derivative of the correction of the condensed microproblem, and neglects the second derivatives of the
 * displacement of the faces. This has two consequences, both proportional to @f$A^+ - A^-@f$ and both of the order
 * of the bending terms: the strain of the point gains the tangential gradient @f$h \nabla_T u^1(y_q)@f$ with
 * @f$u^1(y) = \int_0^y w \, \mathrm{d}y@f$, and the traction changes through the thickness by
 * @f$\partial_y ( \sigma n ) = - h \, \partial_r \sigma e_r@f$, the tangential divergence of the stress. The
 * traction term is present for a homogeneous layer as well, because a surface gradient of the jump produces a
 * tangential divergence of the stress there too. The response is therefore not the one of the equilibrated
 * material, which has neither term, unless the optional last material property (see
 * Marmot::Materials::ThicknessProfile) is zero, which switches the closure off. The derivative of the correction
 * is the one of the converged state of the increment, which is iterated until it does not change, and the tangent
 * neglects its own derivative, which is of the order of the second derivative of the stress.
 *
 * The material properties and the state variables are those of MarmotGradedInterfaceMaterialHypoElastic.
 */
class MarmotGradedXInterfaceMaterialHypoElastic : public MarmotGradedInterfaceMaterialHypoElasticBase {

public:
  using Tensor3d    = MarmotEquilibratedXInterfaceMaterialHypoElastic::Tensor3d;    ///< first order tensor in 3D
  using Tensor33d   = MarmotEquilibratedXInterfaceMaterialHypoElastic::Tensor33d;   ///< second order tensor in 3D
  using Tensor333d  = MarmotEquilibratedXInterfaceMaterialHypoElastic::Tensor333d;  ///< third order tensor in 3D
  using Tensor3333d = MarmotEquilibratedXInterfaceMaterialHypoElastic::Tensor3333d; ///< fourth order tensor in 3D
  using Tensor6d    = MarmotEquilibratedXInterfaceMaterialHypoElastic::Tensor6d;    ///< vector of six components
  using Tensor18d   = MarmotEquilibratedXInterfaceMaterialHypoElastic::Tensor18d;   ///< vector of 18 components

  using IncrementMemory = MarmotEquilibratedXInterfaceMaterialHypoElastic::IncrementMemory; ///< not used, kept for
                                                                                            ///< the element contract
  using State = MarmotEquilibratedXInterfaceMaterialHypoElastic::State;                     ///< response, as in the
                                                                                            ///< equilibrated material
  using Tangents      = MarmotEquilibratedXInterfaceMaterialHypoElastic::Tangents;          ///< block tangent, as there
  using Deformation   = MarmotEquilibratedXInterfaceMaterialHypoElastic::Deformation;       ///< kinematics, as there
  using TimeIncrement = MarmotEquilibratedXInterfaceMaterialHypoElastic::TimeIncrement;     ///< time data, as there

  /**
   * @brief Construct the two-sided graded interface material from a registered hypoelastic bulk material.
   *
   * @param[in] materialName name of the bulk material, registered in the MarmotMaterialHypoElasticFactory
   * @param[in] matProperties material properties, see Marmot::Materials::ThicknessProfile
   * @param[in] nMaterialProperties number of material properties
   * @param[in] materialNumber_ number of the material
   *
   * @throws std::invalid_argument if the properties are inconsistent, or if the bulk material is not registered
   */
  MarmotGradedXInterfaceMaterialHypoElastic( const std::string& materialName,
                                             const double*      matProperties,
                                             int                nMaterialProperties,
                                             int                materialNumber_ )
    : MarmotGradedInterfaceMaterialHypoElasticBase( materialName, matProperties, nMaterialProperties, materialNumber_ )
  {
  }

  /**
   * @brief Solve the microproblem of the layer and compute the generalized stress and the condensed tangent.
   *
   * @param[in,out] state generalized stress and state variables
   * @param[out] tangents condensed block tangent over the jump and the two surface gradients
   * @param[in] deformation increments of the displacements and of the surface gradients of both faces, the
   * normal and the connector between the faces
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

  /**
   * @brief Resultants of the stress through the thickness of the last call of computeStress().
   *
   * The generalized stress of the state is the derivative of the energy of the layer, which has the response of the
   * microproblem and the tangential gradient of the correction in its test function. The resultants are the mean of
   * the stress and its moments @f$ h \int y \sigma \, dy @f$ and @f$ h \int (1 - y) \sigma \, dy @f$, which are
   * what the exact solution of the layer is compared with. They differ from the generalized stress by a divergence
   * along the surface.
   *
   * @return force, moment with respect to the top face and moment with respect to the bottom face
   */
  const Eigen::Matrix< double, 21, 1 >& stressResultant() const { return resultant_; }

private:
  Eigen::Matrix< double, 21, 1 > resultant_ = Eigen::Matrix< double, 21, 1 >::Zero(); ///< resultants of the last call
};
