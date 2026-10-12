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

#include "Marmot/MarmotFastorTensorBasics.h"

namespace Marmot::Materials::InterfaceFiniteStrainHelperFunctions {

  using Marmot::FastorStandardTensors::Tensor18d;
  using Marmot::FastorStandardTensors::Tensor3333d;
  using Marmot::FastorStandardTensors::Tensor333d;
  using Marmot::FastorStandardTensors::Tensor33d;
  using Marmot::FastorStandardTensors::Tensor3d;
  using Marmot::FastorStandardTensors::Tensor6d;
  using Marmot::FastorStandardTensors::Tensor9d;

  /**
   * @brief Geometry of the normal fibre of an interface and the deformation gradient it carries.
   *
   * The reduced normal derivative is read off the connector between the paired faces, and not off the jump divided
   * by the constitutive thickness. With @f$ dX = \ell N + d\tau @f$, the tangential offset @f$ d\tau @f$ is carried
   * by the surface map rather than by the jump, so that the image of the unit normal fibre is
   * @f[
   *   \bar{g} = \frac{1}{\ell} \left( [\![u]\!] - \langle \nabla_s u \rangle \, d\tau \right) ,
   * @f]
   * and the deformation gradient of the layer is
   * @f$ F = I + \langle \nabla_s u \rangle + \bar{g} \otimes N @f$. Leaving the correction out costs the frame
   * indifference of an inclined interface. For a purely normal connector it reduces to @f$ [\![u]\!] / h @f$.
   */
  struct NormalFibre {
    double    ell;                    ///< separation of the paired faces along the normal
    Tensor3d  dTau;                   ///< tangential offset of the connector
    Tensor33d averageSurfaceGradient; ///< average surface gradient of both sides
    Tensor3d  normalGradient;         ///< reduced normal derivative of the displacement
    Tensor33d F;                      ///< deformation gradient of the layer
  };

  /**
   * @brief Compute the normal fibre of an interface from its kinematics.
   *
   * @param[in] U total displacements of the two interface sides, top first
   * @param[in] surfaceGradient total surface gradients of the two interface sides, top first
   * @param[in] normal interface normal, pointing from the first to the second side
   * @param[in] connector reference connector between the paired faces
   * @return the geometry of the fibre and the deformation gradient
   *
   * @throws std::invalid_argument if the paired faces are not separated along the interface normal
   */
  NormalFibre normalFibre( const Tensor6d&  U,
                           const Tensor18d& surfaceGradient,
                           const Tensor3d&  normal,
                           const Tensor3d&  connector );

  /**
   * @brief First Piola-Kirchhoff stress and its derivative from the Kirchhoff stress of a finite-strain material.
   *
   * @param[in] tau Kirchhoff stress
   * @param[in] dTau_dF derivative of the Kirchhoff stress with respect to the deformation gradient
   * @param[in] F deformation gradient
   * @param[out] P first Piola-Kirchhoff stress
   * @param[out] A derivative of the first Piola-Kirchhoff stress with respect to the deformation gradient
   */
  void firstPiolaKirchhoff( const Tensor33d&   tau,
                            const Tensor3333d& dTau_dF,
                            const Tensor33d&   F,
                            Tensor33d&         P,
                            Tensor3333d&       A );

  /**
   * @brief Work conjugates of the reduced kinematics and their derivatives for a given stress state of the layer.
   *
   * The internal work per unit reference midsurface area is @f$ h \, P : \delta F @f$, so that the jump is
   * conjugate to @f$ t = ( h / \ell ) P N @f$ and the average surface gradient to
   * @f$ \Sigma = h P T - ( h / \ell ) ( P N ) \otimes d\tau @f$. The thicknesses @f$ h @f$ and @f$ \ell @f$
   * coincide whenever the faces are separated by the layer; @f$ \ell @f$ is mesh geometry, @f$ h @f$ the
   * constitutive layer.
   */
  struct ReducedResponse {
    Tensor3d    force;         ///< force, conjugate to the jump
    Tensor33d   surfaceStress; ///< surface stress, conjugate to the average surface gradient
    Tensor33d   Q_ik;          ///< derivative of the force with respect to the jump
    Tensor333d  H_ikR;         ///< derivative of the force with respect to the average surface gradient
    Tensor333d  H_iRk;         ///< derivative of the surface stress with respect to the jump
    Tensor3333d A_iRkP;        ///< derivative of the surface stress with respect to the average surface gradient
  };

  /**
   * @brief Reduce a first Piola-Kirchhoff stress, or a derivative of it, to the force and the surface stress.
   *
   * The reduction is linear in the stress, so that it applies to the derivative of the stress with respect to any
   * parameter, for example a nonlocal field, as well.
   *
   * @param[in] P first Piola-Kirchhoff stress of the layer, or its derivative
   * @param[in] normal interface normal
   * @param[in] fibre geometry of the normal fibre
   * @param[in] h constitutive thickness of the layer
   * @param[out] force force, conjugate to the jump
   * @param[out] surfaceStress surface stress, conjugate to the average surface gradient
   */
  void reduceStress( const Tensor33d&   P,
                     const Tensor3d&    normal,
                     const NormalFibre& fibre,
                     double             h,
                     Tensor3d&          force,
                     Tensor33d&         surfaceStress );

  /**
   * @brief Derivatives of a quantity with respect to the jump and to the average surface gradient.
   *
   * With the deformation gradient of the layer @f$ F = I + \langle \nabla_s u \rangle + \bar{g} \otimes N @f$,
   * a quantity @f$ q @f$ with the derivative @f$ G = \partial q / \partial F @f$ changes with the jump as
   * @f$ G N / \ell @f$ and with the average surface gradient as @f$ G T - ( G N ) \otimes d\tau / \ell @f$, where
   * @f$ T @f$ projects on the tangent plane.
   *
   * @param[in] G derivative of the quantity with respect to the deformation gradient
   * @param[in] normal interface normal
   * @param[in] fibre geometry of the normal fibre
   * @param[out] dJump derivative with respect to the jump
   * @param[out] dAverage derivative with respect to the average surface gradient
   */
  void reduceDerivative( const Tensor33d&   G,
                         const Tensor3d&    normal,
                         const NormalFibre& fibre,
                         Tensor3d&          dJump,
                         Tensor33d&         dAverage );

  /**
   * @brief Reduce the stress state of the layer to the force and the surface stress of the interface.
   *
   * @param[in] P first Piola-Kirchhoff stress of the layer
   * @param[in] A derivative of the first Piola-Kirchhoff stress with respect to the deformation gradient
   * @param[in] normal interface normal
   * @param[in] fibre geometry of the normal fibre
   * @param[in] h constitutive thickness of the layer
   * @return the force, the surface stress and the tangents of the interface
   */
  ReducedResponse reduceToInterface( const Tensor33d&   P,
                                     const Tensor3333d& A,
                                     const Tensor3d&    normal,
                                     const NormalFibre& fibre,
                                     double             h );

} // namespace Marmot::Materials::InterfaceFiniteStrainHelperFunctions
