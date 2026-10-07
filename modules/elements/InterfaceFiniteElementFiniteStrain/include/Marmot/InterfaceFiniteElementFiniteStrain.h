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

/**
 * @file InterfaceFiniteElementFiniteStrain.h
 * @brief Interface finite element for finite strains.
 *
 * This file defines the formulation `Marmot::Elements::FiniteStrainDisplacementJumpFormulation` and the element
 * `Marmot::Elements::InterfaceFiniteElementFiniteStrain`, which is the common base of the interface elements
 * (InterfaceFiniteElementBase.h) with this formulation. In contrast to the small-strain element, the interface
 * material (MarmotInterfaceMaterialFiniteStrain) receives the total displacements and the total surface gradients
 * at the end of the increment, and evaluates a bulk material of finite strain for the deformation gradient of the
 * layer.
 */
#pragma once

#include "Marmot/InterfaceFiniteElement.h"
#include "Marmot/InterfaceFiniteElementBase.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotInterfaceMaterialFiniteStrain.h"

#include <Eigen/Dense>

namespace Marmot::Elements {

  /**
   * @struct FiniteStrainDisplacementJumpFormulation
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @brief Formulation of the interface element with displacement-jump kinematics at finite strain.
   *
   * The interface material (`MarmotInterfaceMaterialFiniteStrain`) receives the total displacements of both
   * interface sides, their total surface gradients, the normal of the interface and the connector between the
   * paired faces at the quadrature point, and returns the force and the surface stress together with the tangents.
   * The residual and the tangent are assembled with the jump operator and the operator of the average surface
   * gradient of the quadrature point. The state variables of a quadrature point are the ones of the
   * DisplacementJumpFormulation; the displacements and the surface gradients stored in them are the total ones.
   */
  template < int nDim, int nNodes >
  struct FiniteStrainDisplacementJumpFormulation {

    static constexpr int nTensor = nDim * nDim;

    using Small = DisplacementJumpFormulation< nDim, nNodes >;

    using ForceSized                = typename Small::ForceSized;
    using SurfaceStressSized        = typename Small::SurfaceStressSized;
    using InterfaceDisplSized       = typename Small::InterfaceDisplSized;
    using InterfaceSurfaceGradSized = typename Small::InterfaceSurfaceGradSized;

    using QMatrixSized = typename Small::QMatrixSized;
    using ZMatrixSized = typename Small::ZMatrixSized;
    using HMatrixSized = typename Small::HMatrixSized;

    /// matrix of the derivative of the surface stress with respect to the jump
    using HAverageJumpMatrixSized = Eigen::Matrix< double, nTensor, nDim, Eigen::RowMajor >;

    using Material          = MarmotInterfaceMaterialFiniteStrain; ///< interface material of the element
    using QPStateVarManager = typename Small::QPStateVarManager;   ///< named state variables of a quadrature point

    /// This formulation keeps no data at the quadrature points besides the state variables.
    struct QuadraturePointData {};

    /**
     * @brief Complete the geometry of a quadrature point; nothing is to be done in this formulation.
     *
     * @tparam Element type of the element
     * @tparam QuadraturePoint type of the quadrature point
     */
    template < class Element, class QuadraturePoint >
    static void initializeQuadraturePoint( const Element&, QuadraturePoint& )
    {
    }

    /**
     * @brief Evaluate the material at a quadrature point and add its contributions to the element.
     *
     * Updates the state variables of the quadrature point and adds the contributions of the quadrature point to
     * the residual and to the tangent of the element.
     *
     * @tparam Element type of the element
     * @tparam QuadraturePoint type of the quadrature point
     * @tparam QType type of the total nodal displacements
     * @tparam PeType type of the element residual
     * @tparam KeType type of the element tangent
     * @param[in] element the element, which holds the material
     * @param[in,out] qp quadrature point, whose state variables are updated
     * @param[in] QTotal total nodal displacements of the element at the end of the increment
     * @param[in] dU_GPs increment of the displacements of both interface sides, which this formulation does not need
     * @param[in] dSurface_strain_GPs increment of the surface gradients, which this formulation does not need
     * @param[in] time time at the beginning of the increment
     * @param[in] dT time increment
     * @param[in,out] Pe element residual, to which the contribution is added
     * @param[in,out] Ke element tangent, to which the contribution is added
     */
    template < class Element, class QuadraturePoint, class QType, class PeType, class KeType >
    static void integrate( Element&                         element,
                           QuadraturePoint&                 qp,
                           const QType&                     QTotal,
                           const InterfaceDisplSized&       dU_GPs,
                           const InterfaceSurfaceGradSized& dSurface_strain_GPs,
                           double                           time,
                           double                           dT,
                           PeType&                          Pe,
                           KeType&                          Ke )
    {
      (void)dU_GPs;
      (void)dSurface_strain_GPs;

      using namespace Marmot;

      constexpr int halfSize = nNodes * nDim / 2;

      const auto& Njump    = qp.NmatJump;
      const auto& Bavg     = qp.BmatAverage;
      auto&       material = element.material;

      // total kinematics of both interface sides at the quadrature point, top side first
      const auto QTotalBottom = QTotal.template segment< halfSize >( 0 );
      const auto QTotalTop    = QTotal.template segment< halfSize >( halfSize );

      InterfaceDisplSized totalU_GPs;
      totalU_GPs.template segment< nDim >( 0 )    = qp.NmatSide * QTotalTop;
      totalU_GPs.template segment< nDim >( nDim ) = qp.NmatSide * QTotalBottom;

      InterfaceSurfaceGradSized totalSurfaceGradient_GPs;
      totalSurfaceGradient_GPs.template segment< nTensor >( 0 )       = qp.BmatSide * QTotalTop;
      totalSurfaceGradient_GPs.template segment< nTensor >( nTensor ) = qp.BmatSide * QTotalBottom;

      // The reference connector between the paired faces at this quadrature point. The reduced normal derivative is
      // taken along it, so that the material needs it whole: its normal part is the fibre, its tangential part an
      // offset that the surface map carries.
      const typename Element::VectorDim connector = qp.NmatSide * ( element.getSideCoordinates( 1 ) -
                                                                    element.getSideCoordinates( 0 ) );

      ForceSized              force;
      SurfaceStressSized      surfaceStress;
      QMatrixSized            Q_ik;
      HMatrixSized            H_ikR;
      HAverageJumpMatrixSized H_iRk;
      ZMatrixSized            A_iRkP;

      Material::State         materialState{ Material::Tensor3d( 0.0 ),
                                     Material::Tensor33d( 0.0 ),
                                     qp.managedStateVars->materialStateVars.data() };
      Material::Tangents      materialTangents{ Material::Tensor33d( 0.0 ),
                                           Material::Tensor333d( 0.0 ),
                                           Material::Tensor333d( 0.0 ),
                                           Material::Tensor3333d( 0.0 ) };
      Material::TimeIncrement materialTimeIncrement{ time, dT };

      if constexpr ( nDim == 3 ) {
        Material::Deformation materialDeformation{ Material::Tensor6d( totalU_GPs.data() ),
                                                   Material::Tensor18d( totalSurfaceGradient_GPs.data() ),
                                                   Material::Tensor3d( qp.normal.data() ),
                                                   Material::Tensor3d( connector.data() ) };

        material->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );

        force         = mapEigenToFastor( materialState.force );
        surfaceStress = Eigen::Map< const SurfaceStressSized >( materialState.surfaceStress.data() );
        Q_ik          = mapEigenToFastor( materialTangents.Q_ik );
        H_ikR         = mapEigenToFastor( materialTangents.H_ikR );
        H_iRk         = Eigen::Map< const HAverageJumpMatrixSized >( materialTangents.H_iRk.data() );
        A_iRkP        = mapEigenToFastor( materialTangents.A_iRkP );
      }
      else if constexpr ( nDim == 2 ) {
        /*
         * MarmotInterfaceMaterialFiniteStrain is fixed to 3D Fastor tensor types, so a 2D element embeds its
         * kinematics into the 3D layout and pulls the response back out, see expandSidesTo3D().
         */
        const Eigen::Matrix< double, 6, 1 >  U3d               = expandSidesTo3D< nDim, 1 >( totalU_GPs );
        const Eigen::Matrix< double, 18, 1 > surfaceGradient3d = expandSidesTo3D< nDim, 2 >( totalSurfaceGradient_GPs );

        Material::Deformation materialDeformation{ Material::Tensor6d( U3d.data() ),
                                                   Material::Tensor18d( surfaceGradient3d.data() ),
                                                   expandTo3D( Fastor::Tensor< double, 2 >( qp.normal.data() ) ),
                                                   expandTo3D( Fastor::Tensor< double, 2 >( connector.data() ) ) };

        material->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );

        const Fastor::Tensor< double, 2 >          force2d         = reduceTo2D< U >( materialState.force );
        const Fastor::Tensor< double, 2, 2 >       surfaceStress2d = reduceTo2D< U, U >( materialState.surfaceStress );
        const Fastor::Tensor< double, 2, 2 >       Q2d             = reduceTo2D< U, U >( materialTangents.Q_ik );
        const Fastor::Tensor< double, 2, 2, 2 >    H2d             = reduceTo2D< U, U, U >( materialTangents.H_ikR );
        const Fastor::Tensor< double, 2, 2, 2 >    HAverageJump2d  = reduceTo2D< U, U, U >( materialTangents.H_iRk );
        const Fastor::Tensor< double, 2, 2, 2, 2 > A2d = reduceTo2D< U, U, U, U >( materialTangents.A_iRkP );

        force         = mapEigenToFastor( force2d );
        surfaceStress = Eigen::Map< const SurfaceStressSized >( surfaceStress2d.data() );
        Q_ik          = mapEigenToFastor( Q2d );
        H_ikR         = mapEigenToFastor( H2d );
        H_iRk         = Eigen::Map< const HAverageJumpMatrixSized >( HAverageJump2d.data() );
        A_iRkP        = mapEigenToFastor( A2d );
      }

      qp.managedStateVars->force         = force;
      qp.managedStateVars->surfaceStress = surfaceStress;
      qp.managedStateVars->displacement  = totalU_GPs;
      qp.managedStateVars->surfaceStrain = totalSurfaceGradient_GPs;

      Pe += Njump.transpose() * force * qp.J0xW;
      Pe += Bavg.transpose() * surfaceStress * qp.J0xW;

      Ke += ( Njump.transpose() * Q_ik * Njump + Njump.transpose() * H_ikR * Bavg + Bavg.transpose() * H_iRk * Njump +
              Bavg.transpose() * A_iRkP * Bavg ) *
            qp.J0xW;
    }
  };

  /**
   * @brief Interface finite element with displacement-jump kinematics at finite strain.
   *
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   */
  template < int nDim, int nNodes >
  using InterfaceFiniteElementFiniteStrain = InterfaceFiniteElementBase<
    nDim,
    nNodes,
    FiniteStrainDisplacementJumpFormulation< nDim, nNodes > >;

} // namespace Marmot::Elements
