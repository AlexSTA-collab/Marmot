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

/**
 * @file InterfaceFiniteElement.h
 * @brief Interface finite element formulation for displacement-jump based interface mechanics.
 *
 * This file defines the formulation `Marmot::Elements::DisplacementJumpFormulation` and the element
 * `Marmot::Elements::InterfaceFiniteElement`, which is the common base of the interface elements
 * (InterfaceFiniteElementBase.h) with this formulation. The element evaluates interface traction-like
 * quantities from displacement jumps and average surface kinematics and assembles
 * residual/tangent contributions at quadrature points.
 */
#pragma once

#include "Marmot/InterfaceFiniteElementBase.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotGradedInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotStateVarVectorManager.h"

#include <Eigen/Dense>

namespace Marmot::Elements {

  /**
   * @struct DisplacementJumpFormulation
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @tparam MaterialType Interface material class of the element, which provides the nested types `State`,
   * `Tangents`, `Deformation` and `TimeIncrement` of MarmotInterfaceMaterialHypoElastic.
   * @brief Formulation of the interface element with displacement-jump kinematics.
   *
   * The interface material (`MarmotInterfaceMaterialHypoElastic`) receives the increments of the displacement
   * jump and of the surface gradient of both interface sides and returns the traction and the surface stress
   * together with the tangents Q, Z, H and Y. The formulation uses linearized, small-deformation kinematics.
   * The residual and the tangent are assembled with the jump operator and the operator of the average surface
   * gradient of the quadrature point.
   */
  template < int nDim, int nNodes, class MaterialType = MarmotInterfaceMaterialHypoElastic >
  struct DisplacementJumpFormulation {

    static constexpr int nTensor = nDim * nDim;

    using ForceSized                = Eigen::Matrix< double, nDim, 1 >;
    using SurfaceStressSized        = Eigen::Matrix< double, nTensor, 1 >;
    using InterfaceDisplSized       = Eigen::Matrix< double, 2 * nDim, 1 >;
    using InterfaceSurfaceGradSized = Eigen::Matrix< double, 2 * nTensor, 1 >;

    using QMatrixSized = Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor >;
    using ZMatrixSized = Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor >;
    using HMatrixSized = Eigen::Matrix< double, nDim, nTensor, Eigen::RowMajor >;
    using YMatrixSized = Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor >;

    using Material = MaterialType; ///< interface material of the element

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
     * @brief Named state-variable manager for interface quadrature points.
     */
    class QPStateVarManager : public MarmotStateVarVectorManager {

      /*
       * Persistent state layout for accumulated force, surface stress,
       * displacement, surface strain, and material state variables.
       *
       * The displacement and surface strain entries store the accumulated
       * top/bottom quantities:
       *
       *   displacement   = [u_top, u_bottom]
       *   surface strain = [grad_s u_top, grad_s u_bottom]
       */
      inline const static auto layout = makeLayout( {
        { .name = "force", .length = nDim },
        { .name = "alignment padding", .length = nDim % 2 },
        { .name = "surface stress", .length = nDim * nDim },
        { .name = "displacement", .length = 2 * nDim },
        { .name = "surface strain", .length = 2 * nDim * nDim },
        // For nDim == 3, the material state starts after 40 entries.
        { .name   = "state block alignment padding",
          .length = ( 4 - ( ( nDim + ( nDim % 2 ) + nDim * nDim + 2 * nDim + 2 * nDim * nDim ) % 4 ) ) % 4 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      Eigen::Map< ForceSized >                force;
      Eigen::Map< SurfaceStressSized >        surfaceStress;
      Eigen::Map< InterfaceDisplSized >       displacement;
      Eigen::Map< InterfaceSurfaceGradSized > surfaceStrain;

      Eigen::Map< Eigen::VectorXd > materialStateVars;

      /**
       * @brief Number of state variables of the quadrature point itself, without the material state variables.
       *
       * @return the number of state variables up to the begin of the material state
       */
      static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

      /**
       * @brief Map the state variables of one quadrature point.
       *
       * @param[in] theStateVarVector first state variable of the quadrature point
       * @param[in] nStateVars number of state variables of the quadrature point, including the material state
       */
      QPStateVarManager( double* theStateVarVector, int nStateVars )
        : MarmotStateVarVectorManager( theStateVarVector, layout ),
          force( &find( "force" ) ),
          surfaceStress( &find( "surface stress" ) ),
          displacement( &find( "displacement" ) ),
          surfaceStrain( &find( "surface strain" ) ),
          materialStateVars( &find( "begin of material state" ),
                             nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
      {
      }
    };

    /**
     * @brief Evaluate the material at a quadrature point and add its contributions to the element.
     *
     * Updates the state variables of the quadrature point and adds the contributions of the quadrature point to
     * the residual and to the tangent of the element.
     *
     * @tparam Element type of the element
     * @tparam QuadraturePoint type of the quadrature point
     * @tparam PeType type of the element residual
     * @tparam KeType type of the element tangent
     * @param[in] element the element, which holds the material
     * @param[in,out] qp quadrature point, whose state variables are updated
     * @param[in] QTotal total nodal displacements of the element, which this formulation does not need
     * @param[in] dU_GPs increment of the displacements of both interface sides at the quadrature point
     * @param[in] dSurface_strain_GPs increment of the surface gradient of both interface sides at the quadrature
     * point
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
      (void)QTotal;

      const auto& Njump    = qp.NmatJump;
      const auto& Bavg     = qp.BmatAverage;
      auto&       material = element.material;

      ForceSized         force          = qp.managedStateVars->force;
      SurfaceStressSized surface_stress = qp.managedStateVars->surfaceStress;

      QMatrixSized Q_ij;
      ZMatrixSized Z_ijkl;
      HMatrixSized H_ijk;
      YMatrixSized Y_ijkl;

      Q_ij.setZero();
      Z_ijkl.setZero();
      H_ijk.setZero();
      Y_ijkl.setZero();

      if constexpr ( nDim == 3 ) {
        typename Material::State         materialState{ typename Material::Tensor3d( force.data() ),
                                                typename Material::Tensor33d( surface_stress.data() ),
                                                qp.managedStateVars->materialStateVars.data() };
        typename Material::Tangents      materialTangents{ typename Material::Tensor33d( 0.0 ),
                                                      typename Material::Tensor3333d( 0.0 ),
                                                      typename Material::Tensor333d( 0.0 ),
                                                      typename Material::Tensor3333d( 0.0 ) };
        typename Material::Deformation   materialDeformation{ typename Material::Tensor6d( dU_GPs.data() ),
                                                            typename Material::Tensor18d( dSurface_strain_GPs.data() ),
                                                            typename Material::Tensor3d( qp.normal.data() ) };
        typename Material::TimeIncrement materialTimeIncrement{ time, dT };

        material->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );

        force          = mapEigenToFastor( materialState.force );
        surface_stress = Eigen::Map< const SurfaceStressSized >( materialState.surfaceStress.data() );
        Q_ij           = mapEigenToFastor( materialTangents.Q_ij );
        H_ijk          = mapEigenToFastor( materialTangents.H_ijk );
        Z_ijkl         = mapEigenToFastor( materialTangents.Z_ijkl );
        Y_ijkl         = mapEigenToFastor( materialTangents.Y_ijkl );
      }
      else if constexpr ( nDim == 2 ) {
        using namespace Marmot;

        /*
         * MarmotInterfaceMaterialHypoElastic is fixed to 3D Fastor tensor types, so a 2D element has
         * to embed its quadrature-point state into the 3D layout and pull the response back out.
         * expandTo3D / reduceTo2D do exactly that, and both place the 2D block at indices [0,2),
         * which is the layout the material's tensors expect.
         */
        Fastor::Tensor< double, 3 >    force3d         = expandTo3D( Fastor::Tensor< double, 2 >( force.data() ) );
        Fastor::Tensor< double, 3 >    normal3d        = expandTo3D( Fastor::Tensor< double, 2 >( qp.normal.data() ) );
        Fastor::Tensor< double, 3, 3 > surfaceStress3d = expandTo3D(
          Fastor::Tensor< double, 2, 2 >( surface_stress.data() ) );

        const Eigen::Matrix< double, 6, 1 >  dU3d             = expandSidesTo3D< nDim, 1 >( dU_GPs );
        const Eigen::Matrix< double, 18, 1 > dSurfaceStrain3d = expandSidesTo3D< nDim, 2 >( dSurface_strain_GPs );

        typename Material::State         materialState{ force3d,
                                                surfaceStress3d,
                                                qp.managedStateVars->materialStateVars.data() };
        typename Material::Tangents      materialTangents{ typename Material::Tensor33d( 0.0 ),
                                                      typename Material::Tensor3333d( 0.0 ),
                                                      typename Material::Tensor333d( 0.0 ),
                                                      typename Material::Tensor3333d( 0.0 ) };
        typename Material::Deformation   materialDeformation{ typename Material::Tensor6d( dU3d.data() ),
                                                            typename Material::Tensor18d( dSurfaceStrain3d.data() ),
                                                            normal3d };
        typename Material::TimeIncrement materialTimeIncrement{ time, dT };

        material->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );

        const Fastor::Tensor< double, 2 >          force2d         = reduceTo2D< U >( materialState.force );
        const Fastor::Tensor< double, 2, 2 >       surfaceStress2d = reduceTo2D< U, U >( materialState.surfaceStress );
        const Fastor::Tensor< double, 2, 2 >       Q2d             = reduceTo2D< U, U >( materialTangents.Q_ij );
        const Fastor::Tensor< double, 2, 2, 2 >    H2d             = reduceTo2D< U, U, U >( materialTangents.H_ijk );
        const Fastor::Tensor< double, 2, 2, 2, 2 > Z2d = reduceTo2D< U, U, U, U >( materialTangents.Z_ijkl );
        const Fastor::Tensor< double, 2, 2, 2, 2 > Y2d = reduceTo2D< U, U, U, U >( materialTangents.Y_ijkl );

        force          = mapEigenToFastor( force2d );
        surface_stress = Eigen::Map< const SurfaceStressSized >( surfaceStress2d.data() );
        Q_ij           = mapEigenToFastor( Q2d );
        H_ijk          = mapEigenToFastor( H2d );
        Z_ijkl         = mapEigenToFastor( Z2d );
        Y_ijkl         = mapEigenToFastor( Y2d );
      }

      qp.managedStateVars->force         = force;
      qp.managedStateVars->surfaceStress = surface_stress;
      qp.managedStateVars->displacement += dU_GPs;
      qp.managedStateVars->surfaceStrain += dSurface_strain_GPs;

      Pe += Njump.transpose() * force * qp.J0xW;
      Pe += Bavg.transpose() * surface_stress * qp.J0xW;

      Ke += ( Njump.transpose() * Q_ij * Njump + Bavg.transpose() * Z_ijkl * Bavg + Bavg.transpose() * Y_ijkl * Bavg +
              Njump.transpose() * H_ijk * Bavg + Bavg.transpose() * H_ijk.transpose() * Njump ) *
            qp.J0xW;
    }
  };

  /**
   * @brief Interface finite element with displacement-jump kinematics.
   *
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   */
  template < int nDim, int nNodes >
  using InterfaceFiniteElement = InterfaceFiniteElementBase< nDim,
                                                             nNodes,
                                                             DisplacementJumpFormulation< nDim, nNodes > >;

  /**
   * @brief Interface finite element with displacement-jump kinematics and a layer whose stiffness varies through
   * the thickness.
   *
   * The element has the kinematics and the assembly of InterfaceFiniteElement. Its material is
   * MarmotGradedInterfaceMaterialHypoElastic, which integrates the stiffness through the thickness and condenses
   * the first-order displacement of the microproblem of the layer.
   *
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   */
  template < int nDim, int nNodes >
  using GradedInterfaceFiniteElement = InterfaceFiniteElementBase<
    nDim,
    nNodes,
    DisplacementJumpFormulation< nDim, nNodes, MarmotGradedInterfaceMaterialHypoElastic > >;

} // namespace Marmot::Elements
