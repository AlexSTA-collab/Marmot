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
 * @file GYInterfaceFiniteElementFiniteStrain.h
 * @brief Gradient-enhanced equilibrated two-sided interface finite element at finite strain.
 *
 * This file defines the formulation `Marmot::Elements::GradientEnhancedEquilibratedTwoSidedFormulation` and the
 * element `Marmot::Elements::GYInterfaceFiniteElementFiniteStrain`, which is the common base of the interface
 * elements (InterfaceFiniteElementBase.h) with this formulation. In addition to the displacements, the element
 * carries a nonlocal field, with one degree of freedom on each node of the bottom side, which is the midsurface of the
 * layer.
 */
#pragma once

#include "Marmot/InterfaceFiniteElementBase.h"
#include "Marmot/MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/YInterfaceFiniteElement.h"

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace Marmot::Elements {

  /**
   * @brief Take the entries of a vector over the generalized strain of the material in 3D to the dimension of the
   * element.
   *
   * The vector holds the jump, the top and the bottom surface gradient of the material, as a vector of 21 entries.
   *
   * @param[in] full the vector over the generalized strain in 3D
   * @tparam nDim spatial dimension of the element
   * @return the vector over the generalized strain of the element
   */
  template < int nDim >
  Eigen::Matrix< double, nDim + 2 * nDim * nDim, 1 > reduceGeneralizedStrain( const Fastor::Tensor< double, 21 >& full )
  {
    using namespace Marmot;

    if constexpr ( nDim == 3 ) {
      return Eigen::Map< const Eigen::Matrix< double, 21, 1 > >( full.data() );
    }
    else {
      const Fastor::Tensor< double, 3 >    w( full.data() );
      const Fastor::Tensor< double, 3, 3 > Ap( full.data() + 3 );
      const Fastor::Tensor< double, 3, 3 > Am( full.data() + 12 );

      const Fastor::Tensor< double, 2 >    w2  = reduceTo2D< U >( w );
      const Fastor::Tensor< double, 2, 2 > Ap2 = reduceTo2D< U, U >( Ap );
      const Fastor::Tensor< double, 2, 2 > Am2 = reduceTo2D< U, U >( Am );

      Eigen::Matrix< double, nDim + 2 * nDim * nDim, 1 > reduced;
      reduced.template segment< nDim >( 0 ) = Eigen::Map< const Eigen::Matrix< double, nDim, 1 > >( w2.data() );
      reduced.template segment< nDim * nDim >( nDim ) = Eigen::Map< const Eigen::Matrix< double, nDim * nDim, 1 > >(
        Ap2.data() );
      reduced.template segment< nDim * nDim >(
        nDim + nDim * nDim ) = Eigen::Map< const Eigen::Matrix< double, nDim * nDim, 1 > >( Am2.data() );
      return reduced;
    }
  }

  /**
   * @struct GradientEnhancedEquilibratedTwoSidedFormulation
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @brief Formulation of the interface element with equilibrated two-sided kinematics, at finite strain, with a
   * nonlocal field.
   *
   * The interface material (`MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain`) receives the total
   * displacements and surface gradients of both interface sides, the connector between the faces, and the nonlocal
   * field at the quadrature point. It returns the generalized stress, the local driving field and the nonlocal
   * radius, together with the condensed tangents.
   *
   * The nonlocal field is shared by both faces, which is exact for a nonlocal radius that exceeds the thickness of the
   * layer. Its equation is reduced to the midsurface,
   *
   * @f[
   *   h \left[ \mathbf{N}^\mathsf{T} ( \bar{N} - \bar{L} ) + l^2 \nabla\mathbf{N}^\mathsf{T} \nabla\bar{N} \right] = 0,
   * @f]
   *
   * where @f$ h @f$ is the thickness of the layer, @f$ \bar{L} @f$ the average local field of the two faces and
   * @f$ l @f$ the nonlocal radius. The degrees of freedom of the nonlocal field are placed after the displacements in
   * the vectors of the element.
   */
  template < int nDim, int nNodes >
  struct GradientEnhancedEquilibratedTwoSidedFormulation {

    static constexpr int nTensor         = nDim * nDim;
    static constexpr int nInterfaceNodes = nNodes / 2;
    static constexpr int nDofD           = nNodes * nDim;      ///< degrees of freedom of the displacements
    static constexpr int nDofA           = nInterfaceNodes;    ///< degrees of freedom of the nonlocal field
    static constexpr int nTotalStrain    = nDim + 2 * nTensor; ///< size of the generalized strain (w, A+, A-)

    /// number of degrees of freedom in addition to the displacements
    static constexpr int nExtraDofs = nDofA;

    using ForceSized                = Eigen::Matrix< double, nDim, 1 >;
    using SurfaceStressSized        = Eigen::Matrix< double, nTensor, 1 >;
    using InterfaceDisplSized       = Eigen::Matrix< double, 2 * nDim, 1 >;
    using InterfaceSurfaceGradSized = Eigen::Matrix< double, 2 * nTensor, 1 >;

    using Material = MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain; ///< interface material

    /**
     * @brief Named state-variable manager for interface quadrature points.
     */
    class QPStateVarManager : public MarmotStateVarVectorManager {

      /*
       * Persistent state layout for the generalized force, the two generalized surface stresses, the total
       * displacement, the total surface strain, the local and the nonlocal field, and the material state variables.
       *
       *   displacement   = [u_top, u_bottom]
       *   surface strain = [grad_s u_top, grad_s u_bottom]
       */
      inline const static auto layout = makeLayout( {
        { .name = "generalizedForce", .length = nDim },
        { .name = "surfaceStressPlus", .length = nTensor },
        { .name = "surfaceStressMinus", .length = nTensor },
        { .name = "displacement", .length = 2 * nDim },
        { .name = "surfaceStrain", .length = 2 * nTensor },
        { .name = "localField", .length = 1 },
        { .name = "nonlocalField", .length = 1 },
        // The material state starts at a multiple of four entries.
        { .name   = "state block alignment padding",
          .length = ( 4 - ( ( nDim + 2 * nTensor + 2 * nDim + 2 * nTensor + 2 ) % 4 ) ) % 4 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      Eigen::Map< ForceSized >                generalizedForce;
      Eigen::Map< SurfaceStressSized >        surfaceStressPlus;
      Eigen::Map< SurfaceStressSized >        surfaceStressMinus;
      Eigen::Map< InterfaceDisplSized >       displacement;
      Eigen::Map< InterfaceSurfaceGradSized > surfaceStrain;
      double&                                 localField;
      double&                                 nonlocalField;

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
          generalizedForce( &find( "generalizedForce" ) ),
          surfaceStressPlus( &find( "surfaceStressPlus" ) ),
          surfaceStressMinus( &find( "surfaceStressMinus" ) ),
          displacement( &find( "displacement" ) ),
          surfaceStrain( &find( "surfaceStrain" ) ),
          localField( find( "localField" ) ),
          nonlocalField( find( "nonlocalField" ) ),
          materialStateVars( &find( "begin of material state" ),
                             nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
      {
      }
    };

    /// This formulation keeps the connector between the faces at the quadrature points.
    struct QuadraturePointData : TwoSidedGeometryData< nDim > {};

    /**
     * @brief Complete the geometry of a quadrature point with the connector between the two faces.
     *
     * @tparam Element type of the element
     * @tparam QuadraturePoint type of the quadrature point
     * @param[in] element the element, which holds the coordinates of the nodes
     * @param[in,out] qp quadrature point, whose normal and projections may be reoriented
     *
     * @throws std::invalid_argument if the connector has no positive component along the normal
     */
    template < class Element, class QuadraturePoint >
    static void initializeQuadraturePoint( const Element& element, QuadraturePoint& qp )
    {
      initializeTwoSidedQuadraturePoint( element, qp, "GYInterfaceFiniteElementFiniteStrain" );
    }

    /**
     * @brief Fields on the nodes of the element: the displacement on all nodes, and the nonlocal field on the nodes of
     * the bottom side.
     *
     * @return the names of the fields of every node
     */
    static std::vector< std::vector< std::string > > nodeFields()
    {
      std::vector< std::vector< std::string > > fields( nNodes );

      for ( int i = 0; i < nNodes; ++i ) {
        fields[i].push_back( "displacement" );
        if ( i < nInterfaceNodes )
          fields[i].push_back( "nonlocal damage" );
      }
      return fields;
    }

    /**
     * @brief Permutation from the order of the degrees of freedom of the element to the order node by node.
     *
     * @return for every degree of freedom of the element its index in the order node by node
     */
    static std::vector< int > dofIndicesPermutationPattern()
    {
      std::vector< int > permutation( nDofD + nDofA );

      // a node of the bottom side carries [ u(nDim), abar ] and takes nDim + 1 places, the top nodes carry [ u(nDim) ]
      // after all of the bottom ones
      for ( int A = 0; A < nInterfaceNodes; ++A ) {
        for ( int c = 0; c < nDim; ++c ) {
          permutation[nDim * A + c]                          = ( nDim + 1 ) * A + c;
          permutation[nDim * nInterfaceNodes + nDim * A + c] = ( nDim + 1 ) * nInterfaceNodes + nDim * A + c;
        }
        permutation[nDofD + A] = ( nDim + 1 ) * A + nDim;
      }
      return permutation;
    }

    /**
     * @brief Evaluate the material at a quadrature point and add its contributions to the element.
     *
     * @tparam Element type of the element
     * @tparam QuadraturePoint type of the quadrature point
     * @tparam QType type of the total nodal values
     * @tparam PeType type of the element residual
     * @tparam KeType type of the element tangent
     * @param[in] element the element, which holds the material
     * @param[in,out] qp quadrature point, whose state variables are updated
     * @param[in] QTotal total nodal displacements and nonlocal field of the element at the end of the increment
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

      auto& material = element.material;
      auto& data     = qp.formulationData;

      const InterfaceDisplSized       totalU_GPs               = qp.interpolateToSides( QTotal );
      const InterfaceSurfaceGradSized totalSurfaceGradient_GPs = qp.surfaceGradientsAtSides( QTotal );

      // The nonlocal field at this point, shared by both faces: a nonlocal radius above the thickness makes it flat
      // across the layer, so that there is one value here and no jump of it to carry.
      const Eigen::Matrix< double, nDofA, 1 > nodalNonlocal = QTotal.template segment< nDofA >( nDofD );
      const double                            nonlocalField = ( qp.N * nodalNonlocal )( 0 );

      Material::State state{ { embedVectorIn3D< nDim >( qp.managedStateVars->generalizedForce.data() ),
                               embedTensorIn3D< nDim >( qp.managedStateVars->surfaceStressPlus.data() ),
                               embedTensorIn3D< nDim >( qp.managedStateVars->surfaceStressMinus.data() ),
                               qp.managedStateVars->materialStateVars.data() },
                             0.0,
                             0.0 };

      Material::Tangents tangents{ { Material::Tensor33d( 0.0 ),
                                     Material::Tensor333d( 0.0 ),
                                     Material::Tensor333d( 0.0 ),
                                     Material::Tensor333d( 0.0 ),
                                     Material::Tensor3333d( 0.0 ),
                                     Material::Tensor3333d( 0.0 ),
                                     Material::Tensor333d( 0.0 ),
                                     Material::Tensor3333d( 0.0 ),
                                     Material::Tensor3333d( 0.0 ) },
                                   Material::Tensor21d( 0.0 ),
                                   Material::Tensor21d( 0.0 ),
                                   0.0 };

      // the materials are fixed to three dimensions, a two-dimensional element embeds its kinematics
      const Eigen::Matrix< double, 6, 1 >  U3d               = expandSidesTo3D< nDim, 1 >( totalU_GPs );
      const Eigen::Matrix< double, 18, 1 > surfaceGradient3d = expandSidesTo3D< nDim, 2 >( totalSurfaceGradient_GPs );

      Material::Deformation   deformation{ { Material::Tensor6d( U3d.data() ),
                                             Material::Tensor18d( surfaceGradient3d.data() ),
                                             embedVectorIn3D< nDim >( qp.normal.data() ),
                                             embedVectorIn3D< nDim >( data.separationVector.data() ) },
                                         nonlocalField };
      Material::TimeIncrement timeIncrement{ time, dT };

      material->computeStress( state, tangents, deformation, timeIncrement );

      qp.managedStateVars->displacement  = totalU_GPs;
      qp.managedStateVars->surfaceStrain = totalSurfaceGradient_GPs;
      qp.managedStateVars->localField    = state.localField;
      qp.managedStateVars->nonlocalField = nonlocalField;

      // momentum: the stress and the tangent of the displacements
      assembleTwoSidedResponse< nDim, nNodes >( qp,
                                                state,
                                                tangents,
                                                Pe.template head< nDofD >(),
                                                Ke.template topLeftCorner< nDofD, nDofD >() );

      const auto Bx = generalizedStrainOperator< nDim, nNodes >( qp );

      const Eigen::Matrix< double, nTotalStrain, 1 > dStress_dNonlocal = reduceGeneralizedStrain< nDim >(
        tangents.dGeneralizedStress_dNonlocal );
      const Eigen::Matrix< double, 1, nTotalStrain >
        dLocalField_dStrain = reduceGeneralizedStrain< nDim >( tangents.dLocalField_dGeneralizedStrain ).transpose();

      Ke.template block< nDofD, nDofA >( 0, nDofD ) += Bx.transpose() * dStress_dNonlocal * qp.N * qp.J0xW;

      // The nonlocal equation, reduced to the midsurface,
      //
      //   h [ N^T ( abar - Lbar ) + l^2 gradN^T gradN abar ] ,
      //
      // the Helmholtz equation of the layer without the gradient through the thickness: abar is flat across the layer
      // for a nonlocal radius above the thickness, and the faces pass no flux into the adjacent bodies.
      const double h = material->thickness();
      const double l = state.nonlocalRadius;

      Pe.template tail< nDofA >() += h *
                                     ( qp.N.transpose() * ( nonlocalField - state.localField ) +
                                       l * l * qp.gradN.transpose() * ( qp.gradN * nodalNonlocal ) ) *
                                     qp.J0xW;

      Ke.template block< nDofA, nDofD >( nDofD, 0 ) -= h * qp.N.transpose() * dLocalField_dStrain * Bx * qp.J0xW;

      Ke.template block< nDofA, nDofA >( nDofD, nDofD ) += h *
                                                           ( qp.N.transpose() * qp.N *
                                                               ( 1.0 - tangents.dLocalField_dNonlocal ) +
                                                             l * l * qp.gradN.transpose() * qp.gradN ) *
                                                           qp.J0xW;
    }
  };

  /**
   * @struct GradientEnhancedFaceFieldsFormulation
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @brief Formulation of the equilibrated two-sided interface element at finite strain with a nonlocal field of each
   * face.
   *
   * In contrast to GradientEnhancedEquilibratedTwoSidedFormulation, the faces of the layer have independent nonlocal
   * fields @f$ \bar N^\pm @f$, and so independent damage. The top nodes carry the field of the top face and the
   * bottom nodes the field of the bottom face, after the displacements in the vectors of the element, the bottom
   * side first. The field varies linearly through the thickness @f$ h @f$ of the layer. Each face has its own
   * nonlocal equation on the midsurface, with half of the thickness of the layer, and the derivative through the
   * thickness couples the faces,
   *
   * @f[
   *   \frac{h}{2} \left[ \mathbf{N}^\mathsf{T} ( \bar N^\pm - L^\pm ) + (l^\pm)^2 \nabla\mathbf{N}^\mathsf{T}
   *   \nabla\bar N^\pm \right] \pm \frac{(l^\pm)^2}{h}\, \mathbf{N}^\mathsf{T} ( \bar N^+ - \bar N^- ) = 0,
   * @f]
   *
   * where @f$ L^\pm @f$ is the local driving field of the face. The coupling comes from the derivative of the field
   * through the thickness for the linear profile, and each face has the coefficient of its own nonlocal radius, so
   * that the radii of the faces may differ. It vanishes for a small radius, where the faces are independent, and forces
   * @f$ \bar N^+ = \bar N^- @f$ for a radius that is large against the thickness.
   */
  template < int nDim, int nNodes >
  struct GradientEnhancedFaceFieldsFormulation {

    static constexpr int nTensor         = nDim * nDim;
    static constexpr int nInterfaceNodes = nNodes / 2;
    static constexpr int nDofD           = nNodes * nDim;      ///< degrees of freedom of the displacements
    static constexpr int nDofA           = nInterfaceNodes;    ///< degrees of freedom of the nonlocal field of a face
    static constexpr int nTotalStrain    = nDim + 2 * nTensor; ///< size of the generalized strain (w, A+, A-)

    /// number of degrees of freedom in addition to the displacements: a nonlocal field on every node
    static constexpr int nExtraDofs = 2 * nDofA;

    using ForceSized                = Eigen::Matrix< double, nDim, 1 >;
    using SurfaceStressSized        = Eigen::Matrix< double, nTensor, 1 >;
    using InterfaceDisplSized       = Eigen::Matrix< double, 2 * nDim, 1 >;
    using InterfaceSurfaceGradSized = Eigen::Matrix< double, 2 * nTensor, 1 >;

    using Material = MarmotEquilibratedGradientEnhancedInterfaceMaterialFiniteStrain; ///< interface material

    /**
     * @brief Named state-variable manager for interface quadrature points.
     */
    class QPStateVarManager : public MarmotStateVarVectorManager {

      /*
       * Persistent state layout for the generalized force, the two generalized surface stresses, the total
       * displacement, the total surface strain, the local and the nonlocal field of each face, and the material state
       * variables.
       */
      inline const static auto layout = makeLayout( {
        { .name = "generalizedForce", .length = nDim },
        { .name = "surfaceStressPlus", .length = nTensor },
        { .name = "surfaceStressMinus", .length = nTensor },
        { .name = "displacement", .length = 2 * nDim },
        { .name = "surfaceStrain", .length = 2 * nTensor },
        { .name = "localFieldPlus", .length = 1 },
        { .name = "localFieldMinus", .length = 1 },
        { .name = "nonlocalFieldPlus", .length = 1 },
        { .name = "nonlocalFieldMinus", .length = 1 },
        // The material state starts at a multiple of four entries.
        { .name   = "state block alignment padding",
          .length = ( 4 - ( ( nDim + 2 * nTensor + 2 * nDim + 2 * nTensor + 4 ) % 4 ) ) % 4 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      Eigen::Map< ForceSized >                generalizedForce;
      Eigen::Map< SurfaceStressSized >        surfaceStressPlus;
      Eigen::Map< SurfaceStressSized >        surfaceStressMinus;
      Eigen::Map< InterfaceDisplSized >       displacement;
      Eigen::Map< InterfaceSurfaceGradSized > surfaceStrain;
      double&                                 localFieldPlus;
      double&                                 localFieldMinus;
      double&                                 nonlocalFieldPlus;
      double&                                 nonlocalFieldMinus;

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
          generalizedForce( &find( "generalizedForce" ) ),
          surfaceStressPlus( &find( "surfaceStressPlus" ) ),
          surfaceStressMinus( &find( "surfaceStressMinus" ) ),
          displacement( &find( "displacement" ) ),
          surfaceStrain( &find( "surfaceStrain" ) ),
          localFieldPlus( find( "localFieldPlus" ) ),
          localFieldMinus( find( "localFieldMinus" ) ),
          nonlocalFieldPlus( find( "nonlocalFieldPlus" ) ),
          nonlocalFieldMinus( find( "nonlocalFieldMinus" ) ),
          materialStateVars( &find( "begin of material state" ),
                             nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
      {
      }
    };

    /// This formulation keeps the connector between the faces at the quadrature points.
    struct QuadraturePointData : TwoSidedGeometryData< nDim > {};

    /**
     * @brief Complete the geometry of a quadrature point with the connector between the two faces.
     *
     * @tparam Element type of the element
     * @tparam QuadraturePoint type of the quadrature point
     * @param[in] element the element, which holds the coordinates of the nodes
     * @param[in,out] qp quadrature point, whose normal and projections may be reoriented
     *
     * @throws std::invalid_argument if the connector has no positive component along the normal
     */
    template < class Element, class QuadraturePoint >
    static void initializeQuadraturePoint( const Element& element, QuadraturePoint& qp )
    {
      initializeTwoSidedQuadraturePoint( element, qp, "GY2InterfaceFiniteElementFiniteStrain" );
    }

    /**
     * @brief Fields on the nodes of the element: the displacement and the nonlocal field of its face on every node.
     *
     * @return the names of the fields of every node
     */
    static std::vector< std::vector< std::string > > nodeFields()
    {
      return std::vector< std::vector< std::string > >( nNodes, { "displacement", "nonlocal damage" } );
    }

    /**
     * @brief Permutation from the order of the degrees of freedom of the element to the order node by node.
     *
     * @return for every degree of freedom of the element its index in the order node by node
     */
    static std::vector< int > dofIndicesPermutationPattern()
    {
      std::vector< int > permutation( nDofD + nExtraDofs );

      // every node carries [ u(nDim), N ] and takes nDim + 1 places, the bottom nodes before the top nodes
      for ( int A = 0; A < nInterfaceNodes; ++A ) {
        for ( int c = 0; c < nDim; ++c ) {
          permutation[nDim * A + c]                          = ( nDim + 1 ) * A + c;
          permutation[nDim * nInterfaceNodes + nDim * A + c] = ( nDim + 1 ) * ( nInterfaceNodes + A ) + c;
        }
        permutation[nDofD + A]                   = ( nDim + 1 ) * A + nDim;
        permutation[nDofD + nInterfaceNodes + A] = ( nDim + 1 ) * ( nInterfaceNodes + A ) + nDim;
      }
      return permutation;
    }

    /**
     * @brief Evaluate the material at a quadrature point and add its contributions to the element.
     *
     * @tparam Element type of the element
     * @tparam QuadraturePoint type of the quadrature point
     * @tparam QType type of the total nodal values
     * @tparam PeType type of the element residual
     * @tparam KeType type of the element tangent
     * @param[in] element the element, which holds the material
     * @param[in,out] qp quadrature point, whose state variables are updated
     * @param[in] QTotal total nodal displacements and nonlocal fields of the element at the end of the increment
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

      auto& material = element.material;
      auto& data     = qp.formulationData;

      const InterfaceDisplSized       totalU_GPs               = qp.interpolateToSides( QTotal );
      const InterfaceSurfaceGradSized totalSurfaceGradient_GPs = qp.surfaceGradientsAtSides( QTotal );

      // The nonlocal fields of the faces at this point: the bottom nodes carry the field of the bottom face, the top
      // nodes the one of the top face.
      constexpr int offBottom = nDofD;
      constexpr int offTop    = nDofD + nDofA;

      const Eigen::Matrix< double, nDofA, 1 > nodalNonlocalBottom = QTotal.template segment< nDofA >( offBottom );
      const Eigen::Matrix< double, nDofA, 1 > nodalNonlocalTop    = QTotal.template segment< nDofA >( offTop );
      const double                            nonlocalBottom      = ( qp.N * nodalNonlocalBottom )( 0 );
      const double                            nonlocalTop         = ( qp.N * nodalNonlocalTop )( 0 );

      Material::FaceFieldsState state{ { embedVectorIn3D< nDim >( qp.managedStateVars->generalizedForce.data() ),
                                         embedTensorIn3D< nDim >( qp.managedStateVars->surfaceStressPlus.data() ),
                                         embedTensorIn3D< nDim >( qp.managedStateVars->surfaceStressMinus.data() ),
                                         qp.managedStateVars->materialStateVars.data() },
                                       0.0,
                                       0.0,
                                       0.0,
                                       0.0 };

      Material::FaceFieldsTangents tangents{ { Material::Tensor33d( 0.0 ),
                                               Material::Tensor333d( 0.0 ),
                                               Material::Tensor333d( 0.0 ),
                                               Material::Tensor333d( 0.0 ),
                                               Material::Tensor3333d( 0.0 ),
                                               Material::Tensor3333d( 0.0 ),
                                               Material::Tensor333d( 0.0 ),
                                               Material::Tensor3333d( 0.0 ),
                                               Material::Tensor3333d( 0.0 ) },
                                             Material::Tensor21d( 0.0 ),
                                             Material::Tensor21d( 0.0 ),
                                             Material::Tensor21d( 0.0 ),
                                             Material::Tensor21d( 0.0 ),
                                             0.0,
                                             0.0,
                                             0.0,
                                             0.0 };

      // the materials are fixed to three dimensions, a two-dimensional element embeds its kinematics
      const Eigen::Matrix< double, 6, 1 >  U3d               = expandSidesTo3D< nDim, 1 >( totalU_GPs );
      const Eigen::Matrix< double, 18, 1 > surfaceGradient3d = expandSidesTo3D< nDim, 2 >( totalSurfaceGradient_GPs );

      Material::FaceFieldsDeformation deformation{ { Material::Tensor6d( U3d.data() ),
                                                     Material::Tensor18d( surfaceGradient3d.data() ),
                                                     embedVectorIn3D< nDim >( qp.normal.data() ),
                                                     embedVectorIn3D< nDim >( data.separationVector.data() ) },
                                                   nonlocalTop,
                                                   nonlocalBottom };
      Material::TimeIncrement         timeIncrement{ time, dT };

      material->computeStress( state, tangents, deformation, timeIncrement );

      qp.managedStateVars->displacement       = totalU_GPs;
      qp.managedStateVars->surfaceStrain      = totalSurfaceGradient_GPs;
      qp.managedStateVars->localFieldPlus     = state.localFieldTop;
      qp.managedStateVars->localFieldMinus    = state.localFieldBottom;
      qp.managedStateVars->nonlocalFieldPlus  = nonlocalTop;
      qp.managedStateVars->nonlocalFieldMinus = nonlocalBottom;

      // momentum: the stress and the tangent of the displacements
      assembleTwoSidedResponse< nDim, nNodes >( qp,
                                                state,
                                                tangents,
                                                Pe.template head< nDofD >(),
                                                Ke.template topLeftCorner< nDofD, nDofD >() );

      const auto Bx = generalizedStrainOperator< nDim, nNodes >( qp );

      const Eigen::Matrix< double, nTotalStrain, 1 > dStress_dNonlocalTop = reduceGeneralizedStrain< nDim >(
        tangents.dGeneralizedStress_dNonlocalTop );
      const Eigen::Matrix< double, nTotalStrain, 1 > dStress_dNonlocalBottom = reduceGeneralizedStrain< nDim >(
        tangents.dGeneralizedStress_dNonlocalBottom );
      const Eigen::Matrix< double, 1, nTotalStrain >
        dLocalTop_dStrain = reduceGeneralizedStrain< nDim >( tangents.dLocalFieldTop_dGeneralizedStrain ).transpose();
      const Eigen::Matrix< double, 1, nTotalStrain >
        dLocalBottom_dStrain = reduceGeneralizedStrain< nDim >( tangents.dLocalFieldBottom_dGeneralizedStrain )
                                 .transpose();

      Ke.template block< nDofD, nDofA >( 0, offBottom ) += Bx.transpose() * dStress_dNonlocalBottom * qp.N * qp.J0xW;
      Ke.template block< nDofD, nDofA >( 0, offTop ) += Bx.transpose() * dStress_dNonlocalTop * qp.N * qp.J0xW;

      // The nonlocal equation of a face, reduced to the midsurface,
      //
      //   h/2 [ N^T ( N+- - L+- ) + l+-^2 gradN^T gradN N+- ] ,
      //
      // each face carrying half of the thickness of the layer.
      const double half = 0.5 * material->thickness();

      const double lTop    = state.nonlocalRadiusTop;
      const double lBottom = state.nonlocalRadiusBottom;

      Pe.template segment< nDofA >( offTop ) += half *
                                                ( qp.N.transpose() * ( nonlocalTop - state.localFieldTop ) +
                                                  lTop * lTop * qp.gradN.transpose() *
                                                    ( qp.gradN * nodalNonlocalTop ) ) *
                                                qp.J0xW;
      // The derivative through the thickness couples the faces. Each face has the coefficient of its own nonlocal
      // radius.
      const double couplingTop    = lTop * lTop / ( 2.0 * half );
      const double couplingBottom = lBottom * lBottom / ( 2.0 * half );

      Pe.template segment< nDofA >( offTop ) += couplingTop * qp.N.transpose() * ( nonlocalTop - nonlocalBottom ) *
                                                qp.J0xW;
      Pe.template segment< nDofA >( offBottom ) -= couplingBottom * qp.N.transpose() *
                                                   ( nonlocalTop - nonlocalBottom ) * qp.J0xW;

      Pe.template segment< nDofA >( offBottom ) += half *
                                                   ( qp.N.transpose() * ( nonlocalBottom - state.localFieldBottom ) +
                                                     lBottom * lBottom * qp.gradN.transpose() *
                                                       ( qp.gradN * nodalNonlocalBottom ) ) *
                                                   qp.J0xW;

      // the local fields depend on the displacements and, through the normal-gradient jump, on the field of the
      // other face
      Ke.template block< nDofA, nDofD >( offTop, 0 ) -= half * qp.N.transpose() * dLocalTop_dStrain * Bx * qp.J0xW;
      Ke.template block< nDofA, nDofD >( offBottom, 0 ) -= half * qp.N.transpose() * dLocalBottom_dStrain * Bx *
                                                           qp.J0xW;

      const Eigen::Matrix< double, nDofA, nDofA > NN = qp.N.transpose() * qp.N;

      Ke.template block< nDofA, nDofA >( offTop, offTop ) += half *
                                                             ( NN * ( 1.0 - tangents.dLocalFieldTop_dNonlocalTop ) +
                                                               lTop * lTop * qp.gradN.transpose() * qp.gradN ) *
                                                             qp.J0xW;
      Ke.template block< nDofA, nDofA >( offBottom,
                                         offBottom ) += half *
                                                        ( NN * ( 1.0 - tangents.dLocalFieldBottom_dNonlocalBottom ) +
                                                          lBottom * lBottom * qp.gradN.transpose() * qp.gradN ) *
                                                        qp.J0xW;
      Ke.template block< nDofA, nDofA >( offTop, offTop ) += couplingTop * NN * qp.J0xW;
      Ke.template block< nDofA, nDofA >( offBottom, offBottom ) += couplingBottom * NN * qp.J0xW;
      Ke.template block< nDofA, nDofA >( offTop, offBottom ) -= couplingTop * NN * qp.J0xW;
      Ke.template block< nDofA, nDofA >( offBottom, offTop ) -= couplingBottom * NN * qp.J0xW;

      Ke.template block< nDofA, nDofA >( offTop, offBottom ) -= half * NN * tangents.dLocalFieldTop_dNonlocalBottom *
                                                                qp.J0xW;
      Ke.template block< nDofA, nDofA >( offBottom, offTop ) -= half * NN * tangents.dLocalFieldBottom_dNonlocalTop *
                                                                qp.J0xW;
    }
  };

  /**
   * @brief Gradient-enhanced equilibrated two-sided interface element at finite strain with a nonlocal field of each
   * face.
   *
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   */
  template < int nDim, int nNodes >
  using GY2InterfaceFiniteElementFiniteStrain = InterfaceFiniteElementBase<
    nDim,
    nNodes,
    GradientEnhancedFaceFieldsFormulation< nDim, nNodes > >;

  /**
   * @brief Gradient-enhanced equilibrated two-sided interface element at finite strain.
   *
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   */
  template < int nDim, int nNodes >
  using GYInterfaceFiniteElementFiniteStrain = InterfaceFiniteElementBase<
    nDim,
    nNodes,
    GradientEnhancedEquilibratedTwoSidedFormulation< nDim, nNodes > >;

} // namespace Marmot::Elements
