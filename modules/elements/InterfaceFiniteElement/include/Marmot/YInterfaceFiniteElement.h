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
 * @file YInterfaceFiniteElement.h
 * @brief Equilibrated two-sided interface finite element.
 *
 * This file defines the formulation `Marmot::Elements::EquilibratedTwoSidedFormulation` and the element
 * `Marmot::Elements::YInterfaceFiniteElement`, which is the common base of the interface elements
 * (InterfaceFiniteElementBase.h) with this formulation. The element passes the surface gradients of the top and of
 * the bottom face separately to the equilibrated interface material
 * (MarmotEquilibratedXInterfaceMaterialHypoElastic), which solves a local traction-equilibrium problem and returns
 * a fully populated 3x3 block tangent over the jump and the two surface gradients, including the cross blocks
 * between them. The element assembly therefore uses all nine blocks.
 */
#pragma once

#include "Marmot/InterfaceFiniteElementBase.h"
#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotStateVarVectorManager.h"

#include <Eigen/Dense>

#include <stdexcept>

namespace Marmot::Elements {

  /**
   * @brief Geometry of the paired faces at a quadrature point, which the two-sided formulations keep.
   *
   * @tparam nDim spatial dimension of the element
   */
  template < int nDim >
  struct TwoSidedGeometryData {
    /// connector between the two faces at the quadrature point, from the bottom to the top face
    Eigen::Matrix< double, nDim, 1 > separationVector = Eigen::Matrix< double, nDim, 1 >::Zero();

    /// component of the connector along the normal
    double normalSeparation = 0.0;

    /// component of the connector in the tangent plane
    Eigen::Matrix< double, nDim, 1 > tangentialSeparation = Eigen::Matrix< double, nDim, 1 >::Zero();
  };

  /**
   * @brief Complete the geometry of a quadrature point of a two-sided formulation with the connector of the faces.
   *
   * The connector is the difference of the positions of the top and of the bottom face at the quadrature point. The
   * normal is oriented along the connector, from the bottom to the top face, and the normal and the tangential
   * separation follow from it. Coincident faces have no connector.
   *
   * @tparam Element type of the element
   * @tparam QuadraturePoint type of the quadrature point, whose formulation data derives from TwoSidedGeometryData
   * @param[in] element the element, which holds the coordinates of the nodes
   * @param[in,out] qp quadrature point, whose normal and projections may be reoriented
   * @param[in] formulationName name of the formulation, for the error message
   *
   * @throws std::invalid_argument if the connector has no positive component along the normal
   */
  template < class Element, class QuadraturePoint >
  void initializeTwoSidedQuadraturePoint( const Element& element, QuadraturePoint& qp, const char* formulationName )
  {
    using VectorDim = typename Element::VectorDim;
    using TensorDim = typename Element::TensorDim;

    constexpr double geometryTolerance = 1.0e-12;

    const VectorDim xBottom = qp.NmatSide * element.getSideCoordinates( 0 );
    const VectorDim xTop    = qp.NmatSide * element.getSideCoordinates( 1 );

    auto& data            = qp.formulationData;
    data.separationVector = xTop - xBottom;

    if ( data.separationVector.norm() > geometryTolerance && data.separationVector.dot( qp.normal ) < 0.0 ) {
      qp.normal *= -1.0;
    }

    qp.normalProjection  = qp.normal * qp.normal.transpose();
    qp.tangentProjection = TensorDim::Identity() - qp.normalProjection;

    data.normalSeparation     = data.separationVector.dot( qp.normal );
    data.tangentialSeparation = qp.tangentProjection * data.separationVector;

    if ( data.separationVector.norm() > geometryTolerance && data.normalSeparation <= geometryTolerance ) {
      throw std::invalid_argument( std::string( formulationName ) +
                                   ": paired faces have no positive separation in the interface-normal direction." );
    }
  }

  /**
   * @brief Embed a vector of the element in the 3D layout of the interface materials.
   *
   * @tparam nDim spatial dimension of the element
   * @param[in] data the nDim entries of the vector
   * @return the vector in 3D
   */
  template < int nDim >
  Fastor::Tensor< double, 3 > embedVectorIn3D( const double* data )
  {
    if constexpr ( nDim == 3 )
      return Fastor::Tensor< double, 3 >( data );
    else
      return expandTo3D( Fastor::Tensor< double, nDim >( data ) );
  }

  /**
   * @brief Embed a second-order tensor of the element in the 3D layout of the interface materials.
   *
   * @tparam nDim spatial dimension of the element
   * @param[in] data the nDim x nDim entries of the tensor
   * @return the tensor in 3D
   */
  template < int nDim >
  Fastor::Tensor< double, 3, 3 > embedTensorIn3D( const double* data )
  {
    if constexpr ( nDim == 3 )
      return Fastor::Tensor< double, 3, 3 >( data );
    else
      return expandTo3D( Fastor::Tensor< double, nDim, nDim >( data ) );
  }

  /**
   * @brief Operator that maps the displacements of the element to the generalized strain (w, A+, A-) of a quadrature
   * point.
   *
   * The momentum residual of a quadrature point is the transpose of this operator times the generalized stress, and
   * its tangent is the transpose times the condensed tangent times the operator.
   *
   * @tparam nDim spatial dimension of the element
   * @tparam nNodes number of nodes of the element
   * @tparam QuadraturePoint type of the quadrature point
   * @param[in] qp the quadrature point
   * @return the operator, with the jump, the surface gradient of the top side and the one of the bottom side as rows
   */
  template < int nDim, int nNodes, class QuadraturePoint >
  Eigen::Matrix< double, nDim + 2 * nDim * nDim, nNodes * nDim > generalizedStrainOperator( const QuadraturePoint& qp )
  {
    constexpr int nTensor  = nDim * nDim;
    constexpr int nDofD    = nNodes * nDim;
    constexpr int halfSize = nDofD / 2;

    Eigen::Matrix< double, nDim + 2 * nTensor, nDofD > Bx = Eigen::Matrix< double, nDim + 2 * nTensor, nDofD >::Zero();

    Bx.template block< nDim, nDofD >( 0, 0 )                    = qp.NmatJump;
    Bx.template block< nTensor, halfSize >( nDim, halfSize )    = qp.BmatSide;
    Bx.template block< nTensor, halfSize >( nDim + nTensor, 0 ) = qp.BmatSide;
    return Bx;
  }

  /**
   * @brief Take the response of an equilibrated two-sided interface material, store it in the state of the
   * quadrature point and add the contributions of the quadrature point to the residual and to the tangent.
   *
   * The materials are fixed to three dimensions, so that a two-dimensional element pulls the response back with
   * reduceTo2D. The residual is the internal force, so that the tangent is its derivative.
   *
   * @tparam nDim spatial dimension of the element
   * @tparam nNodes number of nodes of the element
   * @tparam QuadraturePoint type of the quadrature point
   * @tparam MaterialState type of the response of the material
   * @tparam MaterialTangents type of the tangents of the material
   * @tparam PeType type of the element residual
   * @tparam KeType type of the element tangent
   * @param[in,out] qp quadrature point, whose generalized force and surface stresses are updated
   * @param[in] state response of the material
   * @param[in] tangents condensed block tangent of the material
   * @param[in,out] Pe residual of the displacements of the element, to which the contribution is added
   * @param[in,out] Ke tangent of the displacements of the element, to which the contribution is added
   *
   * @details All nine blocks of the tangent are assembled: the cross blocks between the two surface gradients are
   * generally nonzero after the equilibrium condensation and must not be omitted.
   */
  template < int nDim,
             int nNodes,
             class QuadraturePoint,
             class MaterialState,
             class MaterialTangents,
             class PeType,
             class KeType >
  void assembleTwoSidedResponse( QuadraturePoint&        qp,
                                 const MaterialState&    state,
                                 const MaterialTangents& tangents,
                                 PeType&&                Pe,
                                 KeType&&                Ke )
  {
    using namespace Marmot;

    constexpr int nTensor        = nDim * nDim;
    constexpr int sizeLoadVector = nNodes * nDim;
    constexpr int halfSize       = sizeLoadVector / 2;

    using ForceSized         = Eigen::Matrix< double, nDim, 1 >;
    using SurfaceStressSized = Eigen::Matrix< double, nTensor, 1 >;
    using BFullSized         = Eigen::Matrix< double, nTensor, sizeLoadVector >;
    using QwwMatrixSized     = Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor >;
    using QwAMatrixSized     = Eigen::Matrix< double, nDim, nTensor, Eigen::RowMajor >;
    using QAwMatrixSized     = Eigen::Matrix< double, nTensor, nDim, Eigen::RowMajor >;
    using QAAMatrixSized     = Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor >;

    const auto& Njump = qp.NmatJump;
    const auto& Bside = qp.BmatSide;

    // the surface gradients of the top and of the bottom side
    BFullSized BPlusFull                                         = BFullSized::Zero();
    BFullSized BMinusFull                                        = BFullSized::Zero();
    BMinusFull.template block< nTensor, halfSize >( 0, 0 )       = Bside;
    BPlusFull.template block< nTensor, halfSize >( 0, halfSize ) = Bside;

    ForceSized         generalizedForce;
    SurfaceStressSized surfaceStressPlus;
    SurfaceStressSized surfaceStressMinus;

    QwwMatrixSized Q_ww;
    QwAMatrixSized Q_wAp, Q_wAm;
    QAwMatrixSized Q_Apw, Q_Amw;
    QAAMatrixSized Q_ApAp, Q_ApAm, Q_AmAp, Q_AmAm;

    if constexpr ( nDim == 3 ) {
      generalizedForce   = mapEigenToFastor( state.generalizedForce );
      surfaceStressPlus  = Eigen::Map< const SurfaceStressSized >( state.surfaceStressPlus.data() );
      surfaceStressMinus = Eigen::Map< const SurfaceStressSized >( state.surfaceStressMinus.data() );

      Q_ww   = mapEigenToFastor( tangents.Q_ww );
      Q_wAp  = mapEigenToFastor( tangents.Q_wAp );
      Q_wAm  = mapEigenToFastor( tangents.Q_wAm );
      Q_Apw  = Eigen::Map< const QAwMatrixSized >( tangents.Q_Apw.data() );
      Q_Amw  = Eigen::Map< const QAwMatrixSized >( tangents.Q_Amw.data() );
      Q_ApAp = mapEigenToFastor( tangents.Q_ApAp );
      Q_ApAm = mapEigenToFastor( tangents.Q_ApAm );
      Q_AmAp = mapEigenToFastor( tangents.Q_AmAp );
      Q_AmAm = mapEigenToFastor( tangents.Q_AmAm );
    }
    else {
      const Fastor::Tensor< double, 2 >    force2d              = reduceTo2D< U >( state.generalizedForce );
      const Fastor::Tensor< double, 2, 2 > surfaceStressPlus2d  = reduceTo2D< U, U >( state.surfaceStressPlus );
      const Fastor::Tensor< double, 2, 2 > surfaceStressMinus2d = reduceTo2D< U, U >( state.surfaceStressMinus );

      generalizedForce   = mapEigenToFastor( force2d );
      surfaceStressPlus  = Eigen::Map< const SurfaceStressSized >( surfaceStressPlus2d.data() );
      surfaceStressMinus = Eigen::Map< const SurfaceStressSized >( surfaceStressMinus2d.data() );

      const Fastor::Tensor< double, 2, 2 > Qww2d = reduceTo2D< U, U >( tangents.Q_ww );

      const Fastor::Tensor< double, 2, 2, 2 > QwAp2d = reduceTo2D< U, U, U >( tangents.Q_wAp );
      const Fastor::Tensor< double, 2, 2, 2 > QwAm2d = reduceTo2D< U, U, U >( tangents.Q_wAm );
      const Fastor::Tensor< double, 2, 2, 2 > QApw2d = reduceTo2D< U, U, U >( tangents.Q_Apw );
      const Fastor::Tensor< double, 2, 2, 2 > QAmw2d = reduceTo2D< U, U, U >( tangents.Q_Amw );

      const Fastor::Tensor< double, 2, 2, 2, 2 > QApAp2d = reduceTo2D< U, U, U, U >( tangents.Q_ApAp );
      const Fastor::Tensor< double, 2, 2, 2, 2 > QApAm2d = reduceTo2D< U, U, U, U >( tangents.Q_ApAm );
      const Fastor::Tensor< double, 2, 2, 2, 2 > QAmAp2d = reduceTo2D< U, U, U, U >( tangents.Q_AmAp );
      const Fastor::Tensor< double, 2, 2, 2, 2 > QAmAm2d = reduceTo2D< U, U, U, U >( tangents.Q_AmAm );

      Q_ww   = mapEigenToFastor( Qww2d );
      Q_wAp  = mapEigenToFastor( QwAp2d );
      Q_wAm  = mapEigenToFastor( QwAm2d );
      Q_Apw  = Eigen::Map< const QAwMatrixSized >( QApw2d.data() );
      Q_Amw  = Eigen::Map< const QAwMatrixSized >( QAmw2d.data() );
      Q_ApAp = mapEigenToFastor( QApAp2d );
      Q_ApAm = mapEigenToFastor( QApAm2d );
      Q_AmAp = mapEigenToFastor( QAmAp2d );
      Q_AmAm = mapEigenToFastor( QAmAm2d );
    }

    qp.managedStateVars->generalizedForce   = generalizedForce;
    qp.managedStateVars->surfaceStressPlus  = surfaceStressPlus;
    qp.managedStateVars->surfaceStressMinus = surfaceStressMinus;

    Pe += Njump.transpose() * generalizedForce * qp.J0xW;
    Pe += BPlusFull.transpose() * surfaceStressPlus * qp.J0xW;
    Pe += BMinusFull.transpose() * surfaceStressMinus * qp.J0xW;

    Ke += ( Njump.transpose() * Q_ww * Njump + Njump.transpose() * Q_wAp * BPlusFull +
            Njump.transpose() * Q_wAm * BMinusFull + BPlusFull.transpose() * Q_Apw * Njump +
            BPlusFull.transpose() * Q_ApAp * BPlusFull + BPlusFull.transpose() * Q_ApAm * BMinusFull +
            BMinusFull.transpose() * Q_Amw * Njump + BMinusFull.transpose() * Q_AmAp * BPlusFull +
            BMinusFull.transpose() * Q_AmAm * BMinusFull ) *
          qp.J0xW;
  }

  /**
   * @struct EquilibratedTwoSidedFormulation
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @brief Formulation of the interface element with equilibrated two-sided kinematics.
   *
   * The equilibrated interface material receives the increment of the displacement of both interface sides and of
   * the surface gradients of both sides, together with the normal of the interface and the connector between the
   * two faces at the quadrature point. It returns the generalized force, conjugate to the displacement jump, the
   * generalized surface stresses, conjugate to the two surface gradients, and the condensed block tangent. The
   * residual and the tangent are assembled with the jump operator and the operators of the surface gradients of the
   * top and of the bottom side of the quadrature point. The formulation uses linearized, small-deformation
   * kinematics.
   */
  template < int nDim, int nNodes >
  struct EquilibratedTwoSidedFormulation {

    static constexpr int nTensor        = nDim * nDim;
    static constexpr int sizeLoadVector = nNodes * nDim;
    static constexpr int halfSize       = sizeLoadVector / 2;

    using VectorDim = Eigen::Matrix< double, nDim, 1 >;
    using TensorDim = Eigen::Matrix< double, nDim, nDim >;

    using ForceSized                = Eigen::Matrix< double, nDim, 1 >;
    using SurfaceStressSized        = Eigen::Matrix< double, nTensor, 1 >;
    using InterfaceDisplSized       = Eigen::Matrix< double, 2 * nDim, 1 >;
    using InterfaceSurfaceGradSized = Eigen::Matrix< double, 2 * nTensor, 1 >;

    /// operator that maps the displacements of the element to the surface gradient of one side
    using BFullSized = Eigen::Matrix< double, nTensor, sizeLoadVector >;

    using QwwMatrixSized = Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor >;
    using QwAMatrixSized = Eigen::Matrix< double, nDim, nTensor, Eigen::RowMajor >;
    using QAwMatrixSized = Eigen::Matrix< double, nTensor, nDim, Eigen::RowMajor >;
    using QAAMatrixSized = Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor >;

    using Material = MarmotEquilibratedXInterfaceMaterialHypoElastic; ///< interface material of the element

    /**
     * @brief Named state-variable manager for interface quadrature points.
     */
    class QPStateVarManager : public MarmotStateVarVectorManager {

      /*
       * Persistent state layout for the accumulated generalized force, the two generalized surface stresses,
       * displacement, surface strain, and material state variables.
       *
       * The displacement and surface strain entries store the accumulated top/bottom quantities:
       *
       *   displacement   = [u_top, u_bottom]
       *   surface strain = [grad_s u_top, grad_s u_bottom]
       */
      inline const static auto layout = makeLayout( {
        { .name = "generalizedForce", .length = nDim },
        { .name = "surfaceStressPlus", .length = nDim * nDim },
        { .name = "surfaceStressMinus", .length = nDim * nDim },
        { .name = "displacement", .length = 2 * nDim },
        { .name = "surfaceStrain", .length = 2 * nDim * nDim },
        // The material state starts at a multiple of four entries.
        { .name   = "state block alignment padding",
          .length = ( 4 - ( ( nDim + nDim * nDim + nDim * nDim + 2 * nDim + 2 * nDim * nDim ) % 4 ) ) % 4 },
        { .name = "begin of material state", .length = 0 },
      } );

    public:
      Eigen::Map< ForceSized >                generalizedForce;
      Eigen::Map< SurfaceStressSized >        surfaceStressPlus;
      Eigen::Map< SurfaceStressSized >        surfaceStressMinus;
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
          generalizedForce( &find( "generalizedForce" ) ),
          surfaceStressPlus( &find( "surfaceStressPlus" ) ),
          surfaceStressMinus( &find( "surfaceStressMinus" ) ),
          displacement( &find( "displacement" ) ),
          surfaceStrain( &find( "surfaceStrain" ) ),
          materialStateVars( &find( "begin of material state" ),
                             nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
      {
      }
    };

    /**
     * @brief Data of the formulation at a quadrature point, besides the state variables.
     */
    struct QuadraturePointData : TwoSidedGeometryData< nDim > {
      /**
       * @brief Memory of the material that has to survive the iterations of an increment.
       *
       * The state variables are restored at the beginning of every iteration, so this memory is kept outside of
       * them, see MarmotEquilibratedXInterfaceMaterialHypoElastic::IncrementMemory.
       */
      Material::IncrementMemory gaugeMemory;
    };

    /**
     * @brief Complete the geometry of a quadrature point with the connector between the two faces.
     *
     * The connector is the difference of the positions of the top and of the bottom face at the quadrature point. The
     * normal is oriented along the connector, from the bottom to the top face, and the normal and the tangential
     * separation follow from it. Coincident faces have no connector.
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
      initializeTwoSidedQuadraturePoint( element, qp, "YInterfaceFiniteElement" );
    }

    /**
     * @brief Evaluate the material at a quadrature point and add its contributions to the element.
     *
     * Updates the state variables of the quadrature point and adds the contributions of the quadrature point to
     * the residual and to the tangent of the element. The residual is the internal force, so that the tangent is
     * its derivative.
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

      auto& material = element.material;
      auto& data     = qp.formulationData;

      Material::State state{ embedVectorIn3D< nDim >( qp.managedStateVars->generalizedForce.data() ),
                             embedTensorIn3D< nDim >( qp.managedStateVars->surfaceStressPlus.data() ),
                             embedTensorIn3D< nDim >( qp.managedStateVars->surfaceStressMinus.data() ),
                             qp.managedStateVars->materialStateVars.data(),
                             &data.gaugeMemory };

      Material::Tangents tangents{ Material::Tensor33d( 0.0 ),
                                   Material::Tensor333d( 0.0 ),
                                   Material::Tensor333d( 0.0 ),
                                   Material::Tensor333d( 0.0 ),
                                   Material::Tensor3333d( 0.0 ),
                                   Material::Tensor3333d( 0.0 ),
                                   Material::Tensor333d( 0.0 ),
                                   Material::Tensor3333d( 0.0 ),
                                   Material::Tensor3333d( 0.0 ) };

      // the materials are fixed to three dimensions, a two-dimensional element embeds its kinematics
      const Eigen::Matrix< double, 6, 1 >  dU3d             = expandSidesTo3D< nDim, 1 >( dU_GPs );
      const Eigen::Matrix< double, 18, 1 > dSurfaceStrain3d = expandSidesTo3D< nDim, 2 >( dSurface_strain_GPs );

      Material::Deformation   deformation{ Material::Tensor6d( dU3d.data() ),
                                         Material::Tensor18d( dSurfaceStrain3d.data() ),
                                         embedVectorIn3D< nDim >( qp.normal.data() ),
                                         embedVectorIn3D< nDim >( data.separationVector.data() ) };
      Material::TimeIncrement timeIncrement{ time, dT };

      material->computeStress( state, tangents, deformation, timeIncrement );

      qp.managedStateVars->displacement += dU_GPs;
      qp.managedStateVars->surfaceStrain += dSurface_strain_GPs;

      assembleTwoSidedResponse< nDim, nNodes >( qp, state, tangents, Pe, Ke );
    }
  };

  /**
   * @brief Equilibrated two-sided interface finite element.
   *
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   */
  template < int nDim, int nNodes >
  using YInterfaceFiniteElement = InterfaceFiniteElementBase< nDim,
                                                              nNodes,
                                                              EquilibratedTwoSidedFormulation< nDim, nNodes > >;

} // namespace Marmot::Elements
