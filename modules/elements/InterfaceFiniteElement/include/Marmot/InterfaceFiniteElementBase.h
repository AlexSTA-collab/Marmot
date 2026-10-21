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
 * @file InterfaceFiniteElementBase.h
 * @brief Common base of the interface finite elements.
 *
 * This file defines the templated class `Marmot::Elements::InterfaceFiniteElementBase` and its full
 * in-header implementation. It holds everything the interface elements have in common: the geometry
 * of the interface and of its quadrature points, the bookkeeping of the nodes, degrees of freedom and
 * state variables, the assignment of the material and the kinematic increments at the quadrature
 * points. What differs between the interface elements, namely the material, the state variables stored at
 * a quadrature point and the evaluation of the material response together with the assembly of the
 * residual and the tangent, is supplied by a formulation, see InterfaceFiniteElementBase.
 */
#pragma once

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFastorTensorBasics.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotStateVarVectorManager.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Marmot::Elements {

  /**
   * @brief Number of degrees of freedom of the element in addition to the displacements.
   *
   * A formulation may carry further fields, for instance a nonlocal field, on the nodes of the element. It declares
   * the number of these degrees of freedom of the element with a static member `nExtraDofs`, which is placed after the
   * displacements in the vectors of the element; a formulation without the member has none.
   *
   * @tparam Formulation the formulation of the element
   * @return the number of additional degrees of freedom
   */
  template < class Formulation >
  constexpr int numberOfExtraDofs()
  {
    if constexpr ( requires { Formulation::nExtraDofs; } )
      return Formulation::nExtraDofs;
    else
      return 0;
  }

  /**
   * @brief Embed a vector with one block of the interface quantity per side in the 3D layout of the materials.
   *
   * The interface materials are fixed to three dimensions. A two-dimensional element embeds its quantities of both
   * interface sides, which are stored side after side, in the 3D layout: each block is expanded with expandTo3D().
   *
   * @tparam nDim spatial dimension of the element, 2 or 3
   * @tparam rank 1 for a vector per side (jump), 2 for a second-order tensor per side (surface gradient)
   * @param[in] sides the blocks of the bottom and the top side, each holding nDim^rank entries
   * @return the expanded blocks, each holding 3^rank entries
   */
  template < int nDim, int rank >
  Eigen::Matrix< double, 2 * ( rank == 1 ? 3 : 9 ), 1 > expandSidesTo3D(
    const Eigen::Matrix< double, 2 * ( rank == 1 ? nDim : nDim * nDim ), 1 >& sides )
  {
    constexpr int nSide3d = rank == 1 ? 3 : 9;
    constexpr int nSide   = rank == 1 ? nDim : nDim * nDim;

    if constexpr ( nDim == 3 ) {
      return sides;
    }
    else {
      Eigen::Matrix< double, 2 * nSide3d, 1 > expanded = Eigen::Matrix< double, 2 * nSide3d, 1 >::Zero();

      for ( int side = 0; side < 2; ++side ) {
        if constexpr ( rank == 1 ) {
          const Fastor::Tensor< double, 3 > block = expandTo3D(
            Fastor::Tensor< double, nDim >( sides.data() + side * nSide ) );
          expanded.template segment< nSide3d >(
            side * nSide3d ) = Eigen::Map< const Eigen::Matrix< double, nSide3d, 1 > >( block.data() );
        }
        else {
          const Fastor::Tensor< double, 3, 3 > block = expandTo3D(
            Fastor::Tensor< double, nDim, nDim >( sides.data() + side * nSide ) );
          expanded.template segment< nSide3d >(
            side * nSide3d ) = Eigen::Map< const Eigen::Matrix< double, nSide3d, 1 > >( block.data() );
        }
      }
      return expanded;
    }
  }

  /**
   * @class InterfaceFiniteElementBase
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @tparam Formulation Policy that supplies what differs between the interface elements.
   * @brief Common base of the interface finite elements with an interface-material update at each quadrature
   * point.
   *
   * The element combines the geometric interface operators from `MarmotGeometryInterfaceElement<nDim, nNodes>`
   * with an interface material, which is constructed once per element and shared by all quadrature points. At each
   * quadrature point it stores state variables, and it assembles the element residual vector and the algorithmic
   * tangent matrix. Inertia, distributed loads and body forces are not implemented.
   *
   * The @p Formulation has to provide
   * - `Material`: the interface material class, which can be created from a material name, its properties, their
   *   number and the element label, and which provides `setCharacteristicElementLength()`,
   *   `getNumberOfRequiredStateVars()`, `initializeYourself()` and `getStateView()`,
   * - `QPStateVarManager`: the named state variables of a quadrature point. It can be created from a pointer to
   *   them and their number, provides the static function `getNumberOfRequiredStateVarsQuadraturePointOnly()` and
   *   the vector `materialStateVars` of the material state variables, and is a MarmotStateVarVectorManager,
   * - `QuadraturePointData`: the data a formulation keeps at each quadrature point besides the state variables,
   *   for instance additional geometry. It can be default constructed, and may be an empty class,
   * - `initializeQuadraturePoint()`: a static function template that is called for every quadrature point once the
   *   geometry of the point has been evaluated, with the element and the point. It completes the geometry and
   *   fills the QuadraturePointData,
   * - optionally `nExtraDofs`, a static member with the number of degrees of freedom in addition to the displacements,
   *   which are placed after them in the vectors of the element, together with the static functions `nodeFields()` and
   *   `dofIndicesPermutationPattern()` that describe their place on the nodes,
   * - optionally `computeElementKernels()`, a static function template that evaluates the whole element from the
   *   total nodal values and their increment, the time and the time increment, and adds the contributions to the
   *   residual and the tangent. A formulation with an unknown that couples the quadrature points, for instance an
   *   internal degree of freedom that is condensed at the level of the element, has to provide it, and then needs
   *   no `integrate()`,
   * - `integrate()`: a static function template that evaluates the material at a quadrature point from the
   *   total nodal displacements of the element and the increments of the interface displacement and of the surface
   *   gradient, updates the state variables of the
   *   quadrature point, and adds the contributions to the residual and the tangent of the element, see
   *   DisplacementJumpFormulation and EquilibratedTwoSidedFormulation.
   */
  template < int nDim, int nNodes, class Formulation >
  class InterfaceFiniteElementBase : public MarmotElement, public MarmotGeometryInterfaceElement< nDim, nNodes > {

  public:
    /**
     * @brief Section type selector for the interface element.
     */
    enum SectionType {
      Interface,
    };

    static constexpr int nDofPerNodeU    = nDim;
    static constexpr int nInterfaceNodes = nNodes / 2;

    /// number of degrees of freedom of the displacements, the bottom side first
    static constexpr int nDofDisplacement = nNodes * nDim;

    /// number of degrees of freedom of the formulation in addition to the displacements, placed after them
    static constexpr int nExtraDofs = numberOfExtraDofs< Formulation >();

    static constexpr int sizeLoadVector = nDofDisplacement + nExtraDofs;
    static constexpr int nCoordinates   = nNodes * nDim;

    static constexpr int nSideDofs = nInterfaceNodes * nDim;
    static constexpr int nTensor   = nDim * nDim;

    using ParentGeometryElement = MarmotGeometryInterfaceElement< nDim, nNodes >;

    using XiSized              = typename ParentGeometryElement::XiSized;
    using NSized               = typename ParentGeometryElement::NSized;
    using dNdXiSized           = typename ParentGeometryElement::dNdXiSized;
    using SurfaceJacobianSized = typename ParentGeometryElement::SurfaceJacobianSized;
    using MetricSized          = typename ParentGeometryElement::MetricSized;
    using GradSized            = typename ParentGeometryElement::GradSized;

    using VectorDim = typename ParentGeometryElement::VectorDim;
    using TensorDim = typename ParentGeometryElement::TensorDim;

    using NMatrixSized     = typename ParentGeometryElement::NMatrixSized;
    using NJumpMatrixSized = typename ParentGeometryElement::NJumpMatrixSized;

    using BSurfaceSized    = typename ParentGeometryElement::BSurfaceSized;
    using BAvgSurfaceSized = typename ParentGeometryElement::BAvgSurfaceSized;

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;

    using ForceSized                = Eigen::Matrix< double, nDim, 1 >;
    using SurfaceStressSized        = Eigen::Matrix< double, nTensor, 1 >;
    using InterfaceDisplSized       = Eigen::Matrix< double, 2 * nDim, 1 >;
    using InterfaceSurfaceGradSized = Eigen::Matrix< double, 2 * nTensor, 1 >;

    using Material          = typename Formulation::Material;          ///< interface material of the element
    using QPStateVarManager = typename Formulation::QPStateVarManager; ///< named state variables of a quadrature point

    Eigen::Map< const Eigen::VectorXd > elementProperties;

    const int         elLabel;
    const SectionType sectionType;

    /**
     * @brief Container for all per-quadrature-point data.
     */
    struct QuadraturePoint {
      const XiSized xi;
      const double  weight;

      double detJ;
      double sqrtDetG;
      double J0xW;

      NSized               N;
      dNdXiSized           dNdXi;
      SurfaceJacobianSized J;
      MetricSized          G;
      GradSized            gradN;

      VectorDim normal;
      TensorDim normalProjection;
      TensorDim tangentProjection;

      /*
       * One-side operators:
       *
       * NmatSide:
       *   maps one side's nodal displacement vector to u at the interface qp.
       *
       * BmatSide:
       *   maps one side's nodal displacement vector to the surface-gradient
       *   quantity passed to the material.
       */
      NMatrixSized  NmatSide;
      BSurfaceSized BmatSide;

      /*
       * Whole-element operators:
       *
       * NmatJump:
       *   jump u = u_top - u_bottom
       *   NmatJump = [ -Nside , +Nside ]
       *
       * BmatAverage:
       *   grad_s u_avg = 0.5 * (grad_s u_bottom + grad_s u_top)
       *   BmatAverage = 0.5 * [ Bside , Bside ]
       */
      NJumpMatrixSized NmatJump;
      BAvgSurfaceSized BmatAverage;

      std::unique_ptr< QPStateVarManager > managedStateVars;

      /**
       * @brief Characteristic length of this quadrature point, forwarded to the material before each
       * stress evaluation. It varies between quadrature points on a distorted interface, so it cannot
       * live in the shared material object.
       */
      double characteristicLength = 0.0;

      /**
       * @brief Data of the formulation at this quadrature point, besides the state variables.
       */
      typename Formulation::QuadraturePointData formulationData;

      /**
       * @brief Interpolate a field of nodal values to the quadrature point on both interface sides.
       *
       * @tparam Vector type of the nodal vector
       * @param[in] q nodal values of the element, for instance the displacements or their increments, the bottom
       * side first
       * @return the values of the top side followed by the values of the bottom side
       */
      template < class Vector >
      InterfaceDisplSized interpolateToSides( const Vector& q ) const
      {
        constexpr int halfSize = nNodes * nDim / 2;

        InterfaceDisplSized atSides;
        atSides.template segment< nDim >( 0 )    = NmatSide * q.template segment< halfSize >( halfSize );
        atSides.template segment< nDim >( nDim ) = NmatSide * q.template segment< halfSize >( 0 );
        return atSides;
      }

      /**
       * @brief Interpolate the surface gradient of a field of nodal values to the quadrature point on both interface
       * sides.
       *
       * @tparam Vector type of the nodal vector
       * @param[in] q nodal values of the element, for instance the displacements or their increments, the bottom
       * side first
       * @return the surface gradients of the top side followed by the ones of the bottom side
       */
      template < class Vector >
      InterfaceSurfaceGradSized surfaceGradientsAtSides( const Vector& q ) const
      {
        constexpr int halfSize = nNodes * nDim / 2;

        InterfaceSurfaceGradSized atSides;
        atSides.template segment< nTensor >( 0 )       = BmatSide * q.template segment< halfSize >( halfSize );
        atSides.template segment< nTensor >( nTensor ) = BmatSide * q.template segment< halfSize >( 0 );
        return atSides;
      }

      /**
       * @brief Number of non-material state variables stored at this quadrature point.
       */
      int getNumberOfRequiredStateVarsQuadraturePointOnly()
      {
        return QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly();
      }

      /**
       * @brief Total number of state variables including material state.
       *
       * @param nMaterialStateVars Number of state variables the element's material requires. The
       * material is owned by the element, not by the quadrature point, so its size is passed in.
       */
      int getNumberOfRequiredStateVars( int nMaterialStateVars )
      {
        return getNumberOfRequiredStateVarsQuadraturePointOnly() + nMaterialStateVars;
      }

      /**
       * @brief Assign state-variable memory to this quadrature point.
       */
      void assignStateVars( double* stateVars, int nStateVars )
      {
        managedStateVars = std::make_unique< QPStateVarManager >( stateVars, nStateVars );
      }

      /**
       * @brief Construct a quadrature-point data container.
       */
      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ),
          weight( weight ),
          detJ( 0.0 ),
          sqrtDetG( 0.0 ),
          J0xW( 0.0 ),
          N( NSized::Zero() ),
          dNdXi( dNdXiSized::Zero() ),
          J( SurfaceJacobianSized::Zero() ),
          G( MetricSized::Zero() ),
          gradN( GradSized::Zero() ),
          normal( VectorDim::Zero() ),
          normalProjection( TensorDim::Zero() ),
          tangentProjection( TensorDim::Zero() ),
          NmatSide( NMatrixSized::Zero() ),
          BmatSide( BSurfaceSized::Zero() ),
          NmatJump( NJumpMatrixSized::Zero() ),
          BmatAverage( BAvgSurfaceSized::Zero() )
      {
      }
    };

    std::vector< QuadraturePoint > qps;

    /**
     * @brief The interface material, constructed once per element.
     *
     * All quadrature points of an element share one property set, and the material carries no
     * per-quadrature-point state: the state variables are passed in by pointer and the characteristic
     * length is set from QuadraturePoint::characteristicLength before every evaluation. Constructing
     * one instance per quadrature point would repeat the base material's setup work, which is
     * expensive for LinearViscoElasticWiechert.
     */
    std::unique_ptr< Material > material;

    /**
     * @brief Construct an interface finite element.
     * @param[in] elementID Element label.
     * @param[in] integrationType Quadrature rule type.
     * @param[in] sectionType Interface section type.
     */
    InterfaceFiniteElementBase( int                                         elementID,
                                FiniteElement::Quadrature::IntegrationTypes integrationType,
                                SectionType                                 sectionType = SectionType::Interface );

    /**
     * @brief Return number of required state variables per element.
     */
    int getNumberOfRequiredStateVars();

    /**
     * @brief Return nodal primary fields associated with this element.
     */
    std::vector< std::vector< std::string > > getNodeFields();

    /**
     * @brief Return DOF permutation pattern.
     */
    std::vector< int > getDofIndicesPermutationPattern();

    /** @brief Number of nodes of the element, i.e. of both interface sides. */
    int getNNodes() { return nNodes; }

    /** @brief Number of spatial dimensions of the element. */
    int getNSpatialDimensions() { return nDim; }

    /** @brief Number of degrees of freedom of the element. */
    int getNDofPerElement() { return sizeLoadVector; }

    /** @brief Result-file geometry keyword, delegated to the geometry element. */
    std::string getElementShape() { return ParentGeometryElement::getElementShape(); }

    /**
     * @brief Assign element state-variable memory to quadrature points.
     */
    void assignStateVars( double* stateVars, int nStateVars );

    /**
     * @brief Assign element-level properties (e.g., extrusion thickness of the element).
     */
    void assignProperty( const ElementProperties& marmotElementProperty );

    /**
     * @brief Assign a material section to all quadrature points.
     */
    void assignProperty( const MarmotMaterialSection& marmotElementProperty );

    /**
     * @brief Assign a material by name and property array to all quadrature points.
     */
    void assignMaterial( const std::string& materialName, const double* materialProperties, int nMaterialProperties );

    /**
     * @brief Assign nodal coordinates.
     */
    void assignNodeCoordinates( const double* coordinates );

    /**
     * @brief Initialize element geometric operators and quadrature-point geometry data.
     */
    void initializeYourself();

    /**
     * @brief Set initial conditions for supported state categories.
     */
    void setInitialConditions( StateTypes state, const double* values );

    /**
     * @brief Distributed-load routine (currently not implemented).
     */
    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             P,
                                 double*                             K,
                                 const int                           elementFace,
                                 const double*                       load,
                                 const double*                       QTotal,
                                 double                              time,
                                 double                              dT );

    /**
     * @brief Body-force routine (currently not implemented).
     */
    void computeBodyForce( double* P, double* K, const double* load, const double* QTotal, double time, double dT );

    /**
     * @brief Assemble internal force and tangent for one increment.
     *
     * @param QTotal Total displacement vector.
     * @param dQ Incremental displacement.
     * @param Pe Internal force vector (accumulated), same convention as every other element:
     *           Pe = +Pint, so that Ke = +dPe/dQ. External loads are accumulated separately
     *           by computeDistributedLoad() and computeBodyForce() into Pext.
     * @param Ke Tangent stiffness matrix (accumulated).
     * @param time Time data forwarded to materials.
     * @param dT Time increment.
     */
    void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT );

    /**
     * @brief Report that consistent inertia is unsupported for interface elements.
     */
    void computeConsistentInertia( double* M );

    /**
     * @brief Report that lumped inertia is unsupported for interface elements.
     */
    void computeLumpedInertia( double* M );

    /**
     * @brief Access a named state view at a quadrature point.
     */
    StateView getStateView( const std::string& stateName, int qpNumber )
    {
      const auto& qp = qps[qpNumber];

      if ( qp.managedStateVars->contains( stateName ) ) {
        return qp.managedStateVars->getStateView( stateName );
      }

      if ( stateName == "sdv" ) {
        std::cout << __PRETTY_FUNCTION__ << " on 'sdv' is discouraged and deprecated, please use precise state name";
        return { qp.managedStateVars->materialStateVars.data(),
                 static_cast< int >( qp.managedStateVars->materialStateVars.size() ) };
      }

      return material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
    }

    /**
     * @brief Coordinates of the element center (on reference interface side).
     */
    std::vector< double > getCoordinatesAtCenter();

    /**
     * @brief Coordinates of all quadrature points (on reference interface side).
     */
    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    /**
     * @brief Number of quadrature points.
     */
    int getNumberOfQuadraturePoints();
  };

  /**
   * @name Template implementation (header-only)
   * @brief In-header method definitions for `InterfaceFiniteElementBase`.
   */
  ///@{

  template < int nDim, int nNodes, class Formulation >
  InterfaceFiniteElementBase< nDim, nNodes, Formulation >::InterfaceFiniteElementBase(
    int                                         elementID,
    FiniteElement::Quadrature::IntegrationTypes integrationType,
    SectionType                                 sectionType )
    : ParentGeometryElement(),
      elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType )
  {
    const auto qpInfos = FiniteElement::Quadrature::getGaussPointInfo( this->shape, integrationType );

    for ( const auto& qpInfo : qpInfos ) {
      QuadraturePoint qp( qpInfo.xi, qpInfo.weight );
      qps.push_back( std::move( qp ) );
    }
  }

  template < int nDim, int nNodes, class Formulation >
  int InterfaceFiniteElementBase< nDim, nNodes, Formulation >::getNumberOfRequiredStateVars()
  {
    const int nMaterialStateVars = material ? material->getNumberOfRequiredStateVars() : 0;

    return qps[0].getNumberOfRequiredStateVars( nMaterialStateVars ) * qps.size();
  }

  template < int nDim, int nNodes, class Formulation >
  std::vector< std::vector< std::string > > InterfaceFiniteElementBase< nDim, nNodes, Formulation >::getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;

    if ( nodeFields.empty() ) {
      if constexpr ( requires { Formulation::nodeFields(); } ) {
        nodeFields = Formulation::nodeFields();
      }
      else {
        for ( int i = 0; i < nNodes; i++ ) {
          nodeFields.push_back( vector< string >() );
          nodeFields[i].push_back( "displacement" );
        }
      }
    }

    return nodeFields;
  }

  template < int nDim, int nNodes, class Formulation >
  std::vector< int > InterfaceFiniteElementBase< nDim, nNodes, Formulation >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() ) {
      if constexpr ( requires { Formulation::dofIndicesPermutationPattern(); } ) {
        permutationPattern = Formulation::dofIndicesPermutationPattern();
      }
      else {
        for ( int i = 0; i < nNodes * nDim; i++ )
          permutationPattern.push_back( i );
      }
    }

    return permutationPattern;
  }

  template < int nDim, int nNodes, class Formulation >
  void InterfaceFiniteElementBase< nDim, nNodes, Formulation >::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes, class Formulation >
  void InterfaceFiniteElementBase< nDim, nNodes, Formulation >::assignProperty(
    const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes, class Formulation >
  void InterfaceFiniteElementBase< nDim, nNodes, Formulation >::assignProperty( const MarmotMaterialSection& section )
  {
    material = std::make_unique< Material >( section.materialName,
                                             section.materialProperties,
                                             section.nMaterialProperties,
                                             elLabel );
  }

  template < int nDim, int nNodes, class Formulation >
  void InterfaceFiniteElementBase< nDim, nNodes, Formulation >::assignMaterial( const std::string& materialName,
                                                                                const double*      materialProperties,
                                                                                int                nMaterialProperties )
  {
    material = std::make_unique< Material >( materialName, materialProperties, nMaterialProperties, elLabel );
  }

  template < int nDim, int nNodes, class Formulation >
  void InterfaceFiniteElementBase< nDim, nNodes, Formulation >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes, class Formulation >
  void InterfaceFiniteElementBase< nDim, nNodes, Formulation >::initializeYourself()
  {
    const double thickness = elementProperties.size() > 0 ? elementProperties[0] : 1.0;

    for ( QuadraturePoint& qp : qps ) {
      const auto geom = this->evaluateAt( qp.xi, 0 );

      qp.N                 = geom.N;
      qp.dNdXi             = geom.dNdXi;
      qp.J                 = geom.J;
      qp.G                 = geom.G;
      qp.sqrtDetG          = geom.sqrtDetG;
      qp.detJ              = geom.sqrtDetG;
      qp.gradN             = geom.gradN;
      qp.normal            = geom.n;
      qp.normalProjection  = geom.normalProjection;
      qp.tangentProjection = geom.tangentProjection;

      qp.NmatSide    = geom.NmatSide;
      qp.BmatSide    = geom.BmatSide;
      qp.NmatJump    = geom.NmatJump;
      qp.BmatAverage = geom.BmatAverage;

      qp.J0xW = qp.weight * qp.sqrtDetG * thickness;

      if constexpr ( nDim == 3 ) {
        qp.characteristicLength = std::sqrt( qp.sqrtDetG );
      }
      else if constexpr ( nDim == 2 ) {
        qp.characteristicLength = qp.sqrtDetG;
      }

      Formulation::initializeQuadraturePoint( *this, qp );
    }
  }

  template < int nDim, int nNodes, class Formulation >
  void InterfaceFiniteElementBase< nDim, nNodes, Formulation >::computeKernels( const double* QTotal_,
                                                                                const double* dQ_,
                                                                                double*       Pe_,
                                                                                double*       Ke_,
                                                                                double        time,
                                                                                double        dT )
  {
    Eigen::Map< const RhsSized > QTotal( QTotal_ );
    Eigen::Map< const RhsSized > dQ( dQ_ );
    Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
    Eigen::Map< RhsSized >       Pe( Pe_ );

    if constexpr ( requires { Formulation::computeElementKernels( *this, QTotal, dQ, time, dT, Pe, Ke ); } ) {
      Formulation::computeElementKernels( *this, QTotal, dQ, time, dT, Pe, Ke );
    }
    else {
      for ( QuadraturePoint& qp : qps ) {
        const InterfaceDisplSized       dU_GPs              = qp.interpolateToSides( dQ );
        const InterfaceSurfaceGradSized dSurface_strain_GPs = qp.surfaceGradientsAtSides( dQ );

        // the shared material carries no per-quadrature-point state, so the characteristic length of
        // THIS point has to be installed before every evaluation
        material->setCharacteristicElementLength( qp.characteristicLength );

        Formulation::integrate( *this, qp, QTotal, dU_GPs, dSurface_strain_GPs, time, dT, Pe, Ke );
      }
    }
  }

  template < int nDim, int nNodes, class Formulation >
  void InterfaceFiniteElementBase< nDim, nNodes, Formulation >::setInitialConditions( StateTypes    state,
                                                                                      const double* values )
  {
    switch ( state ) {
    case MarmotElement::MarmotMaterialInitialization: {
      for ( QuadraturePoint& qp : qps ) {
        material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                      qp.managedStateVars->materialStateVars.size() );
      }
      break;
    }

    case MarmotElement::MarmotMaterialStateVars: {
      throw std::invalid_argument( "Please use initializeStateVars directly on material" );
    }

    default:
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": invalid initial condition for an interface finite element" );
    }
  }

  template < int nDim, int nNodes, class Formulation >
  void InterfaceFiniteElementBase< nDim, nNodes, Formulation >::computeDistributedLoad(
    MarmotElement::DistributedLoadTypes loadType,
    double*                             P,
    double*                             K,
    const int                           elementFace,
    const double*                       load,
    const double*                       QTotal,
    double                              time,
    double                              dT )
  {
    throw std::invalid_argument(
      MakeString() << __PRETTY_FUNCTION__ << ": distributed loads are not implemented for interface finite elements." );
  }

  template < int nDim, int nNodes, class Formulation >
  void InterfaceFiniteElementBase< nDim, nNodes, Formulation >::computeBodyForce( double*       P,
                                                                                  double*       K,
                                                                                  const double* load,
                                                                                  const double* QTotal,
                                                                                  double        time,
                                                                                  double        dT )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                              << ": body forces are not implemented for interface finite elements." );
  }

  template < int nDim, int nNodes, class Formulation >
  void InterfaceFiniteElementBase< nDim, nNodes, Formulation >::computeConsistentInertia( double* M )
  {
    throw std::runtime_error( MakeString()
                              << __PRETTY_FUNCTION__ << ": inertia is not implemented for interface finite elements." );
  }

  template < int nDim, int nNodes, class Formulation >
  void InterfaceFiniteElementBase< nDim, nNodes, Formulation >::computeLumpedInertia( double* M )
  {
    throw std::runtime_error( MakeString()
                              << __PRETTY_FUNCTION__ << ": inertia is not implemented for interface finite elements." );
  }

  template < int nDim, int nNodes, class Formulation >
  std::vector< double > InterfaceFiniteElementBase< nDim, nNodes, Formulation >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    Eigen::Map< VectorDim > coordsMap( coords.data() );

    const auto centerXi = XiSized::Zero();
    const auto Ncenter  = this->N( centerXi );
    const auto Nmat     = this->NMatrix( Ncenter );

    const auto xSide = this->getSideCoordinates( 0 );

    coordsMap = Nmat * xSide;

    return coords;
  }

  template < int nDim, int nNodes, class Formulation >
  std::vector< std::vector< double > > InterfaceFiniteElementBase< nDim, nNodes, Formulation >::
    getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;

    for ( const auto& qp : qps ) {
      std::vector< double > coords( nDim );

      Eigen::Map< VectorDim > coordsMap( coords.data() );

      const auto xSide = this->getSideCoordinates( 0 );
      coordsMap        = qp.NmatSide * xSide;

      listedCoords.push_back( coords );
    }

    return listedCoords;
  }

  template < int nDim, int nNodes, class Formulation >
  int InterfaceFiniteElementBase< nDim, nNodes, Formulation >::getNumberOfQuadraturePoints()
  {
    return static_cast< int >( qps.size() );
  }

  ///@}

} // namespace Marmot::Elements
