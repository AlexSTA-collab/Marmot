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
 * @file XInterfaceFiniteElement.h
 * @brief Through-thickness-resolved interface finite element (XIQUAD4/XILINE2).
 *
 * This file defines the templated class `Marmot::Elements::XInterfaceFiniteElement`
 * and its full in-header implementation. Unlike `CorrectedInterfaceFiniteElement`,
 * which contracts the material response against the AVERAGE of the top and
 * bottom surface-gradient increments, this element keeps the two sides'
 * surface-gradient B-matrices separate and routes them to
 * `MarmotXInterfaceMaterialHypoElastic`, which evaluates (and does not
 * average away) their individually contracted energies.
 */
#pragma once

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotXInterfaceMaterialHypoElastic.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Marmot::Elements {

  /**
   * @class XInterfaceFiniteElement
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @brief Interface finite element with separately-resolved top/bottom surface kinematics.
   */
  template < int nDim, int nNodes >
  class XInterfaceFiniteElement : public MarmotElement, public MarmotGeometryInterfaceElement< nDim, nNodes > {

  public:
    /**
     * @brief Section type selector for the interface element.
     */
    enum SectionType {
      Interface,
    };

    static constexpr int nDofPerNodeU    = nDim;
    static constexpr int nInterfaceNodes = nNodes / 2;

    static constexpr int sizeLoadVector = nNodes * nDim;
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

    using QMatrixSized = Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor >;
    using ZMatrixSized = Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor >;
    using HMatrixSized = Eigen::Matrix< double, nDim, nTensor, Eigen::RowMajor >;
    using KMatrixSized = Eigen::Matrix< double, nTensor, nDim, Eigen::RowMajor >;

    using Material = MarmotXInterfaceMaterialHypoElastic;

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

      // Actual reference connector between paired lower and upper points:
      //   d = x_top - x_bottom = ell * n + d_tangent.
      // For a zero-thickness interface mesh, d is zero and the material
      // falls back to its constitutive thickness h.
      VectorDim separationVector;
      VectorDim tangentialSeparation;
      double    normalSeparation;

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
       *
       * This element does not use BmatAverage for its own residual/stiffness
       * assembly (it keeps the two sides separate), but retains it for
       * parity with the geometry base class and potential diagnostics.
       */
      NJumpMatrixSized NmatJump;
      BAvgSurfaceSized BmatAverage;

      /**
       * @brief Named state-variable manager for interface quadrature points.
       */
      class QPStateVarManager : public MarmotStateVarVectorManager {

        /*
         * Persistent state layout for accumulated forces, surface stresses,
         * displacement, surface strain, and material state variables.
         *
         * The displacement and surface strain entries store the accumulated
         * top/bottom quantities:
         *
         *   displacement   = [u_top, u_bottom]
         *   surface strain = [grad_s u_top, grad_s u_bottom]
         */
        inline const static auto layout = makeLayout( {
          { .name = "forcePlus", .length = nDim },
          { .name = "forceMinus", .length = nDim },
          { .name = "surfaceStressPlus", .length = nDim * nDim },
          { .name = "surfaceStressMinus", .length = nDim * nDim },
          { .name = "displacement", .length = 2 * nDim },
          { .name = "surfaceStrain", .length = 2 * nDim * nDim },
          { .name   = "state block alignment padding",
            .length = ( 4 - ( ( 2 * nDim + 2 * nDim * nDim + 2 * nDim + 2 * nDim * nDim ) % 4 ) ) % 4 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< ForceSized >                forcePlus;
        Eigen::Map< ForceSized >                forceMinus;
        Eigen::Map< SurfaceStressSized >        surfaceStressPlus;
        Eigen::Map< SurfaceStressSized >        surfaceStressMinus;
        Eigen::Map< InterfaceDisplSized >       displacement;
        Eigen::Map< InterfaceSurfaceGradSized > surfaceStrain;

        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            forcePlus( &find( "forcePlus" ) ),
            forceMinus( &find( "forceMinus" ) ),
            surfaceStressPlus( &find( "surfaceStressPlus" ) ),
            surfaceStressMinus( &find( "surfaceStressMinus" ) ),
            displacement( &find( "displacement" ) ),
            surfaceStrain( &find( "surfaceStrain" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
        {
        }
      };

      std::unique_ptr< QPStateVarManager > managedStateVars;
      std::unique_ptr< Material >          material;

      /**
       * @brief Number of non-material state variables stored at this quadrature point.
       */
      int getNumberOfRequiredStateVarsQuadraturePointOnly()
      {
        return QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly();
      }

      /**
       * @brief Total number of state variables including material state.
       */
      int getNumberOfRequiredStateVars()
      {
        return getNumberOfRequiredStateVarsQuadraturePointOnly() + material->getNumberOfRequiredStateVars();
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
          separationVector( VectorDim::Zero() ),
          tangentialSeparation( VectorDim::Zero() ),
          normalSeparation( 0.0 ),
          NmatSide( NMatrixSized::Zero() ),
          BmatSide( BSurfaceSized::Zero() ),
          NmatJump( NJumpMatrixSized::Zero() ),
          BmatAverage( BAvgSurfaceSized::Zero() )
      {
      }
    };

    std::vector< QuadraturePoint > qps;

    /**
     * @brief Construct an interface finite element.
     * @param[in] elementID Element label.
     * @param[in] integrationType Quadrature rule type.
     * @param[in] sectionType Interface section type.
     */
    XInterfaceFiniteElement( int                                         elementID,
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

    int getNNodes() { return nNodes; }

    int getNSpatialDimensions() { return nDim; }

    int getNDofPerElement() { return sizeLoadVector; }

    /*
     * Important:
     * ParentGeometryElement::getElementShape() returns the computational
     * interface shape, e.g. "iquad4".
     *
     * EnSight/ParaView do not understand "iquad4" as a geometry keyword.
     * Therefore, for output/visualization we map interface elements to valid
     * EnSight element names.
     *
     * The computational element remains an interface element. This only affects
     * the geometry keyword written to result files.
     */
    std::string getElementShape()
    {
      if constexpr ( nDim == 3 && nNodes == 8 ) {
        return "hexa8";
      }
      else if constexpr ( nDim == 2 && nNodes == 4 ) {
        return "bar2";
      }
      else {
        return ParentGeometryElement::getElementShape();
      }
    }

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
     * @brief Assemble residual and tangent for one increment.
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

      return qp.material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
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
   * @brief In-header method definitions for `XInterfaceFiniteElement`.
   */
  ///@{

  template < int nDim, int nNodes >
  XInterfaceFiniteElement< nDim, nNodes >::XInterfaceFiniteElement(
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

  template < int nDim, int nNodes >
  int XInterfaceFiniteElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > XInterfaceFiniteElement< nDim, nNodes >::getNodeFields()
  {
    using namespace std;

    static vector< vector< string > > nodeFields;

    if ( nodeFields.empty() ) {
      for ( int i = 0; i < nNodes; i++ ) {
        nodeFields.push_back( vector< string >() );
        nodeFields[i].push_back( "displacement" );
      }
    }

    return nodeFields;
  }

  template < int nDim, int nNodes >
  std::vector< int > XInterfaceFiniteElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() ) {
      for ( int i = 0; i < nNodes * nDim; i++ )
        permutationPattern.push_back( i );
    }

    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void XInterfaceFiniteElement< nDim, nNodes >::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes >
  void XInterfaceFiniteElement< nDim, nNodes >::assignProperty( const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void XInterfaceFiniteElement< nDim, nNodes >::assignProperty( const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::make_unique< MarmotXInterfaceMaterialHypoElastic >( section.materialName,
                                                                             section.materialProperties,
                                                                             section.nMaterialProperties,
                                                                             elLabel );
    }
  }

  template < int nDim, int nNodes >
  void XInterfaceFiniteElement< nDim, nNodes >::assignMaterial( const std::string& materialName,
                                                                const double*      materialProperties,
                                                                int                nMaterialProperties )
  {
    for ( auto& qp : qps ) {
      qp.material = std::make_unique< MarmotXInterfaceMaterialHypoElastic >( materialName,
                                                                             materialProperties,
                                                                             nMaterialProperties,
                                                                             elLabel );
    }
  }

  template < int nDim, int nNodes >
  void XInterfaceFiniteElement< nDim, nNodes >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void XInterfaceFiniteElement< nDim, nNodes >::initializeYourself()
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

      // Geometry of the actual top--bottom pairing used by NmatJump.
      // The interface normal is oriented from the lower side to the upper side.
      const VectorDim xBottom = qp.NmatSide * this->getSideCoordinates( 0 );
      const VectorDim xTop    = qp.NmatSide * this->getSideCoordinates( 1 );
      qp.separationVector     = xTop - xBottom;

      constexpr double geometryTolerance = 1.0e-12;
      if ( qp.separationVector.norm() > geometryTolerance && qp.separationVector.dot( qp.normal ) < 0.0 ) {
        qp.normal *= -1.0;
      }

      qp.normalProjection     = qp.normal * qp.normal.transpose();
      qp.tangentProjection    = TensorDim::Identity() - qp.normalProjection;
      qp.normalSeparation     = qp.separationVector.dot( qp.normal );
      qp.tangentialSeparation = qp.tangentProjection * qp.separationVector;

      if ( qp.separationVector.norm() > geometryTolerance && qp.normalSeparation <= geometryTolerance ) {
        throw std::invalid_argument(
          "XInterfaceFiniteElement: paired faces have no positive separation in the interface-normal direction." );
      }

      qp.J0xW = qp.weight * qp.sqrtDetG * thickness;

      if ( qp.material ) {
        if constexpr ( nDim == 3 ) {
          qp.material->setCharacteristicElementLength( std::sqrt( qp.sqrtDetG ) );
        }
        else if constexpr ( nDim == 2 ) {
          qp.material->setCharacteristicElementLength( qp.sqrtDetG );
        }
      }
    }
  }

  template < int nDim, int nNodes >
  void XInterfaceFiniteElement< nDim, nNodes >::computeKernels( const double* QTotal_,
                                                                const double* dQ_,
                                                                double*       Pe_,
                                                                double*       Ke_,
                                                                double        time,
                                                                double        dT )
  {
    (void)QTotal_;

    Eigen::Map< const RhsSized > dQ( dQ_ );
    Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
    Eigen::Map< RhsSized >       Pe( Pe_ );

    constexpr int halfSize = nNodes * nDim / 2;

    for ( QuadraturePoint& qp : qps ) {
      const auto& Nside = qp.NmatSide;
      const auto& Bside = qp.BmatSide;
      const auto& Njump = qp.NmatJump;

      // Full-DOF operators giving A+ (top) and A- (bottom) individually,
      // instead of collapsing them into BmatAverage. BmatAverage ==
      // 0.5*(BPlusFull + BMinusFull) by construction, matching
      // MarmotGeometryInterfaceElement::BAverageSurfaceMatrix's own
      // bottom-then-top column convention.
      BAvgSurfaceSized BPlusFull                                   = BAvgSurfaceSized::Zero();
      BAvgSurfaceSized BMinusFull                                  = BAvgSurfaceSized::Zero();
      BMinusFull.template block< nTensor, halfSize >( 0, 0 )       = Bside;
      BPlusFull.template block< nTensor, halfSize >( 0, halfSize ) = Bside;

      const auto dQBottom = dQ.template segment< halfSize >( 0 );
      const auto dQTop    = dQ.template segment< halfSize >( halfSize );

      InterfaceDisplSized dU_GPs;
      dU_GPs.template segment< nDim >( 0 )    = Nside * dQTop;
      dU_GPs.template segment< nDim >( nDim ) = Nside * dQBottom;

      InterfaceSurfaceGradSized dSurface_strain_GPs;
      dSurface_strain_GPs.template segment< nTensor >( 0 )       = Bside * dQTop;
      dSurface_strain_GPs.template segment< nTensor >( nTensor ) = Bside * dQBottom;

      ForceSized         forcePlus          = qp.managedStateVars->forcePlus;
      ForceSized         forceMinus         = qp.managedStateVars->forceMinus;
      SurfaceStressSized surfaceStressPlus  = qp.managedStateVars->surfaceStressPlus;
      SurfaceStressSized surfaceStressMinus = qp.managedStateVars->surfaceStressMinus;

      QMatrixSized Qplus, Qminus;
      ZMatrixSized Zplus, Zminus;
      HMatrixSized Hplus, Hminus;
      KMatrixSized Kplus, Kminus;

      Qplus.setZero();
      Qminus.setZero();
      Zplus.setZero();
      Zminus.setZero();
      Hplus.setZero();
      Hminus.setZero();
      Kplus.setZero();
      Kminus.setZero();

      if constexpr ( nDim == 3 ) {
        Material::State         materialState{ forcePlus.data(),
                                       forceMinus.data(),
                                       surfaceStressPlus.data(),
                                       surfaceStressMinus.data(),
                                       qp.managedStateVars->materialStateVars.data() };
        Material::Tangents      materialTangents{ Qplus.data(),
                                             Qminus.data(),
                                             Hplus.data(),
                                             Hminus.data(),
                                             Kplus.data(),
                                             Kminus.data(),
                                             Zplus.data(),
                                             Zminus.data() };
        Material::Deformation   materialDeformation{ dU_GPs.data(),
                                                   dSurface_strain_GPs.data(),
                                                   qp.normal.data(),
                                                   qp.separationVector.data() };
        Material::TimeIncrement materialTimeIncrement{ time, dT };

        qp.material->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );
      }
      else if constexpr ( nDim == 2 ) {
        Eigen::Vector3d                                forcePlus3d  = Eigen::Vector3d::Zero();
        Eigen::Vector3d                                forceMinus3d = Eigen::Vector3d::Zero();
        Eigen::Matrix< double, 9, 1 >                  surfaceStressPlus3d;
        Eigen::Matrix< double, 9, 1 >                  surfaceStressMinus3d;
        Eigen::Matrix< double, 6, 1 >                  dU3d;
        Eigen::Matrix< double, 18, 1 >                 dSurfaceStrain3d;
        Eigen::Vector3d                                normal3d = Eigen::Vector3d::Zero();
        Eigen::Matrix< double, 3, 3, Eigen::RowMajor > Qplus3d, Qminus3d;
        Eigen::Matrix< double, 9, 9, Eigen::RowMajor > Zplus3d, Zminus3d;
        Eigen::Matrix< double, 3, 9, Eigen::RowMajor > Hplus3d, Hminus3d;
        Eigen::Matrix< double, 9, 3, Eigen::RowMajor > Kplus3d, Kminus3d;
        Eigen::Vector3d                                separation3d = Eigen::Vector3d::Zero();

        surfaceStressPlus3d.setZero();
        surfaceStressMinus3d.setZero();
        dU3d.setZero();
        dSurfaceStrain3d.setZero();
        Qplus3d.setZero();
        Qminus3d.setZero();
        Zplus3d.setZero();
        Zminus3d.setZero();
        Hplus3d.setZero();
        Hminus3d.setZero();
        Kplus3d.setZero();
        Kminus3d.setZero();

        for ( int i = 0; i < nDim; ++i ) {
          forcePlus3d( i )  = forcePlus( i );
          forceMinus3d( i ) = forceMinus( i );
          normal3d( i )     = qp.normal( i );
          separation3d( i ) = qp.separationVector( i );
          dU3d( i )         = dU_GPs( i );
          dU3d( 3 + i )     = dU_GPs( nDim + i );

          for ( int j = 0; j < nDim; ++j ) {
            const int index2d = i * nDim + j;
            const int index3d = i * 3 + j;

            surfaceStressPlus3d( index3d )  = surfaceStressPlus( index2d );
            surfaceStressMinus3d( index3d ) = surfaceStressMinus( index2d );
            dSurfaceStrain3d( index3d )     = dSurface_strain_GPs( index2d );
            dSurfaceStrain3d( 9 + index3d ) = dSurface_strain_GPs( nTensor + index2d );
          }
        }

        Material::State         materialState{ forcePlus3d.data(),
                                       forceMinus3d.data(),
                                       surfaceStressPlus3d.data(),
                                       surfaceStressMinus3d.data(),
                                       qp.managedStateVars->materialStateVars.data() };
        Material::Tangents      materialTangents{ Qplus3d.data(),
                                             Qminus3d.data(),
                                             Hplus3d.data(),
                                             Hminus3d.data(),
                                             Kplus3d.data(),
                                             Kminus3d.data(),
                                             Zplus3d.data(),
                                             Zminus3d.data() };
        Material::Deformation   materialDeformation{ dU3d.data(),
                                                   dSurfaceStrain3d.data(),
                                                   normal3d.data(),
                                                   separation3d.data() };
        Material::TimeIncrement materialTimeIncrement{ time, dT };

        qp.material->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );

        for ( int i = 0; i < nDim; ++i ) {
          forcePlus( i )  = forcePlus3d( i );
          forceMinus( i ) = forceMinus3d( i );
          for ( int k = 0; k < nDim; ++k ) {
            Qplus( i, k )  = Qplus3d( i, k );
            Qminus( i, k ) = Qminus3d( i, k );
          }
        }

        for ( int i = 0; i < nDim; ++i ) {
          for ( int j = 0; j < nDim; ++j ) {
            const int row2d = i * nDim + j;
            const int row3d = i * 3 + j;

            surfaceStressPlus( row2d )  = surfaceStressPlus3d( row3d );
            surfaceStressMinus( row2d ) = surfaceStressMinus3d( row3d );

            for ( int k = 0; k < nDim; ++k ) {
              Kplus( row2d, k )  = Kplus3d( row3d, k );
              Kminus( row2d, k ) = Kminus3d( row3d, k );
            }

            for ( int k = 0; k < nDim; ++k ) {
              for ( int l = 0; l < nDim; ++l ) {
                const int col2d = k * nDim + l;
                const int col3d = k * 3 + l;

                Zplus( row2d, col2d )  = Zplus3d( row3d, col3d );
                Zminus( row2d, col2d ) = Zminus3d( row3d, col3d );
              }
            }
          }
        }

        for ( int i = 0; i < nDim; ++i ) {
          for ( int k = 0; k < nDim; ++k ) {
            for ( int l = 0; l < nDim; ++l ) {
              const int col2d = k * nDim + l;
              const int col3d = k * 3 + l;

              Hplus( i, col2d )  = Hplus3d( i, col3d );
              Hminus( i, col2d ) = Hminus3d( i, col3d );
            }
          }
        }
      }

      qp.managedStateVars->forcePlus          = forcePlus;
      qp.managedStateVars->forceMinus         = forceMinus;
      qp.managedStateVars->surfaceStressPlus  = surfaceStressPlus;
      qp.managedStateVars->surfaceStressMinus = surfaceStressMinus;
      qp.managedStateVars->displacement += dU_GPs;
      qp.managedStateVars->surfaceStrain += dSurface_strain_GPs;

      Pe -= Njump.transpose() * ( forcePlus + forceMinus ) * qp.J0xW;
      Pe -= BPlusFull.transpose() * surfaceStressPlus * qp.J0xW;
      Pe -= BMinusFull.transpose() * surfaceStressMinus * qp.J0xW;

      // No major symmetry is assumed. The cross blocks
      // d(surfaceStressPlus)/d(A-) and d(surfaceStressMinus)/d(A+) are
      // exactly zero (verified in TestMarmotXInterfaceMaterialHypoElastic)
      // and are therefore simply absent from this assembly.
      Ke += ( Njump.transpose() * Qplus * Njump + Njump.transpose() * Qminus * Njump +
              Njump.transpose() * Hplus * BPlusFull + Njump.transpose() * Hminus * BMinusFull +
              BPlusFull.transpose() * Kplus * Njump + BPlusFull.transpose() * Zplus * BPlusFull +
              BMinusFull.transpose() * Kminus * Njump + BMinusFull.transpose() * Zminus * BMinusFull ) *
            qp.J0xW;
    }
  }

  template < int nDim, int nNodes >
  void XInterfaceFiniteElement< nDim, nNodes >::setInitialConditions( StateTypes state, const double* values )
  {
    switch ( state ) {
    case MarmotElement::MarmotMaterialInitialization: {
      for ( QuadraturePoint& qp : qps ) {
        qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                         qp.managedStateVars->materialStateVars.size() );
      }
      break;
    }

    case MarmotElement::MarmotMaterialStateVars: {
      throw std::invalid_argument( "Please use initializeStateVars directly on material" );
    }

    default:
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                                << ": invalid initial condition for XInterfaceFiniteElement" );
    }
  }

  template < int nDim, int nNodes >
  void XInterfaceFiniteElement< nDim, nNodes >::computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                                                        double*                             P,
                                                                        double*                             K,
                                                                        const int                           elementFace,
                                                                        const double*                       load,
                                                                        const double*                       QTotal,
                                                                        double                              time,
                                                                        double                              dT )
  {
    throw std::invalid_argument(
      MakeString() << __PRETTY_FUNCTION__ << ": distributed loads are not implemented for XInterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  void XInterfaceFiniteElement< nDim, nNodes >::computeBodyForce( double*       P,
                                                                  double*       K,
                                                                  const double* load,
                                                                  const double* QTotal,
                                                                  double        time,
                                                                  double        dT )
  {
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                              << ": body forces are not implemented for XInterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  void XInterfaceFiniteElement< nDim, nNodes >::computeConsistentInertia( double* M )
  {
    throw std::runtime_error( MakeString()
                              << __PRETTY_FUNCTION__ << ": inertia is not implemented for XInterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  void XInterfaceFiniteElement< nDim, nNodes >::computeLumpedInertia( double* M )
  {
    throw std::runtime_error( MakeString()
                              << __PRETTY_FUNCTION__ << ": inertia is not implemented for XInterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  std::vector< double > XInterfaceFiniteElement< nDim, nNodes >::getCoordinatesAtCenter()
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

  template < int nDim, int nNodes >
  std::vector< std::vector< double > > XInterfaceFiniteElement< nDim, nNodes >::getCoordinatesAtQuadraturePoints()
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

  template < int nDim, int nNodes >
  int XInterfaceFiniteElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return static_cast< int >( qps.size() );
  }

  ///@}

} // namespace Marmot::Elements
