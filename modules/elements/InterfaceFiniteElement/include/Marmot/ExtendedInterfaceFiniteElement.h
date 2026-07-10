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
 * @file ExtendedInterfaceFiniteElement.h
 * @brief Extended interface finite element formulation for displacement-jump based interface mechanics.
 *
 * This file defines the templated class `Marmot::Elements::ExtendedInterfaceFiniteElement`
 * and its full in-header implementation. The element evaluates interface traction-like
 * quantities from displacement jumps and average surface kinematics and assembles
 * residual/tangent contributions at quadrature points.
 */
#pragma once

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotExtendedInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotStateVarVectorManager.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Marmot::Elements {

  namespace ExtendedInterfaceDebug {

    inline bool enabled()
    {
      return std::getenv( "MARMOT_EI_DEBUG" ) != nullptr ||
             std::getenv( "MARMOT_EI_DEBUG_ELEMENT_ASSEMBLY" ) != nullptr;
    }

    inline int envInt( const char* name, int defaultValue )
    {
      const char* value = std::getenv( name );
      return value == nullptr ? defaultValue : std::atoi( value );
    }

    inline double envDouble( const char* name, double defaultValue )
    {
      const char* value = std::getenv( name );
      return value == nullptr ? defaultValue : std::atof( value );
    }

    inline bool hasEnv( const char* name )
    {
      return std::getenv( name ) != nullptr;
    }

    inline bool coordinateFilterRequested()
    {
      return hasEnv( "MARMOT_EI_DEBUG_X" ) || hasEnv( "MARMOT_EI_DEBUG_Y" ) || hasEnv( "MARMOT_EI_DEBUG_Z" ) ||
             hasEnv( "MARMOT_EI_DEBUG_X_MIN" ) || hasEnv( "MARMOT_EI_DEBUG_X_MAX" ) ||
             hasEnv( "MARMOT_EI_DEBUG_Y_MIN" ) || hasEnv( "MARMOT_EI_DEBUG_Y_MAX" ) ||
             hasEnv( "MARMOT_EI_DEBUG_Z_MIN" ) || hasEnv( "MARMOT_EI_DEBUG_Z_MAX" );
    }

    inline bool axisAllowed( double value, const char* exact, const char* min, const char* max )
    {
      const double tolerance = envDouble( "MARMOT_EI_DEBUG_TOL", 1e-8 );

      if ( hasEnv( exact ) && std::abs( value - envDouble( exact, value ) ) > tolerance )
        return false;

      if ( hasEnv( min ) && value < envDouble( min, value ) - tolerance )
        return false;

      if ( hasEnv( max ) && value > envDouble( max, value ) + tolerance )
        return false;

      return true;
    }

    inline bool coordinatesAllowed( const std::vector< double >& center )
    {
      if ( !coordinateFilterRequested() )
        return true;

      if ( center.size() > 0 &&
           !axisAllowed( center[0], "MARMOT_EI_DEBUG_X", "MARMOT_EI_DEBUG_X_MIN", "MARMOT_EI_DEBUG_X_MAX" ) )
        return false;

      if ( center.size() > 1 &&
           !axisAllowed( center[1], "MARMOT_EI_DEBUG_Y", "MARMOT_EI_DEBUG_Y_MIN", "MARMOT_EI_DEBUG_Y_MAX" ) )
        return false;

      if ( center.size() > 2 &&
           !axisAllowed( center[2], "MARMOT_EI_DEBUG_Z", "MARMOT_EI_DEBUG_Z_MIN", "MARMOT_EI_DEBUG_Z_MAX" ) )
        return false;

      return true;
    }

    inline std::mutex& debugStreamMutex()
    {
      static std::mutex mutex;
      return mutex;
    }

    inline void writeDebugLine( const std::string& line )
    {
      std::lock_guard< std::mutex > lock( debugStreamMutex() );
      std::cerr << line << '\n';
    }

    inline bool shouldPrint( int elementLabel, int qpIndex, const std::vector< double >& center )
    {
      if ( !enabled() )
        return false;

      const char* elementFilter = std::getenv( "MARMOT_EI_DEBUG_ELEMENT" );
      if ( elementFilter != nullptr && std::atoi( elementFilter ) != elementLabel )
        return false;

      const char* qpFilter = std::getenv( "MARMOT_EI_DEBUG_QP" );
      if ( qpFilter != nullptr && std::atoi( qpFilter ) != qpIndex )
        return false;

      if ( !coordinatesAllowed( center ) )
        return false;

      static std::atomic< int > nPrinted{ 0 };
      const int                 maxPrints = envInt( "MARMOT_EI_DEBUG_MAX_CALLS", 50 );
      return nPrinted.fetch_add( 1 ) < maxPrints;
    }

  } // namespace ExtendedInterfaceDebug

  /**
   * @class ExtendedInterfaceFiniteElement
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes.
   * @brief Interface finite element with displacement-jump kinematics.
   *
   * The element combines geometric interface operators from
   * `MarmotGeometryInterfaceElement<nDim, nNodes>` with an interface-material
   * update (`MarmotExtendedInterfaceMaterialHypoElastic`) at each quadrature point.
   * The current formulation uses linearized, small-deformation kinematics.
   * It stores quadrature-point state variables and assembles:
   * - element residual vector,
   * - algorithmic tangent matrix,
   * - zero inertia terms (current formulation).
   */
  template < int nDim, int nNodes >
  class ExtendedInterfaceFiniteElement : public MarmotElement, public MarmotGeometryInterfaceElement< nDim, nNodes > {

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

    using BSurfaceSized     = typename ParentGeometryElement::BSurfaceSized;
    using BAvgSurfaceSized  = typename ParentGeometryElement::BAvgSurfaceSized;
    using BJumpSurfaceSized = typename ParentGeometryElement::BAvgSurfaceSized;

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;

    using ForceSized                = Eigen::Matrix< double, nDim, 1 >;
    using SurfaceStressSized        = Eigen::Matrix< double, nTensor, 1 >;
    using InterfaceDisplSized       = Eigen::Matrix< double, 2 * nDim, 1 >;
    using InterfaceSurfaceGradSized = Eigen::Matrix< double, 2 * nTensor, 1 >;

    using QMatrixSized            = Eigen::Matrix< double, nDim, nDim, Eigen::RowMajor >;
    using ZMatrixSized            = Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor >;
    using HMatrixSized            = Eigen::Matrix< double, nDim, nTensor, Eigen::RowMajor >;
    using YMatrixSized            = Eigen::Matrix< double, nTensor, nTensor, Eigen::RowMajor >;
    using SurfaceJumpUMatrixSized = Eigen::Matrix< double, nTensor, nDim, Eigen::RowMajor >;

    using Material = MarmotExtendedInterfaceMaterialHypoElastic;

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
      NJumpMatrixSized  NmatJump;
      BAvgSurfaceSized  BmatAverage;
      BJumpSurfaceSized BmatJump;

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
          { .name = "surface stress jump", .length = nDim * nDim },
          { .name = "displacement", .length = 2 * nDim },
          { .name = "surface strain", .length = 2 * nDim * nDim },
          // For nDim == 3, the material state starts after 40 entries.
          { .name   = "state block alignment padding",
            .length = ( 4 - ( ( nDim + ( nDim % 2 ) + 2 * nDim * nDim + 2 * nDim + 2 * nDim * nDim ) % 4 ) ) % 4 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< ForceSized >                force;
        Eigen::Map< SurfaceStressSized >        surfaceStress;
        Eigen::Map< SurfaceStressSized >        surfaceStressJump;
        Eigen::Map< InterfaceDisplSized >       displacement;
        Eigen::Map< InterfaceSurfaceGradSized > surfaceStrain;

        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            force( &find( "force" ) ),
            surfaceStress( &find( "surface stress" ) ),
            surfaceStressJump( &find( "surface stress jump" ) ),
            displacement( &find( "displacement" ) ),
            surfaceStrain( &find( "surface strain" ) ),
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
          NmatSide( NMatrixSized::Zero() ),
          BmatSide( BSurfaceSized::Zero() ),
          NmatJump( NJumpMatrixSized::Zero() ),
          BmatAverage( BAvgSurfaceSized::Zero() ),
          BmatJump( BJumpSurfaceSized::Zero() )
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
    ExtendedInterfaceFiniteElement( int                                         elementID,
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
   * @brief In-header method definitions for `ExtendedInterfaceFiniteElement`.
   */
  ///@{

  template < int nDim, int nNodes >
  ExtendedInterfaceFiniteElement< nDim, nNodes >::ExtendedInterfaceFiniteElement(
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
  int ExtendedInterfaceFiniteElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    return qps[0].getNumberOfRequiredStateVars() * qps.size();
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > ExtendedInterfaceFiniteElement< nDim, nNodes >::getNodeFields()
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
  std::vector< int > ExtendedInterfaceFiniteElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    static std::vector< int > permutationPattern;

    if ( permutationPattern.empty() ) {
      for ( int i = 0; i < nNodes * nDim; i++ )
        permutationPattern.push_back( i );
    }

    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void ExtendedInterfaceFiniteElement< nDim, nNodes >::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ ) {
      auto&   qp          = qps[i];
      double* qpStateVars = stateVars + ( i * nQpStateVars );
      qp.assignStateVars( qpStateVars, nQpStateVars );
    }
  }

  template < int nDim, int nNodes >
  void ExtendedInterfaceFiniteElement< nDim, nNodes >::assignProperty( const ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void ExtendedInterfaceFiniteElement< nDim, nNodes >::assignProperty( const MarmotMaterialSection& section )
  {
    for ( auto& qp : qps ) {
      qp.material = std::make_unique< MarmotExtendedInterfaceMaterialHypoElastic >( section.materialName,
                                                                                    section.materialProperties,
                                                                                    section.nMaterialProperties,
                                                                                    elLabel );
    }
  }

  template < int nDim, int nNodes >
  void ExtendedInterfaceFiniteElement< nDim, nNodes >::assignMaterial( const std::string& materialName,
                                                                       const double*      materialProperties,
                                                                       int                nMaterialProperties )
  {
    for ( auto& qp : qps ) {
      qp.material = std::make_unique< MarmotExtendedInterfaceMaterialHypoElastic >( materialName,
                                                                                    materialProperties,
                                                                                    nMaterialProperties,
                                                                                    elLabel );
    }
  }

  template < int nDim, int nNodes >
  void ExtendedInterfaceFiniteElement< nDim, nNodes >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void ExtendedInterfaceFiniteElement< nDim, nNodes >::initializeYourself()
  {
    const double thickness = elementProperties.size() > 0 ? elementProperties[0] : 1.0;

    for ( QuadraturePoint& qp : qps ) {
      const bool fullyProjectedB = ( nDim == 3 );
      const auto geom            = this->evaluateAt( qp.xi, 0, fullyProjectedB );

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
      qp.BmatJump.setZero();
      for ( int row = 0; row < nTensor; ++row ) {
        for ( int col = 0; col < nSideDofs; ++col ) {
          qp.BmatJump( row, col )             = -qp.BmatSide( row, col );
          qp.BmatJump( row, nSideDofs + col ) = qp.BmatSide( row, col );
        }
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
  void ExtendedInterfaceFiniteElement< nDim, nNodes >::computeKernels( const double* QTotal_,
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

    constexpr int halfSize           = nNodes * nDim / 2;
    const auto    debugElementCenter = ExtendedInterfaceDebug::enabled() ? getCoordinatesAtCenter()
                                                                         : std::vector< double >{};

    for ( size_t qpIndex = 0; qpIndex < qps.size(); ++qpIndex ) {
      QuadraturePoint& qp    = qps[qpIndex];
      const auto&      Nside = qp.NmatSide;
      const auto&      Bside = qp.BmatSide;
      const auto&      Njump = qp.NmatJump;
      const auto&      Bavg  = qp.BmatAverage;
      const auto&      Bjump = qp.BmatJump;

      const auto dQBottom = dQ.template segment< halfSize >( 0 );
      const auto dQTop    = dQ.template segment< halfSize >( halfSize );

      InterfaceDisplSized dU_GPs;
      dU_GPs.template segment< nDim >( 0 )    = Nside * dQTop;
      dU_GPs.template segment< nDim >( nDim ) = Nside * dQBottom;

      InterfaceSurfaceGradSized dSurface_strain_GPs;
      dSurface_strain_GPs.template segment< nTensor >( 0 )       = Bside * dQTop;
      dSurface_strain_GPs.template segment< nTensor >( nTensor ) = Bside * dQBottom;

      ForceSized         force               = qp.managedStateVars->force;
      SurfaceStressSized surface_stress      = qp.managedStateVars->surfaceStress;
      SurfaceStressSized surface_stress_jump = qp.managedStateVars->surfaceStressJump;

      QMatrixSized            forceJumpU;
      HMatrixSized            forceAverageSurfaceGradient;
      HMatrixSized            forceJumpSurfaceGradient;
      SurfaceJumpUMatrixSized averageSurfaceStressJumpU;
      ZMatrixSized            averageSurfaceStressAverageSurfaceGradient;
      ZMatrixSized            averageSurfaceStressJumpSurfaceGradient;
      SurfaceJumpUMatrixSized jumpSurfaceStressJumpU;
      ZMatrixSized            jumpSurfaceStressAverageSurfaceGradient;
      ZMatrixSized            jumpSurfaceStressJumpSurfaceGradient;

      forceJumpU.setZero();
      forceAverageSurfaceGradient.setZero();
      forceJumpSurfaceGradient.setZero();
      averageSurfaceStressJumpU.setZero();
      averageSurfaceStressAverageSurfaceGradient.setZero();
      averageSurfaceStressJumpSurfaceGradient.setZero();
      jumpSurfaceStressJumpU.setZero();
      jumpSurfaceStressAverageSurfaceGradient.setZero();
      jumpSurfaceStressJumpSurfaceGradient.setZero();

      const bool debugThisQP = ExtendedInterfaceDebug::shouldPrint( elLabel,
                                                                    static_cast< int >( qpIndex ),
                                                                    debugElementCenter );

      if constexpr ( nDim == 3 ) {
        Material::State         materialState{ force.data(),
                                       surface_stress.data(),
                                       surface_stress_jump.data(),
                                       qp.managedStateVars->materialStateVars.data() };
        Material::Tangents      materialTangents{ forceJumpU.data(),
                                             forceAverageSurfaceGradient.data(),
                                             forceJumpSurfaceGradient.data(),
                                             averageSurfaceStressJumpU.data(),
                                             averageSurfaceStressAverageSurfaceGradient.data(),
                                             averageSurfaceStressJumpSurfaceGradient.data(),
                                             jumpSurfaceStressJumpU.data(),
                                             jumpSurfaceStressAverageSurfaceGradient.data(),
                                             jumpSurfaceStressJumpSurfaceGradient.data() };
        Material::Deformation   materialDeformation{ dU_GPs.data(), dSurface_strain_GPs.data(), qp.normal.data() };
        Material::TimeIncrement materialTimeIncrement{ time, dT };

        qp.material->setDebugOutputForNextCall( debugThisQP );
        qp.material->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );
      }
      else if constexpr ( nDim == 2 ) {
        Eigen::Vector3d                                force3d = Eigen::Vector3d::Zero();
        Eigen::Matrix< double, 9, 1 >                  surfaceStress3d;
        Eigen::Matrix< double, 9, 1 >                  surfaceStressJump3d;
        Eigen::Matrix< double, 6, 1 >                  dU3d;
        Eigen::Matrix< double, 18, 1 >                 dSurfaceStrain3d;
        Eigen::Vector3d                                normal3d = Eigen::Vector3d::Zero();
        Eigen::Matrix< double, 3, 3, Eigen::RowMajor > forceJumpU3d;
        Eigen::Matrix< double, 3, 9, Eigen::RowMajor > forceAverageSurfaceGradient3d;
        Eigen::Matrix< double, 3, 9, Eigen::RowMajor > forceJumpSurfaceGradient3d;
        Eigen::Matrix< double, 9, 3, Eigen::RowMajor > averageSurfaceStressJumpU3d;
        Eigen::Matrix< double, 9, 9, Eigen::RowMajor > averageSurfaceStressAverageSurfaceGradient3d;
        Eigen::Matrix< double, 9, 9, Eigen::RowMajor > averageSurfaceStressJumpSurfaceGradient3d;
        Eigen::Matrix< double, 9, 3, Eigen::RowMajor > jumpSurfaceStressJumpU3d;
        Eigen::Matrix< double, 9, 9, Eigen::RowMajor > jumpSurfaceStressAverageSurfaceGradient3d;
        Eigen::Matrix< double, 9, 9, Eigen::RowMajor > jumpSurfaceStressJumpSurfaceGradient3d;

        surfaceStress3d.setZero();
        surfaceStressJump3d.setZero();
        dU3d.setZero();
        dSurfaceStrain3d.setZero();
        forceJumpU3d.setZero();
        forceAverageSurfaceGradient3d.setZero();
        forceJumpSurfaceGradient3d.setZero();
        averageSurfaceStressJumpU3d.setZero();
        averageSurfaceStressAverageSurfaceGradient3d.setZero();
        averageSurfaceStressJumpSurfaceGradient3d.setZero();
        jumpSurfaceStressJumpU3d.setZero();
        jumpSurfaceStressAverageSurfaceGradient3d.setZero();
        jumpSurfaceStressJumpSurfaceGradient3d.setZero();

        for ( int i = 0; i < nDim; ++i ) {
          force3d( i )  = force( i );
          normal3d( i ) = qp.normal( i );
          dU3d( i )     = dU_GPs( i );
          dU3d( 3 + i ) = dU_GPs( nDim + i );

          for ( int j = 0; j < nDim; ++j ) {
            const int index2d = i * nDim + j;
            const int index3d = i * 3 + j;

            surfaceStress3d( index3d )      = surface_stress( index2d );
            surfaceStressJump3d( index3d )  = surface_stress_jump( index2d );
            dSurfaceStrain3d( index3d )     = dSurface_strain_GPs( index2d );
            dSurfaceStrain3d( 9 + index3d ) = dSurface_strain_GPs( nTensor + index2d );
          }
        }

        Material::State         materialState{ force3d.data(),
                                       surfaceStress3d.data(),
                                       surfaceStressJump3d.data(),
                                       qp.managedStateVars->materialStateVars.data() };
        Material::Tangents      materialTangents{ forceJumpU3d.data(),
                                             forceAverageSurfaceGradient3d.data(),
                                             forceJumpSurfaceGradient3d.data(),
                                             averageSurfaceStressJumpU3d.data(),
                                             averageSurfaceStressAverageSurfaceGradient3d.data(),
                                             averageSurfaceStressJumpSurfaceGradient3d.data(),
                                             jumpSurfaceStressJumpU3d.data(),
                                             jumpSurfaceStressAverageSurfaceGradient3d.data(),
                                             jumpSurfaceStressJumpSurfaceGradient3d.data() };
        Material::Deformation   materialDeformation{ dU3d.data(), dSurfaceStrain3d.data(), normal3d.data() };
        Material::TimeIncrement materialTimeIncrement{ time, dT };

        qp.material->setDebugOutputForNextCall( debugThisQP );
        qp.material->computeStress( materialState, materialTangents, materialDeformation, materialTimeIncrement );

        for ( int i = 0; i < nDim; ++i ) {
          force( i ) = force3d( i );

          for ( int j = 0; j < nDim; ++j ) {
            const int index2d = i * nDim + j;
            const int index3d = i * 3 + j;

            surface_stress( index2d )      = surfaceStress3d( index3d );
            surface_stress_jump( index2d ) = surfaceStressJump3d( index3d );
            forceJumpU( i, j )             = forceJumpU3d( i, j );

            for ( int k = 0; k < nDim; ++k ) {
              const int tensorCol2d = j * nDim + k;
              const int tensorCol3d = j * 3 + k;

              forceAverageSurfaceGradient( i, tensorCol2d ) = forceAverageSurfaceGradient3d( i, tensorCol3d );
              forceJumpSurfaceGradient( i, tensorCol2d )    = forceJumpSurfaceGradient3d( i, tensorCol3d );
              averageSurfaceStressJumpU( tensorCol2d, i )   = averageSurfaceStressJumpU3d( tensorCol3d, i );
              jumpSurfaceStressJumpU( tensorCol2d, i )      = jumpSurfaceStressJumpU3d( tensorCol3d, i );

              for ( int l = 0; l < nDim; ++l ) {
                const int tensorRow2d  = i * nDim + j;
                const int tensorRow3d  = i * 3 + j;
                const int tensorCol2d4 = k * nDim + l;
                const int tensorCol3d4 = k * 3 + l;

                averageSurfaceStressAverageSurfaceGradient( tensorRow2d,
                                                            tensorCol2d4 ) = averageSurfaceStressAverageSurfaceGradient3d( tensorRow3d,
                                                                                                                           tensorCol3d4 );
                averageSurfaceStressJumpSurfaceGradient( tensorRow2d,
                                                         tensorCol2d4 ) = averageSurfaceStressJumpSurfaceGradient3d( tensorRow3d,
                                                                                                                     tensorCol3d4 );
                jumpSurfaceStressAverageSurfaceGradient( tensorRow2d,
                                                         tensorCol2d4 ) = jumpSurfaceStressAverageSurfaceGradient3d( tensorRow3d,
                                                                                                                     tensorCol3d4 );
                jumpSurfaceStressJumpSurfaceGradient( tensorRow2d,
                                                      tensorCol2d4 ) = jumpSurfaceStressJumpSurfaceGradient3d( tensorRow3d,
                                                                                                               tensorCol3d4 );
              }
            }
          }
        }
      }

      qp.managedStateVars->force             = force;
      qp.managedStateVars->surfaceStress     = surface_stress;
      qp.managedStateVars->surfaceStressJump = surface_stress_jump;
      qp.managedStateVars->displacement += dU_GPs;
      qp.managedStateVars->surfaceStrain += dSurface_strain_GPs;

      const RhsSized peContribution = -Njump.transpose() * force * qp.J0xW -
                                      Bavg.transpose() * surface_stress * qp.J0xW -
                                      Bjump.transpose() * surface_stress_jump * qp.J0xW;
      const KeSizedMatrix keContribution = ( Njump.transpose() * forceJumpU * Njump +
                                             Njump.transpose() * forceAverageSurfaceGradient * Bavg +
                                             Njump.transpose() * forceJumpSurfaceGradient * Bjump +
                                             Bavg.transpose() * averageSurfaceStressJumpU * Njump +
                                             Bavg.transpose() * averageSurfaceStressAverageSurfaceGradient * Bavg +
                                             Bavg.transpose() * averageSurfaceStressJumpSurfaceGradient * Bjump +
                                             Bjump.transpose() * jumpSurfaceStressJumpU * Njump +
                                             Bjump.transpose() * jumpSurfaceStressAverageSurfaceGradient * Bavg +
                                             Bjump.transpose() * jumpSurfaceStressJumpSurfaceGradient * Bjump ) *
                                           qp.J0xW;

      Pe += peContribution;
      Ke += keContribution;

      if ( debugThisQP ) {
        const auto displacementJump = dU_GPs.template segment< nDim >( 0 ) - dU_GPs.template segment< nDim >( nDim );
        const auto averageSurfaceGradient = 0.5 * ( dSurface_strain_GPs.template segment< nTensor >( 0 ) +
                                                    dSurface_strain_GPs.template segment< nTensor >( nTensor ) );
        const auto surfaceGradientJump    = dSurface_strain_GPs.template segment< nTensor >( 0 ) -
                                         dSurface_strain_GPs.template segment< nTensor >( nTensor );

        std::ostringstream line;
        line << std::scientific << std::setprecision( 6 ) << "[Marmot EI element] el=" << elLabel << " qp=" << qpIndex
             << " center=";
        for ( size_t i = 0; i < debugElementCenter.size(); ++i ) {
          if ( i > 0 )
            line << ",";
          line << debugElementCenter[i];
        }
        line << " timeOld=" << time << " dT=" << dT << " J0xW=" << qp.J0xW << " normal=" << qp.normal.transpose()
             << " dQInf=" << dQ.template lpNorm< Eigen::Infinity >()
             << " displacementJump=" << displacementJump.transpose()
             << " averageSurfaceGradientNorm=" << averageSurfaceGradient.norm()
             << " surfaceGradientJumpNorm=" << surfaceGradientJump.norm() << " force=" << force.transpose()
             << " surfaceStressNorm=" << surface_stress.norm()
             << " surfaceStressJumpNorm=" << surface_stress_jump.norm()
             << " peContributionInf=" << peContribution.template lpNorm< Eigen::Infinity >()
             << " peGlobalInfAfter=" << Pe.template lpNorm< Eigen::Infinity >()
             << " keContributionInf=" << keContribution.template lpNorm< Eigen::Infinity >()
             << " forceJumpUNorm=" << forceJumpU.norm()
             << " avgStressAvgGradNorm=" << averageSurfaceStressAverageSurfaceGradient.norm()
             << " jumpStressJumpGradNorm=" << jumpSurfaceStressJumpSurfaceGradient.norm();
        ExtendedInterfaceDebug::writeDebugLine( line.str() );
      }
    }
  }

  template < int nDim, int nNodes >
  void ExtendedInterfaceFiniteElement< nDim, nNodes >::setInitialConditions( StateTypes state, const double* values )
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
                                                << ": invalid initial condition for ExtendedInterfaceFiniteElement" );
    }
  }

  template < int nDim, int nNodes >
  void ExtendedInterfaceFiniteElement< nDim, nNodes >::computeDistributedLoad(
    MarmotElement::DistributedLoadTypes loadType,
    double*                             P,
    double*                             K,
    const int                           elementFace,
    const double*                       load,
    const double*                       QTotal,
    double                              time,
    double                              dT )
  {
    throw std::invalid_argument( MakeString()
                                 << __PRETTY_FUNCTION__
                                 << ": distributed loads are not implemented for ExtendedInterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  void ExtendedInterfaceFiniteElement< nDim, nNodes >::computeBodyForce( double*       P,
                                                                         double*       K,
                                                                         const double* load,
                                                                         const double* QTotal,
                                                                         double        time,
                                                                         double        dT )
  {
    throw std::invalid_argument(
      MakeString() << __PRETTY_FUNCTION__ << ": body forces are not implemented for ExtendedInterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  void ExtendedInterfaceFiniteElement< nDim, nNodes >::computeConsistentInertia( double* M )
  {
    throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__
                                           << ": inertia is not implemented for ExtendedInterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  void ExtendedInterfaceFiniteElement< nDim, nNodes >::computeLumpedInertia( double* M )
  {
    throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__
                                           << ": inertia is not implemented for ExtendedInterfaceFiniteElement." );
  }

  template < int nDim, int nNodes >
  std::vector< double > ExtendedInterfaceFiniteElement< nDim, nNodes >::getCoordinatesAtCenter()
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
  std::vector< std::vector< double > > ExtendedInterfaceFiniteElement< nDim,
                                                                       nNodes >::getCoordinatesAtQuadraturePoints()
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
  int ExtendedInterfaceFiniteElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return static_cast< int >( qps.size() );
  }

  ///@}

} // namespace Marmot::Elements
