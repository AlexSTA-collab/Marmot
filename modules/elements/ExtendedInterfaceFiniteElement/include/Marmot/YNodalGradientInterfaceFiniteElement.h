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
 * @file YNodalGradientInterfaceFiniteElement.h
 * @brief Partial-mixed interface finite element with a nodal common-traction
 * field (YIQUAD4_NODALGRADIENT).
 *
 * This element keeps the existing generalized-interface kinematics: the
 * tangential average/jump surface gradients (Abar, DeltaA) stay
 * displacement-derived, and the local normal-gradient jump gamma is still
 * reconstructed from traction equilibrium inside the constitutive kernel.
 * It introduces ONLY two extra globally-shared nodal fields on the midsurface:
 *
 *   g_i : the average normal gradient   <u_{i,r}> N_r   (replaces [u]/h),
 *   t_i : the complementary common traction sigma_ij^+ N_j = sigma_ij^- N_j.
 *
 * This is a partial Hu-Washizu formulation. The multiplier t is NOT
 * eliminated -- it is an independent field -- so the compatibility equation
 * [u] = h g is enforced UNWEIGHTED (never weighted by the material tangent).
 * The resulting element tangent is an indefinite, generally non-symmetric
 * saddle-point system; the g-t block stays invertible whenever the coupling
 * C = h int N_t^T N_g has full rank, even where the material normal stiffness
 * H_gg loses rank in perfect plasticity. Because adjacent interface elements
 * share the same nodal t, the traction field is continuous along the
 * interface -- the intended mechanism for suppressing element-to-element
 * traction oscillations.
 *
 * Element unknowns:  q_e = [ d(24) ; g(12) ; t(12) ]  = 48.
 *   d : displacement, 3 per node, 8 nodes (existing YIQUAD4 order [u-;u+]).
 *   g : average normal gradient, 3 per midsurface node, 4 nodes.
 *   t : common traction,         3 per midsurface node, 4 nodes.
 * The g,t fields are attached to the bottom nodes (0..3) as the midsurface
 * representatives (shared across adjacent interface elements).
 *
 * Residuals (spec section 7; hFac = interface thickness h, dA = surface measure):
 *   R_d = hFac int B_a^T [s_Abar; s_DeltaA] dA + int J^T (N_t p) dA - f_e
 *   R_q = hFac int N_g^T (t^mat - N_t p) dA
 *   R_p = int N_t^T (J d - hFac N_g q) dA
 * Tangent (spec section 10): the nine blocks K_dd,K_dq,K_dp,K_qd,K_qq,K_qp,
 * K_pd,K_pq, and K_pp = 0. Consistent, using the reduced constitutive tangent
 * H = [[H_aa,H_ag],[H_ga,H_gg]] returned by
 * MarmotPartialMixedInterfaceMaterialHypoElastic.
 */
#pragma once

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotPartialMixedInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotStateVarVectorManager.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Marmot::Elements {

  /**
   * @class YNodalGradientInterfaceFiniteElement
   * @brief 3D (QUAD4-based) partial-mixed interface element with nodal
   * average-normal-gradient (g) and common-traction (t) fields.
   */
  class YNodalGradientInterfaceFiniteElement : public MarmotElement, public MarmotGeometryInterfaceElement< 3, 8 > {

  public:
    enum SectionType {
      Interface,
    };

    static constexpr int nDim            = 3;
    static constexpr int nNodes          = 8;
    static constexpr int nInterfaceNodes = 4;
    static constexpr int nTensor         = 9;                     // full 3x3 gradient/stress
    static constexpr int nA              = 18;                    // a = (Abar[9], DeltaA[9])
    static constexpr int nSideDofU       = 12;                    // 4 nodes * 3
    static constexpr int nDofD           = 2 * nSideDofU;         // 24
    static constexpr int nDofG           = 3 * nInterfaceNodes;   // 12
    static constexpr int nDofT           = 3 * nInterfaceNodes;   // 12
    static constexpr int sizeLoadVector  = nDofD + nDofG + nDofT; // 48
    static constexpr int nCoordinates    = nNodes * nDim;

    static constexpr int offD = 0;
    static constexpr int offG = nDofD;         // 24
    static constexpr int offT = nDofD + nDofG; // 36

    using ParentGeometryElement = MarmotGeometryInterfaceElement< 3, 8 >;

    using XiSized              = ParentGeometryElement::XiSized;
    using NSized               = ParentGeometryElement::NSized;
    using dNdXiSized           = ParentGeometryElement::dNdXiSized;
    using SurfaceJacobianSized = ParentGeometryElement::SurfaceJacobianSized;
    using MetricSized          = ParentGeometryElement::MetricSized;
    using GradSized            = ParentGeometryElement::GradSized;

    using VectorDim = ParentGeometryElement::VectorDim;
    using TensorDim = ParentGeometryElement::TensorDim;

    using NMatrixSized     = ParentGeometryElement::NMatrixSized;
    using NJumpMatrixSized = ParentGeometryElement::NJumpMatrixSized; // 3x24 = J

    using BSurfaceSized    = ParentGeometryElement::BSurfaceSized;    // 9x12
    using BAvgSurfaceSized = ParentGeometryElement::BAvgSurfaceSized; // 9x24

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;

    using BaSized  = Eigen::Matrix< double, nA, nDofD >;   // 18x24
    using NgSized  = Eigen::Matrix< double, nDim, nDofG >; // 3x12
    using Sa9Sized = Eigen::Matrix< double, nTensor, 1 >;  // 9

    using Material = MarmotPartialMixedInterfaceMaterialHypoElastic;

    Eigen::Map< const Eigen::VectorXd > elementProperties;

    const int         elLabel;
    const SectionType sectionType;

    struct QuadraturePoint {
      const XiSized xi;
      const double  weight;

      double detJ;
      double sqrtDetG;
      double dA;   // surface measure weight*sqrtDetG (NO thickness factor)
      double hFac; // interface thickness h

      NSized               N;
      dNdXiSized           dNdXi;
      SurfaceJacobianSized J;
      MetricSized          G;
      GradSized            gradN;

      VectorDim normal;
      TensorDim normalProjection;
      TensorDim tangentProjection;

      VectorDim separationVector;

      NMatrixSized  NmatSide;
      BSurfaceSized BmatSide;

      NJumpMatrixSized NmatJump;    // J (3x24)
      BAvgSurfaceSized BmatAverage; // B_Abar (9x24)

      BaSized Ba;                   // [B_Abar ; B_dA] (18x24)
      NgSized Ng;                   // 3x12 (= N_t as well)

      /**
       * @brief Named per-QP state: diagnostic generalized stresses plus the
       * embedded constitutive-kernel state (two bulk states, gamma,
       * sigma+/sigma-). The nodal g,t fields are DOFs, not stored here.
       */
      class QPStateVarManager : public MarmotStateVarVectorManager {

        inline const static auto layout = makeLayout( {
          { .name = "commonTraction", .length = 3 },
          { .name = "surfaceStressAvg", .length = 9 },
          { .name = "surfaceStressJump", .length = 9 },
          { .name = "normalGradientAtQp", .length = 3 },
          { .name = "state block alignment padding", .length = ( 4 - ( ( 3 + 9 + 9 + 3 ) % 4 ) ) % 4 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< Eigen::Matrix< double, 3, 1 > > commonTraction;
        Eigen::Map< Eigen::Matrix< double, 9, 1 > > surfaceStressAvg;
        Eigen::Map< Eigen::Matrix< double, 9, 1 > > surfaceStressJump;
        Eigen::Map< Eigen::Matrix< double, 3, 1 > > normalGradientAtQp;

        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            commonTraction( &find( "commonTraction" ) ),
            surfaceStressAvg( &find( "surfaceStressAvg" ) ),
            surfaceStressJump( &find( "surfaceStressJump" ) ),
            normalGradientAtQp( &find( "normalGradientAtQp" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
        {
        }
      };

      std::unique_ptr< QPStateVarManager > managedStateVars;
      std::unique_ptr< Material >          material;

      int getNumberOfRequiredStateVarsQuadraturePointOnly()
      {
        return QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly();
      }

      int getNumberOfRequiredStateVars()
      {
        return getNumberOfRequiredStateVarsQuadraturePointOnly() + material->getNumberOfRequiredStateVars();
      }

      void assignStateVars( double* stateVars, int nStateVars )
      {
        managedStateVars = std::make_unique< QPStateVarManager >( stateVars, nStateVars );
      }

      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ),
          weight( weight ),
          detJ( 0.0 ),
          sqrtDetG( 0.0 ),
          dA( 0.0 ),
          hFac( 1.0 ),
          N( NSized::Zero() ),
          dNdXi( dNdXiSized::Zero() ),
          J( SurfaceJacobianSized::Zero() ),
          G( MetricSized::Zero() ),
          gradN( GradSized::Zero() ),
          normal( VectorDim::Zero() ),
          normalProjection( TensorDim::Zero() ),
          tangentProjection( TensorDim::Zero() ),
          separationVector( VectorDim::Zero() ),
          NmatSide( NMatrixSized::Zero() ),
          BmatSide( BSurfaceSized::Zero() ),
          NmatJump( NJumpMatrixSized::Zero() ),
          BmatAverage( BAvgSurfaceSized::Zero() ),
          Ba( BaSized::Zero() ),
          Ng( NgSized::Zero() )
      {
      }
    };

    std::vector< QuadraturePoint > qps;

    YNodalGradientInterfaceFiniteElement( int                                         elementID,
                                          FiniteElement::Quadrature::IntegrationTypes integrationType,
                                          SectionType sectionType = SectionType::Interface );

    int getNumberOfRequiredStateVars();

    std::vector< std::vector< std::string > > getNodeFields();

    std::vector< int > getDofIndicesPermutationPattern();

    int getNNodes() { return nNodes; }

    int getNSpatialDimensions() { return nDim; }

    int getNDofPerElement() { return sizeLoadVector; }

    std::string getElementShape() { return "hexa8"; }

    void assignStateVars( double* stateVars, int nStateVars );

    void assignProperty( const ElementProperties& marmotElementProperty );

    void assignProperty( const MarmotMaterialSection& marmotElementProperty );

    void assignMaterial( const std::string& materialName, const double* materialProperties, int nMaterialProperties );

    void assignNodeCoordinates( const double* coordinates );

    void initializeYourself();

    void setInitialConditions( StateTypes state, const double* values );

    void computeDistributedLoad( MarmotElement::DistributedLoadTypes loadType,
                                 double*                             P,
                                 double*                             K,
                                 const int                           elementFace,
                                 const double*                       load,
                                 const double*                       QTotal,
                                 double                              time,
                                 double                              dT );

    void computeBodyForce( double* P, double* K, const double* load, const double* QTotal, double time, double dT );

    void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT );

    void computeConsistentInertia( double* M );

    void computeLumpedInertia( double* M );

    StateView getStateView( const std::string& stateName, int qpNumber )
    {
      const auto& qp = qps[qpNumber];

      if ( qp.managedStateVars->contains( stateName ) ) {
        return qp.managedStateVars->getStateView( stateName );
      }

      return qp.material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
    }

    std::vector< double > getCoordinatesAtCenter();

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    int getNumberOfQuadraturePoints();
  };

} // namespace Marmot::Elements
