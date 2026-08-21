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
 * @file ZInterfaceFiniteElement.h
 * @brief Extended interface element carrying the gradient jump g_k as a
 * globally interpolated field (ZIQUAD4).
 *
 * Drop-in replacement for YIQUAD4 with the SAME node ordering, the same
 * material-property layout (E, nu, h, base-material properties) and the same
 * section input. The only difference in the model file is that the four
 * bottom (midsurface) nodes carry three extra degrees of freedom.
 *
 * Relation to the paper (Stathas & Neuner, interphase paper, section 2--4)
 * -----------------------------------------------------------------------
 * The generalised strain X_M of eq. (43) is UNCHANGED:
 *
 *   X_M = ( [du_i] , du^{S+}_{i,j} , du^{S-}_{i,j} ) ,   M = 1..21 ,
 *
 * and the face reconstruction of eq. (44) is unchanged, including the
 * oblique-connector correction in ell and d^T. What changes is the status of
 * the gradient jump
 *
 *   g_k = [[ du_{k,r} ]] n_r          (eq. 12)
 *
 * of eq. (12). YIQUAD4 fixes it pointwise from the local equilibrium
 * r_i(g) = 0 of eq. (45) and eliminates it by the Schur complement of eq. (47),
 * which requires the averaged acoustic tensor <Q> of eq. (13) to be
 * invertible. This element instead interpolates g_k with the midsurface Q1
 * shape functions and hands the resulting nodal unknowns to the global solve,
 * so that NO inverse of <Q> is ever formed.
 *
 * Which equations of the paper change, and how
 * --------------------------------------------
 *   eq. (10) expanded incremental energy   -- UNCHANGED. It precedes any
 *            statement about g and is what the material kernel evaluates,
 *            one loop over the two faces.
 *   eq. (13) <Q> g = -Phi                  -- NO LONGER SOLVED. It becomes
 *            the variational equation conjugate to the test field ghat.
 *   eq. (16) jump blocks in g              -- this is where the element now
 *            STOPS. The terms carrying ghat_i, which eq. (18) cancels
 *            pointwise, are retained and assembled.
 *   eq. (18) condensed incremental energy  -- DROPPED (contains <Q>^{-1}), and
 *            with it eq. (19) C^eff, eq. (21) A^(1..3), eq. (22) Q^eff,
 *            F^eff, H^eff.
 *   eq. (20) incremental interface energy  -- REPLACED by the unreduced block
 *            form: eq. (16) + eq. (17) + the membrane term of eq. (10).
 *   eq. (23) interface residual            -- EXTENDED. Its derivation drops
 *            the term (h/4)[t_i] ghat_i "since the condensation performed at
 *            the previous step enforced the continuity of the interface
 *            traction". With g independent that premise is gone, the term
 *            survives, and it IS the residual of the new field:
 *
 *              R_g = int_A (h/4) [ t_i ] ghat_i dA .
 *
 *   eqs. (46), (47), (49), (50), (51)      -- REMOVED. No implicit
 *            differentiation, no Schur complement, no gauge direction m^T, no
 *            SVD, no truncated inverse, no frozen per-point binary decision.
 *   eq. (48) degeneracy rate               -- KEPT, as the diagnosis that
 *            motivates all of the above.
 *
 * Element unknowns:  q_e = [ d(24) ; g(12) ] = 36.
 *   d : displacement, 3 per node, 8 nodes, YIQUAD4 order [u-(4 nodes); u+(4 nodes)].
 *   g : gradient jump, 3 per midsurface node, attached to the bottom nodes
 *       0..3 so that it is shared by adjacent interface elements.
 *
 * Residuals (dA = surface measure; the h factors live inside the material):
 *   R_d = int_A ( J^T f + B^{+T} S^+ + B^{-T} S^- ) dA - f_ext
 *   R_g = int_A N_g^T p_z dA ,     p_z = (h/4)[t_i] + (h/4) zeta Q^e dg
 *
 * Tangent: the four blocks K_dd, K_dg, K_gd, K_gg obtained from the material's
 * unreduced (K_xx, K_xz, K_zx, K_zz) by the geometric operator
 * B_X = [ J ; B^+ ; B^- ] (21x24) and N_g (3x12). K_gd is NOT assumed to be
 * K_dg^T: the paper's Phi_i and Psihat_k are independent objects unless the
 * bulk tangent has major symmetry.
 *
 * On the continuity of g: g is a strain-like jump quantity, and the theory
 * does not require it to be continuous along the interface. Interpolating it
 * with C0 shape functions is an added assumption. It is what couples
 * neighbouring stations and therefore what makes an isolated degenerate point
 * harmless, but it also smooths g across a localizing band. See the
 * regularization note in MarmotZInterfaceMaterialHypoElastic for the
 * complementary treatment of a patch that degenerates as a whole.
 */
#pragma once

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotZInterfaceMaterialHypoElastic.h"

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Marmot::Elements {

  /**
   * @class ZInterfaceFiniteElement
   * @brief 3D (QUAD4-based) extended interface element with a nodal
   * gradient-jump field g.
   */
  class ZInterfaceFiniteElement : public MarmotElement, public MarmotGeometryInterfaceElement< 3, 8 > {

  public:
    enum SectionType {
      Interface,
    };

    static constexpr int nDim            = 3;
    static constexpr int nNodes          = 8;
    static constexpr int nInterfaceNodes = 4;
    static constexpr int nTensor         = 9;                   // full 3x3 gradient/stress
    static constexpr int nX              = 21;                  // X_M of eq. (43)
    static constexpr int nZ              = 3;                   // g_k of eq. (12)
    static constexpr int nSideDofU       = 12;                  // 4 nodes * 3
    static constexpr int nDofD           = 2 * nSideDofU;       // 24
    static constexpr int nDofG           = 3 * nInterfaceNodes; // 12
    static constexpr int sizeLoadVector  = nDofD + nDofG;       // 36
    static constexpr int nCoordinates    = nNodes * nDim;

    static constexpr int offD = 0;
    static constexpr int offG = nDofD; // 24

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

    using BXSized = Eigen::Matrix< double, nX, nDofD >;   // 21x24, [ J ; B+ ; B- ]
    using NgSized = Eigen::Matrix< double, nDim, nDofG >; // 3x12

    using Material = MarmotZInterfaceMaterialHypoElastic;

    Eigen::Map< const Eigen::VectorXd > elementProperties;

    const int         elLabel;
    const SectionType sectionType;

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

      VectorDim separationVector;
      VectorDim tangentialSeparation;
      double    normalSeparation;

      NMatrixSized  NmatSide;
      BSurfaceSized BmatSide;

      NJumpMatrixSized NmatJump;    // J (3x24), maps to u+ - u-
      BAvgSurfaceSized BmatAverage; // 9x24

      BXSized BX;                   // [ J ; B+ ; B- ] (21x24)
      NgSized Ng;                   // 3x12

      /**
       * @brief Per-quadrature-point state: the generalized stresses of
       * eq. (23)--(24) plus the traction imbalance that eq. (23) is now
       * allowed to carry, and the embedded two-face constitutive state.
       */
      class QPStateVarManager : public MarmotStateVarVectorManager {

        inline const static auto layout = makeLayout( {
          { .name = "generalizedForce", .length = 3 },
          { .name = "surfaceStressPlus", .length = 9 },
          { .name = "surfaceStressMinus", .length = 9 },
          { .name = "tractionImbalance", .length = 3 },
          { .name = "normalGradientJump", .length = 3 },
          { .name = "state block alignment padding", .length = ( 4 - ( ( 3 + 9 + 9 + 3 + 3 ) % 4 ) ) % 4 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< Eigen::Matrix< double, 3, 1 > > generalizedForce;
        Eigen::Map< Eigen::Matrix< double, 9, 1 > > surfaceStressPlus;
        Eigen::Map< Eigen::Matrix< double, 9, 1 > > surfaceStressMinus;
        Eigen::Map< Eigen::Matrix< double, 3, 1 > > tractionImbalance;
        Eigen::Map< Eigen::Matrix< double, 3, 1 > > normalGradientJump;

        Eigen::Map< Eigen::VectorXd > materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : MarmotStateVarVectorManager( theStateVarVector, layout ),
            generalizedForce( &find( "generalizedForce" ) ),
            surfaceStressPlus( &find( "surfaceStressPlus" ) ),
            surfaceStressMinus( &find( "surfaceStressMinus" ) ),
            tractionImbalance( &find( "tractionImbalance" ) ),
            normalGradientJump( &find( "normalGradientJump" ) ),
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
          BmatAverage( BAvgSurfaceSized::Zero() ),
          BX( BXSized::Zero() ),
          Ng( NgSized::Zero() )
      {
      }
    };

    std::vector< QuadraturePoint > qps;

    ZInterfaceFiniteElement( int                                         elementID,
                             FiniteElement::Quadrature::IntegrationTypes integrationType,
                             SectionType                                 sectionType = SectionType::Interface );

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

    void applyMaterialSettings( QuadraturePoint& qp );

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
