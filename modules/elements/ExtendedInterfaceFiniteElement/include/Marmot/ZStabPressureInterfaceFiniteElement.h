/* ---------------------------------------------------------------------
 *  marmot - MAteRialMOdellingToolbox
 *  Alexandros Stathas alexandros.stathas@boku.ac.at
 *  LGPL 2.1 or later; see LICENSE.md at the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

/**
 * @file ZStabPressureInterfaceFiniteElement.h
 * @brief ZIQUAD4_STABP -- gradient jump as a nodal field, with TWO stabilised
 * mixed pressures.
 *
 * ZIQUAD4 fixes the perfect-plasticity determinacy problem (the averaged
 * acoustic tensor is never inverted) but is measurably worse than YIQUAD4 on
 * the traction jump once the layer yields hard -- 4.9x in the stiff family,
 * unchanged by the gradient-jump regularization zeta, i.e. volumetric locking
 * rather than anything to do with z's determination.
 *
 * A SINGLE mixed pressure cannot fix it. With sigma~ = dev(sigma) - p I and one
 * p shared by both faces, two blocks vanish identically:
 *
 *     d r_p / dz = 0     ( tr sym(z (x) n) = z.n enters the faces as +-z.n/2
 *                          and cancels from the MEAN volumetric strain )
 *     d r_z / dp = 0     ( a shared p cancels from the traction JUMP )
 *
 * so the mode z carries -- the JUMP in volumetric strain -- is invisible to it.
 * This element therefore carries the pressure as a pair, p^pm = pbar +- [p]/2,
 * mirroring the (Abar, DeltaA) and (gbar, z) split the formulation already
 * uses. The two couplings that were missing then exist and are exact:
 *
 *     d r_z     / d[p] = -(h/4) n ,        d r_[pjump] / dz = n .
 *
 * [t] AND [p] ARE RELATED, AND THAT IS THE POINT
 * ----------------------------------------------
 * The normal traction jump is [t].n = n.[dev sigma].n - [p], so enforcing
 * r_z = 0 slaves [p] to n.[dev sigma].n. That does NOT over-determine the pair:
 * the 2x2 system in (z.n, [p]) has
 *
 *     det = (h/4) ( 1 + n.Q_dev.n / K )   ->   h/4
 *
 * as K -> infinity AND the deviatoric acoustic response degenerates -- the
 * double limit in which both the displacement-based and single-pressure forms
 * fail. The off-diagonals -(h/4) and 1 carry the solvability, not the
 * diagonals. Verified to five figures in
 * TestMarmotZStabPressureInterfaceMaterialHypoElastic.
 *
 * The difficulties are then attacked separately and do not interfere:
 *   z tangential   <- tangential r_z, deviatoric acoustic tensor, DEGENERATE
 *                     under perfect plasticity; handled by z being global + zeta
 *   (z.n, [p])     <- normal r_z + r_[pjump], saddle pair, ALWAYS regular
 *
 * STABILISATION
 * -------------
 * (u, pbar) and (z.n, [p]) are both equal-order Q1/Q1 on the same midsurface
 * nodes, so BOTH need stabilising. A Brezzi-Pitkaranta pressure-gradient term
 * is applied to EACH pressure field:
 *
 *   R_pbar = h [ N^T ( <tr eps> + pbar/K ) + (gamma h_e^2 / 2 mu) gradN^T gradN pbar^ ]
 *   R_pjmp = h [ N^T ( [tr eps]  + [p]/K ) + (gamma h_e^2 / 2 mu) gradN^T gradN [p]^  ]
 *
 * The scale is the MESH size h_e, so the term is O(h_e^2) consistent and does
 * not change the converged solution. gamma follows the same override order as
 * YIQUAD4_STABP_MINI: env MARMOT_STABP_GAMMA, then elementProperties[1], then
 * the default.
 *
 * The stabilising mechanism is the MINI displacement bubble, b = (1-xi^2)(1-eta^2),
 * amplitude beta (3), element-local and condensed. It enriches the AVERAGE
 * surface gradient equally on both faces, Abar_ij += beta_i (grad_s b)_j, so it
 * is simply one more contribution to x and couples to everything through the
 * blocks the kernel already returns. R_beta is NONLINEAR in beta (beta enters
 * the plastic return map), so it is solved by an inner Newton with a
 * backtracking line search before condensation, as in YIQUAD4_STABP_MINI.
 *
 * The Brezzi-Pitkaranta pressure-gradient term is retained but DEFAULTS OFF
 * (gamma = 0). On this element it is degenerate: grad_s N is the in-plane
 * gradient of a 2D surface, so gradNp^T gradNp has rank <= 2 per station, and
 * for the pressure JUMP there is no boundary condition to fix its constant
 * mode -- only the 1/K diagonal, which vanishes in the incompressible limit
 * the term exists to serve. MARMOT_STABP_GAMMA or elementProperties[1] turns it
 * on for ablation.
 *
 * Element unknowns: q_e = [ u(24) ; g(12) ; pbar(4) ; [p](4) ] = 44 DOF.
 */
#pragma once

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotZStabPressureInterfaceMaterialHypoElastic.h"

#include <Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Marmot::Elements {

  class ZStabPressureInterfaceFiniteElement : public MarmotElement, public MarmotGeometryInterfaceElement< 3, 8 > {

  public:
    enum SectionType { Interface };

    static constexpr int nDim            = 3;
    static constexpr int nNodes          = 8;
    static constexpr int nInterfaceNodes = 4;
    static constexpr int nTensor         = 9;
    static constexpr int nX              = 21;
    static constexpr int nZ              = 3;
    static constexpr int nSideDofU       = 12;
    static constexpr int nDofU           = 24;
    static constexpr int nDofG           = 12;
    static constexpr int nDofP           = 4;
    static constexpr int sizeLoadVector  = nDofU + nDofG + 2 * nDofP; // 44
    static constexpr int nCoordinates    = nNodes * nDim;

    static constexpr int nDofB = 3; /**< MINI bubble amplitude, condensed */

    static constexpr int offU  = 0;
    static constexpr int offG  = nDofU;         // 24
    static constexpr int offPm = nDofU + nDofG; // 36
    static constexpr int offPj = offPm + nDofP; // 40

    using ParentGeometryElement = MarmotGeometryInterfaceElement< 3, 8 >;
    using XiSized               = ParentGeometryElement::XiSized;
    using NSized                = ParentGeometryElement::NSized;
    using dNdXiSized            = ParentGeometryElement::dNdXiSized;
    using SurfaceJacobianSized  = ParentGeometryElement::SurfaceJacobianSized;
    using MetricSized           = ParentGeometryElement::MetricSized;
    using GradSized             = ParentGeometryElement::GradSized;
    using VectorDim             = ParentGeometryElement::VectorDim;
    using TensorDim             = ParentGeometryElement::TensorDim;
    using NMatrixSized          = ParentGeometryElement::NMatrixSized;
    using NJumpMatrixSized      = ParentGeometryElement::NJumpMatrixSized;
    using BSurfaceSized         = ParentGeometryElement::BSurfaceSized;
    using BAvgSurfaceSized      = ParentGeometryElement::BAvgSurfaceSized;

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;
    using BXSized       = Eigen::Matrix< double, nX, nDofU >;
    using NgSized       = Eigen::Matrix< double, nDim, nDofG >;
    using NpSized       = Eigen::Matrix< double, 1, nDofP >;

    using Material = MarmotZStabPressureInterfaceMaterialHypoElastic;

    static double defaultStabGamma()
    {
      static const double g = [] {
        const char* v = std::getenv( "MARMOT_STABP_GAMMA" );
        return v ? std::atof( v ) : 0.0; // degenerate here; MINI is the stabiliser
      }();
      return g;
    }

    /** Ablation: MARMOT_MINI_BUBBLE=OFF holds beta at zero and skips its
     *  condensation, leaving the element equal-order and unstabilised. */
    static bool bubbleDisabled()
    {
      static const bool d = [] {
        const char* v = std::getenv( "MARMOT_MINI_BUBBLE" );
        return v && std::string( v ) == "OFF";
      }();
      return d;
    }

    Eigen::Map< const Eigen::VectorXd > elementProperties;
    const int                           elLabel;
    const SectionType                   sectionType;
    double                              stabGamma = defaultStabGamma();

    struct QuadraturePoint {
      const XiSized xi;
      const double  weight;
      double        detJ = 0.0, sqrtDetG = 0.0, J0xW = 0.0, meshSize = 0.0;

      NSized                               N;
      dNdXiSized                           dNdXi;
      SurfaceJacobianSized                 J;
      MetricSized                          G;
      GradSized                            gradN;
      VectorDim                            normal, separationVector, tangentialSeparation;
      TensorDim                            normalProjection, tangentProjection;
      double                               normalSeparation = 0.0;
      double                               bubble           = 0.0; /**< b = (1-xi^2)(1-eta^2) */
      VectorDim                            gradBubble;             /**< surface gradient of b */
      NMatrixSized                         NmatSide;
      BSurfaceSized                        BmatSide;
      NJumpMatrixSized                     NmatJump;
      BAvgSurfaceSized                     BmatAverage;
      BXSized                              BX;
      NgSized                              Ng;
      NpSized                              Np;
      Eigen::Matrix< double, nDim, nDofP > gradNp; // surface gradient of the pressure shapes

      class QPStateVarManager : public MarmotStateVarVectorManager {
        inline const static auto layout = makeLayout( {
          { .name = "generalizedForce", .length = 3 },
          { .name = "surfaceStressPlus", .length = 9 },
          { .name = "surfaceStressMinus", .length = 9 },
          { .name = "tractionImbalance", .length = 3 },
          { .name = "normalGradientJump", .length = 3 },
          { .name = "pressureMean", .length = 1 },
          { .name = "pressureJump", .length = 1 },
          { .name = "bubbleAlpha", .length = 3 },
          { .name = "state block alignment padding", .length = ( 4 - ( ( 3 + 9 + 9 + 3 + 3 + 1 + 1 + 3 ) % 4 ) ) % 4 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< Eigen::Matrix< double, 3, 1 > > generalizedForce;
        Eigen::Map< Eigen::Matrix< double, 9, 1 > > surfaceStressPlus;
        Eigen::Map< Eigen::Matrix< double, 9, 1 > > surfaceStressMinus;
        Eigen::Map< Eigen::Matrix< double, 3, 1 > > tractionImbalance;
        Eigen::Map< Eigen::Matrix< double, 3, 1 > > normalGradientJump;
        Eigen::Map< Eigen::Matrix< double, 1, 1 > > pressureMean;
        Eigen::Map< Eigen::Matrix< double, 1, 1 > > pressureJump;
        Eigen::Map< Eigen::Matrix< double, 3, 1 > > bubbleAlpha;
        Eigen::Map< Eigen::VectorXd >               materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

        QPStateVarManager( double* sv, int n )
          : MarmotStateVarVectorManager( sv, layout ),
            generalizedForce( &find( "generalizedForce" ) ),
            surfaceStressPlus( &find( "surfaceStressPlus" ) ),
            surfaceStressMinus( &find( "surfaceStressMinus" ) ),
            tractionImbalance( &find( "tractionImbalance" ) ),
            normalGradientJump( &find( "normalGradientJump" ) ),
            pressureMean( &find( "pressureMean" ) ),
            pressureJump( &find( "pressureJump" ) ),
            bubbleAlpha( &find( "bubbleAlpha" ) ),
            materialStateVars( &find( "begin of material state" ),
                               n - getNumberOfRequiredStateVarsQuadraturePointOnly() )
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
      void assignStateVars( double* sv, int n ) { managedStateVars = std::make_unique< QPStateVarManager >( sv, n ); }

      QuadraturePoint( XiSized xi_, double w )
        : xi( xi_ ),
          weight( w ),
          N( NSized::Zero() ),
          dNdXi( dNdXiSized::Zero() ),
          J( SurfaceJacobianSized::Zero() ),
          G( MetricSized::Zero() ),
          gradN( GradSized::Zero() ),
          normal( VectorDim::Zero() ),
          separationVector( VectorDim::Zero() ),
          tangentialSeparation( VectorDim::Zero() ),
          normalProjection( TensorDim::Zero() ),
          tangentProjection( TensorDim::Zero() ),
          gradBubble( VectorDim::Zero() ),
          NmatSide( NMatrixSized::Zero() ),
          BmatSide( BSurfaceSized::Zero() ),
          NmatJump( NJumpMatrixSized::Zero() ),
          BmatAverage( BAvgSurfaceSized::Zero() ),
          BX( BXSized::Zero() ),
          Ng( NgSized::Zero() ),
          Np( NpSized::Zero() ),
          gradNp( Eigen::Matrix< double, nDim, nDofP >::Zero() )
      {
      }
    };

    std::vector< QuadraturePoint > qps;

    ZStabPressureInterfaceFiniteElement( int                                         elementID,
                                         FiniteElement::Quadrature::IntegrationTypes integrationType,
                                         SectionType sectionType = SectionType::Interface );

    int                                       getNumberOfRequiredStateVars();
    std::vector< std::vector< std::string > > getNodeFields();
    std::vector< int >                        getDofIndicesPermutationPattern();
    int                                       getNNodes() { return nNodes; }
    int                                       getNSpatialDimensions() { return nDim; }
    int                                       getNDofPerElement() { return sizeLoadVector; }
    std::string                               getElementShape() { return "hexa8"; }

    void assignStateVars( double* stateVars, int nStateVars );
    void assignProperty( const ElementProperties& );
    void assignProperty( const MarmotMaterialSection& );
    void assignMaterial( const std::string&, const double*, int );
    void assignNodeCoordinates( const double* );
    void applyMaterialSettings( QuadraturePoint& qp );
    void initializeYourself();
    void setInitialConditions( StateTypes, const double* );
    void computeDistributedLoad( MarmotElement::DistributedLoadTypes,
                                 double*,
                                 double*,
                                 const int,
                                 const double*,
                                 const double*,
                                 double,
                                 double );
    void computeBodyForce( double*, double*, const double*, const double*, double, double );
    void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT );
    void computeConsistentInertia( double* );
    void computeLumpedInertia( double* );

    StateView getStateView( const std::string& stateName, int qpNumber )
    {
      const auto& qp = qps[qpNumber];
      if ( qp.managedStateVars->contains( stateName ) )
        return qp.managedStateVars->getStateView( stateName );
      return qp.material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
    }

    std::vector< double >                getCoordinatesAtCenter();
    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();
    int                                  getNumberOfQuadraturePoints();
  };

} // namespace Marmot::Elements
