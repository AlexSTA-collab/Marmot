/* ---------------------------------------------------------------------
 *  Marmot - WarpingStabPressureMiniInterfaceFiniteElement
 *           (WIQUAD4_STABP_MINI)
 *  Alexandros Stathas alexandros.stathas@boku.ac.at
 * --------------------------------------------------------------------- */

/**
 * @file WarpingStabPressureMiniInterfaceFiniteElement.h
 * @brief Through-thickness-resolved interface element carrying BOTH the
 * symmetric/antisymmetric WARPING micro-fields AND a mixed, stabilised nodal
 * pressure with a MINI displacement bubble (WIQUAD4_STABP_MINI).
 *
 * WHY COMBINE THEM
 * ----------------
 * The two enrichments were each measured against the same resolved full-Cauchy
 * reference (stiff family, angle 10 deg) and each fixes an error the other does
 * not touch:
 *
 *   WIQUAD4 (warping)  displacement-jump error -38%, overall L2 -44% at
 *                      h = 0.01, convergence rate 0.644 -> 0.818; pressure
 *                      oscillation metric 1.400e-1 vs 1.516e-1 for GLIQUAD4,
 *                      i.e. essentially unchanged.
 *   MINI (mixed p)     pressure oscillation metric 1.227e-4, three orders of
 *                      magnitude better, and the only element whose traction
 *                      convergence rate is positive (+0.62 vs -0.51).
 *
 * The reason they compose is that they act on orthogonal parts of the response:
 * warping enriches the through-thickness profile of the DEVIATORIC strain, the
 * mixed pressure replaces the VOLUMETRIC constitutive relation. Neither
 * enrichment touches the face displacements, so both are element-internal and
 * condensable, and their internal blocks can be condensed in ONE Schur
 * complement.
 *
 * UNKNOWNS
 * --------
 *   global    q_e = [ u(24) ; p(4) ] = 28        (+1 DOF per midsurface node,
 *                                                 exactly as YIQUAD4_STABP_MINI)
 *   internal  b   = [ beta_mini(3) ; beta_Ws(9) ; beta_Wa(9) ] = 21, condensed
 *
 * The internal 21 are condensed jointly, not one after the other: the MINI
 * bubble and the warping amplitudes both feed the station strains, so their
 * K_bb blocks couple and a sequential condensation would be wrong.
 *
 * THE TWO INTERNAL FIELDS ARE NOT THE SAME BUBBLE
 * -----------------------------------------------
 * Both expansions use the midsurface bubble b = (1-xi^2)(1-eta^2), but they
 * enrich different things and are NOT redundant:
 *
 *   MINI bubble      ubar += b beta        -> adds beta (x) grad_s b EQUALLY to
 *                                             A+ and A-, i.e. it perturbs the
 *                                             MEAN in-plane displacement. This
 *                                             is the inf-sup bubble: it is what
 *                                             makes the equal-order Q1/Q1
 *                                             displacement/pressure pairing
 *                                             admissible.
 *   warping modes    u += phi_s(zeta) w_s + phi_a(zeta) w_a, each amplitude
 *                                             expanded in { b, b xi, b eta }
 *                                          -> adds grad_s w_s, grad_s w_a as
 *                                             separate generalized strains with
 *                                             their own through-thickness
 *                                             profiles. This is the
 *                                             through-thickness kinematics.
 *
 * A constant MINI bubble amplitude changes the mean gradient; a constant
 * warping amplitude produces Ws = Wa = 0 and is invisible. They span different
 * subspaces of the generalized strain.
 *
 * BOUNDARY-VANISHING IS THE CONSISTENCY CONDITION
 * -----------------------------------------------
 * Every internal mode vanishes on the whole element boundary, so
 * int_Ae grad_s M dA = closed_int_dAe M n ds = 0 and a uniform stress state
 * produces no internal residual (patch test). Measured on the warping element
 * with the non-bubble basis {xi, eta, xi eta}: 21% residual error under uniform
 * loading; with bubbles: 1.4e-17.
 *
 * PRESSURE STABILISATION
 * ----------------------
 * Brezzi-Pitkaranta, identical to YIQUAD4_STABP_MINI:
 *      R_p += h (gamma h_e^2 / 2 mu) gradN^T gradN p^
 * with h_e the MESH size, not a physical length, so the term is O(h_e^2)
 * consistent and does not change the converged solution.
 *
 * The material returns R_p ALREADY multiplied by the interface thickness h, so
 * the whole 40x40 generalized tangent is one consistent block; the element adds
 * only the stabilisation (which carries its own explicit h).
 */
#pragma once

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotWarpingStabPressureInterfaceMaterialHypoElastic.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Marmot::Elements {

  class WarpingStabPressureMiniInterfaceFiniteElement : public MarmotElement,
                                                        public MarmotGeometryInterfaceElement< 3, 8 > {

  public:
    enum SectionType { Interface };

    /** Five Gauss-Lobatto stations, the production rule of the warping family.
     * Fixed rather than templated: the class is a plain (non-template) element,
     * exactly like YIQUAD4_STABP_MINI, so every fixed-size Eigen block below
     * stays non-dependent and needs no `template` disambiguation. */
    static constexpr int nStations = 5;

    static constexpr int nDim            = 3;
    static constexpr int nNodes          = 8;
    static constexpr int nInterfaceNodes = 4;
    static constexpr int nTensor         = 9;
    static constexpr int nDofU           = 24;            // 8 nodes x 3
    static constexpr int nDofP           = 4;             // 4 midsurface nodes x 1
    static constexpr int sizeLoadVector  = nDofU + nDofP; // 28
    static constexpr int nCoordinates    = nNodes * nDim;
    static constexpr int halfU           = nDofU / 2;

    static constexpr int offU = 0;
    static constexpr int offP = nDofU;

    /** Internal, condensed unknowns. */
    static constexpr int nMiniDof   = nDim;                // 3
    static constexpr int nWarpModes = nInterfaceNodes - 1; // 3: { b, b xi, b eta }
    static constexpr int nWarpBlock = nWarpModes * nDim;   // 9 per warping field
    static constexpr int nWarpDof   = 2 * nWarpBlock;      // 18
    static constexpr int nInternal  = nMiniDof + nWarpDof; // 21

    static constexpr int offMini = 0;
    static constexpr int offWarp = nMiniDof;

    /** Columns of the generalized B operator: global unknowns, then internal. */
    static constexpr int nAll = sizeLoadVector + nInternal; // 49

    using Material = MarmotWarpingStabPressureInterfaceMaterialHypoElasticN< nStations >;

    static constexpr int nGen = Material::nGeneralized; // 40

    using ParentGeometryElement = MarmotGeometryInterfaceElement< 3, 8 >;
    using XiSized               = ParentGeometryElement::XiSized;
    using NSized                = ParentGeometryElement::NSized; // 1x4
    using dNdXiSized            = ParentGeometryElement::dNdXiSized;
    using SurfaceJacobianSized  = ParentGeometryElement::SurfaceJacobianSized;
    using MetricSized           = ParentGeometryElement::MetricSized;
    using GradSized             = ParentGeometryElement::GradSized; // 3x4
    using VectorDim             = ParentGeometryElement::VectorDim;
    using TensorDim             = ParentGeometryElement::TensorDim;
    using NMatrixSized          = ParentGeometryElement::NMatrixSized;
    using NJumpMatrixSized      = ParentGeometryElement::NJumpMatrixSized; // 3x24
    using BSurfaceSized         = ParentGeometryElement::BSurfaceSized;    // 9x12
    using BAvgSurfaceSized      = ParentGeometryElement::BAvgSurfaceSized; // 9x24

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;

    using BWarpSized     = Eigen::Matrix< double, nTensor, nWarpBlock >;
    using InternalVector = Eigen::Matrix< double, nInternal, 1 >;
    using InternalMatrix = Eigen::Matrix< double, nInternal, nInternal >;
    using BGenSized      = Eigen::Matrix< double, nGen, nAll >;

    Eigen::Map< const Eigen::VectorXd > elementProperties;
    const int                           elLabel;
    const SectionType                   sectionType;

    struct QuadraturePoint {
      const XiSized xi;
      const double  weight;
      double        sqrtDetG = 0.0, J0xW = 0.0, hElem = 0.0;

      NSized               N;
      dNdXiSized           dNdXi;
      SurfaceJacobianSized J;
      MetricSized          G;
      GradSized            gradN;
      VectorDim            normal, separationVector;
      TensorDim            normalProjection, tangentProjection;
      NMatrixSized         NmatSide;
      BSurfaceSized        BmatSide;
      NJumpMatrixSized     NmatJump;
      BAvgSurfaceSized     BmatAverage;

      VectorDim gradBubble; // surface gradient of b = (1-xi^2)(1-eta^2)
      double    bubble = 0.0;

      /** Surface gradients of the warping modes { b, b xi, b eta }, 3 x nWarpModes. */
      Eigen::Matrix< double, nDim, nWarpModes > gradWarpModes;

      class QPStateVarManager : public MarmotStateVarVectorManager {
        static constexpr int nRaw = nDim           // generalizedForce
                                    + nTensor      // surfaceStressPlus
                                    + nTensor      // surfaceStressMinus
                                    + nTensor      // warpingStressSymmetric
                                    + nTensor      // warpingStressAntisymmetric
                                    + 1            // volumetricResidual
                                    + 1            // internalResidualNorm (diagnostic, qps[0] only)
                                    + nInternal    // internalAmplitudes (element-level, qps[0] only)
                                    + 2 * nDim     // displacement
                                    + 2 * nTensor  // surfaceStrain
                                    + 2 * nTensor; // warpingStrain

        inline const static auto layout = makeLayout( {
          { .name = "generalizedForce", .length = nDim },
          { .name = "surfaceStressPlus", .length = nTensor },
          { .name = "surfaceStressMinus", .length = nTensor },
          { .name = "warpingStressSymmetric", .length = nTensor },
          { .name = "warpingStressAntisymmetric", .length = nTensor },
          { .name = "volumetricResidual", .length = 1 },
          { .name = "internalResidualNorm", .length = 1 },
          { .name = "internalAmplitudes", .length = nInternal },
          { .name = "displacement", .length = 2 * nDim },
          { .name = "surfaceStrain", .length = 2 * nTensor },
          { .name = "warpingStrain", .length = 2 * nTensor },
          { .name = "state block alignment padding", .length = ( 4 - ( nRaw % 4 ) ) % 4 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< Eigen::Matrix< double, nDim, 1 > >        generalizedForce;
        Eigen::Map< Eigen::Matrix< double, nTensor, 1 > >     surfaceStressPlus;
        Eigen::Map< Eigen::Matrix< double, nTensor, 1 > >     surfaceStressMinus;
        Eigen::Map< Eigen::Matrix< double, nTensor, 1 > >     warpingStressSymmetric;
        Eigen::Map< Eigen::Matrix< double, nTensor, 1 > >     warpingStressAntisymmetric;
        Eigen::Map< Eigen::Matrix< double, 1, 1 > >           volumetricResidual;
        Eigen::Map< Eigen::Matrix< double, 1, 1 > >           internalResidualNorm;
        Eigen::Map< Eigen::Matrix< double, nInternal, 1 > >   internalAmplitudes;
        Eigen::Map< Eigen::Matrix< double, 2 * nDim, 1 > >    displacement;
        Eigen::Map< Eigen::Matrix< double, 2 * nTensor, 1 > > surfaceStrain;
        Eigen::Map< Eigen::Matrix< double, 2 * nTensor, 1 > > warpingStrain;
        Eigen::Map< Eigen::VectorXd >                         materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

        QPStateVarManager( double* v, int n )
          : MarmotStateVarVectorManager( v, layout ),
            generalizedForce( &find( "generalizedForce" ) ),
            surfaceStressPlus( &find( "surfaceStressPlus" ) ),
            surfaceStressMinus( &find( "surfaceStressMinus" ) ),
            warpingStressSymmetric( &find( "warpingStressSymmetric" ) ),
            warpingStressAntisymmetric( &find( "warpingStressAntisymmetric" ) ),
            volumetricResidual( &find( "volumetricResidual" ) ),
            internalResidualNorm( &find( "internalResidualNorm" ) ),
            internalAmplitudes( &find( "internalAmplitudes" ) ),
            displacement( &find( "displacement" ) ),
            surfaceStrain( &find( "surfaceStrain" ) ),
            warpingStrain( &find( "warpingStrain" ) ),
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

      QuadraturePoint( XiSized xi_, double w_ ) : xi( xi_ ), weight( w_ )
      {
        N.setZero();
        dNdXi.setZero();
        J.setZero();
        G.setZero();
        gradN.setZero();
        normal.setZero();
        separationVector.setZero();
        normalProjection.setZero();
        tangentProjection.setZero();
        NmatSide.setZero();
        BmatSide.setZero();
        NmatJump.setZero();
        BmatAverage.setZero();
        gradBubble.setZero();
        gradWarpModes.setZero();
      }
    };

    std::vector< QuadraturePoint > qps;

    /** Brezzi-Pitkaranta coefficient; same override chain as YIQUAD4_STABP_MINI. */
    double stabGamma = defaultStabGamma();

    /** Diagnostic switches, mirroring the MARMOT_STABP_GAMMA hook. Values:
     *    "warping" -- suppress the two warping fields  (element becomes a
     *                 five-station MINI)
     *    "mini"    -- suppress the MINI displacement bubble
     * Intended for isolating which enrichment drives a convergence problem;
     * unset in production. */
    static const char* suppressedEnrichment()
    {
      static const char* value = std::getenv( "MARMOT_WSM_SUPPRESS" );
      return value;
    }

    static bool warpingSuppressed()
    {
      const char* v = suppressedEnrichment();
      return v && std::string( v ) == "warping";
    }

    static bool miniBubbleSuppressed()
    {
      const char* v = suppressedEnrichment();
      return v && std::string( v ) == "mini";
    }

    static double defaultStabGamma()
    {
      if ( const char* e = std::getenv( "MARMOT_STABP_GAMMA" ) ) {
        const double v = std::atof( e );
        if ( v > 0.0 ) {
          return v;
        }
      }
      return 0.2;
    }

    WarpingStabPressureMiniInterfaceFiniteElement( int                                         elementID,
                                                   FiniteElement::Quadrature::IntegrationTypes integrationType,
                                                   SectionType sectionType_ = Interface )
      : ParentGeometryElement(),
        elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
        elLabel( elementID ),
        sectionType( sectionType_ )
    {
      for ( const auto& qpi : FiniteElement::Quadrature::getGaussPointInfo( this->shape, integrationType ) ) {
        qps.emplace_back( qpi.xi, qpi.weight );
      }
    }

    int getNumberOfRequiredStateVars() { return qps[0].getNumberOfRequiredStateVars() * (int)qps.size(); }

    std::vector< std::vector< std::string > > getNodeFields()
    {
      static std::vector< std::vector< std::string > > nf;
      if ( nf.empty() ) {
        for ( int i = 0; i < nNodes; i++ ) {
          nf.push_back( { "displacement" } );
          if ( i < nInterfaceNodes ) {
            nf[i].push_back( "interfacePressure" );
          }
        }
      }
      return nf;
    }

    std::vector< int > getDofIndicesPermutationPattern()
    {
      static std::vector< int > perm;
      if ( perm.empty() ) {
        perm.resize( sizeLoadVector );
        // canonical: bottom node A -> 4 slots [u(3), p(1)] at 4A ; top node A -> 3 slots at 16+3A
        for ( int A = 0; A < nInterfaceNodes; A++ ) {
          for ( int c = 0; c < 3; c++ ) {
            perm[offU + 3 * A + c] = 4 * A + c;
          }
          perm[offP + A] = 4 * A + 3;
          for ( int c = 0; c < 3; c++ ) {
            perm[offU + halfU + 3 * A + c] = 16 + 3 * A + c;
          }
        }
      }
      return perm;
    }

    int         getNNodes() { return nNodes; }
    int         getNSpatialDimensions() { return nDim; }
    int         getNDofPerElement() { return sizeLoadVector; }
    std::string getElementShape() { return "hexa8"; }

    void assignStateVars( double* sv, int n )
    {
      const int nq = n / (int)qps.size();
      for ( size_t i = 0; i < qps.size(); i++ ) {
        qps[i].assignStateVars( sv + i * nq, nq );
      }
    }

    void assignProperty( const ElementProperties& p )
    {
      new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( p.elementProperties, p.nElementProperties );
      if ( p.nElementProperties > 1 && p.elementProperties[1] > 0.0 ) {
        stabGamma = p.elementProperties[1];
      }
    }

    void assignProperty( const MarmotMaterialSection& s )
    {
      for ( auto& qp : qps ) {
        qp.material = std::make_unique< Material >( s.materialName,
                                                    s.materialProperties,
                                                    s.nMaterialProperties,
                                                    elLabel );
      }
    }

    void assignMaterial( const std::string& name, const double* props, int nProps )
    {
      for ( auto& qp : qps ) {
        qp.material = std::make_unique< Material >( name, props, nProps, elLabel );
      }
    }

    void assignNodeCoordinates( const double* c ) { ParentGeometryElement::assignNodeCoordinates( c ); }

    /** Parametric derivatives of the warping modes { b, b xi, b eta }, b = (1-xi^2)(1-eta^2). */
    static Eigen::Matrix< double, 2, nWarpModes > warpModeParametricGradients( const XiSized& xi )
    {
      Eigen::Matrix< double, 2, nWarpModes > d = Eigen::Matrix< double, 2, nWarpModes >::Zero();

      const double x = xi( 0 ), e = xi( 1 );
      const double bx = 1.0 - x * x, be = 1.0 - e * e;

      // M0 = b
      d( 0, 0 ) = -2.0 * x * be;
      d( 1, 0 ) = -2.0 * e * bx;

      // M1 = b * xi
      d( 0, 1 ) = be * ( 1.0 - 3.0 * x * x );
      d( 1, 1 ) = -2.0 * x * e * bx;

      // M2 = b * eta
      d( 0, 2 ) = -2.0 * x * e * be;
      d( 1, 2 ) = bx * ( 1.0 - 3.0 * e * e );

      return d;
    }

    void initializeYourself()
    {
      const double thickness = elementProperties.size() > 0 ? elementProperties[0] : 1.0;

      for ( auto& qp : qps ) {
        const auto g         = this->evaluateAt( qp.xi, 0 );
        qp.N                 = g.N;
        qp.dNdXi             = g.dNdXi;
        qp.J                 = g.J;
        qp.G                 = g.G;
        qp.sqrtDetG          = g.sqrtDetG;
        qp.gradN             = g.gradN;
        qp.normal            = g.n;
        qp.normalProjection  = g.normalProjection;
        qp.tangentProjection = g.tangentProjection;
        qp.NmatSide          = g.NmatSide;
        qp.BmatSide          = g.BmatSide;
        qp.NmatJump          = g.NmatJump;
        qp.BmatAverage       = g.BmatAverage;

        const VectorDim xB   = qp.NmatSide * this->getSideCoordinates( 0 );
        const VectorDim xT   = qp.NmatSide * this->getSideCoordinates( 1 );
        qp.separationVector  = xT - xB;
        constexpr double tol = 1.0e-12;
        if ( qp.separationVector.norm() > tol && qp.separationVector.dot( qp.normal ) < 0.0 ) {
          qp.normal *= -1.0;
        }
        qp.normalProjection  = qp.normal * qp.normal.transpose();
        qp.tangentProjection = TensorDim::Identity() - qp.normalProjection;
        if ( qp.separationVector.norm() > tol && qp.separationVector.dot( qp.normal ) <= tol ) {
          throw std::invalid_argument( "WIQUAD4_STABP_MINI: non-positive normal separation." );
        }

        qp.J0xW  = qp.weight * qp.sqrtDetG * thickness;
        qp.hElem = std::sqrt( qp.sqrtDetG ); // in-plane mesh size

        // grad_s M = J G^-1 dM/dxi, the same chain rule the shape-function
        // gradients use.
        const Eigen::Matrix< double, 2, 2 > Ginv = qp.G.inverse();
        {
          const double                  x = qp.xi( 0 ), e = qp.xi( 1 );
          Eigen::Matrix< double, 2, 1 > dbdxi;
          qp.bubble     = ( 1.0 - x * x ) * ( 1.0 - e * e );
          dbdxi( 0 )    = -2.0 * x * ( 1.0 - e * e );
          dbdxi( 1 )    = -2.0 * e * ( 1.0 - x * x );
          qp.gradBubble = qp.J * Ginv * dbdxi;
        }
        qp.gradWarpModes = qp.J * Ginv * warpModeParametricGradients( qp.xi );

        if ( qp.material ) {
          qp.material->setCharacteristicElementLength( qp.hElem );
        }
      }
    }

    void setInitialConditions( StateTypes state, const double* )
    {
      if ( state != MarmotElement::MarmotMaterialInitialization ) {
        throw std::invalid_argument( "WIQUAD4_STABP_MINI: invalid initial condition." );
      }
      for ( auto& qp : qps ) {
        qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                         qp.managedStateVars->materialStateVars.size() );
      }
    }

    /** B operator of one warping amplitude field: (grad_s w)_{ij} = sum_m beta_{m,i} (grad_s M_m)_j. */
    static BWarpSized warpingB( const Eigen::Matrix< double, nDim, nWarpModes >& gradModes )
    {
      BWarpSized B = BWarpSized::Zero();
      for ( int m = 0; m < nWarpModes; ++m ) {
        for ( int i = 0; i < nDim; ++i ) {
          for ( int j = 0; j < nDim; ++j ) {
            B( i * nDim + j, m * nDim + i ) = gradModes( j, m );
          }
        }
      }
      return B;
    }

    /**
     * Remove the trace from a 9-component generalized strain block:
     * vec(X) -> vec( X - (tr X / 3) I ).
     *
     * WHY THE SYMMETRIC WARPING MUST BE DEVIATORIC (measured, not stylistic)
     * ---------------------------------------------------------------------
     * The symmetric warping profile has a NONZERO thickness mean,
     * sum_alpha lambda_alpha phi_s = 2/3, so an unprojected Ws contributes
     * (2/3) tr(Ws) to the thickness-averaged volumetric strain -- the very
     * quantity the mixed pressure constrains. Ws is ELEMENT-LOCAL and
     * condensed, so that hands every element a private, cheap way to absorb an
     * arbitrary element-wise volumetric strain, which is exactly the
     * inter-element coupling the nodal pressure field exists to supply. The
     * checkerboard mode becomes locally free again.
     *
     * Measured on the production benchmark (stiff / angle 10 / h = 0.01,
     * fy = 5): with the unprojected warping the outer Newton limit-cycles from
     * t ~ 0.16 with the pressure correction frozen at ||ddp||inf = 3.3e-2 while
     * BOTH residuals are converged (displacement 1.8e-6, pressure 1.5e-12) --
     * the signature of a null mode, not of a bad step. The internal Newton was
     * verified converged throughout (max |R_b| <= 7e-11), so it is the coupled
     * system, not the condensation. Suppressing the warping entirely
     * (MARMOT_WSM_SUPPRESS=warping) makes the same job run to completion.
     *
     * Projecting Ws onto the traceless part removes the mechanism exactly and
     * CONSISTENTLY -- it is a change of the strain map, applied once in the B
     * operator, so every derived quantity follows automatically:
     *   * tr(Ws) = 0, so the warping drops out of the volumetric residual;
     *   * P^T vec(I) = 0, so the pressure exerts no force on the warping
     *     amplitudes either. The coupling vanishes in both directions.
     * What survives is what the accuracy study actually measured the warping
     * doing: enriching the through-thickness profile of the DEVIATORIC strain.
     *
     * The ANTISYMMETRIC field is deliberately NOT projected: its profile has
     * sum_alpha lambda_alpha phi_a = 0, so it never entered the thickness-
     * averaged volumetric constraint in the first place and carries no such
     * mechanism. Projecting it would discard real kinematics for nothing.
     */
    static Eigen::Matrix< double, nTensor, nTensor > deviatoricProjection()
    {
      Eigen::Matrix< double, nTensor, nTensor > P = Eigen::Matrix< double, nTensor, nTensor >::Identity();
      for ( int i = 0; i < nDim; ++i ) {
        for ( int k = 0; k < nDim; ++k ) {
          P( i * nDim + i, k * nDim + k ) -= 1.0 / 3.0;
        }
      }
      return P;
    }

    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
    {
      Eigen::Map< const RhsSized > QTotal( QTotal_ ), dQ( dQ_ );
      Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
      Eigen::Map< RhsSized >       Pe( Pe_ );

      // The inner Newton on the internal amplitudes re-evaluates the material
      // several times and computeStress COMMITS, so every trial -- including the
      // final, committing one -- must restart from the same incoming state.
      std::vector< std::vector< double > > materialStateIn( qps.size() );
      for ( size_t i = 0; i < qps.size(); i++ ) {
        auto& m            = qps[i].managedStateVars->materialStateVars;
        materialStateIn[i] = std::vector< double >( m.data(), m.data() + m.size() );
      }

      struct Accumulated {
        RhsSized                                           Pe;
        KeSizedMatrix                                      Kdd;
        InternalVector                                     Rb;
        InternalMatrix                                     Kbb;
        Eigen::Matrix< double, nInternal, sizeLoadVector > Kbd;
        Eigen::Matrix< double, sizeLoadVector, nInternal > Kdb;
      };

      auto assembleAll = [&]( const InternalVector& beta, bool commit ) {
        Accumulated acc;
        acc.Pe.setZero();
        acc.Kdd.setZero();
        acc.Rb.setZero();
        acc.Kbb.setZero();
        acc.Kbd.setZero();
        acc.Kdb.setZero();

        for ( size_t iq = 0; iq < qps.size(); iq++ ) {
          auto& qp = qps[iq];

          std::copy( materialStateIn[iq].begin(),
                     materialStateIn[iq].end(),
                     qp.managedStateVars->materialStateVars.data() );

          BAvgSurfaceSized BPlus = BAvgSurfaceSized::Zero(), BMinus = BAvgSurfaceSized::Zero();
          BMinus.block< nTensor, halfU >( 0, 0 )    = qp.BmatSide;
          BPlus.block< nTensor, halfU >( 0, halfU ) = qp.BmatSide;

          // MINI: ubar += b*beta  =>  Abar_ij += beta_i (grad_s b)_j, added
          // equally to A+ and A- because it perturbs the MEAN displacement.
          Eigen::Matrix< double, nTensor, nMiniDof > Bbub = Eigen::Matrix< double, nTensor, nMiniDof >::Zero();
          for ( int i = 0; i < nDim; i++ ) {
            for ( int j = 0; j < nDim; j++ ) {
              Bbub( nDim * i + j, i ) = qp.gradBubble( j );
            }
          }

          BWarpSized       Bw    = warpingB( qp.gradWarpModes );
          const BWarpSized BwDev = deviatoricProjection() * Bw; // see deviatoricProjection()
          if ( warpingSuppressed() ) {
            Bw.setZero();
          }
          if ( miniBubbleSuppressed() ) {
            Bbub.setZero();
          }

          BGenSized Bgen                                                                            = BGenSized::Zero();
          Bgen.block< nDim, nDofU >( Material::offsetJump, offU )                                   = qp.NmatJump;
          Bgen.block< nTensor, nDofU >( Material::offsetSurfacePlus, offU )                         = BPlus;
          Bgen.block< nTensor, nMiniDof >( Material::offsetSurfacePlus, sizeLoadVector + offMini )  = Bbub;
          Bgen.block< nTensor, nDofU >( Material::offsetSurfaceMinus, offU )                        = BMinus;
          Bgen.block< nTensor, nMiniDof >( Material::offsetSurfaceMinus, sizeLoadVector + offMini ) = Bbub;
          Bgen.block< nTensor, nWarpBlock >( Material::offsetWarpingSymmetric,
                                             sizeLoadVector + offWarp )              = warpingSuppressed() ? Bw : BwDev;
          Bgen.block< nTensor, nWarpBlock >( Material::offsetWarpingAntisymmetric,
                                             sizeLoadVector + offWarp + nWarpBlock ) = Bw;
          Bgen.block< 1, nDofP >( Material::offsetPressure, offP )                   = qp.N;

          Eigen::Matrix< double, nAll, 1 > q;
          q.segment< sizeLoadVector >( 0 )         = dQ;
          q.segment< nInternal >( sizeLoadVector ) = beta;

          const Eigen::Matrix< double, nGen, 1 > dGen = Bgen * q;

          Eigen::Matrix< double, 6, 1 >  dU6 = Eigen::Matrix< double, 6, 1 >::Zero();
          Eigen::Matrix< double, 18, 1 > dA  = Eigen::Matrix< double, 18, 1 >::Zero();
          Eigen::Matrix< double, 18, 1 > dW  = Eigen::Matrix< double, 18, 1 >::Zero();

          // The material takes ( u+, u- ) and forms w = u+ - u-; the element
          // already has the jump, so put it entirely on the + side.
          dU6.segment< 3 >( 0 ) = dGen.segment< nDim >( Material::offsetJump );
          dA.segment< 9 >( 0 )  = dGen.segment< nTensor >( Material::offsetSurfacePlus );
          dA.segment< 9 >( 9 )  = dGen.segment< nTensor >( Material::offsetSurfaceMinus );
          dW.segment< 9 >( 0 )  = dGen.segment< nTensor >( Material::offsetWarpingSymmetric );
          dW.segment< 9 >( 9 )  = dGen.segment< nTensor >( Material::offsetWarpingAntisymmetric );

          Eigen::Matrix< double, nGen, 1 > p = Eigen::Matrix< double, nGen, 1 >::Zero();
          Eigen::Matrix< double, nGen, nGen, Eigen::RowMajor >
            K = Eigen::Matrix< double, nGen, nGen, Eigen::RowMajor >::Zero();

          typename Material::State         st( p.data(), qp.managedStateVars->materialStateVars.data() );
          typename Material::Tangents      tg( K.data() );
          typename Material::Deformation   df( dU6.data(),
                                             dA.data(),
                                             dW.data(),
                                             qp.normal.data(),
                                             qp.separationVector.data(),
                                             dGen( Material::offsetPressure ) );
          typename Material::TimeIncrement ti{ time, dT };

          qp.material->computeStress( st, tg, df, ti );

          if ( commit ) {
            qp.managedStateVars->generalizedForce           = p.segment< nDim >( Material::offsetJump );
            qp.managedStateVars->surfaceStressPlus          = p.segment< nTensor >( Material::offsetSurfacePlus );
            qp.managedStateVars->surfaceStressMinus         = p.segment< nTensor >( Material::offsetSurfaceMinus );
            qp.managedStateVars->warpingStressSymmetric     = p.segment< nTensor >( Material::offsetWarpingSymmetric );
            qp.managedStateVars->warpingStressAntisymmetric = p.segment< nTensor >(
              Material::offsetWarpingAntisymmetric );
            qp.managedStateVars->volumetricResidual( 0 ) = p( Material::offsetPressure );
            qp.managedStateVars->displacement.segment< nDim >( 0 ) += dGen.segment< nDim >( Material::offsetJump );
            qp.managedStateVars->surfaceStrain += dGen.segment< 2 * nTensor >( Material::offsetSurfacePlus );
            qp.managedStateVars->warpingStrain += dGen.segment< 2 * nTensor >( Material::offsetWarpingSymmetric );
          }

          const auto Bd = Bgen.block< nGen, sizeLoadVector >( 0, 0 );
          const auto Bb = Bgen.block< nGen, nInternal >( 0, sizeLoadVector );

          acc.Pe -= Bd.transpose() * p * qp.J0xW;
          acc.Rb += Bb.transpose() * p * qp.J0xW;

          acc.Kdd += Bd.transpose() * K * Bd * qp.J0xW;
          acc.Kdb += Bd.transpose() * K * Bb * qp.J0xW;
          acc.Kbd += Bb.transpose() * K * Bd * qp.J0xW;
          acc.Kbb += Bb.transpose() * K * Bb * qp.J0xW;

          // ---- Brezzi-Pitkaranta pressure-gradient stabilisation ----
          // Lives in the element because it needs the MESH size h_e. It is
          // O(h_e^2) consistent, so it does not change the converged solution.
          const double hInterface = qp.material->getInterfaceThickness();
          const double mu         = qp.material->getShearModulus();
          const double tau        = stabGamma * qp.hElem * qp.hElem / ( 2.0 * mu );

          const auto pressureTotal = QTotal.segment< nDofP >( offP );
          acc.Pe.segment< nDofP >( offP ) -= hInterface * tau *
                                             ( qp.gradN.transpose() * ( qp.gradN * pressureTotal ) ) * qp.J0xW;
          acc.Kdd.block< nDofP, nDofP >( offP, offP ) += hInterface * tau * ( qp.gradN.transpose() * qp.gradN ) *
                                                         qp.J0xW;
        }
        return acc;
      };

      // ---- inner Newton on the internal amplitudes ----
      // beta enters the strain and therefore the plastic return mapping, so
      // R_beta is nonlinear in beta. Backtracking line search + exception guard,
      // as required by both parent elements: without it an undamped step at
      // yield onset drives the local traction-equilibrium solve past
      // convergence and the StressUpdateFailed propagates out as a global
      // cutback.
      struct InternalSolution {
        InternalVector beta;
        Accumulated    acc;
        bool           converged;
      };

      auto solveInternal = [&]( const InternalVector& betaInitial ) {
        InternalVector beta      = betaInitial;
        Accumulated    acc       = assembleAll( beta, false );
        const double   tolerance = 1.0e-10 * std::max( 1.0, acc.Pe.norm() );

        for ( int iteration = 0; iteration < 30; iteration++ ) {
          const double residualNorm = acc.Rb.norm();
          if ( residualNorm <= tolerance ) {
            break;
          }

          // Levenberg-Marquardt, not plain Newton. Two things make the internal
          // problem hard enough to need it:
          //   * K_bb has exactly known null directions -- the (grad_s w_a)_{n j}
          //     components of the antisymmetric field, absorbed by the station
          //     gradients because sum_alpha lambda_alpha phi_a = 0, and (after
          //     the deviatoric projection) the volumetric direction of Ws. The
          //     residual is exactly orthogonal to both, so the undamped step is
          //     well defined -- but it is also unbounded along them.
          //   * R_b is only PIECEWISE smooth in beta: each of the
          //     nStations x nQuadraturePoints material points can switch between
          //     elastic and plastic as beta moves, and a plain Newton step then
          //     stalls with a line search that cannot improve. Measured on the
          //     skewed plastic element test: |R_b| frozen at 6.4e-4 instead of
          //     machine zero, which propagates straight into a 14% tangent error.
          // Damping starts at zero, so a well-behaved step costs nothing extra;
          // it is raised only when the line search fails. The CONDENSATION always
          // uses the undamped K_bb, so consistency is untouched.
          const double scale = std::max( 1.0, acc.Kbb.diagonal().cwiseAbs().maxCoeff() );

          bool accepted = false;
          for ( double damping = 0.0; damping <= 1.0 && !accepted;
                damping        = ( damping == 0.0 ? 1.0e-8 : damping * 1.0e2 ) ) {
            InternalMatrix damped = acc.Kbb;
            damped.diagonal().array() += damping * scale;
            const InternalVector dBeta = -damped.completeOrthogonalDecomposition().solve( acc.Rb );

            double alpha = 1.0;
            for ( int lineSearch = 0; lineSearch < 20; lineSearch++ ) {
              try {
                Accumulated trial = assembleAll( beta + alpha * dBeta, false );
                if ( trial.Rb.norm() <= tolerance || trial.Rb.norm() < residualNorm ) {
                  beta += alpha * dBeta;
                  acc      = trial;
                  accepted = true;
                  break;
                }
              }
              catch ( const std::exception& ) {
                // the trial drove the local material solve past convergence:
                // treat exactly like a non-improving step and halve
              }
              alpha *= 0.5;
            }
          }
          if ( !accepted ) {
            break;
          }
        }
        return InternalSolution{ beta, acc, acc.Rb.norm() <= tolerance };
      };

      InternalSolution solution = solveInternal( qps[0].managedStateVars->internalAmplitudes );

      // Warm starting from the previous step's amplitudes is what makes the
      // inner solve cheap, but it can also strand the Newton on the wrong side
      // of a plastic loading/unloading switch: the internal residual is only
      // piecewise smooth in beta, and a stalled line search then leaves R_b at
      // O(1e-3) instead of machine zero. Because the condensed residual
      // correction K_db K_bb^+ R_b is first-order, that stall shows up directly
      // as tangent/residual inconsistency and costs global Newton iterations.
      // Retry from scratch and keep whichever attempt got further -- the same
      // fresh-restart strategy the station material uses for its own local
      // traction-equilibrium solve.
      if ( !solution.converged && qps[0].managedStateVars->internalAmplitudes.norm() > 0.0 ) {
        const InternalSolution fresh = solveInternal( InternalVector::Zero() );
        if ( fresh.converged || fresh.acc.Rb.norm() < solution.acc.Rb.norm() ) {
          solution = fresh;
        }
      }

      const InternalVector beta = solution.beta;
      Accumulated          acc  = assembleAll( beta, true );

      // ---- static condensation of all 21 internal unknowns at once ----
      // Jointly, not sequentially: the MINI bubble and the warping amplitudes
      // both feed the station strains, so their internal blocks couple.
      const auto KbbPinv = acc.Kbb.completeOrthogonalDecomposition().pseudoInverse();

      Pe = acc.Pe + acc.Kdb * ( KbbPinv * acc.Rb );
      Ke = acc.Kdd - acc.Kdb * KbbPinv * acc.Kbd;

      qps[0].managedStateVars->internalAmplitudes = beta;
      // Diagnostic: how far the inner Newton actually got. The condensed
      // residual correction K_db K_bb^+ R_b is only first-order accurate, so a
      // non-converged R_b shows up directly as tangent/residual inconsistency.
      qps[0].managedStateVars->internalResidualNorm( 0 ) = acc.Rb.norm();
    }

    void computeDistributedLoad( MarmotElement::DistributedLoadTypes,
                                 double*,
                                 double*,
                                 const int,
                                 const double*,
                                 const double*,
                                 double,
                                 double )
    {
      throw std::invalid_argument( "WIQUAD4_STABP_MINI: distributed loads not implemented." );
    }
    void computeBodyForce( double*, double*, const double*, const double*, double, double )
    {
      throw std::invalid_argument( "WIQUAD4_STABP_MINI: body forces not implemented." );
    }
    void computeConsistentInertia( double* ) { throw std::runtime_error( "WIQUAD4_STABP_MINI: no inertia." ); }
    void computeLumpedInertia( double* ) { throw std::runtime_error( "WIQUAD4_STABP_MINI: no inertia." ); }

    StateView getStateView( const std::string& name, int qpNumber )
    {
      const auto& qp = qps[qpNumber];
      if ( qp.managedStateVars->contains( name ) ) {
        return qp.managedStateVars->getStateView( name );
      }
      return qp.material->getStateView( name, qp.managedStateVars->materialStateVars.data() );
    }

    std::vector< double > getCoordinatesAtCenter()
    {
      std::vector< double >   c( nDim );
      Eigen::Map< VectorDim > m( c.data() );
      m = this->NMatrix( this->N( XiSized::Zero() ) ) * this->getSideCoordinates( 0 );
      return c;
    }

    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints()
    {
      std::vector< std::vector< double > > out;
      for ( const auto& qp : qps ) {
        std::vector< double >   c( nDim );
        Eigen::Map< VectorDim > m( c.data() );
        m = qp.NmatSide * this->getSideCoordinates( 0 );
        out.push_back( c );
      }
      return out;
    }

    int getNumberOfQuadraturePoints() { return (int)qps.size(); }
  };

} // namespace Marmot::Elements
