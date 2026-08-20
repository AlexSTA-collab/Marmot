/* ---------------------------------------------------------------------
 *  Marmot - YStabPressureMiniInterfaceFiniteElement  (YIQUAD4_STABP_MINI)
 *  Alexandros Stathas alexandros.stathas@boku.ac.at
 * --------------------------------------------------------------------- */

/**
 * @file YStabPressureMiniInterfaceFiniteElement.h
 * @brief Two-station equilibrated interface element with a MIXED, STABILISED
 * nodal pressure (YIQUAD4_STABP_MINI).
 *
 * MOTIVATION (measured, not assumed)
 * ----------------------------------
 * At H/E -> 0 the von Mises interface layer reaches the incompressible limit
 * (K/mu_plastic ~ 1e7). A displacement-only formulation then evaluates
 *      p = K * tr(eps)
 * with tr(eps) ~ 1e-4 against a deviatoric strain of ~0.55 -- four orders of
 * magnitude smaller -- so multiplying by K ~ 3.3e5 amplifies discretisation
 * noise into a spurious element-to-element pressure checkerboard. Measured on
 * the production benchmark: pressure p2p = 37.98 while the von Mises stress is
 * uniform to 5 significant figures (5.055 in all 900 elements) and the
 * deviatoric stress p2p is 1.3e-3. A three-point Poisson-ratio test gave
 * p2p ~ K^0.70 over a 10x range in K with the deviatoric solution bit-identical,
 * confirming pure amplification rather than any mechanical effect.
 *
 * Note what this is NOT: it is not a bending deficiency (the element retains
 * 64% of its elastic bending stiffness when fully yielded, vs 81% for a
 * resolved fiber section) and not a material instability (the second variation
 * of the interface functional is positive definite at every wavenumber and every
 * H; the wrinkle-mode stiffness changes only 21% for a 40000x change in H).
 *
 * FORMULATION
 * -----------
 * The pressure becomes an independent, C0-continuous nodal field (one scalar per
 * midsurface node), so it is a primary unknown instead of K x (near-zero):
 *      sigma~ = dev( sigma(eps) ) - p_h I ,        p_h = N p^
 * Displacement equilibrium is unchanged in structure (it just uses sigma~), and
 * the volumetric constitutive law is enforced weakly with a Brezzi-Pitkaranta
 * pressure-gradient stabilisation:
 *      R_u = W^T f + B+^T S+ + B-^T S-
 *      R_p = h [ N^T ( <tr eps> + p_h/K )  +  (gamma h_e^2 / 2 mu) gradN^T gradN p^ ]
 * The stabilisation scale is the MESH SIZE h_e, not a physical internal length,
 * and the term is consistent (it vanishes as h_e -> 0), so it does not change
 * the converged solution -- it only removes the spurious pressure mode. This is
 * the standard cure for an equal-order (Q1/Q1) displacement/pressure pairing.
 *
 * Element unknowns: q_e = [ u(24) ; p(4) ] = 28 DOF, i.e. only +1 DOF per
 * midsurface node over the standard YIQUAD4.
 *
 * Why an element-local fix cannot work (and why B-bar did nothing): the
 * checkerboard is an element-to-element mode, so any quantity computed inside a
 * single element is blind to it. Mean-dilatation B-bar yields a Q1/P0 pairing,
 * and Q1/P0 is precisely the pairing that ADMITS the checkerboard -- consistent
 * with the measured zero effect of the B-bar variant. Nodal continuity of p is
 * what supplies the inter-element coupling here.
 *
 * RESTRICTION: overriding the constitutive pressure is admissible only because
 * the von Mises yield surface is pressure-INDEPENDENT (see the material header).
 */
#pragma once

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotStabPressureInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotStateVarVectorManager.h"

#include <Eigen/Dense>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Marmot::Elements {

  namespace MiniProfiling {

    /**
     * Env-gated instrumentation for the element-local work of
     * YIQUAD4_STABP_MINI. Enable with MARMOT_MINI_PROFILE=1; the summary is
     * printed at process exit. Disabled it costs one predictable branch per
     * element call, so it can stay in the production header.
     *
     * The clock reads sit around the MATERIAL call and around the small dense
     * solves, so the split between "constitutive" and "everything else the
     * element does" is measured rather than inferred.
     */
    struct Counters {
      std::atomic< long long > elementCalls{ 0 };
      std::atomic< long long > assemblyPasses{ 0 }; // full 4-QP sweeps
      std::atomic< long long > materialCalls{ 0 };  // constitutive evaluations
      std::atomic< long long > newtonIterations{ 0 };
      std::atomic< long long > lineSearchTrials{ 0 };
      std::atomic< long long > maxNewtonIterations{ 0 };
      std::atomic< long long > nanosecondsTotal{ 0 };
      std::atomic< long long > nanosecondsMaterial{ 0 };
      std::atomic< long long > nanosecondsAssembly{ 0 }; // assembly excluding the material
      std::atomic< long long > nanosecondsBubbleSolve{ 0 };
      std::atomic< long long > nanosecondsCondensation{ 0 };

      /** MARMOT_MINI_PROFILE=<path>: the summary is rewritten to that file every
       *  `dumpEvery` element calls. A periodic dump rather than one at exit,
       *  because the library is dlopened from Python and its static
       *  destructors are not reliably run. */
      std::string path = [] {
        const char* v = std::getenv( "MARMOT_MINI_PROFILE" );
        return v ? std::string( v ) : std::string();
      }();
      bool                       enabled   = !path.empty();
      static constexpr long long dumpEvery = 50000;

      ~Counters() { dump(); }

      void dump()
      {
        if ( !enabled || elementCalls.load() == 0 ) {
          return;
        }
        std::FILE* out = std::fopen( path.c_str(), "w" );
        if ( !out ) {
          return;
        }
        const double    toSeconds = 1.0e-9;
        const long long calls     = elementCalls.load();
        std::fprintf( out,
                      "\n==== YIQUAD4_STABP_MINI element-local profile ====\n"
                      "  element calls                : %lld\n"
                      "  assembly passes (4 QP each)  : %lld   (%.2f per element call)\n"
                      "  constitutive evaluations     : %lld   (%.2f per element call)\n"
                      "  bubble Newton iterations     : %lld   (avg %.2f, max %lld per call)\n"
                      "  line-search trials           : %lld   (%.2f per call)\n"
                      "  ---- time ----\n"
                      "  total in computeKernels      : %8.2f s\n"
                      "    constitutive (material)    : %8.2f s  (%5.1f%%)\n"
                      "    assembly excl. material    : %8.2f s  (%5.1f%%)\n"
                      "    3x3 bubble solves          : %8.2f s  (%5.1f%%)\n"
                      "    static condensation        : %8.2f s  (%5.1f%%)\n"
                      "=================================================\n",
                      calls,
                      assemblyPasses.load(),
                      double( assemblyPasses.load() ) / calls,
                      materialCalls.load(),
                      double( materialCalls.load() ) / calls,
                      newtonIterations.load(),
                      double( newtonIterations.load() ) / calls,
                      maxNewtonIterations.load(),
                      lineSearchTrials.load(),
                      double( lineSearchTrials.load() ) / calls,
                      nanosecondsTotal.load() * toSeconds,
                      nanosecondsMaterial.load() * toSeconds,
                      100.0 * nanosecondsMaterial.load() / std::max( 1LL, nanosecondsTotal.load() ),
                      nanosecondsAssembly.load() * toSeconds,
                      100.0 * nanosecondsAssembly.load() / std::max( 1LL, nanosecondsTotal.load() ),
                      nanosecondsBubbleSolve.load() * toSeconds,
                      100.0 * nanosecondsBubbleSolve.load() / std::max( 1LL, nanosecondsTotal.load() ),
                      nanosecondsCondensation.load() * toSeconds,
                      100.0 * nanosecondsCondensation.load() / std::max( 1LL, nanosecondsTotal.load() ) );
        std::fclose( out );
      }
    };

    inline Counters& counters()
    {
      static Counters instance;
      return instance;
    }

    using Clock = std::chrono::steady_clock;

    inline long long since( const Clock::time_point& start )
    {
      return std::chrono::duration_cast< std::chrono::nanoseconds >( Clock::now() - start ).count();
    }

  } // namespace MiniProfiling

  class YStabPressureMiniInterfaceFiniteElement : public MarmotElement, public MarmotGeometryInterfaceElement< 3, 8 > {

  public:
    enum SectionType { Interface };

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

    using Material = MarmotStabPressureInterfaceMaterialHypoElastic;

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
      VectorDim            gradBubble; // surface gradient of b = (1-xi^2)(1-eta^2)
      double               bubble = 0.0;

      class QPStateVarManager : public MarmotStateVarVectorManager {
        inline const static auto layout = makeLayout( {
          { .name = "generalizedForce", .length = 3 },
          { .name = "surfaceStressPlus", .length = 9 },
          { .name = "surfaceStressMinus", .length = 9 },
          { .name = "volumetricResidual", .length = 1 },
          // MINI internal displacement bubble amplitude (3). Only qps[0]'s copy is
          // used -- it is an ELEMENT-level quantity, parked in the per-QP state
          // block because that is the only persistent storage an element gets.
          { .name = "bubbleAlpha", .length = 3 },
          { .name = "state block alignment padding", .length = ( 4 - ( ( 3 + 9 + 9 + 1 + 3 ) % 4 ) ) % 4 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< Eigen::Matrix< double, 3, 1 > > generalizedForce;
        Eigen::Map< Eigen::Matrix< double, 9, 1 > > surfaceStressPlus;
        Eigen::Map< Eigen::Matrix< double, 9, 1 > > surfaceStressMinus;
        Eigen::Map< Eigen::Matrix< double, 1, 1 > > volumetricResidual;
        Eigen::Map< Eigen::Matrix< double, 3, 1 > > bubbleAlpha;
        Eigen::Map< Eigen::VectorXd >               materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

        QPStateVarManager( double* v, int n )
          : MarmotStateVarVectorManager( v, layout ),
            generalizedForce( &find( "generalizedForce" ) ),
            surfaceStressPlus( &find( "surfaceStressPlus" ) ),
            surfaceStressMinus( &find( "surfaceStressMinus" ) ),
            volumetricResidual( &find( "volumetricResidual" ) ),
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
      }
    };

    std::vector< QuadraturePoint > qps;
    // Brezzi-Pitkaranta coefficient. Overridable, in order of precedence:
    //   1. env MARMOT_STABP_GAMMA   (diagnostic sweeps; mirrors the existing
    //      MARMOT_EXT_IFACE_DTAU_SCALE hook in this codebase)
    //   2. elementProperties[1]     (if the section supplies element properties)
    //   3. the default below
    // NOTE: the stabilisation is O(h_e^2) consistent, so gamma does not change the
    // converged solution -- but at a fixed mesh an excessive gamma biases the
    // pressure. Verify the deviatoric response (von Mises stress) is unchanged
    // whenever gamma is raised.
    double stabGamma = defaultStabGamma();

    /**
     * DIAGNOSTIC ABLATION SWITCH. MARMOT_MINI_BUBBLE=OFF disables ONLY the
     * internal displacement bubble: beta is held at zero, its residual/tangent
     * blocks are not condensed, and the element degenerates to the plain
     * equal-order Q1/Q1 displacement-pressure pairing with the SAME pressure
     * interpolation, the SAME pressure equation, the SAME Brezzi-Pitkaranta
     * term, the SAME quadrature and the SAME material. Nothing else changes.
     */
    static bool bubbleDisabled()
    {
      static const bool disabled = [] {
        const char* v = std::getenv( "MARMOT_MINI_BUBBLE" );
        return v && std::string( v ) == "OFF";
      }();
      return disabled;
    }

    /** Local pressure/displacement coupling B_u = d(R_p)/d(u), 4 x 24.
     *  Filled by the last computeKernels call; diagnostic only. */
    Eigen::Matrix< double, nDofP, nDofU > pressureDisplacementCoupling = Eigen::Matrix< double, nDofP, nDofU >::Zero();

    /** Local pressure/bubble coupling B_a = d(R_p)/d(beta), 4 x 3.
     *  Filled by the last computeKernels call; diagnostic only. */
    Eigen::Matrix< double, nDofP, nDim > pressureBubbleCoupling = Eigen::Matrix< double, nDofP, nDim >::Zero();

    // gamma = 0 is admissible and meaningful: it turns the Brezzi-Pitkaranta term
    // off entirely, leaving the MINI bubble as the only stabilising mechanism.
    // Parsed with strtod rather than atof so that a non-numeric value still falls
    // back to the default instead of silently becoming 0.
    static double defaultStabGamma()
    {
      if ( const char* e = std::getenv( "MARMOT_STABP_GAMMA" ) ) {
        char*        end = nullptr;
        const double v   = std::strtod( e, &end );
        if ( end != e && *end == '\0' && v >= 0.0 )
          return v;
      }
      return 0.2;
    }

    YStabPressureMiniInterfaceFiniteElement( int                                         elementID,
                                             FiniteElement::Quadrature::IntegrationTypes integrationType,
                                             SectionType                                 sectionType_ = Interface )
      : ParentGeometryElement(),
        elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
        elLabel( elementID ),
        sectionType( sectionType_ )
    {
      for ( const auto& qpi : FiniteElement::Quadrature::getGaussPointInfo( this->shape, integrationType ) )
        qps.emplace_back( qpi.xi, qpi.weight );
    }

    int getNumberOfRequiredStateVars() { return qps[0].getNumberOfRequiredStateVars() * (int)qps.size(); }

    std::vector< std::vector< std::string > > getNodeFields()
    {
      static std::vector< std::vector< std::string > > nf;
      if ( nf.empty() )
        for ( int i = 0; i < nNodes; i++ ) {
          nf.push_back( { "displacement" } );
          if ( i < nInterfaceNodes )
            nf[i].push_back( "interfacePressure" );
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
          for ( int c = 0; c < 3; c++ )
            perm[offU + 3 * A + c] = 4 * A + c;
          perm[offP + A] = 4 * A + 3;
          for ( int c = 0; c < 3; c++ )
            perm[offU + halfU + 3 * A + c] = 16 + 3 * A + c;
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
      for ( size_t i = 0; i < qps.size(); i++ )
        qps[i].assignStateVars( sv + i * nq, nq );
    }

    void assignProperty( const ElementProperties& p )
    {
      new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( p.elementProperties, p.nElementProperties );
      if ( p.nElementProperties > 1 && p.elementProperties[1] > 0.0 )
        stabGamma = p.elementProperties[1];
    }

    void assignProperty( const MarmotMaterialSection& s )
    {
      for ( auto& qp : qps )
        qp.material = std::make_unique< Material >( s.materialName,
                                                    s.materialProperties,
                                                    s.nMaterialProperties,
                                                    elLabel );
    }

    void assignMaterial( const std::string& name, const double* props, int nProps )
    {
      for ( auto& qp : qps )
        qp.material = std::make_unique< Material >( name, props, nProps, elLabel );
    }

    void assignNodeCoordinates( const double* c ) { ParentGeometryElement::assignNodeCoordinates( c ); }

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
        if ( qp.separationVector.norm() > tol && qp.separationVector.dot( qp.normal ) < 0.0 )
          qp.normal *= -1.0;
        qp.normalProjection  = qp.normal * qp.normal.transpose();
        qp.tangentProjection = TensorDim::Identity() - qp.normalProjection;
        if ( qp.separationVector.norm() > tol && qp.separationVector.dot( qp.normal ) <= tol )
          throw std::invalid_argument( "YIQUAD4_STABP_MINI: non-positive normal separation." );

        qp.J0xW  = qp.weight * qp.sqrtDetG * thickness;
        qp.hElem = std::sqrt( qp.sqrtDetG ); // in-plane mesh size

        // MINI bubble b = (1-xi^2)(1-eta^2): vanishes on the element boundary, so
        // the enrichment is strictly internal and condensable. Its SURFACE gradient
        // is obtained from db/dxi via the same chain rule the shape-function
        // gradients use: grad_s b = sum_A (db/dxi_r) [dXi/dx]_r , which we get by
        // expressing db/dxi in the nodal dN/dxi basis is not possible in general,
        // so use the metric inverse directly:  grad_s b = J * G^-1 * db/dxi.
        {
          const double xi = qp.xi( 0 ), eta = qp.xi( 1 );
          qp.bubble = ( 1.0 - xi * xi ) * ( 1.0 - eta * eta );
          Eigen::Matrix< double, 2, 1 > dbdxi;
          dbdxi( 0 )    = -2.0 * xi * ( 1.0 - eta * eta );
          dbdxi( 1 )    = -2.0 * eta * ( 1.0 - xi * xi );
          qp.gradBubble = qp.J * qp.G.inverse() * dbdxi; // 3x2 * 2x2 * 2x1
        }
        if ( qp.material )
          qp.material->setCharacteristicElementLength( qp.hElem );
      }
    }

    void setInitialConditions( StateTypes state, const double* )
    {
      if ( state != MarmotElement::MarmotMaterialInitialization )
        throw std::invalid_argument( "YIQUAD4_STABP_MINI: invalid initial condition." );
      for ( auto& qp : qps )
        qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                         qp.managedStateVars->materialStateVars.size() );
    }

    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
    {
      auto&      profile      = MiniProfiling::counters();
      const bool profiling    = profile.enabled;
      const auto elementStart = profiling ? MiniProfiling::Clock::now() : MiniProfiling::Clock::time_point{};
      if ( profiling ) {
        const long long call = profile.elementCalls.fetch_add( 1, std::memory_order_relaxed ) + 1;
        if ( call % MiniProfiling::Counters::dumpEvery == 0 ) {
          profile.dump();
        }
      }

      Eigen::Map< const RhsSized > QTotal( QTotal_ ), dQ( dQ_ );
      Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
      Eigen::Map< RhsSized >       Pe( Pe_ );

      using V3  = Eigen::Matrix< double, 3, 1 >;
      using M33 = Eigen::Matrix< double, 3, 3 >;

      // Snapshot the incoming (committed) material state of every QP: the inner
      // bubble Newton re-evaluates the material several times, and computeStress
      // COMMITS, so each trial must start from the same state. Only the final
      // pass commits.
      std::vector< std::vector< double > > matIn( qps.size() );
      for ( size_t i = 0; i < qps.size(); i++ ) {
        auto& m  = qps[i].managedStateVars->materialStateVars;
        matIn[i] = std::vector< double >( m.data(), m.data() + m.size() );
      }

      struct Acc {
        RhsSized                              Pe;
        KeSizedMatrix                         Ke;
        V3                                    Rb;
        M33                                   Kbb;
        Eigen::Matrix< double, 3, nDofU >     Kbu;
        Eigen::Matrix< double, nDofU, 3 >     Kub;
        Eigen::Matrix< double, 3, nDofP >     Kbp;
        Eigen::Matrix< double, nDofP, 3 >     Kpb;
        Eigen::Matrix< double, nDofP, nDofU > Kpu; // diagnostic: B_u

        // Per-QP results carried out of the pass, so that the element state can
        // be written from ANY pass without recomputing it. This is what lets the
        // redundant final "commit" sweep be dropped: it used to re-run all four
        // constitutive updates at a beta that had just been evaluated, and the
        // profile showed it was 1 of every 3.11 assembly passes, i.e. ~32% of all
        // element-local work.
        std::array< Eigen::Matrix< double, 3, 1 >, 8 > qpForce;
        std::array< Eigen::Matrix< double, 9, 1 >, 8 > qpSurfacePlus;
        std::array< Eigen::Matrix< double, 9, 1 >, 8 > qpSurfaceMinus;
        std::array< double, 8 >                        qpVolumetricResidual;
      };

      auto assembleAll = [&]( const V3& beta ) {
        const auto passStart = profiling ? MiniProfiling::Clock::now() : MiniProfiling::Clock::time_point{};
        long long  passMaterialNanoseconds = 0;
        if ( profiling ) {
          profile.assemblyPasses.fetch_add( 1, std::memory_order_relaxed );
        }
        Acc A;
        A.Pe.setZero();
        A.Ke.setZero();
        A.Rb.setZero();
        A.Kbb.setZero();
        A.Kbu.setZero();
        A.Kub.setZero();
        A.Kbp.setZero();
        A.Kpb.setZero();
        A.Kpu.setZero();

        for ( size_t iq = 0; iq < qps.size(); iq++ ) {
          auto& qp  = qps[iq];
          auto& msv = qp.managedStateVars->materialStateVars;
          for ( int k = 0; k < msv.size(); k++ )
            msv[k] = matIn[iq][k];

          const auto& W  = qp.NmatJump;
          const auto& Bs = qp.BmatSide;
          const auto& Np = qp.N;
          const auto& Gp = qp.gradN;

          BAvgSurfaceSized BP = BAvgSurfaceSized::Zero(), BM = BAvgSurfaceSized::Zero();
          BM.block< nTensor, halfU >( 0, 0 )     = Bs;
          BP.block< nTensor, halfU >( 0, halfU ) = Bs;

          const auto dU = dQ.segment< nDofU >( offU );
          const auto dP = dQ.segment< nDofP >( offP );
          const auto Pt = QTotal.segment< nDofP >( offP );

          Eigen::Matrix< double, 6, 1 > dU6;
          dU6.segment< 3 >( 0 ) = qp.NmatSide * dU.segment< halfU >( halfU );
          dU6.segment< 3 >( 3 ) = qp.NmatSide * dU.segment< halfU >( 0 );

          // MINI: ubar += b*beta  =>  Abar_ij += beta_i (grad_s b)_j, added equally
          // to both stations because Abar = (A+ + A-)/2.
          Eigen::Matrix< double, 9, 3 > Bbub = Eigen::Matrix< double, 9, 3 >::Zero();
          for ( int i = 0; i < 3; i++ )
            for ( int j = 0; j < 3; j++ )
              Bbub( 3 * i + j, i ) = qp.gradBubble( j );

          Eigen::Matrix< double, 18, 1 > dSurf;
          dSurf.segment< 9 >( 0 ) = Bs * dU.segment< halfU >( halfU ) + Bbub * beta;
          dSurf.segment< 9 >( 9 ) = Bs * dU.segment< halfU >( 0 ) + Bbub * beta;

          Eigen::Matrix< double, 3, 1 > f  = qp.managedStateVars->generalizedForce;
          Eigen::Matrix< double, 9, 1 > Sp = qp.managedStateVars->surfaceStressPlus;
          Eigen::Matrix< double, 9, 1 > Sm = qp.managedStateVars->surfaceStressMinus;
          double                        rp = 0.0;

          Eigen::Matrix< double, 3, 3, Eigen::RowMajor > Qww;
          Eigen::Matrix< double, 3, 9, Eigen::RowMajor > QwAp, QwAm;
          Eigen::Matrix< double, 9, 3, Eigen::RowMajor > QApw, QAmw;
          Eigen::Matrix< double, 9, 9, Eigen::RowMajor > QApAp, QApAm, QAmAp, QAmAm;
          Eigen::Matrix< double, 3, 1 >                  Qwp;
          Eigen::Matrix< double, 9, 1 >                  QApp, QAmp;
          Eigen::Matrix< double, 1, 3 >                  Qpw;
          Eigen::Matrix< double, 1, 9 >                  QpAp, QpAm;
          double                                         Qpp = 0.0;

          Material::Response      resp{ f.data(), Sp.data(), Sm.data(), &rp };
          Material::Tangents      tg{ Qww.data(),
                                 QwAp.data(),
                                 QwAm.data(),
                                 QApw.data(),
                                 QApAp.data(),
                                 QApAm.data(),
                                 QAmw.data(),
                                 QAmAp.data(),
                                 QAmAm.data(),
                                 Qwp.data(),
                                 QApp.data(),
                                 QAmp.data(),
                                 Qpw.data(),
                                 QpAp.data(),
                                 QpAm.data(),
                                 &Qpp };
          Material::Deformation   def{ dU6.data(),
                                     dSurf.data(),
                                     qp.normal.data(),
                                     qp.separationVector.data(),
                                     ( Np * dP )( 0, 0 ) };
          Material::TimeIncrement tinc{ time, dT };
          if ( profiling ) {
            const auto materialStart = MiniProfiling::Clock::now();
            qp.material->computeStress( msv.data(), resp, tg, def, tinc );
            passMaterialNanoseconds += MiniProfiling::since( materialStart );
            profile.materialCalls.fetch_add( 1, std::memory_order_relaxed );
          }
          else {
            qp.material->computeStress( msv.data(), resp, tg, def, tinc );
          }

          A.qpForce[iq]              = f;
          A.qpSurfacePlus[iq]        = Sp;
          A.qpSurfaceMinus[iq]       = Sm;
          A.qpVolumetricResidual[iq] = rp;

          const double h   = qp.material->getInterfaceThickness();
          const double mu  = qp.material->getShearModulus();
          const double tau = stabGamma * qp.hElem * qp.hElem / ( 2.0 * mu );

          const Eigen::Matrix< double, nDofU, 1 > Ru = W.transpose() * f + BP.transpose() * Sp + BM.transpose() * Sm;
          const Eigen::Matrix< double, nDofP, 1 > Rp = h *
                                                       ( Np.transpose() * rp + tau * ( Gp.transpose() * ( Gp * Pt ) ) );
          A.Pe.segment< nDofU >( offU ) -= Ru * qp.J0xW;
          A.Pe.segment< nDofP >( offP ) -= Rp * qp.J0xW;
          A.Rb += ( Bbub.transpose() * ( Sp + Sm ) ) * qp.J0xW;

          const Eigen::Matrix< double, nDofU, nDofU > Kuu = W.transpose() * Qww * W + W.transpose() * QwAp * BP +
                                                            W.transpose() * QwAm * BM + BP.transpose() * QApw * W +
                                                            BP.transpose() * QApAp * BP + BP.transpose() * QApAm * BM +
                                                            BM.transpose() * QAmw * W + BM.transpose() * QAmAp * BP +
                                                            BM.transpose() * QAmAm * BM;
          const Eigen::Matrix< double, nDofU, nDofP > Kup = ( W.transpose() * Qwp + BP.transpose() * QApp +
                                                              BM.transpose() * QAmp ) *
                                                            Np;
          const Eigen::Matrix< double, nDofP, nDofU > Kpu = h *
                                                            ( Np.transpose() * ( Qpw * W + QpAp * BP + QpAm * BM ) );
          const Eigen::Matrix< double, nDofP, nDofP > Kpp = h * ( Np.transpose() * Qpp * Np +
                                                                  tau * ( Gp.transpose() * Gp ) );

          A.Ke.block< nDofU, nDofU >( offU, offU ) += Kuu * qp.J0xW;
          A.Ke.block< nDofU, nDofP >( offU, offP ) += Kup * qp.J0xW;
          A.Ke.block< nDofP, nDofU >( offP, offU ) += Kpu * qp.J0xW;
          A.Kpu += Kpu * qp.J0xW; // diagnostic copy of the raw (uncondensed) B_u
          A.Ke.block< nDofP, nDofP >( offP, offP ) += Kpp * qp.J0xW;

          const Eigen::Matrix< double, 9, 3 > dSpdb = ( QApAp + QApAm ) * Bbub;
          const Eigen::Matrix< double, 9, 3 > dSmdb = ( QAmAp + QAmAm ) * Bbub;
          Eigen::Matrix< double, 3, 3 >       dfdb;
          for ( int c = 0; c < 3; c++ )
            dfdb.col( c ) = ( QwAp + QwAm ) * Bbub.col( c );

          A.Kbb += ( Bbub.transpose() * ( dSpdb + dSmdb ) ) * qp.J0xW;
          A.Kbu += ( Bbub.transpose() * ( ( QApw + QAmw ) * W + ( QApAp + QAmAp ) * BP + ( QApAm + QAmAm ) * BM ) ) *
                   qp.J0xW;
          A.Kub += ( W.transpose() * dfdb + BP.transpose() * dSpdb + BM.transpose() * dSmdb ) * qp.J0xW;
          A.Kbp += ( Bbub.transpose() * ( QApp + QAmp ) * Np ) * qp.J0xW;
          A.Kpb += ( h * ( Np.transpose() * ( QpAp + QpAm ) * Bbub ) ) * qp.J0xW;
        }
        if ( profiling ) {
          const long long total = MiniProfiling::since( passStart );
          profile.nanosecondsMaterial.fetch_add( passMaterialNanoseconds, std::memory_order_relaxed );
          profile.nanosecondsAssembly.fetch_add( total - passMaterialNanoseconds, std::memory_order_relaxed );
        }
        return A;
      };

      // The material state left in the state vector after a pass belongs to that
      // pass's beta. Snapshot it whenever a pass is ACCEPTED, so a later rejected
      // trial can be undone without re-running the constitutive updates.
      std::vector< std::vector< double > > acceptedMaterialState( qps.size() );

      auto snapshotMaterialState = [&]() {
        for ( size_t i = 0; i < qps.size(); i++ ) {
          auto& m = qps[i].managedStateVars->materialStateVars;
          acceptedMaterialState[i].assign( m.data(), m.data() + m.size() );
        }
      };

      auto commitAccepted = [&]( const Acc& A ) {
        for ( size_t i = 0; i < qps.size(); i++ ) {
          auto& managed = *qps[i].managedStateVars;
          std::copy( acceptedMaterialState[i].begin(),
                     acceptedMaterialState[i].end(),
                     managed.materialStateVars.data() );
          managed.generalizedForce        = A.qpForce[i];
          managed.surfaceStressPlus       = A.qpSurfacePlus[i];
          managed.surfaceStressMinus      = A.qpSurfaceMinus[i];
          managed.volumetricResidual( 0 ) = A.qpVolumetricResidual[i];
        }
      };

      // ---- inner Newton on the internal bubble amplitude: solve R_beta(beta)=0 ----
      // Necessary because beta enters the strain and therefore the plastic
      // return-mapping, so R_beta is NONLINEAR in beta. A single step from beta=0
      // leaves Pe (evaluated at beta=0) inconsistent with the condensed Ke, which
      // shows up directly as a large tangent-vs-FD error on the displacement columns.
      // Backtracking line search + exception guard, mirroring the line search the
      // material's own local gamma-solve already uses. Without it an undamped beta
      // step overshoots at yield onset, drives the strain into a state where the
      // local traction-equilibrium solve cannot converge, and the StressUpdateFailed
      // propagates out as a global cutback (measured: 1 cutback, 1 StressUpdateFailed
      // for MINI vs 0/0 for the same element without the bubble).
      if ( bubbleDisabled() ) {
        // Ablation: no bubble at all. beta stays zero, so it contributes
        // nothing to the strain, and its blocks are NOT condensed -- the
        // element is the plain Q1/Q1 pairing with the identical pressure
        // equation and stabilisation.
        Acc A0 = assembleAll( V3::Zero() );
        snapshotMaterialState();
        commitAccepted( A0 );
        Pe                                   = A0.Pe;
        Ke                                   = A0.Ke;
        qps[0].managedStateVars->bubbleAlpha = V3::Zero();
        pressureDisplacementCoupling         = A0.Kpu;
        pressureBubbleCoupling               = A0.Kpb;
        return;
      }

      V3  beta = qps[0].managedStateVars->bubbleAlpha;
      Acc A    = assembleAll( beta );
      snapshotMaterialState();
      int newtonIterations = 0;
      for ( int it = 0; it < 20; it++ ) {
        const double rn  = A.Rb.norm();
        const double tol = 1.0e-10 * std::max( 1.0, A.Pe.norm() );
        if ( rn <= tol )
          break;
        ++newtonIterations;
        const auto solveStart = profiling ? MiniProfiling::Clock::now() : MiniProfiling::Clock::time_point{};
        const V3   dBeta      = -A.Kbb.fullPivLu().solve( A.Rb );
        if ( profiling ) {
          profile.nanosecondsBubbleSolve.fetch_add( MiniProfiling::since( solveStart ), std::memory_order_relaxed );
        }

        // Backtracking line search on the inner bubble Newton: accept only a step
        // that does not increase ||R_beta||, and halve otherwise. The catch treats
        // a throwing trial (a beta that drove the local traction-equilibrium solve
        // past convergence) exactly like a non-improving step.
        //
        // This is REQUIRED FOR ROBUSTNESS, not for speed. Measured both ways:
        //   * accuracy: irrelevant -- the oscillation result is bit-identical
        //     (p2p 0.03897), i.e. beta reaches the same root either way;
        //   * iterations: no benefit (180 vs 173 on the diagnosis benchmark);
        //   * runtime: it COSTS ~67% there (1953 s vs 1170 s), because it spends
        //     extra material evaluations on trial steps that are valid but merely
        //     non-improving;
        //   * robustness: decisive. Without it, beta can take a full Newton step
        //     that increases the residual, wander into a bad state, and make the
        //     material throw. On the parametric-study mesh (stiff, angle 10,
        //     h=0.01, fy=5) that produced "Element ... requests for a cutback",
        //     minInc exhaustion and a FAILED simulation at t~0.185, on a case the
        //     older EIQUAD4 completes. With the line search the analysis runs.
        // The runtime penalty is worth paying: these jobs are ~50 s, and not
        // completing at all is not a trade.
        double alpha    = 1.0;
        bool   accepted = false;
        for ( int ls = 0; ls < 10; ls++ ) {
          if ( profiling ) {
            profile.lineSearchTrials.fetch_add( 1, std::memory_order_relaxed );
          }
          try {
            Acc trial = assembleAll( beta + alpha * dBeta );
            if ( trial.Rb.norm() <= tol || trial.Rb.norm() < rn ) {
              beta += alpha * dBeta;
              A = trial;
              snapshotMaterialState();
              accepted = true;
              break;
            }
          }
          catch ( const std::exception& ) {
            // the trial beta drove the local material solve past convergence:
            // treat exactly like a non-improving step and halve
          }
          alpha *= 0.5;
        }
        if ( !accepted )
          break; // keep the last good beta; the outer Newton continues from there
      }
      if ( profiling ) {
        profile.newtonIterations.fetch_add( newtonIterations, std::memory_order_relaxed );
        long long previous = profile.maxNewtonIterations.load( std::memory_order_relaxed );
        while ( newtonIterations > previous &&
                !profile.maxNewtonIterations.compare_exchange_weak( previous, newtonIterations ) ) {
        }
      }

      // The accepted pass already evaluated everything at this beta, so the
      // element state is written from it instead of re-running a whole
      // constitutive sweep. The material state is restored from the snapshot,
      // which undoes any rejected trial that ran afterwards.
      //
      // MARMOT_MINI_LEGACY_COMMIT_PASS=1 restores the old behaviour (a full
      // extra sweep at the same beta), kept ONLY so the two can be timed
      // back-to-back under identical machine load.
      static const bool legacyCommitPass = [] {
        const char* v = std::getenv( "MARMOT_MINI_LEGACY_COMMIT_PASS" );
        return v && std::string( v ) == "1";
      }();
      if ( legacyCommitPass ) {
        A = assembleAll( beta );
      }
      commitAccepted( A );

      const auto condensationStart = profiling ? MiniProfiling::Clock::now() : MiniProfiling::Clock::time_point{};

      // One factorisation, three solves -- no explicit K_bb^-1.
      // Spectrum of the bubble block, for comparison with the two-face element's Qbar.
      // grad_s b is a surface gradient, so the contraction sits on tangential slots and
      // lambda_min is expected O(G) with no H-scaling.
      if ( std::getenv( "MARMOT_MINI_KBB_SPECTRUM" ) ) {
        static std::atomic< long > kbbCalls{ 0 };
        const long                 kbbN = ++kbbCalls;
        if ( kbbN % 20000 == 0 ) {
          Eigen::JacobiSVD< Eigen::Matrix3d > kbbSvd( Eigen::Matrix3d( A.Kbb ) );
          const auto                          sv = kbbSvd.singularValues();
          std::printf( "  [kbb] n=%ld  sMax=%.6e sMin=%.6e  sMin/sMax=%.6e\n",
                       kbbN,
                       sv( 0 ),
                       sv( 2 ),
                       sv( 0 ) > 0.0 ? sv( 2 ) / sv( 0 ) : 0.0 );
          std::fflush( stdout );
        }
      }

      const auto                              KbbFactorisation = A.Kbb.fullPivLu();
      const V3                                KbbInvRb         = KbbFactorisation.solve( A.Rb );
      const Eigen::Matrix< double, 3, nDofU > KbbInvKbu        = KbbFactorisation.solve( A.Kbu );
      const Eigen::Matrix< double, 3, nDofP > KbbInvKbp        = KbbFactorisation.solve( A.Kbp );

      Pe = A.Pe;
      Ke = A.Ke;
      // R_beta is ~0 now, so the condensed residual equals Pe; only the tangent
      // needs the Schur complement.
      Pe.segment< nDofU >( offU ) += A.Kub * KbbInvRb;
      Pe.segment< nDofP >( offP ) += A.Kpb * KbbInvRb;
      Ke.block< nDofU, nDofU >( offU, offU ) -= A.Kub * KbbInvKbu;
      Ke.block< nDofU, nDofP >( offU, offP ) -= A.Kub * KbbInvKbp;
      Ke.block< nDofP, nDofU >( offP, offU ) -= A.Kpb * KbbInvKbu;
      Ke.block< nDofP, nDofP >( offP, offP ) -= A.Kpb * KbbInvKbp;

      qps[0].managedStateVars->bubbleAlpha = beta;

      pressureDisplacementCoupling = A.Kpu;
      pressureBubbleCoupling       = A.Kpb;

      if ( profiling ) {
        profile.nanosecondsCondensation.fetch_add( MiniProfiling::since( condensationStart ),
                                                   std::memory_order_relaxed );
        profile.nanosecondsTotal.fetch_add( MiniProfiling::since( elementStart ), std::memory_order_relaxed );
      }
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
      throw std::invalid_argument( "YIQUAD4_STABP_MINI: distributed loads not implemented." );
    }
    void computeBodyForce( double*, double*, const double*, const double*, double, double )
    {
      throw std::invalid_argument( "YIQUAD4_STABP_MINI: body forces not implemented." );
    }
    void computeConsistentInertia( double* ) { throw std::runtime_error( "YIQUAD4_STABP_MINI: no inertia." ); }
    void computeLumpedInertia( double* ) { throw std::runtime_error( "YIQUAD4_STABP_MINI: no inertia." ); }

    StateView getStateView( const std::string& name, int qpNumber )
    {
      const auto& qp = qps[qpNumber];
      if ( qp.managedStateVars->contains( name ) )
        return qp.managedStateVars->getStateView( name );
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
