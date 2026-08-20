/* ---------------------------------------------------------------------
 *  Marmot - WarpingInterfaceFiniteElement  (WIQUAD4 / WILINE2)
 *  Alexandros Stathas alexandros.stathas@boku.ac.at
 * --------------------------------------------------------------------- */

/**
 * @file WarpingInterfaceFiniteElement.h
 * @brief Through-thickness-resolved interface element carrying SYMMETRIC
 * (parabolic) and ANTISYMMETRIC (cubic) WARPING micro-fields, statically
 * condensed (WIQUAD4 / WILINE2).
 *
 * WHAT IT IS
 * ----------
 * Same global unknowns as GLIQUAD4/YIQUAD4 -- nNodes x nDim displacement DOF,
 * nothing else -- so it is a drop-in replacement in an existing mesh: only
 * `type=` changes in the input file. Internally it adds two warping amplitude
 * fields on the midsurface,
 *
 *     u(x_s, zeta) += phi_s(zeta) w_s(x_s) + phi_a(zeta) w_a(x_s),
 *     phi_s = 1 - zeta^2   (even  -> symmetric,     parabolic),
 *     phi_a = zeta - zeta^3 (odd  -> antisymmetric, cubic),
 *
 * which MarmotWarpingInterfaceMaterialHypoElastic turns into two extra
 * generalized strains Ws = grad_s w_s and Wa = grad_s w_a. Both profiles
 * vanish at zeta = +/-1, so the warping never touches the face displacements
 * and therefore never touches the coupling to the surrounding bulk elements.
 * That is exactly what makes the amplitudes element-local and condensable:
 * the assembled global system is bit-for-bit the same size and sparsity as
 * the plain element's.
 *
 * INTERNAL BASIS -- THE MODES MUST BE BUBBLES
 * -------------------------------------------
 * Only grad_s w enters the material (see the material header for the proof
 * that the warping's normal-gradient content is already spanned by the
 * Gauss-Lobatto station gradients g^(alpha) and must NOT be fed in twice).
 * The amplitude fields are expanded in midsurface BUBBLE modes -- functions
 * vanishing on the whole element boundary:
 *
 *     quad midsurface (WIQUAD4): M = { b, b*xi, b*eta },  b = (1-xi^2)(1-eta^2),
 *     line midsurface (WILINE2): M = { b },               b = (1-xi^2),
 *
 * giving nWarpDof = 2 * nWarpModes * nDim internal unknowns per element
 * (18 for WIQUAD4, 4 for WILINE2), element-local and condensable.
 *
 * Boundary-vanishing is the CONSISTENCY condition, not a convenience: an
 * enhanced-strain enrichment must satisfy int_Ae grad_s M dA = 0, otherwise a
 * uniform stress state produces a nonzero internal residual and the element
 * relaxes a state that is already exact. It also kills the in-plane-rotation
 * zero mode, which a general (non-bubble) amplitude field would carry.
 *
 * THE REMAINING, EXACTLY KNOWN NULL SPACE
 * ---------------------------------------
 * One family of directions is still exactly inert, independently of the
 * material: the components (grad_s w_a)_{n j} of the ANTISYMMETRIC field --
 * its normal component varying in-plane. Setting
 * g^(alpha) -= phi_a(xi_alpha) (grad_s w_a)_{n j} cancels the symmetrised
 * strain at EVERY station while preserving the weighted-mean constraint,
 * because sum_alpha lambda_alpha phi_a(xi_alpha) = 0 exactly. Measured on the
 * 6x6 tangential sub-blocks: A+ rank 5, Ws rank 5, Wa rank 3, with the two
 * relaxed Wa directions at 3.6e-15 relative to 1.1e2.
 *
 * Those directions carry exactly zero energy AND exactly zero residual, and the
 * displacement/warping coupling block annihilates them too (the generalized
 * stresses do not change along them, so K_db v = 0 for every null vector v).
 * The Schur complement K_db K_bb^+ K_bd is therefore EXACT with the
 * Moore-Penrose pseudo-inverse -- this is not a regularisation trading accuracy
 * for robustness, it is the correct operator. A rank-revealing complete-
 * orthogonal decomposition is used, which additionally keeps the element
 * well-posed if perfect plasticity drives further directions soft.
 *
 * COST
 * ----
 * The global system is unchanged. The element evaluation is not: each inner
 * Newton iteration re-runs all NStations constitutive updates per quadrature
 * point. This buys through-thickness-nonlinear IN-PLANE strain, which neither
 * the two-station (X/Y) nor the Gauss-Lobatto element can represent -- the
 * latter already has arbitrary transverse-shear warping through its free
 * g^(alpha), but forces the in-plane strain to stay linear across the layer.
 */
#pragma once

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotStateVarVectorManager.h"
#include "Marmot/MarmotWarpingInterfaceMaterialHypoElastic.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Marmot::Elements {

  template < int nDim, int nNodes, int NStations = 5 >
  class WarpingInterfaceFiniteElement : public MarmotElement, public MarmotGeometryInterfaceElement< nDim, nNodes > {

  public:
    enum SectionType { Interface };

    static constexpr int nInterfaceNodes = nNodes / 2;
    static constexpr int nTensor         = nDim * nDim;
    static constexpr int nDofU           = nNodes * nDim;
    static constexpr int sizeLoadVector  = nDofU; // global DOF: displacement only
    static constexpr int nCoordinates    = nNodes * nDim;
    static constexpr int halfU           = nDofU / 2;
    static constexpr int nParam          = nDim - 1;

    /** Non-constant midsurface bilinear modes: {xi, eta, xi*eta} or {xi}. */
    static constexpr int nWarpModes = nInterfaceNodes - 1;
    static constexpr int nWarpBlock = nWarpModes * nDim; // one amplitude field
    static constexpr int nWarpDof   = 2 * nWarpBlock;    // symmetric + antisymmetric

    /** Element-side generalized pair: ( w, A+, A-, Ws, Wa ). */
    static constexpr int nGenElem = nDim + 4 * nTensor;

    using Material = MarmotWarpingInterfaceMaterialHypoElasticN< NStations >;

    static constexpr int nGenMat = Material::nGeneralized; // always 39

    using ParentGeometryElement = MarmotGeometryInterfaceElement< nDim, nNodes >;

    using XiSized              = typename ParentGeometryElement::XiSized;
    using NSized               = typename ParentGeometryElement::NSized;
    using dNdXiSized           = typename ParentGeometryElement::dNdXiSized;
    using SurfaceJacobianSized = typename ParentGeometryElement::SurfaceJacobianSized;
    using MetricSized          = typename ParentGeometryElement::MetricSized;
    using GradSized            = typename ParentGeometryElement::GradSized;
    using VectorDim            = typename ParentGeometryElement::VectorDim;
    using TensorDim            = typename ParentGeometryElement::TensorDim;
    using NMatrixSized         = typename ParentGeometryElement::NMatrixSized;
    using NJumpMatrixSized     = typename ParentGeometryElement::NJumpMatrixSized;
    using BSurfaceSized        = typename ParentGeometryElement::BSurfaceSized;
    using BAvgSurfaceSized     = typename ParentGeometryElement::BAvgSurfaceSized;

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;

    using BWarpSized = Eigen::Matrix< double, nTensor, nWarpBlock >;
    using WarpVector = Eigen::Matrix< double, nWarpDof, 1 >;
    using WarpMatrix = Eigen::Matrix< double, nWarpDof, nWarpDof >;

    /** Generalized B operator over ( displacement , warping amplitudes ). */
    using BGenSized = Eigen::Matrix< double, nGenElem, nDofU + nWarpDof >;

    Eigen::Map< const Eigen::VectorXd > elementProperties;
    const int                           elLabel;
    const SectionType                   sectionType;

    /**
     * Index of each element-side generalized component inside the material's
     * fixed 39-component (3D) layout. The identity for nDim == 3.
     */
    static const std::array< int, nGenElem >& materialIndexMap()
    {
      static const std::array< int, nGenElem > map = [] {
        std::array< int, nGenElem > m{};
        int                         k = 0;
        for ( int i = 0; i < nDim; ++i ) {
          m[k++] = Material::offsetJump + i;
        }
        const int slotOffsets[4] = { Material::offsetSurfacePlus,
                                     Material::offsetSurfaceMinus,
                                     Material::offsetWarpingSymmetric,
                                     Material::offsetWarpingAntisymmetric };
        for ( const int slot : slotOffsets ) {
          for ( int i = 0; i < nDim; ++i ) {
            for ( int j = 0; j < nDim; ++j ) {
              m[k++] = slot + 3 * i + j;
            }
          }
        }
        return m;
      }();
      return map;
    }

    struct QuadraturePoint {
      const XiSized xi;
      const double  weight;

      double sqrtDetG = 0.0;
      double J0xW     = 0.0;
      double hElem    = 0.0;

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

      /** Surface gradients of the non-constant midsurface modes, nDim x nWarpModes. */
      Eigen::Matrix< double, nDim, nWarpModes > gradWarpModes;

      class QPStateVarManager : public MarmotStateVarVectorManager {
        static constexpr int nRaw = nDim           // generalizedForce
                                    + nTensor      // surfaceStressPlus
                                    + nTensor      // surfaceStressMinus
                                    + nTensor      // warpingStressSymmetric
                                    + nTensor      // warpingStressAntisymmetric
                                    + nWarpDof     // warpingAmplitudes (element-level, qps[0] only)
                                    + 2 * nDim     // displacement
                                    + 2 * nTensor  // surfaceStrain
                                    + 2 * nTensor; // warpingStrain

        inline const static auto layout = makeLayout( {
          { .name = "generalizedForce", .length = nDim },
          { .name = "surfaceStressPlus", .length = nTensor },
          { .name = "surfaceStressMinus", .length = nTensor },
          { .name = "warpingStressSymmetric", .length = nTensor },
          { .name = "warpingStressAntisymmetric", .length = nTensor },
          { .name = "warpingAmplitudes", .length = nWarpDof },
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
        Eigen::Map< Eigen::Matrix< double, nWarpDof, 1 > >    warpingAmplitudes;
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
            warpingAmplitudes( &find( "warpingAmplitudes" ) ),
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
        gradWarpModes.setZero();
      }
    };

    std::vector< QuadraturePoint > qps;

    WarpingInterfaceFiniteElement( int                                         elementID,
                                   FiniteElement::Quadrature::IntegrationTypes integrationType,
                                   SectionType                                 sectionType_ = Interface )
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
        }
      }
      return nf;
    }

    std::vector< int > getDofIndicesPermutationPattern()
    {
      static std::vector< int > perm;
      if ( perm.empty() ) {
        perm.resize( sizeLoadVector );
        for ( int i = 0; i < sizeLoadVector; i++ ) {
          perm[i] = i;
        }
      }
      return perm;
    }

    int getNNodes() { return nNodes; }
    int getNSpatialDimensions() { return nDim; }
    int getNDofPerElement() { return sizeLoadVector; }

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

    /**
     * Parametric derivatives of the midsurface warping modes.
     *
     * Every mode is a BUBBLE -- it vanishes on the whole element boundary.
     * That is not a stylistic choice, it is the consistency (patch-test)
     * condition of an enhanced-strain enrichment:
     *
     *     int_Ae grad_s M dA = closed_int_dAe M n ds = 0     for every mode M.
     *
     * If it is violated, a UNIFORM stress state produces a nonzero internal
     * residual -- the enrichment then relaxes a state that is already the exact
     * solution, the element fails the patch test and stops reducing to
     * GLIQUAD4. (Measured with the non-bubble basis {xi, eta, xi*eta}: 21%
     * error in the residual under uniform through-thickness loading.)
     *
     * Vanishing on the element boundary also removes the in-plane-rotation zero
     * mode: grad_s w cannot be pointwise skew for a nonzero combination of
     * bubbles, since the mode gradients are linearly independent as functions.
     */
    static Eigen::Matrix< double, nParam, nWarpModes > warpModeParametricGradients( const XiSized& xi )
    {
      Eigen::Matrix< double, nParam, nWarpModes > d = Eigen::Matrix< double, nParam, nWarpModes >::Zero();

      if constexpr ( nParam == 2 ) {
        const double x = xi( 0 ), e = xi( 1 );
        const double bx = 1.0 - x * x, be = 1.0 - e * e;

        // M0 = b = (1-xi^2)(1-eta^2)
        d( 0, 0 ) = -2.0 * x * be;
        d( 1, 0 ) = -2.0 * e * bx;

        // M1 = b * xi
        d( 0, 1 ) = be * ( 1.0 - 3.0 * x * x );
        d( 1, 1 ) = -2.0 * x * e * bx;

        // M2 = b * eta
        d( 0, 2 ) = -2.0 * x * e * be;
        d( 1, 2 ) = bx * ( 1.0 - 3.0 * e * e );
      }
      else {
        // M0 = b = 1 - xi^2
        d( 0, 0 ) = -2.0 * xi( 0 );
      }
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

        const VectorDim xB  = qp.NmatSide * this->getSideCoordinates( 0 );
        const VectorDim xT  = qp.NmatSide * this->getSideCoordinates( 1 );
        qp.separationVector = xT - xB;

        constexpr double tol = 1.0e-12;
        if ( qp.separationVector.norm() > tol && qp.separationVector.dot( qp.normal ) < 0.0 ) {
          qp.normal *= -1.0;
        }
        qp.normalProjection  = qp.normal * qp.normal.transpose();
        qp.tangentProjection = TensorDim::Identity() - qp.normalProjection;
        if ( qp.separationVector.norm() > tol && qp.separationVector.dot( qp.normal ) <= tol ) {
          throw std::invalid_argument(
            "WarpingInterfaceFiniteElement: paired faces have no positive normal separation." );
        }

        qp.J0xW  = qp.weight * qp.sqrtDetG * thickness;
        qp.hElem = ( nDim == 3 ) ? std::sqrt( qp.sqrtDetG ) : qp.sqrtDetG;

        // Surface gradients of the warping modes, by the same chain rule the
        // MINI bubble uses: grad_s M = J G^-1 dM/dxi.
        qp.gradWarpModes = qp.J * qp.G.inverse() * warpModeParametricGradients( qp.xi );

        if ( qp.material ) {
          qp.material->setCharacteristicElementLength( qp.hElem );
        }
      }
    }

    void setInitialConditions( StateTypes state, const double* )
    {
      if ( state != MarmotElement::MarmotMaterialInitialization ) {
        throw std::invalid_argument( "WarpingInterfaceFiniteElement: invalid initial condition." );
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

    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
    {
      (void)QTotal_;

      Eigen::Map< const RhsSized > dQ( dQ_ );
      Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
      Eigen::Map< RhsSized >       Pe( Pe_ );

      // The inner warping Newton re-evaluates the material several times and
      // computeStress COMMITS, so every trial must restart from the same
      // incoming state. Only the final pass is allowed to commit.
      std::vector< std::vector< double > > materialStateIn( qps.size() );
      for ( size_t i = 0; i < qps.size(); i++ ) {
        auto& m            = qps[i].managedStateVars->materialStateVars;
        materialStateIn[i] = std::vector< double >( m.data(), m.data() + m.size() );
      }

      struct Accumulated {
        RhsSized                                 Pe;
        KeSizedMatrix                            Kdd;
        WarpVector                               Rb;
        WarpMatrix                               Kbb;
        Eigen::Matrix< double, nWarpDof, nDofU > Kbd;
        Eigen::Matrix< double, nDofU, nWarpDof > Kdb;
      };

      const auto& indexMap = materialIndexMap();

      auto assembleAll = [&]( const WarpVector& beta, bool commit ) {
        Accumulated acc;
        acc.Pe.setZero();
        acc.Kdd.setZero();
        acc.Rb.setZero();
        acc.Kbb.setZero();
        acc.Kbd.setZero();
        acc.Kdb.setZero();

        for ( size_t iq = 0; iq < qps.size(); iq++ ) {
          auto& qp = qps[iq];

          // Restore on EVERY pass, including the committing one: computeStress
          // commits, so without this the final pass would start from whatever
          // state the last line-search trial happened to leave behind. That
          // makes the returned residual depend discontinuously on the
          // line-search path and breaks tangent/residual consistency.
          std::copy( materialStateIn[iq].begin(),
                     materialStateIn[iq].end(),
                     qp.managedStateVars->materialStateVars.data() );

          // ---- generalized B over ( displacement , warping ) ----
          BAvgSurfaceSized BPlus = BAvgSurfaceSized::Zero(), BMinus = BAvgSurfaceSized::Zero();
          BMinus.template block< nTensor, halfU >( 0, 0 )    = qp.BmatSide;
          BPlus.template block< nTensor, halfU >( 0, halfU ) = qp.BmatSide;

          const BWarpSized Bw = warpingB( qp.gradWarpModes );

          BGenSized Bgen                                                                       = BGenSized::Zero();
          Bgen.template block< nDim, nDofU >( 0, 0 )                                           = qp.NmatJump;
          Bgen.template block< nTensor, nDofU >( nDim, 0 )                                     = BPlus;
          Bgen.template block< nTensor, nDofU >( nDim + nTensor, 0 )                           = BMinus;
          Bgen.template block< nTensor, nWarpBlock >( nDim + 2 * nTensor, nDofU )              = Bw;
          Bgen.template block< nTensor, nWarpBlock >( nDim + 3 * nTensor, nDofU + nWarpBlock ) = Bw;

          Eigen::Matrix< double, nDofU + nWarpDof, 1 > q;
          q.template segment< nDofU >( 0 )        = dQ;
          q.template segment< nWarpDof >( nDofU ) = beta;

          const Eigen::Matrix< double, nGenElem, 1 > dGen = Bgen * q;

          // ---- scatter into the material's fixed 3D layout ----
          Eigen::Matrix< double, 6, 1 >  dU3d  = Eigen::Matrix< double, 6, 1 >::Zero();
          Eigen::Matrix< double, 18, 1 > dA3d  = Eigen::Matrix< double, 18, 1 >::Zero();
          Eigen::Matrix< double, 18, 1 > dW3d  = Eigen::Matrix< double, 18, 1 >::Zero();
          Eigen::Vector3d                n3d   = Eigen::Vector3d::Zero();
          Eigen::Vector3d                sep3d = Eigen::Vector3d::Zero();

          for ( int i = 0; i < nDim; ++i ) {
            n3d( i )   = qp.normal( i );
            sep3d( i ) = qp.separationVector( i );
          }

          // The material takes ( u+, u- ) and forms w = u+ - u-; the element
          // already has the jump, so put it entirely on the + side.
          for ( int i = 0; i < nDim; ++i ) {
            dU3d( i ) = dGen( i );
          }
          for ( int slot = 0; slot < 4; ++slot ) {
            for ( int i = 0; i < nDim; ++i ) {
              for ( int j = 0; j < nDim; ++j ) {
                const double value  = dGen( nDim + slot * nTensor + i * nDim + j );
                const int    flat3d = 3 * i + j;
                if ( slot == 0 ) {
                  dA3d( flat3d ) = value;
                }
                else if ( slot == 1 ) {
                  dA3d( 9 + flat3d ) = value;
                }
                else if ( slot == 2 ) {
                  dW3d( flat3d ) = value;
                }
                else {
                  dW3d( 9 + flat3d ) = value;
                }
              }
            }
          }

          Eigen::Matrix< double, nGenMat, 1 > p3d = Eigen::Matrix< double, nGenMat, 1 >::Zero();
          Eigen::Matrix< double, nGenMat, nGenMat, Eigen::RowMajor >
            K3d = Eigen::Matrix< double, nGenMat, nGenMat, Eigen::RowMajor >::Zero();

          typename Material::State         st( p3d.data(), qp.managedStateVars->materialStateVars.data() );
          typename Material::Tangents      tg( K3d.data() );
          typename Material::Deformation   df( dU3d.data(), dA3d.data(), dW3d.data(), n3d.data(), sep3d.data() );
          typename Material::TimeIncrement ti{ time, dT };

          qp.material->computeStress( st, tg, df, ti );

          // ---- gather back to the element-side generalized pair ----
          Eigen::Matrix< double, nGenElem, 1 >        p = Eigen::Matrix< double, nGenElem, 1 >::Zero();
          Eigen::Matrix< double, nGenElem, nGenElem > K = Eigen::Matrix< double, nGenElem, nGenElem >::Zero();
          for ( int a = 0; a < nGenElem; ++a ) {
            p( a ) = p3d( indexMap[a] );
            for ( int b = 0; b < nGenElem; ++b ) {
              K( a, b ) = K3d( indexMap[a], indexMap[b] );
            }
          }

          if ( commit ) {
            qp.managedStateVars->generalizedForce           = p.template segment< nDim >( 0 );
            qp.managedStateVars->surfaceStressPlus          = p.template segment< nTensor >( nDim );
            qp.managedStateVars->surfaceStressMinus         = p.template segment< nTensor >( nDim + nTensor );
            qp.managedStateVars->warpingStressSymmetric     = p.template segment< nTensor >( nDim + 2 * nTensor );
            qp.managedStateVars->warpingStressAntisymmetric = p.template segment< nTensor >( nDim + 3 * nTensor );
            qp.managedStateVars->displacement.template segment< nDim >( 0 ) += dGen.template segment< nDim >( 0 );
            qp.managedStateVars->surfaceStrain += dGen.template segment< 2 * nTensor >( nDim );
            qp.managedStateVars->warpingStrain += dGen.template segment< 2 * nTensor >( nDim + 2 * nTensor );
          }

          const auto Bd = Bgen.template block< nGenElem, nDofU >( 0, 0 );
          const auto Bb = Bgen.template block< nGenElem, nWarpDof >( 0, nDofU );

          acc.Pe -= Bd.transpose() * p * qp.J0xW;
          acc.Rb += Bb.transpose() * p * qp.J0xW;

          acc.Kdd += Bd.transpose() * K * Bd * qp.J0xW;
          acc.Kdb += Bd.transpose() * K * Bb * qp.J0xW;
          acc.Kbd += Bb.transpose() * K * Bd * qp.J0xW;
          acc.Kbb += Bb.transpose() * K * Bb * qp.J0xW;
        }
        return acc;
      };

      // ---- inner Newton on the internal warping amplitudes ----
      // beta enters the strain and therefore the plastic return mapping, so
      // R_beta is nonlinear in beta. The backtracking line search + exception
      // guard mirrors the one the MINI element needs: without it an undamped
      // step at yield onset can drive the local traction-equilibrium solve past
      // convergence, and the StressUpdateFailed propagates out as a global
      // cutback.
      WarpVector  beta = qps[0].managedStateVars->warpingAmplitudes;
      Accumulated acc  = assembleAll( beta, false );

      for ( int iteration = 0; iteration < 20; iteration++ ) {
        const double residualNorm = acc.Rb.norm();
        const double tolerance    = 1.0e-10 * std::max( 1.0, acc.Pe.norm() );
        if ( residualNorm <= tolerance ) {
          break;
        }

        // Pseudo-inverse, not an inverse: the internal block has an exactly
        // known, material-independent null space (see the file header). The
        // residual has no component along it, so the step is well defined.
        const WarpVector dBeta = -acc.Kbb.completeOrthogonalDecomposition().solve( acc.Rb );

        double alpha    = 1.0;
        bool   accepted = false;
        for ( int lineSearch = 0; lineSearch < 10; lineSearch++ ) {
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
            // a trial beta drove the local material solve past convergence:
            // treat exactly like a non-improving step and halve
          }
          alpha *= 0.5;
        }
        if ( !accepted ) {
          break; // keep the last good beta; the outer Newton continues from there
        }
      }

      acc = assembleAll( beta, true );

      // ---- static condensation ----
      const auto KbbPinv = acc.Kbb.completeOrthogonalDecomposition().pseudoInverse();

      Pe = acc.Pe + acc.Kdb * ( KbbPinv * acc.Rb );
      Ke = acc.Kdd - acc.Kdb * KbbPinv * acc.Kbd;

      qps[0].managedStateVars->warpingAmplitudes = beta;
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
      throw std::invalid_argument( "WarpingInterfaceFiniteElement: distributed loads not implemented." );
    }

    void computeBodyForce( double*, double*, const double*, const double*, double, double )
    {
      throw std::invalid_argument( "WarpingInterfaceFiniteElement: body forces not implemented." );
    }

    void computeConsistentInertia( double* )
    {
      throw std::runtime_error( "WarpingInterfaceFiniteElement: no inertia." );
    }
    void computeLumpedInertia( double* ) { throw std::runtime_error( "WarpingInterfaceFiniteElement: no inertia." ); }

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
