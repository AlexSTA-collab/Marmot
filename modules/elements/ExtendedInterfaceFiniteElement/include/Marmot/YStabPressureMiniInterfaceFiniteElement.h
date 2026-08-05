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
#include <cmath>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Marmot::Elements {

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

    static double defaultStabGamma()
    {
      if ( const char* e = std::getenv( "MARMOT_STABP_GAMMA" ) ) {
        const double v = std::atof( e );
        if ( v > 0.0 )
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
        RhsSized                          Pe;
        KeSizedMatrix                     Ke;
        V3                                Rb;
        M33                               Kbb;
        Eigen::Matrix< double, 3, nDofU > Kbu;
        Eigen::Matrix< double, nDofU, 3 > Kub;
        Eigen::Matrix< double, 3, nDofP > Kbp;
        Eigen::Matrix< double, nDofP, 3 > Kpb;
      };

      auto assembleAll = [&]( const V3& beta, bool commit ) {
        Acc A;
        A.Pe.setZero();
        A.Ke.setZero();
        A.Rb.setZero();
        A.Kbb.setZero();
        A.Kbu.setZero();
        A.Kub.setZero();
        A.Kbp.setZero();
        A.Kpb.setZero();

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
          qp.material->computeStress( msv.data(), resp, tg, def, tinc );

          if ( commit ) {
            qp.managedStateVars->generalizedForce        = f;
            qp.managedStateVars->surfaceStressPlus       = Sp;
            qp.managedStateVars->surfaceStressMinus      = Sm;
            qp.managedStateVars->volumetricResidual( 0 ) = rp;
          }

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
        return A;
      };

      // ---- inner Newton on the internal bubble amplitude: solve R_beta(beta)=0 ----
      // Necessary because beta enters the strain and therefore the plastic
      // return-mapping, so R_beta is NONLINEAR in beta. A single step from beta=0
      // leaves Pe (evaluated at beta=0) inconsistent with the condensed Ke, which
      // shows up directly as a large tangent-vs-FD error on the displacement columns.
      V3  beta = qps[0].managedStateVars->bubbleAlpha;
      Acc A    = assembleAll( beta, false );
      for ( int it = 0; it < 20; it++ ) {
        const double rn = A.Rb.norm();
        if ( rn <= 1.0e-12 * std::max( 1.0, A.Pe.norm() ) )
          break;
        beta -= A.Kbb.fullPivLu().solve( A.Rb );
        A = assembleAll( beta, false );
      }
      // final pass at the converged beta, committing material state
      A = assembleAll( beta, true );

      const M33 KbbInv = A.Kbb.fullPivLu().inverse();
      Pe               = A.Pe;
      Ke               = A.Ke;
      // R_beta is ~0 now, so the condensed residual equals Pe; only the tangent
      // needs the Schur complement.
      Pe.segment< nDofU >( offU ) += A.Kub * ( KbbInv * A.Rb );
      Pe.segment< nDofP >( offP ) += A.Kpb * ( KbbInv * A.Rb );
      Ke.block< nDofU, nDofU >( offU, offU ) -= A.Kub * KbbInv * A.Kbu;
      Ke.block< nDofU, nDofP >( offU, offP ) -= A.Kub * KbbInv * A.Kbp;
      Ke.block< nDofP, nDofU >( offP, offU ) -= A.Kpb * KbbInv * A.Kbu;
      Ke.block< nDofP, nDofP >( offP, offP ) -= A.Kpb * KbbInv * A.Kbp;

      qps[0].managedStateVars->bubbleAlpha = beta;
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
