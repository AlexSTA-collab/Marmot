/* ---------------------------------------------------------------------
 *  Marmot - BubbleOnlyInterfaceFiniteElement  (GLIQUAD4_BUBBLE)
 *  Alexandros Stathas alexandros.stathas@boku.ac.at
 * --------------------------------------------------------------------- */

/**
 * @file BubbleOnlyInterfaceFiniteElement.h
 * @brief Q1 displacement + MINI displacement bubble, condensed, with NO
 * independent pressure unknown (GLIQUAD4_BUBBLE).
 *
 * WHAT THIS IS FOR
 * ----------------
 * A controlled test of one question: does the MINI displacement BUBBLE ALONE
 * suppress the hydrostatic-stress oscillation, or does suppression require the
 * mixed pressure field it is normally paired with?
 *
 *     u_h = u_Q1 + b a ,        b = (1 - xi^2)(1 - eta^2)
 *
 * exactly the bubble shape function and bubble kinematics of
 * YIQUAD4_STABP_MINI: the amplitude perturbs the MEAN in-plane displacement, so
 * it enters A+ and A- EQUALLY. What is removed relative to that element is the
 * entire mixed-pressure apparatus -- no pressure DOF, no pressure
 * interpolation, no pressure residual, no pressure tangent blocks, no
 * Brezzi-Pitkaranta stabilisation, no pressure condensation. The only local
 * system is
 *
 *     [ K_dd  K_da ]
 *     [ K_ad  K_aa ] ,   condensed to  K_dd - K_da K_aa^-1 K_ad .
 *
 * Hydrostatic pressure is a pure POSTPROCESSING quantity here,
 * p = -tr(sigma)/3 from the station stresses; nothing in the formulation uses
 * it.
 *
 * WHY THIS MATERIAL
 * -----------------
 * It runs MarmotWarpingInterfaceMaterialHypoElasticN<5> with both warping
 * slots driven to ZERO. In that configuration the warping material reproduces
 * MarmotGaussLobattoInterfaceMaterialHypoElasticN<5> BIT-EXACTLY -- same local
 * problem, same generalized stresses, same tangent sub-block (asserted in
 * TestMarmotWarpingInterfaceMaterialHypoElastic, max|dp| = 0, max|dK| = 0).
 * So this element is GLIQUAD4 plus the bubble and nothing else, while getting
 * the material's dense generalized tangent, which makes the bubble column of
 * the generalized B operator trivial to assemble. The five-station Lobatto
 * through-thickness rule and the 2x2 Gauss surface rule are both untouched.
 *
 * The station stresses stay available as stationStress0..4, so the raw
 * material-point pressure diagnostic is identical to the one GLIQUAD4 already
 * provides -- which is what makes the comparison against the existing GLIQUAD4
 * run like-for-like.
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
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Marmot::Elements {

  class BubbleOnlyInterfaceFiniteElement : public MarmotElement, public MarmotGeometryInterfaceElement< 3, 8 > {

  public:
    enum SectionType { Interface };

    static constexpr int nStations       = 5;
    static constexpr int nDim            = 3;
    static constexpr int nNodes          = 8;
    static constexpr int nInterfaceNodes = 4;
    static constexpr int nTensor         = 9;
    static constexpr int nDofU           = 24;
    static constexpr int sizeLoadVector  = nDofU; // displacement only: no pressure DOF
    static constexpr int halfU           = nDofU / 2;

    /** The single internal field: the MINI displacement-bubble amplitude. */
    static constexpr int nBubble = nDim; // 3
    static constexpr int nAll    = nDofU + nBubble;

    using Material = MarmotWarpingInterfaceMaterialHypoElasticN< nStations >;

    static constexpr int nGen = Material::nGeneralized; // 39; the warping slots stay zero

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
    using BubbleVector  = Eigen::Matrix< double, nBubble, 1 >;
    using BubbleMatrix  = Eigen::Matrix< double, nBubble, nBubble >;
    using BGenSized     = Eigen::Matrix< double, nGen, nAll >;

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
      VectorDim            gradBubble;
      double               bubble = 0.0;

      class QPStateVarManager : public MarmotStateVarVectorManager {
        static constexpr int nRaw = nDim + nTensor + nTensor + nBubble + 2 * nDim + 2 * nTensor;

        inline const static auto layout = makeLayout( {
          { .name = "generalizedForce", .length = nDim },
          { .name = "surfaceStressPlus", .length = nTensor },
          { .name = "surfaceStressMinus", .length = nTensor },
          { .name = "bubbleAlpha", .length = nBubble },
          { .name = "displacement", .length = 2 * nDim },
          { .name = "surfaceStrain", .length = 2 * nTensor },
          { .name = "state block alignment padding", .length = ( 4 - ( nRaw % 4 ) ) % 4 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< Eigen::Matrix< double, nDim, 1 > >        generalizedForce;
        Eigen::Map< Eigen::Matrix< double, nTensor, 1 > >     surfaceStressPlus;
        Eigen::Map< Eigen::Matrix< double, nTensor, 1 > >     surfaceStressMinus;
        Eigen::Map< Eigen::Matrix< double, nBubble, 1 > >     bubbleAlpha;
        Eigen::Map< Eigen::Matrix< double, 2 * nDim, 1 > >    displacement;
        Eigen::Map< Eigen::Matrix< double, 2 * nTensor, 1 > > surfaceStrain;
        Eigen::Map< Eigen::VectorXd >                         materialStateVars;

        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

        QPStateVarManager( double* v, int n )
          : MarmotStateVarVectorManager( v, layout ),
            generalizedForce( &find( "generalizedForce" ) ),
            surfaceStressPlus( &find( "surfaceStressPlus" ) ),
            surfaceStressMinus( &find( "surfaceStressMinus" ) ),
            bubbleAlpha( &find( "bubbleAlpha" ) ),
            displacement( &find( "displacement" ) ),
            surfaceStrain( &find( "surfaceStrain" ) ),
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
      }
    };

    std::vector< QuadraturePoint > qps;

    /** Diagnostic ablation, per INSTANCE (not a cached global): the element
     *  reduces to plain GLIQUAD4. Settable from the constructor, or globally
     *  with MARMOT_BUBBLE_ONLY_MODE=nobubble. */
    const bool bubbleDisabledFlag;

    bool bubbleDisabled() const { return bubbleDisabledFlag; }

    static bool bubbleDisabledFromEnvironment()
    {
      const char* v = std::getenv( "MARMOT_BUBBLE_ONLY_MODE" );
      return v && std::string( v ) == "nobubble";
    }

    BubbleOnlyInterfaceFiniteElement( int                                         elementID,
                                      FiniteElement::Quadrature::IntegrationTypes integrationType,
                                      SectionType                                 sectionType_  = Interface,
                                      bool                                        disableBubble = false )
      : ParentGeometryElement(),
        bubbleDisabledFlag( disableBubble || bubbleDisabledFromEnvironment() ),
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
        for ( int i = 0; i < sizeLoadVector; i++ ) {
          perm.push_back( i );
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
          throw std::invalid_argument( "GLIQUAD4_BUBBLE: non-positive normal separation." );
        }

        qp.J0xW  = qp.weight * qp.sqrtDetG * thickness;
        qp.hElem = std::sqrt( qp.sqrtDetG );

        // exactly the MINI bubble and its surface gradient
        const double x = qp.xi( 0 ), e = qp.xi( 1 );
        qp.bubble = ( 1.0 - x * x ) * ( 1.0 - e * e );
        Eigen::Matrix< double, 2, 1 > dbdxi;
        dbdxi( 0 )    = -2.0 * x * ( 1.0 - e * e );
        dbdxi( 1 )    = -2.0 * e * ( 1.0 - x * x );
        qp.gradBubble = qp.J * qp.G.inverse() * dbdxi;

        if ( qp.material ) {
          qp.material->setCharacteristicElementLength( qp.hElem );
        }
      }
    }

    void setInitialConditions( StateTypes state, const double* )
    {
      if ( state != MarmotElement::MarmotMaterialInitialization ) {
        throw std::invalid_argument( "GLIQUAD4_BUBBLE: invalid initial condition." );
      }
      for ( auto& qp : qps ) {
        qp.material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                         qp.managedStateVars->materialStateVars.size() );
      }
    }

    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT )
    {
      (void)QTotal_;

      Eigen::Map< const RhsSized > dQ( dQ_ );
      Eigen::Map< KeSizedMatrix >  Ke( Ke_ );
      Eigen::Map< RhsSized >       Pe( Pe_ );

      std::vector< std::vector< double > > materialStateIn( qps.size() );
      for ( size_t i = 0; i < qps.size(); i++ ) {
        auto& m            = qps[i].managedStateVars->materialStateVars;
        materialStateIn[i] = std::vector< double >( m.data(), m.data() + m.size() );
      }

      struct Accumulated {
        RhsSized                                Pe;
        KeSizedMatrix                           Kdd;
        BubbleVector                            Ra;
        BubbleMatrix                            Kaa;
        Eigen::Matrix< double, nBubble, nDofU > Kad;
        Eigen::Matrix< double, nDofU, nBubble > Kda;
      };

      auto assembleAll = [&]( const BubbleVector& alpha, bool commit ) {
        Accumulated acc;
        acc.Pe.setZero();
        acc.Kdd.setZero();
        acc.Ra.setZero();
        acc.Kaa.setZero();
        acc.Kad.setZero();
        acc.Kda.setZero();

        for ( size_t iq = 0; iq < qps.size(); iq++ ) {
          auto& qp = qps[iq];

          std::copy( materialStateIn[iq].begin(),
                     materialStateIn[iq].end(),
                     qp.managedStateVars->materialStateVars.data() );

          BAvgSurfaceSized BPlus = BAvgSurfaceSized::Zero(), BMinus = BAvgSurfaceSized::Zero();
          BMinus.block< nTensor, halfU >( 0, 0 )    = qp.BmatSide;
          BPlus.block< nTensor, halfU >( 0, halfU ) = qp.BmatSide;

          // MINI bubble: ubar += b*alpha => A±_ij += alpha_i (grad_s b)_j, on BOTH stations
          Eigen::Matrix< double, nTensor, nBubble > Bbub = Eigen::Matrix< double, nTensor, nBubble >::Zero();
          if ( !bubbleDisabled() ) {
            for ( int i = 0; i < nDim; i++ ) {
              for ( int j = 0; j < nDim; j++ ) {
                Bbub( nDim * i + j, i ) = qp.gradBubble( j );
              }
            }
          }

          BGenSized Bgen                                                        = BGenSized::Zero();
          Bgen.block< nDim, nDofU >( Material::offsetJump, 0 )                  = qp.NmatJump;
          Bgen.block< nTensor, nDofU >( Material::offsetSurfacePlus, 0 )        = BPlus;
          Bgen.block< nTensor, nBubble >( Material::offsetSurfacePlus, nDofU )  = Bbub;
          Bgen.block< nTensor, nDofU >( Material::offsetSurfaceMinus, 0 )       = BMinus;
          Bgen.block< nTensor, nBubble >( Material::offsetSurfaceMinus, nDofU ) = Bbub;
          // the warping rows stay ZERO -> the material reduces to Gauss-Lobatto

          Eigen::Matrix< double, nAll, 1 > q;
          q.segment< nDofU >( 0 )       = dQ;
          q.segment< nBubble >( nDofU ) = alpha;

          const Eigen::Matrix< double, nGen, 1 > dGen = Bgen * q;

          Eigen::Matrix< double, 6, 1 >  dU6 = Eigen::Matrix< double, 6, 1 >::Zero();
          Eigen::Matrix< double, 18, 1 > dA  = Eigen::Matrix< double, 18, 1 >::Zero();
          Eigen::Matrix< double, 18, 1 > dW  = Eigen::Matrix< double, 18, 1 >::Zero();

          dU6.segment< 3 >( 0 ) = dGen.segment< nDim >( Material::offsetJump );
          dA.segment< 9 >( 0 )  = dGen.segment< nTensor >( Material::offsetSurfacePlus );
          dA.segment< 9 >( 9 )  = dGen.segment< nTensor >( Material::offsetSurfaceMinus );

          Eigen::Matrix< double, nGen, 1 > p = Eigen::Matrix< double, nGen, 1 >::Zero();
          Eigen::Matrix< double, nGen, nGen, Eigen::RowMajor >
            K = Eigen::Matrix< double, nGen, nGen, Eigen::RowMajor >::Zero();

          typename Material::State         st( p.data(), qp.managedStateVars->materialStateVars.data() );
          typename Material::Tangents      tg( K.data() );
          typename Material::Deformation   df( dU6.data(),
                                             dA.data(),
                                             dW.data(),
                                             qp.normal.data(),
                                             qp.separationVector.data() );
          typename Material::TimeIncrement ti{ time, dT };

          qp.material->computeStress( st, tg, df, ti );

          if ( commit ) {
            qp.managedStateVars->generalizedForce   = p.segment< nDim >( Material::offsetJump );
            qp.managedStateVars->surfaceStressPlus  = p.segment< nTensor >( Material::offsetSurfacePlus );
            qp.managedStateVars->surfaceStressMinus = p.segment< nTensor >( Material::offsetSurfaceMinus );
            qp.managedStateVars->displacement.segment< nDim >( 0 ) += dGen.segment< nDim >( Material::offsetJump );
            qp.managedStateVars->surfaceStrain += dGen.segment< 2 * nTensor >( Material::offsetSurfacePlus );
          }

          const auto Bd = Bgen.block< nGen, nDofU >( 0, 0 );
          const auto Ba = Bgen.block< nGen, nBubble >( 0, nDofU );

          acc.Pe -= Bd.transpose() * p * qp.J0xW;
          acc.Ra += Ba.transpose() * p * qp.J0xW;

          acc.Kdd += Bd.transpose() * K * Bd * qp.J0xW;
          acc.Kda += Bd.transpose() * K * Ba * qp.J0xW;
          acc.Kad += Ba.transpose() * K * Bd * qp.J0xW;
          acc.Kaa += Ba.transpose() * K * Ba * qp.J0xW;
        }
        return acc;
      };

      if ( bubbleDisabled() ) {
        Accumulated acc                      = assembleAll( BubbleVector::Zero(), true );
        Pe                                   = acc.Pe;
        Ke                                   = acc.Kdd;
        qps[0].managedStateVars->bubbleAlpha = BubbleVector::Zero();
        return;
      }

      // ---- local bubble equilibrium: r_a(d, a) = 0 ----
      BubbleVector alpha = qps[0].managedStateVars->bubbleAlpha;
      Accumulated  acc   = assembleAll( alpha, false );

      for ( int iteration = 0; iteration < 20; iteration++ ) {
        const double residualNorm = acc.Ra.norm();
        const double tolerance    = 1.0e-10 * std::max( 1.0, acc.Pe.norm() );
        if ( residualNorm <= tolerance ) {
          break;
        }

        const BubbleVector dAlpha = -acc.Kaa.fullPivLu().solve( acc.Ra );

        double step     = 1.0;
        bool   accepted = false;
        for ( int lineSearch = 0; lineSearch < 20; lineSearch++ ) {
          try {
            Accumulated trial = assembleAll( alpha + step * dAlpha, false );
            if ( trial.Ra.norm() <= tolerance || trial.Ra.norm() < residualNorm ) {
              alpha += step * dAlpha;
              acc      = trial;
              accepted = true;
              break;
            }
          }
          catch ( const std::exception& ) {
          }
          step *= 0.5;
        }
        if ( !accepted ) {
          break;
        }
      }

      acc = assembleAll( alpha, true );

      const BubbleMatrix KaaInverse = acc.Kaa.fullPivLu().inverse();

      Pe = acc.Pe + acc.Kda * ( KaaInverse * acc.Ra );
      Ke = acc.Kdd - acc.Kda * KaaInverse * acc.Kad;

      qps[0].managedStateVars->bubbleAlpha = alpha;
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
      throw std::invalid_argument( "GLIQUAD4_BUBBLE: distributed loads not implemented." );
    }
    void computeBodyForce( double*, double*, const double*, const double*, double, double )
    {
      throw std::invalid_argument( "GLIQUAD4_BUBBLE: body forces not implemented." );
    }
    void computeConsistentInertia( double* ) { throw std::runtime_error( "GLIQUAD4_BUBBLE: no inertia." ); }
    void computeLumpedInertia( double* ) { throw std::runtime_error( "GLIQUAD4_BUBBLE: no inertia." ); }

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
