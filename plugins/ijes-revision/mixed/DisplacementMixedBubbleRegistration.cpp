// C3D8P: eight-node displacement hexahedron with a continuous trilinear nodal pressure and a MINI bubble, the bulk analogue of
// YIQUAD4_STABP_MINI. Unknowns [u(24); p(8)], the amplitude beta of the bubble b = (1-xi^2)(1-eta^2)(1-zeta^2) (3 values, one per
// direction) is condensed element by element.
//
//   sigma~ = P_dev sigma_material(eps) - p m,      m = (1 1 1 0 0 0)^T,  P_dev = I - m m^T / 3
//   R_u = int B^T sigma~,   R_b = int Bb^T sigma~,   R_p = -int N (m^T d_eps + dp / K) ,     K = m^T C m / 9
//
// The stress takes the TOTAL pressure and the volumetric law the INCREMENT (the same split as the interface element). The material
// commits its state in every evaluation, so every trial of the bubble iteration restarts from the state the element was given.
// The mixed stress is written back as the stress of the quadrature point; the hypoelastic material only reads its deviator, so the
// mean stress it carries is irrelevant. Optional Brezzi-Pitkaranta term: environment variable MARMOT_C3D8P_BP = coefficient
// (default 0, off): R_p -= coeff h^2 / (2 mu) int grad N . grad p.
#include "Marmot/DisplacementFiniteElement.h"
#include "Marmot/MarmotElementFactory.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

namespace {
  using Base = Marmot::Elements::DisplacementFiniteElement< 3, 8 >;

  class DisplacementMixedBubbleHex : public Base {
  public:
    static constexpr int nU = 24, nP = 8, nDof = 32;
    using V3   = Eigen::Matrix< double, 3, 1 >;
    using V6   = Eigen::Matrix< double, 6, 1 >;
    using M6   = Eigen::Matrix< double, 6, 6 >;
    using M3   = Eigen::Matrix3d;
    using Rhs  = Eigen::Matrix< double, nDof, 1 >;
    using Kmat = Eigen::Matrix< double, nDof, nDof >;

    struct Pass {
      Rhs                                Pe  = Rhs::Zero();
      Kmat                               Ke  = Kmat::Zero();
      V3                                 Rb  = V3::Zero();
      M3                                 Kbb = M3::Zero();
      Eigen::Matrix< double, 3, nU >     Kbu = Eigen::Matrix< double, 3, nU >::Zero();
      Eigen::Matrix< double, nU, 3 >     Kub = Eigen::Matrix< double, nU, 3 >::Zero();
      Eigen::Matrix< double, 3, nP >     Kbp = Eigen::Matrix< double, 3, nP >::Zero();
      Eigen::Matrix< double, nP, 3 >     Kpb = Eigen::Matrix< double, nP, 3 >::Zero();
    };

    DisplacementMixedBubbleHex( int elementID, Marmot::FiniteElement::Quadrature::IntegrationTypes integrationType, SectionType sectionType )
      : Base( elementID, integrationType, sectionType )
    {
      const char* bp = std::getenv( "MARMOT_C3D8P_BP" );
      bpCoefficient  = bp ? std::atof( bp ) : 0.0;
    }

    int getNDofPerElement() override { return nDof; }

    std::vector< std::vector< std::string > > getNodeFields() override
    {
      std::vector< std::vector< std::string > > fields;
      for ( int node = 0; node < 8; ++node )
        fields.push_back( { "displacement", "interfacePressure" } );
      return fields;
    }

    std::vector< int > getDofIndicesPermutationPattern() override
    {
      std::vector< int > permutation( nDof );
      for ( int node = 0; node < 8; ++node ) {
        for ( int c = 0; c < 3; ++c )
          permutation[3 * node + c] = 4 * node + c;
        permutation[nU + node] = 4 * node + 3;
      }
      return permutation;
    }

    int getNumberOfRequiredStateVars() override { return Base::getNumberOfRequiredStateVars() + 3; }

    void assignStateVars( double* stateVars, int nStateVars ) override
    {
      svBase = stateVars;
      nSv    = nStateVars - 3;
      alpha  = stateVars + nSv;
      Base::assignStateVars( stateVars, nSv );
    }

    void initializeYourself() override
    {
      Base::initializeYourself();
      bubbleB.clear();
      shapeN.clear();
      shapeGrad.clear();
      for ( const auto& qp : qps ) {
        const double xi = qp.xi( 0 ), eta = qp.xi( 1 ), zeta = qp.xi( 2 );
        dNdXiSized dbdxi = dNdXiSized::Zero();
        dbdxi( 0, 0 )    = -2.0 * xi * ( 1 - eta * eta ) * ( 1 - zeta * zeta );
        dbdxi( 1, 0 )    = -2.0 * eta * ( 1 - xi * xi ) * ( 1 - zeta * zeta );
        dbdxi( 2, 0 )    = -2.0 * zeta * ( 1 - xi * xi ) * ( 1 - eta * eta );
        const dNdXiSized    dNdXi = this->dNdXi( qp.xi );
        const JacobianSized J     = this->Jacobian( dNdXi );
        const JacobianSized JInv  = J.inverse();
        const BSized        Bfake = this->B( this->dNdX( dbdxi, JInv ) );
        bubbleB.push_back( Bfake.template leftCols< 3 >() );
        const auto N = this->N( qp.xi );
        Eigen::Matrix< double, nP, 1 > n;
        for ( int i = 0; i < nP; ++i )
          n( i ) = N( i );
        shapeN.push_back( n );
        shapeGrad.push_back( this->dNdX( dNdXi, JInv ) ); // (3 x 8) physical gradients of the shape functions
      }
      // mean element size for the optional pressure stabilisation
      double V = 0;
      for ( const auto& qp : qps )
        V += qp.J0xW;
      hMean = std::cbrt( V );
    }

    void computeKernels( const double* QTotal_, const double* dQ_, double* Pe_, double* Ke_, double time, double dT ) override
    {
      Eigen::Map< const Rhs > QTotal( QTotal_ );
      Eigen::Map< const Rhs > dQ( dQ_ );
      Eigen::Map< Kmat >      Ke( Ke_ );
      Eigen::Map< Rhs >       Pe( Pe_ );

      const std::vector< double > start( svBase, svBase + nSv );
      const V3                    beta0 = Eigen::Map< const V3 >( alpha );

      auto restore = [&]() { std::copy( start.begin(), start.end(), svBase ); };

      auto assemble = [&]( const V3& beta ) -> Pass {
        restore();
        Pass pass;
        const V6 m = ( V6() << 1, 1, 1, 0, 0, 0 ).finished();
        const M6 Pd = M6::Identity() - m * m.transpose() / 3.0;
        const auto Qp  = QTotal.template segment< nP >( nU );
        const auto dPp = dQ.template segment< nP >( nU );
        const auto dQu = dQ.template head< nU >();

        for ( size_t iq = 0; iq < qps.size(); ++iq ) {
          auto&          qp = qps[iq];
          const BSized&  B  = qp.B;
          const auto&    Bb = bubbleB[iq];
          const auto&    N  = shapeN[iq];
          V6             dE = B * dQu + Bb * ( beta - beta0 );

          double elasticEnergyDensity = qp.J0xW > 0.0 ? qp.managedStateVars->elasticStrainEnergy / qp.J0xW : 0.0;
          double dissipation          = qp.J0xW > 0.0 ? qp.managedStateVars->dissipation / qp.J0xW : 0.0;
          MarmotMaterialHypoElastic::state3D  state;
          MarmotMaterialHypoElastic::timeInfo timeInfo;
          state.stress               = qp.managedStateVars->stress;
          state.elasticEnergyDensity = elasticEnergyDensity;
          state.dissipation          = dissipation;
          state.stateVars            = qp.managedStateVars->materialStateVars.data();
          timeInfo.time              = time;
          timeInfo.dT                = dT;
          M6 C;
          qp.material->computeStress( state, C, dE, timeInfo );

          const double p  = N.dot( Qp );
          const double dp = N.dot( dPp );
          const V6     sigmaTilde = Pd * state.stress - p * m;
          const double K          = std::max( ( m.transpose() * C * m )( 0 ) / 9.0, 1e-30 );
          const M6     Cd         = Pd * C;
          const double w          = qp.J0xW;

          qp.managedStateVars->stress              = sigmaTilde;
          qp.managedStateVars->elasticStrainEnergy = state.elasticEnergyDensity * w;
          qp.managedStateVars->dissipation         = state.dissipation * w;
          qp.managedStateVars->totalStrainEnergy   = ( state.elasticEnergyDensity + state.dissipation ) * w;
          qp.managedStateVars->strain += dE;

          const Eigen::Matrix< double, 1, nU > mB = m.transpose() * B;
          const Eigen::Matrix< double, 1, 3 >  mBb = m.transpose() * Bb;

          pass.Pe.template head< nU >() += B.transpose() * sigmaTilde * w;
          pass.Rb += Bb.transpose() * sigmaTilde * w;
          pass.Pe.template tail< nP >() += -N * ( m.dot( dE ) + dp / K ) * w;

          pass.Ke.template block< nU, nU >( 0, 0 ) += B.transpose() * Cd * B * w;
          pass.Ke.template block< nU, nP >( 0, nU ) += -mB.transpose() * N.transpose() * w;
          pass.Ke.template block< nP, nU >( nU, 0 ) += -N * mB * w;
          pass.Ke.template block< nP, nP >( nU, nU ) += -N * N.transpose() / K * w;

          pass.Kub += B.transpose() * Cd * Bb * w;
          pass.Kbu += Bb.transpose() * Cd * B * w;
          pass.Kbb += Bb.transpose() * Cd * Bb * w;
          pass.Kbp += -mBb.transpose() * N.transpose() * w;
          pass.Kpb += -N * mBb * w;

          if ( bpCoefficient > 0.0 ) {
            const double mu    = std::max( ( Cd.template block< 1, 1 >( 3, 3 )( 0 ) ), 1e-30 ); // shear modulus (Voigt 12)
            const double coeff = bpCoefficient * hMean * hMean / ( 2.0 * mu );
            const Eigen::Matrix< double, nP, nP > Gp = shapeGrad[iq].transpose() * shapeGrad[iq];
            pass.Pe.template tail< nP >() += -coeff * Gp * ( Qp + dPp * 0.0 ) * w; // total pressure gradient
            pass.Ke.template block< nP, nP >( nU, nU ) += -coeff * Gp * w;
          }
        }
        return pass;
      };

      V3   beta = beta0;
      Pass pass = assemble( beta );
      for ( int iteration = 0; iteration < 20; ++iteration ) {
        const double residualNorm = pass.Rb.norm();
        const double tolerance    = 1.0e-10 * std::max( 1.0, pass.Pe.norm() );
        if ( residualNorm <= tolerance )
          break;
        const V3 dBeta = -pass.Kbb.fullPivLu().solve( pass.Rb );
        double   a     = 1.0;
        bool     accepted = false;
        for ( int ls = 0; ls < 10; ++ls ) {
          try {
            const Pass trial = assemble( beta + a * dBeta );
            if ( trial.Rb.norm() <= tolerance || trial.Rb.norm() < residualNorm ) {
              beta += a * dBeta;
              pass     = trial;
              accepted = true;
              break;
            }
          }
          catch ( const std::exception& ) {
          }
          a *= 0.5;
        }
        if ( !accepted )
          break;
      }
      pass = assemble( beta ); // leaves the state of the quadrature points at the accepted amplitude

      const auto                     lu        = pass.Kbb.fullPivLu();
      const V3                       KbbInvRb  = lu.solve( pass.Rb );
      const Eigen::Matrix< double, 3, nU > KbbInvKbu = lu.solve( pass.Kbu );
      const Eigen::Matrix< double, 3, nP > KbbInvKbp = lu.solve( pass.Kbp );

      Pe.template head< nU >() += pass.Pe.template head< nU >() - pass.Kub * KbbInvRb;
      Pe.template tail< nP >() += pass.Pe.template tail< nP >() - pass.Kpb * KbbInvRb;
      Ke += pass.Ke;
      Ke.template block< nU, nU >( 0, 0 ) -= pass.Kub * KbbInvKbu;
      Ke.template block< nU, nP >( 0, nU ) -= pass.Kub * KbbInvKbp;
      Ke.template block< nP, nU >( nU, 0 ) -= pass.Kpb * KbbInvKbu;
      Ke.template block< nP, nP >( nU, nU ) -= pass.Kpb * KbbInvKbp;

      Eigen::Map< V3 > alphaMap( alpha );
      alphaMap = beta;
    }

  private:
    double*                                                     svBase = nullptr;
    int                                                         nSv    = 0;
    double*                                                     alpha  = nullptr;
    double                                                      bpCoefficient = 0.0;
    double                                                      hMean         = 0.0;
    std::vector< Eigen::Matrix< double, 6, 3 > >                bubbleB;
    std::vector< Eigen::Matrix< double, nP, 1 > >               shapeN;
    std::vector< Eigen::Matrix< double, 3, nP > >               shapeGrad;
  };

  const static bool C3D8P_isRegistered = MarmotLibrary::MarmotElementFactory::registerElement(
    "C3D8P", []( int elementID ) -> MarmotElement* {
      return new DisplacementMixedBubbleHex( elementID,
                                             Marmot::FiniteElement::Quadrature::IntegrationTypes::FullIntegration,
                                             Base::SectionType::Solid );
    } );
} // namespace
