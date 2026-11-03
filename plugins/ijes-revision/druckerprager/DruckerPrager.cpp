#include "Marmot/DruckerPrager.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotExceptions.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include <cmath>
#include <stdexcept>

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  double DruckerPragerModel::eta( double angleInDegrees )
  {
    const double t = std::tan( angleInDegrees * Constants::Pi / 180. );
    return 3. * t / std::sqrt( 9. + 12. * t * t );
  }

  double DruckerPragerModel::xi( double angleInDegrees )
  {
    const double t = std::tan( angleInDegrees * Constants::Pi / 180. );
    return 3. / std::sqrt( 9. + 12. * t * t );
  }

  double DruckerPragerModel::getDensity( const double* stateVars ) const
  {
    if ( this->nMaterialProperties < 8 )
      throw std::runtime_error( "DruckerPragerModel: no density given (nMaterialProperties < 8)" );
    return this->materialProperties[7];
  }

  DruckerPragerModel::DruckerPragerModel( const double* materialProperties,
                                          const int     nMaterialProperties,
                                          const int     materialLabel )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel )
  {
    if ( nMaterialProperties < 6 )
      throw std::invalid_argument( "DruckerPragerModel: properties are E, nu, c0, phi [deg], psi [deg], H" );
    const double phi = materialProperties[3], psi = materialProperties[4];
    // psi > phi is allowed: phi = 0 < psi is the non-associated, dilatant von Mises model
    if ( !( phi >= 0. && phi < 90. && psi >= 0. && psi < 90. ) )
      throw std::invalid_argument( "DruckerPragerModel: need 0 <= phi, psi < 90 deg" );
    if ( materialProperties[2] <= 0. || materialProperties[5] < 0. )
      throw std::invalid_argument( "DruckerPragerModel: need c0 > 0 and H >= 0" );
    // the tension cap p <= pt must cut the cone below its apex xi c0 / eta, so that the apex is never reached
    if ( nMaterialProperties >= 7 && phi > 0. && !( materialProperties[6] < xi( phi ) * materialProperties[2] / eta( phi ) ) )
      throw std::invalid_argument( "DruckerPragerModel: the tension cap pt must lie below the apex xi c0 / eta" );
    stateLayout.add( "kappa", 1 );
    stateLayout.finalize();
  }

  void DruckerPragerModel::computeStress( state3D&        state,
                                          Matrix6d&       dStress_dStrain,
                                          const Vector6d& dStrain,
                                          const timeInfo& timeInfo ) const
  {
    const double& E  = this->materialProperties[0];
    const double& nu = this->materialProperties[1];
    const double& c0 = this->materialProperties[2];
    const double& H  = this->materialProperties[5];
    const double  et = eta( this->materialProperties[3] ), etb = eta( this->materialProperties[4] );
    const double  x  = xi( this->materialProperties[3] );
    const double  G = E / ( 2. * ( 1. + nu ) ), K = E / ( 3. * ( 1. - 2. * nu ) );

    using namespace ContinuumMechanics::VoigtNotation;
    mVector6d  S( state.stress.data() );
    mMatrix6d  dS_dE( dStress_dStrain.data() );
    const auto Cel = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );

    double&    kappa = stateLayout.getAs< double& >( state.stateVars, "kappa" );
    const auto c     = [&]( double k ) { return c0 + H * k; };

    const Vector6d trial  = S + Cel * dStrain;
    const Vector6d sTrial = IDev * trial;
    const double   pTrial = trial.head< 3 >().sum() / 3.;
    const double   sqJ2   = std::sqrt( Invariants::J2( trial ) );

    const bool   cap = this->nMaterialProperties >= 7; // tension cap p <= pt, associated (volumetric) flow
    const double pt  = cap ? this->materialProperties[6] : 0.;
    const double f1t = sqJ2 + et * pTrial - x * c( kappa );
    const double f2t = cap ? pTrial - pt : -1.;

    if ( f1t <= 0. && f2t <= 0. ) { // elastic
      S     = trial;
      dS_dE = Cel;
      state.elasticEnergyDensity += 0.5 * dStrain.dot( Cel * dStrain );
      return;
    }

    const Vector6d oldStress  = S;
    Matrix6d       IDevStrain = IDev; // maps the engineering strain to the deviatoric strain tensor
    IDevStrain.block< 6, 3 >( 0, 3 ) *= 0.5;
    const Vector6d n = sqJ2 > 0. ? Vector6d( sTrial / ( std::sqrt( 2. ) * sqJ2 ) ) : Vector6d::Zero(); // ||n|| = 1

    // 1. return to the smooth part of the cone: Phi(dGamma) = sqrt(J2t) - G dGamma + eta (pt - K etaBar dGamma) - xi c
    //    is linear in dGamma (linear hardening), so its root is exact; an iteration against an absolute tolerance
    //    fails spuriously once round-off in the trial stress exceeds it (stiff layers, large Newton trial strains).
    //    Valid if it stays on the smooth part (s does not reverse) and below the cap.
    const double dGamma = f1t > 0. ? f1t / ( G + K * et * etb + x * x * H ) : 0.;
    if ( f1t > 0. && sqJ2 - G * dGamma >= 0. && ( !cap || pTrial - K * etb * dGamma <= pt ) ) {
      S = ( 1. - G * dGamma / sqJ2 ) * sTrial + ( pTrial - K * etb * dGamma ) * I;
      kappa += x * dGamma;

      // de Souza Neto Eq. (8.85): consistent tangent of the cone return, engineering shear in the strain
      const double A = 1. / ( G + K * et * etb + x * x * H );
      const double a = G * dGamma / sqJ2; // dGamma / (sqrt2 ||e_d,trial||) of de Souza Neto
      dS_dE          = 2. * G * ( 1. - a ) * IDevStrain + 2. * G * ( a - G * A ) * ( n * n.transpose() ) -
              std::sqrt( 2. ) * G * A * K * ( et * n * I.transpose() + etb * I * n.transpose() ) +
              K * ( 1. - K * et * etb * A ) * ( I * I.transpose() );
    }
    else if ( cap ) {
      // the cap is active (pTrial > pt, or the cone return ends above the cap): p = pt. Below the cap the cone
      // return is always valid (the cap lies below the apex), so pTrial > pt here
      if ( f2t <= 0. )
        throw StressUpdateFailed( "DruckerPragerModel: cap branch reached below the cap" );
      if ( sqJ2 + et * pt - x * c( kappa ) <= 0. ) {
        // 2. cap only: purely volumetric flow, the deviator stays elastic and the stress keeps its shear stiffness
        S     = sTrial + pt * I;
        dS_dE = 2. * G * IDevStrain;
      }
      else {
        // 3. corner of cone and cap: p = pt fixes the pressure, the deviator returns radially on the cone,
        //    sqrt(J2t) - G dG + eta pt - xi c(kappa + xi dG) = 0, and the cap multiplier takes the rest of the
        //    volumetric strain, dLambda = (pTrial - pt) / K - etaBar dG >= 0
        const double dG      = ( sqJ2 + et * pt - x * c( kappa ) ) / ( G + x * x * H );
        const double dLambda = ( pTrial - pt ) / K - etb * dG;
        if ( dG < 0. || dLambda < 0. || sqJ2 - G * dG < 0. )
          throw StressUpdateFailed( "DruckerPragerModel: no admissible return at the corner of cone and cap" );
        S = ( 1. - G * dG / sqJ2 ) * sTrial + pt * I;
        kappa += x * dG;
        const double A = 1. / ( G + x * x * H );
        const double a = G * dG / sqJ2;
        dS_dE          = 2. * G * ( 1. - a ) * IDevStrain + 2. * G * ( a - G * A ) * ( n * n.transpose() );
      }
    }
    else {
      // return to the apex: r(dEpsV) = beta c(kappa + alpha dEpsV) - pt + K dEpsV, s = 0
      // (unreachable for phi = 0: sqrt(J2t) - G dGamma = sqrt(J2t) - G f_t / (G + xi^2 H) > 0 for c > 0)
      if ( et <= 0. )
        throw StressUpdateFailed( "DruckerPragerModel: apex return needs phi > 0" );
      if ( etb <= 0. ) {
        // tension cut-off for psi = 0 without a cap: the non-dilatant flow has no volumetric part and cannot bring
        // the stress back to the apex; the stress goes to the apex of the current cone, the excess volumetric
        // strain opens without hardening (kappa unchanged), and the stress no longer responds to the strain
        S     = ( x * c( kappa ) / et ) * I;
        dS_dE = Matrix6d::Zero();
      }
      else {
        const double alpha = x / etb, beta = x / et;
        const double dEv   = ( pTrial - beta * c( kappa ) ) / ( alpha * beta * H + K ); // r is linear: exact root
        S                  = ( pTrial - K * dEv ) * I;
        kappa += alpha * dEv;
        dS_dE = K * ( 1. - K / ( K + alpha * beta * H ) ) * ( I * I.transpose() );
      }
    }
    // trapezoidal increments of the stress power
    state.elasticEnergyDensity += 0.5 * dStrain.dot( Cel * dStrain );
    state.dissipation += 0.5 * ( S + oldStress ).dot( dStrain ) - 0.5 * dStrain.dot( Cel * dStrain );
  }

} // namespace Marmot::Materials
