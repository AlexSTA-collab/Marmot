#include "Marmot/InterfaceVonMises.h"
#include "Marmot/HaighWestergaard.h"
#include "Marmot/MarmotConstants.h"
#include "Marmot/InterfaceMarmotElasticity.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotVoigt.h"
#include "Marmot/InterfaceVonMisesConstants.h"
#include <iostream>
#include <map>

namespace Marmot::Materials {

  using namespace Eigen;
  using namespace Marmot;

  void InterfaceVonMisesModel::assignStateVars( double* stateVars, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    managedStateVars = std::make_unique< InterfaceVonMisesModelStateVarManager >( stateVars );
    return MarmotMaterialHypoElastic::assignStateVars( stateVars, nStateVars );
  }

  StateView InterfaceVonMisesModel::getStateView( const std::string& stateName )
  {
    return managedStateVars->getStateView( stateName );
  }

  double InterfaceVonMisesModel::getDensity()
  {
    if ( this->nMaterialProperties < 11 )
      throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": No density given! nMaterialProperties < 11" );
    return this->materialProperties[10];
  }

  //double InterfaceVonMisesModel::getNormal()
  //{
  //  if ( this->nMaterialProperties < 12 )
  //    throw std::runtime_error( MakeString() << __PRETTY_FUNCTION__ << ": No normal given! nMaterialProperties < 12" );
  //  return this->materialProperties[11];
  //}

  void InterfaceVonMisesModel::computeStress( double*  force,
                                     double*       surface_stress,
                                     double*       dStress_dStrain,
                                     const double* dU,
                                     const double* dStrain,
                                     const double* normal,
                                     const double* timeOld,
                                     const double  dT,
                                     double&       pNewDT )

  {
    // elasticity parameters
    const double& E_M  = this->materialProperties[0];
    const double& nu_M = this->materialProperties[1];
    const double& E_I  = this->materialProperties[2];
    const double& nu_I = this->materialProperties[3];
    const double& E_0  = this->materialProperties[4];
    const double& nu_0 = this->materialProperties[5];
    
    // plasticity parameters we assume a common yield response but the yield function is not known
    const double& yieldStress      = this->materialProperties[6];
    const double& HLin             = this->materialProperties[7];
    const double& deltaYieldStress = this->materialProperties[8];
    const double& delta            = this->materialProperties[9];

    // map to force, surface stress, displacement, strain, normal and tangent stiffness 
    
    mVector3d  F( force, surface_stress );
    mVector6d  S( force, surface_stress );
    mMatrix6d  dS_dE( dStress_dStrain );
    
    mMatrix15d  dS_dE( dStress_dStrain ); //dimensions x(3+6+6)
    const auto dE = Map< const Vector9d >( dStrain );

    // compute elastic stiffness
    const auto Cel_M = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E_M, nu_M );
    const auto Cel_I = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E_I, nu_I );
    const auto Cel_0 = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E_0, nu_0 );
    
    // Compute characteristic interface moduli, see interface module for geometric properties, during call to the material the normal vector is provided. 
    
    // handle zero strain increment
    if ( dE.isZero( 1e-14 ) ) {
      dS_dE = Cel;
      return;
    }

    // get current hardening variable
    double& kappa = managedStateVars->kappa;

    // isotropic hardening law
    auto fy = [&]( double kappa_ ) {
      return yieldStress + HLin * kappa_ + deltaYieldStress * ( 1. - std::exp( -delta * kappa_ ) );
    };

    // derivative of fy wrt dKappa
    auto dfy_ddKappa = [&]( double kappa_ ) { return HLin + deltaYieldStress * delta * std::exp( -delta * kappa_ ); };

    // yield function
    auto f = [&]( double rho_, double kappa_ ) { return rho_ - Constants::sqrt2_3 * fy( kappa_ ); };

    // compute elastic predictor
    const Vector6d trialStress = S + Cel * dE;

    using namespace ContinuumMechanics::VoigtNotation;
    const double rhoTrial = std::sqrt( 2. * Invariants::J2( trialStress ) );

    if ( f( rhoTrial, kappa ) >= 0.0 ) {
      // plastic step
      const double G = E / ( 2. * ( 1. + nu ) );

      auto g = [&]( double deltaKappa ) {
        return rhoTrial - Constants::sqrt6 * G * deltaKappa - Constants::sqrt2_3 * fy( kappa + deltaKappa );
      };

      // variables for return mapping
      int    counter    = 0;
      double dKappa     = 0;
      double dLambda    = 0;
      double dg_ddKappa = 0;

      // compute return mapping direction
      Vector6d n = ContinuumMechanics::VoigtNotation::IDev * trialStress / rhoTrial;

      while ( std::abs( g( dKappa ) ) > VonMisesConstants::innerNewtonTol ) {

        if ( counter == VonMisesConstants::nMaxInnerNewtonCycles ) {
          pNewDT = 0.5;
          return;
        }
        // compute derivative of g wrt kappa
        dg_ddKappa = -Constants::sqrt6 * G - Constants::sqrt2_3 * dfy_ddKappa( kappa + dKappa );

        // update dKappa and iteration counter
        dKappa -= g( dKappa ) / dg_ddKappa;
        counter += 1;
      }

      dLambda = Constants::sqrt3_2 * dKappa;

      // update material state
      S     = trialStress - 2. * G * dLambda * n;
      kappa = kappa + dKappa;

      // compute consistent tangent in Voigt Notation
      Matrix6d IDevHalfShear = ContinuumMechanics::VoigtNotation::IDev;
      IDevHalfShear.block< 6, 3 >( 0, 3 ) *= 0.5;

      dS_dE = Cel -
              2. * G * ( 1. / ( 1. + dfy_ddKappa( kappa ) / ( 3. * G ) ) - 2. * G * dLambda / rhoTrial ) *
                ( n * n.transpose() ) -
              4. * G * G * dLambda / rhoTrial * IDevHalfShear;
    }
    else {
      // elastic step
      S     = trialStress;
      dS_dE = Cel;
    }
  }

} // namespace Marmot::Materials
