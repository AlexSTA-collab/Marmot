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
    
    mVector6d  F( force );
    mVector18d  S( surface_stress );
    mMatrix21d  dS_dE( dStress_dStrain );
    
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
    //elastic step
    S     = S + Cel * dE;
    dS_dE = Cel;
} // namespace Marmot::Materials
