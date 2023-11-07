#include "Marmot/LinearViscoelasticPowerLaw.h"
#include "Marmot/MarmotElasticity.h"
#include "Marmot/MarmotMaterialHypoElastic.h"
#include "Marmot/MarmotTypedefs.h"
#include "Marmot/MarmotUtility.h"
#include "Marmot/MarmotVoigt.h"
#include "autodiff/forward/real.hpp"
#include <iostream>
#include <map>
#include <string>

using namespace Marmot;
using namespace Eigen;

namespace Marmot::Materials {

  LinearViscoelasticPowerLaw::LinearViscoelasticPowerLaw( const double* materialProperties,
                                                          int           nMaterialProperties,
                                                          int           materialLabel )
    : MarmotMaterialHypoElastic( materialProperties, nMaterialProperties, materialLabel ),
      // clang-format off
      // elastic parameters
      E                                 ( materialProperties[0] ),
      nu                                ( materialProperties[1] ),
      // viscoelastic parameters
      m                                 ( materialProperties[2] ),
      n                                 ( materialProperties[3] ),
      nKelvin                           ( static_cast< size_t > ( materialProperties[4] ) ),
      minTau                            ( materialProperties[5] ),
      timeToDays                        ( materialProperties[6] )
  // clang-format on
  {
    retardationTimes = KelvinChain::generateRetardationTimes( nKelvin, minTau, 10. );

    auto phi_ = [&]( autodiff::Real< powerLawApproximationOrder, double > tau ) { return phi( tau, m, n ); };

    elasticModuli = KelvinChain::computeElasticModuli< powerLawApproximationOrder >( phi_, retardationTimes );
  }

  void LinearViscoelasticPowerLaw::computeStress( double*       stress,
                                                  double*       dStressDDStrain,
                                                  const double* dStrain,
                                                  const double* timeOld,
                                                  const double  dT,
                                                  double&       pNewDT )

  {
    mVector6d nomStress( stress );
    Vector6d  dE( dStrain );
    mMatrix6d C( dStressDDStrain );

    if ( ( dE.array() == 0 ).all() && dT == 0 ) {
      C = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( E, nu );
      return;
    }

    Eigen::Ref< KelvinChain::mapStateVarMatrix > creepStateVars( ( stateVarManager->kelvinStateVars ) );

    const double dTimeDays = dT * timeToDays;

    Matrix6d CelUnitInv = ContinuumMechanics::Elasticity::Isotropic::complianceTensor( 1.0, nu );

    Vector6d creepStrainIncrement = Vector6d::Zero();
    double   creepCompliance      = 0;

    KelvinChain::evaluateKelvinChain( dTimeDays,
                                      elasticModuli,
                                      retardationTimes,
                                      creepStateVars,
                                      creepCompliance,
                                      creepStrainIncrement,
                                      1.0 );

    double effectiveCompliance = 1. / E + creepCompliance;

    C                    = ContinuumMechanics::Elasticity::Isotropic::stiffnessTensor( 1. / effectiveCompliance, nu );
    Vector6d deltaStress = C * ( dE - creepStrainIncrement );
    nomStress            = nomStress + deltaStress;

    KelvinChain::updateStateVarMatrix( dTimeDays,
                                       elasticModuli,
                                       retardationTimes,
                                       creepStateVars,
                                       deltaStress,
                                       CelUnitInv );

    return;
  }

  void LinearViscoelasticPowerLaw::assignStateVars( double* stateVars_, int nStateVars )
  {
    if ( nStateVars < getNumberOfRequiredStateVars() )
      throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__ << ": Not sufficient stateVars!" );

    this->stateVarManager = std::make_unique< LinearViscoelasticPowerLawStateVarManager >( stateVars_, nKelvin );

    MarmotMaterial::assignStateVars( stateVars_, nStateVars );
  }

  StateView LinearViscoelasticPowerLaw::getStateView( const std::string& stateName )
  {
    return stateVarManager->getStateView( stateName );
  }

  int LinearViscoelasticPowerLaw::getNumberOfRequiredStateVars()
  {
    return LinearViscoelasticPowerLawStateVarManager::layout.nRequiredStateVars + nKelvin * 6;
  }
} // namespace Marmot::Materials
