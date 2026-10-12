#include "Marmot/MarmotGradedInterfaceMaterialFiniteStrainBase.h"

#include "Marmot/MarmotInterfaceFiniteStrainHelperFunctions.h"

using namespace Marmot::FastorStandardTensors;

MarmotGradedInterfaceMaterialFiniteStrainBase::BulkResponse MarmotGradedInterfaceMaterialFiniteStrainBase::evaluateBulk(
  int                                              q,
  const double*                                    stateVars,
  const Tensor33d&                                 F,
  const MarmotMaterialFiniteStrain::TimeIncrement& timeIncrement ) const
{
  const auto& point = points[q];

  BulkResponse bulk;
  bulk.state.assign( stateVars + point.offset, stateVars + point.offset + 3 + point.nBaseStateVars );

  MarmotMaterialFiniteStrain::ConstitutiveResponse< 3 > response( Tensor33d( 0.0 ), 0.0, 0.0, bulk.state.data() + 3 );
  MarmotMaterialFiniteStrain::AlgorithmicModuli< 3 >    moduli;
  const MarmotMaterialFiniteStrain::Deformation< 3 >    deformation{ F };

  point.material->computeStress( response, moduli, deformation, timeIncrement );

  Marmot::Materials::InterfaceFiniteStrainHelperFunctions::firstPiolaKirchhoff( response.tau,
                                                                                moduli.dTau_dF,
                                                                                F,
                                                                                bulk.P,
                                                                                bulk.A );
  return bulk;
}
