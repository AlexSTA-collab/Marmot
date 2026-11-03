#pragma once
// Equilibrated two-sided hypoelastic interface material whose two faces carry their own complete property sets, so that
// the two halves of the layer may differ in every parameter of the base material (for instance in the rheology of a
// viscoelastic layer), not only in E and nu as in MarmotEquilibratedXFacePropsInterfaceMaterialHypoElastic.
//
// Properties: [E+, nu+, h, P+ (k entries), E-, nu-, P- (k entries)], k = (N - 5) / 2. The base material of the bottom
// face receives [E-, nu-, P-]. Both faces must require the same number of state variables.
#include "Marmot/MarmotEquilibratedXInterfaceMaterialHypoElastic.h"
#include "Marmot/MarmotMaterialHypoElasticFactory.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

class MarmotEquilibratedXTwoPhaseInterfaceMaterialHypoElastic : public MarmotEquilibratedXInterfaceMaterialHypoElastic {
  std::vector< double > bottomProperties; // kept as a member: the base material keeps the pointer to its properties

  static int topCount( int n )
  {
    return ( n >= 5 && ( n - 5 ) % 2 == 0 ) ? 3 + ( n - 5 ) / 2 : n;
  }

public:
  MarmotEquilibratedXTwoPhaseInterfaceMaterialHypoElastic( const std::string& materialName,
                                                          const double*      matProperties,
                                                          int                nMaterialProperties,
                                                          int                materialNumber_ )
    : MarmotEquilibratedXInterfaceMaterialHypoElastic( materialName,
                                                       matProperties,
                                                       topCount( nMaterialProperties ),
                                                       materialNumber_ )
  {
    if ( nMaterialProperties < 5 || ( nMaterialProperties - 5 ) % 2 != 0 )
      throw std::invalid_argument( "MarmotEquilibratedXTwoPhaseInterfaceMaterialHypoElastic requires "
                                   "[E+, nu+, h, P+(k), E-, nu-, P-(k)], i.e. an odd number of properties >= 5." );
    const int k = ( nMaterialProperties - 5 ) / 2;
    bottomProperties.assign( matProperties + 3 + k, matProperties + nMaterialProperties );
    bottomMaterial = std::unique_ptr< MarmotMaterialHypoElastic >(
      MarmotLibrary::MarmotMaterialHypoElasticFactory::createMaterial( materialName,
                                                                       bottomProperties.data(),
                                                                       static_cast< int >( bottomProperties.size() ),
                                                                       materialNumber ) );
    if ( bottomMaterial->getNumberOfRequiredStateVars() != topMaterial->getNumberOfRequiredStateVars() )
      throw std::invalid_argument( "The two faces of the two-phase interface material must require the same number of "
                                   "state variables." );
  }
};
