#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Fastor/Fastor.h"
#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotInterfaceMaterialHelperFunctions.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>

using namespace Eigen;
using namespace Fastor;

using Tensor1D = Fastor::Tensor< double, 3 >;
using Tensor2D = Fastor::Tensor< double, 3, 3 >;
using Tensor3D = Fastor::Tensor< double, 3, 3, 3 >;
using Tensor4D = Fastor::Tensor< double, 3, 3, 3, 3 >;

void MarmotMaterialHypoElasticInterface::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
}

// void MarmotMaterialHypoElasticInterface::computeStress( Tensor1D&  force,
//                                                         Tensor2D&  surface_stress,
//                                                         Fastor::Tensor<double, 21,21>& dStress_dStrain,
//                                                         const Fastor::Tensor<double, 6,1>& dU,
//                                                         const Fastor::Tensor<double, 18,1>& dSurface_strain,
//                                                         const Tensor1D& normal,
//                                                         const double* timeOld,
//                                                         const double  dT,
//                                                         double&       pNewDT  )
//  {};

void MarmotMaterialHypoElasticInterface::computeStress( double*       force,
                                                        double*       surface_stress,
                                                        double*       dStress_dStrain,
                                                        const double* dU,
                                                        const double* dSurface_strain,
                                                        const double* normal,
                                                        const double* timeOld,
                                                        const double  dT,
                                                        double&       pNewDT )
{
} // namespace Marmot::Materials

namespace MarmotLibrary {
  std::unordered_map< std::string, MarmotMaterialHypoElasticInterfaceFactory::materialFactoryFunction >
    MarmotMaterialHypoElasticInterfaceFactory::materialNameToModelAssociation;

  bool MarmotMaterialHypoElasticInterfaceFactory::registerMaterial( const std::string&      materialName,
                                                                    materialFactoryFunction factoryFunction )
  {
    assert( materialNameToModelAssociation.find( materialName ) == materialNameToModelAssociation.end() );

    materialNameToModelAssociation[materialName] = factoryFunction;

    return true;
  }

  MarmotMaterialHypoElasticInterface* MarmotMaterialHypoElasticInterfaceFactory::createMaterial(
    const std::string& materialName,
    const double*      materialProperties,
    int                nMaterialProperties,
    int                materialNumber )
  {
    try {
      return materialNameToModelAssociation.at(
        materialName )( materialProperties, nMaterialProperties, materialNumber );
    }
    catch ( const std::out_of_range& e ) {
      throw std::invalid_argument( MakeString() << "Invalid material " << materialName << " requested!" );
    }
  }

} // namespace MarmotLibrary
