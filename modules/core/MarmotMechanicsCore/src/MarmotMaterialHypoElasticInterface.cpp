#include "Marmot/MarmotMaterialHypoElasticInterface.h"
#include "Marmot/HughesWinget.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotKinematics.h"
#include "Marmot/MarmotLowerDimensionalStress.h"
#include "Marmot/MarmotMath.h"
#include "Marmot/MarmotTensor.h"
#include "Marmot/MarmotVoigt.h"
#include <iostream>
#include "Fastor/Fastor.h"
#include "Marmot/interface_material_helper_functions.h"

using namespace Eigen;
using namespace Fastor;

using Tensor1D = Fastor::Tensor<double,3>;
using Tensor2D = Fastor::Tensor<double,3,3>;
using Tensor3D = Fastor::Tensor<double,3,3,3>;
using Tensor4D = Fastor::Tensor<double,3,3,3,3 >;

void MarmotMaterialHypoElasticInterface::setCharacteristicElementLength( double length )
{
  characteristicElementLength = length;
}

//void MarmotMaterialHypoElasticInterface::computeStress( Tensor1D&  force,
//                                                        Tensor2D&  surface_stress,
//                                                        Fastor::Tensor<double, 21,21>& dStress_dStrain,
//                                                        const Fastor::Tensor<double, 6,1>& dU,
//                                                        const Fastor::Tensor<double, 18,1>& dSurface_strain,
//                                                        const Tensor1D& normal,
//                                                        const double* timeOld,
//                                                        const double  dT,
//                                                        double&       pNewDT  )
// {};

void MarmotMaterialHypoElasticInterface::computeStress( double*  force,
                                                        double*  surface_stress,
                                                        double* dStress_dStrain,
                                                        const double* dU,
                                                        const double* dSurface_strain,
                                                        const double* normal,
                                                        const double* timeOld,
                                                        const double  dT,
                                                        double&       pNewDT)
{}//namespace Marmot::Materials





